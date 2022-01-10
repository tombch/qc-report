import os
import argparse
import subprocess
import pandas as pd
import time


# This function was written by Heng Li
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


parser = argparse.ArgumentParser()
parser.add_argument('-md', '--metadata', required=True, help='Path to metadata.')
parser.add_argument('-mf', '--multifasta', required=True, help='Path to multifasta.')
parser.add_argument('-s', '--start-date', required=True, help='Minimum sequencing_submission_date.')
parser.add_argument('-e', '--end-date', required=True, help='Maximum sequencing_submission_date.')
parser.add_argument('-r', '--results-dir', help='Directory to save generated report and PNG files (default will be the directory STARTDATE_ENDDATE).')
parser.add_argument('-t', '--temp-dir', default='.', help='Directory to save temporary files (default being the current directory).')
args = parser.parse_args()
start = time.time()


# File prefix and temp file names
prefix = f"{args.start_date}_{args.end_date}"
filtered_metadata_path = f"{args.temp_dir}/{prefix}_metadata.tsv"
filtered_multifasta_path = f"{args.temp_dir}/{prefix}_multifasta.fasta"
swell_filtered_multifasta_path = f"{args.temp_dir}/{prefix}_swell_multifasta.tsv"


# Make results directory if needed
results_dir = args.results_dir
if not results_dir:
    results_dir = prefix
if not os.path.isdir(results_dir):
    subprocess.run(['mkdir', results_dir])


print("Reading metadata...", end=" ", flush=True)
metadata = pd.read_csv(args.metadata, sep="\t", low_memory=False)
metadata["sequencing_submission_date"] = pd.to_datetime(metadata["sequencing_submission_date"], errors="coerce")
print("done.")


print("Filtering metadata by the provided start/end dates...", end=" ", flush=True)
start_date = pd.Timestamp(args.start_date)
end_date = pd.Timestamp(args.end_date)
# Create dictionary that records whether headers are within the date range
header_filter = {}
for i, (index, row) in enumerate(metadata.iterrows()):
    if row['sequencing_submission_date'] >= start_date and row['sequencing_submission_date'] <= end_date:
        header_filter[row['fasta_header']] = True
    else:
        header_filter[row['fasta_header']] = False
# Use a mask to retrieve and store rows within the date range
mask = (metadata['sequencing_submission_date'] >= args.start_date) & (metadata['sequencing_submission_date'] <= args.end_date)
metadata = metadata.loc[mask]
metadata.to_csv(filtered_metadata_path, index=False, sep='\t')
print("done.")


print("Filtering multifasta by the provided start/end dates...", end=" ", flush=True)
with open(args.multifasta) as multifasta, open(filtered_multifasta_path, "w") as filtered_multifasta:
    heng_iter = readfq(multifasta)
    for name, seq, qual in heng_iter:
        header = (name.split('|'))[0]
        if header_filter[header]:
            filtered_multifasta.write(f">{header}\n{seq}\n")
print("done.")


print("Running swell on filtered multifasta...", end=" ", flush=True)
with open(swell_filtered_multifasta_path, "w") as swell_multifasta:
    subprocess.run(['swell', 'separate-fasta', filtered_multifasta_path], stdout=swell_multifasta)
print("done.")


print("Generating report...")
report_dest = f"{results_dir}/{args.start_date}_{args.end_date}_report.pdf"
date_string = f"{args.start_date[6:]}/{args.start_date[4:6]}/{args.start_date[0:4]} - {args.end_date[6:]}/{args.end_date[4:6]}/{args.end_date[0:4]}"
subprocess.run(['Rscript', '-e', f'rmarkdown::render(\'report.Rmd\', output_file = "{report_dest}", params=list(swell_multifasta = \'{swell_filtered_multifasta_path}\', date_filtered_metadata = \'{filtered_metadata_path}\', start_date = \'{args.start_date}\', end_date = \'{args.end_date}\', date_string = \'{date_string}\', results_dir = \'{results_dir}\'))'])
print("done.")


print("Removing temporary files...", end=" ", flush=True)
subprocess.run(['rm', filtered_metadata_path, filtered_multifasta_path, swell_filtered_multifasta_path])
print("done.")


end = time.time()
print(f"Time taken: {int((end - start) // 60)} m {int((end - start) % 60)} s")