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
parser.add_argument('-kt', '--keep-tables', action='store_true', help='If included, will keep the filtered metadata, multifasta and swell files.')
parser.add_argument('-rpdf', '--make-run-pdfs', action='store_true', help='If included, will also generate PDFs for each sequencing org, displaying plots for individual runs.')
args = parser.parse_args()
start = time.time()


# Make results directory (if it doesn't already exist)
date_prefix = f"{args.start_date}_{args.end_date}"
results_dir = date_prefix
if not os.path.isdir(results_dir):
    subprocess.run(['mkdir', results_dir])


sequencing_org_codes = []
filtered_metadata_path = f"{results_dir}/{date_prefix}_metadata.tsv"
if not os.path.exists(filtered_metadata_path):
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
    # Filter metadata
    filtered_metadata = metadata.loc[mask]
    filtered_metadata.to_csv(filtered_metadata_path, index=False, sep='\t')
    # Grab sequencing org codes for potential future use
    sequencing_org_codes = filtered_metadata["sequencing_org_code"].unique()
    print("done.")


filtered_multifasta_path = f"{results_dir}/{date_prefix}_multifasta.fasta"
if not os.path.exists(filtered_multifasta_path):
    print("Filtering multifasta by the provided start/end dates...", end=" ", flush=True)
    with open(args.multifasta) as multifasta, open(filtered_multifasta_path, "w") as filtered_multifasta:
        heng_iter = readfq(multifasta)
        for name, seq, qual in heng_iter:
            header = (name.split('|'))[0]
            if header_filter[header]:
                filtered_multifasta.write(f">{header}\n{seq}\n")
    print("done.")


swell_filtered_multifasta_path = f"{results_dir}/{date_prefix}_swell_multifasta.tsv"
if not os.path.exists(swell_filtered_multifasta_path):
    print("Running swell on filtered multifasta...", end=" ", flush=True)
    with open(swell_filtered_multifasta_path, "w") as swell_multifasta:
        subprocess.run(['swell', 'separate-fasta', filtered_multifasta_path], stdout=swell_multifasta)
    print("done.")


print("Generating report...")
report_dest = f"{results_dir}/{args.start_date}_{args.end_date}_report.pdf"
date_string = f"{args.start_date[6:]}/{args.start_date[4:6]}/{args.start_date[0:4]} - {args.end_date[6:]}/{args.end_date[4:6]}/{args.end_date[0:4]}"
subprocess.run(['Rscript', '-e', f'rmarkdown::render(\'report.Rmd\', output_file = "{report_dest}", params=list(swell_multifasta = \'{swell_filtered_multifasta_path}\', date_filtered_metadata = \'{filtered_metadata_path}\', start_date = \'{args.start_date}\', end_date = \'{args.end_date}\', date_string = \'{date_string}\', results_dir = \'{results_dir}\'))'])
print("done.")


if args.make_run_pdfs:
    print("Generating run PDFs...")
    # Get list of sequencing orgs (if not already obtained)
    if len(sequencing_org_codes) == 0:
        filtered_metadata = pd.read_csv(filtered_metadata_path, sep="\t", low_memory=False)
        sequencing_org_codes = filtered_metadata["sequencing_org_code"].unique()
    # For each code, make a directory (if it doesn't already exist)
    # Then, generate the pdf for that sequencing org
    for code in sequencing_org_codes:
        org_dir = f"{results_dir}/{code}"
        if not os.path.isdir(org_dir):
            subprocess.run(['mkdir', org_dir])
        runs_dest = f"{org_dir}/{args.start_date}_{args.end_date}_{code}.pdf"
        subprocess.run(['Rscript', '-e', f'rmarkdown::render(\'runs.Rmd\', output_file = "{runs_dest}", params=list(code = \'{code}\', swell_multifasta = \'{swell_filtered_multifasta_path}\', date_filtered_metadata = \'{filtered_metadata_path}\', date_prefix = \'{date_prefix}\', date_string = \'{date_string}\', org_dir = \'{org_dir}\'))'])
    # Generate the PDFs
    # subprocess.run(['Rscript', 'runs.R', code, swell_filtered_multifasta_path, filtered_metadata_path, org_dir, date_prefix])
    print("done.")


if not args.keep_tables:
    print("Removing filtered tables...", end=" ", flush=True)
    subprocess.run(['rm', filtered_metadata_path, filtered_multifasta_path, swell_filtered_multifasta_path])
    print("done.")


end = time.time()
print(f"Time taken: {int((end - start) // 60)} m {int((end - start) % 60)} s")