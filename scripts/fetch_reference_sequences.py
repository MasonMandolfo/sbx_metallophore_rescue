#!/usr/bin/env python3
"""
fetch_reference_sequences.py
----------------------------------------------------------
Given a table of best reference accessions, fetch corresponding sequences.
If the reference exists locally, symlink or copy it.
If remote, use NCBI Entrez or wget/rsync to retrieve.
Output: combined reference_sequences.fna
----------------------------------------------------------
"""

import argparse
import pandas as pd
from pathlib import Path
import subprocess
import sys

def fetch_ncbi(accession, outdir):
    """Fetch sequence by accession from NCBI using efetch or wget."""
    outfile = Path(outdir) / f"{accession}.fna"
    cmd = f"datasets download genome accession {accession} --include genome --filename {outfile}.zip"
    subprocess.run(cmd, shell=True, check=True)
    return outfile

def main():
    parser = argparse.ArgumentParser(description="Fetch reference FASTA sequences for confirmation runs.")
    parser.add_argument("--input", required=True, help="TSV of best_reference.tsv")
    parser.add_argument("--out", required=True, help="Output combined reference FASTA")
    args = parser.parse_args()

    ref_table = pd.read_csv(args.input, sep="\t")
    outdir = Path(args.out).parent
    outdir.mkdir(parents=True, exist_ok=True)

    all_fasta_paths = []
    for _, row in ref_table.iterrows():
        acc = str(row.get("best_ref") or row.get("ref") or row[0])
        print(f"Fetching reference {acc}", file=sys.stderr)
        try:
            fasta = fetch_ncbi(acc, outdir)
        except Exception:
            print(f"⚠️ Failed to fetch {acc}, skipping.", file=sys.stderr)
            continue
        all_fasta_paths.append(fasta)

    # Merge into one combined FASTA
    with open(args.out, "w") as out_fh:
        for f in all_fasta_paths:
            subprocess.run(f"unzip -p {f}.zip '*.fna' >> {args.out}", shell=True)
    print(f"✅ Combined references written to {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
