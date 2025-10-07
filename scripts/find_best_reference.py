#!/usr/bin/env python3
"""
find_best_reference.py
----------------------------------------------------------
Identify the best reference genome for metallophore-rescue candidate contigs.
- Uses bin assignments when available (ANI search vs reference DB)
- Falls back to direct nucleotide alignment (minimap2 or blastn)
- Produces two outputs:
    1) best_reference.tsv : one row per contig, best reference accession
    2) best_reference.paf : optional alignment summary
----------------------------------------------------------
"""

import argparse
import pandas as pd
import subprocess
import sys
from pathlib import Path

def run_cmd(cmd):
    """Run a shell command safely."""
    print(f"[CMD] {cmd}", file=sys.stderr)
    subprocess.run(cmd, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description="Find best reference match for metallophore candidates.")
    parser.add_argument("--contigs", required=True, help="FASTA of candidate contigs.")
    parser.add_argument("--bin-map", required=True, help="TSV mapping contig -> bin.")
    parser.add_argument("--ref-db", required=True, help="Path to reference database (FASTA or index).")
    parser.add_argument("--out-best", required=True, help="Output: best reference table.")
    parser.add_argument("--out-align", required=True, help="Output: alignment summary (PAF).")
    parser.add_argument("--min-align-frac", type=float, default=0.5, help="Minimum fraction of contig aligned.")
    args = parser.parse_args()

    contigs = Path(args.contigs)
    bin_map = pd.read_csv(args.bin_map, sep="\t", header=None, names=["contig", "bin"])
    out_best = Path(args.out_best)
    out_align = Path(args.out_align)

    # Try to see if contigs have bin info
    has_bins = bin_map["bin"].notna().any()
    best_hits = []

    if has_bins:
        print("Using bin-level ANI search for matching references...", file=sys.stderr)
        # Example: use skani or fastANI to find nearest reference for each bin
        bins = bin_map["bin"].unique()
        for b in bins:
            tmp_out = Path(f"{b}_ani.tsv")
            cmd = f"skani dist -q bins/{b}.fa -r {args.ref_db} --threads 8 > {tmp_out}"
            run_cmd(cmd)
            # Parse and take top hit
            df = pd.read_csv(tmp_out, sep="\t", header=None)
            if not df.empty:
                ref = df.iloc[df[2].idxmin(), 1]
                ani = 100 - df.iloc[df[2].idxmin(), 2]
                best_hits.append({"bin": b, "best_ref": ref, "ani": ani})
    else:
        print("No bin information found; using direct minimap2 alignment...", file=sys.stderr)
        cmd = f"minimap2 -x asm10 -t 8 {args.ref_db} {contigs} > {out_align}"
        run_cmd(cmd)
        # Simple parser for PAF
        paf = pd.read_csv(out_align, sep="\t", header=None)
        paf = paf[[0,5,9,10]]
        paf.columns = ["contig", "ref", "matches", "aln_len"]
        paf["score"] = paf["matches"] * paf["aln_len"]
        paf = paf.sort_values("score", ascending=False).drop_duplicates("contig")
        paf.to_csv(out_best, sep="\t", index=False)
        sys.exit(0)

    pd.DataFrame(best_hits).to_csv(out_best, sep="\t", index=False)
    print(f"✅ Wrote {out_best} and {out_align}", file=sys.stderr)

if __name__ == "__main__":
    main()
