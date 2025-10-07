#!/usr/bin/env python3
"""
integrate_confirmations.py
----------------------------------------------------------
Integrates discovery (triage), best reference, and confirm/deny evidence.
- Merges antiSMASH + FeGenie confirmations
- Labels contigs as CONFIRM / DENY / UNCERTAIN
Output: metallophore_calls.tsv
----------------------------------------------------------
"""

import argparse
import pandas as pd
from pathlib import Path

def confirm_status(row):
    """Simple confirmation logic."""
    as_hit = row.get("antismash_hit", False)
    fe_hit = row.get("fegenie_hit", False)
    if as_hit and fe_hit:
        return "CONFIRM"
    elif as_hit or fe_hit:
        return "UNCERTAIN"
    else:
        return "DENY"

def main():
    parser = argparse.ArgumentParser(description="Integrate confirmation results for metallophore rescue.")
    parser.add_argument("--triage", required=True, help="Triage candidate TSV")
    parser.add_argument("--bestref", required=True, help="Best reference TSV")
    parser.add_argument("--antismash", required=True, help="antiSMASH confirmation directory")
    parser.add_argument("--fegenie", required=True, help="FeGenie confirmation TSV")
    parser.add_argument("--out", required=True, help="Output TSV with metallophore calls")
    args = parser.parse_args()

    triage = pd.read_csv(args.triage, sep="\t")
    bestref = pd.read_csv(args.bestref, sep="\t")
    fege = pd.read_csv(args.fegenie, sep="\t")

    # Parse antiSMASH results summary file if available
    antismash_summary = Path(args.antismash) / "index.html"
    if antismash_summary.exists():
        anti_df = pd.read_html(antismash_summary)[0]
    else:
        anti_df = pd.DataFrame(columns=["cluster_id", "product"])
    anti_df["antismash_hit"] = anti_df["product"].str.contains("siderophore|nrps|metal", case=False, na=False)

    # Merge evidence
    merged = triage.merge(bestref, how="left", left_on="contig", right_on="contig", suffixes=("", "_ref"))
    merged = merged.merge(fege, how="left", left_on="contig", right_on="contig")
    merged = merged.merge(anti_df, how="left", left_on="contig", right_on="cluster_id")

    merged["fegenie_hit"] = merged.get("Function", "").str.contains("siderophore|iron|zinc", case=False, na=False)
    merged["confirm_status"] = merged.apply(confirm_status, axis=1)

    cols = [
        "sample", "contig", "best_ref", "confirm_status",
        "antismash_hit", "fegenie_hit", "product", "ani"
    ]
    merged[cols].to_csv(args.out, sep="\t", index=False)
    print(f"✅ Wrote {args.out}")

if __name__ == "__main__":
    main()
