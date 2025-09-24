#!/usr/bin/env python3
"""
Parse minimap2 PAF output into tabular format with percent identity and coverage.

Input:  PAF file (from minimap2 -x asm20)
Output: TSV with columns:
    query_id, target_id, percent_identity, query_coverage, aln_length, matches
"""

import sys

if len(sys.argv) < 3:
    print(f"Usage: {sys.argv[0]} input.paf output.tsv", file=sys.stderr)
    sys.exit(1)

infile, outfile = sys.argv[1], sys.argv[2]

with open(infile) as fin, open(outfile, "w") as fout:
    fout.write("query_id\ttarget_id\tpid\tqcov\taln_len\tmatches\n")
    for line in fin:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.strip().split("\t")
        qname, qlen = fields[0], int(fields[1])
        tname = fields[5]
        qstart, qend = int(fields[2]), int(fields[3])
        matches = int(fields[9])       # number of matching bases
        aln_len = int(fields[10])      # alignment block length

        # percent identity (matches / alignment length)
        pid = 100.0 * matches / aln_len if aln_len > 0 else 0.0

        # query coverage (alignment length / query length)
        qcov = 100.0 * aln_len / qlen if qlen > 0 else 0.0

        fout.write(
            f"{qname}\t{tname}\t{pid:.2f}\t{qcov:.2f}\t{aln_len}\t{matches}\n"
        )
