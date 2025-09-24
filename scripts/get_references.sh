#!/usr/bin/env bash
set -euo pipefail

mkdir -p references
cd references/

echo "Downloading MIBiG..."
wget -O references/mibig.fna.gz https://dl.secondarymetabolites.org/mibig/mibig_prot_seqs_3.1.fasta.gz
gunzip -f references/mibig.fna.gz
makeblastdb -in references/mibig.fna -dbtype nucl -out references/mibig_blastdb
minimap2 -d mibig.mmi references/mibig.fna

echo "Downloading UHGG (ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/UHGG_representative_genomes.tar.gz)..."
wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/UHGG_representative_genomes.tar.gz
tar -xvzf UHGG_representative_genomes.tar.gz
cat UHGG_representative_genomes/*.fa > UHGG_representative_genomes.fna
makeblastdb -in UHGG_representative_genomes.fna -dbtype nucl -out UHGG_rep_blastdb
minimap2 -d UHGG_rep.mmi UHGG_representative_genomes.fna
