###############################################################################
# sbx_metallophore_rescue.smk
# Sunbeam extension for metallophore BGC rescue
###############################################################################

# ---------------------------------------------------------------------------
# Extension metadata and logging
# ---------------------------------------------------------------------------
try:
    SBX_METALLOPHORE_RESCUE_VERSION = get_ext_version("sbx_metallophore_rescue")
except NameError:
    SBX_METALLOPHORE_RESCUE_VERSION = "0.0.1"

try:
    logger = get_extension_logger("sbx_metallophore_rescue")
except NameError:
    import logging
    logger = logging.getLogger("sunbeam.pipeline.extensions.sbx_metallophore_rescue")

logger.info("Loading sbx_metallophore_rescue extension...")
logger.info(f"Using version {SBX_METALLOPHORE_RESCUE_VERSION}")

# ---------------------------------------------------------------------------
# Local rules
# ---------------------------------------------------------------------------
localrules: all_metallophores

# ---------------------------------------------------------------------------
# Rule: all_metallophores (pipeline target)
# ---------------------------------------------------------------------------
rule all_metallophores:
    input:
        expand("results/metallophore_rescue/{sample}.validated.tsv", sample=Samples)


###############################################################################
# Rules
###############################################################################

# 1. Run antiSMASH on assemblies
rule run_antismash_assembly:
    input:
        assembly=ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa",
    output:
        directory("results/antismash/{sample}")
    log:
        "logs/run_antismash_assembly_{sample}.log"
        threads: 40
    params:
        datapath=Cfg["sbx_metallophore_rescue"]["antismash_datapath"]
    container:
        "docker://antismash/stanalone:8.0.2"
    shell:
        """
        antismash \
          --databases {params.datapath} \
          --genefinding-tool prodigal-m \
          --output-dir {output} \
          --output-basename {wildcards.sample} \
          --taxon bacteria \
          -c {threads} \
          --cb-knownclusters --cc-mibig --tfbs --fimo --no-enable-genefunctions \
          {input.assembly} \
          > {log} 2>&1
        """

# 2. Parse antiSMASH results
rule parse_antismash_assembly:
    input:
        "results/antismash/{sample}"
    output:
        "results/metallophore_rescue/{sample}.assembly_parsed.tsv"
    log:
        "logs/parse_antismash_assembly_{sample}.log"

# 3. Run FeGenie (optional branch)
rule fegenie_scan:
    input:
        assembly=ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa",
    output:
        "results/metallophore_rescue/{sample}.fegenie.tsv"
    log:
        "logs/fegenie_{sample}.log"

# 4. Select contigs for rescue (merge AntiSMASH + FeGenie evidence)
rule select_rescue_candidates:
    input:
        antismash="results/metallophore_rescue/{sample}.assembly_parsed.tsv",
        fegenie="results/metallophore_rescue/{sample}.fegenie.tsv"
    output:
        "results/metallophore_rescue/{sample}.rescue_candidates.fasta"
    log:
        "logs/select_rescue_candidates_{sample}.log"

# 5. Align candidate contigs to UHGG
rule align_uhgg:
    input:
        contigs="results/metallophore_rescue/{sample}.rescue_candidates.fasta",
        db=Cfg["sbx_metallophore_rescue"]["references"]["uhgg_mmi"]
    output:
        "results/metallophore_rescue/{sample}.uhgg.paf"
    log:
        "logs/align_uhgg_{sample}.log"

# 6. Parse alignments
rule parse_alignments:
    input:
        "results/metallophore_rescue/{sample}.uhgg.paf"
    output:
        "results/metallophore_rescue/{sample}.uhgg_best.tsv"
    log:
        "logs/parse_alignments_{sample}.log"

# 7. Run antiSMASH on rescued MAG regions
rule run_antismash_mag:
    input:
        "results/metallophore_rescue/{sample}.uhgg_best.tsv"
    output:
        directory("results/metallophore_rescue/antismash_mag/{sample}")
    log:
        "logs/run_antismash_mag_{sample}.log"

# 8. Parse MAG antiSMASH results for validation
rule parse_antismash_mag:
    input:
        antismash="results/metallophore_rescue/antismash_mag/{sample}"
    output:
        "results/metallophore_rescue/{sample}.validated.tsv"
    log:
        "logs/parse_antismash_mag_{sample}.log"
###############################################################################
# End of sbx_metallophore_rescue.smk
###############################################################################
