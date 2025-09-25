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
    """
    Top-level target for metallophore BGC rescue workflow.
    """
    input:
        expand("results/metallophore_rescue/{sample}.uhgg.filtered.tsv", sample=Samples)

# ---------------------------------------------------------------------------
# Rule: align_uhgg
# ---------------------------------------------------------------------------
rule align_uhgg:
    """
    Align assembled contigs against UHGG reference with minimap2.
    Produces PAF output for downstream parsing.
    """
    input:
        contigs=ASSEMBLY_FP / "{sample}.fasta",
        db=Cfg["sbx_metallophore_rescue"]["references"]["uhgg_mmi"]
    output:
        paf="results/metallophore_rescue/{sample}.uhgg.paf"
    log:
        LOG_FP / "align_uhgg_{sample}.log"
    benchmark:
        BENCHMARK_FP / "align_uhgg_{sample}.tsv"
    threads: 8
    conda:
        "envs/sbx_metallophore_env.yml"
    shell:
        """
        minimap2 -x asm20 -t {threads} {input.db} {input.contigs} \
          > {output.paf} 2> {log}
        """

# ---------------------------------------------------------------------------
# Rule: parse_alignments
# ---------------------------------------------------------------------------
rule parse_alignments:
    """
    Parse PAF output from minimap2 to extract percent identity and coverage.
    """
    input:
        "results/metallophore_rescue/{sample}.uhgg.paf"
    output:
        "results/metallophore_rescue/{sample}.uhgg.parsed.tsv"
    log:
        LOG_FP / "parse_alignments_{sample}.log"
    conda:
        "envs/sbx_metallophore_env.yml"
    script:
        "scripts/parse_paf.py"

# ---------------------------------------------------------------------------
# Rule: filter_alignments
# ---------------------------------------------------------------------------
rule filter_alignments:
    """
    Filter parsed alignments by identity and query coverage thresholds.
    """
    input:
        "results/metallophore_rescue/{sample}.uhgg.parsed.tsv"
    output:
        "results/metallophore_rescue/{sample}.uhgg.filtered.tsv"
    log:
        LOG_FP / "filter_alignments_{sample}.log"
    params:
        pid_min=Cfg["sbx_metallophore_rescue"]["alignment"]["pid_min"],
        qcov_min=Cfg["sbx_metallophore_rescue"]["alignment"]["qcov_min"]
    conda:
        "envs/sbx_metallophore_env.yml"
    shell:
        """
        awk -v pid={params.pid_min} -v qcov={params.qcov_min} \
          '$2 >= pid && $3 >= qcov {{print}}' {input} > {output} 2> {log}
        """

# ---------------------------------------------------------------------------
# Rule: run_antismash (placeholder)
# ---------------------------------------------------------------------------
rule run_antismash:
    """
    Run antiSMASH on candidate genomes/regions identified by UHGG rescue.
    Placeholder rule — to be expanded.
    """
    input:
        "results/metallophore_rescue/{sample}.uhgg.filtered.tsv"
    output:
        directory("results/metallophore_rescue/antismash/{sample}")
    log:
        LOG_FP / "antismash_{sample}.log"
    conda:
        "envs/sbx_metallophore_env.yml"
    shell:
        """
        echo "antiSMASH would run here using the rescued genome regions." > {log}
        mkdir -p {output}
        """

###############################################################################
# End of sbx_metallophore_rescue.smk
###############################################################################
