###############################################################################
# sbx_metallophore_rescue.smk
# Sunbeam extension for metallophore BGC rescue
###############################################################################

# ---------------------------------------------------------------------------
# Metadata & Logging
# ---------------------------------------------------------------------------
try:
    SBX_METALLOPHORE_RESCUE_VERSION = get_ext_version("sbx_metallophore_rescue")
except NameError:
    SBX_METALLOPHORE_RESCUE_VERSION = "0.1.0"

try:
    logger = get_extension_logger("sbx_metallophore_rescue")
except NameError:
    import logging
    logger = logging.getLogger("sunbeam.pipeline.extensions.sbx_metallophore_rescue")

logger.info(f"Loading sbx_metallophore_rescue v{SBX_METALLOPHORE_RESCUE_VERSION}")

# ---------------------------------------------------------------------------
# Local Rules
# ---------------------------------------------------------------------------
localrules: all_metallophores

# ---------------------------------------------------------------------------
# Final Target
# ---------------------------------------------------------------------------
rule all_metallophores:
    """Final output: validated metallophore BGC calls per sample"""
    input:
        expand("results/metallophore_rescue/{sample}.metallophore_calls.tsv", sample=Samples)

###############################################################################
# PIPELINE STRUCTURE
###############################################################################
# Phase 1: Lenient Discovery
# Phase 2: Rescue by Best Reference (taxonomic match)
# Phase 3: Confirmation/Reannotation
###############################################################################

# ---------------------------------------------------------------------------
# 1️⃣ Run antiSMASH (lenient discovery)
# ---------------------------------------------------------------------------
rule antismash_discovery:
    input:
        assembly = ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa"
    output:
        directory("results/metallophore_rescue/antismash_discovery/{sample}")
    log:
        "logs/antismash_discovery_{sample}.log"
    threads: 32
    params:
        db = Cfg["sbx_metallophore_rescue"]["antismash"]["db"]
    conda:
        "envs/antismash.yml"
    shell:
        r"""
        antismash \
          --databases {params.db} \
          --taxon bacteria \
          --genefinding-tool prodigal-m \
          --cb-knownclusters --cc-mibig --tfbs --fimo --no-enable-genefunctions \
          --clusterhmmer --asf --minimal --fullhmmer \
          -c {threads} \
          {input.assembly} > {log} 2>&1
        """

# ---------------------------------------------------------------------------
# 2️⃣ Parse antiSMASH discovery results
# ---------------------------------------------------------------------------
rule parse_antismash_discovery:
    input:
        "results/metallophore_rescue/antismash_discovery/{sample}"
    output:
        "results/metallophore_rescue/{sample}.antismash_discovery.tsv"
    log:
        "logs/parse_antismash_discovery_{sample}.log"
    conda:
        "envs/sbx_metallophore.yml"
    shell:
        "python {CFG_DIR}/scripts/antismash2csv.py -o {output} {input} > {log} 2>&1"

# ---------------------------------------------------------------------------
# 3️⃣ Run FeGenie (low stringency)
# ---------------------------------------------------------------------------
rule fegenie_discovery:
    input:
        assembly = ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa"
    output:
        "results/metallophore_rescue/{sample}.fegenie_discovery.tsv"
    log:
        "logs/fegenie_discovery_{sample}.log"
    conda:
        "envs/fegenie.yml"
    params:
        stringency = "low"
    shell:
        r"""
        FeGenie.py -bin_dir $(dirname {input.assembly}) \
                   -bin_ext fa \
                   -out results/metallophore_rescue/fegenie_{wildcards.sample} \
                   -stringent {params.stringency} \
                   > {log} 2>&1
        mv results/metallophore_rescue/fegenie_{wildcards.sample}/FeGenie_results_summary.txt {output}
        """

# ---------------------------------------------------------------------------
# 4️⃣ Identify candidate contigs (triage step)
# ---------------------------------------------------------------------------
rule triage_candidates:
    input:
        antismash = "results/metallophore_rescue/{sample}.antismash_discovery.tsv",
        fegenie   = "results/metallophore_rescue/{sample}.fegenie_discovery.tsv"
    output:
        fasta     = "results/metallophore_rescue/{sample}.candidates.fna",
        metadata  = "results/metallophore_rescue/{sample}.candidates.tsv"
    log:
        "logs/triage_candidates_{sample}.log"
    conda:
        "envs/sbx_metallophore.yml"
    shell:
        "python {CFG_DIR}/scripts/triage.py --antismash {input.antismash} --fegenie {input.fegenie} "
        "--out-fasta {output.fasta} --out-tsv {output.metadata} > {log} 2>&1"

# ---------------------------------------------------------------------------
# 5️⃣ Rescue by Best Reference
#     (Prefer bin-based ANI search, else direct nucleotide alignment)
# ---------------------------------------------------------------------------
rule rescue_by_reference:
    input:
        contigs = "results/metallophore_rescue/{sample}.candidates.fna",
        bin_map = "qc/{sample}/refined/{sample}.magscot.refined.contig_to_bin.out",
        ref_db  = Cfg["sbx_metallophore_rescue"]["references"]["refseq"]
    output:
        besthits = "results/metallophore_rescue/{sample}.best_reference.tsv",
        align    = "results/metallophore_rescue/{sample}.best_reference.paf"
    log:
        "logs/rescue_by_reference_{sample}.log"
    threads: 16
    params:
        min_align_frac = 0.5
    conda:
        "envs/mmseqs.yml"
    shell:
        r"""
        python {CFG_DIR}/scripts/find_best_reference.py \
          --contigs {input.contigs} \
          --bin-map {input.bin_map} \
          --ref-db {input.ref_db} \
          --out-best {output.besthits} \
          --out-align {output.align} \
          --min-align-frac {params.min_align_frac} \
          > {log} 2>&1
        """

# ---------------------------------------------------------------------------
# 6️⃣ Fetch rescued references for confirmation
# ---------------------------------------------------------------------------
rule fetch_reference_sequences:
    input:
        "results/metallophore_rescue/{sample}.best_reference.tsv"
    output:
        "results/metallophore_rescue/{sample}.reference_sequences.fna"
    log:
        "logs/fetch_reference_sequences_{sample}.log"
    conda:
        "envs/sbx_metallophore.yml"
    shell:
        "python {CFG_DIR}/scripts/fetch_reference_sequences.py --input {input} --out {output} > {log} 2>&1"

# ---------------------------------------------------------------------------
# 7️⃣ Confirmation antiSMASH + FeGenie pass on rescued reference
# ---------------------------------------------------------------------------
rule antismash_confirm:
    input:
        "results/metallophore_rescue/{sample}.reference_sequences.fna"
    output:
        directory("results/metallophore_rescue/antismash_confirm/{sample}")
    log:
        "logs/antismash_confirm_{sample}.log"
    threads: 32
    params:
        db = Cfg["sbx_metallophore_rescue"]["antismash"]["db"]
    conda:
        "envs/antismash.yml"
    shell:
        r"""
        antismash \
          --databases {params.db} \
          --taxon bacteria \
          --genefinding-tool prodigal-m \
          --cb-knownclusters --cc-mibig --tfbs --fimo --no-enable-genefunctions \
          -c {threads} \
          {input} > {log} 2>&1
        """

rule fegenie_confirm:
    input:
        "results/metallophore_rescue/{sample}.reference_sequences.fna"
    output:
        "results/metallophore_rescue/{sample}.fegenie_confirm.tsv"
    log:
        "logs/fegenie_confirm_{sample}.log"
    conda:
        "envs/fegenie.yml"
    params:
        stringency = "medium"
    shell:
        r"""
        FeGenie.py -bin_dir $(dirname {input}) \
                   -bin_ext fna \
                   -out results/metallophore_rescue/fegenie_confirm_{wildcards.sample} \
                   -stringent {params.stringency} \
                   > {log} 2>&1
        mv results/metallophore_rescue/fegenie_confirm_{wildcards.sample}/FeGenie_results_summary.txt {output}
        """

# ---------------------------------------------------------------------------
# 8️⃣ Integrate confirmation evidence (final metallophore call)
# ---------------------------------------------------------------------------
rule integrate_results:
    input:
        triage  = "results/metallophore_rescue/{sample}.candidates.tsv",
        bestref = "results/metallophore_rescue/{sample}.best_reference.tsv",
        antismash_confirm = directory("results/metallophore_rescue/antismash_confirm/{sample}"),
        fegenie_confirm   = "results/metallophore_rescue/{sample}.fegenie_confirm.tsv"
    output:
        "results/metallophore_rescue/{sample}.metallophore_calls.tsv"
    log:
        "logs/integrate_results_{sample}.log"
    conda:
        "envs/sbx_metallophore.yml"
    shell:
        r"""
        python {CFG_DIR}/scripts/integrate_confirmations.py \
            --triage {input.triage} \
            --bestref {input.bestref} \
            --antismash {input.antismash_confirm} \
            --fegenie {input.fegenie_confirm} \
            --out {output} > {log} 2>&1
        """

###############################################################################
# End of sbx_metallophore_rescue.smk
###############################################################################
