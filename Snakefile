# ==============================================================================
# RKN METAGENOMICS
# taxonomic assignment workflow for root-knot nematode minion sequencing data
# ==============================================================================

configfile: "config.yaml"

# ------------------------------------------------------------------------------
# SAMPLE LIST: load sample information from file
# ------------------------------------------------------------------------------

with open('samples.tsv') as infile:
    sample_list = []
    for line in infile:
        sample_list.append(line.strip())
SAMPLES = sample_list

# ------------------------------------------------------------------------------
# TARGET RULE
# ------------------------------------------------------------------------------

rule all:
    input:
        expand("results/recentrifuge/{sample}.html", sample = SAMPLES)

# ------------------------------------------------------------------------------
# MERGE FASTQ: merge all fastq files per barcode
# ------------------------------------------------------------------------------

rule merge_fastq:
    conda:
        "envs/tax.yaml"
    input:
        reads = config["data_dir"] + "/pass/{SAMPLES}"
    output:
        reads_merged = "results/merged_fastq/{SAMPLES}.merged.fastq"
    shell:
        "cat {input}/*runid*.fastq > {output}"

# ------------------------------------------------------------------------------
# FASTP: qc and trimming of reads
# ------------------------------------------------------------------------------

rule fastp_qc:
    conda:
        "envs/tax.yaml"
    input:
        reads = "results/merged_fastq/{SAMPLES}.merged.fastq"
    output:
        reads_qc = "results/qc/{SAMPLES}.qc.fastq",
        html = "results/qc/fastp_out/{SAMPLES}.fastp.html",
        json = "results/qc/fastp_out/{SAMPLES}.fastp.json",
        failed = "results/qc/failed_out/{SAMPLES}.failed.fastq"

    shell:
        "fastp -i {input.reads} \
        -o {output.reads_qc} \
        -A \
        -f 30 \
        -q 15 \
        -l 2000 \
        --failed_out {output.failed} \
        -h {output.html} \
        -j {output.json}"

# ------------------------------------------------------------------------------
# KRAKEN2: taxonomic assignment
# ------------------------------------------------------------------------------

rule kraken2:
    conda:
        "envs/tax.yaml"
    input:
        reads = "results/qc/{SAMPLES}.qc.fastq"
    output:
        reports = "results/kraken2/reports/{SAMPLES}.txt",
        outputs = "results/kraken2/outputs/{SAMPLES}.tsv"
    threads:
        10
    shell:
        "kraken2 --db {config[kraken2_db]} {input.reads} \
        --minimum-base-quality 10 \
        --use-names \
        --threads {threads} \
        --confidence 0.01 \
        --report {output.reports} \
        --output {output.outputs}"

#  --memory-mapping

# ------------------------------------------------------------------------------
# RECENTRIFUGE: generating taxonomic level assignment html figure
# ------------------------------------------------------------------------------

rule recentrifuge:
    conda:
        "envs/tax.yaml"
    input:
        taxdb = config["taxdump"],
        kraken_output = "results/kraken2/outputs/{SAMPLES}.tsv"
    output:
        "results/recentrifuge/{SAMPLES}.html"
    shell:
        "rcf -n {input.taxdb} -k {input.kraken_output} -o {output}"
