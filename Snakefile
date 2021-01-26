# ==============================================================================
# RKN METAGENOMICS
# taxonomic assignment workflow for root-knot nematode minion sequencing data
# ==============================================================================

configfile: "config.yaml"

# ------------------------------------------------------------------------------
# SAMPLE LIST: load sample information from file
# ------------------------------------------------------------------------------

with open(config['sample_list'], 'r') as infile:
    sample_list = []
    for line in infile:
        sample_list.append(line.strip())
SAMPLES = sample_list

# ------------------------------------------------------------------------------
# TARGET RULE
# ------------------------------------------------------------------------------

rule all:
    input:
        expand("results/kraken2/outputs/{sample}.krk", sample = SAMPLES),
        expand("results/kraken2/reports/{sample}.txt", sample = SAMPLES),
        "results/recentrifuge/run.html",
        "results/recentrifuge/run.xlsx"

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
# NANOFILT: qc and trimming of reads
# ------------------------------------------------------------------------------

rule nanofilt:
    conda:
        "envs/tax.yaml"
    input:
        reads = "results/merged_fastq/{SAMPLES}.merged.fastq"
    output:
        nano_out = "results/qc/{SAMPLES}.qc.fastq"
    shell:
        "NanoFilt {input.reads} \
        -q {config[qual]} \
        -l {config[len]} \
        --headcrop {config[head]} \
        --tailcrop {config[tail]} \
        > {output.nano_out}"

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
        outputs = "results/kraken2/outputs/{SAMPLES}.krk"
    threads:
        10
    shell:
        "kraken2 --db {config[kraken2_db]} {input.reads} \
        --minimum-base-quality {config[minq]} \
        --use-names \
        --threads {threads} \
        --confidence {config[conf]} \
        --report {output.reports} \
        --output {output.outputs}"

# ------------------------------------------------------------------------------
# RECENTRIFUGE: generating taxonomic level assignment html figure
# ------------------------------------------------------------------------------

rule recentrifuge:
    conda:
        "envs/tax.yaml"
    input:
        taxdb = config["taxdump"],
        kraken_output = expand("results/kraken2/outputs/{sample}.krk", sample = SAMPLES)
    params:
        directory("results/kraken2/outputs/")
    output:
        re_html = "results/recentrifuge/run.html",
        re_xlsx = "results/recentrifuge/run.xlsx"
    shell:
        "rcf -a \
        -n {input.taxdb} \
        -k {params} \
        -o {output.re_html}"
