# ==============================================================================
# RKN GENOMIC TAXONOMIC ASSIGNMENT
# taxonomic assignment workflow for root-knot nematode Flongle sequencing data
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
# OPTIONAL INPUT TYPES: guppy or fastq
# ------------------------------------------------------------------------------

if config["input_type"] == "guppy":
    ruleorder: guppy_input > fastq_input
if config["input_type"] == "fastq":
    ruleorder: fastq_input > guppy_input
else:
    print("add input file type to config")

# ------------------------------------------------------------------------------
# TARGET RULE
# ------------------------------------------------------------------------------

rule all:
    input:
        expand("results/kraken2/outputs/{sample}.krk", sample = SAMPLES),
        expand("results/kraken2/reports/{sample}.txt", sample = SAMPLES),
        "results/recentrifuge/" + config['my_experiment'] + ".html",
        "results/recentrifuge/" + config['my_experiment'] + ".xlsx"

# ------------------------------------------------------------------------------
# GUPPY INPUT: merge all fastq files per barcode directory
# ------------------------------------------------------------------------------

rule guppy_input:
    input:
        reads = config["data_dir"] + "/pass/{SAMPLES}"
    output:
        reads_merged = "results/merged_fastq/{SAMPLES}.merged.fastq"
    shell:
        "cat {input.reads}/*runid*.fastq > {output}"

# ------------------------------------------------------------------------------
# FASTQ INPUT:
# ------------------------------------------------------------------------------

rule fastq_input:
    input:
        reads = config["data_dir"] + "/{SAMPLES}.fastq"
    output:
        reads_merged = "results/merged_fastq/{SAMPLES}.merged.fastq"
    shell:
        "cp {input.reads} {output}"

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
# KRAKEN2: taxonomic classification
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
        re_html = "results/recentrifuge/" + config['my_experiment'] + ".html",
        re_xlsx = "results/recentrifuge/" + config['my_experiment'] + ".xlsx"
    shell:
        "rcf -a \
        -n {input.taxdb} \
        -k {params} \
        -o {output.re_html}"
