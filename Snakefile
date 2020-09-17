import os

configfile: "config.yaml"

samples, = glob_wildcards("input/{sample}.fasta")

rule all:
    input:
        reports = expand("results/kraken/{sample}_{size}/{sample}_{size}_{resample}.txt", sample = samples, size = config["size"],resample = range(config["resample"])),
        outputs = expand("results/kraken/{sample}_{size}/{sample}_{size}_{resample}.tsv", sample = samples, size = config["size"],resample = range(config["resample"]))

# CHOP INPUT GENOME INTO FRAGMENTS OF A GIVEN LENGTH

rule chop:
    input:
        expand("input/{sample}.fasta", sample = samples)
    output:
        expand("results/chopped/{sample}_{size}/{sample}_{size}_{resample}.fasta", sample = samples, size = config["size"],resample = range(config["resample"]))
    params:
        outdir = "results/chopped"
    threads:
        10
    script:
        "scripts/chop.py"

# KRAKEN2 FOR TAXONOMIC ASSIGNMENT PER INDIVIDUAL FASTA OUTPUT FROM CHOP

rule kraken:
    input:
        reads = expand("results/chopped/{sample}_{size}/{sample}_{size}_{resample}.fasta", sample = samples, size = config["size"],resample = range(config["resample"]))
    output:
        reports = expand("results/kraken/{sample}_{size}/{sample}_{size}_{resample}.txt", sample = samples, size = config["size"],resample = range(config["resample"])),
        outputs = expand("results/kraken/{sample}_{size}/{sample}_{size}_{resample}.tsv", sample = samples, size = config["size"],resample = range(config["resample"]))
    params:
        db = directory("database"),
        file = expand("{sample}_{size}/{sample}_{size}_{resample}", sample = samples, size = config["size"],resample = range(config["resample"]))
    threads:
        10
    # shell:
    #     "kraken2 --db {params.db} {input.reads} --use-names --memory-mapping --threads 10 --confidence 0.0 --report {output.reports} --output {output.outputs}"
    run:
        for read in params.file:
            infile=("results/chopped/" + read + ".fasta")
            report=("results/kraken/" + read + ".txt")
            output=("results/kraken/" + read + ".tsv")
            shell("kraken2 --db {params.db} {infile} --use-names --memory-mapping --threads 10 --confidence 0.0 --report {report} --output {output}")
