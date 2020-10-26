
samples, = glob_wildcards("01_raw_data/{sample}.fasta")
#this will have to change, needs to be something from the merged barcode files

################################################################################
# MERGE ALL FASTQ FILES IN BARCODE DIRECTORIES

rule merge_files:
    input:
        reads = expand("01_raw_data/{sample}_{size}/{sample}_{size}_{resample}.fasta", sample = samples, size = config["size"],resample = range(config["resample"]))
    output:
        reports = expand("results/kraken/{sample}_{size}/{sample}_{size}_{resample}.txt", sample = samples, size = config["size"],resample = range(config["resample"]))

################################################################################
# QC all reads with fastp, use this to size select too?
# clean out and tidy the data prior to analysis.

rule pastp_qc:
    input:
        reads = xxxxx # wildcard
    output:
        cleaned = yyyyy # wildcard

################################################################################
# KRAKEN2 TAXONOMIC ASSIGNMENT
# will use fastq input
# selected database will influence the outcome

rule kraken2:
    input:
        reads = expand("results/chopped/{sample}_{size}/{sample}_{size}_{resample}.fastq", sample = samples, size = config["size"],resample = range(config["resample"]))
    output:
        reports = expand("results/kraken/{sample}_{size}/{sample}_{size}_{resample}.txt", sample = samples, size = config["size"],resample = range(config["resample"])),
        outputs = expand("results/kraken/{sample}_{size}/{sample}_{size}_{resample}.tsv", sample = samples, size = config["size"],resample = range(config["resample"]))
    params:
        db = directory("database"),
        file = expand("{sample}_{size}/{sample}_{size}_{resample}", sample = samples, size = config["size"],resample = range(config["resample"]))
    threads:
        10
    # shell: # this doesnt work
    #     "kraken2 --db {params.db} {input.reads} --use-names --memory-mapping --threads 10 --confidence 0.0 --report {output.reports} --output {output.outputs}"
    run:
        for read in params.file:
            infile=("results/chopped/" + read + ".fastq")
            report=("results/kraken/" + read + ".txt")
            output=("results/kraken/" + read + ".tsv")
            shell("kraken2 --db {params.db} {infile} --minimum-base-quality 10 --use-names --memory-mapping --threads 10 --confidence 0.01 --report {report} --output {output}")

################################################################################
# KRAKEN2 TAXONOMIC OUTPUT
# wrangle all files into pandas dataframe and reorder to give a simple by sample reads per tax level .tsv file

# needs to be wiggled but look at python-shell_test_scripts/python_scripts_test (making a Kraken2 tsv output)
