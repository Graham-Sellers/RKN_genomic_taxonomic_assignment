import os

configfile: "config.yaml"

samples, = glob_wildcards("input/{sample}.fasta")

rule all:
    input:
        reports = expand("results/kraken/{sample}_{size}/{sample}_{size}_{resample}.txt", sample = samples, size = config["size"],resample = range(config["resample"])),
        outputs = expand("results/kraken/{sample}_{size}/{sample}_{size}_{resample}.tsv", sample = samples, size = config["size"],resample = range(config["resample"]))

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

# directories, samples, = glob_wildcards("results/chopped/{dir}/{sample}.fasta")

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
    #     "kraken2 --db {params.db} {input[0]} --use-names --memory-mapping --threads 10 --confidence 0.0 --report {output[0]} --output {output[1]}"
    run:
        for read in params.file:
            infile=("results/chopped/" + read + ".fasta")
            report=("results/kraken/" + read + ".txt")
            output=("results/kraken/" + read + ".tsv")
            shell("kraken2 --db {params.db} {infile} --use-names --memory-mapping --threads 10 --confidence 0.0 --report {report} --output {output}")


        # print(input.reads[1])
        # print(output.reports[1])
        # print(output.outputs[1])

        # print(len(input.reads))
        # for read in input.reads:
        #     print(read)

#     run:
#         for sample in samples:
#             for i in config["size"]:
#                 print('%s_%s_fasta' %(sample,i))
#
# rule test:
#     run:
#         for i in config["size"]:
#             print(i)
#         for i in config["coverage"]:
#             print(i)
#         for i in config["error"]:
#             print(i)
#
#
# rule all:
#     input:
#         expand("genomes/meloidogyne_genomes/fasta_format/{sample}.fasta", sample = samples)
#
# rule gb2fasta:
#     input:
#         expand("genomes/meloidogyne_genomes/genbank_format/{sample}.gb.gz", sample = samples)
#     output:
#         expand("genomes/meloidogyne_genomes/fasta_format/{sample}.fasta", sample = samples)
#     params:
#         outdir = "genomes/meloidogyne_genomes/fasta_format"
#     threads:
#         10
#     script:
#         "scripts/gb2fasta.py"
#
# rule krakendb:
#     # input:
#     #     expand("genomes/meloidogyne_genomes/fasta_format/{sample}.fasta", sample = samples)
#     output:
#         directory("meloidogyne_db/taxonomy")
#     threads:
#         10
#     params:
#          db = "meloidogyne_db"
#     shell: """
#         kraken2-build --download-taxonomy --db {params.db} --threads {threads}
#         """
#     # run:
#     #     import gzip
    #     from Bio import SeqIO
    #
    #     def gb2fasta(infile, outdir):
    #         print('converting %s.gb.gz' %(infile))
    #         sample=infile.split('/')[-1].split('.gb.gz')[0]
    #         with gzip.open(infile, 'rt') as gb, open(outdir + '/' + sample + '.fasta', 'w') as fasta:
    #             for record in SeqIO.parse(gb, 'genbank'):
    #                 source = [f for f in record.features if f.type == 'source'][0]
    #                 species_name = (' '.join(source.qualifiers['organism'][0].split()[0:2]))
    #                 taxid = source.qualifiers['db_xref'][0].split(':')[1]
    #                 kraken_taxid = ('>' + record.id + '|kraken:taxid|' + taxid + ' ' + species_name)
    #                 fasta.write('%s\n%s\n' %(kraken_taxid, record.seq))
    #
    #     for sample in input:
    #         gb2fasta(sample, params[0])

# rule gb2fasta:
#     output:
#         expand("genomes/meloidogyne_genomes/fasta_format/{sample}.fasta", sample=samples)
#     params:
#         indir="genomes/meloidogyne_genomes/genbank_format",
#         outdir="genomes/meloidogyne_genomes/fasta_format"
#     run:
#         for sample in samples:
#             print('converting %s/%s.gb.gz' %(params.indir,sample))
#             with gzip.open(params.indir+'/'+sample+'.gb.gz','rt') as gb:
#                 with open(params.outdir+'/'+sample+'.fasta','w') as out:
#                     for record in SeqIO.parse(gb,'genbank'):
#                         source = [f for f in record.features if f.type == 'source'][0]
#                         species_name=(' '.join(source.qualifiers['organism'][0].split()[0:2]))
#                         taxid=source.qualifiers['db_xref'][0].split(':')[1]
#                         kraken_taxid=('>'+record.id+'|kraken:taxid|'+taxid+' '+species_name)
#                         out.write('%s\n%s\n' %(kraken_taxid,record.seq))



#
# rule kraken2_taxonomy:
#     input:
#         reads = lambda wildcards: sample_reads[wildcards.samp],
#     output:
#         krak_report = join(outdir, "classification/{samp}.krak.report")
#     params:
#         db = config['database'],
#         paired_string = paired_string
#     threads: kraken_threads
#     resources:
#         mem=kraken_memory,
#         time=6
#     singularity: "shub://bsiranosian/bens_1337_workflows:kraken2"
#     shell: """
#         time kraken2 --db {params.db} --threads {threads} --output {output.krak_report} \
#         --report {output.krak_report} {params.paired_string} {input.reads} --use-names
#         """
