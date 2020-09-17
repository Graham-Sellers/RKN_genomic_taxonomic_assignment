import os
from Bio import SeqIO
import random

def chop(infile, outdir):
    coverage = snakemake.config["coverage"]
    frag_lens = snakemake.config["size"]
    resample = snakemake.config["resample"]
    sample=infile.split("/")[-1].split(".fasta")[0]

    for frag_len in frag_lens:
        print(frag_len)
        with open(infile,"r") as fasta:
            handle = SeqIO.parse(fasta,"fasta")
            seqs=[s for s in handle if len(s) > frag_len]
            total_seqs=(len(seqs))
            for r in range(resample):
                sample_frag = (sample + "_" + str(frag_len))
                with open(outdir + "/" + sample_frag + "/" + sample_frag + "_" + str(r) + ".fasta", "w") as outfile:
                    for i in range(0,coverage):
                        n=random.randrange(0,total_seqs)
                        end=(len(seqs[n])+1)-frag_len
                        start=random.randrange(0,end)
                        stop=start+frag_len
                        wee=seqs[n].seq[start:stop]
                        outfile.write(">sequence_%s\n%s\n" %(i,wee))

for sample in snakemake.input:
    chop(sample, snakemake.params[0])
