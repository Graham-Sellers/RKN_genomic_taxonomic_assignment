import gzip
from Bio import SeqIO

def gb2fasta(infile, outdir):
    print('converting %s.gb.gz' %(infile))
    sample=infile.split('/')[-1].split('.gb.gz')[0]
    with gzip.open(infile, 'rt') as gb, open(outdir + '/' + sample + '.fasta', 'w') as fasta:
        for record in SeqIO.parse(gb, 'genbank'):
            source = [f for f in record.features if f.type == 'source'][0]
            species_name = (' '.join(source.qualifiers['organism'][0].split()[0:2]))
            taxid = source.qualifiers['db_xref'][0].split(':')[1]
            kraken_taxid = ('>' + record.id + '|kraken:taxid|' + taxid + ' ' + species_name)
            fasta.write('%s\n%s\n' %(kraken_taxid, record.seq))

for sample in snakemake.input:
    gb2fasta(sample, snakemake.params[0])
