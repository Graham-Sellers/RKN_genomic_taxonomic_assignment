################################################################################
# GENBANK RECORD GRABBING
# part 1 of database creation

from Bio import Entrez
import urllib.request
import datetime
import os
#
path = '/media/graham/Storage/RKN_snakemake/'
outdir = 'genomes/genbank_format/'
#
# with open(path+'search_terms.txt','r') as sp_list:
#     search_terms=sorted(list(set([line.strip() for line in sp_list])))

search_terms = ['Meloidogyne', 'GCA_000001405.28', 'GCA_012431665.1', 'GCA_002525835.2']
search_list = []
# , 'GCA_000001405.28', 'GCA_012431665.1', 'GCA_002525835.2'

Entrez.email='grahamssellers@gmail.com'
for term in search_terms:
    with Entrez.esearch(db='assembly', term=term, idtype='acc',retmax=100000) as handle:
        record=Entrez.read(handle)
        print('%s:\t\t%s' %(term,record['Count']))
        for r in record['IdList']:
            search_list.append(r)

# print(search_list)

ftp = []
links = []

for rec in search_list:
    with Entrez.esummary(db = 'Assembly', id = rec, report = 'full') as handle:
        seq = Entrez.read(handle, validate=False)
        doc_sum = seq['DocumentSummarySet']['DocumentSummary'][0]
        url = doc_sum['FtpPath_GenBank']
        if url == '':
            continue
        label = os.path.basename(url)
        link = os.path.join(url, label + '_genomic.gbff.gz')
        sp = doc_sum['SpeciesName']
        acc = doc_sum['LastMajorReleaseAccession']
        taxid = doc_sum['SpeciesTaxid']
        ftp.append('%s\t%s\t%s' %(sp, acc, taxid))
        print('%s\t%s' %(len(ftp), sp))
        links.append(link)

for link in links:
        label = os.path.basename(link).replace('.', '_', 1)
        print(label)
        print('start:\t' + datetime.datetime.now().strftime('%X'))
        urllib.request.urlretrieve(link, path + outdir + label)
        print('finish:\t' + datetime.datetime.now().strftime('%X'))

################################################################################
# GENBANK TO FASTA FORMAT FOR KRAKEN 2
# part 2 of database creation

from Bio import SeqIO
import gzip
import os
import datetime

#seq_list=[]
path = '/media/graham/Storage/RKN_snakemake/'
gb_in = 'genomes/genbank_format/'
fa_out = 'genomes/fasta_format/'

files = sorted(os.listdir(path + gb_in))

for file in files:
    file_name = file.split('_genomic.gbff.gz')[0]
    print(file_name)
    print('start:\t' + datetime.datetime.now().strftime('%X'))
    with gzip.open(path + gb_in + file, 'rt') as gb, open(path + fa_out + file_name + '.fasta', 'w') as outfile:
        for record in SeqIO.parse(gb, 'genbank'):
            acc = record.annotations['accessions'][0]
            species_name = record.annotations['organism']
            source = [f for f in record.features if f.type == 'source'][0]
            taxid = source.qualifiers['db_xref'][0].split(':')[1]
            kraken_taxid = ('>' + acc + '|kraken:taxid|' + taxid + ' ' + species_name)
            outfile.write('%s\n%s\n' %(kraken_taxid, record.seq))
    print('finish:\t' + datetime.datetime.now().strftime('%X'))

################################################################################
