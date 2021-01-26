# tests in python
################################################################################
# GB TO FASTA FORMAT for KRAKEN"

from Bio import SeqIO
import gzip
import os
import datetime

#seq_list=[]
path='/media/graham/Storage/nems/'
gb_in='genomes/bacteria_fungi/genbank_format/'
fa_out='genomes/bacteria_fungi/fasta_format/'

files=sorted(os.listdir(path+gb_in))

for file in files:
    file_name=file.split('.gb.')[0]
    print(file_name)
    print('start:\t'+datetime.datetime.now().strftime('%X'))
    with gzip.open(path+gb_in+file,'rt') as gb, open(path+fa_out+file_name+'.fasta','w') as outfile:
        for record in SeqIO.parse(gb,'genbank'):
            source = [f for f in record.features if f.type == 'source'][0]
            species_name=(' '.join(source.qualifiers['organism'][0].split()[0:2]))
            taxid=source.qualifiers['db_xref'][0].split(':')[1]
            kraken_taxid=('>'+record.id+'|kraken:taxid|'+taxid+' '+species_name)
            outfile.write('%s\n%s\n' %(kraken_taxid,record.seq))
    print('finish:\t'+datetime.datetime.now().strftime('%X'))

################################################################################

# path is only here to allow for ease of wrangling from secondary HDD
path='/media/graham/Storage/nems/'

################################################################################
#MERGE ALL FASTQ FILES IN EACH BARCODE DIRECTORY
# merge all fastq outputs per barcode file, rename and stick elsewhere

import shutil
import os
import glob

# path='/media/graham/Storage/nems/'
# inpath=path+'fast5_output/HAC_barcode_qfiltered/pass/'
# outpath=path+'nanopore_test_data/01_nanopore_fastq_HAC_barcode_qfilterpass/'

path='/media/graham/Storage/'
inpath=path+'nanopore_basecalled_data/RKN_lib3_HAC_barcode_bothends_qfiltered/pass/'
outpath=path+'RKN_lib3/01_nanopore_fastq_HAC_barcode_qfilterpass/'

barcodes=os.listdir(inpath)
print(barcodes)

for barcode in barcodes:
    outfilename = outpath + barcode + '.fastq'
    with open(outfilename, 'wb') as outfile:
        for filename in glob.glob(inpath + barcode + '/*.fastq'):
            print(filename)
            if filename == outfilename:
                continue
            with open(filename, 'rb') as readfile:
                shutil.copyfileobj(readfile, outfile)

################################################################################
# SIZE SELECTING READS
# for messing around with different minimum sequence length (pre Kraken2 analysis) see if accuracy increase?

from Bio import SeqIO
import glob
import os

path='/media/graham/Storage/nems/'

size=5000 # change for different minimum sizes

path_to_in=path+'nanopore_test_data/02_nanopore_fasta/'
path_to_out=path+'nanopore_test_data/03_size_selected_fasta/'
read_lengths=[]

files=os.listdir(path_to_in)
files.sort()
for file in files:
    with open(path_to_in+file,'r') as infile, open(path_to_out+file.replace('.fasta','_'+str(size)+'.fasta'),'w') as outfile:
        seqs=[s for s in SeqIO.parse(infile,'fasta') if len(s)>=size]
        #seqs=SeqIO.parse(infile,'fasta')
        #for seq in seqs if len(seq) >= size:
        SeqIO.write(seqs,outfile,'fasta')

################################################################################
# MAKING A KRAKEN2 TSV OUTPUT
# pandas dataframe of readcount per taxominic level for each sample
# human readable and easy to push through R etc?
# based on current rough directory system, will need to be changed in accordance

import pandas as pd
import numpy as np
import glob

path='/media/graham/Storage/nems/'

ma_list = glob.glob(path+'nanopore_test_data/output//HAC_fastq_input/kraken2_standard_qc/*.txt')
ma_list.sort()
whole_bits=[]
for output in ma_list:
    with open(output) as file:
        bits=[]
        bits.append(output)
        for line in file:
            line=line.split('\t')
            if int(line[2])>0:
                bits.append('%s\t%s' %(line[5].strip(),line[2]))
    whole_bits.append(bits)

# make into dataframe:

df=pd.DataFrame()
for i in range(len(whole_bits)):
    sample=whole_bits[i][0].split('/')[-1]
    dfp=pd.DataFrame([x.split('\t') for x in whole_bits[i][1:]],columns=['Taxonomy',sample])
    dfp.index=dfp['Taxonomy']
    dfp=dfp.drop(dfp.columns[0],axis=1)
    df=pd.concat([df,dfp],axis=1,sort=False).fillna(0)
df=df.transpose()

df.to_csv(path+'nanopore_test_data/output/HAC_kraken2_standard_fastq.tsv',sep='\t')
df

# change index names to be input_method_sample:
#
# dictionary={'barcode01':'mig_6',
#             'barcode02': 'mig_9',
#             'barcode03': 'mig_10',
#             'barcode04': 'mig_11',
#             'barcode05': 'mi♀️i_9',
#             'barcode06': 'mi♀️i_10'}
#
# index_rename=[]
# for index in df.index:
#     name=index.split('/')[2:]
#     for key in dictionary:
#         if key in name[2]:
#             name[2]=dictionary[key]
#     name='_'.join(name)
#     index_rename.append(name)
# index_rename
#
# df.index=index_rename
#
# # change order of columns (taxonomy)
#
# cols = list(df.columns.values)
# cols = ['unclassified',
#     'Meloidogyne',
#     'Meloidogyne incognita group',
#     'Meloidogyne incognita',
#     'Meloidogyne arenaria',
#     'Meloidogyne javanica',
#     'Meloidogyne floridensis',
#     'Meloidogyne enterolobii',
#     'Meloidogyne luci',
#     'Meloidogyne hapla',
#     'Meloidogyne graminicola']
# df=df[cols]


################################################################################
# READLENGTHS FOR HISTS IN R
# make pandas dataframe of all reads from each sample
# just for colour, not decided on the actual software for qc visualisation of reads yet.

from Bio import SeqIO
import os
import glob
import pandas as pd


path_to_in = 'results/qc/'
df = pd.DataFrame()

files = glob.glob(path_to_in + '*.fastq')
files.sort()
print(files)

# make into dataframe:

for file in files:
    print(file)
    with open(file, 'r') as infile:
        name = os.path.basename(file).split('.qc')[0]
        seqs = SeqIO.parse(infile, 'fastq')
        read_lengths = []
        for seq in seqs:
            read_lengths.append(len(seq))
        read_lengths.sort()
    dfp = pd.DataFrame(read_lengths, columns = [name])
    df = pd.concat([df, dfp], axis=1, sort=False)

df.to_csv('qc_readlengths.tsv', sep = '\t')

################################################################################
# GENOME ASSEMBLY GRABBING AND DOWNLOADING
# needs fixing

from Bio import Entrez
import urllib.request
import datetime

path='/media/graham/Storage/nems/'

Entrez.email='grahamssellers@gmail.com'

search_list=[]
#search=('Solanum lycopersicum')
term='Meloidogyne hapla'
with Entrez.esearch(db='assembly', term=term, idtype='acc',retmax=100000) as handle:
    record = Entrez.read(handle)
    print(record['Count'])
    for r in record['IdList']:
        search_list.append(r)


################################

GCA_000002995.2_ASM299v2


import os
from Bio import Entrez

term='nematoda'

labels = []
links = []

def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle,validate=False)
    return esummary_record

Entrez.email = "grahamssellers@gmail.com"
handle = Entrez.esearch(db="assembly", term=term, retmax='1000')
record = Entrez.read(handle)
ids = record['IdList']
print (f'found {len(ids)} ids')
links = []
for id in ids:
    #get summary
    summary = get_assembly_summary(id)
    #get ftp link
    url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
    if url == '':
        continue
    label = os.path.basename(url)
    # print(label)
    labels.append(label)
    link = os.path.join(url,label+'_genomic.gbff.gz')
    # print (link)
    links.append(link)


n=1
for shit in links:
    print('%s\t%s' %(n, shit))
    n+=1





'GCA_000147155.1_C_japonica-7.0.1'


if 'GCA_000172435.1_Freeze_1' in labels:
    print('oh yes')
else:
    print('shit, whay int it?')

dip=['i','1','l']

if 'i' in dip:
    print('yep')

ftp=[]
base='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/'

for rec in search_list:
    with Entrez.efetch(db='Assembly',id=rec,report='full',) as handle:
        seq=Entrez.read(handle,validate=False)
        doc_sum=seq['DocumentSummarySet']['DocumentSummary'][0]
        sp=doc_sum['SpeciesName']
        acc=doc_sum['LastMajorReleaseAccession']
        taxid=doc_sum['SpeciesTaxid']
        ftp.append('%s\t%s\t%s' %(taxid,sp,acc))
        print('%s\t%s' %(len(ftp),sp))
        print('start:\t'+datetime.datetime.now().strftime('%X'))
        base='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/'
        url=acc.split('.')[0].replace('_','')
        url=base+'/'.join(url[i:i+3] for i in range(0, len(url), 3))
        url_sub=str(urllib.request.urlopen(url).read()).split()[-1].split('\\')[0]
        url=url+'/'+url_sub+'/'+url_sub+'_genomic.gbff.gz'
        file_name=sp.replace(' ','_').replace('.','')+'.'+acc+'.gb.gz'
        print(url)
        print(doc_sum)
        # urllib.request.urlretrieve(url,path+'genomes/plant_genomes/genbank_format/'+file_name)
        print('finish:\t'+datetime.datetime.now().strftime('%X'))

################################################################################

from Bio import Entrez
import os

term='Solanum lycopersicum'

def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle,validate=False)
    return esummary_record

def get_assemblies(term, download=True, path='assemblies'):
    """Download genbank assemblies for a given search term.
    Args:
        term: search term, usually organism name
        download: whether to download the results
        path: folder to save to
    """

    from Bio import Entrez
    #provide your own mail here
    Entrez.email = "grahamssellers@gmail.com"
    handle = Entrez.esearch(db="assembly", term=term, retmax='200')
    record = Entrez.read(handle)
    ids = record['IdList']
    print (f'found {len(ids)} ids')
    links = []
    for id in ids:
        #get summary
        summary = get_assembly_summary(id)
        #get ftp link
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
        if url == '':
            continue
        label = os.path.basename(url)
        #get the fasta link - change this to get other formats
        link = os.path.join(url,label+'_genomic.gbff.gz')
        # print (link)
        # print(label)
        links.append(link)
#         if download == True:
#             #download link
        urllib.request.urlretrieve(link, f'{label}.fna.gz')
    return links

# ------------------------------------------------------

from Bio import Entrez
import os
import urllib.request

esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
esummary_record = Entrez.read(esummary_handle,validate=False)
term = 'nematoda'
Entrez.email = "grahamssellers@gmail.com"
handle = Entrez.esearch(db="assembly", term=term, retmax='200')
record = Entrez.read(handle)
ids = record['IdList']
print (len(ids))
links = []
for id in ids:
    #get summary
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle,validate=False)
    #get ftp link
    url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
    if url == '':
        continue
    label = os.path.basename(url)
    print(label)
    #get the fasta link - change this to get other formats
    link = os.path.join(url, label + '_genomic.fna.gz')
    print (link)
    links.append(link)
    # if download == True:
        #download link
        # urllib.request.urlretrieve(link, '{label}.fna.gz')
