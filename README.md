# Root-knot nematode metagenomic analysis
### Proof of concept for metagenomic analysis of root-knot nematode sequencing data from Oxford Nanopore Technologies' MinION Flongle platform.

*Kraken 2* taxonomic assignment of guppy basecaller HAC basecalled MinION Flongle sequencing data. This is the final stage of a larger workflow designed for accurate taxonomic identification of individual root-knot nematodes via extraction and sequencing of long read genomic DNA.

---

## The Workflow
*Note: This is a first draft. There will be more to come*

**Inputs:**  
- Basecalled library output from ONT Guppy basecaller (includes pass and fail directories plus sequencing metadata files)  
- Custom Kraken 2 database
- NCBI taxonomy nodes
- Sample sheet (a .tsv file listing barcodes to be analysed)  

### **Steps**  

**1. Merge fastq files**  
Merge all the fastq files inside each barcode directory of the nanopore sequencing data.  

**2. *Nanofilt* QC and trimming**  
Quality control of sequencing data (mean quality of 12), trimming of leading and trailing ends.  

**3. *Kraken 2* taxonomic assignment**  
A small custom database of specific genomes is used for identification of sequences.  

**4. *Recentrifuge* figure generation**  
Create human-readable and interactive html files from the *Kraken 2* outputs

---

### Quickstart

1. Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) (miniconda)

2. git clone this repository.  
run `git clone https://github.com/Graham-Sellers/RKN_genomic_taxonomic_assignment`

3. Install snakemake in your base conda environment.  
run `conda activate base`  
run `conda install -c bioconda -c conda-forge snakemake`

4. Get the data.  
Download test dataset and databases available from [Dropbox](https://www.dropbox.com/sh/5izuwb2ks61xbqg/AACzjETDpjWZh-d8R_qxYzWxa?dl=0)   
Unzip it, there are 3 directories:  
    *RKN_lib3* (MinION basecalled + demultiplexed library)  
    *meloidogyne_tomato_human_sweetpotato_no-mask_db* (Kraken2 database ~ 8Gb)
    *taxdump* (taxonomy nodes for recentrifuge)  
    Move *RKN_lib3* to *data/libraries*.  
Move *meloidogyne_tomato_human_sweetpotato_no-mask_db* and *taxdump* to *data/databases*

5. Make it go.  
run `snakemake --use-conda --cores`
