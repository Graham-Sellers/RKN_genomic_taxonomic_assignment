# Genomic taxonomic assignment of individual root-knot nematodes
### A genomics based species diagnostic tool for root-knot nematodes (RKN; genus *Meloidogyne*) using sequencing data from Oxford Nanopore Technologies' Flongle platform.

*Kraken 2* taxonomic assignment of ONT Guppy basecaller high accuracy (HAC) basecalled Flongle sequencing data. This is the final stage of a larger workflow designed for accurate taxonomic identification of individual RKN via extraction and sequencing of long read genomic DNA.

---

**Input sequencing data preparation**

This workflow is designed to analyse a single Flongle sequencing library prepared with ONT Rapid PCR Barcoding Kit (SQK-RPB004) and basecalled with Guppy GPU high accuracy (HAC) basecaller using the following:  

`guppy_basecaller --input_path path/to/fast5_directory -r --save_path path/to/output_directory --config dna_r9.4.1_450bps_hac.cfg --device cuda:0 --min_qscore 7 --qscore_filtering --barcode_kits "SQK-RPB004" --trim_barcodes --require_barcodes_both_ends`  

This needs to be run on a HPC GPU node or a local machine with GPU computing capability see
[here](https://community.nanoporetech.com/requirements_documents/minion-it-reqs.pdf).  

The result is a basecalled library output directory containing *pass/* and *fail/* directories plus sequencing metadata files.  
This is the basecalled library output directory used in the the workflow.

---

## The Workflow

**Inputs**  
- Basecalled library output directory from ONT Guppy basecaller  
- Kraken 2 database
- NCBI taxonomy nodes (names.dmp, nodes.dmp)
- Sample sheet (a .tsv file listing barcodes/sample names to be analysed, one per line)  

### **Steps**  

**1. Merge fastq files**  
Merge all the fastq files inside each barcode directory of the basecalled library output.  

**2. *Nanofilt* QC and trimming**  
Quality control of sequencing data and read trimming.  

**3. *Kraken 2* taxonomic assignment**  
A small *Meloidogyne* specific database is used for classification of sequences.  

**4. *Recentrifuge* figure generation**  
Create human-readable and interactive html files from the *Kraken 2* outputs

---

### Quick start

1. Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) (miniconda)

2. git clone this repository.  
run `git clone https://github.com/Graham-Sellers/RKN_genomic_taxonomic_assignment`

3. Install snakemake in your base conda environment.  
run `conda activate base`  
run `conda install -c bioconda -c conda-forge snakemake`

4. Get the test data.  
Download *RKN_test_dataset.zip* and *RKN_databases_2021-03.zip* available from [OSF](http://dx.doi.org/10.17605/OSF.IO/VA7S2).  
Unzip them, this will give you three directories and a .tsv file:  
*RKN_db/* (Kraken2 database ~ 6.5 Gb directory containing 3 .k2d files)  
*taxdump/* (nodes.dmp and names.dmp for recentrifuge)  
*RKN_test_lib.tsv* (list of sample names to be processed)  
*RKN_test_lib/* (a Guppy basecalled + demultiplexed library output directory)  
  
5. Put things where they need to be.  
Move *RKN_test_library/* to *data/libraries/*  
Move *RKN_db/* and *taxdump/* to *data/databases/*  
Move RKN_test_samples.tsv to the project main directory (*RKN_genomic_taxonomic_assignment/*)  
*(note: no need to alter config.yaml, it is set for the test dataset)*

6. Make it go.  
run `snakemake --use-conda --cores`

---

### Options for input data

The workflow takes in a Guppy basecaller output directory containing *pass/* and *fail/* directories plus sequencing metadata files (see above). It is recommended to rename the directories in the *pass/* directory to relevant and unique sample names - at present Guppy names them "barcode01", "barcode02" etc. These sample names (if changed) should be reflected in the sample sheet (.tsv file targeted by *"sample_list"* in *config.yaml*). This does however limit the analysis to processing 12 samples at a time.  

**Input data type**

The workflow's *config.yaml* has an input data type option, creatively called *"input_type"*. This is set as *"guppy"* to process a basecalled library output directory from Guppy GPU basecaller. However, by changing the option to *"fastq"* (in *config.yaml*) the workflow then assumes a directory (determined by *"data_dir"* in *config.yaml*) containing fastq files, one per sample (unzipped). This is for those who want to, for example, run the data from the accompanying NCBI BioProject (accession number PRJNA706653) SRA fastq data, without having to basecall raw fast5 MinKNOW output, or multiple libraries (by merging fastqs in each barcoded directory from Guppy basecaller output).
