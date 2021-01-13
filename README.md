# Root-knot nematode metagenomic analysis
### Proof of concept for metagenomic analysis of root-knot nematode sequencing data from Oxford Nanopore Technologies' MinION Flongle platform.

### Quickstart

1. Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) (miniconda)

2. git clone this repository.  
run `git clone https://github.com/Graham-Sellers/RKN_metagenomic_analysis`
    
3. Install snakemake in your base conda environment.  
run `conda activate base`  
run `conda install -c bioconda -c conda-forge snakemake`
    
4. Get the data.  
Download test dataset and databases available at: insert dropbox link.  
Unzip it, there are 3 directories:  
    *RKN_lib3*  
    *meloidogyne_tomato_human_sweetpotato_no-mask_db*  
    *taxdump*
Move *RKN_lib3* to *data/libraries*.  
Move *meloidogyne_tomato_human_sweetpotato_no-mask_db* and *taxdump* to *data/databases*

5. Make it go.  
run `snakemake --use-conda --cores`
