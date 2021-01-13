# Root-knot nematode metagenomic analysis
### Proof of concept for metagenomic analysis of root-knot nematode sequencing data from Oxford Nanopore Technologies' MinION Flongle platform.

## Quickstart

1. Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) (miniconda)

2. git clone this repository  
   run `git clone https://github.com/Graham-Sellers/RKN_metagenomic_analysis`
    
3. Install snakemake in your base conda environment  
   run `conda activate base` (if it isn't already)  
   run `conda install -c bioconda -c conda-forge snakemake`
    
4. Test dataset and databases are available at: insert dropbox link.  
Move RKN_lib3 to data/libraries.  
Move meloidogyne_tomato_human_sweetpotato_no-mask_db and taxdump to data/databases

5. run `snakemake --use-conda --cores`
