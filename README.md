# Root-knot nematode metagenomic analysis
### Proof of concept for metagenomic analysis of root-knot nematode sequencing data from Oxford Nanopore Technologies' MinION Flongle platform.

## Quickstart

1. install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) (miniconda)

2. git clone this repository
    - `git clone https://github.com/davelunt/Tapirs`
    
3. install snakemake in your base conda environment
    - `conda activate base`
    - `conda install -c bioconda -c conda-forge snakemake`
    
4. test dataset and databases are available at: insert drop boxlink. Move RKN_lib3 in data/libraries. Move meloidogyne_tomato_human_sweetpotato_no-mask_db and taxdump into0 data/databases

5. run `snakemake --use-conda --cores`
