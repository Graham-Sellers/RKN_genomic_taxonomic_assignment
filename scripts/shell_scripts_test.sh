# tests in shell
cd ../nems/

cd data/databases/
retaxdump

pwd

ls

kraken2-build --download-taxonomy --db data/databases/master_db

files=$(ls genomes/fasta_format/)
for file in $files
do
kraken2-build --add-to-library genomes/fasta_format/$file --db data/databases/master_db --no-masking
done

kraken2-build --build --db data/databases/master_db

# path is only here to allow for ease of wrangling from secondary HDD
path=/media/graham/Storage/nems/

################################################################################
# FASTQ TO FASTA
# converting fastq to fasta format (just incase this is easier than using straight fastq)

path=/media/graham/Storage/nems/

samples=$(ls $path/nanopore_test_data/01_nanopore_fastq_HAC_barcode_qfilterpass/ | grep '.fastq' | sed 's/.fastq//')

echo -e $samples

for sample in $samples
do
seqkit fq2fa $path/nanopore_test_data/01_nanopore_fastq_HAC_barcode_qfilterpass/${sample}.fastq \
-o $path/nanopore_test_data/02_nanopore_fasta_HAC_barcode_qfilterpass/${sample}.fasta
done

################################################################################
# KRAKEN2 FASTA INPUT
# using fasta input for Kraken2 with a small confidence value (0.01)

path=/media/graham/Storage/nems

mkdir $path/nanopore_test_data/output/HAC_fasta_input/kraken2_standard

files=$(ls $path/nanopore_test_data/02_nanopore_fasta_HAC_barcode_qfilterpass/ | grep '.fastq' | sed 's/.fasta//')

for a in $files
do
kraken2 --db $path/meloidogyne_dbs/meloidogyne_db \
$path/nanopore_test_data/02_nanopore_fasta_HAC_barcode_qfilterpass/${a}.fasta \
--use-names \
--memory-mapping \
--threads 10 \
--confidence 0.01 \
--output $path/nanopore_test_data/output/HAC_fasta_input/kraken2_standard/${a}.tsv \
--report $path/nanopore_test_data/output/HAC_fasta_input/kraken2_standard/${a}.txt
done

################################################################################
# KRAKEN2 FASTQ INPUT
# using fastq input for Kraken2 with a small confidence value (0.01)
# additionally uses --minimum-base-quality 10 (another level of quality control)
# bases/reads below this quality are ignored for taxonomic assignment

path=/media/graham/Storage/nems

mkdir $path/nanopore_test_data/output/HAC_fastq_input/kraken2_standard

files=$(ls $path/nanopore_test_data/01_nanopore_fastq_HAC_barcode_qfilterpass/ | grep '.fastq' | sed 's/.fastq//')

for a in $files
do

db=data/databases/meloidogyne_tomato_human_sweetpotato_no-mask_db

kraken2-inspect --db $db > inspect.txt --threads 10

db=data/databases/meloidogyne_tomato_human_sweetpotato_no-mask_db
kraken2 --db $db \
data/libraries/GCA_004785735.1_ASM478573v1_genomic.fasta \
--use-names \
--memory-mapping \
--threads 10 \
--confidence 0 \
--output report_odd1.tsv \
--report report.odd1.txt

done
ls ~/github/RKN_metagenomic_analysis/data/databases/meloidogyne_tomato_human_sweetpotato_no-mask_db
################################################################################
# KRAKEN2 QC FASTQ INPUT
# using fastq input from fastp for Kraken2 with a small confidence value (0.01)
# additionally uses --minimum-base-quality 10 (another level of quality control)
# bases/reads below this quality are ignored for taxonomic assignment

path=/media/graham/Storage/nems

mkdir $path/nanopore_test_data/output/HAC_fastq_input/kraken2_standard_qc

files=$(ls $path/nanopore_test_data/03_qc/ | grep '.fastq' | sed 's/.fastq//')

for a in $files
do
kraken2 --db $path/meloidogyne_dbs/meloidogyne_db \
$path/nanopore_test_data/03_qc/${a}.fastq \
--minimum-base-quality 10 \
--use-names \
--memory-mapping \
--threads 10 \
--confidence 0.01 \
--output $path/nanopore_test_data/output/HAC_fastq_input/kraken2_standard_qc/${a}.tsv \
--report $path/nanopore_test_data/output/HAC_fastq_input/kraken2_standard_qc/${a}.txt
done


seqkit stats -a
