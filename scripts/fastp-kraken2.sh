#fastp
################################################################################
# FASTP QUALITY CONTROL
# disabled adapter detection, chop front 30 bases, minimum q score 15, minimum length 2000

path=/media/graham/Storage/RKN_lib3
mkdir $path/02_qc
mkdir $path/02_qc/fastp_out/
mkdir $path/02_qc/failed_out/

samples=$(ls $path/01_nanopore_fastq_HAC_barcode_qfilterpass/ | grep '.fastq' | sed 's/.fastq//')

for sample in $samples
do
    fastp -i $path/01_nanopore_fastq_HAC_barcode_qfilterpass/${sample}.fastq \
    -o $path/02_qc/${sample}.fastq \
    -A \
    -f 30 \
    -q 15 \
    -l 2000 \
    --failed_out $path/02_qc/failed_out/${sample}_failed.fastq \
    -h $path/02_qc/fastp_out/${sample}_fastp.html \
    -j $path/02_qc/fastp_out/${sample}_fastp.json
done

################################################################################
# KRAKEN2 QC FASTQ INPUT
# using fastq input from fastp for Kraken2 with a small confidence value (0.01)
# additionally uses --minimum-base-quality 10 (another level of quality control)
# bases/reads below this quality are ignored for taxonomic assignment

path=/media/graham/Storage/RKN_lib3
mkdir $path/Kraken2_output
mkdir $path/Kraken2_output/kraken2_RKN_bothbarcodes_bact
db=/media/graham/Storage/nems/meloidogyne_dbs/meloidogyne_tomato_human_db

files=$(ls $path/02_qc/ | grep '.fastq' | sed 's/.fastq//')

for a in $files
do
kraken2 --db $path/RKN_bact \
$path/02_qc/${a}.fastq \
--minimum-base-quality 10 \
--use-names \
--memory-mapping \
--threads 10 \
--confidence 0.01 \
--output $path/Kraken2_output/kraken2_RKN_bothbarcodes_bact/${a}.tsv \
--report $path/Kraken2_output/kraken2_RKN_bothbarcodes_bact/${a}.txt
done

RKN_bact
$path/meloidogyne_tomato_db \
$path/meloidogyne_tomato_human_db
/media/graham/Storage/nems/meloidogyne_dbs/meloidogyne_tomato_human_db

path=/media/graham/Storage/RKN_lib2

kraken2-inspect --db $path/meloidogyne_tomato_db > $path/inspect.txt

# --minimum-base-quality 12 \
# --minimum-base-quality 15 \

####################################
