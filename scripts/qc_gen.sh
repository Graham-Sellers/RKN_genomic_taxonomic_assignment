
# length and q 0, 7, 15

mkdir qc_gen/fastp
mkdir qc_gen/fastp/fastp_out
mkdir qc_gen/fastp/failed_out
mkdir qc_gen/kraken2
mkdir qc_gen/kraken2/outputs
mkdir qc_gen/kraken2/reports

path=/media/graham/Storage/RKN_snakemake/results/merged_fastq


input=barcode05
qual=(0 15)
len=(0 1000 2000 4000)

for l in ${len[@]}
do
  for q in ${qual[@]}
  do
    fastp -i ${path}/${input}.merged.fastq \
    -o qc_gen/fastp/${input}_${l}_${q}_0_0_0.qc.fastq \
    -A \
    -q $q \
    -l $l \
    --failed_out qc_gen/failed_out/${input}_${l}_${q}_0_0_0.failed \
    -h qc_gen/fastp/fastp_out/${input}_${l}_${q}_0_0_0.html \
    -j qc_gen/fastp/fastp_out/${input}_${l}_${q}_0_0_0.json
  done
done

############

## trim + length @ q15

t=30
len=(0 1000 2000 4000)
q=15
for l in ${len[@]}
do
  fastp -i ${path}/${input}.merged.fastq \
  -o qc_gen/fastp/${input}_${l}_${q}_${t}_0_0.qc.fastq \
  -A \
  -f $t \
  -q $q \
  -l $l \
  --failed_out qc_gen/failed_out/${input}_${l}_${q}_${t}_0_0.failed \
  -h qc_gen/fastp_out/${input}_${l}_${q}_${t}_0_0.html \
  -j qc_gen/fastp_out/${input}_${l}_${q}_${t}_0_0.json
done

# kraken2 base quality 0, confidence 0

database=data/databases/meloidogyne_tomato_human_sweetpotato_no-mask_db
pwd
files=$(ls -1 qc_gen/fastp/ | grep '.qc.fastq' | sed 's/.qc.fastq//')
for file in $files
do
  kraken2 --db $database \
  qc_gen/fastp/${file}.qc.fastq \
  --use-names \
  --threads 10 \
  --confidence 0.0 \
  --report qc_gen/kraken2/reports/${file}.txt \
  --output qc_gen/kraken2/outputs/${file}.krk
done

## confidence

database=data/databases/meloidogyne_tomato_human_sweetpotato_no-mask_db
conf=(0.00 0.01)
files=$(ls -1 qc_gen/fastp/ | grep '_30_0_0.qc.fastq' | sed 's/_0_0.qc.fastq//')
for file in $files
do
  echo $file
  for c in ${conf[@]}
  do
    kraken2 --db $database \
    qc_gen/fastp/${file}_0_0.qc.fastq \
    --minimum-base-quality 10 \
    --use-names \
    --threads 10 \
    --confidence ${c} \
    --report qc_gen/kraken2/reports/${file}_10_${c}.txt \
    --output qc_gen/kraken2/outputs/${file}_10_${c}.krk
  done
done

#####################################
# fastp
database=data/databases/meloidogyne_tomato_human_sweetpotato_no-mask_db
files=$(ls -1 qc_gen/fastp/ | grep '_30_0_0.qc.fastq' | sed 's/_0_0.qc.fastq//')
for file in $files
do
  kraken2 --db $database \
  qc_gen/fastp/${file}_0_0.qc.fastq \
  --minimum-base-quality 7 \
  --use-names \
  --threads 10 \
  --confidence 0.02 \
  --report qc_gen/kraken2/reports/${file}_7_0.02.txt \
  --output qc_gen/kraken2/outputs/${file}_7_0.02.krk
done
