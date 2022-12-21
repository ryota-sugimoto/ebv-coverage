#!/usr/bin/env bash

[ $# == 1 ] || { echo "${0} id" &2>1; exit 1; }

id=${1}
id_url_md5=/home/r-sugimoto/ebv_gwas/ref/id_url_md5
url=$(awk -v id=${id} '$1==id{print $2}' ${id_url_md5})
md5=$(awk -v id=${id} '$1==id{print $3}' ${id_url_md5})

human_ref=/home/r-sugimoto/ebv_gwas/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
human_stats_targets=/home/r-sugimoto/ebv_gwas/ref/bed/GRCh38.target_regions
ebv_ref=/home/r-sugimoto/ebv_gwas/ref/NC_007605.fasta


#download
wget "${url}" &> ${id}.wget.log
cram_file=$(basename ${url})

#check md5
downloaded_md5=$(md5sum ${cram_file} | awk '{print $1}' | tee ${id}.md5)
[[ ${md5} == ${downloaded_md5} ]] \
  || { echo '${id} md5 not much' 2>&1; exit 1; }

#collate
samtools collate --threads 10 --reference ${human_ref} \
  ${cram_file} \
  ${id}.collated &> ${id}.collate.log

#extract fastq
samtools fastq --threads 10 --reference ${human_ref} \
  -1 ${id}.1.fastq.gz \
  -2 ${id}.2.fastq.gz \
  -0 ${id}.unpaired.fastq.gz \
  ${id}.collated.bam &> ${id}.fastq_extraction.log

#map to EBV
bwa mem -t 10 ${ebv_ref} \
  ${id}.1.fastq.gz \
  ${id}.2.fastq.gz 2> ${id}.mem.log \
  | awk '/^@/ || !and($2,0x4)' \
  | samtools sort 2> ${id}.sort.log \
  > ${id}.ebv.bam

#collate
samtools collate --threads 10 ${id}.ebv.bam ${id}.ebv.collated

#fix mate
samtools fixmate -m ${id}.ebv.collated.bam ${id}.ebv.fixed_mates.bam

#mark duplicate
samtools sort ${id}.ebv.fixed_mates.bam > ${id}.ebv.sorted.bam
samtools markdup ${id}.ebv.sorted.bam ${id}.ebv.final.bam

#stats
samtools stats -f 2 -F 1796 -c 1,1000000,1 ${id}.ebv.final.bam \
  > ${id}.ebv.stats
samtools stats -f 2 -F 1796 -c 1,1000000,1 -t ${human_stats_targets} \
  ${cram_file} > ${id}.human.stats

#clean up
rm ${cram_file}
rm ${id}.collated.bam
rm ${id}.1.fastq.gz
rm ${id}.2.fastq.gz
rm ${id}.unpaired.fastq.gz
rm ${id}.wget.log

rm ${id}.ebv.bam
rm ${id}.ebv.collated.bam
rm ${id}.ebv.fixed_mates.bam
rm ${id}.ebv.sorted.bam
