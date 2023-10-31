#!/bin/sh
#
#SBATCH --job-name=RNA-seq
#SBATCH --output=out.txt
#SBATCH --mem=8G
#SBATCH -c 2
#

i=$1
n=$2
ref='ASM972v1'  

mkdir qc
mkdir trimmed
mkdir index
mkdir mapped
mkdir counts
mkdir qc_post
mkdir mutations

  ## initial quality control

fastqc \
raw_data/S0${i}_S${i}_R1_001.fastq.gz \
raw_data/S0${i}_S${i}_R2_001.fastq.gz \
-o qc/
if [${i} == ${n}}
then
  multiqc qc/*_fastqc.zip -o qc/

  ## trim
  
trimmomatic PE -threads 8 \
raw_data/S0${i}_S${i}_R1_001.fastq.gz raw_data/S0${i}_S${i}_R2_001.fastq.gz \
trimmed/S0${i}_fw_p.fastq.gz trimmed/S0${i}_fw_u.fastq.gz \
trimmed/S0${i}_rv_p.fastq.gz trimmed/S0${i}_rv_u.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True  LEADING:33 TRAILING:32 MINLEN:36

	## post trimming qc
 
fastqc \
trimmed/trimmed/S0${i}_fw_p.fastq.gz \
trimmed/S0${i}_rv_p.fastq.gz \
-o qc_post/
if [${i} == ${n}}
then
  srun multiqc qc_post/*_fastqc.zip -o qc/

	## bowtie2
 
if [${i} == 1}
then
  bowtie2-build ref/${ref}.fa index/Syn
bowtie2 --threads 8 -x index/Syn \
-1 trimmed/S0${i}_fw_p.fastq.gz \
-2 trimmed/S0${i}_rv_p.fastq.gz \
-S mapped/S0${i}.sam
samtools sort -o mapped/S0${i}_sorted.sam mapped/S0${i}.sam

  ## quantify
  
htseq-count -s reverse -r pos -a 0 -t gene -i ID  mapped/S0${i}_sorted.sam ref/KO_test.gff > counts/S0${i}.txt

  ## inspect mutations
  
samtools view -b -S -o mapped/S0${i}.bam mapped/S0${i}.sam
samtools sort mapped/S0${i}.bam -o mapped/S0${i}_sorted.bam
samtools index mapped/S0${i}_sorted.bam
bcftools mpileup -f ref/${ref}.fa mapped/S0${i}_sorted.bam > mutations/S0${i}.bcf
bcftools call -v -c mutations/S0${i}.bcf > mutations/S0${i}.vcf
bcftools query -f '%CHROM %POS %DP %REF %ALT %AF1\n' mutations/S0${i}.vcf > mutations/S0${i}.txt

