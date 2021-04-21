#!/bin/bash

#master script for rna-seq analysis

module load singularity 

sample_sheet=($(cat samples.txt))

mkdir -p sequences

for x in ${sample_sheet[@]};do
singularity exec -eC -B "$(pwd)" -H "$(pwd)" rna.sif \
fasterq-dump --split-files -O ./sequences -t ./sequences $x
done

#Quality Control

fastqs=($(find ./sequences -wholename *fastq))
mkdir -p results/quality

singularity exec -eC -B "$(pwd)" -H "$(pwd)" rna.sif \
fastqc -o results/quality ${fastqs[@]}

#Download genome assembly

mkdir genome
wget ftp://ftp.ensembl.org/pub/release-99/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz
mv Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz genome 
gunzip ./genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz

#Download genome annotation
wget ftp://ftp.ensembl.org/pub/release-99/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.99.gtf.gz
mv Saccharomyces_cerevisiae.R64-1-1.99.gtf.gz genome
gunzip ./genome/Saccharomyces_cerevisiae.R64-1-1.99.gtf.gz

#Compress genome assembly & annotation file
mkdir -p ./genome/index

singularity exec -eC -B "$(pwd)" -H "$(pwd)" rna.sif \
STAR --runThreadN 8 --runMode genomeGenerate \
--genomeDir ./genome/index \
--genomeFastaFiles ./genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa \
--sjdbGTFfile ./genome/Saccharomyces_cerevisiae.R64-1-1.99.gtf

#Alignment results

mkdir alignment_results 

S=($(find ./sequences -wholename "*.fastq" -printf '%f\n'|sed "s/_[12].fastq//"|sort|uniq))

for SRR in ${S[@]}; do 
	
	R1="./sequences/${SRR}_1.fastq"
	R2="./sequences/${SRR}_2.fastq"

	singularity exec -eC -B "$(pwd)" -H "$(pwd)" rna.sif \
	STAR --runThreadN 8 --genomeDir ./genome/index \
	--readFilesIn ${R1} ${R2} \
	--outFileNamePrefix ./alignment_results/${SRR}_ \
	--outSAMtype BAM SortedByCoordinate

	#Index BAM file
        singularity exec -eC -B "$(pwd)" -H "$(pwd)" rna.sif \
	samtools index ./alignment_results/${SRR}_Aligned.sortedByCoord.out.bam
done

#Feature Counting 
bamfiles=($(find ./alignment_results -wholename *bam))
mkdir -p counts

singularity exec -eC -B "$(pwd)" -H "$(pwd)" rna.sif \
featureCounts \
-a "./genome/Saccharomyces_cerevisiae.R64-1-1.99.gtf" \
-o ./counts/feature_counts \
-F GTF \
-t exon \
-g gene_id \
--minOverlap 10 \
--largestOverlap \
-s 2 \
-p \
-T 8 \
${bamfiles[@]}

