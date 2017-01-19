#$ -S /bin/sh

# Alignment script for Grid Engine: gets working dir and fastq_id as parameter

bwa="/illumina/pipeline/bwa-0.7.12/bwa"
ref_b37="/illumina/runs/Genome_b37/human_g1k_v37_decoy.fasta"

cd $1

# Alignment using ensembl b37 decoy reference genome (gatk bundle b37)
# $NSLOTS holds the actual number of available CPU cores set by the grid engine
echo "Alignment, BWA mem" 
$bwa mem -M -t $NSLOTS -R "@RG\tID:$2\tSM:$2\tLB:$2\tPL:Illumina" $ref_b37 $2_1.fq.gz $2_2.fq.gz | gzip > $2.sam.gz 
