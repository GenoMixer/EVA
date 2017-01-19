#$ -S /bin/sh

# variant calling script for Grid Engine
# best practices GATK workflow (GATK bundle b37)
# single sample variant calling pipeline using the haplotype caller in gVCF mode 
# paired-end sequencing on exon capture (agilent sure select, all exon v5) using HiSeq2500 Illumina machine

cd $1

# environment variables
   tools_dir="/illumina/pipeline/" 
   ref_dir="/illumina/runs/Genome_b37/"
   fastq_id=$2

# prerequisites tools  
   fastqc="/illumina/pipeline/FastQC/fastqc"
   bwa="/illumina/pipeline/bwa-0.7.12/bwa"
   samtools="/illumina/pipeline/samtools-1.3/samtools"	
   picard="/illumina/pipeline/picard-tools-2.0.1/picard.jar"
   gatk="/illumina/pipeline/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar"
   annovarinput="/illumina/pipeline/annovar/convert2annovar.pl"
   annovar="/illumina/pipeline/annovar/table_annovar.pl"
   annoDB="/illumina/pipeline/annovar/humandb"

# human ref
   ref_b37="/illumina/runs/Genome_b37/human_g1k_v37_decoy.fasta"      

# dbSNP, 1000G, HapMap, Mills 
   phase1_indels="/illumina/runs/Genome_b37/1000G_phase1.indels.b37.vcf"
   mills_indels="/illumina/runs/Genome_b37/Mills_and_1000G_gold_standard.indels.b37.vcf"
   dbsnp_138="/illumina/runs/Genome_b37/dbsnp_138.b37.vcf"
   phase1_snps="/illumina/runs/Genome_b37/1000G_phase1.snps.high_confidence.b37.sites"
   hapmap="/illumina/runs/Genome_b37/hapmap_3.3.b37.vcf"	
   omni="/illumina/runs/Genome_b37/1000G_omni2.5.b37.vcf"

# exon target intervals
   target_intervals="/illumina/runs/Genome_b37/Agilent_AllExon_V5_SureSelect_b37.list"


####################################################################################################################################################################################################################

# Step 1: Converting SAM to BAM

# SamtoBam
   echo "SamtoBam, Samtools"
    $samtools view -bT $ref_b37 -o ${fastq_id}_aln.bam ${fastq_id}.sam.gz
    sleep 10s  

    if [ $? -eq 0 ]; then
    rm ${id}.sam.gz
    fi  

# SortSam
   echo "SortSam, Picard" 
    /illumina/pipeline/jre1.8.0/bin/java -jar $picard SortSam I=${fastq_id}_aln.bam O=${fastq_id}_sorted.bam SO=coordinate CREATE_INDEX=true
    sleep 10s  

    if [ $? -eq 0 ]; then
    rm ${fastq_id}_aln.bam
    fi

    if [ $? -eq 0 ]; then
    rm ${fastq_id}_aln.bai
    fi

####################################################################################################################################################################################################################

## Step 2: Duplicate marking

# Mark Duplicates
   echo "Mark Duplicates, Picard" 
    /illumina/pipeline/jre1.8.0/bin/java -jar $picard MarkDuplicates I=${fastq_id}_sorted.bam O=${fastq_id}_mkdup.bam M=${fastq_id}_markdup_metrics.txt VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true 
    sleep 10s 

    if [ $? -eq 0 ]; then
    rm ${fastq_id}_sorted.bam
    fi

    if [ $? -eq 0 ]; then
    rm ${fastq_id}_sorted.bai
    fi

####################################################################################################################################################################################################################

# Step 3: Local realignment

# Indel Target Creator
   echo "Indel Target Creator, GATK" 
    /illumina/pipeline/jre1.8.0/bin/java -jar $gatk -T RealignerTargetCreator -R $ref_b37 -known $phase1_indels -known $mills_indels -I ${fastq_id}_mkdup.bam -o ${fastq_id}_realn.intervals #-L $target_intervals -ip 10
    sleep 10s 

# Indel Realignment 
   echo "Indel Realignment, GATK" 
    /illumina/pipeline/jre1.8.0/bin/java -jar $gatk -T IndelRealigner -R $ref_b37 -known $phase1_indels -known $mills_indels -targetIntervals ${fastq_id}_realn.intervals -I ${fastq_id}_mkdup.bam -o ${fastq_id}_mkdup_realn.bam -model USE_READS
    sleep 10s

    if [ $? -eq 0 ]; then
    rm ${fastq_id}_mkdup.bam
    fi

    if [ $? -eq 0 ]; then
    rm ${fastq_id}_mkdup.bai
    fi

####################################################################################################################################################################################################################

# Step 4: Base quality score recalibration

# Base Quality Score Recalibration (BQSR)
   echo "BQSR Baserecalibrator, GATK" 
    /illumina/pipeline/jre1.8.0/bin/java -jar $gatk -T BaseRecalibrator -R $ref_b37 -I ${fastq_id}_mkdup_realn.bam -knownSites $dbsnp_138 -knownSites $phase1_indels -knownSites $mills_indels -o ${fastq_id}_mkdup_realn_recal.table
    sleep 10s 

   echo "BQSR Baserecalibrator, GATK" 
    /illumina/pipeline/jre1.8.0/bin/java -jar $gatk -T BaseRecalibrator -R $ref_b37 -I ${fastq_id}_mkdup_realn.bam -knownSites $dbsnp_138 -knownSites $phase1_indels -knownSites $mills_indels -BQSR ${fastq_id}_mkdup_realn_recal.table -o ${fastq_id}_mkdup_realn_after_recal.table
    sleep 10s 

   echo "BQSR Analyze Covariates, GATK" 
    /illumina/pipeline/jre1.8.0/bin/java -jar $gatk -T AnalyzeCovariates -R $ref_b37 -before ${fastq_id}_mkdup_realn_recal.table -after ${fastq_id}_mkdup_realn_after_recal.table -plots ${fastq_id}_recal_plots.pdf
    sleep 10s 

   echo "BQSR Print Reads, GATK" 
    /illumina/pipeline/jre1.8.0/bin/java -jar $gatk -T PrintReads -R $ref_b37 -I ${fastq_id}_mkdup_realn.bam -BQSR ${fastq_id}_mkdup_realn_recal.table -o ${fastq_id}_improved.bam
    sleep 10s 

    if [ $? -eq 0 ]; then
    rm ${fastq_id}_mkdup_realn.bam
    fi

    if [ $? -eq 0 ]; then
    rm ${fastq_id}_mkdup_realn.bai
    fi

####################################################################################################################################################################################################################

# Step 5: Variant calling

# Variant calling by HaplotypeCaller in GVCF mode
   echo "Variant calling by HaplotypeCaller in GVCF mode" 
    /illumina/pipeline/jre1.8.0/bin/java -jar $gatk -T HaplotypeCaller -R $ref_b37 -I ${fastq_id}_improved.bam --emitRefConfidence GVCF --dbsnp $dbsnp_138 -o ${fastq_id}.raw.snps.indels.g.vcf 

echo "Analysis finished: $(date)"

####################################################################################################################################################################################################################


