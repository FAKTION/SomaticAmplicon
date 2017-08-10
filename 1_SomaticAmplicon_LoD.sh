#!/bin/bash
#PBS -m abe
#PBS -M meissnerM1@cardiff.ac.uk
#PBS -N Breast_cancer_pipeline
#PBS -q workq_wgp
#PBS -P PR403
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=12
set -euo pipefail
cd $PBS_O_WORKDIR

#Description: Somatic Amplicon Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_SAMPLE
version="1.6.4"

# Load Modules
module load python/2.7.11-genomics
module load fastqc/0.11.2
module load java/1.8.0_45
module load GATK/3.7
module load picard/2.7.1
module load BWA/0.7.15
module load ampliconrealigner/1.1.1
module load samtools/1.3.1
module load mono/4.4.1
module load softclippcrprimer/1.1.0
module load ActivePerl/5.18
module load pear/0.9.10
module load bedtools/2.26.0
module load vcfparse/1.2.5
module load ensembl-vep/89.6
module load ensemble_tools/86
module load bcftools/1.2
module load htslib/1.2.1

# Directory structure required for pipeline
#
# /data
# └── results
#     └── seqId
#         ├── panel1
#         │   ├── sample1
#         │   ├── sample2
#         │   └── sample3
#         └── panel2
#             ├── sample1
#             ├── sample2
#             └── sample3
#
# Script 1 runs in sample folder, requires fastq files split by lane

#load sample & pipeline variables
. *.variables
. /scratch/mcgmm/Matt_pipeline/data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel".variables

#downsample BAM to p
for p in 0.9 0.75 0.6 0.5 0.4 0.25 0.125 0.0625 0.01; do

    #make folder for dowsampled analysis
    mkdir "$p"
    cd "$p"

    #Downsample BAM
    java -Xmx8g -jar /software/genomics/picard/2.7.1/picard.jar DownsampleSam \
    I=../"$seqId"_"$sampleId".bam \
    O="$seqId"_"$sampleId".bam \
    CREATE_INDEX=true \
    PROBABILITY="$p" \
    STRATEGY=HighAccuracy \
    MAX_RECORDS_IN_RAM=2000000

    ### Variant calling ###

    #make bai alias for Pisces
    ln -s "$seqId"_"$sampleId".bai "$seqId"_"$sampleId".bam.bai

    #extract thick regions
    awk '{print $1"\t"$7"\t"$8}' /scratch/mcgmm/Matt_pipeline/data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed | \
    bedtools merge > "$panel"_ROI_b37_thick.bed

    #Call somatic variants
    mono /scratch/mcgmm/Matt_pipeline/data/db/MiSeqReporter-2.6.3/CallSomaticVariants.exe \
    -B ./"$seqId"_"$sampleId".bam \
    -g /scratch/mcgmm/Matt_pipeline/data/db \
    -f 0.00794 \
    -fo False \
    -b 21 \
    -q 100 \
    -c 50 \
    -s 0.5 \
    -a 20 \
    -F 30 \
    -gVCF False \
    -i false \
    -PhaseSNPs true \
    -MaxPhaseSNPLength 100 \
    -r .

    #fix VCF name
    echo "$sampleId" > name
    bcftools reheader \
    -s name \
    -o "$seqId"_"$sampleId"_fixed.vcf \
    $(echo "$seqId"_"$sampleId" | sed 's/_/-/g')_S999.vcf
    rm name

    #left align and trim variants
    java -Xmx40g -jar /software/genomics/GATK/3.7/GenomeAnalysisTK.jar \
    -T LeftAlignAndTrimVariants \
    -R /scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta \
    -o "$seqId"_"$sampleId"_left_aligned.vcf \
    -V "$seqId"_"$sampleId"_fixed.vcf \
    -L "$panel"_ROI_b37_thick.bed \
    -dt NONE

    #Annotate with GATK contextual information
    java -Xmx40g -jar /software/genomics/GATK/3.7/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R /scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta \
    -I "$seqId"_"$sampleId".bam \
    -V "$seqId"_"$sampleId"_left_aligned.vcf \
    -L "$panel"_ROI_b37_thick.bed \
    -o "$seqId"_"$sampleId"_left_aligned_annotated.vcf \
    -A BaseQualityRankSumTest -A ChromosomeCounts -A MappingQualityRankSumTest -A MappingQualityZero -A RMSMappingQuality \
    -dt NONE

    #Annotate with low complexity region length using mdust
    bcftools annotate \
    -a /scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.mdust.v34.lpad1.bed.gz \
    -c CHROM,FROM,TO,LCRLen \
    -h <(echo '##INFO=<ID=LCRLen,Number=1,Type=Integer,Description="Overlapping mdust low complexity region length (mask cutoff: 34)">') \
    -o "$seqId"_"$sampleId"_lcr.vcf \
    "$seqId"_"$sampleId"_left_aligned_annotated.vcf

    #Filter variants
    java -Xmx40g -jar /software/genomics/GATK/3.7/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R /scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta \
    -V "$seqId"_"$sampleId"_lcr.vcf \
    --filterExpression "LCRLen > 8" \
    --filterName "LowComplexity" \
    --filterExpression "DP < 50" \
    --filterName "LowDP" \
    -L "$panel"_ROI_b37_thick.bed \
    -o "$seqId"_"$sampleId"_filtered.vcf \
    -dt NONE

    #Generate per-base coverage: variant detection sensitivity
    java -Xmx40g -jar /software/genomics/GATK/3.7/GenomeAnalysisTK.jar \
    -T DepthOfCoverage \
    -R /scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta \
    -o "$seqId"_"$sampleId"_DepthOfCoverage \
    -I "$seqId"_"$sampleId".bam \
    -L "$panel"_ROI_b37_thick.bed \
    --countType COUNT_FRAGMENTS \
    --minMappingQuality 20 \
    --minBaseQuality 20 \
    --omitIntervalStatistics \
    -ct "$minimumCoverage" \
    -nt 12 \
    -dt NONE
    
    #tabix index the per-base coverage file
    awk -F'[\t|:]' '{if(NR>1) print $1"\t"$2"\t"$3}' "$seqId"_"$sampleId"_DepthOfCoverage | \
    /software/genomics/htslib/1.2.1/gnu-4.4.7/bin/bgzip > "$seqId"_"$sampleId"_DepthOfCoverage.gz
    /software/genomics/htslib/1.2.1/gnu-4.4.7/bin/tabix -b2 -e2 -s1 "$seqId"_"$sampleId"_DepthOfCoverage.gz

    #move to root folder
    cd ..

done