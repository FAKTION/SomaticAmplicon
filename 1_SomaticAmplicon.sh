#!/bin/bash
#PBS -m abe
#PBS -M meissnerM1@cardiff.ac.uk
#PBS -N Breast_cancer_pipeline
#PBS -q workq_wgp
#PBS -P PR403
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=12
set -euo pipefail
#cd $PBS_O_WORKDIR

#Description: Somatic Amplicon Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_SAMPLE
version="1.6.3"

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
module load ensemble_tools/86
module load bcftools/1.2

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

countQCFlagFails() {
    #count how many core FASTQC tests failed
    grep -E "Basic Statistics|Per base sequence quality|Per tile sequence quality|Per sequence quality scores|Per base N content" "$1" | \
    grep -v ^PASS | \
    grep -v ^WARN | \
    wc -l | \
    sed 's/^[[:space:]]*//g'
}

#load sample & pipeline variables
. *.variables
. /scratch/mcgmm/Matt_pipeline/data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel".variables

### Preprocessing ###

#record FASTQC pass/fail
rawSequenceQuality=PASS

#convert FASTQ to uBAM & add RGIDs
for fastqPair in $(ls "$sampleId"_S*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do

    #parse fastq filenames
    laneId=$(echo "$fastqPair" | cut -d_ -f3)
    read1Fastq=$(ls "$fastqPair"_R1_*fastq.gz)
    read2Fastq=$(ls "$fastqPair"_R2_*fastq.gz)

    #trim adapter
    /software/tools/pyenv/versions/2.7.11-genomics/bin/cutadapt \
    -a "$read1Adapter" \
    -A "$read2Adapter" \
    -m 50 \
    -o "$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    -p "$seqId"_"$sampleId"_"$laneId"_R2.fastq \
    "$read1Fastq" \
    "$read2Fastq"

    #merge overlapping reads
    pear \
    -f "$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    -r "$seqId"_"$sampleId"_"$laneId"_R2.fastq \
    -o "$seqId"_"$sampleId"_"$laneId"_merged.fastq \
    -j 12

    #convert fastq to ubam
    java -Xmx8g -jar /software/genomics/picard/2.7.1/picard.jar FastqToSam \
    F1="$seqId"_"$sampleId"_"$laneId"_merged.fastq.assembled.fastq \
    O="$seqId"_"$sampleId"_"$laneId"_unaligned.bam \
    QUALITY_FORMAT=Standard \
    READ_GROUP_NAME="$seqId"_"$laneId"_"$sampleId" \
    SAMPLE_NAME="$sampleId" \
    LIBRARY_NAME="$worklistId"_"$sampleId"_"$panel" \
    PLATFORM_UNIT="$seqId"_"$laneId" \
    PLATFORM="ILLUMINA" \
    SEQUENCING_CENTER="IMG" \
    PREDICTED_INSERT_SIZE="$expectedInsertSize" \
    SORT_ORDER=queryname \
    MAX_RECORDS_IN_RAM=2000000 \

    #fastqc
    fastqc --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R1.fastq
    fastqc --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R2.fastq

    #check FASTQC output
    if [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R1_fastqc/summary.txt) -gt 0 ] || [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R2_fastqc/summary.txt) -gt 0 ]; then
        rawSequenceQuality=FAIL
    fi

    #clean up
    rm "$seqId"_"$sampleId"_"$laneId"_R1.fastq "$seqId"_"$sampleId"_"$laneId"_R2.fastq "$seqId"_"$sampleId"_"$laneId"_merged.fastq.*

done

#merge lane bams
java -Xmx8g -jar /software/genomics/picard/2.7.1/picard.jar MergeSamFiles \
$(ls "$seqId"_"$sampleId"_*_unaligned.bam | sed 's/^/I=/' | tr '\n' ' ') \
SORT_ORDER=queryname \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT \
USE_THREADING=true \
MAX_RECORDS_IN_RAM=2000000 \
O="$seqId"_"$sampleId"_unaligned.bam

#uBam2fq, map & MergeBamAlignment
java -Xmx8g -jar /software/genomics/picard/2.7.1/picard.jar SamToFastq \
I="$seqId"_"$sampleId"_unaligned.bam \
FASTQ=/dev/stdout \
NON_PF=true \
MAX_RECORDS_IN_RAM=2000000 \
VALIDATION_STRINGENCY=SILENT | \
bwa mem \
-M \
-t 12 \
-p \
/scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta \
/dev/stdin | \
java -Xmx8g -jar /software/genomics/picard/2.7.1/picard.jar MergeBamAlignment \
ATTRIBUTES_TO_RETAIN=X0 \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM="$seqId"_"$sampleId"_unaligned.bam \
OUTPUT="$seqId"_"$sampleId"_aligned.bam \
R=/scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta \
PAIRED_RUN=false \
SORT_ORDER="coordinate" \
IS_BISULFITE_SEQUENCE=false \
ALIGNED_READS_ONLY=false \
CLIP_ADAPTERS=false \
MAX_RECORDS_IN_RAM=2000000 \
MAX_INSERTIONS_OR_DELETIONS=-1 \
UNMAP_CONTAMINANT_READS=false \
CLIP_OVERLAPPING_READS=false \
ALIGNER_PROPER_PAIR_FLAGS=false \
ATTRIBUTES_TO_RETAIN=XS \
INCLUDE_SECONDARY_ALIGNMENTS=true \
CREATE_INDEX=true \

#Realign soft clipped bases
ampliconrealigner.sh \
-I "$seqId"_"$sampleId"_aligned.bam \
-O "$seqId"_"$sampleId"_amplicon_realigned.bam \
-R /scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta \
-T /scratch/mcgmm/Matt_pipeline/data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed

#sort and index BAM
samtools sort -@8 -m8G -o "$seqId"_"$sampleId"_amplicon_realigned_sorted.bam "$seqId"_"$sampleId"_amplicon_realigned.bam
samtools index "$seqId"_"$sampleId"_amplicon_realigned_sorted.bam

#left align indels
java -Xmx40g -jar /software/genomics/GATK/3.7/GenomeAnalysisTK.jar \
-T LeftAlignIndels \
-R /scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta \
-I "$seqId"_"$sampleId"_amplicon_realigned_sorted.bam \
-o "$seqId"_"$sampleId"_amplicon_realigned_left_sorted.bam \
-dt NONE

#Identify regions requiring realignment
java -Xmx40g -jar /software/genomics/GATK/3.7/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta \
-known /scratch/mcgmm/Matt_pipeline/data/db/1000G_phase1.indels.b37.vcf \
-known /scratch/mcgmm/Matt_pipeline/data/db/Mills_and_1000G_gold_standard.indels.b37.vcf \
-known /scratch/mcgmm/Matt_pipeline/data/db/cosmic_78.indels.b37.vcf \
-I "$seqId"_"$sampleId"_amplicon_realigned_left_sorted.bam \
-o "$seqId"_"$sampleId"_indel_realigned.intervals \
-L /scratch/mcgmm/Matt_pipeline/data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip "$padding" \
-nt 12 \
-dt NONE

#Realign around indels
java -Xmx40g -jar /software/genomics/GATK/3.7/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta \
-known /scratch/mcgmm/Matt_pipeline/data/db/1000G_phase1.indels.b37.vcf \
-known /scratch/mcgmm/Matt_pipeline/data/db/Mills_and_1000G_gold_standard.indels.b37.vcf \
-known /scratch/mcgmm/Matt_pipeline/data/db/cosmic_78.indels.b37.vcf \
-targetIntervals "$seqId"_"$sampleId"_indel_realigned.intervals \
--maxReadsForRealignment 500000 \
--maxConsensuses 750 \
--maxReadsForConsensuses 3000 \
--maxReadsInMemory 3750000 \
-LOD 0.4 \
-I "$seqId"_"$sampleId"_amplicon_realigned_left_sorted.bam \
-o "$seqId"_"$sampleId"_indel_realigned.bam \
-dt NONE

#soft clip PCR primers
softclip.sh \
-I "$seqId"_"$sampleId"_indel_realigned.bam \
-O "$seqId"_"$sampleId"_clipped.bam \
-T /scratch/mcgmm/Matt_pipeline/data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed

#sort and index BAM
samtools sort -@8 -m8G -o "$seqId"_"$sampleId"_clipped_sorted.bam "$seqId"_"$sampleId"_clipped.bam
samtools index "$seqId"_"$sampleId"_clipped_sorted.bam

#fix bam tags
java -Xmx8g -jar /software/genomics/picard/2.7.1/picard.jar SetNmMdAndUqTags \
I="$seqId"_"$sampleId"_clipped_sorted.bam \
O="$seqId"_"$sampleId".bam \
CREATE_INDEX=true \
IS_BISULFITE_SEQUENCE=false \
R=/scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta

### Variant calling ###

#make bai alias for Pisces
ln -s "$seqId"_"$sampleId".bai "$seqId"_"$sampleId".bam.bai

#extract thick regions
awk '{print $1"\t"$7"\t"$8}' /scratch/mcgmm/Matt_pipeline/data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed | \
bedtools merge > "$panel"_ROI_b37_thick.bed
bedtools merge > "$panel"_ROI_b37_thick.bed

#Call somatic variants
mono /scratch/mcgmm/Matt_pipeline/data/db/MiSeqReporter-2.6.3/CallSomaticVariants.exe \
-B ./"$seqId"_"$sampleId".bam \
-g /scratch/mcgmm/Matt_pipeline/data/db/gatk/2.8/b37 \    
-f 0.01 \
-fo False \
-b 20 \
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
/software/genomics/bcftools-1.3.1/bcftools annotate \
-a /scratch/mcgmm/Matt_pipeline/data/db/gatk/2.8/b37/human_g1k_v37.mdust.v34.lpad1.bed.gz \  
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

### QC ###

#Convert BED to interval_list for later
java -Xmx8g -jar /software/genomics/picard/2.7.1/picard.jar BedToIntervalList \
I="$panel"_ROI_b37_thick.bed \
O="$panel"_ROI.interval_list \
SD=/scratch/mcgmm/Matt_pipeline/data/db/gatk/2.8/b37/human_g1k_v37.dict  

#HsMetrics: capture & pooling performance
java -Xmx8g -jar /software/genomics/picard/2.7.1/picard.jar CollectHsMetrics \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_hs_metrics.txt \
R=/scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta \
BAIT_INTERVALS="$panel"_ROI.interval_list \
TARGET_INTERVALS="$panel"_ROI.interval_list

#Generate per-base coverage: variant detection sensitivity
java -Xmx40g -jar /software/genomics/GATK/3.7/GenomeAnalysisTK.jar \
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

#Calculate gene (clinical) percentage coverage
java -Xmx8g -jar /scratch/mcgmm/Matt_pipeline/data/db/CoverageCalculator-2.0.2/CoverageCalculator-2.0.2.jar \
"$seqId"_"$sampleId"_DepthOfCoverage \
/scratch/mcgmm/Matt_pipeline/data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_genes.txt \
/scratch/mcgmm/Matt_pipeline/data/db/refseq/ref_GRCh37.p13_top_level.gff3 \
-p5 \
-d"$minimumCoverage" \
> "$seqId"_"$sampleId"_PercentageCoverage.txt

#Gather QC metrics
totalReads=$(head -n8 "$seqId"_"$sampleId"_hs_metrics.txt | tail -n1 | cut -s -f6) #The total number of reads in the SAM or BAM file examine.
pctSelectedBases=$(head -n8 "$seqId"_"$sampleId"_hs_metrics.txt | tail -n1 | cut -s -f19) #On+Near Bait Bases / PF Bases Aligned.
totalTargetedUsableBases=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f2) #total number of usable bases.
meanOnTargetCoverage=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f3) #avg usable coverage
pctTargetBasesCt=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f7) #percentage panel covered with good enough data for variant detection

#Print QC metrics
echo -e "TotalReads\tRawSequenceQuality\tTotalTargetUsableBases\tPctSelectedBases\tPctTargetBasesCt\tMeanOnTargetCoverage" > "$seqId"_"$sampleId"_qc.txt
echo -e "$totalReads\t$rawSequenceQuality\t$totalTargetedUsableBases\t$pctSelectedBases\t$pctTargetBasesCt\t$meanOnTargetCoverage" >> "$seqId"_"$sampleId"_qc.txt

#Add VCF meta data to final VCF
grep '^##' "$seqId"_"$sampleId"_filtered.vcf > "$seqId"_"$sampleId"_filtered_meta.vcf
echo \#\#SAMPLE\=\<ID\="$sampleId",Tissue\=Somatic,WorklistId\="$worklistId",SeqId\="$seqId",Assay\="$panel",PipelineName\=SomaticAmplicon,PipelineVersion\="$version",RawSequenceQuality\="$rawSequenceQuality",TotalReads\="$totalReads",PctSelectedBases\="$pctSelectedBases",MeanOnTargetCoverage\="$meanOnTargetCoverage",PctTargetBasesCt\="$pctTargetBasesCt",TotalTargetedUsableBases\="$totalTargetedUsableBases",RemoteVcfFilePath\=$(find $PWD -type f -name "$seqId"_"$sampleId"_filtered_meta.vcf),RemoteBamFilePath\=$(find $PWD -type f -name "$seqId"_"$sampleId".bam)\> >> "$seqId"_"$sampleId"_filtered_meta.vcf
grep -v '^##' "$seqId"_"$sampleId"_filtered.vcf >> "$seqId"_"$sampleId"_filtered_meta.vcf

#Variant Evaluation
java -Xmx40g -jar /software/genomics/GATK/3.7/GenomeAnalysisTK.jar \
-T VariantEval \
-R /scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta \
-o "$seqId"_"$sampleId"_variant_evaluation.txt \
--eval:"$seqId"_"$sampleId" "$seqId"_"$sampleId"_filtered_meta.vcf \
--comp:omni2.5 /scratch/mcgmm/Matt_pipeline/data/db/gatk/2.8/b37/1000G_omni2.5.b37.vcf \
--comp:hapmap3.3 /scratch/mcgmm/Matt_pipeline/data/db/gatk/2.8/b37/hapmap_3.3.b37.vcf \
--comp:cosmic78 /scratch/mcgmm/Matt_pipeline/data/db/cosmic/b37/cosmic_78.b37.vcf \
-L "$panel"_ROI_b37_thick.bed \
-nt 12 \
-dt NONE

### Reporting ###

#annotate VCF with VEP
variant_effect_predictor.pl \
--verbose \
--no_progress \
--everything \
--fork 12 \
--species homo_sapiens \
--assembly GRCh37 \
--input_file "$seqId"_"$sampleId"_filtered_meta.vcf \
--format vcf \
--output_file "$seqId"_"$sampleId"_filtered_meta_annotated.vcf \
--force_overwrite \
--no_stats \
--cache \
--dir /scratch/mcgmm/Matt_pipeline/data/resources \
--fasta /scratch/mcgmm/Matt_pipeline/data/resources \
--no_intergenic \
--offline \
--cache_version 86 \
--allele_number \
--no_escape \
--shift_hgvs 1 \
--vcf \
--refseq

#check VEP has produced annotated VCF
if [ ! -e "$seqId"_"$sampleId"_filtered_meta_annotated.vcf ]; then
    cp "$seqId"_"$sampleId"_filtered_meta.vcf "$seqId"_"$sampleId"_filtered_meta_annotated.vcf
fi

#index & validate final VCF
java -Xmx40g -jar /software/genomics/GATK/3.7/GenomeAnalysisTK.jar \
-T ValidateVariants \
-R /scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta \          
-V "$seqId"_"$sampleId"_filtered_meta_annotated.vcf \
-dt NONE

#custom coverage reporting
if [ -d /scratch/mcgmm/Matt_pipeline/data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/hotspot_coverage ]; then
    mkdir hotspot_coverage
    echo -e "Target\tSampleId\tAverage\tPercentageAbove$minimumCoverage" > hotspot_coverage/"$seqId"_"$sampleId"_coverage_summary.txt

    for bedFile in $(ls /data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/hotspot_coverage/*.bed); do

        #extract target name
        target=$(basename "$bedFile" | sed 's/\.bed//g')

        #generate per-base coverage
        java -Xmx40g -jar /software/genomics/GATK/3.7/GenomeAnalysisTK.jar \
        -T DepthOfCoverage \
        -R /scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta \
        -o "$seqId"_"$sampleId"_"$target" \
        -I "$seqId"_"$sampleId".bam \
        -L "$bedFile" \
        --countType COUNT_FRAGMENTS \
        --minMappingQuality 20 \
        --minBaseQuality 20 \
        --omitIntervalStatistics \
        --omitLocusTable \
        -ct "$minimumCoverage" \
        -nt 12 \
        -dt NONE

        #extract low depth bases
        awk -v minimumCoverage="$minimumCoverage" '{ if(NR > 1 && $2 < minimumCoverage) {split($1,array,":"); print array[1]"\t"array[2]-1"\t"array[2]} }' "$seqId"_"$sampleId"_"$target" | \
        bedtools merge > hotspot_coverage/"$seqId"_"$sampleId"_"$target"_gaps.bed

        #calculate average coverage
        avg=$(awk '{if (NR > 1) n+= $2} END {print n /(NR-1)}' "$seqId"_"$sampleId"_"$target")

        #count bases above minimumCoverage
        pctAboveThreshold=$(awk -v minimumCoverage="$minimumCoverage" '{if (NR > 1 && $2 >= minimumCoverage) n++} END {print (n /(NR-1)) * 100}' "$seqId"_"$sampleId"_"$target")

        #write summary to file
        echo -e "$target\t$sampleId\t$avg\t$pctAboveThreshold" >> hotspot_coverage/"$seqId"_"$sampleId"_coverage_summary.txt

        rm "$seqId"_"$sampleId"_"$target".sample_statistics
        rm "$seqId"_"$sampleId"_"$target".sample_summary
        rm "$seqId"_"$sampleId"_"$target"

    done
fi

#custom variant reporting
  if [ -d /scratch/mcgmm/Matt_pipeline/data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/hotspot_variants ]; then
    mkdir hotspot_variants

    for bedFile in $(ls /data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/hotspot_variants/*.bed); do

        #extract target name
        target=$(basename "$bedFile" | sed 's/\.bed//g')

        #select variants
        java -Xmx40g -jar /software/genomics/GATK/3.7/GenomeAnalysisTK.jar \
        -T VariantFiltration \
        -R /scratch/mcgmm/Matt_pipeline/data/db/human_g1k_v37.fasta \
        -V "$seqId"_"$sampleId"_filtered_meta_annotated.vcf \
        -L "$bedFile" \
        -o hotspot_variants/"$seqId"_"$sampleId"_"$target"_filtered_meta_annotated.vcf \
        -dt NONE

        #write targeted dataset to table
        vcfparse.sh 
        -V hotspot_variants/"$seqId"_"$sampleId"_"$target"_filtered_meta_annotated.vcf \
        -T /scratch/mcgmm/Matt_pipeline/data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_PreferredTranscripts.txt \
        -C /scratch/mcgmm/Matt_pipeline/data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_KnownVariants.vcf \
        -K

        #move to hotspot_variants
        mv "$seqId"_"$sampleId"_VariantReport.txt hotspot_variants/"$seqId"_"$sampleId"_"$target"_VariantReport.txt

    done
fi

#write full dataset to table
vcfparse.sh \
-V "$seqId"_"$sampleId"_filtered_meta_annotated.vcf \
-T /scratch/mcgmm/Matt_pipeline/data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_PreferredTranscripts.txt \
-C /scratch/mcgmm/Matt_pipeline/data/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_KnownVariants.vcf \
-K

### Clean up ###

#delete unused files
rm "$seqId"_"$sampleId"_*unaligned.bam "$seqId"_"$sampleId"_aligned.bam "$seqId"_"$sampleId"_aligned.bai "$seqId"_"$sampleId"_amplicon_realigned.bam
rm "$seqId"_"$sampleId"_amplicon_realigned_sorted.bam "$seqId"_"$sampleId"_amplicon_realigned_sorted.bam.bai "$seqId"_"$sampleId"_indel_realigned.intervals
rm "$seqId"_"$sampleId"_clipped.bam "$seqId"_"$sampleId"_clipped_sorted.bam "$seqId"_"$sampleId"_clipped_sorted.bam.bai "$panel"_ROI.interval_list "$panel"_ROI_b37_thick.bed
rm "$seqId"_"$sampleId"_left_aligned.vcf "$seqId"_"$sampleId"_left_aligned.vcf.idx "$seqId"_"$sampleId".bam.bai "$seqId"_"$sampleId"_amplicon_realigned_left_sorted.bam
rm "$seqId"_"$sampleId"_amplicon_realigned_left_sorted.bai "$seqId"_"$sampleId"_filtered_meta.vcf "$seqId"_"$sampleId"_filtered_meta.vcf.idx "$seqId"_"$sampleId"_filtered.vcf
rm "$seqId"_"$sampleId"_filtered.vcf.idx "$seqId"_"$sampleId"_fixed.vcf "$seqId"_"$sampleId"_fixed.vcf.idx "$seqId"_"$sampleId"_indel_realigned.bam "$seqId"_"$sampleId"_indel_realigned.bai
rm "$seqId"_"$sampleId"_*_fastqc.zip "$seqId"_"$sampleId"_lcr.vcf "$seqId"_"$sampleId"_lcr.vcf.idx "$seqId"_"$sampleId"_left_aligned_annotated.vcf "$seqId"_"$sampleId"_left_aligned_annotated.vcf.idx
rm $(echo "$seqId"_"$sampleId" | sed 's/_/-/g')_S999.vcf
rm -r VariantCallingLogs