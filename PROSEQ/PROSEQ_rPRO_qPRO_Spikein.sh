#!/usr/bin/bash

# *******************************************
# Date: 202411 Wangshaoxuan
# Description: Pipeline for PRO-seq
# 1. Spike-in: With or without spike-in
# 2. qPRO-seq protocol (duplex 6 nt UMI) or rPRO-seq protocol (5'/read2 single n nt UMI)
# 3. Alignment: Bowtie2 (there is no consensus on which aligner works best with PRO-seq)
# 4. Peak calling: dREG, PINTS and HOMER (NRSA)
# 5. QC metrics modified from TG Scott et.al 2022 bioRxiv and RF Sigauke et al. 2023 bioRxiv
# *******************************************

# usage
# nohup ./xxx.sh /chenfeilab/Gaux/rawDataBackup/xxx/xxxxxx_PRO-seq sampleInfo.txt 241108_PROSEQ hg38 N rPRO 6 all &> ./241105_PROseq.log &

# sampleInfo.txt file
# XXX_001 XXX_PRO-seq_220308_DLD1-PNUTS-dTAG-3h-rep1
# XXX_002 XXX_PRO-seq_220319_DLD1-PNUTS-dTAG-3h-rep2

usage() {
    echo "Usage: $0 <rawDataRawDir> <sampleInfo> <runInfo> <spikeIn> <expRef> <libType> <umiLen> <callPeakMethod>"
    echo "  rawDataRawDir: raw data directory"
    echo "  sampleInfo: space separated sample information file"
    echo "  runInfo: run information, such as date, owner, etc."
    echo "  spikeIn: if spike-in PRO-seq library, 'Y' or 'N' "
    echo "  expRef: experiment genome reference, 'hg38', 'mm10' "
    echo "  libType: PRO-seq library type, 'qPRO' or 'rPRO' "
    echo "  umiLen: the length (n nt) of UMI, between 6 and 12"
    echo "  callPeakMethod: peak calling method, 'dREG', 'PINTS', 'HOMER' or 'all' "

    exit 1
}

# args
rawDataRawDir=$1
sampleInfo=$2
runInfo=$3
spikeIn=$4
expRef=$5
libType=$6
umiLen=$7
callPeakMethod=$8

# check args
if [ $# -ne 8 ]; then
    echo "Error: Invalid number of arguments!"
    usage
fi
if [[ "$rawDataRawDir" == *'*'* ]]; then
    baseDir="$rawDataRawDir"
    while [[ "$baseDir" == *'*'* ]]; do
        baseDir="${baseDir%/*}"
    done
    if [[ -d "$baseDir" && $(find $rawDataRawDir 2>/dev/null | head -n 1) ]]; then
        echo "Using $rawDataRawDir as input data."
    else
        echo "Error: Invalid or empty rawDataRawDir directory!"
	    usage
    fi
else
    if [[ -d "$rawDataRawDir" && $(ls -A "$rawDataRawDir") ]]; then
        echo "Using $rawDataRawDir as input data."
    else
        echo "Error: Invalid or empty rawDataRawDir directory!"
	    usage
    fi
fi
if [[ ! -s "$sampleInfo" ]];then
    echo "Error: $sampleInfo is not a valid file!"
    usage
fi
if [[ "$spikeIn" != "Y" && "$spikeIn" != "N" ]]; then
    echo "Error: Invalid spikeIn value! Must be 'Y' or 'N'."
    usage
fi
if [[ "$expRef" != "hg38" && "$expRef" != "hg19" && "$expRef" != "mm10" ]]; then
    echo "Error: Invalid expRef value! Must be 'hg38', 'hg19', or 'mm10'."
    usage
fi
if [[ "$libType" != "rPRO" && "$libType" != "qPRO" ]]; then
    echo "Error: Invalid libType value! Must be 'rPRO' or 'qPRO'."
    usage
fi
if [[ ! ($umiLen -ge 6 && $umiLen -le 12) ]]; then
    echo "Error: Invalid umiLen value! Must between 6 and 12."
    usage
fi
valid_methods=("dREG" "PINTS" "HOMER" "all")
if [[ ! " ${valid_methods[@]} " =~ " $callPeakMethod " ]]; then
    echo "Error: Invalid callPeakMethod value! Must be 'dREG', 'PINTS', 'HOMER' or 'all'."
    usage
fi
echo -e "PRO-seq parameters: \nrawDataRawDir: $rawDataRawDir\nsampleInfo: $sampleInfo\nrunInfo: $runInfo\nspikeIn: $spikeIn\nexpRef: $expRef\nlibType: $libType\numiLen: $umiLen\ncallPeakMethod: $callPeakMethod"

# library struceture (DNA): 
# qPRO: { 5' - [ P7 -- i7 (6 base variable) -- GTGACTGGAGTT -- CCTTGGCACCCGAGAATTCCA -- UMI (6 nt) -- C ] -- [insert] -- [ T -- UMI (6nt) -- GATCGTCGGACTGTAGAACTCTGAAC -- i5 (8 base variable) -- P5 ] -- 3' }
# rPRO: { 5' - [ P7 -- i7 (6 base variable) -- GTGACTGGAGTT -- CCTTGGCACCCGAGAATTCCA -- UMI (n nt) -- C ] -- [insert] -- [ GATCGTCGGACTGTAGAACTCTGAAC -- i5 (8 base variable) -- P5 ] -- 3' }
# R1 adaptor:  TGGAATTCTCGGGTGCCAAGG -- AACTCCAGTCAC (reverse complement of GTGACTGGAGTT -- CCTTGGCACCCGAGAATTCCA)
# R2 adaptor : GATCGTCGGACTGTAGAACTCTGAAC -- i5 (GGAACTTA) -- TGTAGATCTCGGTGGTCGCCGTATCATT
ADAPTOR_R1="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC" 
ADAPTOR_R2="GATCGTCGGACTGTAGAACTCTGAAC"
MAPQ=10

# output dirs
workDir=`pwd`
rawDataDir=${workDir}/00_raw_data
rawQcDir=${workDir}/01-1_raw_data_fastqc
trimDir=${workDir}/01-2_cutadapt_trimmed_data
trimFastpDir=${workDir}/01-3_fastp_labelUMI
rmrRNADir=${workDir}/01-4_remove_rRNA_reads
cleanQcDir=${workDir}/01-5_clean_data_fastqc
map2ExpDir=${workDir}/02-1_experimental_genome_alignment
map2SpkDir=${workDir}/02-2_spikein_genome_alignment

# log dirs
logDir=${workDir}/logs
trimLogDir=${logDir}/cutadapt_trim_log
trimFastpLogDir=${logDir}/fastp_label_UMI_log
rmrRNALogDir=${logDir}/bowtie2_remove_rRNA_log
map2ExpLogDir=${logDir}/bowtie2_experimental_genome_alignment_log
map2SpkLogDir=${logDir}/bowtie2_spikein_genome_alignment_log

# QC dirs
insertLenDir=${logDir}/flash_estimate_insert_length


# all tools in proseq env update at Oct 2024
# # fastqc version 0.12.1
fastqc='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/bin/fastqc'
# # cutadapt version 4.9
cutadapt='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/bin/cutadapt'
# # multiqc version 1.25.1
multiqc='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/bin/multiqc'
# # fastp version 0.23.4
fastp='/share/home/Grape/software/install_pkg/miniconda3/envs/proseq/bin/fastp'
# # bowtie2 version 2.5.4
bowtie2='/share/home/Grape/software/install_pkg/miniconda3/envs/proseq/bin/bowtie2'
# # umi_tools version 1.1.6
umi_tools='/share/home/Grape/software/install_pkg/miniconda3/envs/proseq/bin/umi_tools'
# # samtools version 1.21
samtools='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/download/samtools/bin/samtools'
# # deepTools version 3.5.5
bamCoverage='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/bin/bamCoverage'
# # dREG version

# # PINTS verson 1.1.14
pints_caller='/share/home/Grape/software/install_pkg/miniconda3/envs/proseq/bin/pints_caller'
# # HOMER

# # flash version 1.2.11
flash='/share/home/Grape/software/install_pkg/FLASH-1.2.11/flash'

# # R script
Rscript='/opt/ohpc/pub/apps/R/4.2.2/bin/Rscript'
get_degradation_ratio='/gpfs/chenfeilab/Gaux/PanCaNascent/02_Alignment/get_degradation_ratio.R'


# PROseq bowtie2 index
hg38Index='/share/home/Grape/genome/Homo_sapiens/bowtie2_index/hg38/hg38'
hg19Index='/share/home/Blueberry/reference/index/bowtie2/hg19/hg19'
mm10Index='/share/home/Grape/genome/Mus_musculus/bowtie2_index/mm10/mm10'
rDNAIndexHuman='/share/home/Grape/genome/Homo_sapiens/rDNA/rDNA_bowtie2_index/human_rDNA'
rDNAIndexMouse='/share/home/Grape/genome/Mus_musculus/rDNA/rDNA_bowtie2_index/mouse_rDNA'

# genome ref
if [[ $spikeIn == 'Y' ]]; then
    case $expRef in
        "hg38")
            expBowtie2Index=${hg38Index}
            spkBowtie2Index=${mm10Index}
            rDNABowtie2Index=${rDNAIndexHuman}
            exp_info="hg38"
            spike_info="mm10"
            ;;
        "hg19")
            expBowtie2Index=${hg19Index}
            spkBowtie2Index=${mm10Index}
            rDNABowtie2Index=${rDNAIndexHuman}
            exp_info="hg19"
            spike_info="mm10"
            ;;
        "mm10")
            expBowtie2Index=${mm10Index}
            spkBowtie2Index=${hg38Index}
            rDNABowtie2Index=${rDNAIndexMouse}
            exp_info="mm10"
            spike_info="hg38"
            ;;
        *)
            exit 1
            ;;
    esac
else
    spike_info=''
    case $expRef in
        "hg38")
            expBowtie2Index=${hg38Index}
            rDNABowtie2Index=${rDNAIndexHuman}
            exp_info="hg38"
            ;;
        "hg19")
            expBowtie2Index=${hg19Index}
            rDNABowtie2Index=${rDNAIndexHuman}
            exp_info="hg19"
            ;;
        "mm10")
            expBowtie2Index=${mm10Index}
            rDNABowtie2Index=${rDNAIndexMouse}
            exp_info="mm10"
            ;;
        *)
            exit 1
            ;;
    esac
fi


echo -e "\n***************************\nPRO-seq processing at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
echo -e "experimental genome is: ${exp_info} \nspike-in genome is: ${spike_info}"

####
# Step0: link raw data and rename
# WARNING: File name must have same suffix pattern!
echo -e "\n***************************\nRenaming files at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
mkdir $logDir
if [[ ! -d $rawDataDir ]]; then
	mkdir -p $rawDataDir
	cat $sampleInfo | while read file; do
        arr=($file)
        old=${arr[0]}
        new=${arr[1]}
        full=`ls ${rawDataRawDir}/* | grep $old`
        suffix=$(echo -e "$full" | rev | sed -e '1{h;d;}' -e 'G;s,\(.*\).*\n\1.*,\1,;h;$!d' | rev)
        r1=`ls ${rawDataRawDir}/*1${suffix} | grep $old`
        r2=`ls ${rawDataRawDir}/*2${suffix} | grep $old`
        ln -s $r1 ${rawDataDir}/${new}_R1.fq.gz
        ln -s $r2 ${rawDataDir}/${new}_R2.fq.gz
        echo -e "${arr[@]} $suffix" >> ${logDir}/sampleInfo_new.txt
	done
    echo -e "Finish renaming"
else
	:
fi


####
# Step1: QC
echo -e "\n***************************\nTrimming and QC at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

# Step1.1: raw data fastqc
if [[ ! -d $rawQcDir ]]; then
    mkdir -p $rawQcDir
	$fastqc -t 48 --memory 1024 ${rawDataDir}/*.fq.gz -o $rawQcDir &> /dev/null
fi
if [[ ! -s ${rawQcDir}/rawdata_multiqc.html ]] || [[ ! -d ${rawQcDir}/rawdata_multiqc_data ]]; then
	$multiqc -f -n ${runInfo}_rawdata_multiqc -o $rawQcDir $rawQcDir
fi

# Step1.2: trim adapter and low quality sequence (cutadapt)
# Use cutadapt but not trim_galore to clip more adapter

if [[ ! -d $trimDir ]]; then
    mkdir -p $trimDir
	mkdir -p $trimLogDir
    echo -e "\n***************************\nStart run cutadapt at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
	for r1 in `ls ${rawDataDir}/*_R1.fq.gz`; do
		r2=${r1/R1.fq.gz/R2.fq.gz}
		sampleName=$(basename ${r1%_R1.fq.gz})
        $cutadapt \
            -a $ADAPTOR_R1 -a "TGGAATTCTCGG" \
            -A $ADAPTOR_R2 -A "GATCGTCGGACT" \
            --times 2 \
            --cores 48 \
            --quality-cutoff 20 \
            --minimum-length 15 \
            --error-rate 0.15 \
            --overlap 5 \
            --pair-filter=any \
            --compression-level=1 \
            --action=trim \
            --nextseq-trim=20 \
            -o ${trimDir}/${sampleName}_trimmed_R1.fq.gz -p ${trimDir}/${sampleName}_trimmed_R2.fq.gz \
            $r1 $r2 &> ${trimLogDir}/${sampleName}_cutadapt.log
        echo -e "Finish cutadapt trim for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
	done
    echo -e "\nFinish run cutadapt for all samples!"
fi


# QC 1: estimate insert size by FLASH

if [[ ! -d $insertLenDir ]]; then
    mkdir -p $insertLenDir
    echo -e "\n***************************\nStart estimate insert size at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    if [ ! -s ${insertLenDir}/*.hist ]; then
        for r1 in `ls ${trimDir}/*_trimmed_R1.fq.gz`; do
            r2=${r1/R1.fq.gz/R2.fq.gz}
            sampleName=$(basename ${r1%_trimmed_R1.fq.gz})
            # eastimate insert size by flash
            $flash --quiet --thread 36 --compress-prog=pigz --compress-prog-args='--processes 36' --suffix=gz -d $insertLenDir -o $sampleName $r1 $r2
            echo -e "Finish estimate insert for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
        done
    fi
    # plot insert size distribution and degradation ratio like
    echo -e "\nFinish for all samples, remove intermediate files and plot insert size distributions...\n"
    ls ${insertLenDir}/* | grep 'fastq.gz$' | while read file; do unlink $file; done
    $Rscript $get_degradation_ratio $insertLenDir $runInfo $libType $umiLen
fi


# Step1.3: UMI labeling and clipping (fastp)

fastp_label_umi() {
    local read1=$1               # -i
    local read2=$2               # -I
    local out_prefix=$3          # -o and -O
    local adapter_r1=$4          # --adapter_sequence
    local adapter_r2=$5          # --adapter_sequence_r2
    local umi_loc=$6             # --umi_loc
    local log_file_prefix=$7     # --html and logs
    local extra_params=$8        

    # remove R2 umi "C" (and R1 umi "A"(only for qPRO)) by skip 1 base
    # Use 15bp as threshold
    $fastp \
        -i $read1 -I $read2 \
        -o ${out_prefix}_fastp_R1.fq.gz -O ${out_prefix}_fastp_R2.fq.gz \
        --adapter_sequence $adapter_r1 --adapter_sequence_r2 $adapter_r2 \
        --length_required 15 \
        --trim_poly_g \
        --umi \
        --umi_prefix="UMI" \
        --umi_skip=1 \
        --umi_loc=$umi_loc \
        --umi_len=$umiLen \
        --html ${log_file_prefix}_fastp.html \
        --json ${log_file_prefix}_fastp.json \
        --thread 24 \
        --correction \
        --overlap_len_require 15 \
        $extra_params 2> ${log_file_prefix}_fastp.log

    # echo -e "$read1\t$read2\t$out_prefix\t$adapter_r1\t$adapter_r2\t$umi_loc\t$log_file_prefix\t$extra_params"
}

# qPRO have 6 bp R1 and R2 duplex UMI
# rPRO only have R2 n bp UMI+C on the prefix/head of R2 (single UMI)
# sequencing read-through always happen on single UMI rPRO library, so here we clip R1 reverse complementary R2 UMI+G (n+1 bp) by `--trim_tail1` and `--trim_tail2` to keep consistent

if [[ ! -d $trimFastpDir ]]; then
    mkdir -p ${trimFastpDir}
    mkdir -p ${trimFastpLogDir}
    echo -e "\n***************************\nLabeling UMI at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    if [ ! -s ${trimFastpDir}/*_fastp_R1.fq.gz ]; then
        for r1 in `ls ${trimDir}/*_trimmed_R1.fq.gz`; do
            r2=${r1/R1.fq.gz/R2.fq.gz}
            sampleName=$(basename ${r1%_trimmed_R1.fq.gz})
            outPrefix=${trimFastpDir}/${sampleName}
            logName=${trimFastpLogDir}/${sampleName}
            if [[ $libType == 'qPRO' ]]; then
                # always 6nt duplex UMI
                fastp_label_umi $r1 $r2 $outPrefix $ADAPTOR_R1 $ADAPTOR_R2 "per_read" $logName ""
            elif [[ $libType == 'rPRO' ]]; then
                # clip R1 reverse complementary R2 UMI+G (n+1 bp)
                fastp_label_umi $r1 $r2 $outPrefix $ADAPTOR_R1 $ADAPTOR_R2 "read2" $logName "--trim_tail1=$(echo "${umiLen}+1" | bc -l) --trim_tail2=0"
            else
                exit 1
            fi
            echo -e "Finish label UMI for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
        done
    fi
fi

# Step1.4: remove rRNA reads
run_bowtie2_alignment() {

    local genome_dir=$1          # -x
    local read1=$2               # -1 and -2
    local read2=$3
    local extra_params=$4
    local log_file_prefix=$5

    echo -e "Run bowtie2 with:\n-x $genome_dir -1 $read1 -2 $read2 --threads 48 --local $extra_params\n" > ${log_file_prefix}_bowtie2.log
    
    # use local mode to allow soft clip
    $bowtie2 -x $genome_dir \
        -1 $read1 -2 $read2 \
        --threads 48 \
        --local \
        $extra_params \
        2>> ${log_file_prefix}_bowtie2.log

}

# remove rRNA reads
if [[ ! -d $rmrRNADir ]]; then
    mkdir -p ${rmrRNADir}
    mkdir -p ${rmrRNALogDir}
    echo -e "\n***************************\nRemoving rRNA reads at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    for r1 in `ls ${trimFastpDir}/*_fastp_R1.fq.gz`; do
        r2=${r1/fastp_R1.fq.gz/fastp_R2.fq.gz}
        sampleName=$(basename ${r1/_fastp_R1.fq.gz/})
        logName=${rmrRNALogDir}/${sampleName}_rmrRNA
        if [ ! -s ${rmrRNALogDir}/${sampleName}_rmrRNA_bowtie2.log ]; then
            # We are interested primarily in quickly identifying and removing any reads that have a valid alignment to the serial alignment genome (-k 1 parameter)
            run_bowtie2_alignment $rDNABowtie2Index $r1 $r2 "--fast-local -N 1 -L 18 -k 1 --un-conc-gz ${rmrRNADir}/${sampleName}_rmrRNA_R%.fq.gz" $logName > /dev/null
        fi
        echo -e "Finish remove rRNA reads for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
    done
fi

# Step1.5: clean data fastqc
if [[ ! -d $cleanQcDir ]]; then
	mkdir -p $cleanQcDir
    echo -e "\n***************************\nSummary QC at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
	$fastqc -t 48 --memory 1024 ${rmrRNADir}/*.fq.gz -o $cleanQcDir &> /dev/null
fi
if [[ ! -s ${cleanQcDir}/cleanData_multiqc.html ]] || [[ ! -d ${cleanQcDir}/cleanData_multiqc_data ]]; then
	$multiqc -f -n ${runInfo}_cleanData_multiqc -o $cleanQcDir $trimLogDir $trimFastpLogDir $rmrRNALogDir $cleanQcDir 
fi

# Step2: Alignment
# Step2.1: align to experimental
# !! Set bowtie2 seed length to 18 (-L 18) to rescue more reads for highly degraded libarary (like tumor sample...)
# peak calling methods like dREG needs enriched reads to detect the transcriptional peaks, so we need more reads: https://github.com/Danko-Lab/proseq2.0/tree/master?tab=readme-ov-file#notes-for-dreg-users

echo -e "\n***************************\nAligning to experimental at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
if [[ ! -d $map2ExpDir ]]; then
    mkdir -p $map2ExpDir
	mkdir -p $map2ExpLogDir
    if [[ ! `ls ${map2ExpDir}/*_${exp_info}.bam 2> /dev/null` ]]; then
        for r1 in `ls ${rmrRNADir}/*_rmrRNA_R1.fq.gz`;do
            r2=${r1/R1.fq.gz/R2.fq.gz}
            sampleName=$(basename ${r1%_rmrRNA_R1.fq.gz})
            logName=${map2ExpLogDir}/${sampleName}_${exp_info}
            if [ ! -s ${map2ExpLogDir}/${sampleName}_${exp_info}_bowtie2.log ]; then
                run_bowtie2_alignment $expBowtie2Index $r1 $r2 "--very-sensitive-local -X 750 -L 18" $logName | \
                $samtools view -@ 48 -bS -q ${MAPQ} -f 2 | \
                $samtools sort -@ 48 -o ${map2ExpDir}/${sampleName}_${exp_info}.bam 2> /dev/null
                $samtools index -@ 48 ${map2ExpDir}/${sampleName}_${exp_info}.bam
            fi
            echo -e "Finish align to ${exp_info} for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
        done
    fi
fi

# Step2.2: align to spike-in
if [[ $spikeIn == 'Y' ]];then
    echo -e "\n***************************\nAligning to spike-in at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    if [[ ! -d $map2SpkDir ]];then
        mkdir -p $map2SpkDir
        mkdir -p $map2SpkLogDir
        if [[ ! `ls ${map2SpkDir}/*_${spike_info}.bam 2> /dev/null` ]]; then
            for r1 in `ls ${rmrRNADir}/*_rmrRNA_R1.fq.gz`;do
                r2=${r1/R1.fq.gz/R2.fq.gz}
                sampleName=$(basename ${r1%_rmrRNA_R1.fq.gz})
                logName=${map2SpkLogDir}/${sampleName}_${spike_info}
                if [ ! -s ${map2SpkLogDir}/${sampleName}_${spike_info}_bowtie2.log ]; then
                    run_bowtie2_alignment $spkBowtie2Index $r1 $r2 "--very-sensitive-local -X 750 -L 18" $logName | \
                    $samtools view -@ 48 -bS -q ${MAPQ} -f 2 | \
                    $samtools sort -@ 48 -o ${map2SpkDir}/${sampleName}_${spike_info}.bam 2> /dev/null
                    $samtools index -@ 48 ${map2SpkDir}/${sampleName}_${spike_info}.bam
                fi
                echo -e "Finish align to ${spike_info} for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
            done
        fi
    fi
fi
