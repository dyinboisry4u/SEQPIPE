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
# nohup ./xxx.sh /chenfeilab/Gaux/rawDataBackup/xxx/xxxxxx_PRO-seq sampleInfo.txt 241108_PROSEQ N hg38 rPRO 6 none &> 241105_PROseq.log &
# nohup ./xxx.sh "/chenfeilab/Pomelo/try/SXY/24rawdata/*/*" sampleInfo.txt test_SXY51to60_qPRO Y mm10 qPRO 6 all &> test_SXY51to60_qPRO.log &

# sampleInfo.txt file
# XXX_001 XXX_PRO-seq_220308_DLD1-PNUTS-dTAG-3h-rep1
# XXX_002 XXX_PRO-seq_220319_DLD1-PNUTS-dTAG-3h-rep2

usage() {
    echo "Usage: $0 <rawDataRawDir> <sampleInfo> <runInfo> <spikeIn> <expRef> <libType> <umiLen> <identifyTRE> <noNormBw>"
    echo "  rawDataRawDir: raw data directory"
    echo "  sampleInfo: space separated sample information file"
    echo "  runInfo: a description for the run, e.g., 'PROseq_241130' "
    echo "  spikeIn: if spike-in PRO-seq library, 'Y' or 'N' "
    echo "  expRef: experiment genome reference, 'hg38', 'hg19', 'mm10' "
    echo "  libType: Library type, 'qPRO', 'rPRO' or 'PROcap' "
    echo "  umiLen: the length (n nt) of UMI, between 6 and 12"
    echo "  identifyTRE: peak calling method, 'dREG', 'PINTS', 'HOMER', 'all' or 'none' "
    echo "  noNormBw: if get no normalized bigwig (only for low sequencing depth sample to merge), 'Y' or 'N' "

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
identifyTRE=$8
noNormBw=$9

# check args
if [ $# -ne 9 ]; then
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
if [[ "$libType" != "rPRO" && "$libType" != "qPRO" && "$libType" != "PROcap" ]]; then
    echo "Error: Invalid libType value! Must be 'rPRO', 'qPRO' or 'PROcap'."
    usage
fi
if [[ ! ($umiLen -ge 6 && $umiLen -le 12) ]]; then
    echo "Error: Invalid umiLen value! Must between 6 and 12."
    usage
fi
valid_methods=("dREG" "PINTS" "HOMER" "all" "none")
if [[ ! " ${valid_methods[@]} " =~ " $identifyTRE " ]]; then
    echo "Error: Invalid identifyTRE value! Must be 'dREG', 'PINTS', 'HOMER', 'all' or 'none'."
    usage
fi
if [[ "$noNormBw" != "Y" && "$noNormBw" != "N" ]]; then
    echo "Error: Invalid noNormBw value! Must be 'Y' or 'N'."
    usage
fi
echo -e "PRO-seq parameters: \nrawDataRawDir: $rawDataRawDir\nsampleInfo: $sampleInfo\nrunInfo: $runInfo\nspikeIn: $spikeIn\nexpRef: $expRef\nlibType: $libType\numiLen: $umiLen\nidentifyTRE: $identifyTRE\nnoNormBw: $noNormBw"

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
rmdupExpDir=${workDir}/02-3_remove_duplicates_experimental
rmdupSpkDir=${workDir}/02-4_remove_duplicates_spikein
spkScaledBwDir=${workDir}/03_spikein_scaled_bw_track
cpmScaledBwDir=${workDir}/03_cpm_normalized_bw_track
noScaledBwDir=${workDir}/03_no_normalized_bw_track
identifyTREDir=${workDir}/04_identifyTRE
pintsDir=${workDir}/04_identifyTRE/PINTS_TRE
dREGDir=${workDir}/04_identifyTRE/dREG_TRE

# log dirs
logDir=${workDir}/logs
trimLogDir=${logDir}/cutadapt_trim_log
trimFastpLogDir=${logDir}/fastp_label_UMI_log
rmrRNALogDir=${logDir}/bowtie2_remove_rRNA_log
map2ExpLogDir=${logDir}/bowtie2_experimental_genome_alignment_log
map2SpkLogDir=${logDir}/bowtie2_spikein_genome_alignment_log
rmdupExpLogDir=${logDir}/umitools_rmdup_experimental_log
rmdupSpkLogDir=${logDir}/umitools_rmdup_spikein_log
summaryDir=${logDir}/alignment_summary_${runInfo}.txt
spkScaledBwLogDir=${logDir}/bamCoverage_spikein_scaled_bw_log
cpmScaledBwLogDir=${logDir}/bamCoverage_cpm_normalized_bw_log
noScaledBwLogDir=${logDir}/bamCoverage_no_normalized_bw_log

# QC dirs
insertLenDir=${logDir}/flash_estimate_insert_length

# call TRE log
pintsLogDir=${logDir}/PINTS_call_TRE_log
dREGLogDir=${logDir}/dREG_call_TRE_log

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
# # dREG version Feb 17, 2024 (https://github.com/Danko-Lab/dREG/commit/4a1643a69c36b64ce2b8ba7fb76701e744692e4a)
dREG='/share/home/Grape/software/install_pkg/dREG-master/run_dREG.bsh'
# dREG SVR model: (ftp://cbsuftp.tc.cornell.edu/danko/hub/dreg.models/asvm.gdm.6.6M.20170828.rdata)
dREG_model='/share/home/Grape/software/install_pkg/dREG-master/model/asvm.gdm.6.6M.20170828.rdata'
# # R packages in rDeps.R and Rgtsvm(https://github.com/Danko-Lab/Rgtsvm) to accelerate SVM with GPU, please please install it first!
# # bedtools version 2.31.1 (dREG)
# # bedops version 2.4.41 (dREG) (bedops; bedmap; sort-bed) (https://bedops.readthedocs.io/en/latest/index.html)
# # tabix (dREG)
# # bedGraphToBigWig version 2.10 (dREG) (https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig)

# # PINTS verson 1.1.14
pints_caller='/share/home/Grape/software/install_pkg/miniconda3/envs/proseq/bin/pints_caller'
# # HOMER

# # flash version 1.2.11
flash='/share/home/Grape/software/install_pkg/FLASH-1.2.11/flash'

# # R script
Rscript='/opt/ohpc/pub/apps/R/4.2.2/bin/Rscript'
get_degradation_ratio='/share/home/Grape/SEQPIPE/PROSEQ/get_degradation_ratio.R'


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
if [[ ! -d $logDir ]]; then
    mkdir $logDir
fi

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
        echo -e "${arr[@]} $suffix" >> ${logDir}/sampleInfo_new_${runInfo}.txt
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
if [[ ! -s ${rawQcDir}/rawdata_multiqc_${runInfo}.html ]] || [[ ! -d ${rawQcDir}/rawdata_multiqc_${runInfo}_data ]]; then
    $multiqc -f -n rawdata_multiqc_${runInfo} -o $rawQcDir $rawQcDir
fi

# Step1.2: trim adapter and low quality sequence (cutadapt)
# Use cutadapt but not trim_galore to clip more adapter

if [[ ! -d $trimDir ]]; then
    mkdir -p $trimDir
    mkdir -p $trimLogDir
fi
echo -e "\n***************************\nStart run cutadapt at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
for r1 in `ls ${rawDataDir}/*_R1.fq.gz`; do
    r2=${r1/R1.fq.gz/R2.fq.gz}
    sampleName=$(basename ${r1%_R1.fq.gz})
    if [ ! -s ${trimDir}/${sampleName}_trimmed_R2.fq.gz ]; then
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
    fi
    echo -e "Finish cutadapt trim for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
done
echo -e "\nFinish run cutadapt for all samples!"


# QC 1: estimate insert size (FLASH)

if [[ $libType != 'PROcap' ]]; then
    if [[ ! -d $insertLenDir ]]; then
        mkdir -p $insertLenDir
    fi
    echo -e "\n***************************\nStart estimate insert size at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    for r1 in `ls ${trimDir}/*_trimmed_R1.fq.gz`; do
        r2=${r1/R1.fq.gz/R2.fq.gz}
        sampleName=$(basename ${r1%_trimmed_R1.fq.gz})
        if [ ! -s ${insertLenDir}/${sampleName}.hist ]; then
            # eastimate insert size by flash
            $flash --quiet --thread 36 --max-overlap 150 --compress-prog=pigz --compress-prog-args='--processes 36' --suffix=gz -d $insertLenDir -o $sampleName $r1 $r2
        fi
        echo -e "Finish estimate insert for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
    done
    # plot insert size distribution and degradation ratio
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
        --thread 16 \
        --correction \
        --overlap_len_require 15 \
        $extra_params 2> ${log_file_prefix}_fastp.log

}


# qPRO have 6 nt R1 and R2 duplex UMI
# rPRO only have R2 n nt UMI+C on the head of R2 (single UMI)
# sequencing read-through always happen on single UMI rPRO library, clip R1 reverse complementary R2 UMI+G (n+1 nt)

if [[ ! -d $trimFastpDir ]]; then
    mkdir -p ${trimFastpDir}
    mkdir -p ${trimFastpLogDir}
fi
echo -e "\n***************************\nLabeling UMI at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
for r1 in `ls ${trimDir}/*_trimmed_R1.fq.gz`; do
    r2=${r1/R1.fq.gz/R2.fq.gz}
    sampleName=$(basename ${r1%_trimmed_R1.fq.gz})
    outPrefix=${trimFastpDir}/${sampleName}
    logName=${trimFastpLogDir}/${sampleName}
    if [ ! -s ${trimFastpDir}/${sampleName}_fastp_R2.fq.gz ]; then
        if [[ $libType == 'qPRO' || $libType == 'PROcap' ]]; then
            # always 6nt duplex UMI
            fastp_label_umi $r1 $r2 $outPrefix $ADAPTOR_R1 $ADAPTOR_R2 "per_read" $logName ""
        elif [[ $libType == 'rPRO' ]]; then
            # clip R1 reverse complementary R2 UMI+G (n+1 nt)
            fastp_label_umi $r1 $r2 $outPrefix $ADAPTOR_R1 $ADAPTOR_R2 "read2" $logName "--trim_tail1=$(echo "${umiLen}+1" | bc -l) --trim_tail2=0"
        else
            exit 1
        fi
    fi
    echo -e "Finish label UMI for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
done


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
fi
echo -e "\n***************************\nRemoving rRNA reads at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
for r1 in `ls ${trimFastpDir}/*_fastp_R1.fq.gz`; do
    r2=${r1/fastp_R1.fq.gz/fastp_R2.fq.gz}
    sampleName=$(basename ${r1/_fastp_R1.fq.gz/})
    logName=${rmrRNALogDir}/${sampleName}_rmrRNA
    if [ ! -s ${rmrRNALogDir}/${sampleName}_rmrRNA_bowtie2.log ]; then
        # we are interested primarily in quickly identifying and removing any reads that have a valid alignment to the serial alignment genome (-k 1)
        run_bowtie2_alignment $rDNABowtie2Index $r1 $r2 "--fast-local -k 1 --un-conc-gz ${rmrRNADir}/${sampleName}_rmrRNA_R%.fq.gz" $logName > /dev/null
    fi
    echo -e "Finish remove rRNA reads for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
done

# Step1.5: clean data fastqc

if [[ ! -d $cleanQcDir ]]; then
    mkdir -p $cleanQcDir
    echo -e "\n***************************\nSummary QC at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    $fastqc -t 48 --memory 1024 ${rmrRNADir}/*.fq.gz -o $cleanQcDir &> /dev/null
fi
if [[ ! -s ${cleanQcDir}/cleanData_multiqc_${runInfo}.html ]] || [[ ! -d ${cleanQcDir}/cleanData_multiqc_${runInfo}_data ]]; then
    $multiqc -f -n cleanData_multiqc_${runInfo} -o $cleanQcDir $trimLogDir $trimFastpLogDir $rmrRNALogDir $cleanQcDir 
fi

# Step2: Alignment
# Step2.1: align to experimental
# NOTE: To rescue more degradation inset (<20bp, highly degraded libarary like tumor sample...) we set a more sensitive/loose parameters, like change bowtie2 seed length and min score function, although it could cause more multimapper
# in old pipeline, <22bp reads can not be aligned, cause too short seed length and too strict score threshold
# peak calling methods like dREG needs enriched reads to detect the transcriptional peaks, we need more reads: https://github.com/Danko-Lab/proseq2.0/tree/master?tab=readme-ov-file#notes-for-dreg-users

echo -e "\n***************************\nAligning to experimental genome at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
if [[ ! -d $map2ExpDir ]]; then
    mkdir -p $map2ExpDir
    mkdir -p $map2ExpLogDir
fi
# if [[ ! `ls ${map2ExpDir}/*_${exp_info}.bam 2> /dev/null` ]]; then
#     echo -e ""
# fi
for r1 in `ls ${rmrRNADir}/*_rmrRNA_R1.fq.gz`;do
    r2=${r1/R1.fq.gz/R2.fq.gz}
    sampleName=$(basename ${r1%_rmrRNA_R1.fq.gz})
    logName=${map2ExpLogDir}/${sampleName}_${exp_info}
    if [ ! -s ${map2ExpLogDir}/${sampleName}_${exp_info}_bowtie2.log ]; then
        if [[ $libType == 'qPRO' || $libType == 'PROcap' ]]; then
            # remove fastp separator
            run_bowtie2_alignment $expBowtie2Index $r1 $r2 "--very-sensitive-local -X 750 -L 18 --score-min G,20,6" $logName | \
            $samtools view -@ 48 -Sh -q ${MAPQ} -f 2 | \
            sed 's/:UMI_\([A-Z]*\)_\([A-Z]*\)/:UMI_\1\2/g' | \
            $samtools sort -@ 48 -o ${map2ExpDir}/${sampleName}_${exp_info}.bam 2> /dev/null
        elif [[ $libType == 'rPRO' ]]; then
            run_bowtie2_alignment $expBowtie2Index $r1 $r2 "--very-sensitive-local -X 750 -L 18 --score-min G,20,6" $logName | \
            $samtools view -@ 48 -bS -q ${MAPQ} -f 2 | \
            $samtools sort -@ 48 -o ${map2ExpDir}/${sampleName}_${exp_info}.bam 2> /dev/null
        else
            exit 1
        fi
        $samtools index -@ 48 ${map2ExpDir}/${sampleName}_${exp_info}.bam
    fi
    echo -e "Finish align to ${exp_info} for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
done

# Step2.2: align to spike-in

if [[ $spikeIn == 'Y' ]];then
    echo -e "\n***************************\nAligning to spike-in genome at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    if [[ ! -d $map2SpkDir ]];then
        mkdir -p $map2SpkDir
        mkdir -p $map2SpkLogDir
    fi
    # if [[ ! `ls ${map2SpkDir}/*_${spike_info}.bam 2> /dev/null` ]]; then
    #     echo -e ""
    # fi
    for r1 in `ls ${rmrRNADir}/*_rmrRNA_R1.fq.gz`;do
        r2=${r1/R1.fq.gz/R2.fq.gz}
        sampleName=$(basename ${r1%_rmrRNA_R1.fq.gz})
        logName=${map2SpkLogDir}/${sampleName}_${spike_info}
        if [ ! -s ${map2SpkLogDir}/${sampleName}_${spike_info}_bowtie2.log ]; then
            if [[ $libType == 'qPRO' || $libType == 'PROcap' ]]; then
                run_bowtie2_alignment $spkBowtie2Index $r1 $r2 "--very-sensitive-local -X 750 -L 18 --score-min G,20,6" $logName | \
                $samtools view -@ 48 -Sh -q ${MAPQ} -f 2 | \
                sed 's/:UMI_\([A-Z]*\)_\([A-Z]*\)/:UMI_\1\2/g' | \
                $samtools sort -@ 48 -o ${map2SpkDir}/${sampleName}_${spike_info}.bam 2> /dev/null
            elif [[ $libType == 'rPRO' ]]; then
                run_bowtie2_alignment $spkBowtie2Index $r1 $r2 "--very-sensitive-local -X 750 -L 18 --score-min G,20,6" $logName | \
                $samtools view -@ 48 -bS -q ${MAPQ} -f 2 | \
                $samtools sort -@ 48 -o ${map2SpkDir}/${sampleName}_${spike_info}.bam 2> /dev/null
            else
                exit 1
            fi
            $samtools index -@ 48 ${map2SpkDir}/${sampleName}_${spike_info}.bam
        fi
        echo -e "Finish align to ${spike_info} for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
    done
fi

##
## Ongoing: add combined index (experimental + spike-in) and then align to save time for small genome (like dm6)
##

# Step2.3 remove duplicates for experimental bam

echo -e "\n***************************\nDeduplicating for experimental bam at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
if [[ ! -d $rmdupExpDir ]]; then
    mkdir -p $rmdupExpDir
    mkdir -p $rmdupExpLogDir
fi
for bam in `ls ${map2ExpDir}/*_${exp_info}.bam`; do
    fileName=$(basename ${bam%.bam})
    logName=${rmdupExpLogDir}/${fileName}_umitools.log
    rmdupBam=${rmdupExpDir}/${fileName}_rmdup.bam
    rmdupBamStat=${rmdupExpDir}/${fileName}_rmdup.flagstat
    if [ ! -s $rmdupBam ]; then
        $umi_tools dedup --paired --umi-separator=":UMI_" -I $bam -S $rmdupBam &> $logName
        $samtools index -@ 24 $rmdupBam
        $samtools flagstat -@ 24 $rmdupBam > $rmdupBamStat
    fi
    echo -e "Finish deduplicate for ${fileName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
done

# Step2.4 remove duplicates for spike-in bam

if [[ $spikeIn == 'Y' ]]; then
    echo -e "\n***************************\nDeduplicating for spike-in bam at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    if [[ ! -d $rmdupSpkDir ]]; then
        mkdir -p $rmdupSpkDir
        mkdir -p $rmdupSpkLogDir
    fi
    for bam in `ls ${map2SpkDir}/*_${spike_info}.bam`; do
        fileName=$(basename ${bam%.bam})
        logName=${rmdupSpkLogDir}/${fileName}_umitools.log
        rmdupBam=${rmdupSpkDir}/${fileName}_rmdup.bam
        rmdupBamStat=${rmdupSpkDir}/${fileName}_rmdup.flagstat
        if [ ! -s $rmdupBam ]; then
            $umi_tools dedup --paired --umi-separator=":UMI_" -I $bam -S $rmdupBam &> $logName
            $samtools index -@ 24 $rmdupBam
            $samtools flagstat -@ 24 $rmdupBam > $rmdupBamStat
        fi
        echo -e "Finish deduplicate for ${fileName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
    done
fi

# Step3: calculate spike-in scale factor
# For PRO-seq, scale factor is CPM of spike-in mapped reads

if [[ ! -d $summaryDir ]]; then
    echo -e "\n***************************\nCalculating alignment summary at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    echo -e "sample_name\tall_reads\trRNA_reads\trRNA_mapping_ratio\tclean_reads\t${exp_info}_reads\t${exp_info}_mapping_ratio\t${exp_info}_qc_reads\t${exp_info}_qc_ratio\t${spike_info}_qc_reads\t${spike_info}_qc_ratio\t${spike_info}_unique_reads\t${spike_info}_unique_ratio\tscale_qc_factor\tscale_unique_factor" > $summaryDir
    basename -a $(ls ${map2ExpLogDir}/*_${exp_info}_bowtie2.log) | while read file; do
        sampleName=${file/_${exp_info}_bowtie2.log/}
        # allReads: after QC reads
        riboLog=${rmrRNALogDir}/${sampleName}_rmrRNA_bowtie2.log
        allReads=`grep 'reads; of these' $riboLog | awk '{print $1*2}'`
        riboReads=`cat $riboLog | perl -ne 'print "$1\t" while /(\d+) \(/g' | awk '{print 2*($3+$4+$5)+$7+$8}'`
        riboMapRatio=`grep 'overall alignment rate' $riboLog | awk -F '%' '{print $1}'`
        # cleanReads: after remove rRNA reads
        # # exp bowtie2 map reads (aligned concordantly 1 times + aligned concordantly >1 times + aligned discordantly 1 time + singletons)
        expLog=${map2ExpLogDir}/${file}
        cleanReads=`grep 'reads; of these' $expLog | awk '{print $1*2}'`
    	expReads=`cat $expLog | perl -ne 'print "$1\t" while /(\d+) \(/g' | awk '{print 2*($3+$4+$5)+$7+$8}'`
		expMapRatio=`grep 'overall alignment rate' $expLog | awk -F '%' '{print $1}'`
        # remove singleton and >MAPQ reads
        # expQcReads and spkQcRatio: divide by cleanReads
        expStat=${rmdupExpDir}/${sampleName}_${exp_info}_rmdup.flagstat
        expQcReads=`echo "2*$(cat $expStat | grep "read1" | cut -d " " -f 1)" | bc -l`
        expQcRatio=`printf "%.2f\n" $(echo "100*${expQcReads}/${cleanReads}" | bc -l)`
        if [[ $spikeIn == 'Y' ]]; then
            spkStat=${rmdupSpkDir}/${sampleName}_${spike_info}_rmdup.flagstat
            spkQcReads=`echo "2*$(cat $spkStat | grep "read1" | cut -d " " -f 1)" | bc -l`
            spkQcRatio=`printf "%.2f\n" $(echo "100*${spkQcReads}/${cleanReads}" | bc -l)`
            # Unique map to spkin but not map to exp (to avoid homogenous)
            spkUniqueReads=`echo "2*$(comm --check-order -23 <(samtools view -@ 12 ${rmdupSpkDir}/${sampleName}_${spike_info}_rmdup.bam | cut -f 1 | sort -S48G --parallel=24 | uniq) \
                <(samtools view -@ 12 ${rmdupExpDir}/${sampleName}_${exp_info}_rmdup.bam | cut -f 1 | sort -S62G --parallel=36 | uniq) | wc -l)" | bc -l`
            spkUniqueRatio=`printf "%.2f\n" $(echo "100*${spkUniqueReads}/${cleanReads}" | bc -l)`
            scaleQcFactor=`echo "1000000/${spkQcReads}" | bc -l`
            scaleUniqueFactor=`echo "1000000/${spkUniqueReads}" | bc -l`
        fi
        echo -e ${sampleName}"\t"${allReads}"\t"${riboReads}"\t"${riboMapRatio}"\t"${cleanReads}"\t"${expReads}"\t"${expMapRatio}"\t"${expQcReads}"\t"${expQcRatio}"\t"${spkQcReads}"\t"${spkQcRatio}"\t"${spkUniqueReads}"\t"${spkUniqueRatio}"\t"${scaleQcFactor}"\t"${scaleUniqueFactor} >> $summaryDir
    done
    echo -e "Finish calculate alignment summary at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
fi

# Step4: get bw track

run_bamCoverage() {

    local input=$1                # --bam
    local output=$2               # --outFileName
    local extract_flag=$3         # --samFlagInclude
    local norm=$4                 # --normalizeUsing
    local scale_factor=$5         # --scaleFactor
    local extra_params=$6
    local log_file_prefix=$7

    $bamCoverage \
        --skipNonCoveredRegions \
        --binSize 1 \
        --numberOfProcessors 48 \
        --bam $input --outFileName $output \
        --samFlagInclude $extract_flag \
        --normalizeUsing $norm \
        --scaleFactor $scale_factor \
        $extra_params \
        &> ${log_file_prefix}_bamCoverage.log

}

# Step4.1: get track for pol2 position
# Use CPM but not RPKM, cause: https://www.biostars.org/p/9474318/
if [[ $spikeIn == 'Y' ]]; then
    echo -e "\n***************************\nGetting spikein normalized track at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    if [[ ! -d $spkScaledBwDir ]]; then
        mkdir -p $spkScaledBwDir
        mkdir -p $spkScaledBwLogDir
    fi
    cat ${logDir}/alignment_summary_${runInfo}.txt | sed '1d' | while read line; do
        arr=($line)
        sampleName=${arr[0]}
        # if use unique: scaleFactor=${arr[14]}
        scaleFactor=${arr[13]}
        bam=${rmdupExpDir}/${sampleName}_${exp_info}_rmdup.bam
        if [[ $libType != 'PROcap' ]]; then
            sb_fwd=${spkScaledBwDir}/${sampleName}_singlebase_fwd.bw; sb_rev=${spkScaledBwDir}/${sampleName}_singlebase_rev.bw; sb_rev_minus=${spkScaledBwDir}/${sampleName}_singlebase_rev_minus.bw
            fl_fwd=${spkScaledBwDir}/${sampleName}_fulllength_fwd.bw; fl_rev=${spkScaledBwDir}/${sampleName}_fulllength_rev.bw; fl_rev_minus=${spkScaledBwDir}/${sampleName}_fulllength_rev_minus.bw
            if [[ ! -s $sb_fwd || ! -s $fl_fwd ]]; then
                echo "Generating file: $sb_fwd and $fl_fwd"
                run_bamCoverage $bam $sb_fwd 83 "None" $scaleFactor "--Offset 1" ${spkScaledBwLogDir}/${sampleName}_singlebase_fwd
                run_bamCoverage $bam $fl_fwd 83 "None" $scaleFactor "" ${spkScaledBwLogDir}/${sampleName}_fulllength_fwd
            fi
            if [[ ! -s $sb_rev || ! -s $fl_rev ]]; then
                echo "Generating file: $sb_rev and $fl_rev"
                run_bamCoverage $bam $sb_rev 99 "None" $scaleFactor "--Offset 1" ${spkScaledBwLogDir}/${sampleName}_singlebase_rev
                run_bamCoverage $bam $fl_rev 99 "None" $scaleFactor "" ${spkScaledBwLogDir}/${sampleName}_fulllength_rev
            fi
            if [[ ! -s $sb_rev_minus || ! -s $fl_rev_minus ]]; then
                echo "Generating file: $sb_rev_minus and $fl_rev_minus"
                run_bamCoverage $bam $sb_rev_minus 99 "None" "-$scaleFactor" "--Offset 1" ${spkScaledBwLogDir}/${sampleName}_singlebase_rev_minus
                run_bamCoverage $bam $fl_rev_minus 99 "None" "-$scaleFactor" "" ${spkScaledBwLogDir}/${sampleName}_fulllength_rev_minus
            fi
        else
            sb_fwd=${spkScaledBwDir}/${sampleName}_singlebase_fwd.bw; sb_rev=${spkScaledBwDir}/${sampleName}_singlebase_rev.bw; sb_rev_minus=${spkScaledBwDir}/${sampleName}_singlebase_rev_minus.bw
            if [[ ! -s $sb_fwd || ! -s $sb_rev || ! -s $sb_rev_minus ]]; then
                echo "Generating file: $sb_fwd, $sb_rev and $sb_rev_minus"
                run_bamCoverage $bam $sb_fwd 163 "None" $scaleFactor "--Offset 1" ${spkScaledBwLogDir}/${sampleName}_singlebase_fwd
                run_bamCoverage $bam $sb_minus 147 "None" $scaleFactor "--Offset 1" ${spkScaledBwLogDir}/${sampleName}_singlebase_rev
                run_bamCoverage $bam $sb_rev_minus 147 "None" "-$scaleFactor" "--Offset 1" ${spkScaledBwLogDir}/${sampleName}_singlebase_rev_minus
            fi
        fi
        echo -e "Finish get track for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
    done
else
    echo -e "\n***************************\nGetting CPM normalized track at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    if [[ ! -d $cpmScaledBwDir ]]; then
        mkdir -p $cpmScaledBwDir
        mkdir -p $cpmScaledBwLogDir
    fi
    ls ${rmdupExpDir}/*_${exp_info}_rmdup.bam | while read sample; do
        sampleName=$(basename ${sample%_${exp_info}_rmdup.bam})
        if [[ $libType != 'PROcap' ]]; then
            sb_fwd=${cpmScaledBwDir}/${sampleName}_singlebase_fwd.bw; sb_rev=${cpmScaledBwDir}/${sampleName}_singlebase_rev.bw; sb_rev_minus=${cpmScaledBwDir}/${sampleName}_singlebase_rev_minus.bw
            fl_fwd=${cpmScaledBwDir}/${sampleName}_fulllength_fwd.bw; fl_rev=${cpmScaledBwDir}/${sampleName}_fulllength_rev.bw; fl_rev_minus=${cpmScaledBwDir}/${sampleName}_fulllength_rev_minus.bw
            if [[ ! -s $sb_fwd || ! -s $fl_fwd ]]; then
                echo "Generating file: $sb_fwd and $fl_fwd"
                run_bamCoverage $sample $sb_fwd 83 "CPM" 1 "--Offset 1" ${cpmScaledBwLogDir}/${sampleName}_singlebase_fwd
                run_bamCoverage $sample $fl_fwd 83 "CPM" 1 "" ${cpmScaledBwLogDir}/${sampleName}_fulllength_fwd
            fi
            if [[ ! -s $sb_rev || ! -s $fl_rev ]]; then
                echo "Generating file: $sb_rev and $fl_rev"
                run_bamCoverage $sample $sb_rev 99 "CPM" 1 "--Offset 1" ${cpmScaledBwLogDir}/${sampleName}_singlebase_rev
                run_bamCoverage $sample $fl_rev 99 "CPM" 1 "" ${cpmScaledBwLogDir}/${sampleName}_fulllength_rev
            fi
            if [[ ! -s $sb_rev_minus || ! -s $fl_rev_minus ]]; then
                echo "Generating file: $sb_rev_minus and $fl_rev_minus"
                run_bamCoverage $sample $sb_rev_minus 99 "CPM" "-1" "--Offset 1" ${cpmScaledBwLogDir}/${sampleName}_singlebase_rev_minus
                run_bamCoverage $sample $fl_rev_minus 99 "CPM" "-1" "" ${cpmScaledBwLogDir}/${sampleName}_fulllength_rev_minus
            fi
        else
            sb_fwd=${cpmScaledBwDir}/${sampleName}_singlebase_fwd.bw; sb_rev=${cpmScaledBwDir}/${sampleName}_singlebase_rev.bw; sb_rev_minus=${cpmScaledBwDir}/${sampleName}_singlebase_rev_minus.bw
            if [[ ! -s $sb_fwd || ! -s $sb_rev || ! -s $sb_rev_minus ]]; then
                echo "Generating file: $sb_fwd, $sb_rev and $sb_rev_minus"
                run_bamCoverage $sample $sb_fwd 163 "CPM" 1 "--Offset 1" ${cpmScaledBwLogDir}/${sampleName}_singlebase_fwd
                run_bamCoverage $sample $sb_minus 147 "CPM" 1 "--Offset 1" ${cpmScaledBwLogDir}/${sampleName}_singlebase_rev
                run_bamCoverage $sample $sb_rev_minus 147 "CPM" "-1" "--Offset 1" ${cpmScaledBwLogDir}/${sampleName}_singlebase_rev_minus
            fi
        fi
        echo -e "Finish get track for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
    done
fi

# Step4.2: get track for dREG input

if [[ $identifyTRE == 'all' || $identifyTRE == 'dREG' || $noNormBw == 'Y' ]]; then
    echo -e "\n***************************\nGetting no normalized track at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    if [[ ! -d $noScaledBwDir ]]; then
    	mkdir -p $noScaledBwDir
        mkdir -p $noScaledBwLogDir
    fi
    ls ${rmdupExpDir}/*_${exp_info}_rmdup.bam | while read sample; do
        sampleName=$(basename ${sample%_${exp_info}_rmdup.bam})
        sb_fwd=${noScaledBwDir}/${sampleName}_none_normalized_singlebase_fwd.bw; sb_rev_minus=${noScaledBwDir}/${sampleName}_none_normalized_singlebase_rev_minus.bw
        if [[ ! -s $sb_fwd ]]; then
            echo "Generating file: $sb_fwd"
            run_bamCoverage $sample $sb_fwd 83 "None" 1 "--Offset 1" ${noScaledBwLogDir}/${sampleName}_none_normalized_singlebase_fwd
            # if get TSS position
            # run_bamCoverage $sample $sb_fwd 163 "None" 1 "--Offset 1" ${noScaledBwLogDir}/${sampleName}_none_normalized_singlebase_fwd
        fi
        if [[ ! -s $sb_rev_minus ]]; then
            echo "Generating file: $sb_rev_minus"
            run_bamCoverage $sample $sb_rev_minus 99 "None" "-1" "--Offset 1" ${noScaledBwLogDir}/${sampleName}_none_normalized_singlebase_rev_minus
            # if get TSS position
            # run_bamCoverage $sample $sb_rev_minus 147 "None" "-1" "--Offset 1" ${noScaledBwLogDir}/${sampleName}_none_normalized_singlebase_rev_minus
        fi
        echo -e "Finish get track for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
    done
fi

# Step5: identify TRE

# Step5.1 PINTS
pints_call_peak() {

    local bam=$1               # --bam-file
    local out_dir=$2           # --save-to
    local file_prefix=$3       # --file-prefix
    local log_file_prefix=$4

    # PINTS need TSS position
    $pints_caller \
        --exp-type R2_5 \
        --thread 48 \
        --mapq-threshold $MAPQ \
        --bam-file $bam \
        --save-to $out_dir \
        --file-prefix $file_prefix \
        &> ${log_file_prefix}_PINTS.log
}

# Step5.2 dREG

dREG_call_peak() {

    local fwd_bw=$1               # plus_strand.bw
    local rev_minus_bw=$2         # minus_strand.bw (rev_minus)
    local out_prefix=$3
    local log_file_prefix=$4
    
    # bash run_dREG.bsh plus_strand.bw minus_strand.bw out_prefix asvm.RData CPU_cores GPU_id
    $dREG $fwd_bw $rev_minus_bw $out_prefix $dREG_model 28 0 &> ${log_file_prefix}_dREG.log
}


if [[ $identifyTRE != 'none' ]]; then
    echo -e "\n***************************\nIdentifying TRE at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    if [[ ! -d $identifyTREDir ]]; then
        mkdir -p $identifyTREDir
    fi
    # PINTS
    if [[ $identifyTRE == 'all' || $identifyTRE == 'PINTS' ]]; then
        echo -e "\n***************************\nCalling TRE with PINTS at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
        echo -e "PINTS can take multiple bams as input and gain some insight from the replicates. If you have biological replicates, you may need to merge them to call peak: --bam-file test_rep1.bam test_rep2.bam test_rep3.bam\n\n"
        if [[ ! -d $pintsDir ]]; then
            mkdir -p $pintsDir
            mkdir -p $pintsLogDir
        fi
        ls ${rmdupExpDir}/*_${exp_info}_rmdup.bam | while read bam; do
            sampleName=$(basename ${bam%_${exp_info}_rmdup.bam})
            sampleOutDir=${pintsDir}/${sampleName}
            sampleLogPrefix=${pintsLogDir}/${sampleName}
            if [[ ! -d $sampleOutDir ]]; then
                mkdir -p $sampleOutDir
            fi
            if [[ ! -s ${sampleOutDir}/${sampleName}_1_divergent_peaks.bed ]]; then
                pints_call_peak $bam $sampleOutDir $sampleName $sampleLogPrefix
            fi
            echo -e "Finish run pints_caller for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
        done
    fi
    # dREG
    if [[ $identifyTRE == 'all' || $identifyTRE == 'dREG' ]]; then
        echo -e "\n***************************\nCalling TRE with dREG at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
        echo -e "To a reasonable statistical power for discovering TREs, dREG need as few as ~40M uniquely mappable reads. If the sequencing depth is not enough and you have biological replicates, you may need to merge them to call peak\n\n"
        if [[ ! -d $dREGDir ]]; then
            mkdir -p $dREGDir
            mkdir -p $dREGLogDir
        fi
        ls ${noScaledBwDir}/*_none_normalized_singlebase_fwd.bw | while read plus_bw; do
            sampleName=$(basename ${plus_bw%_none_normalized_singlebase_fwd.bw})
            rev_minus_bw=${plus_bw/fwd.bw/rev_minus.bw}
            sampleOutDir=${dREGDir}/${sampleName}
            sampleLogPrefix=${dREGLogDir}/${sampleName}
            if [[ ! -d $sampleOutDir ]]; then
                mkdir -p $sampleOutDir
            fi
            if [[ ! -s ${sampleOutDir}/${sampleName}.dREG.peak.score.bed.gz ]]; then
                ssh -Tqn clg005 "$(typeset -f dREG_call_peak); dREG='$dREG'; dREG_model='$dREG_model'; dREG_call_peak '$plus_bw' '$rev_minus_bw' '${sampleOutDir}/${sampleName}' '$sampleLogPrefix'" &> /dev/null
            fi
            echo -e "Finish run dREG for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
        done
    fi
fi

echo -e "\n***************************\nFinish PRO-seq processing at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"