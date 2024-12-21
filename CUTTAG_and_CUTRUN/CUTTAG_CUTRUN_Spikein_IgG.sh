#!/usr/bin/bash

# *******************************************
# Date: 202412 Wangshaoxuan (shaoxuanwang@hotmail.com)
# Description: Pipeline for CUT&TAG and CUT&RUN
# 1. Spike-in: With or without spike-in
# 2. Spike-in strategy: DNA or nuclei
# 2. Library: CUT&TAG and CUT&RUN
# 3. Negative/Background control: IgG or not (CUT&RUN and CUT&Tag use IgG but not input, Henikoff 2019)
# 4. Alignment: Bowtie2
# 5. Black/Suspect list: Ogata JD, 2023
# 6. Peak calling: SEACR and MACS3
# 7. QC metrics: ENCODE and Zheng Y, 2020
# *******************************************

# usage
# nohup ./xxx.sh /chenfeilab/Gaux/rawDataBackup/xxx/xxxxxx_PRO-seq sampleInfo.txt 241108_PROSEQ N hg38 rPRO 6 none &> 241105_PROseq.log &
# nohup ./xxx.sh "/chenfeilab/Pomelo/try/SXY/24rawdata/*/*" sampleInfo.txt test_SXY51to60_qPRO Y mm10 qPRO 6 all &> test_SXY51to60_qPRO.log &

# sampleInfo.txt file
# rawName newName controlName
# XXX_001 XXX_CUTTAG_H3K4me3_DLD1_DMSO_rep1 XXX_CUTTAG_IgG_DLD1_DMSO_rep1
# XXX_002 XXX_CUTTAG_H3K27ac_DLD1_DMSO_rep1 XXX_CUTTAG_IgG_DLD1_DMSO_rep1
# XXX_003 XXX_CUTTAG_H3K4me3_DLD1_PNUTS_dTAG_rep1 XXX_CUTTAG_IgG_DLD1_PNUTS_dTAG_rep1
# XXX_004 XXX_CUTTAG_H3K27ac_DLD1_PNUTS_dTAG_rep1 XXX_CUTTAG_IgG_DLD1_PNUTS_dTAG_rep1
# XXX_005 XXX_CUTTAG_IgG_DLD1_DMSO_rep1 IgG
# XXX_006 XXX_CUTTAG_IgG_DLD1_PNUTS_dTAG_rep1 IgG
# XXX_007 XXX_CUTTAG_H3K4me3_293T_DMSO_rep1 none
# XXX_008 XXX_CUTTAG_H3K27ac_293T_DMSO_rep1 none


usage() {

    echo "Usage: $0 <rawDataRawDir> <sampleInfo> <runInfo> <spikeIn> <expRef> <spkRef> <spkStrategy> <libType> <rmDup> <controlIgG> <callPeak> <peakType>"

    echo "  rawDataRawDir: raw data directory. If path contains wildcards, enclose the path in quotes to prevent bash wildcard expansion "
    echo "  sampleInfo: space separated sample information file, must have 2 (rawName newName) columns or 3 (rawName newName controlName) columns (if controlIgG is 'Y') "
    echo "  runInfo: a description for the run, e.g., 'CUTTAG_H3K4me3_241230' "
    echo "  spikeIn: if spike-in CUTTAG library, 'Y' or 'N' "
    echo "  expRef: experiment genome reference, 'hg38', 'hg19', 'mm10' "
    echo "  spkRef: spike-in genome reference, 'dm6', 'k12', 'hg38', 'mm10', 'none' "
    echo "  spkStrategy: spike-in normalization strategy, 'DNA', 'nuclei' or 'none'. DNA refers to a fixed DNA sequence, as included in many kits. (However, if the spike-in is randomly fragmented DNA, please set it as 'nuclei'). "
    echo "  libType: library type, 'CUTTAG' or 'CUTRUN' "
    echo "  rmDup: remove duplicates or just mark them, 'remove' or 'mark' "
    echo "  controlIgG: if IgG as negative control, 'Y' or 'N' "
    echo "  callPeak: peak calling method, 'MACS3', 'SEACR' or 'none'"
    echo "  peakType: peak type, 'narrow' or 'broad' "

    exit 1
}

# args
rawDataRawDir=$1
sampleInfo=$2
runInfo=$3
spikeIn=$4
expRef=$5
spkRef=$6
spkStrategy=$7
libType=$8
rmDup=$9
controlIgG=${10}
callPeak=${11}
peakType=${12}

# check args
if [ $# -ne 12 ]; then
    echo "Error: Invalid number of arguments!"
    usage
fi
if [[ "$rawDataRawDir" == *'*'* ]]; then
    baseDir="$rawDataRawDir"
    while [[ "$baseDir" == *'*'* ]]; do
        baseDir="${baseDir%/*}"
    done
    if [[ -d "$baseDir" && $(find $rawDataRawDir 2>/dev/null | head -n 1) ]]; then
        echo -e "Using $rawDataRawDir as input data.\n"
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
if [[ $spikeIn == 'Y' ]]; then
    if [[ "$spkRef" != "dm6" && "$spkRef" != "k12" && "$spkRef" != "hg38" && "$spkRef" != "mm10" ]]; then
        echo "Error: Invalid spkRef value! Must be 'dm6', 'k12', 'hg38' or 'mm10'."
        usage
    fi
    if [[ "$spkStrategy" != "DNA" && "$spkStrategy" != 'nuclei' ]]; then
        echo "Error: Invalid spkStrategy value! Must be 'DNA' or 'nuclei'."
        usage
    fi
else
    if [[ "$spkRef" != "none" ]]; then
        echo "Error: spkRef must be 'none' if spikeIn is 'N'."
        usage
    fi
    if [[ "$spkStrategy" != "none" ]]; then
        echo "Error: spkStrategy must be 'none' if spikeIn is 'N'."
        usage
    fi
fi
if [[ "$libType" != "CUTTAG" && "$libType" != "CUTRUN" ]]; then
    echo "Error: Invalid libType value! Must be 'CUTTAG' or 'CUTRUN'."
    usage
fi
if [[ "$rmDup" != "remove" && "$rmDup" != "mark" ]]; then
    echo "Error: Invalid rmDup value! Must be 'remove' or 'mark'."
    usage
fi
if [[ "$controlIgG" != "Y" && "$controlIgG" != "N" ]]; then
    echo "Error: Invalid controlIgG value! Must be 'Y' or 'N'."
    usage
fi
if [[ $controlIgG == 'Y' ]]; then
    # check 'IgG' character
    if [[ ! $(cat $sampleInfo | awk '$2 ~ /IgG/') ]]; then
        echo "Error: There is no 'IgG' character in your sampleInfo file, please check if you have the correct IgG name or set controlIgG to 'N'."
        usage
    fi
    # check column
    if awk '{if (NF != 3) {exit 1}}' $sampleInfo; then
        :
    else
        echo "Error: Invaild sampleInfo file format, must have exactly 3 columns when exist IgG control."
        exit 1
    fi
else
    if [[ $(cat $sampleInfo | awk '$2 ~ /IgG/') ]]; then
        echo "Warning: It seems you have IgG in this run, but controlIgG is set to 'N'. If you do not intend to use IgG, you may ignore this message."
    fi
fi
if [[ "$callPeak" != "MACS3" && "$callPeak" != "SEACR" && "$callPeak" != "none" ]]; then
    echo "Error: Invalid callPeak value! Must be 'MACS3', 'SEACR' or 'none'."
    usage
fi
if [[ "$peakType" != "narrow" && "$peakType" != "broad" ]]; then
    echo "Error: Invalid peakType value! Must be 'narrow' or 'broad'."
    usage
fi

# print
echo -e "${libType/CUT/CUT&} parameters: \nrawDataRawDir: $rawDataRawDir\nsampleInfo: $sampleInfo\nrunInfo: $runInfo\nspikeIn: $spikeIn\nexpRef: $expRef\nspkRef: $spkRef\nspkStrategy: $spkStrategy\nlibType: $libType\nrmDup: $rmDup\ncontrolIgG: $controlIgG\ncallPeak: $callPeak\npeakType: $peakType"
# MAPQ
MAPQ=30


# output dirs
workDir=`pwd`
rawDataDir=${workDir}/00_raw_data

rawQcDir=${workDir}/01-1_raw_data_fastqc
trimDir=${workDir}/01-2_trimmed_data
cleanQcDir=${workDir}/01-3_trimmed_data_fastqc

map2ExpDir=${workDir}/02-1_experimental_genome_alignment
map2SpkDir=${workDir}/02-2_spikein_genome_alignment
markDupExpDir=${workDir}/02-3_experimental_mark_duplicates
markDupSpkDir=${workDir}/02-4_spikein_mark_duplicates
filterExpBamDir=${workDir}/02-5_experimental_bam_filter
filterSpkBamDir=${workDir}/02-6_spikein_bam_filter

spkScaledBwDir=${workDir}/03_spikein_scaled_bw_track
cpmScaledBwDir=${workDir}/03_cpm_normalized_bw_track

macs3Dir=${workDir}/04_peak_calling_MACS3
seacrDir=${workDir}/04_peak_calling_SEACR


# log dirs
logDir=${workDir}/logs
trimLogDir=${logDir}/trim_galore_log

map2ExpLogDir=${logDir}/bowtie2_experimental_genome_alignment_log
map2SpkLogDir=${logDir}/bowtie2_spikein_genome_alignment_log
markDupExpLogDir=${logDir}/picard_MarkDuplicates_experimental_log
markDupSpkLogDir=${logDir}/picard_MarkDuplicates_spikein_log

summaryDir=${logDir}/alignment_summary_${runInfo}.txt
advSummaryDir=${logDir}/advanced_summary_${runInfo}.txt

spkScaledBwLogDir=${logDir}/bamCoverage_spikein_scaled_bw_log
cpmScaledBwLogDir=${logDir}/bamCoverage_cpm_normalized_bw_log

# call peak log
macs3LogDir=${logDir}/MACS3_call_peak_log
seacrLogDir=${logDir}/SEACR_call_peak_log

# QC dir
insertSizeDir=${logDir}/insert_size_distribution

# all tools in cutrun env update at Oct 2024
# # fastqc version 0.12.1
fastqc='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/bin/fastqc'
# # multiqc version 1.25.1
multiqc='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/bin/multiqc'
# # trim_galore version 0.6.10 with cutadapt version 4.9
trim_galore='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/bin/trim_galore'
# # cutadapt version 4.9
cutadapt='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/bin/cutadapt'
# # bowtie2 version 2.5.4
bowtie2='/share/home/Grape/software/install_pkg/miniconda3/envs/proseq/bin/bowtie2'
# # samtools version 1.21
samtools='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/download/samtools/bin/samtools'
# # picard version 2.25.0 (use old version, we do not need add read group for MarkDuplicates)
picard='/share/home/Grape/software/install_pkg/miniconda3/envs/mainenv/bin/picard'
# # deepTools version 3.5.5
bamCoverage='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/bin/bamCoverage'
alignmentSieve='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/bin/alignmentSieve'
# # bedtools version 2.31.1
bedtools='/share/home/Grape/software/install_pkg/miniconda3/envs/proseq/bin/bedtools'
# # bedops version 2.4.41
bedops='/share/home/Grape/software/bin/bedops'
# # MACS version 3.0.2
macs3='/share/home/Grape/software/install_pkg/miniconda3/envs/chipseq/bin/macs3'
# # SEACR version 1.4-beta.2 (https://github.com/FredHutch/SEACR/archive/refs/tags/v1.4-beta.2.tar.gz)
# # SEACR requires R and Bedtools to be available in your path
seacr='/share/home/Grape/software/install_pkg/SEACR-1.4-beta.2/SEACR_1.4.sh'

# # R script
Rscript='/opt/ohpc/pub/apps/R/4.2.2/bin/Rscript'
plot_insert_size_distribution='/share/home/Grape/SEQPIPE/CUTTAG_and_CUTRUN/plot_insert_size_distribution.R'

# Bowtie2 index
hg38Index='/share/home/Grape/genome/Homo_sapiens/bowtie2_index/hg38/hg38'
hg19Index='/share/home/Blueberry/reference/index/bowtie2/hg19/hg19'
mm10Index='/share/home/Grape/genome/Mus_musculus/bowtie2_index/mm10/mm10'
dm6Index='/share/home/Grape/genome/Drosophila_melanogaster/bowtie2_index/dm6/dm6'
k12Index='/share/home/Grape/genome/Escherichia_coli/bowtie2_index/ASM584v2/ASM584v2'

# genome size
hg38gSize=2913022398
hg19gSize=2864785220
mm10gSize=2652783500

# black list
# Jonathan 2023: hg38 we use Kundaje 2020, hg19 and mm10 we use Boyle v2
hg38Blacklist='/share/home/Grape/genome/Homo_sapiens/blacklist/ENCFF356LFX.bed'
hg19Blacklist='/share/home/Grape/genome/Homo_sapiens/blacklist/hg19-blacklist.v2.bed'
mm10Blacklist='/share/home/Grape/genome/Mus_musculus/blacklist/mm10-blacklist.v2.bed'
# CUT&RUN blacklist (only for CUT&RUN)
hg38BlacklistCR='/share/home/Grape/genome/Homo_sapiens/blacklist/hg38_Nordin_CUTRUN_blacklist.bed'
mm10BlacklistCR='/share/home/Grape/genome/Mus_musculus/blacklist/mm10_Nordin_CUTRUN_blacklist.bed'

# reference genome
if [[ $(echo $expRef | sed 's/[0-9]//g') == $(echo $spkRef | sed 's/[0-9]//g') ]]; then
    echo "Error: expRef and spkRef must from different organism genomes!"
    exit 1
fi

if [[ $spikeIn == 'Y' ]]; then
    spkIndex="${spkRef}Index"
    spkBowtie2Index=${!spkIndex}
    spike_info=${spkRef}
else
    spike_info=''
fi

expIndex="${expRef}Index"
expBowtie2Index=${!expIndex}
exp_info=${expRef}
gSize="${expRef}gSize"
expGenomeSize=${!gSize}
blacklist="${expRef}Blacklist"
expBlacklist=${!blacklist}
echo -e "$exp_info blacklist is: $expBlacklist"

# Start
echo -e "\n***************************\n${libType/CUT/CUT&} processing at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
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

# Step1.2: trim adapter and low quality sequence (trim_galore)

nohup_number=0
if [[ ! -d $trimDir ]]; then
    mkdir -p $trimDir
    mkdir -p $trimLogDir
    for r1 in `ls ${rawDataDir}/*_R1.fq.gz`; do
        r2=${r1/R1.fq.gz/R2.fq.gz}
        sampleName=$(basename ${r1%_R1.fq.gz})
        $trim_galore --phred33 -j 4 -o $trimDir -q 25 --length 15 -e 0.1 --stringency 4 \
        --path_to_cutadapt $cutadapt \
        --paired $r1 $r2 &> ${trimLogDir}/${sampleName}_cutadapt.log &
        nohup_number=`echo $nohup_number+8 | bc`
        if [[ $nohup_number -eq 48 ]]; then
            wait
            nohup_number=0
        fi
    done
    wait
    echo -e "\n\nFinish trimming at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n"
fi

# Step1.3: trimmed data qc

if [[ ! -d $cleanQcDir ]]; then
    mkdir -p $cleanQcDir
    $fastqc -t 48 --memory 1024 ${trimDir}/*.fq.gz -o $cleanQcDir &> /dev/null
fi
if [[ ! -s ${cleanQcDir}/cleanData_multiqc_${runInfo}.html ]] || [[ ! -d ${cleanQcDir}/cleanData_multiqc_${runInfo}_data ]]; then
    $multiqc -f -n cleanData_multiqc_${runInfo} -o $cleanQcDir $trimLogDir $cleanQcDir 
fi

####
# Step2: Alignment and QC

run_bowtie2_alignment() {

    local genome_dir=$1          # -x
    local read1=$2               # -1 and -2
    local read2=$3
    local extra_params=$4
    local log_file_prefix=$5
    local out_prefix=$6

    echo -e "Run bowtie2 with:\n-x $genome_dir -1 $read1 -2 $read2 --threads 48 $extra_params\n" > ${log_file_prefix}_bowtie2.log

    $bowtie2 -x $genome_dir \
        -1 $read1 -2 $read2 \
        --threads 48 \
        $extra_params \
        2>> ${log_file_prefix}_bowtie2.log | \
        $samtools sort -@ 48 -o ${out_prefix}.bam 2> /dev/null

}

# Step2.1: align to experimental

echo -e "\n***************************\nAligning to experimental genome ($exp_info) at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
if [[ ! -d $map2ExpDir ]]; then
    mkdir -p $map2ExpDir
    mkdir -p $map2ExpLogDir
fi

for r1 in `ls ${trimDir}/*_1.fq.gz`; do
    r2=${r1/R1_val_1.fq.gz/R2_val_2.fq.gz}
    sampleName=$(basename ${r1%_R1_val_1.fq.gz})
    logName=${map2ExpLogDir}/${sampleName}_${exp_info}
    outName=${map2ExpDir}/${sampleName}_${exp_info}
    if [ ! -s ${logName}_bowtie2.log ]; then
        # for histone modification CUTTAG, set -X 1000 (may be set 2000 is also good), do not set --dovetail
        if [[ $libType == 'CUTTAG' ]]; then
            run_bowtie2_alignment $expBowtie2Index $r1 $r2 "--local --very-sensitive-local --no-mixed --no-discordant -I 15 -X 1000" $logName $outName
        # set --dovetail for short CUTRUN fragments (adapters may could not be removed completely..)
        # for TF CUT&TAG -X 700 is enough
        # --soft-clipped-unmapped-tlen for get --dovetail fragment size
        elif [[ $libType == 'CUTRUN' ]]; then
            run_bowtie2_alignment $expBowtie2Index $r1 $r2 "--local --very-sensitive-local --dovetail --soft-clipped-unmapped-tlen --no-mixed --no-discordant -I 15 -X 700" $logName $outName
        else
            exit 1
        fi
        $samtools index -@ 48 ${outName}.bam
    fi
    echo -e "Finish align to ${exp_info} for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
done

# Step2.2: align to spike-in

if [[ $spikeIn == 'Y' ]]; then
    echo -e "\n***************************\nAligning to spike-in genome ($spike_info) at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    if [[ ! -d $map2SpkDir ]]; then
        mkdir -p $map2SpkDir
        mkdir -p $map2SpkLogDir
    fi
    for r1 in `ls ${trimDir}/*_1.fq.gz`; do
        r2=${r1/R1_val_1.fq.gz/R2_val_2.fq.gz}
        sampleName=$(basename ${r1%_R1_val_1.fq.gz})
        logName=${map2SpkLogDir}/${sampleName}_${spike_info}
        outName=${map2SpkDir}/${sampleName}_${spike_info}
        # the alignment options is depends on the normalization strategy of spike-in
        if [ ! -s ${logName}_bowtie2.log ]; then
            if [[ $spkStrategy == 'DNA' ]]; then
                run_bowtie2_alignment $spkBowtie2Index $r1 $r2 "--end-to-end --very-sensitive --no-mixed --no-discordant" $logName $outName
            else
                if [[ $libType == 'CUTTAG' ]]; then
                    run_bowtie2_alignment $spkBowtie2Index $r1 $r2 "--local --very-sensitive-local --no-mixed --no-discordant -I 15 -X 1000" $logName $outName
                elif [[ $libType == 'CUTRUN' ]]; then
                    run_bowtie2_alignment $spkBowtie2Index $r1 $r2 "--local --very-sensitive-local --dovetail --soft-clipped-unmapped-tlen --no-mixed --no-discordant -I 15 -X 700" $logName $outName
                else
                    exit 1
                fi
            fi
            $samtools index -@ 48 ${outName}.bam
        fi
        echo -e "Finish align to ${spike_info} for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
    done
fi

# Step2.3 mark or remove duplicates for experimental bam

if [[ ! -d $markDupExpDir ]]; then
    mkdir -p $markDupExpDir
    mkdir -p $markDupExpLogDir
fi

echo -e "\n***************************\nMarking or removing duplicates for experimental bam ($exp_info) at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

for bam in `ls ${map2ExpDir}/*_${exp_info}.bam`; do
    fileName=$(basename ${bam%.bam})
    markDupBam=${markDupExpDir}/${fileName}.markdup.bam
    markDupMetric=${markDupExpDir}/${fileName}.markdup.metrics
    if [[ ! -s $markDupMetric ]]; then
        $picard MarkDuplicates \
            --ASSUME_SORT_ORDER coordinate \
            --VALIDATION_STRINGENCY LENIENT \
            --REMOVE_DUPLICATES false \
            --INPUT $bam \
            --OUTPUT $markDupBam \
            --METRICS_FILE $markDupMetric &> ${markDupExpLogDir}/${fileName}_markdup.log
    fi
    echo -e "Finish mark exp duplicates for ${fileName%_$exp_info} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"

    # remove IgG duplicates (IgG control samples always have high level duplicates)
    if [[ $controlIgG == 'Y' && $fileName =~ .*IgG.* ]]; then
        rmdupBam=${markDupExpDir}/${fileName}.rmdup.bam
        if [[ ! -s $rmdupBam ]]; then
            $samtools view -@ 48 -bSh -F 1024 $markDupBam > $rmdupBam
            echo -e "Finish remove IgG duplicates for ${fileName%_$exp_info}, removing marked bam..."
            unlink $markDupBam
        fi
    fi

    # if remove exp duplicates
    if [[ $rmDup == 'remove' ]]; then
        rmdupBam=${markDupExpDir}/${fileName}.rmdup.bam
        if [[ ! -s $rmdupBam ]]; then
            $samtools view -@ 48 -bSh -F 1024 $markDupBam > $rmdupBam
            echo -e "Finish remove exp duplicates for ${fileName%_$exp_info}, removing marked bam..."
            unlink $markDupBam
        fi
    fi
done


# Step2.4 remove duplicates for spike-in bam

if [[ $spikeIn == 'Y' && $rmDup == 'remove' ]]; then
    if [[ $spkStrategy == 'DNA' ]]; then
        echo -e "For DNA spike-in, SKIP remove duplicates!\n"
    else
        echo -e "\n***************************\nDeduplicating for spike-in bam ($spike_info) at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
        if [[ ! -d $markDupSpkDir ]]; then
            mkdir -p $markDupSpkDir
            mkdir -p $markDupSpkLogDir
        fi
        # just remove!
        for bam in `ls ${map2SpkDir}/*_${spike_info}.bam`; do
            fileName=$(basename ${bam%.bam})
            rmDupBam=${markDupSpkDir}/${fileName}.rmdup.bam
            rmDupMetric=${markDupSpkDir}/${fileName}.markdup.metrics
            if [[ ! -s $rmDupMetric ]]; then
                $picard MarkDuplicates \
                    --ASSUME_SORT_ORDER coordinate \
                    --VALIDATION_STRINGENCY LENIENT \
                    --REMOVE_DUPLICATES true \
                    --INPUT $bam \
                    --OUTPUT $rmDupBam \
                    --METRICS_FILE $rmDupMetric &> ${markDupSpkLogDir}/${fileName}_markdup.log
            fi
            echo -e "Finish remove spk duplicates for ${fileName%_$spike_info} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
        done
    fi
fi

# Step2.5 filter exp bam and shift reads for C&T

if [[ ! -d $filterExpBamDir ]]; then
    mkdir -p $filterExpBamDir
fi

echo -e "\n***************************\nFilter for experimental bam ($exp_info) at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

ls ${markDupExpDir}/*_${exp_info}.*.bam | while read bam; do 
    fileName=$(basename ${bam%.*.bam})
    cleanBam=${filterExpBamDir}/${fileName}.clean.bam
    cleanBamStat=${filterExpBamDir}/${fileName}.clean.flagstat
    if [[ ! -s $cleanBam ]]; then
        echo -e "Filter ${fileName%_$exp_info}..."
        $samtools view -@ 48 -bS -q ${MAPQ} -f 2 -e 'rname != "chrM" && rname != "chrMT"' $bam | $samtools sort -@ 48 -O BAM -o $cleanBam 2> /dev/null
        $samtools index -@ 48 $cleanBam
        $samtools flagstat -@ 48 $cleanBam > $cleanBamStat
    fi
    # shift
    if [[ $libType == 'CUTTAG' ]]; then
        tempBam=${filterExpBamDir}/${fileName}.tmp.bam
        shiftBam=${filterExpBamDir}/${fileName}.clean.shift.bam
        if [[ ! -s $shiftBam ]]; then
            echo -e "Shifting reads for ${fileName%_$exp_info}..."
            $alignmentSieve --numberOfProcessors 48 --ATACshift --bam $cleanBam -o $tempBam
            $samtools sort -@ 48 $tempBam -o $shiftBam 2> /dev/null
            $samtools index -@ 48 $shiftBam
            unlink $tempBam
        fi
    fi
done

echo -e "Finish all experimental filter\n"

# summary
$multiqc -f -n exp_alignment_multiqc_${runInfo} -o $filterExpBamDir $map2ExpDir $markDupExpDir $filterExpBamDir

# Step2.6 filter spk bam

if [[ $spikeIn == 'Y' ]]; then
    echo -e "\n***************************\nFilter for spike-in bam ($spike_info) at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    if [[ ! -d $filterSpkBamDir ]]; then
        mkdir -p $filterSpkBamDir
    fi
    ls ${map2SpkDir}/*_${spike_info}.bam | while read file; do
        fileName=$(basename ${file%.bam})
        cleanBam=${filterSpkBamDir}/${fileName}.clean.bam
        cleanBamStat=${filterSpkBamDir}/${fileName}.clean.flagstat
        if [[ $spkStrategy == 'DNA' || $rmDup == 'mark' ]]; then
            bam=${file}
        else
            bam=${markDupSpkDir}/$(basename ${file%.bam}).rmdup.bam
        fi
        # filter
        if [[ ! -s $cleanBam ]]; then
            echo -e "Filter $bam..."
            $samtools view -@ 48 -bS -q ${MAPQ} -f 2 -e 'rname != "chrM" && rname != "chrMT"' $bam | $samtools sort -@ 48 -O BAM -o $cleanBam 2> /dev/null
            $samtools index -@ 48 $cleanBam
            $samtools flagstat -@ 48 $cleanBam > $cleanBamStat
        fi
    done
    echo -e "Finish all spike-in filter\n"
    # summary
    $multiqc -f -n spk_alignment_multiqc_${runInfo} -o $filterSpkBamDir $map2SpkDir $filterSpkBamDir
fi


# Step3: calculate spike-in scale factor
# For C&T and C&R, scale factor is "C / fragments mapped to spike-in genome" (the C/const is depend on the percentage of spike-in reads)

if [[ ! -s $summaryDir ]]; then
    echo -e "\n***************************\nCalculating alignment summary at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    echo -e "sample_name\tall_reads\t${exp_info}_reads\t${exp_info}_mapping_ratio\t${exp_info}_qc_reads\t${exp_info}_qc_ratio\t${spike_info}_qc_reads\t${spike_info}_qc_ratio\t${spike_info}_unique_reads\t${spike_info}_unique_ratio\tscale_qc_factor\tscale_unique_factor" > $summaryDir
    basename -a $(ls ${map2ExpLogDir}/*_${exp_info}_bowtie2.log) | while read file; do
        sampleName=${file/_${exp_info}_bowtie2.log/}
        expLog=${map2ExpLogDir}/${file}
        # allReads: after QC (alignment input) reads
        allReads=`grep 'reads; of these' $expLog | awk '{print $1*2}'`
        # exp bowtie2 map reads (aligned concordantly 1 times + aligned concordantly >1 times + aligned discordantly 1 time + singletons)
    	expReads=`cat $expLog | perl -ne 'print "$1\t" while /(\d+) \(/g' | awk '{print 2*($3+$4+$5)+$7+$8}'`
		expMapRatio=`grep 'overall alignment rate' $expLog | awk -F '%' '{print $1}'`
        # filter out mitochondria reads, singleton and >MAPQ reads (and duplicates reads (if rmDup is 'remove'))
        expStat=${filterExpBamDir}/${sampleName}_${exp_info}.clean.flagstat
        expQcReads=`echo "2*$(cat $expStat | grep "read1" | cut -d " " -f 1)" | bc -l`
        expQcRatio=`printf "%.4f\n" $(echo "100*${expQcReads}/${allReads}" | bc -l)`
        if [[ $spikeIn == 'Y' ]]; then
            spkStat=${filterSpkBamDir}/${sampleName}_${spike_info}.clean.flagstat
            spkQcReads=`echo "2*$(cat $spkStat | grep "read1" | cut -d " " -f 1)" | bc -l`
            spkQcRatio=`printf "%.4f\n" $(echo "100*${spkQcReads}/${allReads}" | bc -l)`
            # Unique map to spkin but not map to exp (to avoid homogenous)
            spkUniqueReads=`echo "2*$(comm --check-order -23 <($samtools view -@ 12 ${filterSpkBamDir}/${sampleName}_${spike_info}.clean.bam | cut -f 1 | sort -S48G --parallel=24 | uniq) \
                <($samtools view -@ 12 ${filterExpBamDir}/${sampleName}_${exp_info}.clean.bam | cut -f 1 | sort -S62G --parallel=36 | uniq) | wc -l)" | bc -l`
            spkUniqueRatio=`printf "%.4f\n" $(echo "100*${spkUniqueReads}/${allReads}" | bc -l)`
            if [[ $spkStrategy == 'DNA' ]]; then
                const=10000
            else
                const=1000000
            fi
            scaleQcFactor=`echo "${const}/${spkQcReads}" | bc -l`
            scaleUniqueFactor=`echo "${const}/${spkUniqueReads}" | bc -l`
        fi
        echo -e ${sampleName}"\t"${allReads}"\t"${expReads}"\t"${expMapRatio}"\t"${expQcReads}"\t"${expQcRatio}"\t"${spkQcReads}"\t"${spkQcRatio}"\t"${spkUniqueReads}"\t"${spkUniqueRatio}"\t"${scaleQcFactor}"\t"${scaleUniqueFactor} >> $summaryDir
        if [[ $spikeIn == 'N' ]]; then
            cut -f 1-6 < $summaryDir > temp.txt && mv temp.txt $summaryDir
        fi        
    done
    echo -e "Finish calculate alignment summary at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
fi

# Step3.2: get insert size distribution

if [[ ! -d $insertSizeDir ]]; then
    mkdir -p $insertSizeDir
fi

echo -e "\n***************************\nPlot insert size distribution at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

cat ${logDir}/alignment_summary_${runInfo}.txt | sed '1d' | while read line; do
    arr=($line)
    sampleName=${arr[0]}
    txtFile=${insertSizeDir}/${sampleName}_insert_size.txt
    if [[ ! -s $txtFile ]]; then
        if [[ $libType == 'CUTTAG' ]]; then
            bam=${filterExpBamDir}/${sampleName}_${exp_info}.clean.shift.bam
        else
            bam=${filterExpBamDir}/${sampleName}_${exp_info}.clean.bam
        fi
        $samtools view -@ 48 -f 66 $bam | cut -f 9 | awk '{print sqrt($0^2)}' | sort -S48G --parallel=24 | uniq -c | sort -b -k2,2n | sed -e 's/^[ \\t]*//' > $txtFile
    fi
done

# plot
echo -e "Finish get insert size, plot...\n"
$Rscript $plot_insert_size_distribution $insertSizeDir $runInfo $libType

# Step4: get bw track

run_bamCoverage() {

    local input=$1                # --bam
    local output=$2               # --outFileName
    local black_list=$3           # --blackListFileName
    local norm=$4                 # --normalizeUsing
    local scale_factor=$5         # --scaleFactor
    local log_file_prefix=$6

    $bamCoverage \
        --skipNonCoveredRegions \
        --binSize 1 \
        --numberOfProcessors 48 \
        --bam $input --outFileName $output \
        --blackListFileName $black_list \
        --normalizeUsing $norm \
        --scaleFactor $scale_factor \
        &> ${log_file_prefix}_bamCoverage.log

}


if [[ $spikeIn == 'Y' ]]; then
    echo -e "\n***************************\nGetting spikein normalized track at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    if [[ ! -d $spkScaledBwDir ]]; then
        mkdir -p $spkScaledBwDir
        mkdir -p $spkScaledBwLogDir
    fi
    cat ${logDir}/alignment_summary_${runInfo}.txt | sed '1d' | while read line; do
        arr=($line)
        sampleName=${arr[0]}
        # if use unique: scaleFactor=${arr[11]}
        scaleFactor=${arr[10]}
        bw=${spkScaledBwDir}/${sampleName}.bw
        bwLogPrefix=${spkScaledBwLogDir}/${sampleName}
        if [[ ! -s $bw ]]; then
            if [[ $libType == 'CUTTAG' ]]; then
                bam=${filterExpBamDir}/${sampleName}_${exp_info}.clean.shift.bam
            else
                bam=${filterExpBamDir}/${sampleName}_${exp_info}.clean.bam
            fi
            run_bamCoverage $bam $bw $expBlacklist "None" $scaleFactor $bwLogPrefix
        fi
        echo -e "Finish get track for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
    done
else
    echo -e "\n***************************\nGetting CPM normalized track at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    if [[ ! -d $cpmScaledBwDir ]]; then
        mkdir -p $cpmScaledBwDir
        mkdir -p $cpmScaledBwLogDir
    fi
    cat ${logDir}/alignment_summary_${runInfo}.txt | sed '1d' | while read line; do
        arr=($line)
        sampleName=${arr[0]}
        bw=${cpmScaledBwDir}/${sampleName}.bw
        bwLogPrefix=${cpmScaledBwLogDir}/${sampleName}
        if [[ ! -s $bw ]]; then
            if [[ $libType == 'CUTTAG' ]]; then
                bam=${filterExpBamDir}/${sampleName}_${exp_info}.clean.shift.bam
            else
                bam=${filterExpBamDir}/${sampleName}_${exp_info}.clean.bam
            fi
            run_bamCoverage $bam $bw $expBlacklist "CPM" 1 $bwLogPrefix
        fi
        echo -e "Finish get track for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
    done 
fi


# Step5: peak calling

# Step5.1 MACS3
macs3_call_peak() {

    local bam=$1               # --treatment
    local out_dir=$2           # --outdir
    local file_prefix=$3       # --name
    local extra_params=$4
    local log_file_prefix=$5

    $macs3 callpeak \
        --format BAMPE \
        --qvalue 0.05 \
        --keep-dup all \
        --bdg \
        --treatment $bam \
        --gsize $expGenomeSize \
        --outdir $out_dir \
        --name $file_prefix \
        $extra_params \
        2> ${log_file_prefix}_macs3.log
}

# Step5.2 SEACR
seacr_call_peak() {

    local bdg=$1                  # --bedgraph, -b
    local norm=$2                 # --normalize, -n
    local out_prefix=$3           # --output, -o
    local extra_params=$4
    local log_file_prefix=$5
    
    $seacr \
    --bedgraph $bdg \
    --normalize $norm \
    --mode stringent \
    --output $out_prefix \
    --extension 0.1 \
    $extra_params \
    2> ${log_file_prefix}_SEACR.log
}

# call peak

if [[ $callPeak == 'MACS3' ]]; then
    echo -e "\n***************************\nCalling peak with $callPeak at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    echo -e "For lots of TFs and some histone modifications in a specific cell line or condition, there may not be a best-practice set of parameters. So it may require trying different parameters, you could set callPeak to "none" to skip this step.\n"
    # MACS3
    if [[ ! -d $macs3Dir ]]; then
        mkdir -p $macs3Dir
        mkdir -p $macs3LogDir
    fi
    cat $sampleInfo | while read line; do
        arr=($line)
        sampleName=${arr[1]}
        # call peak with control
        if [[ $controlIgG == 'Y' ]]; then
            # prepare input and log
            controlName=${arr[2]}
            if [[ $controlName != 'none' && $controlName != 'IgG' ]]; then
                peakLogPrefix=${macs3LogDir}/${sampleName}_${peakType}_with_control
                if [[ $libType == 'CUTTAG' ]]; then
                    treatmentBam=${filterExpBamDir}/${sampleName}_${exp_info}.clean.shift.bam
                    controlBam=${filterExpBamDir}/${controlName}_${exp_info}.clean.shift.bam
                else
                    treatmentBam=${filterExpBamDir}/${sampleName}_${exp_info}.clean.bam
                    controlBam=${filterExpBamDir}/${controlName}_${exp_info}.clean.bam
                fi
                echo -e "${sampleName} peak calling with control: ${controlName}..."
            elif [[ $controlName == 'none' ]]; then
                peakLogPrefix=${macs3LogDir}/${sampleName}_${peakType}_without_control
                if [[ $libType == 'CUTTAG' ]]; then
                    treatmentBam=${filterExpBamDir}/${sampleName}_${exp_info}.clean.shift.bam
                else
                    treatmentBam=${filterExpBamDir}/${sampleName}_${exp_info}.clean.bam
                fi
                echo -e "${sampleName} peak calling without control..."
            else
                echo -e "${sampleName} is control or invalid file, skip..."
                break
            fi
            # output
            if [[ $peakType == 'narrow' ]]; then
                peakBed=${macs3Dir}/${sampleName}_peaks.narrowPeak
            else
                peakBed=${macs3Dir}/${sampleName}_peaks.broadPeak
            fi
            # call peak
            if [[ ! -s $peakBed ]] ; then
                if [[ $controlName != 'none' ]]; then
                    if [[ $peakType == 'narrow' ]]; then
                        macs3_call_peak $treatmentBam $macs3Dir $sampleName "--control $controlBam" $peakLogPrefix
                    else
                        macs3_call_peak $treatmentBam $macs3Dir $sampleName "--control $controlBam --broad --broad-cutoff 0.1" $peakLogPrefix
                    fi
                else
                    if [[ $peakType == 'narrow' ]]; then
                        macs3_call_peak $treatmentBam $macs3Dir $sampleName "" $peakLogPrefix
                    else
                        macs3_call_peak $treatmentBam $macs3Dir $sampleName "--nolambda --broad --broad-cutoff 0.1" $peakLogPrefix
                    fi
                fi
            fi
        else
            # input
            echo -e "${sampleName} peak calling without control..."
            peakLogPrefix=${macs3LogDir}/${sampleName}_${peakType}_without_control
            if [[ $libType == 'CUTTAG' ]]; then
                treatmentBam=${filterExpBamDir}/${sampleName}_${exp_info}.clean.shift.bam
            else
                treatmentBam=${filterExpBamDir}/${sampleName}_${exp_info}.clean.bam
            fi
            # output
            if [[ $peakType == 'narrow' ]]; then
                peakBed=${macs3Dir}/${sampleName}_peaks.narrowPeak
            else
                peakBed=${macs3Dir}/${sampleName}_peaks.broadPeak
            fi
            # call peak
            if [[ ! -s $peakBed ]] ; then
                if [[ $peakType == 'narrow' ]]; then
                    macs3_call_peak $treatmentBam $macs3Dir $sampleName "" $peakLogPrefix
                else
                    macs3_call_peak $treatmentBam $macs3Dir $sampleName "--nolambda --broad --broad-cutoff 0.1" $peakLogPrefix
                fi
            fi
        fi
        # filter black list peak
        if [[ $peakType == 'narrow' ]]; then
            filterPeakBed=${macs3Dir}/${sampleName}_filter_peaks.narrowPeak
        else
            filterPeakBed=${macs3Dir}/${sampleName}_filter_peaks.broadPeak
        fi
        if [[ ! -s $filterPeakBed ]] ; then
            bedops -n 1 $peakBed $expBlacklist > $filterPeakBed
        fi
        echo -e "Finish call ${peakType} peak for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
    done
elif [[ $callPeak == 'SEACR' ]]; then
    :
else
    echo -e "\nSkip peak calling\n"
fi


# Step5.3: advanced QC metrics
# Note: --local always has lower "unique mapped ratio" than --end-to-end mode, especially for short reads. It could be less strict than the >50%-70% (some people expected).

if [[ ! -s $advSummaryDir ]]; then
    echo -e "\n***************************\nCalculating advanced summary at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    echo -e "sample_name\t${exp_info}_mapping_ratio\t${exp_info}_unique_mapping_ratio\tmt_mapping_ratio\t${exp_info}_dup_ratio\tpeak_counts\tFRiP" > $advSummaryDir
    cat $summaryDir | sed '1d' | while read line; do
        sampleName=`echo $line | awk '{print $1}'`
        # exp genome mapped: basic metrics column 4
        expMapRatio=`echo $line | awk '{print $4}'`
        # exp genome unique mapped
        bowtie2Log=${map2ExpLogDir}/${sampleName}_${exp_info}_bowtie2.log
        expUniqMapRatio=`cat $bowtie2Log | grep 'aligned concordantly exactly' | perl -ne 'print $1 if /\(([\d\.]+)%\)/'`
        # mt ratio
        expBam=${map2ExpDir}/${sampleName}_${exp_info}.bam
        expReads=`echo $line | awk '{print $3}'`
        mtMapReads=`samtools idxstats -@ 48 $expBam | grep 'chrM' | awk '{sum+=$3} END {print sum}'`
        mtMapRatio=`printf "%.2f\n" $(echo "100*${mtMapReads}/${expReads}" | bc -l)`
        # duplication ratio
        markDupMetric=${markDupExpDir}/${sampleName}_${exp_info}.markdup.metrics
        dupRatio=`printf "%.2f\n" $(echo "100*$(cat $markDupMetric | grep -A 2 '## METRICS CLASS' | awk -F '\t' '{print $9}' | sed '1,2d')" | bc -l)`
        # peak counts
        if [[ $peakType == 'narrow' ]]; then
            peakBed=${macs3Dir}/${sampleName}_peaks.narrowPeak
            nCol=11
        else
            peakBed=${macs3Dir}/${sampleName}_peaks.broadPeak
            nCol=10
        fi
        peakCounts=`wc -l $peakBed`
        # FRiP (fraction of reads in peaks), peak already sorted
        filterExpBam=${filterExpBamDir}/${sampleName}_${exp_info}.clean.bam
        # peakReads=`bedtools intersect -u -ubam -nonamecheck -a $filterExpBam -b $peakBed | samtools view -@ 24 -c`
        peakReads=`bedtools intersect -c -nonamecheck -a $peakBed -b $filterExpBam | awk -v col="$nCol" '{sum+=$col} END {print sum}'`
        expQcReads=`echo $line | awk '{print $5}'`
        FRiP=`printf "%.2f\n" $(echo "100*${peakReads}/${expQcReads}" | bc -l)`
        echo -e ${sampleName}"\t"${expMapRatio}"\t"${expUniqMapRatio}"\t"${mtMapRatio}"\t"${dupRatio}"\t"${peakCounts}"\t"${FRiP} >> $advSummaryDir
    done
    echo -e "Finish calculate advanced metrics summary at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
fi

echo -e "\n***************************\nFinish ${libType/CUT/CUT&} processing at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"