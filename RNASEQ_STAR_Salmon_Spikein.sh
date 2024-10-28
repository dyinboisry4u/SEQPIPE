#!/bin/bash

# *******************************************
# Date: 202410 Wangshaoxuan
# Description: Pipeline for RNA-Seq
# 1. Spike-in: With/or not with spike-in
# 2. Alignment: STAR 
# 3. Quantification: featureCounts or Salmon
# *******************************************

# usage
# nohup ./xxx.sh /chenfeilab/Gaux/rawDataBackup/xxx/xxxxxx_RNA-Seq ./sampleInfo.txt N hg38 Salmon &> ./241018_RNASeq.log &

# sampleInfo.txt file
# XXX_001 XXX_RNA-Seq_220308_Hela-rRNA-depletion-rep1
# XXX_002 XXX_RNA-Seq_220319_Hela-rRNA-depletion-rep2

usage() {
    echo "Usage: $0 <rawDataRawDir> <sampleInfo> <spikeIn> <expRef> <quantMethod>"
    echo "  rawDataRawDir: raw data directory"
    echo "  sampleInfo: space separated sample information file"
    echo "  spikeIn: if spike-in RNA-Seq, 'Y' or 'N' "
    echo "  expRef: experiment genome reference, 'hg38', 'mm10' "
    echo "  quantMethod: quantification method, 'featureCounts' or 'Salmon' "
    exit 1
}

# args 
rawDataRawDir=$1
sampleInfo=$2
spikeIn=$3
expRef=$4
quantMethod=$5

# check args
if [ $# -ne 5 ]; then
    echo "Error: Invalid number of arguments!"
    usage
fi
if [[ -d "$rawDataRawDir" && $(ls -A "$rawDataRawDir") ]]; then
    echo "Using $rawDataRawDir as input data."
else
    echo "Error: Invalid or empty rawDataRawDir directory!"
	usage
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
if [[ "$quantMethod" != "featureCounts" && "$quantMethod" != "Salmon" ]]; then
    echo "Error: Invalid quantMethod value! Must be 'featureCounts' or 'Salmon'."
    usage
fi

echo -e "RNA-seq parameters: \nrawDataRawDir: $rawDataRawDir\nsampleInfo: $sampleInfo\nspikeIn: $spikeIn\nexpRef: $expRef\nquantMethod: $quantMethod"

# output dirs
workDir=`pwd`
rawDataDir=${workDir}/00_raw_data
rawQcDir=${workDir}/01-1_raw_data_fastqc
trimDir=${workDir}/01-2_trimmed_data
cleanQcDir=${workDir}/01-3_trimmed_data_fastqc
map2ExpDir=${workDir}/02-1_experimental_genome_alignment
map2SpkDir=${workDir}/02-2_spikein_genome_alignment
filterExpBamDir=${workDir}/02-3_experimental_bam_qc
filterSpkBamDir=${workDir}/02-4_spikein_bam_qc
spkScaledBwDir=${workDir}/03_spikein_scaled_bw_track
cpmScaledBwDir=${workDir}/03_cpm_scaled_bw_track
quantificationDir=${workDir}/04_quantification

# log dirs
logDir=${workDir}/logs
trimLogDir=${logDir}/trim_galore_log
map2ExpLogDir=${logDir}/experimental_genome_alignment_STAR_log
map2SpkLogDir=${logDir}/spikein_genome_alignment_STAR_log
markDupLogDir=${logDir}/picard_MarkDuplicates_log
summaryDir=${logDir}/alignment_summary.txt
# filterExpLogDir=${logDir}/exp_filter
# filterSpkLogDir=${logDir}/spk_filter
spkScaledBwLogDir=${logDir}/spikein_scaled_bw_bamCoverage_log
cpmScaledBwLogDir=${logDir}/cpm_scaled_bw_bamCoverage_log
featureCountsLogDir=${logDir}/quantification_featureCounts_log
salmonLogDir=${logDir}/quantification_Salmon_log

# all tools in rnaseq env update at Oct 2024
# # fastqc version 0.12.1
fastqc='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/bin/fastqc'
# # multiqc version 1.25.1
multiqc='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/bin/multiqc'
# # trim_galore version 0.6.10 with cutadapt version 4.9
trim_galore='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/bin/trim_galore'
# # STAR version 2.7.11b
STAR='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/download_bin/STAR'
# # samtools version 1.21
samtools='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/download/samtools/bin/samtools'
# # Salmon version 1.10.0
salmon='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/download_bin/salmon'
# # picard version 3.3.0; after 3.0.0, picard need java 17
# picard='java -jar /share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/download/picard.jar'
# # Use old version, we do not need add read group for MarkDuplicates
# https://gatk.broadinstitute.org/hc/en-us/community/posts/24452455341467-Error-in-MarkDuplicates
# https://gatk.broadinstitute.org/hc/en-us/community/posts/14923447538715-Best-practice-with-read-groups
# https://www.biostars.org/p/9479469/
picard='/share/home/Grape/software/install_pkg/miniconda3/envs/mainenv/bin/picard'
# # featureCounts version 2.0.7
featureCounts='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/download/subread-2.0.7-Linux-x86_64/bin/featureCounts'
# # deepTools version 3.5.5
bamCoverage='/share/home/Grape/software/install_pkg/miniconda3/envs/rnaseq/bin/bamCoverage'

# deal ref genome: 
# hg38 gencode v39; mm10 gencode vM25; hg19 gencode v19?
hg38Index='/share/home/Grape/genome/summary/hg38_starIndex'
hg19Index='/share/home/Blueberry/reference/index/star/hg19_star_index'
mm10Index='/share/home/Grape/genome/summary/mm10_starIndex'
hg38Anno='/share/home/Grape/genome/Homo_sapiens/genecode/primary_assembly_GRCh38_v39/gencode.v39.primary_assembly.annotation.gtf'
hg19Anno='/share/home/Grape/genome/Homo_sapiens/genecode/gencode.v19.annotation.gtf'
mm10Anno='/share/home/Grape/genome/Mus_musculus/gencode/release_M25/gencode.vM25.primary_assembly.annotation.gtf'
hg38TranscriptFa='/share/home/Grape/genome/Homo_sapiens/genecode/primary_assembly_GRCh38_v39/gencode.v39.transcripts.gffread.fa'
hg19TranscriptFa=''
mm10TranscriptFa='/share/home/Grape/genome/Mus_musculus/gencode/release_M25/gencode.vM25.transcripts.gffread.fa'

if [[ $spikeIn == 'Y' ]]; then
    case $expRef in
        "hg38")
            expStarIndex=${hg38Index}
            spkStarIndex=${mm10Index}
            expAnno=${hg38Anno}
            expFa=${hg38TranscriptFa}
            exp_info="hg38"
            spike_info="mm10"
            ;;
        "hg19")
            expStarIndex=${hg19Index}
            spkStarIndex=${mm10Index}
            expAnno=${hg19Anno}
            expFa=${hg19TranscriptFa}
            exp_info="hg19"
            spike_info="mm10"
            ;;
        "mm10")
            expStarIndex=${mm10Index}
            spkStarIndex=${hg38Index}
            expAnno=${mm10Anno}
            expFa=${mm10TranscriptFa}
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
            expStarIndex=${hg38Index}
            expAnno=${hg38Anno}
            expFa=${hg38TranscriptFa}
            exp_info="hg38"
            ;;
        "hg19")
            expStarIndex=${hg19Index}
            expAnno=${hg19Anno}
            expFa=${hg19TranscriptFa}
            exp_info="hg19"
            ;;
        "mm10")
            expStarIndex=${mm10Index}
            expAnno=${mm10Anno}
            expFa=${mm10TranscriptFa}
            exp_info="mm10"
            ;;
        *)
            exit 1
            ;;
    esac
fi

echo -e "\n***************************\nRNA-seq processing at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
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
# Step1: fastq QC
echo -e "\n***************************\nTrimming and QC at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
# Step1.1: raw data qc
if [[ ! -d $rawQcDir ]]; then
    mkdir -p $rawQcDir
	$fastqc -t 48 --memory 1024 ${rawDataDir}/*.fq.gz -o $rawQcDir &> /dev/null
fi
if [[ ! -s ${rawQcDir}/rawdata_multiqc.html ]] || [[ ! -d ${rawQcDir}/rawdata_multiqc_data ]]; then
	$multiqc -f -n rawdata_multiqc -o $rawQcDir $rawQcDir
fi

# Step1.2: trim adapter and low quality sequence (trim_galore)
if [[ ! -d $trimDir ]]; then
    mkdir -p $trimDir
	mkdir -p $trimLogDir
	for r1 in `ls ${rawDataDir}/*_R1.fq.gz`; do
		r2=${r1/R1.fq.gz/R2.fq.gz}
		sampleName=$(basename ${r1%_R1.fq.gz})
		$trim_galore --phred33 -j 4 -o $trimDir -q 25 --length 30 -e 0.1 --stringency 4 \
		--paired $r1 $r2 &> ${trimLogDir}/${sampleName}_cutadapt.log &
	done
    echo -e "\n\nFinish trimming for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n"
fi
wait

# Step1.3: trimmed data qc
if [[ ! -d $cleanQcDir ]]; then
	mkdir -p $cleanQcDir
	$fastqc -t 48 --memory 1024 ${trimDir}/*.fq.gz -o $cleanQcDir &> /dev/null
fi
if [[ ! -s ${cleanQcDir}/cleanData_multiqc.html ]] || [[ ! -d ${cleanQcDir}/cleanData_multiqc_data ]]; then
	$multiqc -f -n cleanData_multiqc -o $cleanQcDir $trimLogDir $cleanQcDir 
fi

####
# Step2: Alignment and QC
# args ref by two famous and widely used pipeline: 1) Encode long RNA-Seq 2) nf-core RNA-Seq
# # For 2*150 reads --outFilterMismatchNmax 3 (equal as --outFilterMismatchNoverReadLmax 0.01) looks too stringent, although average sequencing error rate <0.01, but we want to allow for some SNPs.

run_star_alignment() {
    local genome_dir=$1          # --genomeDir
    local read1=$2               # --readFilesIn
    local read2=$3               
    local out_prefix=$4          # --outFileNamePrefix
    local n_multimap=$5          # --outFilterMultimapNmax
    local quant_mode=$6          # --quantMode
    local extra_params=$7        
    local log_file=$8          
    # local threadN=${9:-48}
    # local bamThreadN=${10:-48}

    $STAR --genomeDir $genome_dir \
        --readFilesIn $read1 $read2 \
        --outFileNamePrefix $out_prefix \
        --runThreadN 48 \
        --runRNGseed 777 \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 48 \
        --twopassMode Basic \
        --outFilterMismatchNmax 5 \
        --outFilterMultimapNmax $n_multimap \
        --quantMode $quant_mode \
        $extra_params \
        &> $log_file
}

# align to experimental
echo -e "\n***************************\nAligning to experimental at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
if [[ ! -d $map2ExpDir ]]; then
    mkdir -p $map2ExpDir
	mkdir -p $map2ExpLogDir
    if [[ ! `ls ${map2ExpDir}/*_${exp_info}_Aligned.sortedByCoord.out.bam 2> /dev/null` ]]; then
        for r1 in `ls ${trimDir}/*_1.fq.gz`;do
            r2=${r1/R1_val_1.fq.gz/R2_val_2.fq.gz}
            sampleName=$(basename ${r1%_R1_val_1.fq.gz})
            outPrefix=${map2ExpDir}/${sampleName}_${exp_info}_
            logName=${map2ExpLogDir}/${sampleName}_${exp_info}_STAR.log
            if [[ $quantMethod == 'featureCounts' ]]; then
                # for old lab pipeline, we always used uniquely mapped reads for gene expression, but in new pipe we conside the multimappers: see also https://www.biostars.org/p/9476989/
                run_star_alignment $expStarIndex $r1 $r2 $outPrefix 1 "GeneCounts" "" $logName
            elif [[ $quantMethod == 'Salmon' ]]; then
                # ignore '--outSAMstrandField intronMotif' which suggested for non-stranded RNA-seq
                # WARNING: after STAR 2.7.11b '--quantTranscriptomeBan Singleend' was substituted by '--quantTranscriptomeSAMoutput BanSingleEnd'. see https://github.com/alexdobin/star/releases
                run_star_alignment $expStarIndex $r1 $r2 $outPrefix 20 "TranscriptomeSAM GeneCounts" "--quantTranscriptomeSAMoutput BanSingleEnd --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outSAMattributes NH HI AS NM MD" $logName
            else
                exit 1
            fi
            echo -e "Finish alignment for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
        done
    fi
fi
# align to spike-in
if [[ $spikeIn == 'Y' ]];then
    echo -e "\n***************************\nAligning to spike-in at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    if [[ ! -d $map2SpkDir ]];then
        mkdir -p $map2SpkDir
        mkdir -p $map2SpkLogDir
        if [[ ! `ls ${map2SpkDir}/*_${spike_info}_Aligned.sortedByCoord.out.bam 2> /dev/null` ]]; then
            for r1 in `ls ${trimDir}/*_1.fq.gz`;do
                r2=${r1/R1_val_1.fq.gz/R2_val_2.fq.gz}
                sampleName=$(basename ${r1%_R1_val_1.fq.gz})
                outPrefix=${map2SpkDir}/${sampleName}_${spike_info}_
                logName=${map2SpkLogDir}/${sampleName}_${spike_info}_STAR.log
                if [[ $quantMethod == 'featureCounts' ]]; then
                    run_star_alignment $spkStarIndex $r1 $r2 $outPrefix 1 "GeneCounts" "" $logName
                elif [[ $quantMethod == 'Salmon' ]]; then
                    run_star_alignment $spkStarIndex $r1 $r2 $outPrefix 20 "TranscriptomeSAM GeneCounts" "--quantTranscriptomeSAMoutput BanSingleEnd --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outSAMattributes NH HI AS NM MD" $logName
                else
                    exit 1
                fi
                echo -e "Finish alignment for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
            done
        fi
    fi
fi

# DO NOT REMOVE RUP FOR RNA-SEQ
REMOVE_DUP="false"

if [[ ! -d $filterExpBamDir ]];then
    mkdir -p $filterExpBamDir
	mkdir -p $markDupLogDir
	if [[ ! `ls ${filterExpBamDir}/*_${exp_info}.markdup.metrics 2> /dev/null` ]]; then
        echo -e "\n***************************\nMarking duplicate for experimental at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
		ls ${map2ExpDir}/*_${exp_info}_Aligned.sortedByCoord.out.bam | while read bam; do

            # do not rmdup but mark them: https://nf-co.re/rnaseq/3.2/docs/output
            # we do not want the secondary and supplementary reads to be duplicate marked as well as the primary reads, so we do not get queryname sorted file
            # Mark duplicate at genome level (for Aligned.sortedByCoord.out.bam) but not at transcriptome level (for Aligned.toTranscriptome.out.bam): https://groups.google.com/g/rna-star/c/fTLo7vkJhWg
            sampleName=$(basename ${bam%_${exp_info}_Aligned.sortedByCoord.out.bam})
            markDupBam=${filterExpBamDir}/${sampleName}_${exp_info}.markdup.bam
            markDupMetric=${filterExpBamDir}/${sampleName}_${exp_info}.markdup.metrics

			$picard MarkDuplicates \
                --ASSUME_SORT_ORDER coordinate \
                --VALIDATION_STRINGENCY LENIENT \
                --REMOVE_DUPLICATES $REMOVE_DUP \
			    --INPUT $bam \
			    --OUTPUT $markDupBam \
			    --METRICS_FILE $markDupMetric &> ${markDupLogDir}/${sampleName}_${exp_info}_markdup.log
            
            filteredGenomeBam=${filterExpBamDir}/${sampleName}_${exp_info}.markdup.paired.bam
            $samtools view -@ 48 -bS -f 2 -o $filteredGenomeBam $markDupBam
            $samtools index -@ 48 $filteredGenomeBam
            bamGenomeStat=${filterExpBamDir}/${sampleName}_${exp_info}.genome.flagstat
            $samtools flagstat -@ 48 $filteredGenomeBam > $bamGenomeStat

            if [[ $quantMethod == 'featureCounts' ]]; then
                filteredBam=$filteredGenomeBam
            # do not sort by coordinates for transcriptome
            elif [[ $quantMethod == 'Salmon' ]]; then
                transcriptomeBam=${map2ExpDir}/${sampleName}_${exp_info}_Aligned.toTranscriptome.out.bam
                filteredBam=${filterExpBamDir}/${sampleName}_${exp_info}_transcriptome.paired.bam
                $samtools view -@ 48 -bS -f 2 -o $filteredBam $transcriptomeBam
                bamTranscriptomeStat=${filterExpBamDir}/${sampleName}_${exp_info}.transcriptome.flagstat
                $samtools flagstat -@ 48 $filteredBam > $bamTranscriptomeStat
            else
                exit 1
            fi
            echo -e "Finish MarkDuplicates, filter singletons and get statistics(flagstat) of $(basename ${filteredBam}) for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
		done
        $multiqc -f -n exp_alignment_multiqc -o $filterExpBamDir -x "*STARpass1" $filterExpBamDir $map2ExpDir
	fi
fi
# spikeIn for genome
if [[ $spikeIn == 'Y' ]];then
    mkdir -p $filterSpkBamDir
    if [[ ! `ls ${filterSpkBamDir}/*_${spike_info}.markdup.metrics 2> /dev/null` ]]; then
        echo -e "\n***************************\nMarking duplicate for spike-in at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
        ls ${map2SpkDir}/*_${spike_info}_Aligned.sortedByCoord.out.bam | while read bam; do
            sampleName=$(basename ${bam%_${spike_info}_Aligned.sortedByCoord.out.bam})
            markDupBam=${filterSpkBamDir}/${sampleName}_${spike_info}.markdup.bam
            markDupMetric=${filterSpkBamDir}/${sampleName}_${spike_info}.markdup.metrics

            $picard MarkDuplicates \
                --ASSUME_SORT_ORDER coordinate \
                --VALIDATION_STRINGENCY LENIENT \
                --REMOVE_DUPLICATES $REMOVE_DUP \
                --INPUT $bam \
                --OUTPUT $markDupBam \
                --METRICS_FILE $markDupMetric &> ${markDupLogDir}/${sampleName}_${spike_info}_markdup.log
            
            filteredGenomeBam=${filterSpkBamDir}/${sampleName}_${spike_info}.markdup.paired.bam
            # -q 30 (any value >3 and < 255) for only keep unique reads for scalefactor calculate!
            $samtools view -@ 48 -bS -q 30 -f 2 -o $filteredGenomeBam $markDupBam
            $samtools index -@ 48 $filteredGenomeBam
            bamGenomeStat=${filterSpkBamDir}/${sampleName}_${spike_info}.genome.flagstat
            $samtools flagstat -@ 48 $filteredGenomeBam > $bamGenomeStat
            echo -e "Finish MarkDuplicates, filter singletons and get statistics(flagstat) of $(basename ${filteredGenomeBam}) for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
        done
        $multiqc -f -n spk_alignment_multiqc -o $filterSpkBamDir -x "*STARpass1" $filterSpkBamDir $map2SpkDir
    fi
fi

# Step3: calculate spike-in scale factor
# # For RNA-Seq, scale factor is CPM of spike-in mapped reads
if [[ ! -d $summaryDir ]]; then
    echo -e "\n***************************\nCalculating alignment summary at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    echo -e "sample_name\tall_reads\t${exp_info}_reads\t${exp_info}_mapping_ratio\t${exp_info}_qc_reads\t${exp_info}_qc_ratio\t${spike_info}_qc_reads\t${spike_info}_qc_ratio\t${spike_info}_unique_reads\t${spike_info}_unique_ratio\tscale_qc_factor\tscale_unique_factor" > $summaryDir
    basename -a $(ls ${map2ExpDir}/*_${exp_info}_Log.final.out) | while read file; do
        sampleName=${file/_${exp_info}_Log.final.out/}
        allReads=`grep "Number of input reads" ${map2ExpDir}/${file} | awk '{print $NF*2}'`
        # unique read and multi-mapping reads (pass STAR filter(multiple loci) and not pass(to too many loci))
    	expReads=`grep -e "Number.* mapped" -e "mapped.*number" ${map2ExpDir}/${file} | cut -f 2 | awk '{sum+=$0};END{print sum*2}'`
		expMapRatio=`grep " mapped.*%" ${map2ExpDir}/${file} | cut -f 2 | awk '{sum+=$0};END{print sum}'`
        # remove singleton and not pass STAR multi-mapping filter reads
        expQcReads=`echo "2*$(cat ${filterExpBamDir}/${sampleName}_${exp_info}.genome.flagstat | grep "read1" | cut -d " " -f 1)" | bc -l`
        expQcRatio=`printf "%.2f\n" $(echo "100*${expQcReads}/${allReads}" | bc -l)`
        spkQcReads=`echo "2*$(cat ${filterSpkBamDir}/${sampleName}_${spike_info}.genome.flagstat | grep "read1" | cut -d " " -f 1)" | bc -l`
        spkQcRatio=`printf "%.2f\n" $(echo "100*${spkQcReads}/${allReads}" | bc -l)`
        # Unique map to spkin but not map to exp (to avoid homogenous)
        spkUniqueReads=`echo "2*$(comm --check-order -23 <(samtools view -@ 12 ${map2SpkDir}/${sampleName}_${spike_info}_Aligned.sortedByCoord.out.bam | cut -f 1 | sort -S48G --parallel=24 | uniq) \
            <(samtools view -@ 12 ${map2ExpDir}/${sampleName}_${exp_info}_Aligned.sortedByCoord.out.bam | cut -f 1 | sort -S62G --parallel=36 | uniq) | wc -l)" | bc -l`
        spkUniqueRatio=`printf "%.2f\n" $(echo "100*${spkUniqueReads}/${allReads}" | bc -l)`
        scaleQcFactor=`echo "1000000/${spkQcReads}" | bc -l`
        scaleUniqueFactor=`echo "1000000/${spkUniqueReads}" | bc -l`
        echo -e ${sampleName}"\t"${allReads}"\t"${expReads}"\t"${expMapRatio}"\t"${expQcReads}"\t"${expQcRatio}"\t"${spkQcReads}"\t"${spkQcRatio}"\t"${spkUniqueReads}"\t"${spkUniqueRatio}"\t"${scaleQcFactor}"\t"${scaleUniqueFactor} >> $summaryDir
    done
    echo -e "Finish calculate alignment summary at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
fi

# Step4: get bw track
if [[ $spikeIn == 'Y' ]];then
	if [[ ! -d $spkScaledBwDir ]];then
		mkdir -p $spkScaledBwDir
		mkdir -p $spkScaledBwLogDir
        echo -e "\n***************************\nGetting spikein normalized track at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
		cat ${logDir}/alignment_summary.txt | sed '1d' | while read line;do
			arr=($line)
			sampleName=${arr[0]}
            # if use unique: scaleFactor=${arr[11]}
			scaleFactor=${arr[10]}
			bam=${filterExpBamDir}/${sampleName}_${exp_info}.markdup.paired.bam
			if [[ ! -s "${spkScaledBwDir}/${sampleName}_ds.bw" ]];then
                echo "Generating file: ${spkScaledBwDir}/${sampleName}_ds.bw"
        		$bamCoverage --bam $bam --skipNonCoveredRegions --binSize 1 --numberOfProcessors 48 --outFileName ${spkScaledBwDir}/${sampleName}_ds.bw --scaleFactor $scaleFactor --normalizeUsing None &> ${spkScaledBwLogDir}/${sampleName}_ds.log
			fi
            if [[ ! -s "${spkScaledBwDir}/${sampleName}_fwd.bw" ]];then
                echo "Generating file: ${spkScaledBwDir}/${sampleName}_fwd.bw"
                $bamCoverage --bam $bam --skipNonCoveredRegions --filterRNAstrand forward --binSize 1 --numberOfProcessors 48 --outFileName ${spkScaledBwDir}/${sampleName}_fwd.bw --scaleFactor $scaleFactor --normalizeUsing None &> ${spkScaledBwLogDir}/${sampleName}_fwd.log
            fi
            if [[ ! -s "${spkScaledBwDir}/${sampleName}_rev.bw" ]];then
                echo "Generating file: ${spkScaledBwDir}/${sampleName}_rev.bw"
                $bamCoverage --bam $bam --skipNonCoveredRegions --filterRNAstrand reverse --binSize 1 --numberOfProcessors 48 --outFileName ${spkScaledBwDir}/${sampleName}_rev.bw --scaleFactor $scaleFactor --normalizeUsing None &> ${spkScaledBwLogDir}/${sampleName}_rev.log
            fi
            echo -e "Finish get track for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
		done
    fi
else
	if [[ ! -d $cpmScaledBwDir ]];then
		mkdir -p $cpmScaledBwDir
		mkdir -p $cpmScaledBwLogDir
        echo -e "\n***************************\nGetting CPM normalized track at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
		ls ${filterExpBamDir}/*_${exp_info}.markdup.paired.bam | while read sample;do
			sampleName=$(basename ${sample%_${exp_info}.markdup.paired.bam})
			if [[ ! -s "${cpmScaledBwDir}/${sampleName}_ds.bw" ]];then
                echo "Generating file: ${cpmScaledBwDir}/${sampleName}_ds.bw"
				$bamCoverage --bam $sample --skipNonCoveredRegions --binSize 1 --numberOfProcessors 48 --outFileName ${cpmScaledBwDir}/${sampleName}_ds.bw --scaleFactor 1 --normalizeUsing CPM &> ${cpmScaledBwLogDir}/${sampleName}_ds.log
			fi
            if [[ ! -s "${cpmScaledBwDir}/${sampleName}_fwd.bw" ]];then
                echo "Generating file: ${cpmScaledBwDir}/${sampleName}_fwd.bw"
                $bamCoverage --bam $bam --skipNonCoveredRegions --filterRNAstrand forward --binSize 1 --numberOfProcessors 48 --outFileName ${cpmScaledBwDir}/${sampleName}_fwd.bw --scaleFactor 1 --normalizeUsing CPM &> ${cpmScaledBwLogDir}/${sampleName}_fwd.log
            fi
            if [[ ! -s "${cpmScaledBwDir}/${sampleName}_rev.bw" ]];then
                echo "Generating file: ${cpmScaledBwDir}/${sampleName}_rev.bw"
                $bamCoverage --bam $bam --skipNonCoveredRegions --filterRNAstrand reverse --binSize 1 --numberOfProcessors 48 --outFileName ${cpmScaledBwDir}/${sampleName}_rev.bw --scaleFactor 1 --normalizeUsing CPM &> ${cpmScaledBwLogDir}/${sampleName}_rev.log
            fi
            echo -e "Finish get track for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
		done
	fi
fi

# Step5: Quantification
# # since featureCounts v2.0.2, both -p and --countReadPairs need to be used to explicitly specify count read pairs (fragments) instead of reads
# -s 2 to set reversely stranded library
if [[ ! -d $quantificationDir ]];then
    mkdir -p $quantificationDir
    echo -e "\n***************************\nDoing quantification with ${quantMethod} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
    if [[ $quantMethod == 'featureCounts' ]]; then
        mkdir -p $featureCountsLogDir
        ls ${filterExpBamDir}/*_${exp_info}.markdup.paired.bam | while read bam; do
		    sampleName=$(basename ${bam%_${exp_info}.markdup.paired.bam})
            if [[ ! -s "${quantificationDir}/${sampleName}_featureCounts.counts" ]]; then
                $featureCounts \
                -T 48 \
                -a $expAnno \
                -F GTF \
                -t exon \
                -g gene_id \
                -s 2 \
                -o ${quantificationDir}/${sampleName}_featureCounts.counts \
                -p --countReadPairs \
                ${bam} &> ${featureCountsLogDir}/${sampleName}_featureCounts.log
            fi
            echo -e "Finish quantification for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
	    done
    elif [[ $quantMethod == 'Salmon' ]]; then
        mkdir -p $salmonLogDir
        ls ${filterExpBamDir}/*_${exp_info}_transcriptome.paired.bam | while read bam; do
            sampleName=$(basename ${bam%_${exp_info}_transcriptome.paired.bam})
            if [[ ! -d "${quantificationDir}/${sampleName}_salmon_quant" ]]; then
                $salmon quant \
                -t $expFa \
                --libType A \
                -a $bam \
                --threads 48 \
                --gcBias --seqBias \
                -o ${quantificationDir}/${sampleName}_salmon_quant &> ${salmonLogDir}/${sampleName}_salmon.log
            fi
            echo -e "Finish quantification for ${sampleName} at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)"
        done
    else
        exit 1
    fi
fi
