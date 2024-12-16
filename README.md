# SEQPIPE: A collection of bash pipelines for sequencing data upstream analysis.

A collection of bash pipelines for NGS data upstream analysis from the [ChenLab](https://chenf-lab.fudan.edu.cn/)

**NGS data:**<br>
- [RNA-seq](#RNA-seq)
- [PRO-seq/cap](#PRO-seqcap)
- [TT-seq](#TT-seq)
- [ATAC-seq](#ATAC-seq)
- [ChIP-seq](#ChIP-seq)
- [CUT&Tag/RUN](#cuttagrun)
...

Correspondence: ShaoxuanWang@hotmail.com

## Reference genome and annotation
> **hg38**: [gencode release 39](https://www.gencodegenes.org/human/release_39.html)
```bash
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.primary_assembly.genome.fa.gz
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.primary_assembly.annotation.gtf.gz
# PE150 STAR index
STAR --runMode genomeGenerate --runThreadN 48 --genomeDir ./PE150_gencode_v39 --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile gencode.v39.primary_assembly.annotation.gtf --sjdbOverhang 149 --limitGenomeGenerateRAM 549755813888
# transcriptome fasta
gffread -w gencode.v39.transcripts.gffread.fa -g GRCh38.primary_assembly.genome.fa gencode.v39.primary_assembly.annotation.gtf
# Bowtie2 index
bowtie2-build -f GRCh38.primary_assembly.genome.fa --threads 48 ./hg38
# Bowtie2 human rDNA: https://github.com/databio/ref_decoy/blob/master/human_rDNA.fa.gz
bowtie2-build -f human_rDNA.fa.gz --threads 12 ./human_rDNA
```
> **mm10**: [gencode release M25](https://www.gencodegenes.org/mouse/release_M25.html)
```bash
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz
# PE150 STAR index
STAR --runMode genomeGenerate --runThreadN 48 --genomeDir ./PE150_gencode_m25 --genomeFastaFiles GRCm38.primary_assembly.genome.fa --sjdbGTFfile gencode.vM25.primary_assembly.annotation.gtf --sjdbOverhang 149 --limitGenomeGenerateRAM 549755813888
# transcriptome fasta
gffread -w gencode.vM25.transcripts.gffread.fa -g GRCm38.primary_assembly.genome.fa gencode.vM25.primary_assembly.annotation.gtf
# Bowtie2 index
bowtie2-build -f GRCm38.primary_assembly.genome.fa --threads 48 ./mm10
# Bowtie2 mouse rDNA: https://github.com/databio/ref_decoy/blob/master/mouse_rDNA.fa.gz
bowtie2-build -f mouse_rDNA.fa.gz --threads 12 ./mouse_rDNA
```

## RNA-seq
RNA-seq data with (or without) spike-in normalizaton. <br>
Alignment: [STAR](https://github.com/alexdobin/STAR) <br>
Quantification: [Salmon](https://github.com/COMBINE-lab/salmon) or [featureCounts](https://subread.sourceforge.net/featureCounts.html)

> Run RNA-seq pipe
```
bash RNASEQ_STAR_Salmon_Spikein.sh <rawDataRawDir> <sampleInfo> <runInfo> <spikeIn> <expRef> <quantMethod> <bw>

<rawDataRawDir>: raw data directory
<sampleInfo>: space separated sample information file
<runInfo>: a description for the run, e.g., 'RNAseq_241120'
<spikeIn>: if spike-in RNA-Seq, 'Y' or 'N'
<expRef>: experiment genome reference, 'hg38', 'mm10'
<quantMethod>: quantification method, 'featureCounts' or 'Salmon'
<bw>: if bigwig for RNA-Seq, 'Y', or 'N'
```
For example,
```bash
nohup bash RNASEQ_STAR_Salmon_Spikein.sh /chenfeilab/Gaux/rawDataBackup/test/241018_RNA-Seq sampleInfo.txt 241018_RNASEQ N hg38 Salmon N &> 241018_RNASeq.log &
# If your file path contain wildcards, please enclose the path in double quotes to prevent bash wildcard expansion.
nohup bash RNASEQ_STAR_Salmon_Spikein.sh "/chenfeilab/Gaux/rawDataBackup/*/*" sampleInfo.txt 241111_RNASEQ Y hg38 featureCounts N &> 241111_RNASeq.log &
```

## PRO-seq/cap
PRO-seq data with (or without) spike-in normalizaton for [qPRO-seq protocol](https://www.biorxiv.org/content/10.1101/2020.05.18.102277v1.full) and [rPRO-seq protocol](https://www.biorxiv.org/content/10.1101/2024.05.08.593182v1).<br>
PRO-cap is [PRO-cap protocol](https://www.nature.com/articles/nprot.2016.086) but a library struceture same as qPRO-seq.<br>
Alignment: [Bowtie2](https://github.com/BenLangmead/bowtie2) <br>
Peak Calling: [dREG](https://github.com/Danko-Lab/dREG) and [PINTS](https://github.com/hyulab/PINTS)

> Run PRO-seq pipe
```
bash PROSEQ_rPRO_qPRO_Spikein.sh <rawDataRawDir> <sampleInfo> <runInfo> <spikeIn> <expRef> <libType> <umiLen> <identifyTRE> <noNormBw>

<rawDataRawDir>: raw data directory
<sampleInfo>: space separated sample information file
<runInfo>: a description for the run, e.g., 'PROseq_241130'
<spikeIn>: if spike-in PRO-seq library, 'Y' or 'N'
<expRef>: experiment genome reference, 'hg38', 'hg19' or 'mm10'
<libType>: PRO-seq library type, 'qPRO', 'rPRO' or 'PROcap'
<umiLen>: the length (n nt) of UMI, between 6 and 12
<identifyTRE>: peak calling method, 'dREG', 'PINTS', 'all' or 'none'
<noNormBw>: if get no normalized bigwig (only for low sequencing depth sample to merge), 'Y' or 'N'
```
For example,
```bash
nohup bash PROSEQ_rPRO_qPRO_Spikein.sh /chenfeilab/Gaux/rawDataBackup/xxx/xxxxxx_PRO-seq sampleInfo.txt 241108_PROSEQ N hg38 rPRO 6 all Y &> 241105_rPROseq.log &
nohup bash PROSEQ_rPRO_qPRO_Spikein.sh "/chenfeilab/Pomelo/try/SXY/24rawdata/*/*" sampleInfo.txt SXY51to60_qPRO Y mm10 qPRO 6 none N &> SXY51to60_qPROseq.log &
nohup bash PROSEQ_rPRO_qPRO_Spikein.sh /chenfeilab/Gaux/rawDataBackup/xxx/xxxxxx_PRO-cap sampleInfo.txt 241130_PROCAP N hg38 PROcap 6 none N &> 241130_PROcap.log &
```

## TT-seq
...
...


## ATAC-seq
...
...


## ChIP-seq
...
...


## CUT&Tag/RUN
CUT&Tag and CUT&RUN data with (or without) spike-in normalizaton and IgG negative control.<br>
Spike-in strategy: fixed DNA or nuclei (randomly DNA) <br>
Alignment: [Bowtie2](https://github.com/BenLangmead/bowtie2) <br>
Peak calling: [MACS3](https://github.com/macs3-project/MACS) and [SEACR](https://github.com/FredHutch/SEACR)

> Run CUT&Tag CUT&RUN pipe
```
bash CUTTAG_CUTRUN_Spikein_IgG.sh <rawDataRawDir> <sampleInfo> <runInfo> <spikeIn> <expRef> <spkRef> <spkStrategy> <libType> <rmDup> <controlIgG> <callPeak> <peakType>

<rawDataRawDir>: raw data directory. If path contains wildcards, enclose the path in quotes to prevent bash wildcard expansion
<sampleInfo>: space separated sample information file, must have 2 (rawName newName) columns or 3 (rawName newName controlName) columns (if controlIgG is 'Y')
<runInfo>: a description for the run, e.g., 'CUTTAG_H3K4me3_241230' 
<spikeIn>: if spike-in CUTTAG library, 'Y' or 'N'
<expRef>: experiment genome reference, 'hg38', 'hg19', 'mm10'
<spkRef>: spike-in genome reference, 'dm6', 'k12', 'hg38', 'mm10', 'none'
<spkStrategy>: spike-in normalization strategy, 'DNA', 'nuclei' or 'none'. DNA refers to a fixed DNA sequence, as included in many kits. (However, if the spike-in is randomly fragmented DNA, please set it as 'nuclei')
<libType>: library type, 'CUTTAG' or 'CUTRUN'
<rmDup>: remove duplicates or just mark them, 'remove' or 'mark'
<controlIgG>: if IgG as negative control, 'Y' or 'N'
<callPeak>: peak calling method, 'MACS3', 'SEACR' or 'none'
<peakType>: peak type, 'narrow' or 'broad'
```