# SEQPIPE: A collection of bash pipelines for sequencing data upstream analysis.

A collection of bash pipelines for NGS data upstream analysis from the [ChenLab](https://chenf-lab.fudan.edu.cn/)

**NGS data:**<br>
- [RNA-seq](#RNA-seq)
- [PRO-seq](#PRO-seq)
- [TT-seq](#TT-seq)
- [ATAC-seq](#ATAC-seq)
- [ChIP-seq](#ChIP-seq)
- [CUT&Tag](#CUT&Tag)
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
```

## RNA-seq
RNA-seq data with (or without) spike-in normalizaton. <br>
Alignment: [STAR](https://github.com/alexdobin/STAR) <br>
Quantification: [Salmon](https://github.com/COMBINE-lab/salmon) or featureCounts(https://subread.sourceforge.net/featureCounts.html)

> Run RNA-seq pipe
```bash
nohup bash RNASEQ_STAR_Salmon_Spikein.sh /chenfeilab/Gaux/rawDataBackup/test/241018_RNA-Seq sampleInfo.txt N hg38 Salmon &> ./241018_RNASeq.log &
```

## PRO-seq
PRO-seq data with(or without) spike-in normalizaton for [qPRO-seq protocol](https://www.biorxiv.org/content/10.1101/2020.05.18.102277v1.full) and [rPRO-seq protocol](https://www.biorxiv.org/content/10.1101/2024.05.08.593182v1). <br>
Alignment: [Bowtie2](https://github.com/BenLangmead/bowtie2) <br>
Peak Calling: [dREG](https://github.com/Danko-Lab/dREG) and [PINTS](https://github.com/hyulab/PINTS)

> Run PRO-seq pipe
```bash
nohup bash PROSEQ...
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

## CUT&Tag
...
...
