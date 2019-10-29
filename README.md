This repository contains data access, source code and demo of the custom computational analyses developed as part of our paper **Sequence preference of DNMT3A/B dictates genome-wide DNA methylation signature** *(under review)*


## Data

All the raw data from the E.coli bisulfite sequencing experiments have been deposited in the ArrayExpress database at EMBL-EBI under accession number [E-MTAB-8406](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8406) (*currently only accesible to referees - login details included in the cover letter*)

Sample name | Organism | Treatment | Number of libraries
------------|----------|-----------|--------------------
DNMT3A-Abcam-30min | unmethylated Escherichia coli str. K-12 substr. MG1655 | human DNMT3A	(30 min) | 4 
DNMT3A-Abcam-120min | unmethylated Escherichia coli str. K-12 substr. MG1655 | human DNMT3A	(120 min) | 4 
DNMT3A-Abcam-240min | unmethylated Escherichia coli str. K-12 substr. MG1655 | human DNMT3A	(240 min) | 4
DNMT3B-Abcam-30min | unmethylated Escherichia coli str. K-12 substr. MG1655 | human DNMT3B	(30 min) | 4
DNMT3B-Abcam-120min | unmethylated Escherichia coli str. K-12 substr. MG1655 | human DNMT3B	(120 min) | 4
DNMT3B-Abcam-240min | unmethylated Escherichia coli str. K-12 substr. MG1655 | human DNMT3B	(240 min) | 4
MSssI-10min | unmethylated Escherichia coli str. K-12 substr. MG1655 | bacterial MSssI (10 min) | 4
MSssI-30min | unmethylated Escherichia coli str. K-12 substr. MG1655 | bacterial MSssI (30 min) | 4
MSssI-240min | unmethylated Escherichia coli str. K-12 substr. MG1655 | bacterial MSssI (240 min) | 4
Ecoli-unmethylated-DNA | unmethylated Escherichia coli str. K-12 substr. MG1655 | none | 4	



## Code

### Contents

- [System requirements and installation](README.md#system-requirements-and-installation)
- [Instructions for use and demo](README.md#instructions-for-use-and-demo)
  - [Sequence context diversity in reference genomes](README.md#sequence-context-diversity-in-reference-genomes)
  - [Quality check](README.md#quality-check)
  - [Trimming](README.md#trimming)
  - [Prepare reference genome](README.md#prepare-reference-genome)
  - [Alignment](README.md#alignment)
  - [Deduplication](README.md#deduplication)
  - [Extract methylation](README.md#extract-methylation)
  - [Sequence context plots](README.md#sequence-context-plots)


### System requirements and installation

Software, for installation details of the individual tools follow the links:

- [EMBOSS v6.6.0.0](http://emboss.sourceforge.net/)
- [fastaRegexFinder.py v0.1.1](https://github.com/dariober/bioinformatics-cafe/tree/master/fastaRegexFinder)
- Standard Unix tools: cut, awk, sort, uniq, grep, tr, echo ...
- [bedtools v2.27.0](http://bedtools.readthedocs.io/en/latest/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [bismark v0.19.0](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
- [tableCat.py](https://github.com/dariober/bioinformatics-cafe/blob/master/tableCat/tableCat.py)
- [R v3.3.2](https://www.r-project.org/). Libraries:
  - [data.table v1.10.4](https://cran.r-project.org/web/packages/data.table/index.html)
  - [ggplot2 v2.2.1](http://ggplot2.org/)

Operating system:

  - CentOS Linux release 7.3.1611 (OS used during code development)
  - [slurm job scheduling system v19.05.0](https://slurm.schedmd.com/quickstart.html)


### Instructions for use and demo

#### Sequence context diversity in reference genomes

Calculating all available sequence contexts in the lambda, E.coli and human genomes as follows:

```bash
# define the reference genome
ref=genome.fa # either lamba, E.coli or human

# reverse complement the reference genome
revseq $ref -reverse -complement -notag -outseq genome_revcom.fa

# generate index for reference genomes
samtools faidx $ref
samtools faidx genome_revcom.fa


##########
## NCGN ##
##########

# get all CGs on the forward strand
fastaRegexFinder.py -f $ref -r CG --noreverse | cut -f1-3 | bedtools sort -i > cg_plus.bed

# get all CGs on the reverse strand
fastaRegexFinder.py -f genome_revcom.fa -r CG --noreverse | cut -f1-3 | bedtools sort -i > cg_minus.bed

# calculate all CG contexts in both forward and reverse strands
for b in {1..8}
do
  l=`echo $b*2+2 | bc`
  c=`echo "4^($b*2)" | bc`
  bedtools slop -i cg_plus.bed -g $ref.fai -b $b | bedtools getfasta -fi $ref -bed - | awk 'NR % 2 == 0' |  tr '[:lower:]' '[:upper:]' | grep -v "N" | grep -x ".\{$l\}" | sort | uniq > cg_plus.$b
  bedtools slop -i cg_minus.bed -g genome_revcom.fa.fai -b $b | bedtools getfasta -fi genome_revcom.fa -bed - | awk 'NR % 2 == 0' |  tr '[:lower:]' '[:upper:]' | grep -v "N" | grep -x ".\{$l\}" | sort | uniq > cg_minus.$b
  m=`cat cg_plus.$b cg_minus.$b | sort | uniq -c | wc -l`
  rm cg_plus.$b cg_minus.$b
  echo "CpG +/-" $b "bp:" $m "out of" $c "combinations (" `echo "scale=2; 100*$m/$c" | bc` "%)"
done


##########
## NCAN ##
##########

# get all CAs on the forward strand
fastaRegexFinder.py -f $ref -r CA --noreverse | cut -f1-3 | bedtools sort -i > ca_plus.bed

# get all CAs on the reverse strand
fastaRegexFinder.py -f genome_revcom.fa -r CA --noreverse | cut -f1-3 | bedtools sort -i > ca_minus.bed

# calculate all CA contexts in both forward and reverse strands
for b in {1..8}
do
  l=`echo $b*2+2 | bc`
  c=`echo "4^($b*2)" | bc`
  bedtools slop -i ca_plus.bed -g $ref.fai -b $b | bedtools getfasta -fi $ref -bed - | awk 'NR % 2 == 0' |  tr '[:lower:]' '[:upper:]' | grep -v "N" | grep -x ".\{$l\}" | sort | uniq > ca_plus.$b
  bedtools slop -i ca_minus.bed -g genome_revcom.fa.fai -b $b | bedtools getfasta -fi genome_revcom.fa -bed - | awk 'NR % 2 == 0' |  tr '[:lower:]' '[:upper:]' | grep -v "N" | grep -x ".\{$l\}" | sort | uniq > ca_minus.$b
  m=`cat ca_plus.$b ca_minus.$b | awk 'NR % 2 == 0' | grep -x ".\{$l\}" | sort | uniq -c | wc -l`
  rm ca_plus.$b ca_minus.$b
  echo "CpA +/-" $b "bp:" $m "out of" $c "combinations (" `echo "scale=2; 100*$m/$c" | bc` "%)"
done


#########
## NCN ##
#########

# get all Cs on the forward strand
fastaRegexFinder.py -f $ref -r C --noreverse | cut -f1-3 | bedtools sort -i > c_plus.bed

# get all Cs on the reverse strand
fastaRegexFinder.py -f genome_revcom.fa -r C --noreverse | cut -f1-3 | bedtools sort -i > c_minus.bed

# calculate all C contexts in both forward and reverse strands
for b in {1..8}
do
  l=`echo $b*2+1 | bc`
  c=`echo "4^($b*2)" | bc`
  bedtools slop -i c_plus.bed -g $ref.fai -b $b | bedtools getfasta -fi $ref -bed - | awk 'NR % 2 == 0' |  tr '[:lower:]' '[:upper:]' | grep -v "N" | grep -x ".\{$l\}" | sort | uniq > c_plus.$b
  bedtools slop -i c_minus.bed -g genome_revcom.fa.fai -b $b | bedtools getfasta -fi genome_revcom.fa -bed - | awk 'NR % 2 == 0' |  tr '[:lower:]' '[:upper:]' | grep -v "N" | grep -x ".\{$l\}" | sort | uniq > c_minus.$b
  m=`cat c_plus.$b c_minus.$b | awk 'NR % 2 == 0' | grep -x ".\{$l\}" | sort | uniq -c | wc -l`
  rm c_plus.$b c_minus.$b
  echo "C +/-" $b "bp:" $m "out of" $c "combinations (" `echo "scale=2; 100*$m/$c" | bc` "%)"
done
```


#### Quality check

`*.fastq.gz` files to be downloaded from [E-MTAB-8406](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8406) and located in the `~/fastq` folder:

```bash
cd ~/fastq

mkdir ../fastqc

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../fastqc/$bname.log --mem 4G --wrap "fastqc --noextract --nogroup -q -o ../fastqc $fq"
done
```


#### Trimming

```bash
cd ~/fastq

mkdir ../fastq_trimmed

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../fastq_trimmed/$bname.log --mem 4G --wrap "cutadapt -a AGATCGGAAGAGC -m 10 -u 6 -o ../fastq_trimmed/$fq $fq > ../fastq_trimmed/$bname.txt"
done
```


#### Prepare reference genome

Genome and annotations for `Escherichia coli str. K-12 substr. MG1655` obtained from [Ensembl Bacteria release 42](http://bacteria.ensembl.org/Escherichia_coli_str_k_12_substr_mg1655/Info/Index) and located in `~/reference`

```bash
cd ~/reference
bismark_genome_preparation .
```


#### Alignment

```bash
cd ~/fastq_trimmed

mkdir ../bismark

ref=../reference

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../bismark/$bname.log --mem 32G --wrap "bismark --non_directional --unmapped -o ../bismark -p 4 $ref $fq"
done
```


#### Deduplication

```bash
cd ~/bismark

for id in `ls *_R1_001_bismark_bt2.bam | cut -d "_" -f1-3 | sort | uniq`
do
  bams=`echo ${id}_L00[1-4]_R1_001_bismark_bt2.bam`
  # echo $id, $bams
  sbatch -J $id -o $id.dedup.log --mem 8GB --wrap "deduplicate_bismark -s --output_dir . --bam --multiple $bams"
done
```

Sort and index deduplicated bam files:

```bash
cd ~/bismark

for bam in *.deduplicated.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o $bname.sorted.log --mem 8G --wrap "samtools sort -@ 20 -T ~/tmp/$bname -o $bname.sorted.bam $bam && \
  samtools index $bname.sorted.bam"
done
```

Check depth:

```bash
cd ~/bismark

for bam in *.deduplicated.sorted.bam
do
  echo $bam, `samtools depth $bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'`
done
```


#### Extract methylation

```bash
cd ~/bismark

mkdir ../methylation

ref=../reference

for bam in *.deduplicated.bam
do
  bname=${bam%_L001_R1_001_bismark_bt2.multiple.deduplicated.bam}
  # echo $bname, $ref, $bam
  sbatch -J $bname -o ../methylation/$bname.log --mem 8GB --wrap "bismark_methylation_extractor -s --comprehensive -o ../methylation --gzip --parallel 8 --bedGraph -CX_context --cytosine_report --genome_folder $ref $bam"
done
```

Uncompress resulting `*.CX_report.txt.gz` files:

```bash
cd ~/methylation
pigz -d *.CX_report.txt.gz
```


#### Sequence context plots

Obtain sequence context and methylation files:

```bash
srun --mem 32G --pty /usr/bin/bash

## Sequence context file
cd ~/reference

ref=Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa

fastaRegexFinder.py -f $ref -r C | \
cut -f2-3,6 | \
awk -v OFS="\t" '{print "Chromosome", $1, $2, ".", ".", $3}' | \
bedtools slop -i - -g $ref.fai -b 6 | \
bedtools getfasta -fi $ref -bed - -bedOut -s | \
awk 'length($7) == 13' | \
awk -v OFS="\t" '{print $1, $2+6, $3-6, $4, $5, $6, $7}' > Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.context.bed

wc -l Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.context.bed # 2357525 Cs with +/- 6bp sequence context in E.coli both in forward and reverse strands

## Methylation files
cd ~/methylation

for report in *.CX_report.txt
do
  bname=${report%.CX_report.txt.CX_report.txt}
  nohup awk -v OFS="\t" '{print $1, $2-1, $2, $4, $5, $3}' $report | \
  bedtools intersect -a - -b ../reference/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.context.bed -loj -s | \
  awk -v OFS="\t" '{if($7 ~ /^Chr/) print $1, $2, $3, $4, $5, $6, $13}' > $bname.context.txt &
done
```

Plotting (compact boxplots):

```r
library(data.table)
library(ggplot2)

# Set width
options(width = 300)

# Load data
data <- fread("tableCat.py -i ~/methylation/*.context.txt -r _L001_R1_001_bismark_bt2.multiple.deduplicated.context.txt")
setnames(data, c("ref", "start", "end", "cnt_met", "cnt_unmet", "strand", "context", "library"))

data[, library := sapply(data$library, function(x) unlist(strsplit(x, "_"))[1])]

# Calculate percentage methylation
data[, pct_met := round(100 * cnt_met/(cnt_met+cnt_unmet), 2)]

# Filter with cnt_met + cnt_unmet > 5
data <- data[cnt_met + cnt_unmet > 5]


#######
# XCY #
#######
data_xcy <- copy(data)

# Extract XCY context
data_xcy[, context_xcy := as.vector(sapply(data_xcy$context, function(x) paste(unlist(strsplit(x, ""))[6:8], collapse = "")))]

# load function for number of observations
give.n <- function(x){
  return(c(y = 102, label = length(x)))
}

# Plot for each library
for (l in unique(data_xcy$library)){
  print(l)
  # Normal compact boxplot
  gg <- ggplot(data = data_xcy[library == l & !is.na(pct_met)], aes(x = reorder(context_xcy, pct_met, FUN = median), y = pct_met)) +
  geom_boxplot(outlier.shape=NA) +
  ylab(expression("% Methylation")) +
  xlab(expression("")) +
  ggtitle(l) +
  theme_bw() +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
  coord_flip(ylim = c(0, 100))
  ggsave(paste("~/figures/", l, "_xcy_compact.pdf",sep = ""))
  # Numbered compact boxplot
  gg <- ggplot(data = data_xcy[library == l & !is.na(pct_met)], aes(x = reorder(context_xcy, pct_met, FUN = median), y = pct_met)) +
  geom_boxplot(outlier.shape=NA) +
  ylab(expression("% Methylation")) +
  xlab(expression("")) +
  ggtitle(l) +
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text", hjust = 0) +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
  coord_flip(ylim = c(0, 110)) +
  scale_y_continuous(breaks=c(0, 25, 50, 75, 100))
  ggsave(paste("~/figures/", l, "_xcy_compact_number.pdf",sep = ""))
}

rm(data_xcy)


########
# XCGY #
########
data_xcgy <- copy(data[grepl("......CG.....", data$context)])

# Extract XCGY context
data_xcgy[, context_xcgy := as.vector(sapply(data_xcgy$context, function(x) paste(unlist(strsplit(x, ""))[6:9], collapse = "")))]

# load function for number of observations
give.n <- function(x){
  return(c(y = 102, label = length(x)))
}

# Plot for each library
for (l in unique(data_xcgy$library)){
  print(l)
  # Normal compact boxplot
  gg <- ggplot(data = data_xcgy[library == l & !is.na(pct_met)], aes(x = reorder(context_xcgy, pct_met, FUN = median), y = pct_met)) +
  geom_boxplot(outlier.shape=NA) +
  ylab(expression("% Methylation")) +
  xlab(expression("")) +
  ggtitle(l) +
  theme_bw() +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
  coord_flip(ylim = c(0, 100))
  ggsave(paste("~/figures/", l, "_xcgy_compact.pdf", sep = ""))
  # Numbered compact boxplot
  gg <- ggplot(data = data_xcgy[library == l & !is.na(pct_met)], aes(x = reorder(context_xcgy, pct_met, FUN = median), y = pct_met)) +
  geom_boxplot(outlier.shape=NA) +
  ylab(expression("% Methylation")) +
  xlab(expression("")) +
  ggtitle(l) +
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text", hjust = 0) +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
  coord_flip(ylim = c(0, 110)) +
  scale_y_continuous(breaks=c(0, 25, 50, 75, 100))
  ggsave(paste("~/figures/", l, "_xcgy_compact_number.pdf", sep = ""))
}

rm(data_xcgy)


##########
# WXCGYZ #
##########
data_wxcgyz <- copy(data[grepl("......CG.....", data$context)])

# Extract WXCGYZ context
data_wxcgyz[, context_wxcgyz := as.vector(sapply(data_wxcgyz$context, function(x) paste(unlist(strsplit(x, ""))[5:10], collapse = "")))]

# load function for number of observations
give.n <- function(x){
  return(c(y = 102, label = length(x)))
}

# Plot for each library
for (l in unique(data_wxcgyz$library)){
  print(l)
  # Normal compact boxplot
  gg <- ggplot(data = data_wxcgyz[library == l & !is.na(pct_met)], aes(x = reorder(context_wxcgyz, pct_met, FUN = median), y = pct_met)) +
  geom_boxplot(outlier.shape=NA) +
  ylab(expression("% Methylation")) +
  xlab(expression("")) +
  ggtitle(l) +
  theme_bw() +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
  coord_flip(ylim = c(0, 100))
  ggsave(paste("~/figures/", l, "_wxcgyz_compact.pdf",sep = ""), height = 140, units= 'cm', limitsize = FALSE)
  # Numbered compact boxplot
  gg <- ggplot(data = data_wxcgyz[library == l & !is.na(pct_met)], aes(x = reorder(context_wxcgyz, pct_met, FUN = median), y = pct_met)) +
  geom_boxplot(outlier.shape=NA) +
  ylab(expression("% Methylation")) +
  xlab(expression("")) +
  ggtitle(l) +
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text", hjust = 0) +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
  coord_flip(ylim = c(0, 110)) +
  scale_y_continuous(breaks=c(0, 25, 50, 75, 100))
  ggsave(paste("~/figures/", l, "_wxcgyz_compact_number.pdf",sep = ""), height = 140, units= 'cm', limitsize = FALSE)
}

rm(data_wxcgyz)
```
