This repository contains data access, source code and a small demo for the custom computational analyses developed as part of our paper **Sequence preference of DNMT3A/B dictates genome-wide DNA methylation signature** *(under review)*


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
- [Instructions for use](README.md#instructions-for-use)
  - [Sequence context diversity in reference genomes](README.md#sequence-context-diversity-in-reference-genomes)
  - [Quality check](README.md#quality-check)
  - [Trimming](README.md#trimming)
  - [Prepare reference genome](README.md#prepare-reference-genome)
  - [Alignment](README.md#alignment)
  - [Deduplication](README.md#deduplication)
  - [Extract methylation](README.md#extract-methylation)
  - [Sequence context plots](README.md#sequence-context-plots)
  
- [Small demo](README.md#small-demo)


### System requirements and installation

Software, for installation details for the individual tools follow the links:

- [EMBOSS v6.6.0.0](http://emboss.sourceforge.net/)
- [fastaRegexFinder.py v0.1.1](https://github.com/dariober/bioinformatics-cafe/tree/master/fastaRegexFinder)
- Standard Unix tools: cut, awk, sort, uniq, grep, tr, echo ...
- [bedtools v2.27.0](http://bedtools.readthedocs.io/en/latest/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [bismark v0.19.0](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
- [R v3.3.2](https://www.r-project.org/). Libraries:
  - [data.table v1.10.4](https://cran.r-project.org/web/packages/data.table/index.html)
  - [ggplot2 v2.2.1](http://ggplot2.org/)

Operating system:

  - CentOS Linux release 7.3.1611 (OS used during code development)
  - [slurm job scheduling system v19.05.0](https://slurm.schedmd.com/quickstart.html)


### Instructions for use

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

`*.fastq.gz` files to be downloaded from [E-MTAB-8406](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8406) and located in `~/fastq` folder:

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
data <- fread("tableCat.py -i /scratchb/sblab/martin03/repository/20181213_DNMT_preference/data/20190326/methylation/*.context.txt -r _L001_R1_001_bismark_bt2.multiple.deduplicated.context.txt")
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

# Explore XCY context for e.g. library DNMT3B-Abcam-240min
data_xcy[library == "DNMT3B-Abcam-240min", .(.N, pct_met = median(pct_met, na.rm=TRUE)), by = .(context_xcy)][order(-pct_met)]
#    context_xcy      N pct_met
# 1:         ACG 141173   57.14
# 2:         TCG 138123   50.00
# 3:         GCG 223988   45.00
# 4:         CCG 167518   37.50
# 5:         ACA 111348   15.00
# 6:         GCA 184419    7.69
# 7:         ACT  95133    7.14
# 8:         GCT 155339    6.67
# 9:         TCA 160172    6.25
#10:         CCA 165274    4.17
#11:         ACC 143463    4.17
#12:         TCT 107414    3.45
#13:         GCC 178590    2.86
#14:         CCC  90656    0.00
#15:         TCC 107317    0.00
#16:         CCT  96499    0.00


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
  ggsave(paste("/Users/martin03/github/sblab-bioinformatics/projects/20181213_DNMT_preference/figures/20190329_", l, "_xcy_compact.pdf",sep = ""))
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
  ggsave(paste("/Users/martin03/github/sblab-bioinformatics/projects/20181213_DNMT_preference/figures/20190329_", l, "_xcy_compact_number.pdf",sep = ""))
}

rm(data_xcy)


########
# XCGY #
########
data_xcgy <- copy(data[grepl("......CG.....", data$context)])

# Extract XCGY context
data_xcgy[, context_xcgy := as.vector(sapply(data_xcgy$context, function(x) paste(unlist(strsplit(x, ""))[6:9], collapse = "")))]


# Explore XCGY context for e.g. libraries DNMT3A-Abcam-240min and DNMT3B-Abcam-240min
## DNMT3A-Abcam-240min
data_xcgy[library == "DNMT3A-Abcam-240min", .(.N, pct_met = median(pct_met, na.rm=TRUE)), by = .(context_xcgy)][order(-pct_met)]
#    context_xcgy     N pct_met
# 1:         CCGC 55726   46.15
# 2:         ACGC 51048   45.45
# 3:         TCGC 49777   45.45
# 4:         GCGC 68921   40.00
# 5:         ACGT 28263   37.50
# 6:         CCGT 35097   36.36
# 7:         TCGT 28071   35.29
# 8:         CCGA 31141   34.48
# 9:         CCGG 46820   32.43
#10:         GCGT 51058   32.00
#11:         ACGA 28064   31.03
#12:         ACGG 35123   29.41
#13:         TCGG 31161   29.17
#14:         GCGA 49775   28.00
#15:         TCGA 30174   25.93
#16:         GCGG 55758   25.00

## DNMT3B-Abcam-240min
data_xcgy[library == "DNMT3B-Abcam-240min", .(.N, pct_met = median(pct_met, na.rm=TRUE)), by = .(context_xcgy)][order(-pct_met)]
#    context_xcgy     N pct_met
# 1:         TCGG 30986   67.50
# 2:         ACGG 34830   66.67
# 3:         TCGA 29912   62.50
# 4:         ACGA 27759   61.54
# 5:         GCGG 55410   57.14
# 6:         CCGG 46522   56.00
# 7:         CCGA 30850   55.56
# 8:         GCGA 49406   54.55
# 9:         ACGC 50560   53.85
#10:         ACGT 28024   45.00
#11:         TCGC 49357   42.86
#12:         GCGC 68432   36.11
#13:         GCGT 50740   33.33
#14:         TCGT 27868   30.00
#15:         CCGC 55270   23.81
#16:         CCGT 34876   18.18


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
  ggsave(paste("/Users/martin03/github/sblab-bioinformatics/projects/20181213_DNMT_preference/figures/20190329_", l, "_xcgy_compact.pdf", sep = ""))
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
  ggsave(paste("/Users/martin03/github/sblab-bioinformatics/projects/20181213_DNMT_preference/figures/20190329_", l, "_xcgy_compact_number.pdf", sep = ""))
}

rm(data_xcgy)


##########
# WXCGYZ #
##########
data_wxcgyz <- copy(data[grepl("......CG.....", data$context)])

# Extract WXCGYZ context
data_wxcgyz[, context_wxcgyz := as.vector(sapply(data_wxcgyz$context, function(x) paste(unlist(strsplit(x, ""))[5:10], collapse = "")))]


# Explore WXCGYZ context for e.g. libraries DNMT3A-Abcam-240min and DNMT3B-Abcam-240min
## DNMT3A-Abcam-240min
data_wxcgyz[library == "DNMT3A-Abcam-240min", .(.N, pct_met = median(pct_met, na.rm=TRUE)), by = .(context_wxcgyz)][order(-pct_met)]
#     context_wxcgyz    N pct_met
#  1:         TACGCC 3237   66.67
#  2:         TACGTC 1237   65.71
#  3:         TCCGCC 3212   65.22
#  4:         TCCGTC 1557   62.50
#  5:         CACGCC 3964   62.50
# ---                            
#252:         GGCGTG 3961   15.15
#253:         GGCGGG 2991   15.00
#254:         AGCGTG 2956   14.71
#255:         GTCGGG 1559   14.29
#256:         AGCGGG 2595   12.90

## DNMT3B-Abcam-240min
data_wxcgyz[library == "DNMT3B-Abcam-240min", .(.N, pct_met = median(pct_met, na.rm=TRUE)), by = .(context_wxcgyz)][order(-pct_met)]
#     context_wxcgyz    N pct_met
#  1:         GTCGGC 2723   73.91
#  2:         ATCGGC 4304   73.68
#  3:         TACGGC 2376   73.33
#  4:         TGCGGC 4289   73.08
#  5:         TTCGGC 3150   72.73
# ---                            
#252:         GCCGCA 4275   11.11
#253:         GCCGTT 4190   11.11
#254:         GCCGCG 4077   10.34
#255:         ACCGTG 2522    9.09
#256:         GCCGTG 2635    8.33


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
  ggsave(paste("/Users/martin03/github/sblab-bioinformatics/projects/20181213_DNMT_preference/figures/20190329_", l, "_wxcgyz_compact.pdf",sep = ""), height = 140, units= 'cm', limitsize = FALSE)
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
  ggsave(paste("/Users/martin03/github/sblab-bioinformatics/projects/20181213_DNMT_preference/figures/20190329_", l, "_wxcgyz_compact_number.pdf",sep = ""), height = 140, units= 'cm', limitsize = FALSE)
}

rm(data_wxcgyz)


############
# UWXCGYZV #
############
data_uwxcgyzv <- copy(data[grepl("......CG.....", data$context)])

# Extract UWXCGYZV context
data_uwxcgyzv[, context_uwxcgyzv := as.vector(sapply(data_uwxcgyzv$context, function(x) paste(unlist(strsplit(x, ""))[4:11], collapse = "")))]


# Explore UWXCGYZV context for e.g. libraries DNMT3A-Abcam-240min and DNMT3B-Abcam-240min
## DNMT3A-Abcam-240min
data_uwxcgyzv[library == "DNMT3A-Abcam-240min", .(.N, pct_met = median(pct_met, na.rm=TRUE)), by = .(context_uwxcgyzv)][order(-pct_met)]
#      context_uwxcgyzv   N pct_met
#   1:         GTACGTCA  99  71.880
#   2:         CTACGTCA  85  70.370
#   3:         TTACGTCA 170  69.440
#   4:         GTCCGCCA 203  68.750
#   5:         TTACGCCA 382  68.750
#  ---                             
#4092:         AATCGGAG  42   7.155
#4093:         TAGCGGAG  46   6.965
#4094:         AATCGGGG  66   6.905
#4095:         AATCGGGT 128   6.670
#4096:         TAGCGGGG  53   6.450


## DNMT3B-Abcam-240min
data_uwxcgyzv[library == "DNMT3B-Abcam-240min", .(.N, pct_met = median(pct_met, na.rm=TRUE)), by = .(context_uwxcgyzv)][order(-pct_met)]
#      context_uwxcgyzv   N pct_met
#   1:         GTACGGGG  40  76.920
#   2:         AGTCGGCT  70  76.235
#   3:         TTACGGCG 221  76.190
#   4:         ATGCGGCG 538  76.000
#   5:         CATCGGCG 424  75.930
#  ---                             
#4092:         GGCCGCGA 132   5.260
#4093:         GGCCGTGG 110   5.130
#4094:         AACCGCGA 219   5.000
#4095:         TACCGTGG 202   4.760
#4096:         AACCGTGG 244   3.850


# Plot for each library
for (l in unique(data_uwxcgyzv$library)){
  print(l)
  # Normal compact boxplot
  gg <- ggplot(data = data_uwxcgyzv[library == l & !is.na(pct_met)], aes(x = reorder(context_uwxcgyzv, pct_met, FUN = median), y = pct_met)) +
  geom_boxplot(outlier.shape=NA) +
  ylab(expression("% Methylation")) +
  xlab(expression("")) +
  ggtitle(l) +
  theme_bw() +
  theme(axis.title = element_text(size=16), axis.text.x = element_text(size=16, color = "black"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text = element_text(size=16, color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
  coord_flip(ylim = c(0, 100))
  ggsave(paste("/Users/martin03/github/sblab-bioinformatics/projects/20181213_DNMT_preference/figures/20190329_", l, "_uwxcgyzv_compact.pdf",sep = ""), height = 140, units= 'cm', limitsize = FALSE)
}

rm(data_uwxcgyzv)
```

Normal compact boxplots:

XCY | XCGY | WXCGYZ | UWXCGYZV
----|------|--------|---------
[DNMT3A-Abcam-30min](figures/20190329_DNMT3A-Abcam-30min_xcy_compact.pdf) | [DNMT3A-Abcam-30min](figures/20190329_DNMT3A-Abcam-30min_xcgy_compact.pdf) | [DNMT3A-Abcam-30min](figures/20190329_DNMT3A-Abcam-30min_wxcgyz_compact.pdf) | [DNMT3A-Abcam-30min](figures/20190329_DNMT3A-Abcam-30min_uwxcgyzv_compact.pdf)
[DNMT3A-Abcam-120min](figures/20190329_DNMT3A-Abcam-120min_xcy_compact.pdf) | [DNMT3A-Abcam-120min](figures/20190329_DNMT3A-Abcam-120min_xcgy_compact.pdf) | [DNMT3A-Abcam-120min](figures/20190329_DNMT3A-Abcam-120min_wxcgyz_compact.pdf) | [DNMT3A-Abcam-120min](figures/20190329_DNMT3A-Abcam-120min_uwxcgyzv_compact.pdf)
[DNMT3A-Abcam-240min](figures/20190329_DNMT3A-Abcam-240min_xcy_compact.pdf) | [DNMT3A-Abcam-240min](figures/20190329_DNMT3A-Abcam-240min_xcgy_compact.pdf) | [DNMT3A-Abcam-240min](figures/20190329_DNMT3A-Abcam-240min_wxcgyz_compact.pdf) | [DNMT3A-Abcam-240min](figures/20190329_DNMT3A-Abcam-240min_uwxcgyzv_compact.pdf)
[DNMT3B-Abcam-30min](figures/20190329_DNMT3B-Abcam-30min_xcy_compact.pdf) | [DNMT3B-Abcam-30min](figures/20190329_DNMT3B-Abcam-30min_xcgy_compact.pdf) | [DNMT3B-Abcam-30min](figures/20190329_DNMT3B-Abcam-30min_wxcgyz_compact.pdf) | [DNMT3B-Abcam-30min](figures/20190329_DNMT3B-Abcam-30min_uwxcgyzv_compact.pdf)
[DNMT3B-Abcam-120min](figures/20190329_DNMT3B-Abcam-120min_xcy_compact.pdf) | [DNMT3B-Abcam-120min](figures/20190329_DNMT3B-Abcam-120min_xcgy_compact.pdf) | [DNMT3B-Abcam-120min](figures/20190329_DNMT3B-Abcam-120min_wxcgyz_compact.pdf) | [DNMT3B-Abcam-120min](figures/20190329_DNMT3B-Abcam-120min_uwxcgyzv_compact.pdf)
[DNMT3B-Abcam-240min](figures/20190329_DNMT3B-Abcam-240min_xcy_compact.pdf) | [DNMT3B-Abcam-240min](figures/20190329_DNMT3B-Abcam-240min_xcgy_compact.pdf) | [DNMT3B-Abcam-240min](figures/20190329_DNMT3B-Abcam-240min_wxcgyz_compact.pdf) | [DNMT3B-Abcam-240min](figures/20190329_DNMT3B-Abcam-240min_uwxcgyzv_compact.pdf)
[MSssI-10min](figures/20190329_MSssI-10min_xcy_compact.pdf) | [MSssI-10min](figures/20190329_MSssI-10min_xcgy_compact.pdf) | [MSssI-10min](figures/20190329_MSssI-10min_wxcgyz_compact.pdf) | [MSssI-10min](figures/20190329_MSssI-10min_uwxcgyzv_compact.pdf)
[MSssI-30min](figures/20190329_MSssI-30min_xcy_compact.pdf) | [MSssI-30min](figures/20190329_MSssI-30min_xcgy_compact.pdf) | [MSssI-30min](figures/20190329_MSssI-30min_wxcgyz_compact.pdf) | [MSssI-30min](figures/20190329_MSssI-30min_uwxcgyzv_compact.pdf)
[MSssI-240min](figures/20190329_MSssI-240min_xcy_compact.pdf) | [MSssI-240min](figures/20190329_MSssI-240min_xcgy_compact.pdf) | [MSssI-240min](figures/20190329_MSssI-240min_wxcgyz_compact.pdf) | [MSssI-240min](figures/20190329_MSssI-240min_uwxcgyzv_compact.pdf)

Numbered compact boxplots:

XCY | XCGY | WXCGYZ
----|------|-------
[DNMT3A-Abcam-30min](figures/20190329_DNMT3A-Abcam-30min_xcy_compact_number.pdf) | [DNMT3A-Abcam-30min](figures/20190329_DNMT3A-Abcam-30min_xcgy_compact_number.pdf) | [DNMT3A-Abcam-30min](figures/20190329_DNMT3A-Abcam-30min_wxcgyz_compact_number.pdf)
[DNMT3A-Abcam-120min](figures/20190329_DNMT3A-Abcam-120min_xcy_compact_number.pdf) | [DNMT3A-Abcam-120min](figures/20190329_DNMT3A-Abcam-120min_xcgy_compact_number.pdf) | [DNMT3A-Abcam-120min](figures/20190329_DNMT3A-Abcam-120min_wxcgyz_compact_number.pdf)
[DNMT3A-Abcam-240min](figures/20190329_DNMT3A-Abcam-240min_xcy_compact_number.pdf) | [DNMT3A-Abcam-240min](figures/20190329_DNMT3A-Abcam-240min_xcgy_compact_number.pdf) | [DNMT3A-Abcam-240min](figures/20190329_DNMT3A-Abcam-240min_wxcgyz_compact_number.pdf)
[DNMT3B-Abcam-30min](figures/20190329_DNMT3B-Abcam-30min_xcy_compact_number.pdf) | [DNMT3B-Abcam-30min](figures/20190329_DNMT3B-Abcam-30min_xcgy_compact_number.pdf) | [DNMT3B-Abcam-30min](figures/20190329_DNMT3B-Abcam-30min_wxcgyz_compact_number.pdf)
[DNMT3B-Abcam-120min](figures/20190329_DNMT3B-Abcam-120min_xcy_compact_number.pdf) | [DNMT3B-Abcam-120min](figures/20190329_DNMT3B-Abcam-120min_xcgy_compact_number.pdf) | [DNMT3B-Abcam-120min](figures/20190329_DNMT3B-Abcam-120min_wxcgyz_compact_number.pdf)
[DNMT3B-Abcam-240min](figures/20190329_DNMT3B-Abcam-240min_xcy_compact_number.pdf) | [DNMT3B-Abcam-240min](figures/20190329_DNMT3B-Abcam-240min_xcgy_compact_number.pdf) | [DNMT3B-Abcam-240min](figures/20190329_DNMT3B-Abcam-240min_wxcgyz_compact_number.pdf)
[MSssI-10min](figures/20190329_MSssI-10min_xcy_compact_number.pdf) | [MSssI-10min](figures/20190329_MSssI-10min_xcgy_compact_number.pdf) | [MSssI-10min](figures/20190329_MSssI-10min_wxcgyz_compact_number.pdf)
[MSssI-30min](figures/20190329_MSssI-30min_xcy_compact_number.pdf) | [MSssI-30min](figures/20190329_MSssI-30min_xcgy_compact_number.pdf) | [MSssI-30min](figures/20190329_MSssI-30min_wxcgyz_compact_number.pdf)
[MSssI-240min](figures/20190329_MSssI-240min_xcy_compact_number.pdf) | [MSssI-240min](figures/20190329_MSssI-240min_xcgy_compact_number.pdf) | [MSssI-240min](figures/20190329_MSssI-240min_wxcgyz_compact_number.pdf)



