## Contents

- [Quality check](ecoli_bs.md#quality-check)
- [Trimming](ecoli_bs.md#trimming)
- [Prepare reference genome](ecoli_bs.md#prepare-reference-genome)
- [Alignment](ecoli_bs.md#alignment)
- [Deduplication](ecoli_bs.md#deduplication)
- [Extract methylation](ecoli_bs.md#extract-methylation)
- [Sequence context plots](ecoli_bs.md#sequence-context-plots)


### Quality check

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


### Trimming

```bash
cd ~/fastq

mkdir ../fastq_trimmed

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../fastq_trimmed/$bname.log --mem 4G --wrap "cutadapt -a AGATCGGAAGAGC -m 10 -u 6 -o ../fastq_trimmed/$fq $fq > ../fastq_trimmed/$bname.txt"
done
```


### Prepare reference genome

Genome and annotations for `Escherichia coli str. K-12 substr. MG1655` obtained from [Ensembl Bacteria release 42](http://bacteria.ensembl.org/Escherichia_coli_str_k_12_substr_mg1655/Info/Index) and located in `~/reference`

```bash
cd ~/reference
bismark_genome_preparation .
```


### Alignment

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


### Deduplication

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


### Extract methylation

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


### Sequence context plots

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


########
# XCAY #
########
data_xcay <- copy(data[grepl("......CA.....", data$context)])

# Extract XCAY context
data_xcay[, context_xcay := as.vector(sapply(data_xcay$context, function(x) paste(unlist(strsplit(x, ""))[6:9], collapse = "")))]

# load function for number of observations
give.n <- function(x){
  return(c(y = 102, label = length(x)))
}

# Plot for each library
for (l in unique(data_xcay$library)){
  print(l)
  # Normal compact boxplot
  gg <- ggplot(data = data_xcay[library == l & !is.na(pct_met)], aes(x = reorder(context_xcay, pct_met, FUN = median), y = pct_met)) +
  geom_boxplot(outlier.shape=NA) +
  ylab(expression("% Methylation")) +
  xlab(expression("")) +
  ggtitle(l) +
  theme_bw() +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
  coord_flip(ylim = c(0, 100))
  ggsave(paste("~/figures/", l, "_xcay_compact.pdf", sep = ""))
  # Numbered compact boxplot
  gg <- ggplot(data = data_xcay[library == l & !is.na(pct_met)], aes(x = reorder(context_xcay, pct_met, FUN = median), y = pct_met)) +
  geom_boxplot(outlier.shape=NA) +
  ylab(expression("% Methylation")) +
  xlab(expression("")) +
  ggtitle(l) +
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text", hjust = 0) +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
  coord_flip(ylim = c(0, 110)) +
  scale_y_continuous(breaks=c(0, 25, 50, 75, 100))
  ggsave(paste("~/figures/", l, "_xcay_compact_number.pdf", sep = ""))
}

rm(data_xcay)


#######
# CAY #
#######
data_cay <- copy(data[grepl("......CA.....", data$context)])

# Extract CAY context
data_cay[, context_cay := as.vector(sapply(data_cay$context, function(x) paste(unlist(strsplit(x, ""))[7:9], collapse = "")))]

# load function for number of observations
give.n <- function(x){
  return(c(y = 52, label = length(x)))
}

# Plot for each library
for (l in unique(data_cay$library)){
  print(l)
  # Normal compact boxplot
  gg <- ggplot(data = data_cay[library == l & !is.na(pct_met)], aes(x = reorder(context_cay, pct_met, FUN = median), y = pct_met)) +
  geom_boxplot(outlier.shape=NA) +
  ylab(expression("% Methylation")) +
  xlab(expression("")) +
  ggtitle(l) +
  theme_bw() +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
  coord_flip(ylim = c(0, 50))
  ggsave(paste("~/figures/", l, "_cay_compact.pdf", sep = ""), height = 4, width = 4)
  # Numbered compact boxplot
  gg <- ggplot(data = data_cay[library == l & !is.na(pct_met)], aes(x = reorder(context_cay, pct_met, FUN = median), y = pct_met)) +
  geom_boxplot(outlier.shape=NA) +
  ylab(expression("% Methylation")) +
  xlab(expression("")) +
  ggtitle(l) +
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text", hjust = 0) +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
  coord_flip(ylim = c(0, 60)) +
  scale_y_continuous(breaks=c(0, 25, 50, 75, 100))
  ggsave(paste("~/figures/", l, "_cay_compact_number.pdf", sep = ""), height = 4, width = 4)
}


###########
# CAC/CAG #
###########
cac_cag <- dcast(data_cay[, .(pct_met = median(pct_met, na.rm=TRUE)), by = .(library, context_cay)], library ~ context_cay, value.var = "pct_met")
cac_cag[, CAC_CAG := CAC/CAG]
cac_cag$library <- factor(cac_cag$library, levels=c("DNMT3A-Abcam-30min", "DNMT3A-Abcam-120min", "DNMT3A-Abcam-240min", "DNMT3B-Abcam-30min", "DNMT3B-Abcam-120min", "DNMT3B-Abcam-240min", "MSssI-10min", "MSssI-30min", "MSssI-240min", "Ecoli-unmethylated-DNA"))

gg <- ggplot(cac_cag, aes(x=library, y=CAC_CAG)) +
geom_bar(stat="identity", color="black", position=position_dodge(), alpha = 0.5) +
theme_classic() +
ylab(expression("[mCAC/CAC]/[mCAG/CAG]")) +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 12, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black"), plot.margin = margin(0.5, 0.5, 0.5, 2, "cm")) +
coord_cartesian(ylim = c(0, 1.5))

ggsave("~/figures/CAC_CAG.pdf", width = 5, height = 7)
```

Plotting (custom logos/heatmaps):

```r
#cd ~/methylation
#R

library(data.table)
library(ggplot2)
library(ggseqlogo)
library(grid)
library(gridExtra)

# Set width
options(width = 300)

# Load data
data <- fread("tableCat.py -i ~/methylation/*.context.txt -r _L001_R1_001_bismark_bt2.multiple.deduplicated.context.txt")
setnames(data, c("ref", "start", "end", "cnt_met", "cnt_unmet", "strand", "context", "library"))

data[, library := sapply(data$library, function(x) unlist(strsplit(x, "_"))[1])]

# Calculate rounded percentage methylation
data[, pct_met := round(100 * cnt_met/(cnt_met+cnt_unmet), 2)]

# Filter with cnt_met + cnt_unmet > 5
data <- data[cnt_met + cnt_unmet > 5]


###########
# -5_C_+5 #
###########
data_c <- copy(data)

# Extract -5_C_+5 context
data_c[, context_c := as.vector(sapply(data_c$context, function(x) paste(unlist(strsplit(x, ""))[2:12], collapse = "")))]

# Count distributions
gg <- ggplot(data_c[!is.na(pct_met)], aes(x=cnt_met+cnt_unmet)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("depth") +
ylab("frequency") +
facet_wrap(~ factor(library, levels = c("DNMT3A-Abcam-30min", "DNMT3A-Abcam-120min", "DNMT3A-Abcam-240min", "DNMT3B-Abcam-30min", "DNMT3B-Abcam-120min", "DNMT3B-Abcam-240min", "MSssI-10min", "MSssI-30min", "MSssI-240min", "Ecoli-unmethylated-DNA")), ncol = 3) +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"))

ggsave('~/figures/count_distribution_c.pdf', width = 24, height = 24, units = 'cm')

# % methylation distributions
gg <- ggplot(data_c[!is.na(pct_met) & (cnt_met+cnt_unmet) > 9], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency") +
facet_wrap(~ factor(library, levels = c("DNMT3A-Abcam-30min", "DNMT3A-Abcam-120min", "DNMT3A-Abcam-240min", "DNMT3B-Abcam-30min", "DNMT3B-Abcam-120min", "DNMT3B-Abcam-240min", "MSssI-10min", "MSssI-30min", "MSssI-240min", "Ecoli-unmethylated-DNA")), ncol = 3, scales = "free_y") +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black")) +
coord_cartesian(xlim = c(0, 100))

ggsave('~/figures/pctmet_distribution_c.pdf', width = 28, height = 24, units = 'cm')


# Plot for each library
for (l in unique(data_c$library)){
  print(l)
  # Obtain filtered and order matrix
  data_c_filter <- data_c[library == l & !is.na(pct_met) & (cnt_met+cnt_unmet) > 9][order(-pct_met)]
  # Obtain top and bottom sequences
  seqs_top <- data_c_filter[1:round(nrow(data_c_filter)/10)]$context_c
  seqs_bottom <- data_c_filter[(nrow(data_c_filter)-round(nrow(data_c_filter)/10)):nrow(data_c_filter)]$context_c
  # Plot logo
  print("- bits top 10%")
  gg_bits_top <- ggplot() +
  geom_logo(seqs_top, method = 'bits', rev_stack_order = T) +
  ggtitle("top 10%") +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 11, 1), labels =  seq(-5, 5, 1)) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  print("- prob top 10%")
  gg_prob_top <- ggplot() +
  geom_logo(seqs_top, method = 'prob', rev_stack_order = T) +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 11, 1), labels = seq(-5, 5, 1)) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  print("- bits bottom 10%")
  gg_bits_bottom <- ggplot() +
  geom_logo(seqs_bottom, method = 'bits', rev_stack_order = T) +
  ggtitle("bottom 10%") +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 11, 1), labels =  seq(-5, 5, 1)) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  print("- prob bottom 10%")
  gg_prob_bottom <- ggplot() +
  geom_logo(seqs_bottom, method = 'prob', rev_stack_order = T) +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 11, 1), labels = seq(-5, 5, 1)) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  gg <- arrangeGrob(gg_bits_top, gg_bits_bottom, gg_prob_top, gg_prob_bottom, ncol = 2, top = textGrob(l, gp=gpar(fontface="bold", fontsize=20)))
  ggsave(file = paste("~/figures/", l, "_ggseqlogo_c_topbottom10pct.pdf",sep = ""), gg, width = 10)
}

rm(data_c)


############
# -5_CG_+5 #
############
data_cg <- copy(data[grepl("......CG.....", data$context)])

# Extract -5_CG_+5 context
data_cg[, context_cg := as.vector(sapply(data_cg$context, function(x) paste(unlist(strsplit(x, ""))[2:13], collapse = "")))]

# Count distributions
gg <- ggplot(data_cg[!is.na(pct_met)], aes(x=cnt_met+cnt_unmet)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("depth") +
ylab("frequency") +
facet_wrap(~ factor(library, levels = c("DNMT3A-Abcam-30min", "DNMT3A-Abcam-120min", "DNMT3A-Abcam-240min", "DNMT3B-Abcam-30min", "DNMT3B-Abcam-120min", "DNMT3B-Abcam-240min", "MSssI-10min", "MSssI-30min", "MSssI-240min", "Ecoli-unmethylated-DNA")), ncol = 3) +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"))

ggsave('~/figures/count_distribution_cg.pdf', width = 24, height = 24, units = 'cm')

# % methylation distributions
gg <- ggplot(data_cg[!is.na(pct_met) & (cnt_met+cnt_unmet) > 9], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency") +
facet_wrap(~ factor(library, levels = c("DNMT3A-Abcam-30min", "DNMT3A-Abcam-120min", "DNMT3A-Abcam-240min", "DNMT3B-Abcam-30min", "DNMT3B-Abcam-120min", "DNMT3B-Abcam-240min", "MSssI-10min", "MSssI-30min", "MSssI-240min", "Ecoli-unmethylated-DNA")), ncol = 3, scales = "free_y") +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black")) +
coord_cartesian(xlim = c(0, 100))

ggsave('~/figures/pctmet_distribution_cg.pdf', width = 28, height = 24, units = 'cm')


# Plot for each library
for (l in unique(data_cg$library)){
  print(l)
  # Obtain filtered and order matrix
  data_cg_filter <- data_cg[library == l & !is.na(pct_met) & (cnt_met+cnt_unmet) > 9][order(-pct_met)]
  # Obtain top and bottom sequences
  seqs_top <- data_cg_filter[1:round(nrow(data_cg_filter)/10)]$context_cg
  seqs_bottom <- data_cg_filter[(nrow(data_cg_filter)-round(nrow(data_cg_filter)/10)):nrow(data_cg_filter)]$context_cg
  # Plot logo
  print("- bits top 10%")
  gg_bits_top <- ggplot() +
  geom_logo(seqs_top, method = 'bits', rev_stack_order = T) +
  ggtitle("top 10%") +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  print("- prob top 10%")
  gg_prob_top <- ggplot() +
  geom_logo(seqs_top, method = 'prob', rev_stack_order = T) +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  print("- bits bottom 10%")
  gg_bits_bottom <- ggplot() +
  geom_logo(seqs_bottom, method = 'bits', rev_stack_order = T) +
  ggtitle("bottom 10%") +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  print("- prob bottom 10%")
  gg_prob_bottom <- ggplot() +
  geom_logo(seqs_bottom, method = 'prob', rev_stack_order = T) +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  gg <- arrangeGrob(gg_bits_top, gg_bits_bottom, gg_prob_top, gg_prob_bottom, ncol = 2, top = textGrob(l, gp=gpar(fontface="bold", fontsize=20)))
  ggsave(file = paste("~/figures/", l, "_ggseqlogo_cg_topbottom10pct.pdf",sep = ""), gg, width = 10)
}

rm(data_cg)


############
# -5_CA_+5 #
############
data_ca <- copy(data[grepl("......CA.....", data$context)])

# Extract -5_CA_+5 context
data_ca[, context_ca := as.vector(sapply(data_ca$context, function(x) paste(unlist(strsplit(x, ""))[2:13], collapse = "")))]

# Count distributions
gg <- ggplot(data_ca[!is.na(pct_met)], aes(x=cnt_met+cnt_unmet)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("depth") +
ylab("frequency") +
facet_wrap(~ factor(library, levels = c("DNMT3A-Abcam-30min", "DNMT3A-Abcam-120min", "DNMT3A-Abcam-240min", "DNMT3B-Abcam-30min", "DNMT3B-Abcam-120min", "DNMT3B-Abcam-240min", "MSssI-10min", "MSssI-30min", "MSssI-240min", "Ecoli-unmethylated-DNA")), ncol = 3) +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"))

ggsave('~/figures/count_distribution_ca.pdf', width = 24, height = 24, units = 'cm')

# % methylation distributions
gg <- ggplot(data_ca[!is.na(pct_met) & (cnt_met+cnt_unmet) > 9], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency") +
facet_wrap(~ factor(library, levels = c("DNMT3A-Abcam-30min", "DNMT3A-Abcam-120min", "DNMT3A-Abcam-240min", "DNMT3B-Abcam-30min", "DNMT3B-Abcam-120min", "DNMT3B-Abcam-240min", "MSssI-10min", "MSssI-30min", "MSssI-240min", "Ecoli-unmethylated-DNA")), ncol = 3, scales = "free_y") +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black")) +
coord_cartesian(xlim = c(0, 100))

ggsave('~/figures/pctmet_distribution_ca.pdf', width = 28, height = 24, units = 'cm')


# Plot for each library
for (l in unique(data_ca$library)){
  print(l)
  # Obtain filtered and order matrix
  data_ca_filter <- data_ca[library == l & !is.na(pct_met) & (cnt_met+cnt_unmet) > 9][order(-pct_met)]
  # Obtain top and bottom sequences
  seqs_top <- data_ca_filter[1:round(nrow(data_ca_filter)/10)]$context_ca
  seqs_bottom <- data_ca_filter[(nrow(data_ca_filter)-round(nrow(data_ca_filter)/10)):nrow(data_ca_filter)]$context_ca
  # Plot logo
  print("- bits top 10%")
  gg_bits_top <- ggplot() +
  geom_logo(seqs_top, method = 'bits', rev_stack_order = T) +
  ggtitle("top 10%") +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  print("- prob top 10%")
  gg_prob_top <- ggplot() +
  geom_logo(seqs_top, method = 'prob', rev_stack_order = T) +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  print("- bits bottom 10%")
  gg_bits_bottom <- ggplot() +
  geom_logo(seqs_bottom, method = 'bits', rev_stack_order = T) +
  ggtitle("bottom 10%") +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  print("- prob bottom 10%")
  gg_prob_bottom <- ggplot() +
  geom_logo(seqs_bottom, method = 'prob', rev_stack_order = T) +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  gg <- arrangeGrob(gg_bits_top, gg_bits_bottom, gg_prob_top, gg_prob_bottom, ncol = 2, top = textGrob(l, gp=gpar(fontface="bold", fontsize=20)))
  ggsave(file = paste("~/figures/", l, "_ggseqlogo_ca_topbottom10pct.pdf",sep = ""), gg, width = 10)
}

rm(data_ca)


############
# -5_CC_+5 #
############
data_cc <- copy(data[grepl("......CC.....", data$context)])

# Extract -5_CC_+5 context
data_cc[, context_cc := as.vector(sapply(data_cc$context, function(x) paste(unlist(strsplit(x, ""))[2:13], collapse = "")))]

# Count distributions
gg <- ggplot(data_cc[!is.na(pct_met)], aes(x=cnt_met+cnt_unmet)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("depth") +
ylab("frequency") +
facet_wrap(~ factor(library, levels = c("DNMT3A-Abcam-30min", "DNMT3A-Abcam-120min", "DNMT3A-Abcam-240min", "DNMT3B-Abcam-30min", "DNMT3B-Abcam-120min", "DNMT3B-Abcam-240min", "MSssI-10min", "MSssI-30min", "MSssI-240min", "Ecoli-unmethylated-DNA")), ncol = 3) +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"))

ggsave('~/figures/count_distribution_cc.pdf', width = 24, height = 24, units = 'cm')

# % methylation distributions
gg <- ggplot(data_cc[!is.na(pct_met) & (cnt_met+cnt_unmet) > 9], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency") +
facet_wrap(~ factor(library, levels = c("DNMT3A-Abcam-30min", "DNMT3A-Abcam-120min", "DNMT3A-Abcam-240min", "DNMT3B-Abcam-30min", "DNMT3B-Abcam-120min", "DNMT3B-Abcam-240min", "MSssI-10min", "MSssI-30min", "MSssI-240min", "Ecoli-unmethylated-DNA")), ncol = 3, scales = "free_y") +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black")) +
coord_cartesian(xlim = c(0, 100))

ggsave('~/figures/pctmet_distribution_cc.pdf', width = 28, height = 24, units = 'cm')


# Plot for each library
for (l in unique(data_cc$library)){
  print(l)
  # Obtain filtered and order matrix
  data_cc_filter <- data_cc[library == l & !is.na(pct_met) & (cnt_met+cnt_unmet) > 9][order(-pct_met)]
  # Obtain top and bottom sequences
  seqs_top <- data_cc_filter[1:round(nrow(data_cc_filter)/10)]$context_cc
  seqs_bottom <- data_cc_filter[(nrow(data_cc_filter)-round(nrow(data_cc_filter)/10)):nrow(data_cc_filter)]$context_cc
  # Plot logo
  print("- bits top 10%")
  gg_bits_top <- ggplot() +
  geom_logo(seqs_top, method = 'bits', rev_stack_order = T) +
  ggtitle("top 10%") +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  print("- prob top 10%")
  gg_prob_top <- ggplot() +
  geom_logo(seqs_top, method = 'prob', rev_stack_order = T) +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  print("- bits bottom 10%")
  gg_bits_bottom <- ggplot() +
  geom_logo(seqs_bottom, method = 'bits', rev_stack_order = T) +
  ggtitle("bottom 10%") +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  print("- prob bottom 10%")
  gg_prob_bottom <- ggplot() +
  geom_logo(seqs_bottom, method = 'prob', rev_stack_order = T) +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  gg <- arrangeGrob(gg_bits_top, gg_bits_bottom, gg_prob_top, gg_prob_bottom, ncol = 2, top = textGrob(l, gp=gpar(fontface="bold", fontsize=20)))
  ggsave(file = paste("~/figures/", l, "_ggseqlogo_cc_topbottom10pct.pdf",sep = ""), gg, width = 10)
}

rm(data_cc)


############
# -5_CT_+5 #
############
data_ct <- copy(data[grepl("......CT.....", data$context)])

# Extract -5_CT_+5 context
data_ct[, context_ct := as.vector(sapply(data_ct$context, function(x) paste(unlist(strsplit(x, ""))[2:13], collapse = "")))]

# Count distributions
gg <- ggplot(data_ct[!is.na(pct_met)], aes(x=cnt_met+cnt_unmet)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("depth") +
ylab("frequency") +
facet_wrap(~ factor(library, levels = c("DNMT3A-Abcam-30min", "DNMT3A-Abcam-120min", "DNMT3A-Abcam-240min", "DNMT3B-Abcam-30min", "DNMT3B-Abcam-120min", "DNMT3B-Abcam-240min", "MSssI-10min", "MSssI-30min", "MSssI-240min", "Ecoli-unmethylated-DNA")), ncol = 3) +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"))

ggsave('~/figures/count_distribution_ct.pdf', width = 24, height = 24, units = 'cm')

# % methylation distributions
gg <- ggplot(data_ct[!is.na(pct_met) & (cnt_met+cnt_unmet) > 9], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency") +
facet_wrap(~ factor(library, levels = c("DNMT3A-Abcam-30min", "DNMT3A-Abcam-120min", "DNMT3A-Abcam-240min", "DNMT3B-Abcam-30min", "DNMT3B-Abcam-120min", "DNMT3B-Abcam-240min", "MSssI-10min", "MSssI-30min", "MSssI-240min", "Ecoli-unmethylated-DNA")), ncol = 3, scales = "free_y") +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black")) +
coord_cartesian(xlim = c(0, 100))

ggsave('~/figures/pctmet_distribution_ct.pdf', width = 28, height = 24, units = 'cm')


# Plot for each library
for (l in unique(data_ct$library)){
  print(l)
  # Obtain filtered and order matrix
  data_ct_filter <- data_ct[library == l & !is.na(pct_met) & (cnt_met+cnt_unmet) > 9][order(-pct_met)]
  # Obtain top and bottom sequences
  seqs_top <- data_ct_filter[1:round(nrow(data_ct_filter)/10)]$context_ct
  seqs_bottom <- data_ct_filter[(nrow(data_ct_filter)-round(nrow(data_ct_filter)/10)):nrow(data_ct_filter)]$context_ct
  # Plot logo
  print("- bits top 10%")
  gg_bits_top <- ggplot() +
  geom_logo(seqs_top, method = 'bits', rev_stack_order = T) +
  ggtitle("top 10%") +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  print("- prob top 10%")
  gg_prob_top <- ggplot() +
  geom_logo(seqs_top, method = 'prob', rev_stack_order = T) +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  print("- bits bottom 10%")
  gg_bits_bottom <- ggplot() +
  geom_logo(seqs_bottom, method = 'bits', rev_stack_order = T) +
  ggtitle("bottom 10%") +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  print("- prob bottom 10%")
  gg_prob_bottom <- ggplot() +
  geom_logo(seqs_bottom, method = 'prob', rev_stack_order = T) +
  theme_logo() +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = c(seq(-5, 0, 1), seq(0, 5, 1))) +
  theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
  gg <- arrangeGrob(gg_bits_top, gg_bits_bottom, gg_prob_top, gg_prob_bottom, ncol = 2, top = textGrob(l, gp=gpar(fontface="bold", fontsize=20)))
  ggsave(file = paste("~/figures/", l, "_ggseqlogo_ct_topbottom10pct.pdf",sep = ""), gg, width = 10)
}

rm(data_ct)
```
