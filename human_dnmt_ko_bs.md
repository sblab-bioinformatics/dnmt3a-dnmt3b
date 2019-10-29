
WT:

- [HUES64 WGBS WT replicate 1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1112840)
  - SRX259085: SRR1067557 SRR1067560 SRR1067564 SRR1067567 SRR1067578
- [HUES64 WGBS WT replicate 2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1112841)
  - SRX259086: SRR1067565 SRR1067563 SRR1067562

DNMT KO:

- [HUES64 WGBS DNMT3A KO Early](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1545002)
  - SRX759479: SRR1652272
- [HUES64 WGBS DNMT3A KO Late](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1545005)
  - SRX759482: SRR1652275
- [HUES64 WGBS DNMT3B KO Early](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1545003)
  - SRX759480: SRR1652273
- [HUES64 WGBS DNMT3B KO Late](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1545006)
  - SRX759483: SRR1652276
- [HUES64 WGBS DNMT3A/3B Double KO Early](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1545004)
  - SRX759481: SRR1652274
- [HUES64 WGBS DNMT3A/3B Double KO Late](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1545007)
  - SRX759484: SRR1652277



## Contents

- [Obtain fastq files](human_dnmt_ko_bs.md#obtain-fastq-files)
- [Renaming](human_dnmt_ko_bs.md#renaming)
- [Quality check](human_dnmt_ko_bs.md#quality-check)
- [Trimming](human_dnmt_ko_bs.md#trimming)
- [Prepare reference genome](human_dnmt_ko_bs.md#prepare-reference-genome)
- [Alignment](human_dnmt_ko_bs.md#alignment)
- [Deduplication](human_dnmt_ko_bs.md#deduplication)
- [Sort and index deduplicated bam files, extract chr1, index and sort by name](human_dnmt_ko_bs.md#sort-and-index-deduplicated-bam-files-extract-chr1-index-and-sort-by-name)
- [Methylation chr1](human_dnmt_ko_bs.md#methylation-chr1)
- [Sequence context plots chr1](human_dnmt_ko_bs.md#sequence-context-plots-chr1)



### Obtain fastq files

```bash
cd ~

mkdir -p fastq && cd fastq

# WT
for i in SRR1067557 SRR1067560 SRR1067564 SRR1067567 SRR1067578 SRR1067565 SRR1067563 SRR1067562
do
  sbatch -J $i -o $i.log --mem 16G --wrap "fastq-dump --outdir . --gzip --split-files -F $i"
done

# KO
for i in SRR1652272 SRR1652275 SRR1652273 SRR1652276 SRR1652274 SRR1652277
do
  sbatch -J $i -o $i.log --mem 16G --wrap "fastq-dump --outdir . --gzip --split-files -F $i"
done

tail *.log # fine
```


### Renaming

- SRR1067557: HUES64_WT_rep1_1
- SRR1067560: HUES64_WT_rep1_2
- SRR1067564: HUES64_WT_rep1_3
- SRR1067567: HUES64_WT_rep1_4
- SRR1067578: HUES64_WT_rep1_5
- SRR1067565: HUES64_WT_rep2_1
- SRR1067563: HUES64_WT_rep2_2
- SRR1067562: HUES64_WT_rep2_3
- SRR1652272: HUES64_DNMT3AKOEarly_rep1_1
- SRR1652275: HUES64_DNMT3AKOLate_rep1_1
- SRR1652273: HUES64_DNMT3BKOEarly_rep1_1
- SRR1652276: HUES64_DNMT3BKOLate_rep1_1
- SRR1652274: HUES64_DNMT3A3BDoubleKOEarly_rep1_1
- SRR1652277: HUES64_DNMT3A3BDoubleKOLate_rep1_1

```bash
cd ~/fastq

# SRR1067557: HUES64_WT_rep1_1
mv SRR1067557_1.fastq.gz HUES64_WT_rep1_1_SRR1067557_1.fastq.gz
mv SRR1067557_2.fastq.gz HUES64_WT_rep1_1_SRR1067557_2.fastq.gz
# SRR1067560: HUES64_WT_rep1_2
mv SRR1067560_1.fastq.gz HUES64_WT_rep1_2_SRR1067560_1.fastq.gz
mv SRR1067560_2.fastq.gz HUES64_WT_rep1_2_SRR1067560_2.fastq.gz
# SRR1067564: HUES64_WT_rep1_3
mv SRR1067564_1.fastq.gz HUES64_WT_rep1_3_SRR1067564_1.fastq.gz
mv SRR1067564_2.fastq.gz HUES64_WT_rep1_3_SRR1067564_2.fastq.gz
# SRR1067567: HUES64_WT_rep1_4
mv SRR1067567_1.fastq.gz HUES64_WT_rep1_4_SRR1067567_1.fastq.gz
mv SRR1067567_2.fastq.gz HUES64_WT_rep1_4_SRR1067567_2.fastq.gz
# SRR1067578: HUES64_WT_rep1_5
mv SRR1067578_1.fastq.gz HUES64_WT_rep1_5_SRR1067578_1.fastq.gz
mv SRR1067578_2.fastq.gz HUES64_WT_rep1_5_SRR1067578_2.fastq.gz
# SRR1067565: HUES64_WT_rep2_1
mv SRR1067565_1.fastq.gz HUES64_WT_rep2_1_SRR1067565_1.fastq.gz
mv SRR1067565_2.fastq.gz HUES64_WT_rep2_1_SRR1067565_2.fastq.gz
# SRR1067563: HUES64_WT_rep2_2
mv SRR1067563_1.fastq.gz HUES64_WT_rep2_2_SRR1067563_1.fastq.gz
mv SRR1067563_2.fastq.gz HUES64_WT_rep2_2_SRR1067563_2.fastq.gz
# SRR1067562: HUES64_WT_rep2_3
mv SRR1067562_1.fastq.gz HUES64_WT_rep2_3_SRR1067562_1.fastq.gz
mv SRR1067562_2.fastq.gz HUES64_WT_rep2_3_SRR1067562_2.fastq.gz
# SRR1652272: HUES64_DNMT3AKOEarly_rep1_1
mv SRR1652272_1.fastq.gz HUES64_DNMT3AKOEarly_rep1_1_SRR1652272_1.fastq.gz
mv SRR1652272_2.fastq.gz HUES64_DNMT3AKOEarly_rep1_1_SRR1652272_2.fastq.gz
# SRR1652275: HUES64_DNMT3AKOLate_rep1_1
mv SRR1652275_1.fastq.gz HUES64_DNMT3AKOLate_rep1_1_SRR1652275_1.fastq.gz
mv SRR1652275_2.fastq.gz HUES64_DNMT3AKOLate_rep1_1_SRR1652275_2.fastq.gz
# SRR1652273: HUES64_DNMT3BKOEarly_rep1_1
mv SRR1652273_1.fastq.gz HUES64_DNMT3BKOEarly_rep1_1_SRR1652273_1.fastq.gz
mv SRR1652273_2.fastq.gz HUES64_DNMT3BKOEarly_rep1_1_SRR1652273_2.fastq.gz
# SRR1652276: HUES64_DNMT3BKOLate_rep1_1
mv SRR1652276_1.fastq.gz HUES64_DNMT3BKOLate_rep1_1_SRR1652276_1.fastq.gz
mv SRR1652276_2.fastq.gz HUES64_DNMT3BKOLate_rep1_1_SRR1652276_2.fastq.gz
# SRR1652274: HUES64_DNMT3A3BDoubleKOEarly_rep1_1
mv SRR1652274_1.fastq.gz HUES64_DNMT3A3BDoubleKOEarly_rep1_1_SRR1652274_1.fastq.gz
mv SRR1652274_2.fastq.gz HUES64_DNMT3A3BDoubleKOEarly_rep1_1_SRR1652274_2.fastq.gz
# SRR1652277: HUES64_DNMT3A3BDoubleKOLate_rep1_1
mv SRR1652277_1.fastq.gz HUES64_DNMT3A3BDoubleKOLate_rep1_1_SRR1652277_1.fastq.gz
mv SRR1652277_2.fastq.gz HUES64_DNMT3A3BDoubleKOLate_rep1_1_SRR1652277_2.fastq.gz
```



### Quality check

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

for fq1 in *_1.fastq.gz
do
  fq2=${fq1/_1./_2.}
  bname=${fq1%_1.fastq.gz}
  #echo $fq1, $fq2, $bname
  sbatch -J $bname -o ../fastq_trimmed/$bname.log --mem 16G --wrap "cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 10 -q 20 -o ../fastq_trimmed/$fq1 -p ../fastq_trimmed/$fq2 $fq1 $fq2 > ../fastq_trimmed/$bname.txt"
done
```



### Prepare reference genome

```bash
cd ~

mkdir reference && cd reference

cp ../../GRCh38.p12.genome.fa .

bismark_genome_preparation .
```



### Alignment

```bash
cd ~/fastq_trimmed

mkdir ../bismark

ref=../reference

# KO libraries
for fq1 in *KO*_1.fastq.gz
do
  fq2=${fq1/_1./_2.}
  bname=${fq1%_1.fastq.gz}
  sbatch -J $bname -o ../bismark/$bname.log --mem 64G --wrap "bismark -o ../bismark --parallel 4 $ref -1 $fq1 -2 $fq2"
done

# WT libraries
for fq1 in *WT*_1.fastq.gz
do
  fq2=${fq1/_1./_2.}
  bname=${fq1%_1.fastq.gz}
  sbatch -J $bname -o ../bismark/$bname.log --mem 64G --wrap "bismark -o ../bismark --parallel 4 $ref -1 $fq1 -2 $fq2"
done
```



### Deduplication

```bash
cd ~/bismark

# KO libraries
for bam in *KO*_1_bismark_bt2_pe.bam
do
  bname=${bam%_1_bismark_bt2_pe.bam}
  sbatch -J $bname -o $bname.dedup.log --mem 32GB --wrap "deduplicate_bismark -p --output_dir . --bam $bam"
done

# WT libraries
for id in `ls *WT*_1_bismark_bt2_pe.bam | cut -d "_" -f1-3 | sort | uniq`
do
  bams=`echo ${id}_*_*_1_bismark_bt2_pe.bam`
  sbatch -J $id -o $id.dedup.log --mem 96GB --wrap "deduplicate_bismark -p --output_dir . --bam --multiple $bams"
done
```



### Sort and index deduplicated bam files, extract chr1, index and sort by name

```bash
cd ~/bismark

# KO libraries
for bam in *KO*.deduplicated.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o $bname.log --mem 64G --wrap "samtools sort -@ 20 -T ~/tmp/$bname -o $bname.coordsort.bam $bam && \
  samtools index $bname.coordsort.bam && \
  samtools view -@ 20 $bname.coordsort.bam chr1 -b -o $bname.coordsort.chr1.bam && \
  samtools index $bname.coordsort.chr1.bam && \
  samtools sort -@ 20 -T ~/tmp/$bname.coordsort.chr1 -o $bname.coordsort.chr1.namesort.bam $bname.coordsort.chr1.bam -n && \
  rm $bname.coordsort.bam* $bname.coordsort.chr1.bam*"
done

# WT libraries
for bam in *WT*.deduplicated.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o $bname.log --mem 64G --wrap "samtools sort -@ 20 -T ~/tmp/$bname -o $bname.coordsort.bam $bam && \
  samtools index $bname.coordsort.bam && \
  samtools view -@ 20 $bname.coordsort.bam chr1 -b -o $bname.coordsort.chr1.bam && \
  samtools index $bname.coordsort.chr1.bam && \
  samtools sort -@ 20 -T ~/tmp/$bname.coordsort.chr1 -o $bname.coordsort.chr1.namesort.bam $bname.coordsort.chr1.bam -n && \
  rm $bname.coordsort.bam* $bname.coordsort.chr1.bam*"
done
```



### Methylation chr1

```bash
cd ~/bismark

mkdir ../methylation

ref=../reference

# KO libraries
for bam in *KO*.deduplicated.coordsort.chr1.namesort.bam
do
  bname=${bam%.bam}
  # echo $bname, $ref, $bam
  sbatch -J $bname -o ../methylation/$bname.log --mem 64GB --wrap "bismark_methylation_extractor -p --comprehensive -o ../methylation --gzip --parallel 8 --bedGraph -CX_context --cytosine_report --genome_folder $ref $bam"
done

# WT libraries
for bam in *WT*.deduplicated.coordsort.chr1.namesort.bam
do
  bname=${bam%.bam}
  # echo $bname, $ref, $bam
  sbatch -J $bname -o ../methylation/$bname.log --mem 64GB --wrap "bismark_methylation_extractor -p --comprehensive -o ../methylation --gzip --parallel 8 --bedGraph -CX_context --cytosine_report --genome_folder $ref $bam"
done
```



### Sequence context plots chr1

Obtain sequence context and methylation files:

```bash
## Create chr1 reference
cd ~/reference
samtools faidx GRCh38.p12.genome.fa
samtools faidx GRCh38.p12.genome.fa chr1 > GRCh38.p12.genome.chr1.fa
samtools faidx GRCh38.p12.genome.chr1.fa


## Sequence context file
cd ~/reference

ref=GRCh38.p12.genome.chr1.fa

fastaRegexFinder.py -f $ref -r C | \
cut -f1-3,6 | \
awk -v OFS="\t" '{print $1, $2, $3, ".", ".", $4}' | \
bedtools slop -i - -g $ref.fai -b 6 | \
bedtools getfasta -fi $ref -bed - -bedOut -s | \
awk 'length($7) == 13' | \
awk -v OFS="\t" '{print $1, $2+6, $3-6, $4, $5, $6, $7}' > GRCh38.p12.genome.chr1.context.bed


## Methylation files
cd ~/methylation

# KO libraries
for report in *KO*.CX_report.txt.gz
do
  bname=${report%.CX_report.txt.CX_report.txt.gz}
  nohup zcat $report | \
  grep -P '^chr1\t' | \
  awk -v OFS='\t' '{if ($4+$5 > 5) print $1, $2-1, $2, $4, $5, $3}' | \
  bedtools intersect -a - -b ../reference/GRCh38.p12.genome.chr1.context.bed -loj -s -sorted | \
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $13}' > $bname.context.txt &
done

# WT libraries
for report in *WT*.CX_report.txt.gz
do
  bname=${report%.CX_report.txt.CX_report.txt.gz}
  nohup zcat $report | \
  grep -P '^chr1\t' | \
  awk -v OFS='\t' '{if ($4+$5 > 5) print $1, $2-1, $2, $4, $5, $3}' | \
  bedtools intersect -a - -b ../reference/GRCh38.p12.genome.chr1.context.bed -loj -s -sorted | \
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $13}' > $bname.context.txt &
done


## Concatenate context files
cd ~/methylation
tableCat.py -i *.context.txt | awk '{split($8,a,"_"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"a[2]"_"a[3]}' > HUES64.context.txt
```

Plotting (compact boxplots):

```r
library(data.table)
library(ggplot2)

# Set width
options(width = 300)

# Load data
data <- fread("~/methylation/HUES64.context.txt")
setnames(data, c("chr", "start", "end", "cnt_met", "cnt_unmet", "strand", "context", "library"))

# Calculate percentage methylation
data[, pct_met := round(100 * cnt_met/(cnt_met+cnt_unmet), 2)]


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
data_xcy <- data_xcy[!context_xcy %like% "N"]

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
data_xcgy <- data_xcgy[!context_xcgy %like% "N"]

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
data_wxcgyz <- data_wxcgyz[!context_wxcgyz %like% "N"]

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

# Plot for each library
data_uwxcgyzv <- data_uwxcgyzv[!context_uwxcgyzv %like% "N"]

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
  ggsave(paste("~/figures/", l, "_uwxcgyzv_compact.pdf",sep = ""), height = 140, units= 'cm', limitsize = FALSE)
}

rm(data_uwxcgyzv)
```

Plotting (custom logos/heatmaps):

```r
library(data.table)
library(ggplot2)
library(ggseqlogo)
library(grid)
library(gridExtra)

# Set width
options(width = 300)

# Load data, data was filtered to have more than 5 reads coverage
data <- fread("~/methylation/HUES64.context.txt")
setnames(data, c("chr", "start", "end", "cnt_met", "cnt_unmet", "strand", "context", "library"))

# Calculate percentage methylation
data[, pct_met := round(100 * cnt_met/(cnt_met+cnt_unmet), 2)]


############
# -4_CG_+4 #
############

data_cg <- copy(data[grepl("......CG.....", data$context)])

# Extract -4_CG_+4 context
data_cg[, context_cg := as.vector(sapply(data_cg$context, function(x) paste(unlist(strsplit(x, ""))[3:12], collapse = "")))]

# remove context with N
data_cg <- data_cg[!context_cg %like% "N"]

# % methylation distributions
gg <- ggplot(data_cg[!is.na(pct_met)], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency") +
facet_wrap(~ factor(library, levels = c("WT_rep1", "WT_rep2", "DNMT3AKOEarly_rep1", "DNMT3AKOLate_rep1", "DNMT3BKOEarly_rep1", "DNMT3BKOLate_rep1", "DNMT3A3BDoubleKOEarly_rep1", "DNMT3A3BDoubleKOLate_rep1")), ncol = 2, scales = "free_y") +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black")) +
coord_cartesian(xlim = c(0, 100))
ggsave('~/figures/pctmet_distribution_cg.pdf')

# log y axis
gg <- ggplot(data_cg[!is.na(pct_met)], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency (log10)") +
facet_wrap(~ factor(library, levels = c("WT_rep1", "WT_rep2", "DNMT3AKOEarly_rep1", "DNMT3AKOLate_rep1", "DNMT3BKOEarly_rep1", "DNMT3BKOLate_rep1", "DNMT3A3BDoubleKOEarly_rep1", "DNMT3A3BDoubleKOLate_rep1")), ncol = 2, scales = "free_y") +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black")) +
scale_y_continuous(trans='log10') +
coord_cartesian(xlim = c(0, 100))
ggsave('~/figures/pctmet_distribution_cg_log10.pdf')

# plot genome average presence
for (l in unique(data_cg$library)){
    print(l)
    seqs <- data_cg$context_cg
    # Plot logo
    print("- bits")
    gg_bits <- ggplot() +
    geom_logo(seqs, method = 'bits') +
    ggtitle(paste0("CG coverage")) +
    theme_logo() +
    scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
    theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
    print("- prob")
    gg_prob <- ggplot() +
    geom_logo(seqs, method = 'prob', rev_stack_order = T) +
    theme_logo() +
    scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
    theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
    ggsave(file = paste("/~figures/", l, "_cov_CG.pdf",sep = ""), arrangeGrob(gg_bits, gg_prob, ncol = 1, top = textGrob(l, gp=gpar(fontface="bold", fontsize=25))), width = 6)
}

# Plot for each library for top and bottom 1/5/10%
for (l in unique(data_cg$library)){
    print(l)
    # Obtain filtered and order matrix
    data_cg_filter <- data_cg[library == l & !is.na(pct_met)][order(-pct_met)]
    # Obtain top and bottom sequences
    for (p in c(1,5,10)){
        seqs_top <- data_cg_filter[1:round(nrow(data_cg_filter)/(100/p))]$context_cg
        seqs_bottom <- data_cg_filter[(nrow(data_cg_filter)-round(nrow(data_cg_filter)/(100/p))):nrow(data_cg_filter)]$context_cg
        # Plot logo
        print(paste0("- bits top ", p, "%"))
        gg_bits_top <- ggplot() +
        geom_logo(seqs_top, method = 'bits') +
        ggtitle(paste0("top ", p, "%")) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- prob top ", p, "%"))
        gg_prob_top <- ggplot() +
        geom_logo(seqs_top, method = 'prob', rev_stack_order = T) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- bits bottom ", p, "%"))
        gg_bits_bottom <- ggplot() +
        geom_logo(seqs_bottom, method = 'bits') +
        ggtitle(paste0("bottom ", p, "%")) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- prob bottom ", p, "%"))
        gg_prob_bottom <- ggplot() +
        geom_logo(seqs_bottom, method = 'prob', rev_stack_order = T) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        ggsave(file = paste("~/figures/", l, "_topbottom", p, "pct_CG.pdf",sep = ""), arrangeGrob(gg_bits_top, gg_bits_bottom, gg_prob_top, gg_prob_bottom, ncol = 2, top = textGrob(l, gp=gpar(fontface="bold", fontsize=25))), width = 12)
    }
}

# Plot for each library for top and bottom 100, 500, 1000, 5000, 10000
for (l in unique(data_cg$library)){
    print(l)
    # Obtain filtered and order matrix
    data_cg_filter <- data_cg[library == l & !is.na(pct_met)][order(-pct_met)]
    # Obtain top and bottom sequences
    for (p in c(100, 500, 1000, 5000, 10000)){
        seqs_top <- data_cg_filter[1:p]$context_cg
        seqs_bottom <- data_cg_filter[(nrow(data_cg_filter)-p):nrow(data_cg_filter)]$context_cg
        # Plot logo
        print(paste0("- bits top ", p))
        gg_bits_top <- ggplot() +
        geom_logo(seqs_top, method = 'bits') +
        ggtitle(paste0("top ", p)) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- prob top ", p))
        gg_prob_top <- ggplot() +
        geom_logo(seqs_top, method = 'prob', rev_stack_order = T) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- bits bottom ", p))
        gg_bits_bottom <- ggplot() +
        geom_logo(seqs_bottom, method = 'bits') +
        ggtitle(paste0("bottom ", p)) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- prob bottom ", p))
        gg_prob_bottom <- ggplot() +
        geom_logo(seqs_bottom, method = 'prob', rev_stack_order = T) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        ggsave(file = paste("~/figures/", l, "_topbottom", p, "_CG.pdf",sep = ""), arrangeGrob(gg_bits_top, gg_bits_bottom, gg_prob_top, gg_prob_bottom, ncol = 2, top = textGrob(l, gp=gpar(fontface="bold", fontsize=25))), width = 12)
    }
}

rm(data_cg)


############
# -4_CA_+4 #
############
data_ca <- copy(data[grepl("......CA.....", data$context)])

# Extract -4_CA_+4 context
data_ca[, context_ca := as.vector(sapply(data_ca$context, function(x) paste(unlist(strsplit(x, ""))[3:12], collapse = "")))]

# remove context with N
data_ca <- data_ca[!context_ca %like% "N"]

# % methylation distributions of non-zero
gg <- ggplot(data_ca[pct_met != 0 & !is.na(pct_met)], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency") +
facet_wrap(~ factor(library, levels = c("WT_rep1", "WT_rep2", "DNMT3AKOEarly_rep1", "DNMT3AKOLate_rep1", "DNMT3BKOEarly_rep1", "DNMT3BKOLate_rep1", "DNMT3A3BDoubleKOEarly_rep1", "DNMT3A3BDoubleKOLate_rep1")), nrow = 3, scales = "free_y") +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black")) +
coord_cartesian(xlim = c(0, 100))
ggsave('~/figures/pctmet_distribution_ca_nonzero.pdf')

# log10 y axis
gg <- ggplot(data_ca[pct_met != 0 & !is.na(pct_met)], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency (log10)") +
facet_wrap(~ factor(library, levels = c("WT_rep1", "WT_rep2", "DNMT3AKOEarly_rep1", "DNMT3AKOLate_rep1", "DNMT3BKOEarly_rep1", "DNMT3BKOLate_rep1", "DNMT3A3BDoubleKOEarly_rep1", "DNMT3A3BDoubleKOLate_rep1")), nrow = 3, scales = "free_y") +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black")) +
scale_y_continuous(trans='log10') +
coord_cartesian(xlim = c(0, 100))
ggsave('~/figures/pctmet_distribution_ca_nonzero_log10.pdf')

# plot genome average presence
for (l in unique(data_ca$library)){
    print(l)
    seqs <- data_ca$context_ca
    # Plot logo
    print("- bits")
    gg_bits <- ggplot() +
    geom_logo(seqs, method = 'bits') +
    ggtitle(paste0("CA coverage")) +
    theme_logo() +
    scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
    theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
    print("- prob")
    gg_prob <- ggplot() +
    geom_logo(seqs, method = 'prob', rev_stack_order = T) +
    theme_logo() +
    scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
    theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
    ggsave(file = paste("~/figures/", l, "_cov_CA.pdf",sep = ""), arrangeGrob(gg_bits, gg_prob, ncol = 1, top = textGrob(l, gp=gpar(fontface="bold", fontsize=20))), width = 6)
}

# Plot for each library for top and bottom 1/5/10%
for (l in unique(data_ca$library)){
    print(l)
    # Obtain filtered and order matrix
    data_ca_filter <- data_ca[library == l & !is.na(pct_met)][order(-pct_met)]
    # Obtain top and bottom sequences
    for (p in c(1,5,10)){
        seqs_top <- data_ca_filter[1:round(nrow(data_ca_filter)/(100/p))]$context_ca
        seqs_bottom <- data_ca_filter[(nrow(data_ca_filter)-round(nrow(data_ca_filter)/(100/p))):nrow(data_ca_filter)]$context_ca
        # Plot logo
        print(paste0("- bits top ", p, "%"))
        gg_bits_top <- ggplot() +
        geom_logo(seqs_top, method = 'bits') +
        ggtitle(paste0("top ", p, "%")) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- prob top ", p, "%"))
        gg_prob_top <- ggplot() +
        geom_logo(seqs_top, method = 'prob', rev_stack_order = T) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- bits bottom ", p, "%"))
        gg_bits_bottom <- ggplot() +
        geom_logo(seqs_bottom, method = 'bits') +
        ggtitle(paste0("bottom ", p, "%")) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- prob bottom ", p, "%"))
        gg_prob_bottom <- ggplot() +
        geom_logo(seqs_bottom, method = 'prob', rev_stack_order = T) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        ggsave(file = paste("~/figures/", l, "_topbottom", p, "pct_CA.pdf",sep = ""), arrangeGrob(gg_bits_top, gg_bits_bottom, gg_prob_top, gg_prob_bottom, ncol = 2, top = textGrob(l, gp=gpar(fontface="bold", fontsize=20))), width = 12)
    }
}

# Plot for each library for top and bottom 100, 500, 1000, 5000, 10000
for (l in unique(data_ca$library)){
    print(l)
    # Obtain filtered and order matrix
    data_ca_filter <- data_ca[library == l & !is.na(pct_met)][order(-pct_met)]
    # Obtain top and bottom sequences
    for (p in c(100, 500, 1000, 5000, 10000)){
        seqs_top <- data_ca_filter[1:p]$context_ca
        seqs_bottom <- data_ca_filter[(nrow(data_ca_filter)-p):nrow(data_ca_filter)]$context_ca
        # Plot logo
        print(paste0("- bits top ", p))
        gg_bits_top <- ggplot() +
        geom_logo(seqs_top, method = 'bits') +
        ggtitle(paste0("top ", p)) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- prob top ", p))
        gg_prob_top <- ggplot() +
        geom_logo(seqs_top, method = 'prob', rev_stack_order = T) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- bits bottom ", p))
        gg_bits_bottom <- ggplot() +
        geom_logo(seqs_bottom, method = 'bits') +
        ggtitle(paste0("bottom ", p)) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- prob bottom ", p))
        gg_prob_bottom <- ggplot() +
        geom_logo(seqs_bottom, method = 'prob', rev_stack_order = T) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        ggsave(file = paste("~/figures/", l, "_topbottom", p, "_CA.pdf",sep = ""), arrangeGrob(gg_bits_top, gg_bits_bottom, gg_prob_top, gg_prob_bottom, ncol = 2, top = textGrob(l, gp=gpar(fontface="bold", fontsize=20))), width = 12)
    }
}


################
# -4_nonCpG_+4 #
################

data_nonCpG <- copy(data[!grepl("......CG.....", data$context)])

# Extract -4_nonCpG_+4 context
data_nonCpG[, context_nonCpG := as.vector(sapply(data_nonCpG$context, function(x) paste(unlist(strsplit(x, ""))[3:12], collapse = "")))]

# remove context with N
data_nonCpG <- data_nonCpG[!context_nonCpG %like% "N"]

# % methylation distributions of non-zero
gg <- ggplot(data_nonCpG[pct_met != 0 & !is.na(pct_met)], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency") +
facet_wrap(~ factor(library, levels = c("WT_rep1", "WT_rep2", "DNMT3AKOEarly_rep1", "DNMT3AKOLate_rep1", "DNMT3BKOEarly_rep1", "DNMT3BKOLate_rep1", "DNMT3A3BDoubleKOEarly_rep1", "DNMT3A3BDoubleKOLate_rep1")), nrow = 3, scales = "free_y") +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black")) +
coord_cartesian(xlim = c(0, 100))
ggsave('~/figures/pctmet_distribution_ca_nonzero.pdf')
# log10 y axis
gg <- ggplot(data_nonCpG[pct_met != 0 & !is.na(pct_met)], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency (log10)") +
facet_wrap(~ factor(library, levels = c("WT_rep1", "WT_rep2", "DNMT3AKOEarly_rep1", "DNMT3AKOLate_rep1", "DNMT3BKOEarly_rep1", "DNMT3BKOLate_rep1", "DNMT3A3BDoubleKOEarly_rep1", "DNMT3A3BDoubleKOLate_rep1")), nrow = 3, scales = "free_y") +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black")) +
scale_y_continuous(trans='log10') +
coord_cartesian(xlim = c(0, 100))
ggsave('~/figures/pctmet_distribution_ca_nonzero_log10.pdf')

# plot genome average presence
for (l in unique(data_nonCpG$library)){
    print(l)
    seqs <- data_nonCpG$context_nonCpG
    # Plot logo
    print("- bits")
    gg_bits <- ggplot() +
    geom_logo(seqs, method = 'bits') +
    ggtitle(paste0("nonCpG coverage")) +
    theme_logo() +
    scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
    theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
    print("- prob")
    gg_prob <- ggplot() +
    geom_logo(seqs, method = 'prob', rev_stack_order = T) +
    theme_logo() +
    scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
    theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
    ggsave(file = paste("~/figures/", l, "_cov_nonCpG.pdf",sep = ""), arrangeGrob(gg_bits, gg_prob, ncol = 1, top = textGrob(l, gp=gpar(fontface="bold", fontsize=20))), width = 6)
}

# Plot for each library for top and bottom 1/5/10%
for (l in unique(data_nonCpG$library)){
    print(l)
    # Obtain filtered and order matrix
    data_nonCpG_filter <- data_nonCpG[library == l & !is.na(pct_met)][order(-pct_met)]
    # Obtain top and bottom sequences
    for (p in c(1,5,10)){
        seqs_top <- data_nonCpG_filter[1:round(nrow(data_nonCpG_filter)/(100/p))]$context_nonCpG
        seqs_bottom <- data_nonCpG_filter[(nrow(data_nonCpG_filter)-round(nrow(data_nonCpG_filter)/(100/p))):nrow(data_nonCpG_filter)]$context_nonCpG
        # Plot logo
        print(paste0("- bits top ", p, "%"))
        gg_bits_top <- ggplot() +
        geom_logo(seqs_top, method = 'bits') +
        ggtitle(paste0("top ", p, "%")) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- prob top ", p, "%"))
        gg_prob_top <- ggplot() +
        geom_logo(seqs_top, method = 'prob', rev_stack_order = T) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- bits bottom ", p, "%"))
        gg_bits_bottom <- ggplot() +
        geom_logo(seqs_bottom, method = 'bits') +
        ggtitle(paste0("bottom ", p, "%")) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- prob bottom ", p, "%"))
        gg_prob_bottom <- ggplot() +
        geom_logo(seqs_bottom, method = 'prob', rev_stack_order = T) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        ggsave(file = paste("~/figures/", l, "_topbottom", p, "pct_nonCpG.pdf",sep = ""), arrangeGrob(gg_bits_top, gg_bits_bottom, gg_prob_top, gg_prob_bottom, ncol = 2, top = textGrob(l, gp=gpar(fontface="bold", fontsize=20))), width = 12)
    }
}

# Plot for each library fot top and bottom 100, 500, 1000, 5000, 10000
for (l in unique(data_nonCpG$library)){
    print(l)
    # Obtain filtered and order matrix
    data_nonCpG_filter <- data_nonCpG[library == l & !is.na(pct_met)][order(-pct_met)]
    # Obtain top and bottom sequences
    for (p in c(100, 500, 1000, 5000, 10000)){
        seqs_top <- data_nonCpG_filter[1:p]$context_nonCpG
        seqs_bottom <- data_nonCpG_filter[(nrow(data_nonCpG_filter)-p):nrow(data_nonCpG_filter)]$context_nonCpG
        # Plot logo
        print(paste0("- bits top ", p))
        gg_bits_top <- ggplot() +
        geom_logo(seqs_top, method = 'bits') +
        ggtitle(paste0("top ", p)) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- prob top ", p))
        gg_prob_top <- ggplot() +
        geom_logo(seqs_top, method = 'prob', rev_stack_order = T) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- bits bottom ", p))
        gg_bits_bottom <- ggplot() +
        geom_logo(seqs_bottom, method = 'bits') +
        ggtitle(paste0("bottom ", p)) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        print(paste0("- prob bottom ", p))
        gg_prob_bottom <- ggplot() +
        geom_logo(seqs_bottom, method = 'prob', rev_stack_order = T) +
        theme_logo() +
        scale_x_continuous(breaks = seq(1, 10, 1), labels = c(seq(-4, 0, 1), seq(0, 4, 1))) +
        theme(axis.line.y = element_line(color="black"), axis.title = element_text(size=20), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(size=18, color = "black"), plot.title = element_text(size=20, hjust = 0.5))
        ggsave(file = paste("~/figures/", l, "_topbottom", p, "_nonCpG.pdf",sep = ""), arrangeGrob(gg_bits_top, gg_bits_bottom, gg_prob_top, gg_prob_bottom, ncol = 2, top = textGrob(l, gp=gpar(fontface="bold", fontsize=20))), width = 12)
    }
}


###########
# CAC/CAG #
###########
data_cay <- copy(data[grepl("......CA.....", data$context)])

# Extract CAY context
data_cay[, context_cay := as.vector(sapply(data_cay$context, function(x) paste(unlist(strsplit(x, ""))[7:9], collapse = "")))]

# remove context with N
data_cay <- data_cay[!context_cay %like% "N"]

cac_cag <- dcast(data_cay[, .(.N, pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(library, context_cay)], library ~ context_cay, value.var = "pct_met_mean")
cac_cag[, `:=`(CAC_CAG = CAC/CAG, CAY_CAR = (CAC + CAT)/(CAG + CAA))]
cac_cag
#                       library  CAA  CAC  CAG  CAT   CAC_CAG   CAY_CAR
# 1: DNMT3A3BDoubleKOEarly_rep1 1.62 1.69 1.73 1.63 0.9768786 0.9910448
# 2:  DNMT3A3BDoubleKOLate_rep1 0.13 0.16 0.14 0.18 1.1428571 1.2592593
# 3:         DNMT3AKOEarly_rep1 2.39 2.02 3.35 1.85 0.6029851 0.6742160
# 4:          DNMT3AKOLate_rep1 0.85 0.47 1.52 0.36 0.3092105 0.3502110
# 5:         DNMT3BKOEarly_rep1 1.75 2.40 1.97 1.84 1.2182741 1.1397849
# 6:          DNMT3BKOLate_rep1 0.21 0.93 0.30 0.30 3.1000000 2.4117647
# 7:                    WT_rep1 1.34 1.91 2.18 0.86 0.8761468 0.7869318
# 8:                    WT_rep2 1.16 1.74 1.79 0.84 0.9720670 0.8745763


###########
# CGC/CGG #
###########
data_cgy <- copy(data[grepl("......CG.....", data$context)])

# Extract CGY context
data_cgy[, context_cgy := as.vector(sapply(data_cgy$context, function(x) paste(unlist(strsplit(x, ""))[7:9], collapse = "")))]

# remove context with N
data_cgy <- data_cgy[!context_cgy %like% "N"]

cgc_cgg <- dcast(data_cgy[, .(.N, pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(library, context_cgy)], library ~ context_cgy, value.var = "pct_met_mean")
cgc_cgg[, `:=`(CGC_CGG = CGC/CGG, CAY_CAR = (CGC + CGT)/(CGG + CGA))]
cgc_cgg
#                       library   CGA   CGC   CGG   CGT   CGC_CGG   CAY_CAR
# 1: DNMT3A3BDoubleKOEarly_rep1 79.77 77.84 75.51 74.65 1.0308568 0.9820325
# 2:  DNMT3A3BDoubleKOLate_rep1 65.25 60.02 56.91 59.38 1.0546477 0.9774067
# 3:         DNMT3AKOEarly_rep1 85.35 82.46 83.37 82.51 0.9890848 0.9777738
# 4:          DNMT3AKOLate_rep1 79.50 68.61 71.58 78.35 0.9585080 0.9727297
# 5:         DNMT3BKOEarly_rep1 83.61 84.96 83.25 83.43 1.0205405 1.0091694
# 6:          DNMT3BKOLate_rep1 78.67 73.36 72.83 79.36 1.0072772 1.0080528
# 7:                    WT_rep1 84.35 78.02 78.31 83.29 0.9962968 0.9917005
# 8:                    WT_rep2 84.24 77.20 78.04 83.75 0.9892363 0.9918043
```

