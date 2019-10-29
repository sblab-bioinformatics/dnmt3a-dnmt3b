
TKO (paired-end): [GSE57411](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57411):

  - TKO_DNMT3A2_bisulfite_hiseq2000, [GSM1382253](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1382253)
    - SRX535245: SRR8139275
  - TKO_DNMT3B1_bisulfite_hiseq2000, [GSM1382256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1382256)
    - SRX535248: SRR8139276

WT (single-end): [GSE30202](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30202):

  - MouseES_BisSeq_HiSeq, [GSM748786](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM748786)
    - SRX080191: SRR299053 SRR299054 SRR299055



## Contents

- [Obtain fastq files](mouse_bs.md#obtain-fastq-files)
- [Renaming](mouse_bs.md#renaming)
- [Quality check](mouse_bs.md#quality-check)
- [Trimming](mouse_bs.md#trimming)
- [Prepare reference genome](mouse_bs.md#prepare-reference-genome)
- [Alignment](mouse_bs.md#alignment)
- [Deduplication](mouse_bs.md#deduplication)
- [Sort and index deduplicated bam files, extract chr1, index and sort by name](mouse_bs.md#sort-and-index-deduplicated-bam-files-extract-chr1-index-and-sort-by-name)
- [Extract methylation chr1](mouse_bs.md#extract-methylation-chr1)
- [Sequence context plots chr1](mouse_bs.md#sequence-context-plots-chr1)



### Obtain fastq files

```bash
cd ~

mkdir -p fastq_bsseq && cd fastq_bsseq

# TKO
for i in SRR8139275 SRR8139276
do
  sbatch -J $i -o $i.log --mem 16G --wrap "fastq-dump --outdir . --gzip --split-files -F $i"
done

# WT
for i in SRR299053 SRR299054 SRR299055
do
  sbatch -J $i -o $i.log --mem 16G --wrap "fastq-dump --outdir . --gzip --split-files -F $i"
done
```


### Renaming

```bash
cd ~/fastq_bsseq

# TKO
mv SRR8139275_1.fastq.gz TKO_DNMT3A2_bisulfite_hiseq2000_SRR8139275_1.fastq.gz
mv SRR8139275_2.fastq.gz TKO_DNMT3A2_bisulfite_hiseq2000_SRR8139275_2.fastq.gz
mv SRR8139276_1.fastq.gz TKO_DNMT3B1_bisulfite_hiseq2000_SRR8139276_1.fastq.gz
mv SRR8139276_2.fastq.gz TKO_DNMT3B1_bisulfite_hiseq2000_SRR8139276_2.fastq.gz

# WT
mv SRR299053_1.fastq.gz MouseES_BisSeq_HiSeq_SRR299053_1.fastq.gz
mv SRR299054_1.fastq.gz MouseES_BisSeq_HiSeq_SRR299054_1.fastq.gz
mv SRR299055_1.fastq.gz MouseES_BisSeq_HiSeq_SRR299055_1.fastq.gz
```


### Quality check

```bash
cd ~/fastq_bsseq

mkdir ../fastqc_bsseq

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../fastqc_bsseq/$bname.log --mem 4G --wrap "fastqc --noextract --nogroup -q -o ../fastqc_bsseq $fq"
done
```


### Trimming

```bash
cd ~/fastq_bsseq

mkdir ../fastq_bsseq_trimmed

# TKO
for fq1 in TKO_*_1.fastq.gz
do
  fq2=${fq1/_1./_2.}
  bname=${fq1%_1.fastq.gz}
  #echo $fq1, $fq2, $bname
  sbatch -J $bname -o ../fastq_bsseq_trimmed/$bname.log --mem 16G --wrap "cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 10 -q 20 -o ../fastq_bsseq_trimmed/$fq1 -p ../fastq_bsseq_trimmed/$fq2 $fq1 $fq2 > ../fastq_bsseq_trimmed/$bname.txt"
done

# WT
cd ~/fastq_bsseq_trimmed

for fq in MouseES_BisSeq_HiSeq_*.fastq.gz
do
  bname=${fq%_1.fastq.gz}
  sbatch -J $bname -o ../fastq_bsseq_trimmed/$bname.log --mem 16G --wrap "cutadapt -a AGATCGGAAGAGC -m 10 -q 20 -o ../fastq_bsseq_trimmed/$fq $fq > ../fastq_bsseq_trimmed/$bname.txt"
done
```


### Prepare reference genome

```bash
cd ~

mkdir reference && cd reference

cp ../../GRCm38.p6.genome.fa .

bismark_genome_preparation .
```


### Alignment

```bash
cd ~/fastq_bsseq_trimmed

mkdir ../bismark

ref=../reference


# TKO
for fq1 in TKO_*_1.fastq.gz
do
  fq2=${fq1/_1./_2.}
  bname=${fq1%_1.fastq.gz}
  sbatch -J $bname -o ../bismark/$bname.log --mem 64G --wrap "bismark -o ../bismark --parallel 4 $ref -1 $fq1 -2 $fq2"
done


# WT
for fq in MouseES_BisSeq_HiSeq_*_1.fastq.gz
do
  bname=${fq%_1.fastq.gz}
  sbatch -J $bname -o ../bismark/$bname.log --mem 64G --wrap "bismark -o ../bismark --parallel 4 $ref $fq"
done
```


### Deduplication

```bash
cd ~/bismark

# TKO
for bam in TKO*_1_bismark_bt2_pe.bam
do
  bname=${bam%_1_bismark_bt2_pe.bam}
  sbatch -J $bname -o $bname.dedup.log --mem 32GB --wrap "deduplicate_bismark -p --output_dir . --bam $bam"
done

# WT
for id in `ls MouseES_BisSeq_HiSeq_*_1_bismark_bt2.bam | cut -d "_" -f1-3 | sort | uniq`
do
  bams=`echo ${id}_*_1_bismark_bt2.bam`
  sbatch -J $id -o $id.dedup.log --mem 32GB --wrap "deduplicate_bismark -s --output_dir . --bam --multiple $bams"
done
```


### Sort and index deduplicated bam files, extract chr1, index and sort by name

```bash
cd ~/bismark

# TKO
for bam in TKO*.deduplicated.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o $bname.log --mem 64G --wrap "samtools sort -@ 20 -T ~/tmp/$bname -o $bname.coordsort.bam $bam && \
  samtools index $bname.coordsort.bam && \
  samtools view -@ 20 $bname.coordsort.bam chr1 -b -o $bname.coordsort.chr1.bam && \
  samtools index $bname.coordsort.chr1.bam && \
  samtools sort -@ 20 -T ~/tmp/$bname.coordsort.chr1 -o $bname.coordsort.chr1.namesort.bam $bname.coordsort.chr1.bam -n && \
  rm $bname.coordsort.bam* $bname.coordsort.chr1.bam*"
done

# WT
for bam in MouseES_BisSeq_HiSeq*.deduplicated.bam
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


### Extract methylation chr1

```bash
cd ~/bismark

mkdir ../methylation

ref=../reference

# TKO
for bam in TKO*.deduplicated.coordsort.chr1.namesort.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o ../methylation/$bname.log --mem 64GB --wrap "bismark_methylation_extractor -p --comprehensive -o ../methylation --gzip --parallel 8 --bedGraph -CX_context --cytosine_report --genome_folder $ref $bam"
done

# WT
for bam in MouseES_BisSeq_HiSeq*.deduplicated.coordsort.chr1.namesort.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o ../methylation/$bname.log --mem 64GB --wrap "bismark_methylation_extractor -s --comprehensive -o ../methylation --gzip --parallel 8 --bedGraph -CX_context --cytosine_report --genome_folder $ref $bam"
done
```


### Sequence context plots chr1

Obtain sequence context and methylation files:

```bash
cd ~
mkdir reference2 && cd reference2

samtools faidx ../reference/GRCm38.p6.genome.fa
samtools faidx ../reference/GRCm38.p6.genome.fa chr1 > GRCm38.p6.genome.chr1.fa
samtools faidx GRCm38.p6.genome.chr1.fa


## Sequence context file
cd ~/reference2

ref=GRCm38.p6.genome.chr1.fa

fastaRegexFinder.py -f $ref -r C | \
cut -f1-3,6 | \
awk -v OFS="\t" '{print $1, $2, $3, ".", ".", $4}' | \
bedtools slop -i - -g $ref.fai -b 6 | \
bedtools getfasta -fi $ref -bed - -bedOut -s | \
awk 'length($7) == 13' | \
awk -v OFS="\t" '{print $1, $2+6, $3-6, $4, $5, $6, $7}' > GRCm38.p6.genome.chr1.context.bed


## Methylation files
cd ~/methylation

# TKO
for report in TKO*.CX_report.txt.gz
do
  bname=${report%.CX_report.txt.CX_report.txt.gz}
  nohup zcat $report | \
  grep -P '^chr1\t' | \
  awk -v OFS='\t' '{if ($4+$5 > 5) print $1, $2-1, $2, $4, $5, $3}' | \
  bedtools intersect -a - -b ../reference2/GRCm38.p6.genome.chr1.context.bed -loj -s -sorted | \
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $13}' > $bname.context.txt &
done

# WT
for report in MouseES_BisSeq_HiSeq*.CX_report.txt.gz
do
  bname=${report%.CX_report.txt.CX_report.txt.gz}
  nohup zcat $report | \
  grep -P '^chr1\t' | \
  awk -v OFS='\t' '{if ($4+$5 > 5) print $1, $2-1, $2, $4, $5, $3}' | \
  bedtools intersect -a - -b ../reference2/GRCm38.p6.genome.chr1.context.bed -loj -s -sorted | \
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $13}' > $bname.context.txt &
done

## Concatenate context files
cd ~/methylation

tableCat.py -i *.context.txt | awk '{split($8,a,"_"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"a[1]"_"a[2]"_"a[3]"_"a[4]}' > TKO_DNMT3A2_DNMT3B1_MouseES.context.txt
```

Plotting:

```r
library(data.table)
library(ggplot2)
library(ggseqlogo)
library(grid)
library(gridExtra)

# Set width
options(width = 300)

# Load data, data was filtered to have more than 5 reads coverage
data <- fread("TKO_DNMT3A2_DNMT3B1_MouseES.context.txt")
setnames(data, c("chr", "start", "end", "cnt_met", "cnt_unmet", "strand", "context", "library"))

# Calculate percentage methylation
data[, pct_met := round(100 * cnt_met/(cnt_met+cnt_unmet), 2)]

# Rename library fields
data[library == "MouseES_BisSeq_HiSeq_SRR299053", library := "WT"]
data[library == "TKO_DNMT3A2_bisulfite_hiseq2000", library := "TKO-DNMT3A2"]
data[library == "TKO_DNMT3B1_bisulfite_hiseq2000", library := "TKO-DNMT3B1"]


############ 
# -4_CG_+4 #
############
data_cg <- copy(data[grepl("......CG.....", data$context)])

# Extract -4_CG_+4 context
data_cg[, context_cg := as.vector(sapply(data_cg$context, function(x) paste(unlist(strsplit(x, ""))[3:12], collapse = "")))]

# remove context with N
data_cg <- data_cg[!context_cg %like% "N"]

table(data_cg$library)
# TKO-DNMT3A2 TKO-DNMT3B1          WT
#      169336      112897      746877

# % methylation distributions
gg <- ggplot(data_cg[!is.na(pct_met)], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency") +
facet_wrap(~ factor(library, levels = c("WT", "TKO-DNMT3A2", "TKO-DNMT3B1")), nrow = 3, scales = "free_y") +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black")) +
coord_cartesian(xlim = c(0, 100))
ggsave('~/figures/pctmet_distribution_cg.pdf')

# log y axis
gg <- ggplot(data_cg[!is.na(pct_met)], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency (log10)") +
facet_wrap(~ factor(library, levels = c("WT", "TKO-DNMT3A2", "TKO-DNMT3B1")), nrow = 3, scales = "free_y") +
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
    ggsave(file = paste("~/figures/", l, "_cov_CG.pdf",sep = ""), arrangeGrob(gg_bits, gg_prob, ncol = 1, top = textGrob(l, gp=gpar(fontface="bold", fontsize=25))), width = 6)
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

# Plot for each library for top and bottom 100, 500, 1000, 2000, 5000
for (l in unique(data_cg$library)){
    print(l)
    # Obtain filtered and order matrix
    data_cg_filter <- data_cg[library == l & !is.na(pct_met)][order(-pct_met)]
    # Obtain top and bottom sequences
    for (p in c(100, 500, 1000, 2000, 5000)){
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

table(data_ca$library)

# % methylation distributions of non-zero
gg <- ggplot(data_ca[pct_met != 0 & !is.na(pct_met)], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency (log10)") +
facet_wrap(~ factor(library, levels = c("WT", "TKO-DNMT3A2", "TKO-DNMT3B1")), nrow = 3, scales = "free_y") +
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

# Plot for each library fot top and bottom 1/5/10%
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

# Plot for each library fot top and bottom 100, 500, 1000, 2000, 5000
for (l in unique(data_ca$library)){
    print(l)
    # Obtain filtered and order matrix
    data_ca_filter <- data_ca[library == l & !is.na(pct_met)][order(-pct_met)]
    # Obtain top and bottom sequences
    for (p in c(100, 500, 1000, 2000, 5000)){
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

table(data_nonCpG$library)

# % methylation distributions of non-zero
gg <- ggplot(data_nonCpG[pct_met != 0 & !is.na(pct_met)], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency") +
facet_wrap(~ factor(library, levels = c("WT", "TKO-DNMT3A2", "TKO-DNMT3B1")), nrow = 3, scales = "free_y") +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black")) +
coord_cartesian(xlim = c(0, 100))
ggsave('~/figures/pctmet_distribution_nonCpG_nonzero.pdf')

# log y axis
gg <- ggplot(data_nonCpG[pct_met != 0 & !is.na(pct_met)], aes(x=pct_met)) +
geom_histogram(binwidth=1) +
theme_bw() +
xlab("% methylation") +
ylab("frequency (log10)") +
facet_wrap(~ factor(library, levels = c("WT", "TKO-DNMT3A2", "TKO-DNMT3B1")), nrow = 3, scales = "free_y") +
theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black")) +
scale_y_continuous(trans='log10') +
coord_cartesian(xlim = c(0, 100))
ggsave('~/figures/pctmet_distribution_nonCpG_nonzero_log10.pdf')

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

# Plot for each library fot top and bottom 1/5/10%
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

# Plot for each library fot top and bottom 100, 500, 1000, 2000, 5000
for (l in unique(data_nonCpG$library)){
    print(l)
    # Obtain filtered and order matrix
    data_nonCpG_filter <- data_nonCpG[library == l & !is.na(pct_met)][order(-pct_met)]
    # Obtain top and bottom sequences
    for (p in c(100, 500, 1000, 2000, 5000)){
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

table(data_cay$library)
# TKO-DNMT3A2 TKO-DNMT3B1          WT
#     1245759      806670     6363697

cac_cag <- dcast(data_cay[, .(.N, pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(library, context_cay)], library ~ context_cay, value.var = "pct_met_mean")
cac_cag[, `:=`(CAC_CAG = CAC/CAG, CAY_CAR = (CAC + CAT)/(CAG + CAA))]
cac_cag
#        library  CAA  CAC  CAG  CAT   CAC_CAG   CAY_CAR
# 1: TKO-DNMT3A2 0.13 0.61 0.22 0.19 2.7727273 2.2857143
# 2: TKO-DNMT3B1 0.15 0.08 0.23 0.07 0.3478261 0.3947368
# 3:          WT 2.54 3.87 3.23 2.21 1.1981424 1.0537262


###########
# CGC/CGG #
###########

data_cgy <- copy(data[grepl("......CG.....", data$context)])

# Extract CGY context
data_cgy[, context_cgy := as.vector(sapply(data_cgy$context, function(x) paste(unlist(strsplit(x, ""))[7:9], collapse = "")))]

# remove context with N
data_cgy <- data_cgy[!context_cgy %like% "N"]

table(data_cgy$library)
# TKO-DNMT3A2 TKO-DNMT3B1          WT
#      169336      112897      746877

cgc_cgg <- dcast(data_cgy[, .(.N, pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(library, context_cgy)], library ~ context_cgy, value.var = "pct_met_mean")
cgc_cgg[, `:=`(CGC_CGG = CGC/CGG, CAY_CAR = (CGC + CGT)/(CGG + CGA))]
cgc_cgg
#        library   CGA   CGC   CGG   CGT   CGC_CGG   CAY_CAR
# 1: TKO-DNMT3A2  5.83 10.15  4.53  7.91 2.2406181 1.7432432
# 2: TKO-DNMT3B1  3.90  1.31  3.97  1.02 0.3299748 0.2960610
# 3:          WT 83.12 79.26 81.08 83.17 0.9775530 0.9892205
```
