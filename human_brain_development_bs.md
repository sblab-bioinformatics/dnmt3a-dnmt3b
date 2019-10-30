### Download datasets

```sh
cd ~
mkdir data
cd data

# get brain chr1 data from GSE47966 and HUES6 from GSM1173778
# MethylC-Seq_hs_fc_fetal
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1163nnn/GSM1163695/suppl/GSM1163695_allC.MethylC-Seq_hs_fc_fetal.chr1.txt.txt.gz
mv GSM1163695_allC.MethylC-Seq_hs_fc_fetal.chr1.txt.txt.gz GSM1163695_allC.MethylC-Seq_hs_fc_fetal.chr1.txt.gz

# MethylC-Seq_hs_mfg_12yr
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1164nnn/GSM1164630/suppl/GSM1164630_allC.MethylC-Seq_hs_mfg_12yr.chr1.txt.gz

# MethylC-Seq_hs_mfg_16yr
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1164nnn/GSM1164631/suppl/GSM1164631_allC.MethylC-Seq_hs_mfg_16yr.chr1.txt.gz

# MethylC-Seq_hs_mfg_25yr
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1164nnn/GSM1164632/suppl/GSM1164632_allC.MethylC-Seq_hs_mfg_25yr.chr1.txt.gz

# MethylC-Seq_hs_mfg_5yr
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1166nnn/GSM1166274/suppl/GSM1166274_allC.MethylC-Seq_hs_mfg_5yr.chr1.txt.gz

# MethylC-Seq_hs_mfg_35do
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1167nnn/GSM1167004/suppl/GSM1167004_allC.MethylC-Seq_hs_mfg_35do.chr1.txt.gz

# MethylC-Seq_hs_mfg_2yr
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1167nnn/GSM1167005/suppl/GSM1167005_allC.MethylC-Seq_hs_mfg_2yr.chr1.txt.gz

# MethylC-Seq_hs_fc_64yr
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173772/suppl/GSM1173772_allC.MethylC-Seq_hs_fc_64yr.chr1.txt.gz

# MethylC-Seq_hs_fc_female_53yr_NeuN_pos
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173773/suppl/GSM1173773_allC.MethylC-Seq_hs_fc_female_53yr_NeuN_pos.chr1.txt.gz

# MethylC-Seq_hs_fc_female_53yr_NeuN_neg
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173774/suppl/GSM1173774_allC.MethylC-Seq_hs_fc_female_53yr_NeuN_neg.chr1.txt.gz

# MethylC-Seq_hs_fc_male_55yr_tissue
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173775/suppl/GSM1173775_allC.MethylC-Seq_hs_fc_male_55yr_tissue.chr1.txt.gz

# MethylC-Seq_hs_fc_male_55yr_NeuN_pos
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173776/suppl/GSM1173776_allC.MethylC-Seq_hs_fc_male_55yr_NeuN_pos.chr1.txt.gz

# MethylC-Seq_hs_fc_male_55yr_NeuN_neg
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173777/suppl/GSM1173777_allC.MethylC-Seq_hs_fc_male_55yr_NeuN_neg.chr1.txt.gz

# MethylC-Seq_hs_hues6
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173778/suppl/GSM1173778_allC.MethylC-Seq_hs_hues6.chr1.txt.gz

```

### prepare hg18 reference and methylation context

```sh
srun --mem 96G --pty /usr/bin/bash


## Create chr1 reference
cd ~
mkdir reference && cd reference

ref=hg18_genome.fa

samtools faidx $ref
samtools faidx $ref chr1 > hg18_genome.chr1.fa
samtools faidx hg18_genome.chr1.fa

ref=hg18_genome.chr1.fa

fastaRegexFinder.py -f $ref -r C | \
cut -f1-3,6 | \
awk -v OFS="\t" '{print $1, $2, $3, ".", ".", $4}' | \
bedtools slop -i - -g $ref.fai -b 6 | \
bedtools getfasta -fi $ref -bed - -bedOut -s | \
awk 'length($7) == 13' | \
awk -v OFS="\t" '{print $1, $2+6, $3-6, $4, $5, $6, toupper($7)}' > hg18_genome.chr1.context.bed

wc -l hg18_genome.chr1.context.bed 
# 93921242 Cs, with +/- 6bp sequence context in chr1 both in forward and reverse strands


## Sequence context file
cd ~/data

for report in *.chr1.txt.gz
do
  bname=${report%.chr1.txt.gz}
  nohup zcat $report | \
  awk -v OFS='\t' '{if ($6 > 10) print "chr1", $2, $2+1, $5, $6, $3}' | \
  bedtools intersect -a - -b ../reference/hg18_genome.chr1.context.bed -loj -s -sorted | \
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $13}' > $bname.context.txt &
done

wc -l *.context.txt 
#    26666299 GSM1163695_allC.MethylC-Seq_hs_fc_fetal.context.txt
#    43770632 GSM1164630_allC.MethylC-Seq_hs_mfg_12yr.context.txt
#    44567083 GSM1164631_allC.MethylC-Seq_hs_mfg_16yr.context.txt
#    41767596 GSM1164632_allC.MethylC-Seq_hs_mfg_25yr.context.txt
#    42911516 GSM1166274_allC.MethylC-Seq_hs_mfg_5yr.context.txt
#    41535200 GSM1167004_allC.MethylC-Seq_hs_mfg_35do.context.txt
#    45409169 GSM1167005_allC.MethylC-Seq_hs_mfg_2yr.context.txt
#    11292867 GSM1173772_allC.MethylC-Seq_hs_fc_64yr.context.txt
#    32893028 GSM1173773_allC.MethylC-Seq_hs_fc_female_53yr_NeuN_pos.context.txt
#    13858409 GSM1173774_allC.MethylC-Seq_hs_fc_female_53yr_NeuN_neg.context.txt
#       60001 GSM1173775_allC.MethylC-Seq_hs_fc_male_55yr_tissue.context.txt
#    26296266 GSM1173776_allC.MethylC-Seq_hs_fc_male_55yr_NeuN_pos.context.txt
#    21802900 GSM1173777_allC.MethylC-Seq_hs_fc_male_55yr_NeuN_neg.context.txt
#     7366782 GSM1173778_allC.MethylC-Seq_hs_hues6.context.txt

for report in GSM1173775_allC.MethylC-Seq_hs_fc_male_55yr_tissue.chr1.txt.gz
do
  bname=${report%.chr1.txt.gz}
  nohup zcat $report | \
  awk -v OFS='\t' '{if ($6 > 5) print "chr1", $2, $2+1, $5, $6, $3}' | \
  bedtools intersect -a - -b ../reference/hg18_genome.chr1.context.bed -loj -s -sorted | \
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $13}' > $bname.cov5.context.txt &
done

wc -l GSM1173775_allC.MethylC-Seq_hs_fc_male_55yr_tissue*
#    60001 GSM1173775_allC.MethylC-Seq_hs_fc_male_55yr_tissue.context.txt
#  1367135 GSM1173775_allC.MethylC-Seq_hs_fc_male_55yr_tissue.cov5.context.txt

# Concatenate context files
cd ~/data

### mfg
tableCat.py -i GSM*mfg*.context.txt -r .context.txt | awk '{split($8,a,"_"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"a[3]"_"a[4]"_"a[5]}' > mfg.context.txt &

### fc
tableCat.py -i GSM*fc*.context.txt -r .context.txt | awk '{split($8,a,"_"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"a[3]"_"a[4]"_"a[5]"_"a[6]"_"a[7]"_"a[8]}' > fc.context.txt &

### hues6
cp GSM1173778_allC.MethylC-Seq_hs_hues6.context.txt hues6.context.txt

```

### load mfg data in R

```sh
library(data.table)
library(ggplot2)

# Set width
options(width = 300)

# Load mfg data
data <- fread("~/data/mfg.context.txt")
setnames(data, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand", "context", "library"))

# Calculate percentage methylation
data[, pct_met := round(100 * cnt_met/(cnt_tot), 2)]

table(data$library)
# hs_mfg_12yr hs_mfg_16yr hs_mfg_25yr  hs_mfg_2yr hs_mfg_35do  hs_mfg_5yr
#    43770632    44567083    41767596    45409169    41535200    42911516
```


#### extract CGN methylation

```sh
data_cgy <- copy(data[grepl("......CG.....", data$context)])

# remove context with N
# The like operator is a simple wrapper for grep(..., value=TRUE)
data_cgy <- data_cgy[!context_cgy %like% "N"]


# Extract CGY context
data_cgy[, context_cgy := as.vector(sapply(data_cgy$context, function(x) paste(unlist(strsplit(x, ""))[7:9], collapse = "")))]


# Explore CGY context
## hs_mfg_25yr
data_cgy[library == "hs_mfg_25yr", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cgy)][order(-pct_met_median)]
#    context_cgy      N pct_met_median pct_met_mean
# 1:         CGC 273064          94.44        88.45
# 2:         CGG 361337          92.31        85.41
# 3:         CGA 386166          91.67        84.75
# 4:         CGT 454532          90.91        83.20

## hs_mfg_35do
data_cgy[library == "hs_mfg_35do", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cgy)][order(-pct_met_median)]
#    context_cgy      N pct_met_median pct_met_mean
# 1:         CGC 271625          94.44        87.56
# 2:         CGG 367447          92.31        83.90
# 3:         CGA 394385          91.67        82.61
# 4:         CGT 453847          90.00        80.68

############ CGC/CGG ##############
cgc_cgg <- dcast(data_cgy[, .(pct_met = mean(pct_met, na.rm=TRUE)), by = .(library, context_cgy)], library ~ context_cgy, value.var = "pct_met")
cgc_cgg[, CGC_CGG := CGC/CGG]
cgc_cgg
#        library      CGA      CGC      CGG      CGT  CGC_CGG
# 1: hs_mfg_12yr 84.44367 88.31599 85.07937 82.92465 1.038042
# 2: hs_mfg_16yr 84.59878 88.17560 85.02975 83.12611 1.036997
# 3: hs_mfg_25yr 84.75081 88.44712 85.40665 83.19899 1.035600
# 4:  hs_mfg_2yr 83.83509 87.89283 84.80865 82.29094 1.036366
# 5: hs_mfg_35do 82.60919 87.55776 83.90150 80.67508 1.043578
# 6:  hs_mfg_5yr 84.26716 88.39173 85.12829 82.69841 1.038336

cgc_cgg_median <- dcast(data_cgy[, .(pct_met = median(pct_met, na.rm=TRUE)), by = .(library, context_cgy)], library ~ context_cgy, value.var = "pct_met")
cgc_cgg_median[, CGC_CGG := CGC/CGG]
cgc_cgg_median
#        library   CGA   CGC   CGG   CGT  CGC_CGG
# 1: hs_mfg_12yr 91.67 94.44 92.31 90.91 1.023074
# 2: hs_mfg_16yr 91.67 94.44 92.31 90.91 1.023074
# 3: hs_mfg_25yr 91.67 94.44 92.31 90.91 1.023074
# 4:  hs_mfg_2yr 91.67 94.12 92.31 90.91 1.019608
# 5: hs_mfg_35do 91.67 94.44 92.31 90.00 1.023074
# 6:  hs_mfg_5yr 91.67 94.44 92.31 90.91 1.023074

gg <- ggplot(cgc_cgg, aes(x=library, y=CGC_CGG)) +
geom_bar(stat="identity", color="black", position=position_dodge(), width= 0.75, alpha = 0.5) +
theme_classic() +
ylab(expression("[mCGC/CGC]/[mCGG/CGG]")) +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 12, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black"), aspect.ratio = 5/3) +
coord_cartesian(ylim = c(0, 1.2))
ggsave("~/data/brain_cgn/mfg_CGC_CGG.pdf")


rm(cgc_cgg)
rm(cgc_cgg_median)
rm(data_cgy)
```

#### extract CAN methylation

```sh
data_cay <- copy(data[grepl("......CA.....", data$context)])

# Extract CAY context
data_cay[, context_cay := as.vector(sapply(data_cay$context, function(x) paste(unlist(strsplit(x, ""))[7:9], collapse = "")))]

# remove context with N
# The like operator is a simple wrapper for grep(..., value=TRUE)
data_cay <- data_cay[!context_cay %like% "N"]

# Explore CAY context
## hs_mfg_35do
data_cay[library == "hs_mfg_35do", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
# 1:         CAC 2918099              0         2.77
# 2:         CAT 4438671              0         1.14
# 3:         CAG 4306944              0         1.10
# 4:         CAA 4942410              0         0.87

## hs_mfg_2yr
data_cay[library == "hs_mfg_2yr", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
# 1:         CAC 3188307              0         5.98
# 2:         CAG 4689169              0         2.21
# 3:         CAT 4772060              0         2.18
# 4:         CAA 5180649              0         1.56

## hs_mfg_5yr
data_cay[library == "hs_mfg_5yr", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
# 1:         CAC 3019585              0         5.54
# 2:         CAG 4418025              0         2.12
# 3:         CAT 4568533              0         2.09
# 4:         CAA 4998206              0         1.53

## hs_mfg_12yr
data_cay[library == "hs_mfg_12yr", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
# 1:         CAC 3076610              0         6.82
# 2:         CAG 4477913              0         2.62
# 3:         CAT 4653505              0         2.51
# 4:         CAA 5028765              0         1.81

## hs_mfg_16yr
data_cay[library == "hs_mfg_16yr", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
1:         CAC 3130041           4.76         7.78
2:         CAG 4591086           0.00         2.99
3:         CAT 4680165           0.00         2.87
4:         CAA 5054790           0.00         2.08

## hs_mfg_25yr
data_cay[library == "hs_mfg_25yr", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
# 1:         CAC 2937173              5         7.87
# 2:         CAG 4289554              0         3.06
# 3:         CAT 4443736              0         2.94
# 4:         CAA 4833980              0         2.14


############ CAC/CAG ##############
cac_cag <- dcast(data_cay[, .(pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(library, context_cay)], library ~ context_cay, value.var = "pct_met_mean")
cac_cag[, CAC_CAG := CAC/CAG][order(CAC_CAG)]
#        library  CAA  CAC  CAG  CAT  CAC_CAG
# 1: hs_mfg_35do 0.87 2.77 1.10 1.14 2.518182
# 2: hs_mfg_25yr 2.14 7.87 3.06 2.94 2.571895
# 3: hs_mfg_16yr 2.08 7.78 2.99 2.87 2.602007
# 4: hs_mfg_12yr 1.81 6.82 2.62 2.51 2.603053
# 5:  hs_mfg_5yr 1.53 5.54 2.12 2.09 2.613208
# 6:  hs_mfg_2yr 1.56 5.98 2.21 2.18 2.705882

# set the order of libraries
cac_cag$library <- factor(cac_cag$library, level = c("hs_mfg_35do", "hs_mfg_2yr", "hs_mfg_5yr", "hs_mfg_12yr", "hs_mfg_16yr", "hs_mfg_25yr"))

gg <- ggplot(cac_cag, aes(x=library, y=CAC_CAG)) +
geom_bar(stat="identity", color="black", position=position_dodge(), width= 0.75) +
theme_bw() +
ylab(expression("[mCAC/CAC]/[mCAG/CAG]")) +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 12, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black"), aspect.ratio = 5/3) +
coord_cartesian(ylim = c(0, 3))
ggsave("~/data/brain_can/mfg_CAC_CAG.pdf")


########### plot CAN changes through age ###########
cac_cag_melt <- melt(cac_cag, id.vars = "library", variable.name = "CAN", value.name = "pct_met")

# set the order of libraries
cac_cag_melt$library <- factor(cac_cag_melt$library, level = c("hs_mfg_35do", "hs_mfg_2yr", "hs_mfg_5yr", "hs_mfg_12yr", "hs_mfg_16yr", "hs_mfg_25yr"))

write.table(cac_cag_melt, file = "~/data/brain_can/mfg_CAC_CAG.txt", row.names = FALSE, quote = FALSE, sep = "\t")

gg <- ggplot(cac_cag_melt, aes(x=library, y= pct_met, color = CAN, group = CAN)) +
geom_line(aes(linetype = CAN), lwd = 1) +
scale_linetype_manual(values = c(1,1,1,1,2)) +
theme_bw() +
ylab(expression("[mCAC/CAC]/[mCAG/CAG]")) + 
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 12, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black"), aspect.ratio = 5/3) +
coord_cartesian(ylim = c(0, 8))
ggsave("~/data/brain_can/mfg_CAC_CAG_age.pdf")

########################################################################
########### plot CAN changes through age (fc datasets added) ###########
# combine and edit mfg_CAC_CAG.txt and fc.CAC_CAG.txt locally
cac_cag_age <- fread("~/data/brain_can/mfg_fc_CAC_CAG.txt")

# set the order of libraries by ranking CAC_CAG
cac_cag_age$library <- factor(cac_cag_age$library, level = c("hs_fc_fetal___", "hs_mfg_35do", "hs_mfg_2yr", "hs_mfg_5yr", "hs_mfg_12yr", "hs_mfg_16yr", "hs_mfg_25yr", "hs_fc_male_55yr_tissue_", "hs_fc_64yr___"))


gg <- ggplot(cac_cag_age, aes(x=library, y= pct_met, color = CAN, group = CAN)) +
geom_line(aes(linetype = CAN), lwd = 1) +
scale_linetype_manual(values = c(1,1,2,1,1)) +
theme_bw() +
ylab(expression("[mCAC/CAC]/[mCAG/CAG]")) + 
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 12, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black"), aspect.ratio = 5/3) +
coord_cartesian(ylim = c(0, 8))
ggsave("~/data/brain_can/mfg_fc_CAC_CAG_age.pdf")

rm(cac_cag_age)
rm(cac_cag_melt)
rm(cac_cag)
rm(data_cay)
```

### load fc data in R

```sh
library(data.table)
library(ggplot2)

# Set width
options(width = 300)

# Load fc data
data <- fread("~/data/fc.context.txt")
setnames(data, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand", "context", "library"))

# Calculate percentage methylation
data[, pct_met := round(100 * cnt_met/(cnt_tot), 2)]

table(data$library)
# hs_fc_64yr___ hs_fc_female_53yr_NeuN_neg hs_fc_female_53yr_NeuN_pos             hs_fc_fetal___   hs_fc_male_55yr_NeuN_neg   hs_fc_male_55yr_NeuN_pos    hs_fc_male_55yr_tissue_
#     11292867                   13858409                   32893028                   26666299                   21802900                   26296266                      60001
```


#### extract CGN methylation

```sh
data_cgy <- copy(data[grepl("......CG.....", data$context)])

# Extract CGY context
data_cgy[, context_cgy := as.vector(sapply(data_cgy$context, function(x) paste(unlist(strsplit(x, ""))[7:9], collapse = "")))]

# remove context with N
# The like operator is a simple wrapper for grep(..., value=TRUE)
data_cgy <- data_cgy[!context_cgy %like% "N"]

# Explore CGY context
## hs_fc_fetal___
data_cgy[library == "hs_fc_fetal___", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cgy)][order(-pct_met_median)]
#    context_cgy      N pct_met_median pct_met_mean
# 1:         CGC 434942          94.44        76.27
# 2:         CGA 394745          92.86        79.48
# 3:         CGG 574789          92.31        74.51
# 4:         CGT 450613          91.67        78.28

## hs_fc_female_53yr_NeuN_neg
data_cgy[library == "hs_fc_female_53yr_NeuN_neg", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cgy)][order(-pct_met_mean)]
#    context_cgy      N pct_met_median pct_met_mean
# 1:         CGC  64225         100.00        88.79
# 2:         CGA 116345          92.86        85.56
# 3:         CGG  89495          92.86        85.49
# 4:         CGT 142204          91.67        83.78

## hs_fc_female_53yr_NeuN_pos
data_cgy[library == "hs_fc_female_53yr_NeuN_pos", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cgy)][order(-pct_met_mean)]
#    context_cgy      N pct_met_median pct_met_mean
# 1:         CGA 409743          95.74        86.27
# 2:         CGT 441069          95.12        85.67
# 3:         CGC 338063          98.25        85.40
# 4:         CGG 508144          95.35        84.31

############ CGC/CGG ##############
cgc_cgg <- dcast(data_cgy[, .(pct_met = mean(pct_met, na.rm=TRUE)), by = .(library, context_cgy)], library ~ context_cgy, value.var = "pct_met")
cgc_cgg[, CGC_CGG := CGC/CGG]
cgc_cgg
#                       library      CGA      CGC      CGG      CGT   CGC_CGG
# 1:              hs_fc_64yr___ 78.43547 77.30880 75.74430 77.48458 1.0206550
# 2: hs_fc_female_53yr_NeuN_neg 85.56226 88.79166 85.49066 83.78242 1.0386125
# 3: hs_fc_female_53yr_NeuN_pos 86.26595 85.40203 84.30661 85.67425 1.0129933
# 4:             hs_fc_fetal___ 79.47817 76.26736 74.51160 78.28011 1.0235636
# 5:   hs_fc_male_55yr_NeuN_neg 85.40804 89.12862 85.73595 83.73348 1.0395711
# 6:   hs_fc_male_55yr_NeuN_pos 88.20676 90.29616 88.34535 87.48740 1.0220817
# 7:    hs_fc_male_55yr_tissue_ 70.09863 65.99386 71.86425 71.34957 0.9183127

cgc_cgg_median <- dcast(data_cgy[, .(pct_met = median(pct_met, na.rm=TRUE)), by = .(library, context_cgy)], library ~ context_cgy, value.var = "pct_met")
cgc_cgg_median[, CGC_CGG := CGC/CGG]
cgc_cgg_median
#                       library    CGA    CGC    CGG   CGT   CGC_CGG
# 1:              hs_fc_64yr___  90.91  91.67  90.91 87.50 1.0083599
# 2: hs_fc_female_53yr_NeuN_neg  92.86 100.00  92.86 91.67 1.0768899
# 3: hs_fc_female_53yr_NeuN_pos  95.74  98.25  95.35 95.12 1.0304143
# 4:             hs_fc_fetal___  92.86  94.44  92.31 91.67 1.0230744
# 5:   hs_fc_male_55yr_NeuN_neg  92.86 100.00  93.33 92.31 1.0714668
# 6:   hs_fc_male_55yr_NeuN_pos 100.00 100.00 100.00 95.45 1.0000000
# 7:    hs_fc_male_55yr_tissue_  81.82  78.57  84.03 83.33 0.9350232

gg <- ggplot(cgc_cgg, aes(x=library, y=CGC_CGG)) +
geom_bar(stat="identity", color="black", position=position_dodge(), width= 0.75, alpha = 0.5) +
theme_classic() +
ylab(expression("[mCGC/CGC]/[mCGG/CGG]")) +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 12, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black"), aspect.ratio = 5/3) +
coord_cartesian(ylim = c(0, 1.2))
ggsave("~/data/brain_cgn/fc_CGC_CGG.pdf")

rm(cgc_cgg)
rm(cgc_cgg_median)
rm(data_cgy)
```

#### extract CAN methylation

```sh
data_cay <- copy(data[grepl("......CA.....", data$context)])

# Extract CAY context
data_cay[, context_cay := as.vector(sapply(data_cay$context, function(x) paste(unlist(strsplit(x, ""))[7:9], collapse = "")))]

# remove context with N
# The like operator is a simple wrapper for grep(..., value=TRUE)
data_cay <- data_cay[!context_cay %like% "N"]

# Explore CAY context
## hs_fc_fetal___
data_cay[library == "hs_fc_fetal___", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
#    context_cay       N pct_met_median pct_met_mean
# 1:         CAC 1639142              0         1.02
# 2:         CAT 2083628              0         0.73
# 3:         CAG 3515084              0         0.69
# 4:         CAA 2217011              0         0.65

## hs_fc_female_53yr_NeuN_neg
data_cay[library == "hs_fc_female_53yr_NeuN_neg", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
#    context_cay       N pct_met_median pct_met_mean
# 1:         CAC  912872              0         3.16
# 2:         CAG 1376090              0         1.24
# 3:         CAT 1662552              0         1.20
# 4:         CAA 1814025              0         0.90

## hs_fc_female_53yr_NeuN_pos
data_cay[library == "hs_fc_female_53yr_NeuN_pos", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
#    context_cay       N pct_met_median pct_met_mean
# 1:         CAC 2253256          15.79        26.30
# 2:         CAT 2673212           3.70        10.04
# 3:         CAG 4200759           4.55         9.89
# 4:         CAA 3006805           1.82         7.74

## hs_fc_male_55yr_NeuN_neg
data_cay[library == "hs_fc_male_55yr_NeuN_neg", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
#    context_cay       N pct_met_median pct_met_mean
# 1:         CAC 1461378              0         4.07
# 2:         CAG 2165033              0         1.56
# 3:         CAT 2535517              0         1.47
# 4:         CAA 2727061              0         1.08

## hs_fc_male_55yr_NeuN_pos
data_cay[library == "hs_fc_male_55yr_NeuN_pos", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
#    context_cay       N pct_met_median pct_met_mean
# 1:         CAC 1784413          14.29        23.55
# 2:         CAG 2643565           0.00         9.33
# 3:         CAT 2966071           0.00         8.84
# 4:         CAA 3159227           0.00         6.36

## hs_fc_male_55yr_tissue_
data_cay[library == "hs_fc_male_55yr_tissue_", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
#    context_cay    N pct_met_median pct_met_mean
# 1:         CAC 3374              0         3.37
# 2:         CAG 9017              0         1.70
# 3:         CAT 5552              0         1.43
# 4:         CAA 5323              0         1.26

## hs_fc_64yr___
data_cay[library == "hs_fc_64yr___", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
#    context_cay       N pct_met_median pct_met_mean
# 1:         CAC  686876              0         7.13
# 2:         CAG 1483803              0         2.71
# 3:         CAT  861497              0         2.54
# 4:         CAA  901394              0         1.97


############ CAC/CAG ##############
cac_cag <- dcast(data_cay[, .(pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(library, context_cay)], library ~ context_cay, value.var = "pct_met_mean")
cac_cag[, CAC_CAG := CAC/CAG][order(CAC_CAG)]
#                       library  CAA   CAC  CAG   CAT  CAC_CAG
# 1:             hs_fc_fetal___ 0.65  1.02 0.69  0.73 1.478261
# 2:    hs_fc_male_55yr_tissue_ 1.26  3.37 1.70  1.43 1.982353
# 3:   hs_fc_male_55yr_NeuN_pos 6.36 23.55 9.33  8.84 2.524116
# 4: hs_fc_female_53yr_NeuN_neg 0.90  3.16 1.24  1.20 2.548387
# 5:   hs_fc_male_55yr_NeuN_neg 1.08  4.07 1.56  1.47 2.608974
# 6:              hs_fc_64yr___ 1.97  7.13 2.71  2.54 2.630996
# 7: hs_fc_female_53yr_NeuN_pos 7.74 26.30 9.89 10.04 2.659252

# set the order of libraries by ranking CAC_CAG
cac_cag$library <- factor(cac_cag$library, level = cac_cag$library[order(cac_cag$CAC_CAG)])

gg <- ggplot(cac_cag, aes(x=library, y=CAC_CAG)) +
geom_bar(stat="identity", color="black", position=position_dodge(), width= 0.75) +
theme_bw() +
ylab(expression("[mCAC/CAC]/[mCAG/CAG]")) +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 12, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black"), aspect.ratio = 5/3) +
coord_cartesian(ylim = c(0, 3))
ggsave("~/data/brain_can/fc_CAC_CAG.pdf")


########### plot CAN changes through age ###########
cac_cag_melt <- melt(cac_cag, id.vars = "library", variable.name = "CAN", value.name = "pct_met")

# set the order of libraries by ranking CAC_CAG
cac_cag_melt$library <- factor(cac_cag_melt$library, level = c("hs_fc_fetal___", "hs_fc_female_53yr_NeuN_neg", "hs_fc_female_53yr_NeuN_pos", "hs_fc_male_55yr_NeuN_neg", "hs_fc_male_55yr_NeuN_pos", "hs_fc_male_55yr_tissue_", "hs_fc_64yr___"))

write.table(cac_cag_melt, file = "~/data/brain_can/fc_CAC_CAG.txt", row.names = FALSE, quote = FALSE, sep = "\t")

gg <- ggplot(cac_cag_melt, aes(x=library, y= pct_met, color = CAN, group = CAN)) +
geom_line(aes(linetype = CAN), lwd = 1) +
scale_linetype_manual(values = c(1,1,1,1,2)) +
theme_bw() +
ylab(expression("[mCAC/CAC]/[mCAG/CAG]")) + 
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 12, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black"), aspect.ratio = 5/3) +
coord_cartesian(ylim = c(0, 30))
ggsave("~/data/brain_can/fc_CAC_CAG_age.pdf")

rm(cac_cag_melt)
rm(cac_cag)
rm(data_cay)
rm(data)
```