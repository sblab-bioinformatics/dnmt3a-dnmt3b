
## download dataset from GEO

Global epigenomic reconfiguration during mammalian brain development. Science 2013 Aug 9;341(6146):1237905.
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47966>

```sh
cd ~
mkdir -p 20190722_mouse_brain_development/data
cd ~/20190722_mouse_brain_development/data

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173779/suppl/GSM1173779_allC.MethylC-Seq_mm_fc_fetal.chr19.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173780/suppl/GSM1173780_allC.MethylC-Seq_mm_fc_1wk.chr19.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173781/suppl/GSM1173781_allC.MethylC-Seq_mm_fc_2wk.chr19.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173782/suppl/GSM1173782_allC.MethylC-Seq_mm_fc_4wk.chr19.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173783/suppl/GSM1173783_allC.MethylC-Seq_mm_fc_6wk.chr19.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173784/suppl/GSM1173784_allC.MethylC-Seq_mm_fc_10wk.chr19.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173785/suppl/GSM1173785_allC.MethylC-Seq_mm_fc_22mo.chr19.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173786/suppl/GSM1173786_allC.MethylC-Seq_mm_fc_male_7wk_neun_pos.chr19.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173787/suppl/GSM1173787_allC.MethylC-Seq_mm_fc_male_7wk_neun_neg.chr19.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173788/suppl/GSM1173788_allC.MethylC-Seq_mm_fc_female_6wk_neun_pos.chr19.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173789/suppl/GSM1173789_allC.MethylC-Seq_mm_fc_female_6wk_neun_neg.chr19.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173790/suppl/GSM1173790_allC.MethylC-Seq_mm_fc_female_12mo_neun_pos.chr19.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173791/suppl/GSM1173791_allC.MethylC-Seq_mm_fc_female_12mo_neun_neg.chr19.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173792/suppl/GSM1173792_allC.MethylC-Seq_mm_fc_glia_S100b_pos.chr19.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173793/suppl/GSM1173793_allC.MethylC-Seq_mm_fc_tet2ko.chr19.txt.gz &

```

### prepare mm9 reference and methylation context

```sh
srun --mem 96G --pty /usr/bin/bash

## Create mm9.chr19 reference
cd ~/20190722_mouse_brain_development
mkdir reference
cd reference

ref=genome.fa

samtools faidx $ref chr19 > mm9.chr19.fa
samtools faidx mm9.chr19.fa

ref=mm9.chr19.fa

fastaRegexFinder.py -f $ref -r C | \
cut -f1-3,6 | \
awk -v OFS="\t" '{print $1, $2, $3, ".", ".", $4}' | \
bedtools slop -i - -g $ref.fai -b 6 | \
bedtools getfasta -fi $ref -bed - -bedOut -s | \
awk 'length($7) == 13' | \
awk -v OFS="\t" '{print $1, $2+6, $3-6, $4, $5, $6, toupper($7)}' > mm9.chr19.context.bed &

wc -l mm9.chr19.context.bed
# 24845749 mm9.chr19.context.bed (C±6)
# 24845751 mm9.chr19.context.bed (C±3, overwritten)
```

## extract context information around C (±6)

```sh
cd ~/20190722_mouse_brain_development/data

bed=~/20190722_mouse_brain_development/reference/mm9.chr19.context.bed

for report in *.chr19.txt.gz
do
  bname=${report%.txt.gz}
  # echo $bname $report
  zcat "$report" | \
  awk -v OFS='\t' '{if ($6 > 9) print "chr19", $2, $2+1, $5, $6, $3}' | \
  bedtools intersect -a - -b "$bed" -loj -s -sorted | \
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $13}' > $bname.context.txt &
done

wc -l *.context.txt 
#    5251961 GSM1173779_allC.MethylC-Seq_mm_fc_fetal.context.txt
#    5241652 GSM1173780_allC.MethylC-Seq_mm_fc_1wk.context.txt
#    9740227 GSM1173781_allC.MethylC-Seq_mm_fc_2wk.context.txt
#   10579190 GSM1173782_allC.MethylC-Seq_mm_fc_4wk.context.txt
#    9977493 GSM1173783_allC.MethylC-Seq_mm_fc_6wk.context.txt
#    7002418 GSM1173784_allC.MethylC-Seq_mm_fc_10wk.context.txt
#   11869853 GSM1173785_allC.MethylC-Seq_mm_fc_22mo.context.txt
#    8396847 GSM1173786_allC.MethylC-Seq_mm_fc_male_7wk_neun_pos.context.txt
#   11556716 GSM1173787_allC.MethylC-Seq_mm_fc_male_7wk_neun_neg.context.txt
#    3115075 GSM1173788_allC.MethylC-Seq_mm_fc_female_6wk_neun_pos.context.txt
#     217580 GSM1173789_allC.MethylC-Seq_mm_fc_female_6wk_neun_neg.context.txt
#    5006539 GSM1173790_allC.MethylC-Seq_mm_fc_female_12mo_neun_pos.context.txt
#    2785635 GSM1173791_allC.MethylC-Seq_mm_fc_female_12mo_neun_neg.context.txt
#   12282691 GSM1173792_allC.MethylC-Seq_mm_fc_glia_S100b_pos.context.txt
#   12751741 GSM1173793_allC.MethylC-Seq_mm_fc_tet2ko.context.txt
#     547297 GSM1173794_allC.TAB-Seq_mm_fc_fetal.context.txt
#    5976088 GSM1173795_allC.TAB-Seq_mm_fc_6wk.context.txt

## Error occured processing GSM1173794_allC.TAB-Seq_mm_fc_fetal.chr19.txt.gz
## Error: Sorted input specified, but the file - has the following out of order record
## chr19   8387233 8387234 0       20      +

zgrep -C 3 -w "8387233" GSM1173794_allC.TAB-Seq_mm_fc_fetal.chr19.txt.gz
# 19	8387230	-	CCC	0	3
# 19	8387232	+	CCC	0	20
# 19	8387234	+	CNN	0	20
# 19	8387233	+	CNN	0	20 <<<*
# 19	8387236	-	CTN	0	3
# 19	8387243	-	CAA	0	3
# 19	8387244	-	CCA	0	2

# reprocess GSM1173794_allC.TAB-Seq_mm_fc_fetal.chr19.txt.gz in unsorted mode
for report in GSM1173794_allC.TAB-Seq_mm_fc_fetal.chr19.txt.gz
do
  bname=${report%.txt.gz}
  # echo $bname $report
  zcat "$report" | \
  sort -k2,2n | \
  awk -v OFS='\t' '{if ($6 > 9) print "chr19", $2, $2+1, $5, $6, $3}' | \
  bedtools intersect -a - -b "$bed" -loj -s -sorted | \
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $13}' > $bname.context.txt &
done

# reprocess GSM1173789_allC.MethylC-Seq_mm_fc_female_6wk_neun_neg.chr19.txt.gz
for report in GSM1173789_allC.MethylC-Seq_mm_fc_female_6wk_neun_neg.chr19.txt.gz
do
  bname=${report%.txt.gz}
  # echo $bname $report
  zcat "$report" | \
  awk -v OFS='\t' '{if ($6 > 5) print "chr19", $2, $2+1, $5, $6, $3}' | \
  bedtools intersect -a - -b "$bed" -loj -s -sorted | \
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $13}' > $bname.context.txt &
done

wc -l *.context.txt 
#     5251961 GSM1173779_allC.MethylC-Seq_mm_fc_fetal.chr19.context.txt
#     5241652 GSM1173780_allC.MethylC-Seq_mm_fc_1wk.chr19.context.txt
#     9740227 GSM1173781_allC.MethylC-Seq_mm_fc_2wk.chr19.context.txt
#    10579190 GSM1173782_allC.MethylC-Seq_mm_fc_4wk.chr19.context.txt
#     9977493 GSM1173783_allC.MethylC-Seq_mm_fc_6wk.chr19.context.txt
#     7002418 GSM1173784_allC.MethylC-Seq_mm_fc_10wk.chr19.context.txt
#    11869853 GSM1173785_allC.MethylC-Seq_mm_fc_22mo.chr19.context.txt
#     8396847 GSM1173786_allC.MethylC-Seq_mm_fc_male_7wk_neun_pos.chr19.context.txt
#    11556716 GSM1173787_allC.MethylC-Seq_mm_fc_male_7wk_neun_neg.chr19.context.txt
#     3115075 GSM1173788_allC.MethylC-Seq_mm_fc_female_6wk_neun_pos.chr19.context.txt
#     3078045 GSM1173789_allC.MethylC-Seq_mm_fc_female_6wk_neun_neg.chr19.context.txt
#     5006539 GSM1173790_allC.MethylC-Seq_mm_fc_female_12mo_neun_pos.chr19.context.txt
#     2785635 GSM1173791_allC.MethylC-Seq_mm_fc_female_12mo_neun_neg.chr19.context.txt
#    12282691 GSM1173792_allC.MethylC-Seq_mm_fc_glia_S100b_pos.chr19.context.txt
#    12751741 GSM1173793_allC.MethylC-Seq_mm_fc_tet2ko.chr19.context.txt
#     5954244 GSM1173794_allC.TAB-Seq_mm_fc_fetal.chr19.context.txt
#     5976088 GSM1173795_allC.TAB-Seq_mm_fc_6wk.chr19.context.txt

```

# look at NCAN context

### load MethylC-Seq data in R

```R
R
library(data.table)
library(ggplot2)

# Set width
options(width = 300)

# Load MethylC-Seq data
data <- fread("~/20190722_mouse_brain_development/data/mouse_brain_MethylC-Seq.chr19.context.txt")
setnames(data, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand", "context", "library"))

# Calculate percentage methylation
data[, pct_met := round(100 * cnt_met/(cnt_tot), 2)]

table(data$library)
#                 fc_10wk                  fc_1wk                 fc_22mo                  fc_2wk                  fc_4wk                  fc_6wk fc_female_12mo_neun_neg fc_female_12mo_neun_pos  fc_female_6wk_neun_neg  fc_female_6wk_neun_pos                fc_fetal       fc_glia_S100b_pos
#                 7002418                 5241652                11869853                 9740227                10579190                 9977493                 2785635                 5006539                 3078045                 3115075                 5251961                12282691
#    fc_male_7wk_neun_neg    fc_male_7wk_neun_pos               fc_tet2ko
#                11556716                 8396847                12751741
```


# extract CAN methylation

```r
data_ncan <- fread("~/20190722_mouse_brain_development/data/20190723_mouse_brain_chr19_NCAN.txt")

# Extract cay context
data_ncan[, context_cay := as.vector(sapply(data_ncan$context_xcay, function(x) paste(unlist(strsplit(x, ""))[2:4], collapse = "")))]

data_cay <- data_ncan[, .(.N, pct_met_median = median(pct_met_median, na.rm=TRUE), pct_met_mean = round(mean(pct_met_mean, na.rm=TRUE), 2)), by = .(library, context_cay)]

# Explore cay context
data_cay[library == "fc_22mo", ][order(-pct_met_mean)]
#    library context_cay N pct_met_median pct_met_mean
# 1: fc_22mo         CAC 4                     0              7.27
# 2: fc_22mo         CAG 4                     0              2.51
# 3: fc_22mo         CAT 4                     0              2.50
# 4: fc_22mo         CAA 4                     0              1.85

data_cay[library == "fc_fetal", ][order(-pct_met_mean)]
#     library context_cay N pct_met_median pct_met_mean
# 1: fc_fetal         CAC 4                     0              0.86
# 2: fc_fetal         CAG 4                     0              0.64
# 3: fc_fetal         CAT 4                     0              0.61
# 4: fc_fetal         CAA 4                     0              0.61


# set the order of libraries
data_cay$library <- factor(data_cay$library, levels = c("fc_fetal", "fc_1wk", "fc_2wk", "fc_4wk", "fc_6wk","fc_10wk", "fc_22mo", "fc_male_7wk_neun_pos", "fc_male_7wk_neun_neg",  "fc_female_6wk_neun_pos", "fc_female_6wk_neun_neg", "fc_female_12mo_neun_pos", "fc_female_12mo_neun_neg", "fc_glia_S100b_pos", "fc_tet2ko"))

data_cay$context_cay <- factor(data_cay$context_cay, levels = data_cay[, .(.N, pct_met_mean = round(mean(pct_met_mean, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]$context_cay)
## levels: "CAC" "CAG" "CAT" "CAA"

# plot
gg <- ggplot(data_cay, aes(x = library, y = pct_met_mean, color = context_cay, group = context_cay)) +
geom_point(size = 1) +
geom_line() +
theme_bw() +
ylab(expression("mCAN/CAN")) + 
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=12), axis.text.y = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, size = 10, color = "black", hjust = 1), legend.text = element_text(size = 12, color = "black")) +
coord_cartesian(ylim = c(0, 15))
ggsave("~/figures/CAN_age.pdf", width = 25, height = 20, units = "cm")
```

### neun_pos & neun_neg & tet2ko  

```r
data_ncan <- fread("~/20190722_mouse_brain_development/data/20190723_mouse_brain_chr19_NCAN.txt")

# Extract cay context
data_ncan[, context_cay := as.vector(sapply(data_ncan$context_xcay, function(x) paste(unlist(strsplit(x, ""))[2:4], collapse = "")))]

data_cay <- data_ncan[, .(.N, pct_met_median = median(pct_met_median, na.rm=TRUE), pct_met_mean = round(mean(pct_met_mean, na.rm=TRUE), 2)), by = .(library, context_cay)]

# extract cell-type-specific rows
can_celltype <- data_cay[grepl("neun_pos|neun_neg|fc_tet2ko|fc_glia_S100b_pos", data_cay$library), ]

table(can_celltype$library)

can_celltype$library <- factor(can_celltype$library, levels = c("fc_male_7wk_neun_pos", "fc_male_7wk_neun_neg",  "fc_female_6wk_neun_pos", "fc_female_6wk_neun_neg", "fc_female_12mo_neun_pos", "fc_female_12mo_neun_neg", "fc_glia_S100b_pos", "fc_tet2ko"))

can_celltype$context_cay <- factor(can_celltype$context_cay, levels = data_cay[, .(.N, pct_met_mean = round(mean(pct_met_mean, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]$context_cay)
## levels: "CAC" "CAG" "CAT" "CAA"

# plot point and line
gg <- ggplot(can_celltype, aes(x = library, y = pct_met_mean, color = context_cay, group = context_cay)) +
geom_point(size = 1.5) +
geom_line() +
theme_bw() +
ylab(expression("mCAN/CAN")) + 
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=12), axis.text.y = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, size = 10, color = "black", hjust = 1), legend.text = element_text(size = 12, color = "black")) +
coord_cartesian(ylim = c(0, 15))
ggsave("~/figures/CAN_age_celltype_line.pdf", width = 25, height = 20, units = "cm")


# plot bar
gg <- ggplot(can_celltype, aes(x = library, y = pct_met_mean, color = context_cay, group = context_cay)) +
geom_bar(aes(fill = context_cay), stat="identity", position=position_dodge()) +
theme_bw() +
ylab(expression("mCAN/CAN")) + 
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=12), axis.text.y = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, size = 10, color = "black", hjust = 1), legend.text = element_text(size = 12, color = "black")) +
coord_cartesian(ylim = c(0, 15))
ggsave("~/figures/CAN_age_celltype_bar.pdf", width = 25, height = 20, units = "cm")

# plot stacked bar
gg <- ggplot(can_celltype, aes(x = library, y = pct_met_mean, colour = context_cay, group = context_cay)) +
geom_bar(aes(fill = context_cay), stat="identity") +
geom_text(aes(label = paste0(pct_met_mean,"%")), size=4, colour = "black", position = position_stack(vjust=0.5)) +
theme_bw() +
ylab(expression("mCAN/CAN")) + 
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=12), axis.text.y = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, size = 10, color = "black", hjust = 1), legend.text = element_text(size = 12, color = "black")) +
coord_cartesian(ylim = c(0, 25))
ggsave("~/figures/CAN_age_celltype_stackbar.pdf", width = 25, height = 20, units = "cm")

######################### plot Percent stacked bar ######################### 
# amounts of CAA,CAC,CAG,CAT on both strands of mm10.chr19 are in similar range
# See: /Users/mao01/Google Drive/PROGRAM/20190724_Trinucleotide_Frenquency
# Following plot is in approximation assuming that: CAA ~ CAC ~ CAG ~ CAT

# plot Percent stacked bar
gg <- ggplot(can_celltype, aes(x = library, y = pct_met_mean, colour = context_cay, group = context_cay)) +
geom_bar(aes(fill = context_cay), stat="identity", position = "fill") +
geom_text(aes(label = paste0(pct_met_mean,"%")), size=4, colour = "black", position = position_fill(vjust=0.5)) +
theme_bw() +
ylab(expression("[mCAN/CAN]/sum([mCAN/CAN])")) + 
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=12), axis.text.y = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, size = 10, color = "black", hjust = 1), legend.text = element_text(size = 12, color = "black")) +
coord_cartesian(ylim = c(0, 1))
ggsave("~/figures/CAN_age_celltype_stackbar_percent.pdf", width = 25, height = 20, units = "cm")
```


### plot CAC/CAG


```r
data_ncan <- fread("~/20190722_mouse_brain_development/data/20190723_mouse_brain_chr19_NCAN.txt")

# Extract cay context
data_ncan[, context_cay := as.vector(sapply(data_ncan$context_xcay, function(x) paste(unlist(strsplit(x, ""))[2:4], collapse = "")))]

data_cay <- data_ncan[, .(.N, pct_met_median = median(pct_met_median, na.rm=TRUE), pct_met_mean = round(mean(pct_met_mean, na.rm=TRUE), 2)), by = .(library, context_cay)]

cac_cag <- dcast(data_cay, library ~ context_cay, value.var = "pct_met_mean")
cac_cag[, `:=`(CAC_CAG = CAC/CAG, CAY_CAR = (CAC + CAT)/(CAG + CAA))]
cac_cag[order(CAC_CAG)]
#                     library  CAA   CAC  CAG  CAT  CAC_CAG  CAY_CAR
#  1:                fc_fetal 0.61  0.86 0.64 0.61 1.343750 1.176000
#  2:                  fc_1wk 0.68  1.14 0.74 0.73 1.540541 1.316901
#  3:  fc_female_6wk_neun_neg 0.82  2.84 1.13 1.07 2.513274 2.005128
#  4:                  fc_2wk 1.20  3.70 1.46 1.47 2.534247 1.943609
#  5:    fc_male_7wk_neun_neg 0.93  3.48 1.32 1.27 2.636364 2.111111
#  6: fc_female_12mo_neun_neg 0.85  3.17 1.19 1.14 2.663866 2.112745
#  7:                  fc_4wk 2.29  8.01 2.88 2.98 2.781250 2.125725
#  8:       fc_glia_S100b_pos 0.91  3.31 1.17 1.17 2.829060 2.153846
#  9:                 fc_22mo 1.85  7.27 2.51 2.50 2.896414 2.240826
# 10:                  fc_6wk 2.09  7.96 2.74 2.79 2.905109 2.225673
# 11:                 fc_10wk 1.67  6.80 2.30 2.23 2.956522 2.274559
# 12: fc_female_12mo_neun_pos 2.51 11.29 3.78 3.72 2.986772 2.386328
# 13:    fc_male_7wk_neun_pos 2.52 11.37 3.80 3.79 2.992105 2.398734
# 14:  fc_female_6wk_neun_pos 2.51 11.42 3.79 3.72 3.013193 2.403175
# 15:               fc_tet2ko 1.55  6.91 2.20 2.23 3.140909 2.437333

# transform from wide to long format
cac_cag_melt <- melt(cac_cag[, c(1,6,7)], value.var = "library")

# set levels
cac_cag_melt$library <- factor(cac_cag_melt$library, levels = c("fc_fetal", "fc_1wk", "fc_2wk", "fc_4wk", "fc_6wk","fc_female_6wk_neun_pos", "fc_female_6wk_neun_neg", "fc_male_7wk_neun_pos", "fc_male_7wk_neun_neg", "fc_10wk", "fc_female_12mo_neun_pos", "fc_female_12mo_neun_neg", "fc_22mo", "fc_glia_S100b_pos", "fc_tet2ko"))

## plot CAC_CAG, CAY_CAR side by side, ranked by age
gg <- ggplot(cac_cag_melt, aes(x = library, y = value, fill = variable)) +
geom_bar(stat="identity", color="black", position=position_dodge(), width= 0.6, alpha = 0.5) +
theme_bw() +
ylab(expression("[mCAC/CAC]/[mCAG/CAG]")) +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 12, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
coord_cartesian(ylim = c(0, 3.5))
ggsave("~/figures/CAC_CAG_CAY_CAR_age.pdf")
```
