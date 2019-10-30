
## download dataset from GEO

Global epigenomic reconfiguration during mammalian brain development. Science 2013 Aug 9;341(6146):1237905.
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47966>

```sh
cd ~
mkdir data
cd ~/data

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173779/suppl/GSM1173779_allC.MethylC-Seq_mm_fc_fetal.chr1.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173780/suppl/GSM1173780_allC.MethylC-Seq_mm_fc_1wk.chr1.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173781/suppl/GSM1173781_allC.MethylC-Seq_mm_fc_2wk.chr1.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173782/suppl/GSM1173782_allC.MethylC-Seq_mm_fc_4wk.chr1.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173783/suppl/GSM1173783_allC.MethylC-Seq_mm_fc_6wk.chr1.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173784/suppl/GSM1173784_allC.MethylC-Seq_mm_fc_10wk.chr1.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173785/suppl/GSM1173785_allC.MethylC-Seq_mm_fc_22mo.chr1.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173786/suppl/GSM1173786_allC.MethylC-Seq_mm_fc_male_7wk_neun_pos.chr1.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173787/suppl/GSM1173787_allC.MethylC-Seq_mm_fc_male_7wk_neun_neg.chr1.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173788/suppl/GSM1173788_allC.MethylC-Seq_mm_fc_female_6wk_neun_pos.chr1.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173789/suppl/GSM1173789_allC.MethylC-Seq_mm_fc_female_6wk_neun_neg.chr1.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173790/suppl/GSM1173790_allC.MethylC-Seq_mm_fc_female_12mo_neun_pos.chr1.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173791/suppl/GSM1173791_allC.MethylC-Seq_mm_fc_female_12mo_neun_neg.chr1.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173792/suppl/GSM1173792_allC.MethylC-Seq_mm_fc_glia_S100b_pos.chr1.txt.gz &

curl -sSOL ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1173nnn/GSM1173793/suppl/GSM1173793_allC.MethylC-Seq_mm_fc_tet2ko.chr1.txt.gz &

```

### prepare mm9 reference and methylation context

```sh
srun --mem 16G --pty /usr/bin/bash

## Create mm9.chr1 reference
cd ~
mkdir reference
cd reference

ref=genome.fa

samtools faidx $ref chr1 > mm9.chr1.fa
samtools faidx mm9.chr1.fa

ref=mm9.chr1.fa

fastaRegexFinder.py -f $ref -r C | \
cut -f1-3,6 | \
awk -v OFS="\t" '{print $1, $2, $3, ".", ".", $4}' | \
bedtools slop -i - -g $ref.fai -b 6 | \
bedtools getfasta -fi $ref -bed - -bedOut -s | \
awk 'length($7) == 13' | \
awk -v OFS="\t" '{print $1, $2+6, $3-6, $4, $5, $6, toupper($7)}' > mm9.chr1.context.bed &
```

## extract context information around C (±6)

```sh
cd ~/data

bed=~/reference/mm9.chr1.context.bed

for report in *.chr1.txt.gz
do
  bname=${report%.txt.gz}
  # echo $bname $report
  zcat "$report" | \
  awk -v OFS='\t' '{if ($6 > 9) print "chr1", $2, $2+1, $5, $6, $3}' | \
  bedtools intersect -a - -b "$bed" -loj -s -sorted | \
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $13}' > $bname.context.txt &
done

## Error occured processing GSM1173794_allC.TAB-Seq_mm_fc_fetal.chr1.txt.gz
## Error: Sorted input specified, but the file - has the following out of order record
## chr1   8387233 8387234 0       20      +

zgrep -C 3 -w "8387233" GSM1173794_allC.TAB-Seq_mm_fc_fetal.chr1.txt.gz
# 19	8387230	-	CCC	0	3
# 19	8387232	+	CCC	0	20
# 19	8387234	+	CNN	0	20
# 19	8387233	+	CNN	0	20 <<<*
# 19	8387236	-	CTN	0	3
# 19	8387243	-	CAA	0	3
# 19	8387244	-	CCA	0	2

# reprocess GSM1173794_allC.TAB-Seq_mm_fc_fetal.chr1.txt.gz in unsorted mode
for report in GSM1173794_allC.TAB-Seq_mm_fc_fetal.chr1.txt.gz
do
  bname=${report%.txt.gz}
  # echo $bname $report
  zcat "$report" | \
  sort -k2,2n | \
  awk -v OFS='\t' '{if ($6 > 9) print "chr1", $2, $2+1, $5, $6, $3}' | \
  bedtools intersect -a - -b "$bed" -loj -s -sorted | \
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $13}' > $bname.context.txt &
done

# reprocess GSM1173789_allC.MethylC-Seq_mm_fc_female_6wk_neun_neg.chr1.txt.gz
for report in GSM1173789_allC.MethylC-Seq_mm_fc_female_6wk_neun_neg.chr1.txt.gz
do
  bname=${report%.txt.gz}
  # echo $bname $report
  zcat "$report" | \
  awk -v OFS='\t' '{if ($6 > 5) print "chr1", $2, $2+1, $5, $6, $3}' | \
  bedtools intersect -a - -b "$bed" -loj -s -sorted | \
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $13}' > $bname.context.txt &
done

# Concatenate MethylC-Seq context files
tableCat.py -i GSM*MethylC-Seq*.context.txt -r .context.txt | awk '{split($8,a,"."); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"substr(a[2],16)}' > mouse_brain_MethylC-Seq.chr1.context.txt &

```


### load MethylC-Seq data in R

```r
R
library(data.table)
library(ggplot2)

# Set width
options(width = 300)

# Load MethylC-Seq data
data <- fread("~/data/mouse_brain_MethylC-Seq.chr1.context.txt")
setnames(data, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand", "context", "library"))

# Calculate percentage methylation
data[, pct_met := round(100 * cnt_met/(cnt_tot), 2)]

table(data$library)
#                 fc_10wk                  fc_1wk                 fc_22mo                  fc_2wk                  fc_4wk                  fc_6wk fc_female_12mo_neun_neg fc_female_12mo_neun_pos  fc_female_6wk_neun_neg  fc_female_6wk_neun_pos                fc_fetal       fc_glia_S100b_pos
#                 7002418                 5241652                11869853                 9740227                10579190                 9977493                 2785635                 5006539                 3078045                 3115075                 5251961                12282691
#    fc_male_7wk_neun_neg    fc_male_7wk_neun_pos               fc_tet2ko
#                11556716                 8396847                12751741

```

### extract NCAN methylation

```r
data_xcay <- copy(data[grepl("...CA..", data$context)])

# Extract XCAY context
data_xcay[, context_xcay := as.vector(sapply(data_xcay$context, function(x) paste(unlist(strsplit(x, ""))[3:6], collapse = "")))]

# remove context with N
data_xcay <- data_xcay[!context_xcay %like% "N"]

# Explore XCAY context
## hs_mfg_35do
data_xcay[library == "fc_female_6wk_neun_pos", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_xcay)][order(-pct_met_mean)]
#     context_xcay      N pct_met_median pct_met_mean
#  1:         ACAC  83627           7.69        14.04
#  2:         GCAC  39201           0.00        11.80
#  3:         CCAC  43212           0.00        10.39
#  4:         TCAC  57923           0.00         9.46
#  5:         ACAT 121200           0.00         5.28
#  6:         ACAG 113703           0.00         5.04
#  7:         CCAG  69927           0.00         4.18
#  8:         GCAT  65551           0.00         3.56
#  9:         CCAT  65718           0.00         3.30
# 10:         GCAG  65265           0.00         3.28
# 11:         ACAA 146929           0.00         2.99
# 12:         CCAA  77883           0.00         2.94
# 13:         TCAT  96305           0.00         2.75
# 14:         TCAG  88321           0.00         2.67
# 15:         GCAA  73170           0.00         2.44
# 16:         TCAA 106481           0.00         1.66

data_xcay[library == "fc_female_6wk_neun_neg", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_xcay)][order(-pct_met_mean)]
#     context_xcay      N pct_met_median pct_met_mean
#  1:         ACAC  82485              0         3.31
#  2:         GCAC  38800              0         2.94
#  3:         CCAC  44621              0         2.84
#  4:         TCAC  56699              0         2.25
#  5:         CCAG  71849              0         1.41
#  6:         ACAT 115594              0         1.29
#  7:         ACAG 110897              0         1.27
#  8:         CCAT  64548              0         1.18
#  9:         CCAA  77722              0         1.06
# 10:         GCAG  64700              0         1.00
# 11:         GCAT  62830              0         0.97
# 12:         TCAT  91151              0         0.85
# 13:         ACAA 143508              0         0.84
# 14:         TCAG  87027              0         0.83
# 15:         GCAA  71617              0         0.76
# 16:         TCAA 102889              0         0.64

data_xcay[library == "fc_glia_S100b_pos", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_xcay)][order(-pct_met_mean)]
#     context_xcay      N pct_met_median pct_met_mean
#  1:         ACAC 202200              0         4.13
#  2:         CCAC 191513              0         3.36
#  3:         GCAC 194909              0         3.04
#  4:         TCAC 197143              0         2.72
#  5:         CCAG 389346              0         1.45
#  6:         ACAT 238559              0         1.44
#  7:         ACAG 401350              0         1.39
#  8:         CCAT 228700              0         1.31
#  9:         CCAA 213828              0         1.22
# 10:         GCAT 241180              0         1.00
# 11:         ACAA 244744              0         0.97
# 12:         GCAG 384087              0         0.96
# 13:         TCAT 239173              0         0.93
# 14:         TCAG 383810              0         0.89
# 15:         GCAA 238627              0         0.76
# 16:         TCAA 230140              0         0.70

data_xcay[library == "fc_10wk", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_xcay)][order(-pct_met_mean)]
#     context_xcay      N pct_met_median pct_met_mean
#  1:         ACAC 123250              0         8.47
#  2:         GCAC 111177              0         6.66
#  3:         CCAC  97357              0         6.36
#  4:         TCAC 113510              0         5.69
#  5:         ACAT 159171              0         2.99
#  6:         ACAG 252741              0         2.93
#  7:         CCAG 219311              0         2.68
#  8:         CCAT 135679              0         2.16
#  9:         CCAA 136279              0         2.08
# 10:         GCAT 152239              0         2.06
# 11:         ACAA 168767              0         1.93
# 12:         GCAG 221975              0         1.90
# 13:         TCAT 153769              0         1.71
# 14:         TCAG 233682              0         1.68
# 15:         GCAA 153658              0         1.50
# 16:         TCAA 154335              0         1.17

data_xcay[library == "fc_tet2ko", .(.N, pct_met_median = median(pct_met, na.rm=TRUE), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_xcay)][order(-pct_met_mean)]
#     context_xcay      N pct_met_median pct_met_mean
#  1:         ACAC 259235              0         8.72
#  2:         GCAC 189461              0         6.77
#  3:         CCAC 204200              0         6.30
#  4:         TCAC 228004              0         5.84
#  5:         ACAT 329817              0         3.15
#  6:         ACAG 418428              0         2.92
#  7:         CCAG 351499              0         2.48
#  8:         GCAT 260217              0         2.06
#  9:         CCAT 266753              0         2.02
# 10:         ACAA 326498              0         1.92
# 11:         GCAG 335665              0         1.85
# 12:         CCAA 245618              0         1.84
# 13:         TCAT 308822              0         1.67
# 14:         TCAG 379586              0         1.56
# 15:         GCAA 250634              0         1.39
# 16:         TCAA 287240              0         1.04

# save data_xcay table
fwrite(data_xcay[, .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(library, context_xcay)], file = "~/data/mouse_brain_chr1_NCAN.txt", sep = "\t", row.names=FALSE, quote=FALSE)

rm(data_xcay)
```

# extract CAN methylation

```r
data_ncan <- fread("~/data/mouse_brain_chr1_NCAN.txt")


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


### plot CAC/CAG


```r
data_ncan <- fread("~/data/mouse_brain_chr1_NCAN.txt")

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
