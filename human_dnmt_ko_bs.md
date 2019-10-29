
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
