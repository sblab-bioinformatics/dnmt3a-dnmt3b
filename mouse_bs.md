
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

nohup bismark_genome_preparation . &
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


#### Deduplication

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
cd /scratchb/sblab/martin03/repository/20181213_DNMT_preference/data/20190624/bismark

# TKO
for bam in TKO*.deduplicated.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o $bname.log --mem 64G --wrap "samtools sort -@ 20 -T /scratchb/sblab/martin03/tmp/$bname -o $bname.coordsort.bam $bam && \
  samtools index $bname.coordsort.bam && \
  samtools view -@ 20 $bname.coordsort.bam chr1 -b -o $bname.coordsort.chr1.bam && \
  samtools index $bname.coordsort.chr1.bam && \
  samtools sort -@ 20 -T /scratchb/sblab/martin03/tmp/$bname.coordsort.chr1 -o $bname.coordsort.chr1.namesort.bam $bname.coordsort.chr1.bam -n && \
  rm $bname.coordsort.bam* $bname.coordsort.chr1.bam*"
done

tail TKO*.deduplicated.log


# WT
for bam in MouseES_BisSeq_HiSeq*.deduplicated.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o $bname.log --mem 64G --wrap "samtools sort -@ 20 -T /scratchb/sblab/martin03/tmp/$bname -o $bname.coordsort.bam $bam && \
  samtools index $bname.coordsort.bam && \
  samtools view -@ 20 $bname.coordsort.bam chr1 -b -o $bname.coordsort.chr1.bam && \
  samtools index $bname.coordsort.chr1.bam && \
  samtools sort -@ 20 -T /scratchb/sblab/martin03/tmp/$bname.coordsort.chr1 -o $bname.coordsort.chr1.namesort.bam $bname.coordsort.chr1.bam -n && \
  rm $bname.coordsort.bam* $bname.coordsort.chr1.bam*"
done

tail MouseES_BisSeq_HiSeq*.deduplicated.log
```





## Methylation chr1

```bash
cd /scratchb/sblab/martin03/repository/20181213_DNMT_preference/data/20190624/bismark

mkdir ../methylation

ref=../reference

# TKO
for bam in TKO*.deduplicated.coordsort.chr1.namesort.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o ../methylation/$bname.log --mem 64GB --wrap "bismark_methylation_extractor -p --comprehensive -o ../methylation --gzip --parallel 8 --bedGraph -CX_context --cytosine_report --genome_folder $ref $bam"
done

tail ../methylation/TKO*.deduplicated.coordsort.chr1.namesort.log

grep "C methylated in CpG context" ../methylation/TKO*.log
# ../methylation/TKO_DNMT3A2_bisulfite_hiseq2000_SRR8139275_1_bismark_bt2_pe.deduplicated.coordsort.chr1.namesort.log:C methylated in CpG context:	7.1%
# ../methylation/TKO_DNMT3B1_bisulfite_hiseq2000_SRR8139276_1_bismark_bt2_pe.deduplicated.coordsort.chr1.namesort.log:C methylated in CpG context:	2.8%

grep "C methylated in CH. context" ../methylation/TKO*.log
# ../methylation/TKO_DNMT3A2_bisulfite_hiseq2000_SRR8139275_1_bismark_bt2_pe.deduplicated.coordsort.chr1.namesort.log:C methylated in CHG context:	0.1%
# ../methylation/TKO_DNMT3A2_bisulfite_hiseq2000_SRR8139275_1_bismark_bt2_pe.deduplicated.coordsort.chr1.namesort.log:C methylated in CHH context:	0.1%
# ../methylation/TKO_DNMT3B1_bisulfite_hiseq2000_SRR8139276_1_bismark_bt2_pe.deduplicated.coordsort.chr1.namesort.log:C methylated in CHG context:	0.2%
# ../methylation/TKO_DNMT3B1_bisulfite_hiseq2000_SRR8139276_1_bismark_bt2_pe.deduplicated.coordsort.chr1.namesort.log:C methylated in CHH context:	0.1%


# WT
for bam in MouseES_BisSeq_HiSeq*.deduplicated.coordsort.chr1.namesort.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o ../methylation/$bname.log --mem 64GB --wrap "bismark_methylation_extractor -s --comprehensive -o ../methylation --gzip --parallel 8 --bedGraph -CX_context --cytosine_report --genome_folder $ref $bam"
done

tail ../methylation/MouseES_BisSeq_HiSeq*.deduplicated.coordsort.chr1.namesort.log

grep "C methylated in CpG context" ../methylation/MouseES_BisSeq_HiSeq*.log
# C methylated in CpG context:	81.7%

grep "C methylated in CH. context" ../methylation/MouseES_BisSeq_HiSeq*.log
# C methylated in CHG context:	2.5%
# C methylated in CHH context:	2.3%
```





## Sequence context plots chr1

### Obtain sequence context and methylation files

```bash
srun --mem 96G --pty /usr/bin/bash


## Create chr1 reference
cd /scratcha/sblab/martin03/repository/20181213_DNMT_preference/data/20190624
mkdir reference2 && cd reference2

samtools faidx ../reference/GRCm38.p6.genome.fa
samtools faidx ../reference/GRCm38.p6.genome.fa chr1 > GRCm38.p6.genome.chr1.fa
samtools faidx GRCm38.p6.genome.chr1.fa


## Sequence context file
cd /scratcha/sblab/martin03/repository/20181213_DNMT_preference/data/20190624/reference2

ref=GRCm38.p6.genome.chr1.fa

fastaRegexFinder.py -f $ref -r C | \
cut -f1-3,6 | \
awk -v OFS="\t" '{print $1, $2, $3, ".", ".", $4}' | \
bedtools slop -i - -g $ref.fai -b 6 | \
bedtools getfasta -fi $ref -bed - -bedOut -s | \
awk 'length($7) == 13' | \
awk -v OFS="\t" '{print $1, $2+6, $3-6, $4, $5, $6, $7}' > GRCm38.p6.genome.chr1.context.bed

wc -l GRCm38.p6.genome.chr1.context.bed # 78962721 Cs with +/- 6bp sequence context in chr1 both in forward and reverse strands


## Methylation files
cd /scratcha/sblab/martin03/repository/20181213_DNMT_preference/data/20190624/methylation

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

wc -l *.context.txt
# 17100124 MouseES_BisSeq_HiSeq_SRR299053_1_bismark_bt2.multiple.deduplicated.coordsort.chr1.namesort.context.txt
#  3159295 TKO_DNMT3A2_bisulfite_hiseq2000_SRR8139275_1_bismark_bt2_pe.deduplicated.coordsort.chr1.namesort.context.txt
#  2047331 TKO_DNMT3B1_bisulfite_hiseq2000_SRR8139276_1_bismark_bt2_pe.deduplicated.coordsort.chr1.namesort.context.txt


## Concatenate context files
cd /scratchb/sblab/martin03/repository/20181213_DNMT_preference/data/20190624/methylation

tableCat.py -i *.context.txt | awk '{split($8,a,"_"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"a[1]"_"a[2]"_"a[3]"_"a[4]}' > TKO_DNMT3A2_DNMT3B1_MouseES.context.txt

exit
```
