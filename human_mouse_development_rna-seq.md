
## Contents

- [Download fastq files](human_mouse_development_rna-seq.md#download-fastq-files)
- [Rename fastq files](human_mouse_development_rna-seq.md#rename-fastq-files)
- [Quality check](human_mouse_development_rna-seq.md#quality-check)
- [Trim illumina adapters and quality trimming](human_mouse_development_rna-seq.md#trim-illumina-adapters-and-quality-trimming)
- [Merge fastq files](human_mouse_development_rna-seq.md#merge-fastq-files)
- [Alignment](human_mouse_development_rna-seq.md#alignment)
- [Obtain tables](human_mouse_development_rna-seq.md#obtain-tables)



### Download fastq files

Go to [GSE47966](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47966). Click on `SRP026048` and select `Source:RNA`. Click on `Send to:`, then `Run Selector`, then `Go`. Finally click on the `RunInfo Table` button, which will download the file `SraRunTable.txt` to the Desktop. Copy it to the cluster:

```bash
cd ~/geo

mkdir fastq && cd fastq

for run_id in `cut -f8 ../SraRunTable.txt | tail -n +2`
do
  sbatch -J $run_id -o $run_id.log --mem 4G --wrap "prefetch $run_id && \
  fastq-dump --outdir . --gzip --split-files ~/sra/$run_id.sra && \
  rm ~/sra/$run_id.sra"
done
```



### Rename fastq files

```bash
cd ~/geo/fastq

# SRR921934, SRR921935: RNA-Seq_hs_mfg_35do
for filename in SRR92193{4,5}_{1,2}.fastq.gz; do mv "$filename" "RNA-Seq_hs_mfg_35do_$filename"; done;

# SRR921936, SRR921937: RNA-Seq_hs_mfg_2yr
for filename in SRR92193{6,7}_{1,2}.fastq.gz; do mv "$filename" "RNA-Seq_hs_mfg_2yr_$filename"; done;

# SRR921938, SRR921939: RNA-Seq_hs_mfg_5yr
for filename in SRR92193{8,9}_{1,2}.fastq.gz; do mv "$filename" "RNA-Seq_hs_mfg_5yr_$filename"; done;

# SRR921940, SRR921941: RNA-Seq_hs_mfg_12yr
for filename in SRR92194{0,1}_{1,2}.fastq.gz; do mv "$filename" "RNA-Seq_hs_mfg_12yr_$filename"; done;

# SRR921942, SRR921943: RNA-Seq_hs_mfg_16yr
for filename in SRR92194{2,3}_{1,2}.fastq.gz; do mv "$filename" "RNA-Seq_hs_mfg_16yr_$filename"; done;

# SRR921944, SRR921945: RNA-Seq_hs_mfg_25yr
for filename in SRR92194{4,5}_{1,2}.fastq.gz; do mv "$filename" "RNA-Seq_hs_mfg_25yr_$filename"; done;

# SRR921946, SRR921947, SRR921948, SRR921949: RNA-Seq_hs_fc_53yr
for filename in SRR92194{6,7,8,9}_1.fastq.gz; do mv "$filename" "RNA-Seq_hs_fc_53yr_$filename"; done;

# SRR921950, SRR921951, SRR921952, SRR921953, SRR921954: RNA-Seq_mm_fc_fetal_rep1
for filename in SRR92195{0,1,2,3,4}_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_fetal_rep1_$filename"; done;

# SRR921955, SRR921956, SRR921957, SRR921958: RNA-Seq_mm_fc_fetal_rep2
for filename in SRR92195{5,6,7,8}_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_fetal_rep2_$filename"; done;

# SRR921959, SRR921960, SRR921961, SRR921962: RNA-Seq_mm_fc_1wk_rep1a
for filename in SRR9219{59,60,61,62}_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_1wk_rep1a_$filename"; done;

# SRR921963: RNA-Seq_mm_fc_1wk_rep1b
for filename in SRR921963_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_1wk_rep1b_$filename"; done;

# SRR921964: RNA-Seq_mm_fc_1wk_rep2
for filename in SRR921964_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_1wk_rep2_$filename"; done;

# SRR921965: RNA-Seq_mm_fc_1wk_rep3
for filename in SRR921965_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_1wk_rep3_$filename"; done;

# SRR921966, SRR921967, SRR921968, SRR921969, SRR921970: RNA-Seq_mm_fc_2wk_rep1a
for filename in SRR9219{66,67,68,69,70}_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_2wk_rep1a_$filename"; done;

# SRR921971: RNA-Seq_mm_fc_2wk_rep1b
for filename in SRR921971_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_2wk_rep1b_$filename"; done;

# SRR921972: RNA-Seq_mm_fc_2wk_rep2
for filename in SRR921972_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_2wk_rep2_$filename"; done;

# SRR921973: RNA-Seq_mm_fc_2wk_rep3
for filename in SRR921973_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_2wk_rep3_$filename"; done;

# SRR921974, SRR921975, SRR921976, SRR921977, SRR921978: RNA-Seq_mm_fc_4wk_rep1
for filename in SRR92197{4,5,6,7,8}_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_4wk_rep1_$filename"; done;

# SRR921979, SRR921980, SRR921981, SRR921982: RNA-Seq_mm_fc_4wk_rep2
for filename in SRR9219{79,80,81,82}_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_4wk_rep2_$filename"; done;

# SRR921983, SRR921984, SRR921985, SRR921986, SRR921987: RNA-Seq_mm_fc_6wk_rep1a
for filename in SRR92198{3,4,5,6,7}_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_6wk_rep1a_$filename"; done;

# SRR921988: RNA-Seq_mm_fc_6wk_rep1b
for filename in SRR921988_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_6wk_rep1b_$filename"; done;

# SRR921989: RNA-Seq_mm_fc_6wk_rep2
for filename in SRR921989_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_6wk_rep2_$filename"; done;

# SRR921990: RNA-Seq_mm_fc_6wk_rep3
for filename in SRR921990_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_6wk_rep3_$filename"; done;

# SRR921991, SRR921992, SRR921993, SRR921994, SRR921995: RNA-Seq_mm_fc_10wk_rep1a
for filename in SRR92199{1,2,3,4,5}_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_10wk_rep1a_$filename"; done;

# SRR921996: RNA-Seq_mm_fc_10wk_rep1b
for filename in SRR921996_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_10wk_rep1b_$filename"; done;

# SRR921997: RNA-Seq_mm_fc_10wk_rep2
for filename in SRR921997_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_10wk_rep2_$filename"; done;

# SRR921998: RNA-Seq_mm_fc_10wk_rep3
for filename in SRR921998_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_10wk_rep3_$filename"; done;

# SRR921999, SRR922000, SRR922001, SRR922002: RNA-Seq_mm_fc_22mo
for filename in SRR92{1999,2000,2001,2002}_1.fastq.gz; do mv "$filename" "RNA-Seq_mm_fc_22mo_$filename"; done;
```



### Quality check

```bash
cd ~/geo/fastq

mkdir ../fastqc

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../fastqc/$bname.log --mem 4G --wrap "fastqc --noextract --nogroup -q -o ../fastqc $fq"
done
```



### Trim illumina adapters and quality trimming

```bash
cd ~/geo/fastq

mkdir ../fastq_trimmed

# human paired-end
for fq1 in RNA-Seq_hs_mfg_*_1.fastq.gz
do
  fq2=${fq1/_1./_2.}
  bname=${fq1%_1.fastq.gz}
  #echo $fq1, $fq2, $bname
  sbatch -J $bname -o ../fastq_trimmed/$bname.log --mem 4G --wrap "cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 15 -q 20 -o ../fastq_trimmed/$fq1 -p ../fastq_trimmed/$fq2 $fq1 $fq2 > ../fastq_trimmed/$bname.txt"
done

# human single-end
for fq in RNA-Seq_hs_fc_*_1.fastq.gz
do
  bname=${fq%_1.fastq.gz}
  #echo $fq1, $fq2, $bname
  sbatch -J $bname -o ../fastq_trimmed/$bname.log --mem 4G --wrap "cutadapt -a AGATCGGAAGAGC -m 15 -q 20 -o ../fastq_trimmed/$fq $fq > ../fastq_trimmed/$bname.txt"
done

# mouse single-end
for fq in RNA-Seq_mm_fc_*_1.fastq.gz
do
  bname=${fq%_1.fastq.gz}
  #echo $fq1, $fq2, $bname
  sbatch -J $bname -o ../fastq_trimmed/$bname.log --mem 4G --wrap "cutadapt -a AGATCGGAAGAGC -m 15 -q 20 -o ../fastq_trimmed/$fq $fq > ../fastq_trimmed/$bname.txt"
done
```



### Merge fastq files

```bash
cd ~/geo/fastq_trimmed

mkdir ../fastq_trimmed_merged

# SRR921934, SRR921935: RNA-Seq_hs_mfg_35do
nohup cat RNA-Seq_hs_mfg_35do_SRR921934_1.fastq.gz RNA-Seq_hs_mfg_35do_SRR921935_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_hs_mfg_35do_1.fastq.gz &
nohup cat RNA-Seq_hs_mfg_35do_SRR921934_2.fastq.gz RNA-Seq_hs_mfg_35do_SRR921935_2.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_hs_mfg_35do_2.fastq.gz &

# SRR921936, SRR921937: RNA-Seq_hs_mfg_2yr
nohup cat RNA-Seq_hs_mfg_2yr_SRR921936_1.fastq.gz RNA-Seq_hs_mfg_2yr_SRR921937_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_hs_mfg_2yr_1.fastq.gz &
nohup cat RNA-Seq_hs_mfg_2yr_SRR921936_2.fastq.gz RNA-Seq_hs_mfg_2yr_SRR921937_2.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_hs_mfg_2yr_2.fastq.gz &

# SRR921938, SRR921939: RNA-Seq_hs_mfg_5yr
nohup cat RNA-Seq_hs_mfg_5yr_SRR921938_1.fastq.gz RNA-Seq_hs_mfg_5yr_SRR921939_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_hs_mfg_5yr_1.fastq.gz &
nohup cat RNA-Seq_hs_mfg_5yr_SRR921938_2.fastq.gz RNA-Seq_hs_mfg_5yr_SRR921939_2.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_hs_mfg_5yr_2.fastq.gz &

# SRR921940, SRR921941: RNA-Seq_hs_mfg_12yr
nohup cat RNA-Seq_hs_mfg_12yr_SRR921940_1.fastq.gz RNA-Seq_hs_mfg_12yr_SRR921941_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_hs_mfg_12yr_1.fastq.gz &
nohup cat RNA-Seq_hs_mfg_12yr_SRR921940_2.fastq.gz RNA-Seq_hs_mfg_12yr_SRR921941_2.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_hs_mfg_12yr_2.fastq.gz &

# SRR921942, SRR921943: RNA-Seq_hs_mfg_16yr
nohup cat RNA-Seq_hs_mfg_16yr_SRR921942_1.fastq.gz RNA-Seq_hs_mfg_16yr_SRR921943_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_hs_mfg_16yr_1.fastq.gz &
nohup cat RNA-Seq_hs_mfg_16yr_SRR921942_2.fastq.gz RNA-Seq_hs_mfg_16yr_SRR921943_2.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_hs_mfg_16yr_2.fastq.gz &

# SRR921944, SRR921945: RNA-Seq_hs_mfg_25yr
nohup cat RNA-Seq_hs_mfg_25yr_SRR921944_1.fastq.gz RNA-Seq_hs_mfg_25yr_SRR921945_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_hs_mfg_25yr_1.fastq.gz &
nohup cat RNA-Seq_hs_mfg_25yr_SRR921944_2.fastq.gz RNA-Seq_hs_mfg_25yr_SRR921945_2.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_hs_mfg_25yr_2.fastq.gz &

# SRR921946, SRR921947, SRR921948, SRR921949: RNA-Seq_hs_fc_53yr
nohup cat RNA-Seq_hs_fc_53yr_SRR92194{6,7,8,9}_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_hs_fc_53yr_1.fastq.gz &

# SRR921950, SRR921951, SRR921952, SRR921953, SRR921954: RNA-Seq_mm_fc_fetal_rep1
nohup cat RNA-Seq_mm_fc_fetal_rep1_SRR92195{0,1,2,3,4}_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_mm_fc_fetal_rep1_1.fastq.gz &

# SRR921955, SRR921956, SRR921957, SRR921958: RNA-Seq_mm_fc_fetal_rep2
nohup cat RNA-Seq_mm_fc_fetal_rep2_SRR92195{5,6,7,8}_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_mm_fc_fetal_rep2_1.fastq.gz &

# SRR921959, SRR921960, SRR921961, SRR921962: RNA-Seq_mm_fc_1wk_rep1a
nohup cat RNA-Seq_mm_fc_1wk_rep1a_SRR9219{59,60,61,62}_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_mm_fc_1wk_rep1a_1.fastq.gz &

# SRR921963: RNA-Seq_mm_fc_1wk_rep1b
cp RNA-Seq_mm_fc_1wk_rep1b_SRR921963_1.fastq.gz ../fastq_trimmed_merged/RNA-Seq_mm_fc_1wk_rep1b_1.fastq.gz

# SRR921964: RNA-Seq_mm_fc_1wk_rep2
cp RNA-Seq_mm_fc_1wk_rep2_SRR921964_1.fastq.gz ../fastq_trimmed_merged/RNA-Seq_mm_fc_1wk_rep2_1.fastq.gz

# SRR921965: RNA-Seq_mm_fc_1wk_rep3
cp RNA-Seq_mm_fc_1wk_rep3_SRR921965_1.fastq.gz ../fastq_trimmed_merged/RNA-Seq_mm_fc_1wk_rep3_1.fastq.gz

# SRR921966, SRR921967, SRR921968, SRR921969, SRR921970: RNA-Seq_mm_fc_2wk_rep1a
nohup cat RNA-Seq_mm_fc_2wk_rep1a_SRR9219{66,67,68,69,70}_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_mm_fc_2wk_rep1a_1.fastq.gz &

# SRR921971: RNA-Seq_mm_fc_2wk_rep1b
cp RNA-Seq_mm_fc_2wk_rep1b_SRR921971_1.fastq.gz ../fastq_trimmed_merged/RNA-Seq_mm_fc_2wk_rep1b_1.fastq.gz

# SRR921972: RNA-Seq_mm_fc_2wk_rep2
cp RNA-Seq_mm_fc_2wk_rep2_SRR921972_1.fastq.gz ../fastq_trimmed_merged/RNA-Seq_mm_fc_2wk_rep2_1.fastq.gz

# SRR921973: RNA-Seq_mm_fc_2wk_rep3
cp RNA-Seq_mm_fc_2wk_rep3_SRR921973_1.fastq.gz ../fastq_trimmed_merged/RNA-Seq_mm_fc_2wk_rep3_1.fastq.gz

# SRR921974, SRR921975, SRR921976, SRR921977, SRR921978: RNA-Seq_mm_fc_4wk_rep1
nohup cat RNA-Seq_mm_fc_4wk_rep1_SRR92197{4,5,6,7,8}_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_mm_fc_4wk_rep1_1.fastq.gz &

# SRR921979, SRR921980, SRR921981, SRR921982: RNA-Seq_mm_fc_4wk_rep2
nohup cat RNA-Seq_mm_fc_4wk_rep2_SRR9219{79,80,81,82}_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_mm_fc_4wk_rep2_1.fastq.gz &

# SRR921983, SRR921984, SRR921985, SRR921986, SRR921987: RNA-Seq_mm_fc_6wk_rep1a
nohup cat RNA-Seq_mm_fc_6wk_rep1a_SRR92198{3,4,5,6,7}_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_mm_fc_6wk_rep1a_1.fastq.gz &

# SRR921988: RNA-Seq_mm_fc_6wk_rep1b
cp RNA-Seq_mm_fc_6wk_rep1b_SRR921988_1.fastq.gz ../fastq_trimmed_merged/RNA-Seq_mm_fc_6wk_rep1b_1.fastq.gz

# SRR921989: RNA-Seq_mm_fc_6wk_rep2
cp RNA-Seq_mm_fc_6wk_rep2_SRR921989_1.fastq.gz ../fastq_trimmed_merged/RNA-Seq_mm_fc_6wk_rep2_1.fastq.gz

# SRR921990: RNA-Seq_mm_fc_6wk_rep3
cp RNA-Seq_mm_fc_6wk_rep3_SRR921990_1.fastq.gz ../fastq_trimmed_merged/RNA-Seq_mm_fc_6wk_rep3_1.fastq.gz

# SRR921991, SRR921992, SRR921993, SRR921994, SRR921995: RNA-Seq_mm_fc_10wk_rep1a
nohup cat RNA-Seq_mm_fc_10wk_rep1a_SRR92199{1,2,3,4,5}_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_mm_fc_10wk_rep1a_1.fastq.gz &

# SRR921996: RNA-Seq_mm_fc_10wk_rep1b
cp RNA-Seq_mm_fc_10wk_rep1b_SRR921996_1.fastq.gz ../fastq_trimmed_merged/RNA-Seq_mm_fc_10wk_rep1b_1.fastq.gz

# SRR921997: RNA-Seq_mm_fc_10wk_rep2
cp RNA-Seq_mm_fc_10wk_rep2_SRR921997_1.fastq.gz ../fastq_trimmed_merged/RNA-Seq_mm_fc_10wk_rep2_1.fastq.gz

# SRR921998: RNA-Seq_mm_fc_10wk_rep3
cp RNA-Seq_mm_fc_10wk_rep3_SRR921998_1.fastq.gz ../fastq_trimmed_merged/RNA-Seq_mm_fc_10wk_rep3_1.fastq.gz

# SRR921999, SRR922000, SRR922001, SRR922002: RNA-Seq_mm_fc_22mo
nohup cat RNA-Seq_mm_fc_22mo_SRR92{1999,2000,2001,2002}_1.fastq.gz > ../fastq_trimmed_merged/RNA-Seq_mm_fc_22mo_1.fastq.gz &
```



### Alignment

#### Prepare rsem reference sequences

```bash
mkdir rsem_bowtie2_basic

# bowtie2
rsem-prepare-reference --gtf ../gtf/gencode.vM18.annotation.gtf \
--bowtie2 \
-p 20 \
../fasta/GRCm38.p6.genome.fa \
GRCm38.p6.genome
```


#### Align

```bash
cd ~/fastq_trimmed_merged

mkdir ../rsem

bwtidx_human=../rsem_bowtie2_basic/GRCh38.p12.genome
bwtidx_mouse=../rsem_bowtie2_basic/GRCm38.p6.genome


# human paired-end
for fq1 in RNA-Seq_hs_mfg_*_1.fastq.gz
do
  fq2=${fq1/_1./_2.}
  bname=${fq1%_1.fastq.gz}
  #echo $fq1, $fq2, $bname
  sbatch -J $bname -o ../rsem/$bname.log --mem 32G --wrap "rsem-calculate-expression \
  --paired-end \
  -p 20 \
  --bowtie2 \
  --append-names \
  $fq1 \
  $fq2 \
  $bwtidx_human \
  ../rsem/$bname"
done


# human single-end
for fq in RNA-Seq_hs_fc_*_1.fastq.gz
do
  bname=${fq%_1.fastq.gz}
  #echo $fq1, $fq2, $bname
  sbatch -J $bname -o ../rsem/$bname.log --mem 32G --wrap "rsem-calculate-expression \
  -p 20 \
  --bowtie2 \
  --append-names \
  $fq \
  $bwtidx_human \
  ../rsem/$bname"
done


# mouse single-end
for fq in RNA-Seq_mm_fc_*_1.fastq.gz
do
  bname=${fq%_1.fastq.gz}
  #echo $fq1, $fq2, $bname
  sbatch -J $bname -o ../rsem/$bname.log --mem 32G --wrap "rsem-calculate-expression \
  -p 20 \
  --bowtie2 \
  --append-names \
  $fq \
  $bwtidx_mouse \
  ../rsem/$bname"
done
```



### Obtain tables

- human genes
- human transcripts
- mouse genes
- mouse transcripts

```python
# cd ~/geo/rsem
import os

#########
# human #
#########
hs_genes = []

for f in os.listdir("."):
  if f.startswith("RNA-Seq_hs") & f.endswith(".genes.results"):
    bname=f.replace(".genes.results", "")
    ifile=open(f, "r")
    ilines=ifile.readlines()
    ifile.close()
    for line in ilines:
      fields = line.split()
      gene_id = fields[0]
      ensembl_id = gene_id.split("_")[0].split(".")[0]
      gene_name = gene_id.split("_")[1]
      length = fields[2]
      tpm = fields[5]
      fpkm = fields[6]
      if ("DNMT1" in gene_id) or ("DNMT3A" in gene_id) or ("DNMT3B" in gene_id) or ("DNMT3L" in gene_id) or ("TET1" in gene_id) or ("TET2" in gene_id) or ("TET3" in gene_id):
        hs_genes.append([bname, gene_name, ensembl_id, length, tpm, fpkm])


ofile = open("~/geo/human_genes.txt", "w")
ofile.write("\t".join(["library_id", "gene_name", "ensembl_id", "length_average", "TPM", "FPKM"]) + "\n")
ofile.write("\n".join(["\t".join(i) for i in hs_genes]) + "\n")
ofile.close()


hs_transcripts = []

for f in os.listdir("."):
  if f.startswith("RNA-Seq_hs") & f.endswith(".isoforms.results"):
    bname=f.replace(".isoforms.results", "")
    ifile=open(f, "r")
    ilines=ifile.readlines()
    ifile.close()
    for line in ilines:
      fields = line.split()
      transcript_id = fields[0]
      transcript_ensembl_id = transcript_id.split("_")[0].split(".")[0]
      transcript_name = transcript_id.split("_")[1]
      gene_id = fields[1]
      gene_ensembl_id = gene_id.split("_")[0].split(".")[0]
      gene_name = gene_id.split("_")[1]
      length = fields[2]
      tpm = fields[5]
      fpkm = fields[6]
      isopct = fields[7]
      if ("DNMT1" in gene_id) or ("DNMT3A" in gene_id) or ("DNMT3B" in gene_id) or ("DNMT3L" in gene_id) or ("TET1" in gene_id) or ("TET2" in gene_id) or ("TET3" in gene_id):
        hs_transcripts.append([bname, transcript_name, transcript_ensembl_id, gene_name, gene_ensembl_id, length, tpm, fpkm, isopct])


ofile = open("~/geo/human_transcripts.txt", "w")
ofile.write("\t".join(["library_id", "transcript_name", "transcript_ensembl_id", "gene_name", "gene_ensembl_id", "length", "TPM", "FPKM", "IsoPct"]) + "\n")
ofile.write("\n".join(["\t".join(i) for i in hs_transcripts]) + "\n")
ofile.close()


#########
# mouse #
#########
mm_genes = []

for f in os.listdir("."):
  if f.startswith("RNA-Seq_mm") & f.endswith(".genes.results"):
    bname=f.replace(".genes.results", "")
    ifile=open(f, "r")
    ilines=ifile.readlines()
    ifile.close()
    for line in ilines:
      fields = line.split()
      gene_id = fields[0]
      ensembl_id = gene_id.split("_")[0].split(".")[0]
      gene_name = gene_id.split("_")[1]
      length = fields[2]
      tpm = fields[5]
      fpkm = fields[6]
      if (gene_name == "Dnmt1") or (gene_name == "Dnmt3a") or (gene_name == "Dnmt3b") or (gene_name == "Dnmt3l") or (gene_name == "Tet1") or (gene_name == "Tet2") or (gene_name == "Tet3"):
        mm_genes.append([bname, gene_name, ensembl_id, length, tpm, fpkm])


ofile = open("~/geo/mouse_genes.txt", "w")
ofile.write("\t".join(["library_id", "gene_name", "ensembl_id", "length_average", "TPM", "FPKM"]) + "\n")
ofile.write("\n".join(["\t".join(i) for i in mm_genes]) + "\n")
ofile.close()


mm_transcripts = []

for f in os.listdir("."):
  if f.startswith("RNA-Seq_mm") & f.endswith(".isoforms.results"):
    bname=f.replace(".isoforms.results", "")
    ifile=open(f, "r")
    ilines=ifile.readlines()
    ifile.close()
    for line in ilines:
      fields = line.split()
      transcript_id = fields[0]
      transcript_ensembl_id = transcript_id.split("_")[0].split(".")[0]
      transcript_name = transcript_id.split("_")[1]
      gene_id = fields[1]
      gene_ensembl_id = gene_id.split("_")[0].split(".")[0]
      gene_name = gene_id.split("_")[1]
      length = fields[2]
      tpm = fields[5]
      fpkm = fields[6]
      isopct = fields[7]
      if (gene_name == "Dnmt1") or (gene_name == "Dnmt3a") or (gene_name == "Dnmt3b") or (gene_name == "Dnmt3l") or (gene_name == "Tet1") or (gene_name == "Tet2") or (gene_name == "Tet3"):
        mm_transcripts.append([bname, transcript_name, transcript_ensembl_id, gene_name, gene_ensembl_id, length, tpm, fpkm, isopct])


ofile = open("~/geo/mouse_transcripts.txt", "w")
ofile.write("\t".join(["library_id", "transcript_name", "transcript_ensembl_id", "gene_name", "gene_ensembl_id", "length", "TPM", "FPKM", "IsoPct"]) + "\n")
ofile.write("\n".join(["\t".join(i) for i in mm_transcripts]) + "\n")
ofile.close()
```

Same steps as above were performed for RNA-seq files from ENCODE:

- ENCFF920CNZ,ENCFF320FJX: RNA-Seq_mm_forebrain_10.5days_rep1
- ENCFF663SNC,ENCFF528EVC: RNA-Seq_mm_forebrain_10.5days_rep2
- ENCFF294JRP,ENCFF920QAY: RNA-Seq_mm_forebrain_12.5days_rep1
- ENCFF700OLU,ENCFF203BWA: RNA-Seq_mm_forebrain_12.5days_rep2
- ENCFF931IVO: RNA-Seq_mm_forebrain_16.5days_rep1
- ENCFF114DRT: RNA-Seq_mm_forebrain_16.5days_rep2
- ENCFF179JEC: RNA-Seq_mm_forebrain_15.5days_rep1
- ENCFF891HIX: RNA-Seq_mm_forebrain_15.5days_rep2
- ENCFF235DNM: RNA-Seq_mm_forebrain_13.5days_rep1
- ENCFF959PSX: RNA-Seq_mm_forebrain_13.5days_rep2
- ENCFF037JQC,ENCFF358MFI: RNA-Seq_mm_forebrain_0days_rep1
- ENCFF447EXU,ENCFF458NWF: RNA-Seq_mm_forebrain_0days_rep2
- ENCFF329ACL: RNA-Seq_mm_forebrain_11.5days_rep1
- ENCFF896COV,ENCFF251LNG: RNA-Seq_mm_forebrain_11.5days_rep2
- ENCFF270GKY,ENCFF460TCF: RNA-Seq_mm_forebrain_14.5days_rep1
- ENCFF126IRS,ENCFF748SRJ: RNA-Seq_mm_forebrain_14.5days_rep2
