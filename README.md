This repository contains data access, source code and demo of the custom computational analyses developed as part of our paper **Sequence preference of DNMT3A/B dictates genome-wide DNA methylation signature** *(under review)*


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

To access the publicly available human and mouse bisulfite and RNA-seq datasets - have a look at the Supplementary table


## Code

### Contents

- [System requirements and installation](README.md#system-requirements-and-installation)
- [Instructions for use and demo](README.md#instructions-for-use-and-demo)
  - [Sequence context diversity in reference genomes](README.md#sequence-context-diversity-in-reference-genomes)
  - [*E.coli* bisulfite sequencing](ecoli_bs.md)
  - [Mouse TKO DNMT bisulfite sequencing](mouse_tko_dnmt_bs.md)
  - [Human DNMT KO bisulfite sequencing](human_dnmt_ko_bs.md)
  - [Human tissues bisulfite sequencing](human_tissues_bs.md)
  - [Mouse tissues bisulfite sequencing](mouse_tissue_bs.md)
  - [Mouse brain development bisulfite sequencing](mouse_brain_development_bs.md)
  - [Human and mouse development bisulfite RNA-seq](human_mouse_development_rna-seq.md)



### System requirements and installation

Software, for installation details of the individual tools follow the links:

- [SRA toolkit.v2.8.0](https://github.com/ncbi/sra-tools)
- [EMBOSS v6.6.0.0](http://emboss.sourceforge.net/)
- [fastaRegexFinder.py v0.1.1](https://github.com/dariober/bioinformatics-cafe/tree/master/fastaRegexFinder)
- Standard Unix tools: cut, awk, sort, uniq, grep, tr, echo ...
- [bedtools v2.27.0](http://bedtools.readthedocs.io/en/latest/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [bismark v0.19.0](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
- [tableCat.py](https://github.com/dariober/bioinformatics-cafe/blob/master/tableCat/tableCat.py)
- [R v3.3.2](https://www.r-project.org/). Libraries:
  - [data.table v1.10.4](https://cran.r-project.org/web/packages/data.table/index.html)
  - [ggplot2 v2.2.1](http://ggplot2.org/)
  - [ggseqlogo v0.1](https://cran.rstudio.com/web/packages/ggseqlogo/index.html)
  - [gridExtra v2.2.1](https://cran.r-project.org/web/packages/gridExtra/index.html)
- [Python v2.7.12](https://www.python.org/)

Operating system:

  - CentOS Linux release 7.3.1611 (OS used during code development)
  - [slurm job scheduling system v19.05.0](https://slurm.schedmd.com/quickstart.html)


### Instructions for use and demo

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
