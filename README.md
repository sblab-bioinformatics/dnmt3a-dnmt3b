This repository contains data access and computational analyses for the methods developed in our DNMT3A/B sequence preference paper (*under review*)


## Data

All the raw data from the E.coli bisulfite sequencing experiments have been deposited in the ArrayExpress database at EMBL-EBI under accession number [E-MTAB-8406](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8406) (*currently only accesible to referees*)

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



## Code

### Contents

- [Software required](README.md#software-required)
- [Sequence context diversity in reference genomes](README.md#sequence-context-diversity-in-reference-genomes)



### Software required

- [EMBOSS v6.6.0.0](http://emboss.sourceforge.net/)
- [fastaRegexFinder.py v0.1.1](https://github.com/dariober/bioinformatics-cafe/tree/master/fastaRegexFinder)
- Standard Unix tools: cut, awk, sort, uniq, grep, tr, echo ...
- [bedtools v2.27.0](http://bedtools.readthedocs.io/en/latest/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)



### Sequence context diversity in reference genomes

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

