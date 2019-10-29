
[63 datasets available in ENCODE](https://www.encodeproject.org/matrix/?type=Experiment&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.classification=tissue&assay_slims=DNA+methylation&assay_title=WGBS&files.file_type=bed+bedMethyl)


## Contents

- [Download bedmethyl files](human_tissues_bs.md#download-bedmethyl-files)
- [Renaming](human_tissues_bs.md#renaming)
- [Extract chr20 info](human_tissues_bs.md#extract-chr20-info)
- [Extract context information around C](human_tissues_bs.md#extract-context-information-around-C)
- [Plotting](human_tissues_bs.md#plotting)



### Download bedmethyl files

```bash
cd ~

mkdir data && cd data

# downlaod files.txt from above ENCODE address
xargs -L 1 -P 12 curl -sSOL < files.txt &
```


### Renaming

```bash
mv ENCFF068IHF.bed.gz AdiposeTissue_ENCFF068IHF_CHH.bed.gz
mv ENCFF103DNU.bed.gz AdiposeTissue_ENCFF103DNU_CpG.bed.gz
mv ENCFF112SWY.bed.gz AdiposeTissue_ENCFF112SWY_CHH.bed.gz
mv ENCFF152EGF.bed.gz AdiposeTissue_ENCFF152EGF_CHG.bed.gz
mv ENCFF318AMC.bed.gz AdiposeTissue_ENCFF318AMC_CpG.bed.gz
mv ENCFF477GKI.bed.gz AdiposeTissue_ENCFF477GKI_CpG.bed.gz
mv ENCFF532BRS.bed.gz AdiposeTissue_ENCFF532BRS_CHG.bed.gz
mv ENCFF789XRU.bed.gz AdiposeTissue_ENCFF789XRU_CHG.bed.gz
mv ENCFF810FOG.bed.gz AdiposeTissue_ENCFF810FOG_CHH.bed.gz
mv ENCFF018AWT.bed.gz AdrenalGland_ENCFF018AWT_CHG.bed.gz
mv ENCFF106BER.bed.gz AdrenalGland_ENCFF106BER_CHH.bed.gz
mv ENCFF116KYV.bed.gz AdrenalGland_ENCFF116KYV_CHG.bed.gz
mv ENCFF210XTE.bed.gz AdrenalGland_ENCFF210XTE_CpG.bed.gz
mv ENCFF216DJL.bed.gz AdrenalGland_ENCFF216DJL_CpG.bed.gz
mv ENCFF282HNG.bed.gz AdrenalGland_ENCFF282HNG_CHG.bed.gz
mv ENCFF427PVN.bed.gz AdrenalGland_ENCFF427PVN_CHG.bed.gz
mv ENCFF483MGW.bed.gz AdrenalGland_ENCFF483MGW_CHH.bed.gz
mv ENCFF485BWP.bed.gz AdrenalGland_ENCFF485BWP_CHG.bed.gz
mv ENCFF552DWF.bed.gz AdrenalGland_ENCFF552DWF_CHH.bed.gz
mv ENCFF618WAT.bed.gz AdrenalGland_ENCFF618WAT_CpG.bed.gz
mv ENCFF631CXH.bed.gz AdrenalGland_ENCFF631CXH_CHH.bed.gz
mv ENCFF826PSA.bed.gz AdrenalGland_ENCFF826PSA_CpG.bed.gz
mv ENCFF913ZNZ.bed.gz AdrenalGland_ENCFF913ZNZ_CpG.bed.gz
mv ENCFF987LCQ.bed.gz AdrenalGland_ENCFF987LCQ_CHH.bed.gz
mv ENCFF017JID.bed.gz Aorta_ENCFF017JID_CHG.bed.gz
mv ENCFF223PYT.bed.gz Aorta_ENCFF223PYT_CHH.bed.gz
mv ENCFF553HJV.bed.gz Aorta_ENCFF553HJV_CpG.bed.gz
mv ENCFF245HSZ.bed.gz BodyOfPancreas_ENCFF245HSZ_CHH.bed.gz
mv ENCFF359JGR.bed.gz BodyOfPancreas_ENCFF359JGR_CHH.bed.gz
mv ENCFF699RBP.bed.gz BodyOfPancreas_ENCFF699RBP_CpG.bed.gz
mv ENCFF722JPQ.bed.gz BodyOfPancreas_ENCFF722JPQ_CHG.bed.gz
mv ENCFF748MTS.bed.gz BodyOfPancreas_ENCFF748MTS_CpG.bed.gz
mv ENCFF835YHZ.bed.gz BodyOfPancreas_ENCFF835YHZ_CHG.bed.gz
mv ENCFF250XVE.bed.gz Esophagus_ENCFF250XVE_CHG.bed.gz
mv ENCFF510EMT.bed.gz Esophagus_ENCFF510EMT_CpG.bed.gz
mv ENCFF546PSU.bed.gz Esophagus_ENCFF546PSU_CHG.bed.gz
mv ENCFF591PNU.bed.gz Esophagus_ENCFF591PNU_CHH.bed.gz
mv ENCFF625GVK.bed.gz Esophagus_ENCFF625GVK_CpG.bed.gz
mv ENCFF877IOX.bed.gz Esophagus_ENCFF877IOX_CHH.bed.gz
mv ENCFF560SMW.bed.gz Heart_ENCFF560SMW_CpG.bed.gz
mv ENCFF734QDJ.bed.gz Heart_ENCFF734QDJ_CHH.bed.gz
mv ENCFF994KKA.bed.gz Heart_ENCFF994KKA_CHG.bed.gz
mv ENCFF259ZEJ.bed.gz HeartLeftVentricle_ENCFF259ZEJ_CHH.bed.gz
mv ENCFF536RSX.bed.gz HeartLeftVentricle_ENCFF536RSX_CpG.bed.gz
mv ENCFF587RNQ.bed.gz HeartLeftVentricle_ENCFF587RNQ_CHH.bed.gz
mv ENCFF684JHX.bed.gz HeartLeftVentricle_ENCFF684JHX_CpG.bed.gz
mv ENCFF833HAZ.bed.gz HeartLeftVentricle_ENCFF833HAZ_CHG.bed.gz
mv ENCFF881IGR.bed.gz HeartLeftVentricle_ENCFF881IGR_CHG.bed.gz
mv ENCFF188WGW.bed.gz HeartRightVentricle_ENCFF188WGW_CHH.bed.gz
mv ENCFF375AGR.bed.gz HeartRightVentricle_ENCFF375AGR_CHG.bed.gz
mv ENCFF513ITC.bed.gz HeartRightVentricle_ENCFF513ITC_CpG.bed.gz
mv ENCFF566BGW.bed.gz HeartRightVentricle_ENCFF566BGW_CHG.bed.gz
mv ENCFF739YHS.bed.gz HeartRightVentricle_ENCFF739YHS_CHH.bed.gz
mv ENCFF831OYO.bed.gz HeartRightVentricle_ENCFF831OYO_CpG.bed.gz
mv ENCFF390EXZ.bed.gz LargeIntestine_ENCFF390EXZ_CHH.bed.gz
mv ENCFF411YBY.bed.gz LargeIntestine_ENCFF411YBY_CHG.bed.gz
mv ENCFF923CZC.bed.gz LargeIntestine_ENCFF923CZC_CpG.bed.gz
mv ENCFF071OFT.bed.gz LowerLegSkin_ENCFF071OFT_CHG.bed.gz
mv ENCFF121VIX.bed.gz LowerLegSkin_ENCFF121VIX_CpG.bed.gz
mv ENCFF177AGZ.bed.gz LowerLegSkin_ENCFF177AGZ_CHG.bed.gz
mv ENCFF219GCQ.bed.gz LowerLegSkin_ENCFF219GCQ_CpG.bed.gz
mv ENCFF523EPC.bed.gz LowerLegSkin_ENCFF523EPC_CHH.bed.gz
mv ENCFF673WUB.bed.gz LowerLegSkin_ENCFF673WUB_CHH.bed.gz
mv ENCFF039JFT.bed.gz Lung_ENCFF039JFT_CpG.bed.gz
mv ENCFF299AJA.bed.gz Lung_ENCFF299AJA_CHG.bed.gz
mv ENCFF477AUC.bed.gz Lung_ENCFF477AUC_CpG.bed.gz
mv ENCFF563DHU.bed.gz Lung_ENCFF563DHU_CHH.bed.gz
mv ENCFF721VBY.bed.gz Lung_ENCFF721VBY_CHH.bed.gz
mv ENCFF982BAD.bed.gz Lung_ENCFF982BAD_CHG.bed.gz
mv ENCFF030NLF.bed.gz MuscleOfLeg_ENCFF030NLF_CHH.bed.gz
mv ENCFF588ETU.bed.gz MuscleOfLeg_ENCFF588ETU_CpG.bed.gz
mv ENCFF717QGR.bed.gz MuscleOfLeg_ENCFF717QGR_CHG.bed.gz
mv ENCFF419NZF.bed.gz MuscleOfTrunk_ENCFF419NZF_CHH.bed.gz
mv ENCFF623YLP.bed.gz MuscleOfTrunk_ENCFF623YLP_CHG.bed.gz
mv ENCFF645AZF.bed.gz MuscleOfTrunk_ENCFF645AZF_CpG.bed.gz
mv ENCFF052DAU.bed.gz Ovary_ENCFF052DAU_CHG.bed.gz
mv ENCFF189WPY.bed.gz Ovary_ENCFF189WPY_CpG.bed.gz
mv ENCFF247ILV.bed.gz Ovary_ENCFF247ILV_CpG.bed.gz
mv ENCFF303ZGP.bed.gz Ovary_ENCFF303ZGP_CpG.bed.gz
mv ENCFF425UGD.bed.gz Ovary_ENCFF425UGD_CHH.bed.gz
mv ENCFF676ALA.bed.gz Ovary_ENCFF676ALA_CHH.bed.gz
mv ENCFF702EVG.bed.gz Ovary_ENCFF702EVG_CHG.bed.gz
mv ENCFF833ORQ.bed.gz Ovary_ENCFF833ORQ_CHG.bed.gz
mv ENCFF874VKO.bed.gz Ovary_ENCFF874VKO_CHH.bed.gz
mv ENCFF167RDF.bed.gz Pancreas_ENCFF167RDF_CHH.bed.gz
mv ENCFF340WWH.bed.gz Pancreas_ENCFF340WWH_CHH.bed.gz
mv ENCFF396EJH.bed.gz Pancreas_ENCFF396EJH_CHG.bed.gz
mv ENCFF500DKA.bed.gz Pancreas_ENCFF500DKA_CpG.bed.gz
mv ENCFF661ORW.bed.gz Pancreas_ENCFF661ORW_CHG.bed.gz
mv ENCFF763RUE.bed.gz Pancreas_ENCFF763RUE_CpG.bed.gz
mv ENCFF284SPU.bed.gz Placenta_ENCFF284SPU_CHH.bed.gz
mv ENCFF437OKM.bed.gz Placenta_ENCFF437OKM_CpG.bed.gz
mv ENCFF579LUK.bed.gz Placenta_ENCFF579LUK_CHG.bed.gz
mv ENCFF027KTR.bed.gz ProstateGland_ENCFF027KTR_CpG.bed.gz
mv ENCFF091WVA.bed.gz ProstateGland_ENCFF091WVA_CHG.bed.gz
mv ENCFF250DBS.bed.gz ProstateGland_ENCFF250DBS_CHH.bed.gz
mv ENCFF121ZES.bed.gz PsoasMuscle_ENCFF121ZES_CpG.bed.gz
mv ENCFF195FUI.bed.gz PsoasMuscle_ENCFF195FUI_CHG.bed.gz
mv ENCFF498IKD.bed.gz PsoasMuscle_ENCFF498IKD_CHG.bed.gz
mv ENCFF806ZTQ.bed.gz PsoasMuscle_ENCFF806ZTQ_CHH.bed.gz
mv ENCFF913UZU.bed.gz PsoasMuscle_ENCFF913UZU_CpG.bed.gz
mv ENCFF927XBL.bed.gz PsoasMuscle_ENCFF927XBL_CHH.bed.gz
mv ENCFF110AZO.bed.gz RightCardiacAtrium_ENCFF110AZO_CpG.bed.gz
mv ENCFF209DKJ.bed.gz RightCardiacAtrium_ENCFF209DKJ_CHH.bed.gz
mv ENCFF937SZG.bed.gz RightCardiacAtrium_ENCFF937SZG_CHG.bed.gz
mv ENCFF157POM.bed.gz SigmoidColon_ENCFF157POM_CpG.bed.gz
mv ENCFF336YNE.bed.gz SigmoidColon_ENCFF336YNE_CHG.bed.gz
mv ENCFF417ERA.bed.gz SigmoidColon_ENCFF417ERA_CHG.bed.gz
mv ENCFF455TQO.bed.gz SigmoidColon_ENCFF455TQO_CpG.bed.gz
mv ENCFF519PPJ.bed.gz SigmoidColon_ENCFF519PPJ_CHH.bed.gz
mv ENCFF576QRO.bed.gz SigmoidColon_ENCFF576QRO_CHH.bed.gz
mv ENCFF122LEF.bed.gz SmallIntestine_ENCFF122LEF_CpG.bed.gz
mv ENCFF152DFK.bed.gz SmallIntestine_ENCFF152DFK_CHH.bed.gz
mv ENCFF166NAS.bed.gz SmallIntestine_ENCFF166NAS_CHH.bed.gz
mv ENCFF204VWK.bed.gz SmallIntestine_ENCFF204VWK_CHG.bed.gz
mv ENCFF241AQC.bed.gz SmallIntestine_ENCFF241AQC_CpG.bed.gz
mv ENCFF266NGW.bed.gz SmallIntestine_ENCFF266NGW_CpG.bed.gz
mv ENCFF312ABE.bed.gz SmallIntestine_ENCFF312ABE_CHH.bed.gz
mv ENCFF500MJZ.bed.gz SmallIntestine_ENCFF500MJZ_CHG.bed.gz
mv ENCFF521DHD.bed.gz SmallIntestine_ENCFF521DHD_CpG.bed.gz
mv ENCFF684YXF.bed.gz SmallIntestine_ENCFF684YXF_CHG.bed.gz
mv ENCFF763TBC.bed.gz SmallIntestine_ENCFF763TBC_CHH.bed.gz
mv ENCFF971IXP.bed.gz SmallIntestine_ENCFF971IXP_CHG.bed.gz
mv ENCFF164EAU.bed.gz SpinalCord_ENCFF164EAU_CpG.bed.gz
mv ENCFF242KMS.bed.gz SpinalCord_ENCFF242KMS_CHG.bed.gz
mv ENCFF406RVP.bed.gz SpinalCord_ENCFF406RVP_CHH.bed.gz
mv ENCFF083GEQ.bed.gz Spleen_ENCFF083GEQ_CHH.bed.gz
mv ENCFF178PXB.bed.gz Spleen_ENCFF178PXB_CHG.bed.gz
mv ENCFF200MJQ.bed.gz Spleen_ENCFF200MJQ_CpG.bed.gz
mv ENCFF204WSG.bed.gz Spleen_ENCFF204WSG_CHG.bed.gz
mv ENCFF333OHK.bed.gz Spleen_ENCFF333OHK_CpG.bed.gz
mv ENCFF333VOO.bed.gz Spleen_ENCFF333VOO_CHH.bed.gz
mv ENCFF441NPK.bed.gz Spleen_ENCFF441NPK_CHG.bed.gz
mv ENCFF526PFA.bed.gz Spleen_ENCFF526PFA_CpG.bed.gz
mv ENCFF550FZT.bed.gz Spleen_ENCFF550FZT_CpG.bed.gz
mv ENCFF698NXR.bed.gz Spleen_ENCFF698NXR_CHH.bed.gz
mv ENCFF730NQT.bed.gz Spleen_ENCFF730NQT_CpG.bed.gz
mv ENCFF736IHW.bed.gz Spleen_ENCFF736IHW_CHG.bed.gz
mv ENCFF781KFP.bed.gz Spleen_ENCFF781KFP_CHG.bed.gz
mv ENCFF818IZZ.bed.gz Spleen_ENCFF818IZZ_CHH.bed.gz
mv ENCFF852FGI.bed.gz Spleen_ENCFF852FGI_CHH.bed.gz
mv ENCFF026NBR.bed.gz Stomach_ENCFF026NBR_CHG.bed.gz
mv ENCFF071MCV.bed.gz Stomach_ENCFF071MCV_CHG.bed.gz
mv ENCFF217SOG.bed.gz Stomach_ENCFF217SOG_CHG.bed.gz
mv ENCFF274AEJ.bed.gz Stomach_ENCFF274AEJ_CHH.bed.gz
mv ENCFF279AZR.bed.gz Stomach_ENCFF279AZR_CHH.bed.gz
mv ENCFF284GTX.bed.gz Stomach_ENCFF284GTX_CHH.bed.gz
mv ENCFF404YOB.bed.gz Stomach_ENCFF404YOB_CHH.bed.gz
mv ENCFF435SPL.bed.gz Stomach_ENCFF435SPL_CpG.bed.gz
mv ENCFF489CEV.bed.gz Stomach_ENCFF489CEV_CpG.bed.gz
mv ENCFF497YOO.bed.gz Stomach_ENCFF497YOO_CpG.bed.gz
mv ENCFF520WQW.bed.gz Stomach_ENCFF520WQW_CHG.bed.gz
mv ENCFF534RNT.bed.gz Stomach_ENCFF534RNT_CpG.bed.gz
mv ENCFF578VJG.bed.gz Stomach_ENCFF578VJG_CHG.bed.gz
mv ENCFF674JTR.bed.gz Stomach_ENCFF674JTR_CHG.bed.gz
mv ENCFF699TSO.bed.gz Stomach_ENCFF699TSO_CHH.bed.gz
mv ENCFF811QOG.bed.gz Stomach_ENCFF811QOG_CpG.bed.gz
mv ENCFF844EFX.bed.gz Stomach_ENCFF844EFX_CpG.bed.gz
mv ENCFF968GHE.bed.gz Stomach_ENCFF968GHE_CHH.bed.gz
mv ENCFF038JFQ.bed.gz Testis_ENCFF038JFQ_CHH.bed.gz
mv ENCFF497RCO.bed.gz Testis_ENCFF497RCO_CHG.bed.gz
mv ENCFF507JBR.bed.gz Testis_ENCFF507JBR_CHG.bed.gz
mv ENCFF638QVP.bed.gz Testis_ENCFF638QVP_CpG.bed.gz
mv ENCFF715DMX.bed.gz Testis_ENCFF715DMX_CpG.bed.gz
mv ENCFF949IHY.bed.gz Testis_ENCFF949IHY_CHH.bed.gz
mv ENCFF168HTX.bed.gz Thymus_ENCFF168HTX_CpG.bed.gz
mv ENCFF382ZHD.bed.gz Thymus_ENCFF382ZHD_CHH.bed.gz
mv ENCFF392XPZ.bed.gz Thymus_ENCFF392XPZ_CpG.bed.gz
mv ENCFF430MOL.bed.gz Thymus_ENCFF430MOL_CHG.bed.gz
mv ENCFF481RBU.bed.gz Thymus_ENCFF481RBU_CHG.bed.gz
mv ENCFF614LCK.bed.gz Thymus_ENCFF614LCK_CHH.bed.gz
mv ENCFF222RJM.bed.gz ThyroidGland_ENCFF222RJM_CHH.bed.gz
mv ENCFF223LJW.bed.gz ThyroidGland_ENCFF223LJW_CpG.bed.gz
mv ENCFF497IYX.bed.gz ThyroidGland_ENCFF497IYX_CpG.bed.gz
mv ENCFF633YWU.bed.gz ThyroidGland_ENCFF633YWU_CHH.bed.gz
mv ENCFF726QXC.bed.gz ThyroidGland_ENCFF726QXC_CHG.bed.gz
mv ENCFF867USJ.bed.gz ThyroidGland_ENCFF867USJ_CHG.bed.gz
mv ENCFF439JNT.bed.gz TibialNerve_ENCFF439JNT_CHH.bed.gz
mv ENCFF699KTW.bed.gz TibialNerve_ENCFF699KTW_CpG.bed.gz
mv ENCFF759UUY.bed.gz TibialNerve_ENCFF759UUY_CHG.bed.gz
mv ENCFF843SYR.bed.gz TibialNerve_ENCFF843SYR_CpG.bed.gz
mv ENCFF900SVP.bed.gz TibialNerve_ENCFF900SVP_CHG.bed.gz
mv ENCFF945DJY.bed.gz TibialNerve_ENCFF945DJY_CHH.bed.gz
mv ENCFF311XOU.bed.gz UpperLobeOfLeftLung_ENCFF311XOU_CHG.bed.gz
mv ENCFF724LUJ.bed.gz UpperLobeOfLeftLung_ENCFF724LUJ_CHH.bed.gz
mv ENCFF733EFJ.bed.gz UpperLobeOfLeftLung_ENCFF733EFJ_CpG.bed.gz
mv ENCFF842MHJ.bed.gz UpperLobeOfLeftLung_ENCFF842MHJ_CpG.bed.gz
mv ENCFF894SZK.bed.gz UpperLobeOfLeftLung_ENCFF894SZK_CHG.bed.gz
mv ENCFF936JAC.bed.gz UpperLobeOfLeftLung_ENCFF936JAC_CHH.bed.gz
```


### Extract chr20 info

```sh
cd ~/data

# add double quotes around variables to prevent issues caused by space in name
for bed in *.bed.gz
do
    bname=${bed%.bed.gz}
    echo "$bname"
    zcat "$bed" |
    awk -vOFS="\t" '{if ($1 == "chr20" && $10 > 9) print $1, $2, $3, $11, $10, $6}' > $bname.chr20.bed &
done
```


### Extract context information around C

```sh
cd ~/data

ref=genome.fa

for bed in *.bed
do
    bname=${bed%.bed}
    echo "$bname"
    awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $6}' "$bed" | \
    bedtools slop -i - -g "$ref".fai -b 3 | \
    bedtools getfasta -bedOut -s -fi "$ref" -bed - | \
    awk 'length($7) == 7' | \
    awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $6, toupper($7)}' > "$bname".context.txt &
done


# concatenate all CpG context files
tableCat.py -i *ENCFF*CpG*.context.txt -r .context.txt | awk '{split($8,a,"_"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"a[1]}' > human_tissue.CpG.chr20.context.txt &

# concatenate all CHG context files
tableCat.py -i *ENCFF*CHG*.context.txt -r .context.txt | awk '{split($8,a,"_"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"a[1]}' > human_tissue.CHG.chr20.context.txt &

# concatenate all CHH context files
tableCat.py -i *ENCFF*CHH*.context.txt -r .context.txt | awk '{split($8,a,"_"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"a[1]}' > human_tissue.CHH.chr20.context.txt &

# remove unneeded columns to save space and memory for later analysis
tableCat.py -i human_tissue.CH*.chr20.context.txt -r .context.txt | awk '{print $7"\t"$4"\t"$8}' > human_tissue.nonCpG.chr20.context.txt &

# only keep methylated cytosine for later analysis
cat human_tissue.{CpG,CHG,CHH}.chr20.context.txt | awk '{if ($4 > 0) print $7"\t"$4"\t"$8}' > human_tissue.methylated.chr20.context.txt &
```


### Plotting

```r
library(data.table)
library(ggplot2)

# Set width
options(width = 300)

# Load mfg data
data <- fread("~/data/human_tissue.nonCpG.chr20.context.txt")
setnames(data, c("context", "pct_met", "human_tissue"))

##########
## NCAN ##
##########

data_xcay <- copy(data[grepl("...CA..", data$context)])

# Extract XCAY context
data_xcay[, context_xcay := as.vector(sapply(data_xcay$context, function(x) paste(unlist(strsplit(x, ""))[3:6], collapse = "")))]

# remove context with N
# The like operator is a simple wrapper for grep(..., value=TRUE)
data_xcay <- data_xcay[!context_xcay %like% "N"]

# Explore xcay context
## Placenta
data_xcay[human_tissue == "Placenta", .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_xcay)][order(-pct_met_mean)]
#     context_xcay      N pct_met_median pct_met_mean
#  1:         ACAC 167424              0         1.13
#  2:         GCAC 210339              0         1.00
#  3:         CCAC 239483              0         0.87
#  4:         TCAC 220790              0         0.83
#  5:         ACAG 348028              0         0.46
#  6:         ACAT 195958              0         0.46
#  7:         CCAT 213981              0         0.44
#  8:         CCAG 455262              0         0.42
#  9:         GCAT 207934              0         0.41
# 10:         GCAG 401864              0         0.39
# 11:         TCAT 204828              0         0.38
# 12:         TCAG 340228              0         0.35
# 13:         ACAA 207536              0         0.35
# 14:         CCAA 226424              0         0.35
# 15:         GCAA 229779              0         0.31
# 16:         TCAA 213859              0         0.29

data_xcay[human_tissue == "Ovary", .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_xcay)][order(-pct_met_mean)]
#     context_xcay       N pct_met_median pct_met_mean
#  1:         CCAC 1303192              0         1.40
#  2:         ACAC 1267714              0         1.16
#  3:         GCAC  977705              0         1.12
#  4:         TCAC 1339092              0         0.97
#  5:         CCAT 1387012              0         0.93
#  6:         CCAG 1955426              0         0.91
#  7:         CCAA 1444538              0         0.84
#  8:         GCAG 1582975              0         0.64
#  9:         ACAG 1768887              0         0.63
# 10:         GCAT 1107819              0         0.60
# 11:         ACAT 1615422              0         0.58
# 12:         TCAG 1652051              0         0.57
# 13:         TCAT 1542378              0         0.57
# 14:         GCAA 1244459              0         0.51
# 15:         TCAA 1569985              0         0.51
# 16:         ACAA 1732149              0         0.50

data_xcay[human_tissue == "Testis", .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_xcay)][order(-pct_met_mean)]
#     context_xcay       N pct_met_median pct_met_mean
#  1:         CCAC  632218              0         1.90
#  2:         ACAC  690722              0         1.77
#  3:         GCAC  514525              0         1.66
#  4:         TCAC  698022              0         1.38
#  5:         CCAG  997462              0         1.18
#  6:         CCAT  720127              0         1.11
#  7:         CCAA  786157              0         1.01
#  8:         ACAG  974180              0         0.93
#  9:         GCAG  844753              0         0.83
# 10:         ACAT  926875              0         0.79
# 11:         TCAG  877999              0         0.75
# 12:         GCAT  606178              0         0.74
# 13:         TCAT  847043              0         0.71
# 14:         ACAA 1030360              0         0.65
# 15:         GCAA  705297              0         0.64
# 16:         TCAA  892987              0         0.62

data_xcay[human_tissue == "MuscleOfLeg", .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_xcay)][order(-pct_met_mean)]
#     context_xcay      N pct_met_median pct_met_mean
#  1:         ACAC 120223              0         0.75
#  2:         GCAC 165182              0         0.67
#  3:         CCAC 179895              0         0.64
#  4:         TCAC 158229              0         0.58
#  5:         ACAT 131497              0         0.50
#  6:         GCAT 154290              0         0.49
#  7:         CCAT 152309              0         0.49
#  8:         ACAG 266853              0         0.45
#  9:         TCAT 136828              0         0.45
# 10:         GCAG 333380              0         0.44
# 11:         CCAG 360694              0         0.42
# 12:         ACAA 142010              0         0.42
# 13:         CCAA 164978              0         0.42
# 14:         TCAG 256965              0         0.41
# 15:         GCAA 171466              0         0.40
# 16:         TCAA 144547              0         0.38

data_xcay[human_tissue == "TibialNerve", .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_xcay)][order(-pct_met_mean)]
#     context_xcay       N pct_met_median pct_met_mean
#  1:         CCAC  677255              0         1.40
#  2:         ACAC  738431              0         1.06
#  3:         GCAC  556274              0         1.02
#  4:         CCAT  777090              0         1.02
#  5:         CCAG 1081143              0         0.99
#  6:         TCAC  751810              0         0.93
#  7:         CCAA  846421              0         0.92
#  8:         GCAG  915704              0         0.70
#  9:         ACAG 1043585              0         0.65
# 10:         TCAG  947128              0         0.63
# 11:         GCAT  654677              0         0.62
# 12:         TCAT  912008              0         0.61
# 13:         ACAT  989413              0         0.59
# 14:         GCAA  753829              0         0.55
# 15:         TCAA  952370              0         0.55
# 16:         ACAA 1091771              0         0.52

data_xcay[human_tissue == "TibialNerve", .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_xcay)][order(-pct_met_mean)]
#     context_xcay       N pct_met_median pct_met_mean
#  1:         CCAC  677255              0         1.40
#  2:         ACAC  738431              0         1.06
#  3:         GCAC  556274              0         1.02
#  4:         CCAT  777090              0         1.02
#  5:         CCAG 1081143              0         0.99
#  6:         TCAC  751810              0         0.93
#  7:         CCAA  846421              0         0.92
#  8:         GCAG  915704              0         0.70
#  9:         ACAG 1043585              0         0.65
# 10:         TCAG  947128              0         0.63
# 11:         GCAT  654677              0         0.62
# 12:         TCAT  912008              0         0.61
# 13:         ACAT  989413              0         0.59
# 14:         GCAA  753829              0         0.55
# 15:         TCAA  952370              0         0.55
# 16:         ACAA 1091771              0         0.52

# save data_cay table
fwrite(data_xcay[, .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(human_tissue, context_xcay)][order(-pct_met_mean)], file = "~/data/human_tissue_chr20_NCAN.txt", sep = "\t", row.names=FALSE, quote=FALSE)

rm(data_xcay)


#########
## CAN ##
#########

data_cay <- copy(data[grepl("...CA..", data$context)])

# Extract cay context
data_cay[, context_cay := as.vector(sapply(data_cay$context, function(x) paste(unlist(strsplit(x, ""))[4:6], collapse = "")))]

# remove context with N
# The like operator is a simple wrapper for grep(..., value=TRUE)
data_cay <- data_cay[!context_cay %like% "N"]

# Explore cay context
## Placenta
data_cay[human_tissue == "Placenta", .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
#    context_cay       N pct_met_median pct_met_mean
# 1:         CAC  838036              0         0.94
# 2:         CAT  822701              0         0.42
# 3:         CAG 1545382              0         0.40
# 4:         CAA  877598              0         0.33

data_cay[human_tissue == "Ovary", .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
#    context_cay       N pct_met_median pct_met_mean
# 1:         CAC 4887703              0         1.16
# 2:         CAG 6959339              0         0.70
# 3:         CAT 5652631              0         0.67
# 4:         CAA 5991131              0         0.59

data_cay[human_tissue == "Testis", .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
#    context_cay       N pct_met_median pct_met_mean
# 1:         CAC 2535487              0         1.67
# 2:         CAG 3694394              0         0.93
# 3:         CAT 3100224              0         0.83
# 4:         CAA 3414801              0         0.72

# save data_cay table
fwrite(data_cay[, .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(human_tissue, context_cay)][order(-pct_met_mean)], file = "/scratcha/sblab/mao01/20190718_human_tissue_wgbs/20190719_human_tissue_chr20_CAN.txt", sep = "\t", row.names=FALSE, quote=FALSE)


############ CAC/CAG ##############
cac_cag <- dcast(data_cay[, .(pct_met = mean(pct_met, na.rm=TRUE)), by = .(human_tissue, context_cay)], human_tissue ~ context_cay, value.var = "pct_met")
cac_cag[, CAC_CAG := CAC/CAG]
cac_cag[order(CAC_CAG)]
#            human_tissue       CAA       CAC       CAG       CAT  CAC_CAG
#  1:              Thymus 0.5006638 0.7483201 0.5629262 0.5835693 1.329340
#  2:      LargeIntestine 1.0393552 1.6054828 1.1867700 1.1262631 1.352817
#  3: UpperLobeOfLeftLung 0.6175940 1.0536052 0.7358940 0.6867047 1.431735
#  4:              Spleen 0.5603811 0.9408853 0.6456746 0.6328499 1.457213
#  5:         TibialNerve 0.6273569 1.0981346 0.7466446 0.6993104 1.470759
#  6:        LowerLegSkin 0.6063790 1.0696134 0.7226043 0.6824536 1.480220
#  7:        ThyroidGland 0.6118996 1.1048872 0.7390538 0.6928680 1.495002
#  8:         MuscleOfLeg 0.4036703 0.6540498 0.4289231 0.4845684 1.524865
#  9:       ProstateGland 0.6011349 1.1068201 0.7258137 0.6862457 1.524937
# 10:               Heart 0.4991144 0.8181888 0.5302931 0.5839153 1.542899
# 11:        SigmoidColon 0.5383444 0.9668448 0.6204463 0.6183566 1.558305
# 12:      SmallIntestine 0.5636567 1.0123029 0.6484970 0.6314459 1.560999
# 13:                Lung 0.5193606 0.9784718 0.6126233 0.6010943 1.597183
# 14:             Stomach 0.5363041 0.9805265 0.6118826 0.6149256 1.602475
# 15:           Esophagus 0.5342156 1.0000482 0.6228066 0.6119235 1.605712
# 16:      BodyOfPancreas 0.6200443 1.1967990 0.7444648 0.6883698 1.607597
# 17:               Ovary 0.5853568 1.1616641 0.6982288 0.6681697 1.663730
# 18:       AdiposeTissue 0.5500266 1.1080573 0.6508966 0.6377314 1.702355
# 19:               Aorta 0.5395777 1.1220032 0.6433499 0.6239357 1.744001
# 20:        AdrenalGland 0.5982425 1.2703552 0.7112889 0.6909096 1.785991
# 21:              Testis 0.7214128 1.6736785 0.9336625 0.8327337 1.792595
# 22:       MuscleOfTrunk 0.3810520 0.7218410 0.4011979 0.4726336 1.799214
# 23:  RightCardiacAtrium 0.5417952 1.1838345 0.6519282 0.6388541 1.815897
# 24:            Pancreas 0.5507324 1.2232407 0.6673877 0.6447651 1.832879
# 25: HeartRightVentricle 0.6384832 1.7130347 0.8182735 0.7728564 2.093474
# 26:  HeartLeftVentricle 0.6125030 1.6595298 0.7917043 0.7404093 2.096148
# 27:         PsoasMuscle 0.6591718 2.1122773 0.9290245 0.8439796 2.273651
# 28:            Placenta 0.3266051 0.9442267 0.4046631 0.4218252 2.333365
# 29:          SpinalCord 0.4840304 1.5761166 0.5796134 0.6592147 2.719255
#            human_tissue       CAA       CAC       CAG       CAT  CAC_CAG

cac_cag$human_tissue <- factor(cac_cag$human_tissue, level = cac_cag$human_tissue[order(cac_cag$CAC_CAG)])

gg <- ggplot(cac_cag, aes(x=human_tissue, y=CAC_CAG)) +
geom_bar(stat="identity", color="black", position=position_dodge(), width= 0.6, alpha = 0.5) +
theme_bw() +
ylab(expression("[mCAC/CAC]/[mCAG/CAG]")) +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 12, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
coord_cartesian(ylim = c(0, 3))
ggsave("~/data/human_tissue_chr20_CAC_CAG.pdf")


############ CAY_CAR ##############
cac_cag <- dcast(data_cay[, .(pct_met = mean(pct_met, na.rm=TRUE)), by = .(human_tissue, context_cay)], human_tissue ~ context_cay, value.var = "pct_met")
cac_cag[, `:=`(CAC_CAG = CAC/CAG, CAY_CAR = (CAC + CAT)/(CAG + CAA))]
cac_cag[order(CAC_CAG)]
#           human_tissue       CAA       CAC       CAG       CAT  CAC_CAG  CAY_CAR
#  1:              Thymus 0.5006638 0.7483201 0.5629262 0.5835693 1.329340 1.252258
#  2:      LargeIntestine 1.0393552 1.6054828 1.1867700 1.1262631 1.352817 1.227130
#  3: UpperLobeOfLeftLung 0.6175940 1.0536052 0.7358940 0.6867047 1.431735 1.285796
#  4:              Spleen 0.5603811 0.9408853 0.6456746 0.6328499 1.457213 1.304861
#  5:         TibialNerve 0.6273569 1.0981346 0.7466446 0.6993104 1.470759 1.308183
#  6:        LowerLegSkin 0.6063790 1.0696134 0.7226043 0.6824536 1.480220 1.318351
#  7:        ThyroidGland 0.6118996 1.1048872 0.7390538 0.6928680 1.495002 1.330731
#  8:         MuscleOfLeg 0.4036703 0.6540498 0.4289231 0.4845684 1.524865 1.367556
#  9:       ProstateGland 0.6011349 1.1068201 0.7258137 0.6862457 1.524937 1.351270
# 10:               Heart 0.4991144 0.8181888 0.5302931 0.5839153 1.542899 1.362050
# 11:        SigmoidColon 0.5383444 0.9668448 0.6204463 0.6183566 1.558305 1.367979
# 12:      SmallIntestine 0.5636567 1.0123029 0.6484970 0.6314459 1.560999 1.356056
# 13:                Lung 0.5193606 0.9784718 0.6126233 0.6010943 1.597183 1.395396
# 14:             Stomach 0.5363041 0.9805265 0.6118826 0.6149256 1.602475 1.389541
# 15:           Esophagus 0.5342156 1.0000482 0.6228066 0.6119235 1.605712 1.393207
# 16:      BodyOfPancreas 0.6200443 1.1967990 0.7444648 0.6883698 1.607597 1.381573
# 17:               Ovary 0.5853568 1.1616641 0.6982288 0.6681697 1.663730 1.425564
# 18:       AdiposeTissue 0.5500266 1.1080573 0.6508966 0.6377314 1.702355 1.453705
# 19:               Aorta 0.5395777 1.1220032 0.6433499 0.6239357 1.744001 1.475947
# 20:        AdrenalGland 0.5982425 1.2703552 0.7112889 0.6909096 1.785991 1.497684
# 21:              Testis 0.7214128 1.6736785 0.9336625 0.8327337 1.792595 1.514380
# 22:       MuscleOfTrunk 0.3810520 0.7218410 0.4011979 0.4726336 1.799214 1.526973
# 23:  RightCardiacAtrium 0.5417952 1.1838345 0.6519282 0.6388541 1.815897 1.526894
# 24:            Pancreas 0.5507324 1.2232407 0.6673877 0.6447651 1.832879 1.533515
# 25: HeartRightVentricle 0.6384832 1.7130347 0.8182735 0.7728564 2.093474 1.706456
# 26:  HeartLeftVentricle 0.6125030 1.6595298 0.7917043 0.7404093 2.096148 1.709106
# 27:         PsoasMuscle 0.6591718 2.1122773 0.9290245 0.8439796 2.273651 1.861393
# 28:            Placenta 0.3266051 0.9442267 0.4046631 0.4218252 2.333365 1.868059
# 29:          SpinalCord 0.4840304 1.5761166 0.5796134 0.6592147 2.719255 2.101579
#            human_tissue       CAA       CAC       CAG       CAT  CAC_CAG  CAY_CAR


cac_cag$human_tissue <- factor(cac_cag$human_tissue, level = cac_cag$human_tissue[order(cac_cag$CAY_CAR)])

gg <- ggplot(cac_cag, aes(x=human_tissue, y=CAY_CAR)) +
geom_bar(stat="identity", color="black", position=position_dodge(), width= 0.6, alpha = 0.5) +
theme_bw() +
ylab(expression("[mCAC/CAC+mCAT/CAT]/[mCAG/CAG+mCAA/CAA]")) +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 12, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
coord_cartesian(ylim = c(0, 2.5))
ggsave("~/data/human_tissue_chr20_CAY_CAR.pdf")

rm(cac_cag)
rm(data_cay)
```
