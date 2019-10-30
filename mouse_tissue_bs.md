[72 datasets available from ENCODE]
<https://www.encodeproject.org/matrix/?type=Experiment&status=released&assay_slims=DNA+methylation&assay_title=WGBS&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&biosample_ontology.classification=tissue&files.file_type=bed+bedMethyl>

## Download bedmethyl files

```sh
cd ~
mkdir -p ~/data
cd ~/data
# downlaod files.txt from above ENCODE address
xargs -L 1 -P 12 curl -sSOL < files.txt &
```

## rename bedmethyl files

```sh
cd ~/data

mv ENCFF006QDP.bed.gz EmbryonicFacialProminence_E10.5_ENCFF006QDP_1_CHH.bed.gz
mv ENCFF050FNM.bed.gz EmbryonicFacialProminence_E10.5_ENCFF050FNM_2_CHH.bed.gz
mv ENCFF095VLL.bed.gz EmbryonicFacialProminence_E10.5_ENCFF095VLL_2_CpG.bed.gz
mv ENCFF198NON.bed.gz EmbryonicFacialProminence_E10.5_ENCFF198NON_2_CHG.bed.gz
mv ENCFF722PJA.bed.gz EmbryonicFacialProminence_E10.5_ENCFF722PJA_1_CHG.bed.gz
mv ENCFF886MNA.bed.gz EmbryonicFacialProminence_E10.5_ENCFF886MNA_1_CpG.bed.gz
mv ENCFF083XWZ.bed.gz EmbryonicFacialProminence_E11.5_ENCFF083XWZ_2_CpG.bed.gz
mv ENCFF306DKT.bed.gz EmbryonicFacialProminence_E11.5_ENCFF306DKT_1_CpG.bed.gz
mv ENCFF590MRS.bed.gz EmbryonicFacialProminence_E11.5_ENCFF590MRS_1_CHG.bed.gz
mv ENCFF754BON.bed.gz EmbryonicFacialProminence_E11.5_ENCFF754BON_2_CHH.bed.gz
mv ENCFF884JVJ.bed.gz EmbryonicFacialProminence_E11.5_ENCFF884JVJ_1_CHH.bed.gz
mv ENCFF997OLQ.bed.gz EmbryonicFacialProminence_E11.5_ENCFF997OLQ_2_CHG.bed.gz
mv ENCFF063YAO.bed.gz EmbryonicFacialProminence_E14.5_ENCFF063YAO_1_CHH.bed.gz
mv ENCFF686ZQQ.bed.gz EmbryonicFacialProminence_E14.5_ENCFF686ZQQ_2_CpG.bed.gz
mv ENCFF709AXF.bed.gz EmbryonicFacialProminence_E14.5_ENCFF709AXF_2_CHH.bed.gz
mv ENCFF744VPT.bed.gz EmbryonicFacialProminence_E14.5_ENCFF744VPT_1_CHG.bed.gz
mv ENCFF757MSK.bed.gz EmbryonicFacialProminence_E14.5_ENCFF757MSK_2_CHG.bed.gz
mv ENCFF989UNH.bed.gz EmbryonicFacialProminence_E14.5_ENCFF989UNH_1_CpG.bed.gz
mv ENCFF365XZL.bed.gz Forebrain_E10.5_ENCFF365XZL_2_CpG.bed.gz
mv ENCFF369TZO.bed.gz Forebrain_E10.5_ENCFF369TZO_1_CHG.bed.gz
mv ENCFF400MQF.bed.gz Forebrain_E10.5_ENCFF400MQF_2_CHH.bed.gz
mv ENCFF506SUF.bed.gz Forebrain_E10.5_ENCFF506SUF_2_CHG.bed.gz
mv ENCFF590OFM.bed.gz Forebrain_E10.5_ENCFF590OFM_1_CHH.bed.gz
mv ENCFF977BKF.bed.gz Forebrain_E10.5_ENCFF977BKF_1_CpG.bed.gz
mv ENCFF144FZS.bed.gz Forebrain_E11.5_ENCFF144FZS_2_CpG.bed.gz
mv ENCFF170AZF.bed.gz Forebrain_E11.5_ENCFF170AZF_1_CHH.bed.gz
mv ENCFF211VME.bed.gz Forebrain_E11.5_ENCFF211VME_2_CHG.bed.gz
mv ENCFF359TMU.bed.gz Forebrain_E11.5_ENCFF359TMU_2_CHH.bed.gz
mv ENCFF485KXN.bed.gz Forebrain_E11.5_ENCFF485KXN_1_CHG.bed.gz
mv ENCFF530SDQ.bed.gz Forebrain_E11.5_ENCFF530SDQ_1_CpG.bed.gz
mv ENCFF078XSR.bed.gz Forebrain_E12.5_ENCFF078XSR_2_CpG.bed.gz
mv ENCFF140TCV.bed.gz Forebrain_E12.5_ENCFF140TCV_1_CpG.bed.gz
mv ENCFF318YVL.bed.gz Forebrain_E12.5_ENCFF318YVL_1_CHH.bed.gz
mv ENCFF326BGA.bed.gz Forebrain_E12.5_ENCFF326BGA_1_CpG.bed.gz
mv ENCFF355LMF.bed.gz Forebrain_E12.5_ENCFF355LMF_1_CHG.bed.gz
mv ENCFF356GTX.bed.gz Forebrain_E12.5_ENCFF356GTX_2_CpG.bed.gz
mv ENCFF386BDW.bed.gz Forebrain_E12.5_ENCFF386BDW_1_CHH.bed.gz
mv ENCFF619YYK.bed.gz Forebrain_E12.5_ENCFF619YYK_2_CHH.bed.gz
mv ENCFF728ZXA.bed.gz Forebrain_E12.5_ENCFF728ZXA_2_CHG.bed.gz
mv ENCFF861QZG.bed.gz Forebrain_E12.5_ENCFF861QZG_2_CHG.bed.gz
mv ENCFF880YON.bed.gz Forebrain_E12.5_ENCFF880YON_1_CHG.bed.gz
mv ENCFF984JGN.bed.gz Forebrain_E12.5_ENCFF984JGN_2_CHH.bed.gz
mv ENCFF098IMN.bed.gz Forebrain_E14.5_ENCFF098IMN_2_CHH.bed.gz
mv ENCFF287NGW.bed.gz Forebrain_E14.5_ENCFF287NGW_2_CHG.bed.gz
mv ENCFF321BQX.bed.gz Forebrain_E14.5_ENCFF321BQX_2_CpG.bed.gz
mv ENCFF811RWO.bed.gz Forebrain_E14.5_ENCFF811RWO_1_CHG.bed.gz
mv ENCFF873SBV.bed.gz Forebrain_E14.5_ENCFF873SBV_1_CpG.bed.gz
mv ENCFF913UYX.bed.gz Forebrain_E14.5_ENCFF913UYX_1_CHH.bed.gz
mv ENCFF146KIN.bed.gz Forebrain_E15.5_ENCFF146KIN_2_CpG.bed.gz
mv ENCFF163LYC.bed.gz Forebrain_E15.5_ENCFF163LYC_1_CHH.bed.gz
mv ENCFF197VSH.bed.gz Forebrain_E15.5_ENCFF197VSH_1_CHG.bed.gz
mv ENCFF277ADP.bed.gz Forebrain_E15.5_ENCFF277ADP_2_CHG.bed.gz
mv ENCFF310XUL.bed.gz Forebrain_E15.5_ENCFF310XUL_1_CpG.bed.gz
mv ENCFF626JUP.bed.gz Forebrain_E15.5_ENCFF626JUP_2_CHH.bed.gz
mv ENCFF072NIR.bed.gz Forebrain_P0_ENCFF072NIR_2_CHH.bed.gz
mv ENCFF184ZTA.bed.gz Forebrain_P0_ENCFF184ZTA_1_CHG.bed.gz
mv ENCFF446IJS.bed.gz Forebrain_P0_ENCFF446IJS_1_CpG.bed.gz
mv ENCFF783ZKP.bed.gz Forebrain_P0_ENCFF783ZKP_2_CpG.bed.gz
mv ENCFF848QLE.bed.gz Forebrain_P0_ENCFF848QLE_1_CHH.bed.gz
mv ENCFF890QOU.bed.gz Forebrain_P0_ENCFF890QOU_2_CHG.bed.gz
mv ENCFF167OJH.bed.gz Heart_E10.5_ENCFF167OJH_1_CHG.bed.gz
mv ENCFF327MVH.bed.gz Heart_E10.5_ENCFF327MVH_1_CHH.bed.gz
mv ENCFF428AXW.bed.gz Heart_E10.5_ENCFF428AXW_1_CpG.bed.gz
mv ENCFF677YTO.bed.gz Heart_E10.5_ENCFF677YTO_2_CHG.bed.gz
mv ENCFF751DLO.bed.gz Heart_E10.5_ENCFF751DLO_2_CpG.bed.gz
mv ENCFF945JPE.bed.gz Heart_E10.5_ENCFF945JPE_2_CHH.bed.gz
mv ENCFF008KDS.bed.gz Heart_E11.5_ENCFF008KDS_2_CpG.bed.gz
mv ENCFF107UNJ.bed.gz Heart_E11.5_ENCFF107UNJ_1_CpG.bed.gz
mv ENCFF466ROG.bed.gz Heart_E11.5_ENCFF466ROG_1_CHH.bed.gz
mv ENCFF527YFM.bed.gz Heart_E11.5_ENCFF527YFM_2_CHH.bed.gz
mv ENCFF745JIA.bed.gz Heart_E11.5_ENCFF745JIA_1_CHG.bed.gz
mv ENCFF764TIB.bed.gz Heart_E11.5_ENCFF764TIB_2_CHG.bed.gz
mv ENCFF075MVM.bed.gz Heart_E14.5_ENCFF075MVM_2_CHG.bed.gz
mv ENCFF203TRR.bed.gz Heart_E14.5_ENCFF203TRR_2_CpG.bed.gz
mv ENCFF384GAO.bed.gz Heart_E14.5_ENCFF384GAO_2_CHH.bed.gz
mv ENCFF624VWF.bed.gz Heart_E14.5_ENCFF624VWF_1_CHH.bed.gz
mv ENCFF925EPJ.bed.gz Heart_E14.5_ENCFF925EPJ_1_CpG.bed.gz
mv ENCFF997BOO.bed.gz Heart_E14.5_ENCFF997BOO_1_CHG.bed.gz
mv ENCFF050DUZ.bed.gz Heart_E15.5_ENCFF050DUZ_2_CHH.bed.gz
mv ENCFF116FND.bed.gz Heart_E15.5_ENCFF116FND_2_CHG.bed.gz
mv ENCFF698PML.bed.gz Heart_E15.5_ENCFF698PML_1_CHG.bed.gz
mv ENCFF951TJQ.bed.gz Heart_E15.5_ENCFF951TJQ_1_CpG.bed.gz
mv ENCFF963WCA.bed.gz Heart_E15.5_ENCFF963WCA_2_CpG.bed.gz
mv ENCFF966QQO.bed.gz Heart_E15.5_ENCFF966QQO_1_CHH.bed.gz
mv ENCFF091PCT.bed.gz Heart_P0_ENCFF091PCT_2_CHH.bed.gz
mv ENCFF095ULX.bed.gz Heart_P0_ENCFF095ULX_2_CHG.bed.gz
mv ENCFF096RYX.bed.gz Heart_P0_ENCFF096RYX_2_CpG.bed.gz
mv ENCFF467UEZ.bed.gz Heart_P0_ENCFF467UEZ_1_CpG.bed.gz
mv ENCFF641HYV.bed.gz Heart_P0_ENCFF641HYV_1_CHH.bed.gz
mv ENCFF800JTN.bed.gz Heart_P0_ENCFF800JTN_1_CHG.bed.gz
mv ENCFF193XZY.bed.gz Hindbrain_E14.5_ENCFF193XZY_1_CHG.bed.gz
mv ENCFF335KGS.bed.gz Hindbrain_E14.5_ENCFF335KGS_2_CHG.bed.gz
mv ENCFF362ENF.bed.gz Hindbrain_E14.5_ENCFF362ENF_2_CpG.bed.gz
mv ENCFF424YVZ.bed.gz Hindbrain_E14.5_ENCFF424YVZ_2_CHH.bed.gz
mv ENCFF434IOR.bed.gz Hindbrain_E14.5_ENCFF434IOR_1_CpG.bed.gz
mv ENCFF536XTB.bed.gz Hindbrain_E14.5_ENCFF536XTB_1_CHH.bed.gz
mv ENCFF084TEF.bed.gz Hindbrain_E15.5_ENCFF084TEF_2_CHH.bed.gz
mv ENCFF268GBN.bed.gz Hindbrain_E15.5_ENCFF268GBN_1_CHG.bed.gz
mv ENCFF279YQH.bed.gz Hindbrain_E15.5_ENCFF279YQH_2_CHG.bed.gz
mv ENCFF764EVM.bed.gz Hindbrain_E15.5_ENCFF764EVM_1_CHH.bed.gz
mv ENCFF941BTB.bed.gz Hindbrain_E15.5_ENCFF941BTB_2_CpG.bed.gz
mv ENCFF949ZVJ.bed.gz Hindbrain_E15.5_ENCFF949ZVJ_1_CpG.bed.gz
mv ENCFF163RKZ.bed.gz Hindbrain_E16.5_ENCFF163RKZ_1_CHH.bed.gz
mv ENCFF239RWZ.bed.gz Hindbrain_E16.5_ENCFF239RWZ_2_CHG.bed.gz
mv ENCFF407TGX.bed.gz Hindbrain_E16.5_ENCFF407TGX_2_CpG.bed.gz
mv ENCFF567XKD.bed.gz Hindbrain_E16.5_ENCFF567XKD_1_CHG.bed.gz
mv ENCFF835UNN.bed.gz Hindbrain_E16.5_ENCFF835UNN_2_CHH.bed.gz
mv ENCFF994CUN.bed.gz Hindbrain_E16.5_ENCFF994CUN_1_CpG.bed.gz
mv ENCFF303GFV.bed.gz Hindbrain_P0_ENCFF303GFV_1_CpG.bed.gz
mv ENCFF433NPI.bed.gz Hindbrain_P0_ENCFF433NPI_2_CpG.bed.gz
mv ENCFF487WLW.bed.gz Hindbrain_P0_ENCFF487WLW_1_CHH.bed.gz
mv ENCFF849GKH.bed.gz Hindbrain_P0_ENCFF849GKH_2_CHH.bed.gz
mv ENCFF861OKE.bed.gz Hindbrain_P0_ENCFF861OKE_1_CHG.bed.gz
mv ENCFF959XMN.bed.gz Hindbrain_P0_ENCFF959XMN_2_CHG.bed.gz
mv ENCFF082AEW.bed.gz Intestine_E15.5_ENCFF082AEW_1_CHH.bed.gz
mv ENCFF488REA.bed.gz Intestine_E15.5_ENCFF488REA_1_CpG.bed.gz
mv ENCFF578MUQ.bed.gz Intestine_E15.5_ENCFF578MUQ_2_CHH.bed.gz
mv ENCFF711PYT.bed.gz Intestine_E15.5_ENCFF711PYT_1_CHG.bed.gz
mv ENCFF874BVN.bed.gz Intestine_E15.5_ENCFF874BVN_2_CHG.bed.gz
mv ENCFF963EBD.bed.gz Intestine_E15.5_ENCFF963EBD_2_CpG.bed.gz
mv ENCFF331NBR.bed.gz Intestine_E16.5_ENCFF331NBR_1_CHG.bed.gz
mv ENCFF690ADB.bed.gz Intestine_E16.5_ENCFF690ADB_2_CpG.bed.gz
mv ENCFF888TCY.bed.gz Intestine_E16.5_ENCFF888TCY_1_CHH.bed.gz
mv ENCFF895RRK.bed.gz Intestine_E16.5_ENCFF895RRK_1_CpG.bed.gz
mv ENCFF899HXB.bed.gz Intestine_E16.5_ENCFF899HXB_2_CHG.bed.gz
mv ENCFF991NMN.bed.gz Intestine_E16.5_ENCFF991NMN_2_CHH.bed.gz
mv ENCFF036IYN.bed.gz Intestine_P0_ENCFF036IYN_2_CHG.bed.gz
mv ENCFF101IXI.bed.gz Intestine_P0_ENCFF101IXI_1_CHH.bed.gz
mv ENCFF498OHJ.bed.gz Intestine_P0_ENCFF498OHJ_2_CHH.bed.gz
mv ENCFF644XMA.bed.gz Intestine_P0_ENCFF644XMA_1_CpG.bed.gz
mv ENCFF714DLC.bed.gz Intestine_P0_ENCFF714DLC_1_CHG.bed.gz
mv ENCFF968JQK.bed.gz Intestine_P0_ENCFF968JQK_2_CpG.bed.gz
mv ENCFF078LCR.bed.gz Kidney_E14.5_ENCFF078LCR_1_CpG.bed.gz
mv ENCFF140ZSS.bed.gz Kidney_E14.5_ENCFF140ZSS_2_CHH.bed.gz
mv ENCFF558IPJ.bed.gz Kidney_E14.5_ENCFF558IPJ_1_CHG.bed.gz
mv ENCFF616VND.bed.gz Kidney_E14.5_ENCFF616VND_2_CHG.bed.gz
mv ENCFF783ZZK.bed.gz Kidney_E14.5_ENCFF783ZZK_2_CpG.bed.gz
mv ENCFF891CVJ.bed.gz Kidney_E14.5_ENCFF891CVJ_1_CHH.bed.gz
mv ENCFF098GFU.bed.gz Kidney_E16.5_ENCFF098GFU_1_CHG.bed.gz
mv ENCFF210JXM.bed.gz Kidney_E16.5_ENCFF210JXM_2_CHG.bed.gz
mv ENCFF407XCC.bed.gz Kidney_E16.5_ENCFF407XCC_1_CHH.bed.gz
mv ENCFF444WNU.bed.gz Kidney_E16.5_ENCFF444WNU_2_CpG.bed.gz
mv ENCFF478YVH.bed.gz Kidney_E16.5_ENCFF478YVH_1_CpG.bed.gz
mv ENCFF689XHV.bed.gz Kidney_E16.5_ENCFF689XHV_2_CHH.bed.gz
mv ENCFF055XUI.bed.gz Kidney_P0_ENCFF055XUI_1_CHH.bed.gz
mv ENCFF226ISA.bed.gz Kidney_P0_ENCFF226ISA_1_CpG.bed.gz
mv ENCFF459BQO.bed.gz Kidney_P0_ENCFF459BQO_1_CHG.bed.gz
mv ENCFF549SKO.bed.gz Kidney_P0_ENCFF549SKO_2_CHH.bed.gz
mv ENCFF745SMP.bed.gz Kidney_P0_ENCFF745SMP_2_CpG.bed.gz
mv ENCFF924VGN.bed.gz Kidney_P0_ENCFF924VGN_2_CHG.bed.gz
mv ENCFF230BWH.bed.gz Limb_E10.5_ENCFF230BWH_1_CHG.bed.gz
mv ENCFF245SIR.bed.gz Limb_E10.5_ENCFF245SIR_1_CHH.bed.gz
mv ENCFF572FTW.bed.gz Limb_E10.5_ENCFF572FTW_1_CpG.bed.gz
mv ENCFF423RKY.bed.gz Limb_E11.5_ENCFF423RKY_2_CHH.bed.gz
mv ENCFF567KXJ.bed.gz Limb_E11.5_ENCFF567KXJ_1_CHH.bed.gz
mv ENCFF644JJD.bed.gz Limb_E11.5_ENCFF644JJD_2_CpG.bed.gz
mv ENCFF789FQF.bed.gz Limb_E11.5_ENCFF789FQF_2_CHG.bed.gz
mv ENCFF800JCD.bed.gz Limb_E11.5_ENCFF800JCD_1_CpG.bed.gz
mv ENCFF899PXM.bed.gz Limb_E11.5_ENCFF899PXM_1_CHG.bed.gz
mv ENCFF038IVS.bed.gz Limb_E14.5_ENCFF038IVS_2_CpG.bed.gz
mv ENCFF147UAB.bed.gz Limb_E14.5_ENCFF147UAB_1_CHH.bed.gz
mv ENCFF154RZC.bed.gz Limb_E14.5_ENCFF154RZC_2_CHG.bed.gz
mv ENCFF176VLZ.bed.gz Limb_E14.5_ENCFF176VLZ_1_CpG.bed.gz
mv ENCFF488NDY.bed.gz Limb_E14.5_ENCFF488NDY_2_CHH.bed.gz
mv ENCFF989SSB.bed.gz Limb_E14.5_ENCFF989SSB_1_CHG.bed.gz
mv ENCFF342HSU.bed.gz Limb_E15.5_ENCFF342HSU_2_CHG.bed.gz
mv ENCFF567HPM.bed.gz Limb_E15.5_ENCFF567HPM_1_CHH.bed.gz
mv ENCFF613ASH.bed.gz Limb_E15.5_ENCFF613ASH_2_CpG.bed.gz
mv ENCFF724AJN.bed.gz Limb_E15.5_ENCFF724AJN_1_CpG.bed.gz
mv ENCFF727ANI.bed.gz Limb_E15.5_ENCFF727ANI_2_CHH.bed.gz
mv ENCFF792NCZ.bed.gz Limb_E15.5_ENCFF792NCZ_1_CHG.bed.gz
mv ENCFF306MCB.bed.gz Liver_E11.5_ENCFF306MCB_2_CpG.bed.gz
mv ENCFF567UEW.bed.gz Liver_E11.5_ENCFF567UEW_2_CHH.bed.gz
mv ENCFF772YDD.bed.gz Liver_E11.5_ENCFF772YDD_1_CHH.bed.gz
mv ENCFF900FQC.bed.gz Liver_E11.5_ENCFF900FQC_1_CHG.bed.gz
mv ENCFF903TZX.bed.gz Liver_E11.5_ENCFF903TZX_1_CpG.bed.gz
mv ENCFF904WEW.bed.gz Liver_E11.5_ENCFF904WEW_2_CHG.bed.gz
mv ENCFF079JQZ.bed.gz Liver_E14.5_ENCFF079JQZ_2_CpG.bed.gz
mv ENCFF288JXO.bed.gz Liver_E14.5_ENCFF288JXO_1_CpG.bed.gz
mv ENCFF322SWN.bed.gz Liver_E14.5_ENCFF322SWN_1_CHH.bed.gz
mv ENCFF324ZRK.bed.gz Liver_E14.5_ENCFF324ZRK_2_CHG.bed.gz
mv ENCFF633ZCK.bed.gz Liver_E14.5_ENCFF633ZCK_1_CHG.bed.gz
mv ENCFF925DUI.bed.gz Liver_E14.5_ENCFF925DUI_2_CHH.bed.gz
mv ENCFF036LQW.bed.gz Liver_E15.5_ENCFF036LQW_1_CHG.bed.gz
mv ENCFF354TYJ.bed.gz Liver_E15.5_ENCFF354TYJ_1_CHH.bed.gz
mv ENCFF371FSG.bed.gz Liver_E15.5_ENCFF371FSG_2_CpG.bed.gz
mv ENCFF553AFX.bed.gz Liver_E15.5_ENCFF553AFX_1_CpG.bed.gz
mv ENCFF756ITF.bed.gz Liver_E15.5_ENCFF756ITF_2_CHH.bed.gz
mv ENCFF758ZHO.bed.gz Liver_E15.5_ENCFF758ZHO_2_CHG.bed.gz
mv ENCFF243ELJ.bed.gz Liver_E16.5_ENCFF243ELJ_1_CpG.bed.gz
mv ENCFF363NVM.bed.gz Liver_E16.5_ENCFF363NVM_2_CpG.bed.gz
mv ENCFF834POI.bed.gz Liver_E16.5_ENCFF834POI_1_CHH.bed.gz
mv ENCFF849TAC.bed.gz Liver_E16.5_ENCFF849TAC_1_CHG.bed.gz
mv ENCFF852HTJ.bed.gz Liver_E16.5_ENCFF852HTJ_2_CHH.bed.gz
mv ENCFF963ROS.bed.gz Liver_E16.5_ENCFF963ROS_2_CHG.bed.gz
mv ENCFF054AYU.bed.gz Lung_E14.5_ENCFF054AYU_2_CpG.bed.gz
mv ENCFF070FVI.bed.gz Lung_E14.5_ENCFF070FVI_1_CpG.bed.gz
mv ENCFF167VJT.bed.gz Lung_E14.5_ENCFF167VJT_1_CHG.bed.gz
mv ENCFF291LSY.bed.gz Lung_E14.5_ENCFF291LSY_2_CHH.bed.gz
mv ENCFF488VTG.bed.gz Lung_E14.5_ENCFF488VTG_2_CHG.bed.gz
mv ENCFF815VIM.bed.gz Lung_E14.5_ENCFF815VIM_1_CHH.bed.gz
mv ENCFF117RRX.bed.gz Lung_E15.5_ENCFF117RRX_1_CHG.bed.gz
mv ENCFF138NBE.bed.gz Lung_E15.5_ENCFF138NBE_1_CpG.bed.gz
mv ENCFF417IKO.bed.gz Lung_E15.5_ENCFF417IKO_1_CHH.bed.gz
mv ENCFF589YIF.bed.gz Lung_E15.5_ENCFF589YIF_2_CHH.bed.gz
mv ENCFF775PQC.bed.gz Lung_E15.5_ENCFF775PQC_2_CHG.bed.gz
mv ENCFF922AMS.bed.gz Lung_E15.5_ENCFF922AMS_2_CpG.bed.gz
mv ENCFF054SJJ.bed.gz Lung_E16.5_ENCFF054SJJ_1_CHG.bed.gz
mv ENCFF222XYK.bed.gz Lung_E16.5_ENCFF222XYK_2_CpG.bed.gz
mv ENCFF359CFG.bed.gz Lung_E16.5_ENCFF359CFG_1_CpG.bed.gz
mv ENCFF412QLO.bed.gz Lung_E16.5_ENCFF412QLO_1_CHH.bed.gz
mv ENCFF448MFG.bed.gz Lung_E16.5_ENCFF448MFG_2_CHH.bed.gz
mv ENCFF891UGL.bed.gz Lung_E16.5_ENCFF891UGL_2_CHG.bed.gz
mv ENCFF086MTL.bed.gz Lung_P0_ENCFF086MTL_1_CHH.bed.gz
mv ENCFF441XUW.bed.gz Lung_P0_ENCFF441XUW_2_CHH.bed.gz
mv ENCFF555BBK.bed.gz Lung_P0_ENCFF555BBK_1_CHG.bed.gz
mv ENCFF555SRI.bed.gz Lung_P0_ENCFF555SRI_2_CpG.bed.gz
mv ENCFF920WHY.bed.gz Lung_P0_ENCFF920WHY_1_CpG.bed.gz
mv ENCFF933BGO.bed.gz Lung_P0_ENCFF933BGO_2_CHG.bed.gz
mv ENCFF022CRO.bed.gz Midbrain_E10.5_ENCFF022CRO_2_CHG.bed.gz
mv ENCFF318LFN.bed.gz Midbrain_E10.5_ENCFF318LFN_1_CHH.bed.gz
mv ENCFF515VOU.bed.gz Midbrain_E10.5_ENCFF515VOU_2_CHH.bed.gz
mv ENCFF672VOP.bed.gz Midbrain_E10.5_ENCFF672VOP_1_CpG.bed.gz
mv ENCFF735TII.bed.gz Midbrain_E10.5_ENCFF735TII_1_CHG.bed.gz
mv ENCFF995WMY.bed.gz Midbrain_E10.5_ENCFF995WMY_2_CpG.bed.gz
mv ENCFF061VBW.bed.gz Midbrain_E11.5_ENCFF061VBW_1_CHH.bed.gz
mv ENCFF584DMO.bed.gz Midbrain_E11.5_ENCFF584DMO_1_CHG.bed.gz
mv ENCFF611XAY.bed.gz Midbrain_E11.5_ENCFF611XAY_1_CpG.bed.gz
mv ENCFF750DYD.bed.gz Midbrain_E11.5_ENCFF750DYD_2_CpG.bed.gz
mv ENCFF982XEC.bed.gz Midbrain_E11.5_ENCFF982XEC_2_CHH.bed.gz
mv ENCFF983JGH.bed.gz Midbrain_E11.5_ENCFF983JGH_2_CHG.bed.gz
mv ENCFF004ZOY.bed.gz Midbrain_E14.5_ENCFF004ZOY_2_CHG.bed.gz
mv ENCFF575NQV.bed.gz Midbrain_E14.5_ENCFF575NQV_2_CHH.bed.gz
mv ENCFF673VYW.bed.gz Midbrain_E14.5_ENCFF673VYW_1_CHG.bed.gz
mv ENCFF793ZHN.bed.gz Midbrain_E14.5_ENCFF793ZHN_2_CpG.bed.gz
mv ENCFF885FQW.bed.gz Midbrain_E14.5_ENCFF885FQW_1_CHH.bed.gz
mv ENCFF897NOL.bed.gz Midbrain_E14.5_ENCFF897NOL_1_CpG.bed.gz
mv ENCFF320GGU.bed.gz Midbrain_E15.5_ENCFF320GGU_1_CHH.bed.gz
mv ENCFF426XRG.bed.gz Midbrain_E15.5_ENCFF426XRG_2_CpG.bed.gz
mv ENCFF647NFJ.bed.gz Midbrain_E15.5_ENCFF647NFJ_1_CHG.bed.gz
mv ENCFF688RHG.bed.gz Midbrain_E15.5_ENCFF688RHG_1_CpG.bed.gz
mv ENCFF960ADL.bed.gz Midbrain_E15.5_ENCFF960ADL_2_CHH.bed.gz
mv ENCFF985GEB.bed.gz Midbrain_E15.5_ENCFF985GEB_2_CHG.bed.gz
mv ENCFF370QDC.bed.gz Midbrain_P0_ENCFF370QDC_1_CHH.bed.gz
mv ENCFF504HKC.bed.gz Midbrain_P0_ENCFF504HKC_1_CpG.bed.gz
mv ENCFF597NPZ.bed.gz Midbrain_P0_ENCFF597NPZ_2_CHH.bed.gz
mv ENCFF747LEB.bed.gz Midbrain_P0_ENCFF747LEB_1_CHG.bed.gz
mv ENCFF933NAQ.bed.gz Midbrain_P0_ENCFF933NAQ_2_CHG.bed.gz
mv ENCFF966UJE.bed.gz Midbrain_P0_ENCFF966UJE_2_CpG.bed.gz
mv ENCFF572VAP.bed.gz NeuralTube_E11.5_ENCFF572VAP_1_CHG.bed.gz
mv ENCFF793ATY.bed.gz NeuralTube_E11.5_ENCFF793ATY_2_CHH.bed.gz
mv ENCFF859JHQ.bed.gz NeuralTube_E11.5_ENCFF859JHQ_1_CpG.bed.gz
mv ENCFF870TZB.bed.gz NeuralTube_E11.5_ENCFF870TZB_2_CHG.bed.gz
mv ENCFF933JRY.bed.gz NeuralTube_E11.5_ENCFF933JRY_1_CHH.bed.gz
mv ENCFF943HMN.bed.gz NeuralTube_E11.5_ENCFF943HMN_2_CpG.bed.gz
mv ENCFF152IRW.bed.gz NeuralTube_E14.5_ENCFF152IRW_2_CpG.bed.gz
mv ENCFF197WQW.bed.gz NeuralTube_E14.5_ENCFF197WQW_2_CHH.bed.gz
mv ENCFF237MDB.bed.gz NeuralTube_E14.5_ENCFF237MDB_1_CHG.bed.gz
mv ENCFF622LPY.bed.gz NeuralTube_E14.5_ENCFF622LPY_2_CHG.bed.gz
mv ENCFF937GFF.bed.gz NeuralTube_E14.5_ENCFF937GFF_1_CpG.bed.gz
mv ENCFF944KCR.bed.gz NeuralTube_E14.5_ENCFF944KCR_1_CHH.bed.gz
mv ENCFF042EKD.bed.gz NeuralTube_E15.5_ENCFF042EKD_1_CHH.bed.gz
mv ENCFF114BEJ.bed.gz NeuralTube_E15.5_ENCFF114BEJ_2_CHG.bed.gz
mv ENCFF115JNG.bed.gz NeuralTube_E15.5_ENCFF115JNG_1_CpG.bed.gz
mv ENCFF192NIG.bed.gz NeuralTube_E15.5_ENCFF192NIG_2_CpG.bed.gz
mv ENCFF364IIS.bed.gz NeuralTube_E15.5_ENCFF364IIS_1_CHG.bed.gz
mv ENCFF736MNI.bed.gz NeuralTube_E15.5_ENCFF736MNI_2_CHH.bed.gz
mv ENCFF176GYO.bed.gz Stomach_E15.5_ENCFF176GYO_2_CHH.bed.gz
mv ENCFF284YYM.bed.gz Stomach_E15.5_ENCFF284YYM_1_CHH.bed.gz
mv ENCFF301QSZ.bed.gz Stomach_E15.5_ENCFF301QSZ_2_CpG.bed.gz
mv ENCFF455WPB.bed.gz Stomach_E15.5_ENCFF455WPB_2_CHG.bed.gz
mv ENCFF748HYK.bed.gz Stomach_E15.5_ENCFF748HYK_1_CHG.bed.gz
mv ENCFF880SXA.bed.gz Stomach_E15.5_ENCFF880SXA_1_CpG.bed.gz
mv ENCFF051OZQ.bed.gz Stomach_P0_ENCFF051OZQ_1_CpG.bed.gz
mv ENCFF184VZB.bed.gz Stomach_P0_ENCFF184VZB_2_CHH.bed.gz
mv ENCFF349RJN.bed.gz Stomach_P0_ENCFF349RJN_1_CHG.bed.gz
mv ENCFF361ACS.bed.gz Stomach_P0_ENCFF361ACS_2_CHG.bed.gz
mv ENCFF562RGL.bed.gz Stomach_P0_ENCFF562RGL_2_CpG.bed.gz
mv ENCFF782FDT.bed.gz Stomach_P0_ENCFF782FDT_1_CHH.bed.gz

mv ENCFF004RDJ.bed.gz EmbryonicFacialProminence_E12.5_ENCFF004RDJ_2_CHG.bed.gz
mv ENCFF132PIM.bed.gz EmbryonicFacialProminence_E12.5_ENCFF132PIM_2_CHH.bed.gz
mv ENCFF497NGX.bed.gz EmbryonicFacialProminence_E12.5_ENCFF497NGX_2_CpG.bed.gz
mv ENCFF698HGV.bed.gz EmbryonicFacialProminence_E12.5_ENCFF698HGV_1_CpG.bed.gz
mv ENCFF836OZV.bed.gz EmbryonicFacialProminence_E12.5_ENCFF836OZV_1_CHG.bed.gz
mv ENCFF881TBC.bed.gz EmbryonicFacialProminence_E12.5_ENCFF881TBC_1_CHH.bed.gz
mv ENCFF176SNG.bed.gz EmbryonicFacialProminence_E13.5_ENCFF176SNG_1_CHG.bed.gz
mv ENCFF235LSX.bed.gz EmbryonicFacialProminence_E13.5_ENCFF235LSX_2_CpG.bed.gz
mv ENCFF311GQV.bed.gz EmbryonicFacialProminence_E13.5_ENCFF311GQV_1_CpG.bed.gz
mv ENCFF495ICZ.bed.gz EmbryonicFacialProminence_E13.5_ENCFF495ICZ_2_CHH.bed.gz
mv ENCFF901KED.bed.gz EmbryonicFacialProminence_E13.5_ENCFF901KED_1_CHH.bed.gz
mv ENCFF915RPR.bed.gz EmbryonicFacialProminence_E13.5_ENCFF915RPR_2_CHG.bed.gz
mv ENCFF028HZV.bed.gz EmbryonicFacialProminence_E15.5_ENCFF028HZV_1_CHH.bed.gz
mv ENCFF329QCJ.bed.gz EmbryonicFacialProminence_E15.5_ENCFF329QCJ_1_CHG.bed.gz
mv ENCFF396QSA.bed.gz EmbryonicFacialProminence_E15.5_ENCFF396QSA_2_CpG.bed.gz
mv ENCFF559ZAX.bed.gz EmbryonicFacialProminence_E15.5_ENCFF559ZAX_2_CHH.bed.gz
mv ENCFF581CMO.bed.gz EmbryonicFacialProminence_E15.5_ENCFF581CMO_1_CpG.bed.gz
mv ENCFF820XFX.bed.gz EmbryonicFacialProminence_E15.5_ENCFF820XFX_2_CHG.bed.gz
mv ENCFF060VXB.bed.gz Forebrain_E16.5_ENCFF060VXB_1_CHG.bed.gz
mv ENCFF098QRE.bed.gz Forebrain_E16.5_ENCFF098QRE_2_CpG.bed.gz
mv ENCFF533JUK.bed.gz Forebrain_E16.5_ENCFF533JUK_2_CHG.bed.gz
mv ENCFF592SYK.bed.gz Forebrain_E16.5_ENCFF592SYK_1_CpG.bed.gz
mv ENCFF831GBW.bed.gz Forebrain_E16.5_ENCFF831GBW_2_CHH.bed.gz
mv ENCFF972VAG.bed.gz Forebrain_E16.5_ENCFF972VAG_1_CHH.bed.gz
mv ENCFF002TAP.bed.gz Heart_E12.5_ENCFF002TAP_2_CHG.bed.gz
mv ENCFF179OQV.bed.gz Heart_E12.5_ENCFF179OQV_2_CHH.bed.gz
mv ENCFF464PLY.bed.gz Heart_E12.5_ENCFF464PLY_2_CpG.bed.gz
mv ENCFF519VEW.bed.gz Heart_E12.5_ENCFF519VEW_1_CHG.bed.gz
mv ENCFF792LFL.bed.gz Heart_E12.5_ENCFF792LFL_1_CHH.bed.gz
mv ENCFF816XTJ.bed.gz Heart_E12.5_ENCFF816XTJ_1_CpG.bed.gz
mv ENCFF114EIZ.bed.gz Heart_E13.5_ENCFF114EIZ_1_CpG.bed.gz
mv ENCFF391LHX.bed.gz Heart_E13.5_ENCFF391LHX_2_CHH.bed.gz
mv ENCFF711PFS.bed.gz Heart_E13.5_ENCFF711PFS_2_CpG.bed.gz
mv ENCFF800YTE.bed.gz Heart_E13.5_ENCFF800YTE_1_CHH.bed.gz
mv ENCFF832CEM.bed.gz Heart_E13.5_ENCFF832CEM_1_CHG.bed.gz
mv ENCFF897PZN.bed.gz Heart_E13.5_ENCFF897PZN_2_CHG.bed.gz
mv ENCFF171RUF.bed.gz Heart_E16.5_ENCFF171RUF_2_CHG.bed.gz
mv ENCFF367TSY.bed.gz Heart_E16.5_ENCFF367TSY_1_CHG.bed.gz
mv ENCFF574RDQ.bed.gz Heart_E16.5_ENCFF574RDQ_1_CpG.bed.gz
mv ENCFF633ZDN.bed.gz Heart_E16.5_ENCFF633ZDN_2_CpG.bed.gz
mv ENCFF908XGD.bed.gz Heart_E16.5_ENCFF908XGD_2_CHH.bed.gz
mv ENCFF942RGX.bed.gz Heart_E16.5_ENCFF942RGX_1_CHH.bed.gz
mv ENCFF424GSL.bed.gz Hindbrain_E10.5_ENCFF424GSL_1_CHG.bed.gz
mv ENCFF474FXF.bed.gz Hindbrain_E10.5_ENCFF474FXF_2_CHG.bed.gz
mv ENCFF507ENU.bed.gz Hindbrain_E10.5_ENCFF507ENU_1_CpG.bed.gz
mv ENCFF550IDR.bed.gz Hindbrain_E10.5_ENCFF550IDR_2_CpG.bed.gz
mv ENCFF742BEA.bed.gz Hindbrain_E10.5_ENCFF742BEA_2_CHH.bed.gz
mv ENCFF755NTZ.bed.gz Hindbrain_E10.5_ENCFF755NTZ_1_CHH.bed.gz
mv ENCFF069YRI.bed.gz Hindbrain_E11.5_ENCFF069YRI_1_CHH.bed.gz
mv ENCFF159VIE.bed.gz Hindbrain_E11.5_ENCFF159VIE_2_CpG.bed.gz
mv ENCFF223ZGM.bed.gz Hindbrain_E11.5_ENCFF223ZGM_2_CHG.bed.gz
mv ENCFF246BSY.bed.gz Hindbrain_E11.5_ENCFF246BSY_2_CHH.bed.gz
mv ENCFF279GQT.bed.gz Hindbrain_E11.5_ENCFF279GQT_1_CpG.bed.gz
mv ENCFF440PTK.bed.gz Hindbrain_E11.5_ENCFF440PTK_1_CHG.bed.gz
mv ENCFF172LOQ.bed.gz Hindbrain_E12.5_ENCFF172LOQ_1_CpG.bed.gz
mv ENCFF189KVN.bed.gz Hindbrain_E12.5_ENCFF189KVN_2_CHG.bed.gz
mv ENCFF476XYT.bed.gz Hindbrain_E12.5_ENCFF476XYT_1_CHH.bed.gz
mv ENCFF599BKR.bed.gz Hindbrain_E12.5_ENCFF599BKR_2_CpG.bed.gz
mv ENCFF769XAN.bed.gz Hindbrain_E12.5_ENCFF769XAN_2_CHH.bed.gz
mv ENCFF941FFK.bed.gz Hindbrain_E12.5_ENCFF941FFK_1_CHG.bed.gz
mv ENCFF283KKF.bed.gz Hindbrain_E13.5_ENCFF283KKF_1_CHH.bed.gz
mv ENCFF467JWP.bed.gz Hindbrain_E13.5_ENCFF467JWP_1_CpG.bed.gz
mv ENCFF548NAI.bed.gz Hindbrain_E13.5_ENCFF548NAI_1_CHG.bed.gz
mv ENCFF585WMY.bed.gz Hindbrain_E13.5_ENCFF585WMY_2_CHH.bed.gz
mv ENCFF894POI.bed.gz Hindbrain_E13.5_ENCFF894POI_2_CpG.bed.gz
mv ENCFF923CEU.bed.gz Hindbrain_E13.5_ENCFF923CEU_2_CHG.bed.gz
mv ENCFF068EHZ.bed.gz Intestine_E14.5_ENCFF068EHZ_2_CHH.bed.gz
mv ENCFF232RQE.bed.gz Intestine_E14.5_ENCFF232RQE_2_CpG.bed.gz
mv ENCFF550FRR.bed.gz Intestine_E14.5_ENCFF550FRR_2_CHG.bed.gz
mv ENCFF550JOH.bed.gz Intestine_E14.5_ENCFF550JOH_1_CpG.bed.gz
mv ENCFF597EAD.bed.gz Intestine_E14.5_ENCFF597EAD_1_CHG.bed.gz
mv ENCFF811VVC.bed.gz Intestine_E14.5_ENCFF811VVC_1_CHH.bed.gz
mv ENCFF228AMN.bed.gz Kidney_E15.5_ENCFF228AMN_2_CHH.bed.gz
mv ENCFF253NMT.bed.gz Kidney_E15.5_ENCFF253NMT_1_CHH.bed.gz
mv ENCFF354NPD.bed.gz Kidney_E15.5_ENCFF354NPD_2_CHG.bed.gz
mv ENCFF364RJZ.bed.gz Kidney_E15.5_ENCFF364RJZ_2_CpG.bed.gz
mv ENCFF399WHI.bed.gz Kidney_E15.5_ENCFF399WHI_1_CHG.bed.gz
mv ENCFF420PMS.bed.gz Kidney_E15.5_ENCFF420PMS_1_CpG.bed.gz
mv ENCFF558LRK.bed.gz Limb_E10.5_ENCFF558LRK_2_CpG.bed.gz
mv ENCFF607KLP.bed.gz Limb_E10.5_ENCFF607KLP_2_CHG.bed.gz
mv ENCFF925AWP.bed.gz Limb_E10.5_ENCFF925AWP_2_CHH.bed.gz
mv ENCFF080DSY.bed.gz Limb_E12.5_ENCFF080DSY_2_CpG.bed.gz
mv ENCFF103DRR.bed.gz Limb_E12.5_ENCFF103DRR_1_CpG.bed.gz
mv ENCFF116OZB.bed.gz Limb_E12.5_ENCFF116OZB_1_CHH.bed.gz
mv ENCFF221AHX.bed.gz Limb_E12.5_ENCFF221AHX_2_CHG.bed.gz
mv ENCFF723SRT.bed.gz Limb_E12.5_ENCFF723SRT_1_CHG.bed.gz
mv ENCFF944MOE.bed.gz Limb_E12.5_ENCFF944MOE_2_CHH.bed.gz
mv ENCFF043XQY.bed.gz Limb_E13.5_ENCFF043XQY_1_CHG.bed.gz
mv ENCFF213NJE.bed.gz Limb_E13.5_ENCFF213NJE_2_CHH.bed.gz
mv ENCFF404MCS.bed.gz Limb_E13.5_ENCFF404MCS_1_CHH.bed.gz
mv ENCFF488GSF.bed.gz Limb_E13.5_ENCFF488GSF_2_CHG.bed.gz
mv ENCFF556ENA.bed.gz Limb_E13.5_ENCFF556ENA_1_CpG.bed.gz
mv ENCFF834PXN.bed.gz Limb_E13.5_ENCFF834PXN_2_CpG.bed.gz
mv ENCFF080NOK.bed.gz Liver_E12.5_ENCFF080NOK_1_CHH.bed.gz
mv ENCFF134ECD.bed.gz Liver_E12.5_ENCFF134ECD_2_CpG.bed.gz
mv ENCFF312AJC.bed.gz Liver_E12.5_ENCFF312AJC_2_CHG.bed.gz
mv ENCFF368QJJ.bed.gz Liver_E12.5_ENCFF368QJJ_1_CHG.bed.gz
mv ENCFF624JYN.bed.gz Liver_E12.5_ENCFF624JYN_1_CpG.bed.gz
mv ENCFF869DDJ.bed.gz Liver_E12.5_ENCFF869DDJ_2_CHH.bed.gz
mv ENCFF178WUX.bed.gz Liver_E13.5_ENCFF178WUX_2_CpG.bed.gz
mv ENCFF195RVI.bed.gz Liver_E13.5_ENCFF195RVI_1_CHG.bed.gz
mv ENCFF516VJL.bed.gz Liver_E13.5_ENCFF516VJL_1_CpG.bed.gz
mv ENCFF658PEL.bed.gz Liver_E13.5_ENCFF658PEL_1_CHH.bed.gz
mv ENCFF672FLQ.bed.gz Liver_E13.5_ENCFF672FLQ_2_CHH.bed.gz
mv ENCFF771OGD.bed.gz Liver_E13.5_ENCFF771OGD_2_CHG.bed.gz
mv ENCFF195LEN.bed.gz Liver_P0_ENCFF195LEN_2_CpG.bed.gz
mv ENCFF221DYF.bed.gz Liver_P0_ENCFF221DYF_2_CHG.bed.gz
mv ENCFF347TCD.bed.gz Liver_P0_ENCFF347TCD_1_CHH.bed.gz
mv ENCFF448ZJB.bed.gz Liver_P0_ENCFF448ZJB_2_CHH.bed.gz
mv ENCFF770AUO.bed.gz Liver_P0_ENCFF770AUO_1_CpG.bed.gz
mv ENCFF969GPP.bed.gz Liver_P0_ENCFF969GPP_1_CHG.bed.gz
mv ENCFF390NTN.bed.gz Midbrain_E12.5_ENCFF390NTN_1_CpG.bed.gz
mv ENCFF393VFF.bed.gz Midbrain_E12.5_ENCFF393VFF_2_CpG.bed.gz
mv ENCFF431WCK.bed.gz Midbrain_E12.5_ENCFF431WCK_2_CHH.bed.gz
mv ENCFF557RBS.bed.gz Midbrain_E12.5_ENCFF557RBS_2_CHG.bed.gz
mv ENCFF787CPG.bed.gz Midbrain_E12.5_ENCFF787CPG_1_CHH.bed.gz
mv ENCFF801KFY.bed.gz Midbrain_E12.5_ENCFF801KFY_1_CHG.bed.gz
mv ENCFF118LMO.bed.gz Midbrain_E13.5_ENCFF118LMO_2_CHH.bed.gz
mv ENCFF326GKU.bed.gz Midbrain_E13.5_ENCFF326GKU_1_CpG.bed.gz
mv ENCFF400NZS.bed.gz Midbrain_E13.5_ENCFF400NZS_2_CHG.bed.gz
mv ENCFF403IYQ.bed.gz Midbrain_E13.5_ENCFF403IYQ_1_CHH.bed.gz
mv ENCFF537UAE.bed.gz Midbrain_E13.5_ENCFF537UAE_1_CHG.bed.gz
mv ENCFF646TSJ.bed.gz Midbrain_E13.5_ENCFF646TSJ_2_CpG.bed.gz
mv ENCFF320XZV.bed.gz Midbrain_E16.5_ENCFF320XZV_2_CpG.bed.gz
mv ENCFF353LYL.bed.gz Midbrain_E16.5_ENCFF353LYL_2_CHG.bed.gz
mv ENCFF379CTK.bed.gz Midbrain_E16.5_ENCFF379CTK_2_CHH.bed.gz
mv ENCFF461HQE.bed.gz Midbrain_E16.5_ENCFF461HQE_1_CHH.bed.gz
mv ENCFF796POX.bed.gz Midbrain_E16.5_ENCFF796POX_1_CpG.bed.gz
mv ENCFF875BEE.bed.gz Midbrain_E16.5_ENCFF875BEE_1_CHG.bed.gz
mv ENCFF222CLY.bed.gz NeuralTube_E12.5_ENCFF222CLY_2_CHH.bed.gz
mv ENCFF231BLL.bed.gz NeuralTube_E12.5_ENCFF231BLL_1_CpG.bed.gz
mv ENCFF369BXZ.bed.gz NeuralTube_E12.5_ENCFF369BXZ_1_CHH.bed.gz
mv ENCFF697TML.bed.gz NeuralTube_E12.5_ENCFF697TML_2_CpG.bed.gz
mv ENCFF749BME.bed.gz NeuralTube_E12.5_ENCFF749BME_2_CHG.bed.gz
mv ENCFF864WRU.bed.gz NeuralTube_E12.5_ENCFF864WRU_1_CHG.bed.gz
mv ENCFF101JDR.bed.gz NeuralTube_E13.5_ENCFF101JDR_2_CpG.bed.gz
mv ENCFF152ZUF.bed.gz NeuralTube_E13.5_ENCFF152ZUF_1_CpG.bed.gz
mv ENCFF165GRG.bed.gz NeuralTube_E13.5_ENCFF165GRG_2_CHG.bed.gz
mv ENCFF548WLP.bed.gz NeuralTube_E13.5_ENCFF548WLP_2_CHH.bed.gz
mv ENCFF658VBN.bed.gz NeuralTube_E13.5_ENCFF658VBN_1_CHH.bed.gz
mv ENCFF693QZX.bed.gz NeuralTube_E13.5_ENCFF693QZX_1_CHG.bed.gz
mv ENCFF029NDO.bed.gz Stomach_E14.5_ENCFF029NDO_2_CpG.bed.gz
mv ENCFF369DBO.bed.gz Stomach_E14.5_ENCFF369DBO_2_CHG.bed.gz
mv ENCFF389YPY.bed.gz Stomach_E14.5_ENCFF389YPY_1_CHG.bed.gz
mv ENCFF720OAC.bed.gz Stomach_E14.5_ENCFF720OAC_1_CpG.bed.gz
mv ENCFF760PKX.bed.gz Stomach_E14.5_ENCFF760PKX_1_CHH.bed.gz
mv ENCFF840GNZ.bed.gz Stomach_E14.5_ENCFF840GNZ_2_CHH.bed.gz
mv ENCFF163CYY.bed.gz Stomach_E16.5_ENCFF163CYY_1_CHG.bed.gz
mv ENCFF230NVJ.bed.gz Stomach_E16.5_ENCFF230NVJ_2_CHG.bed.gz
mv ENCFF260QXS.bed.gz Stomach_E16.5_ENCFF260QXS_1_CHH.bed.gz
mv ENCFF532SZJ.bed.gz Stomach_E16.5_ENCFF532SZJ_1_CpG.bed.gz
mv ENCFF659YFO.bed.gz Stomach_E16.5_ENCFF659YFO_2_CHH.bed.gz
mv ENCFF803ATF.bed.gz Stomach_E16.5_ENCFF803ATF_2_CpG.bed.gz

```

## extract chr19 info

```sh
cd ~/data

for bed in *.bed.gz
do
    bname=${bed%.bed.gz}
    echo $bed $bname
    zcat $bed | awk -vOFS="\t" '{if ($1 == "chr19" && $10 > 9) print $1, $2, $3, $11, $10, $6 }' > $bname.chr19.bed &
done

wc -l *.bed
#    16776602 EmbryonicFacialProminence_E10.5_ENCFF006QDP_1_CHH.chr19.bed
#    15432189 EmbryonicFacialProminence_E10.5_ENCFF050FNM_2_CHH.chr19.bed
#      911351 EmbryonicFacialProminence_E10.5_ENCFF095VLL_2_CpG.chr19.bed
#     4637824 EmbryonicFacialProminence_E10.5_ENCFF198NON_2_CHG.chr19.bed
#     4990059 EmbryonicFacialProminence_E10.5_ENCFF722PJA_1_CHG.chr19.bed
#     1004429 EmbryonicFacialProminence_E10.5_ENCFF886MNA_1_CpG.chr19.bed
#      940055 EmbryonicFacialProminence_E11.5_ENCFF083XWZ_2_CpG.chr19.bed
#      964416 EmbryonicFacialProminence_E11.5_ENCFF306DKT_1_CpG.chr19.bed
#     4834729 EmbryonicFacialProminence_E11.5_ENCFF590MRS_1_CHG.chr19.bed
#    16076330 EmbryonicFacialProminence_E11.5_ENCFF754BON_2_CHH.chr19.bed
#    16376632 EmbryonicFacialProminence_E11.5_ENCFF884JVJ_1_CHH.chr19.bed
#     4745384 EmbryonicFacialProminence_E11.5_ENCFF997OLQ_2_CHG.chr19.bed
#    12017242 EmbryonicFacialProminence_E14.5_ENCFF063YAO_1_CHH.chr19.bed
#      815174 EmbryonicFacialProminence_E14.5_ENCFF686ZQQ_2_CpG.chr19.bed
#    14350517 EmbryonicFacialProminence_E14.5_ENCFF709AXF_2_CHH.chr19.bed
#     3528743 EmbryonicFacialProminence_E14.5_ENCFF744VPT_1_CHG.chr19.bed
#     4241469 EmbryonicFacialProminence_E14.5_ENCFF757MSK_2_CHG.chr19.bed
#      670539 EmbryonicFacialProminence_E14.5_ENCFF989UNH_1_CpG.chr19.bed
#      919351 Forebrain_E10.5_ENCFF365XZL_2_CpG.chr19.bed
#     4628308 Forebrain_E10.5_ENCFF369TZO_1_CHG.chr19.bed
#    15512425 Forebrain_E10.5_ENCFF400MQF_2_CHH.chr19.bed
#     4662702 Forebrain_E10.5_ENCFF506SUF_2_CHG.chr19.bed
#    15390436 Forebrain_E10.5_ENCFF590OFM_1_CHH.chr19.bed
#      909495 Forebrain_E10.5_ENCFF977BKF_1_CpG.chr19.bed

```

## prepare mm10 reference

```sh
cd ~/ref
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz

tar -zxvf Mus_musculus_UCSC_mm10.tar.gz

# genome.fa.fai IS older than genome.fa, generate new index
cd ~/ref/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/
samtools faidx genome.fa
```

## extract context information around C (±3)

```sh
cd ~/data

ref=~/ref/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa

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
tableCat.py -i *ENCFF*CpG*.context.txt -r .context.txt | awk '{split($8,a,"_"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"a[1]"\t"a[2]"\t"a[4]}' > mouse_tissue.CpG.chr19.context.txt &
# takes 6.5G of space

# concatenate all CHG context files
tableCat.py -i *ENCFF*CHG*.context.txt -r .context.txt | awk '{split($8,a,"_"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"a[1]"\t"a[2]"\t"a[4]}' > mouse_tissue.CHG.chr19.context.txt &
# take 33G of space

# concatenate all CHH context files
tableCat.py -i *ENCFF*CHH*.context.txt -r .context.txt | awk '{split($8,a,"_"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"a[1]"\t"a[2]"\t"a[4]}' > mouse_tissue.CHH.chr19.context.txt &
# take 110G of space

# remove unneeded columns to save space and memory for later analysis
# concatenate all context files, final size 70G
tableCat.py -i mouse_tissue.CH*.chr19.context.txt -r .context.txt | awk '{print $7"\t"$4"\t"$8"\t"$9"\t"$10}' > mouse_tissue.nonCpG.chr19.context.txt &

# only keep methylated cytosine for later analysis
# concatenate all context files, final size 8.3G
cat mouse_tissue.{CpG,CHG,CHH}.chr19.context.txt | awk '{if ($4 > 0) print $7"\t"$4"\t"$8"\t"$9"\t"$10}' > mouse_tissue.methylated.chr19.context.txt &

# extract NCAN context, 24G
awk -v OFS="\t" '{if ($1 ~ "...CA..") print substr($1,3,4)"\t"$2"\t"$3"\t"$4"\t"$5}' mouse_tissue.nonCpG.chr19.context.txt > mouse_tissue.ncan.chr19.context.txt &

# extract CAN context, 23G
awk -v OFS="\t" '{if ($1 ~ ".CA.") print substr($1,2,3)"\t"$2"\t"$3"\t"$4"\t"$5}' mouse_tissue.ncan.chr19.context.txt > mouse_tissue.can.chr19.context.txt &


wc -l mouse_tissue.*
#    628983453 mouse_tissue.CHG.chr19.context.txt
#   2111756248 mouse_tissue.CHH.chr19.context.txt
#    121951594 mouse_tissue.CpG.chr19.context.txt
#    319594031 mouse_tissue.methylated.chr19.context.txt
##   11.2% are methylated
#   2740739701 mouse_tissue.nonCpG.chr19.context.txt
#   1046548152 mouse_tissue.ncan.chr19.context.txt

```

## load datasets in R

```r
R
library(data.table)
library(ggplot2)

# Set width
options(width = 300)

# Load mfg data
# data <- fread("~/data/mouse_tissue.nonCpG.chr19.context.txt")
# # nrow larger than current 2^31 limit
# # See ?"Memory-limits" for more details.

data <- fread("~/data/mouse_tissue.ncan.chr19.context.txt")

setnames(data, c("context", "pct_met", "mouse_tissue", "age", "rep"))

table(data$mouse_tissue)
# EmbryonicFacialProminence                 Forebrain                     Heart                 Hindbrain                 Intestine                    Kidney                      Limb                     Liver                      Lung                  Midbrain                NeuralTube
#                  88650271                 116612512                 116017875                 114687386                  59614115                  58367569                  86677707                 102851963                  57196468                 115279796                  70454952
#                   Stomach
#                  60137538

table(data$age)
#     E10.5     E11.5     E12.5     E13.5     E14.5     E15.5     E16.5        P0
#  91116854 117884145 128076443  97822849 171322970 175816157 134221615 130287119
```


## extract CAN methylation

```r
data_cay <- data

# Extract cay context
data_cay[, context_cay := as.vector(sapply(data_cay$context, function(x) paste(unlist(strsplit(x, ""))[2:4], collapse = "")))]

# remove context with N
# The like operator is a simple wrapper for grep(..., value=TRUE)
data_cay <- data_cay[!context_cay %like% "N"]

# Explore cay context
## Forebrain
data_cay[mouse_tissue == "Forebrain" & age == "E10.5", .(.N, pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
#    context_cay       N pct_met_mean
# 1:         CAC 3150820         0.56
# 2:         CAG 4445654         0.42
# 3:         CAT 3746466         0.38
# 4:         CAA 3832375         0.36

data_cay[mouse_tissue == "Forebrain" & age == "P0", .(.N, pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay)][order(-pct_met_mean)]
#    context_cay       N pct_met_mean
# 1:         CAC 3294836         1.07
# 2:         CAG 4394884         0.72
# 3:         CAT 3778649         0.67
# 4:         CAA 3836885         0.61

data_cay[mouse_tissue == "Forebrain", .(.N, pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cay, age)][order(-pct_met_mean)]
#     context_cay   age       N pct_met_mean
#  1:         CAC    P0 3294836         1.07
#  2:         CAC E14.5 3065323         0.78
#  3:         CAC E12.5 5752464         0.74
#  4:         CAG    P0 4394884         0.72
#  5:         CAC E15.5 2989962         0.72
#  6:         CAC E11.5 3081536         0.71
#  7:         CAT    P0 3778649         0.67
#  8:         CAA    P0 3836885         0.61
#  9:         CAC E16.5 3004688         0.60
# 10:         CAG E14.5 4187631         0.57
# 11:         CAC E10.5 3150820         0.56
# 12:         CAG E12.5 8224290         0.55
# 13:         CAG E11.5 4194668         0.54
# 14:         CAT E14.5 3615956         0.54
# 15:         CAG E15.5 4168500         0.52
# 16:         CAT E12.5 7013917         0.51
# 17:         CAA E14.5 3699725         0.51
# 18:         CAT E11.5 3612488         0.50
# 19:         CAA E12.5 7248591         0.49
# 20:         CAT E15.5 3584545         0.48
# 21:         CAA E11.5 3695172         0.47
# 22:         CAA E15.5 3679847         0.45
# 23:         CAG E10.5 4445654         0.42
# 24:         CAG E16.5 4105870         0.42
# 25:         CAT E16.5 3554377         0.39
# 26:         CAT E10.5 3746466         0.38
# 27:         CAA E10.5 3832375         0.36
# 28:         CAA E16.5 3652393         0.36

# save data_cay table
fwrite(data_cay[, .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(mouse_tissue, context_cay, age, rep)][order(-pct_met_mean)], file = "~/mouse_tissue_chr19_CAN.txt", sep = "\t", row.names=FALSE, quote=FALSE)
```

### plot CAC/CAG, CAY_CAR

```r
cac_cag <- dcast(data_cay[, .(pct_met = mean(pct_met, na.rm=TRUE)), by = .(mouse_tissue, context_cay)], mouse_tissue ~ context_cay, value.var = "pct_met")
cac_cag[, `:=`(CAC_CAG = CAC/CAG, CAY_CAR = (CAC + CAT)/(CAG + CAA))]
cac_cag[order(CAC_CAG)]
#                  mouse_tissue       CAA       CAC       CAG       CAT  CAC_CAG  CAY_CAR
#  1:                     Liver 0.4836600 0.7454694 0.5406185 0.5091728 1.378920 1.224903
#  2:                 Forebrain 0.4676460 0.7421644 0.5357661 0.4982549 1.385240 1.236201
#  3:                      Limb 0.4766360 0.7566539 0.5423338 0.5051103 1.395181 1.238274
#  4: EmbryonicFacialProminence 0.4686779 0.7550125 0.5363783 0.4974769 1.407612 1.246189
#  5:                NeuralTube 0.4741638 0.7672549 0.5430646 0.5044901 1.412824 1.250206
#  6:                  Midbrain 0.4880899 0.8211095 0.5607379 0.5234554 1.464338 1.281969
#  7:                      Lung 0.4959886 0.8295126 0.5644742 0.5293456 1.469532 1.281382
#  8:                    Kidney 0.4572306 0.7864982 0.5265964 0.4898921 1.493550 1.297373
#  9:                 Hindbrain 0.4947994 0.8617753 0.5746118 0.5343253 1.499752 1.305485
# 10:                 Intestine 0.4732612 0.8397022 0.5489881 0.5092102 1.529545 1.319553
# 11:                     Heart 0.4688009 0.8723488 0.5529133 0.5080850 1.577731 1.351096
# 12:                   Stomach 0.4900126 0.9078036 0.5712848 0.5317115 1.589056 1.356373

## plot CAC_CAG
cac_cag$mouse_tissue <- factor(cac_cag$mouse_tissue, level = cac_cag$mouse_tissue[order(cac_cag$CAC_CAG)])

gg <- ggplot(cac_cag, aes(x=mouse_tissue, y=CAC_CAG)) +
geom_bar(stat="identity", color="black", position=position_dodge(), width= 0.6, alpha = 0.5) +
theme_bw() +
ylab(expression("[mCAC/CAC]/[mCAG/CAG]")) +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 12, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
coord_cartesian(ylim = c(0, 1.75))
ggsave("~/mouse_tissue_chr19_CAC_CAG.pdf")

## plot CAY_CAR
cac_cag$mouse_tissue <- factor(cac_cag$mouse_tissue, level = cac_cag$mouse_tissue[order(cac_cag$CAC_CAG)])

gg <- ggplot(cac_cag, aes(x=mouse_tissue, y=CAY_CAR)) +
geom_bar(stat="identity", color="black", position=position_dodge(), width= 0.6, alpha = 0.5) +
theme_bw() +
ylab(expression("[mCAC/CAC+mCAT/CAT]/[mCAG/CAG+mCAA/CAA]")) +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 12, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
coord_cartesian(ylim = c(0, 1.75))
ggsave("~/mouse_tissue_chr19_CAY_CAR.pdf")


## plot CAC_CAG, CAY_CAR side by side
cac_cag_melt <- melt(cac_cag[, c(1,6,7)], value.var = "mouse_tissue")

cac_cag_melt$mouse_tissue <- factor(cac_cag_melt$mouse_tissue, level = cac_cag$mouse_tissue[order(cac_cag$CAC_CAG)])

gg <- ggplot(cac_cag_melt, aes(x = mouse_tissue, y = value, fill = variable)) +
geom_bar(stat="identity", color="black", position=position_dodge(), width= 0.6, alpha = 0.5) +
theme_bw() +
ylab(expression("[mCAC/CAC]/[mCAG/CAG]")) +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 12, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
coord_cartesian(ylim = c(0, 1.75))
ggsave("~/mouse_tissue_chr19_CAC_CAG_CAR_CAR.pdf")

rm(cac_cag_melt)
rm(cac_cag)
```

### plot CAN changes through age

```r
can_age <- data_cay[, .(pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(mouse_tissue, context_cay, age, rep)]

# save can_age table
fwrite(can_age, file = "~/mouse_tissue_chr19_CAN_age.txt", sep = "\t", row.names=FALSE, quote=FALSE)

## load data
# R
# library(data.table)
# library(ggplot2)
# options(width = 300)
# can_age <- fread("~/mouse_tissue_chr19_CAN_age.txt")

# caculate mean,sd in replicates
library(plyr)

can_age_rep <- ddply(can_age, c("mouse_tissue", "context_cay", "age"), summarise, rep_N = length(rep), rep_mean = mean(pct_met_mean),rep_sd = sd(pct_met_mean), rep_se = rep_sd/sqrt(rep_N))
# after ddply,  can_age_rep is just a data frame, while can_age is data table and data frame, setDT could turn  can_age_rep into data table. different between DT & DF: https://stackoverflow.com/questions/13618488/what-you-can-do-with-a-data-frame-that-you-cant-with-a-data-table
setDT(can_age_rep)

# set the order of libraries
can_age_rep$age <- factor( can_age_rep$age, levels = c("E10.5", "E11.5", "E12.5", "E13.5", "E14.5","E15.5", "E16.5", "P0"))

can_age_rep$mouse_tissue <- factor( can_age_rep$mouse_tissue, levels = c("Forebrain", "Midbrain", "Hindbrain", "NeuralTube", "Heart","Liver", "EmbryonicFacialProminence", "Limb", "Kidney",  "Intestine", "Stomach", "Lung"))

can_age_rep$context_cay <- factor(can_age_rep$context_cay, levels = can_age_rep[, .(.N, rep_mean_mean = round(mean(rep_mean, na.rm=TRUE), 2)), by = .(context_cay)][order(-rep_mean_mean)]$context_cay)
## levels = "CAC" "CAG" "CAT" "CAA"

# plot
gg <- ggplot( can_age_rep, aes(x = age, y = rep_mean, color = context_cay, group = context_cay)) +
geom_point() +
geom_line() +
geom_errorbar(aes(ymin = rep_mean - rep_sd, ymax = rep_mean + rep_sd), width = 0.2) +
# scale_linetype_manual(values = c(1,1,1,1,1)) +
theme_bw() +
facet_wrap(~ mouse_tissue, ncol = 4) +
ylab(expression("CA methylation (%)")) + 
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=12), axis.text.y = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, size = 10, color = "black", hjust = 1), legend.text = element_text(size = 12, color = "black")) +
coord_cartesian(ylim = c(0, 1.5))
ggsave("~/mouse_tissue_chr19_CAN_age.pdf", width = 25, height = 20, units = "cm")

# plot brain (add colour)
gg <- ggplot( can_age_rep[mouse_tissue %in% c("Forebrain", "Midbrain", "Hindbrain")], aes(x = age, y = rep_mean, color = context_cay, group = context_cay)) +
geom_point() +
geom_line() +
## levels = "CAC" "CAG" "CAT" "CAA"
scale_colour_manual(values = c('blue', 'orange', 'red', 'darkgreen')) +
geom_errorbar(aes(ymin = rep_mean - rep_sd, ymax = rep_mean + rep_sd), width = 0.2) +
# scale_linetype_manual(values = c(1,1,1,1,1)) +
theme_bw() +
facet_wrap(~ mouse_tissue, ncol = 3) +
ylab(expression("CA methylation (%)")) + 
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=12), axis.text.y = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, size = 10, color = "black", hjust = 1), legend.text = element_text(size = 12), legend.position="top") +
coord_cartesian(ylim = c(0, 1.5))
ggsave("~/mouse_tissue_chr19_CAN_age.pdf", width = 25, height = 20, units = "cm")

# plot others (add colour)
gg <- ggplot( can_age_rep[mouse_tissue %in% c("NeuralTube", "Heart","Liver", "EmbryonicFacialProminence", "Limb", "Kidney",  "Intestine", "Stomach", "Lung")], aes(x = age, y = rep_mean, color = context_cay, group = context_cay)) +
geom_point() +
geom_line() +
## levels = "CAC" "CAG" "CAT" "CAA"
scale_colour_manual(values = c('blue', 'orange', 'red', 'darkgreen')) +
geom_errorbar(aes(ymin = rep_mean - rep_sd, ymax = rep_mean + rep_sd), width = 0.2) +
# scale_linetype_manual(values = c(1,1,1,1,1)) +
theme_bw() +
facet_wrap(~ mouse_tissue, ncol = 3) +
ylab(expression("mCAN/CAN")) + 
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=12), axis.text.y = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, size = 10, color = "black", hjust = 1), legend.text = element_text(size = 12), legend.position="top") +
coord_cartesian(ylim = c(0, 1.5))
ggsave("~/mouse_tissue_chr19_CAN_age.pdf", width = 25, height = 20, units = "cm")


########### plot CAC/CAG, CAY_CAR ratios through age ###########
can_age <- data_cay[, .(pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(mouse_tissue, context_cay, age, rep)]

cac_cag_age <- dcast(can_age, mouse_tissue + age + rep ~ context_cay, value.var = "pct_met_mean")
cac_cag_age[, `:=`(CAC_CAG = CAC/CAG, CAY_CAR = (CAC + CAT)/(CAG + CAA))]

# convert wide to long format
cac_cag_age_melt <- melt(cac_cag_age[, c(1:3,8,9)], id.vars = c("mouse_tissue", "age", "rep"))

# caculate mean,sd in replicates for CAC_CAG, CAY_CAR
library(plyr)

cac_cag_age_melt_rep <- ddply(cac_cag_age_melt, c("mouse_tissue", "age", "variable"), summarise, rep_N = length(rep), rep_mean = mean(value), rep_sd = sd(value), rep_se = rep_sd/sqrt(rep_N))

setDT(cac_cag_age_melt_rep)

# set the order of libraries
cac_cag_age_melt_rep$age <- factor(cac_cag_age_melt_rep$age, levels = c("E10.5", "E11.5", "E12.5", "E13.5", "E14.5","E15.5", "E16.5", "P0"))

cac_cag_age_melt_rep$mouse_tissue <- factor(cac_cag_age_melt_rep$mouse_tissue, levels = c("Forebrain", "Midbrain", "Hindbrain", "NeuralTube", "Heart","Liver", "EmbryonicFacialProminence", "Limb", "Kidney",  "Intestine", "Stomach", "Lung"))

# plot CAC_CAG, CAY_CAR ratios through age
gg <- ggplot(cac_cag_age_melt_rep, aes(x = age, y = rep_mean, color = variable, group = variable)) +
geom_point() +
geom_line() +
geom_errorbar(aes(ymin = rep_mean - rep_sd, ymax = rep_mean + rep_sd), width = 0.2) +
# scale_linetype_manual(values = c(1,1,1,1,1)) +
theme_bw() +
facet_wrap(~ mouse_tissue, ncol = 4) +
ylab(expression("[mCAC/CAC]/[mCAG/CAG]")) + 
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=12), axis.text.y = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, size = 10, color = "black", hjust = 1), legend.text = element_text(size = 12, color = "black")) +
coord_cartesian(ylim = c(0, 2))
ggsave("~/20190722_mouse_tissue_chr19_CAC_CAG_CAY_CAR_age.pdf", width = 25, height = 20, units = "cm")

rm(cac_cag_age_melt_rep)
rm(cac_cag_age_melt)
rm(can_age)
rm(can_age_rep)
rm(data_cay)
rm(data)

```



## save and delete files

```sh
cd ~/data

# remove bedmethyl & bed files
rm *bed*
rm *ENCFF*.context.txt

# compress context files
pigz *.context.txt &
```


# extract CGN methylation


```r
cd ~/data
unpigz mouse_tissue.CpG.chr19.context.txt

R
library(data.table)
library(ggplot2)

# Set width
options(width = 300)

data <- fread("~/data/mouse_tissue.CpG.chr19.context.txt")

setnames(data, c("chr", "start", "end", "pct_met", "cnt_tot", "strand", "context", "mouse_tissue", "age", "rep"))

table(data$mouse_tissue)


data_cgy <- data

# Extract cgy context
data_cgy[, context_cgy := as.vector(sapply(data_cgy$context, function(x) paste(unlist(strsplit(x, ""))[4:6], collapse = "")))]

# remove context with N
# The like operator is a simple wrapper for grep(..., value=TRUE)
data_cgy <- data_cgy[!context_cgy %like% "N"]

# Explore cgy context
## Forebrain
data_cgy[mouse_tissue == "Forebrain" & age == "E10.5", .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cgy)][order(-pct_met_mean)]
#    context_cgy      N pct_met_median pct_met_mean
# 1:         CGA 446322             89        81.40
# 2:         CGT 522594             88        80.97
# 3:         CGC 376848             91        79.86
# 4:         CGG 483082             90        79.76

data_cgy[mouse_tissue == "Forebrain" & age == "P0", .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(context_cgy)][order(-pct_met_mean)]
#    context_cgy      N pct_met_median pct_met_mean
# 1:         CGC 363881             91        79.83
# 2:         CGA 436008             87        77.68
# 3:         CGG 464827             89        77.50
# 4:         CGT 514761             86        76.90

data_cgy[mouse_tissue == "Forebrain", .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE))), by = .(context_cgy, age)][order(age)]
#     context_cgy   age      N pct_met_median pct_met_mean
#  1:         CGA E10.5 446322             89           81
#  2:         CGT E10.5 522594             88           81
#  3:         CGG E10.5 483082             90           80
#  4:         CGC E10.5 376848             91           80
#  5:         CGA E11.5 415358             88           81
#  6:         CGT E11.5 490255             87           80
#  7:         CGG E11.5 439509             90           80
#  8:         CGC E11.5 340322             91           81
#  9:         CGA E12.5 819991             88           80
# 10:         CGT E12.5 958346             87           79
# 11:         CGG E12.5 860216             90           80
# 12:         CGC E12.5 658677             91           81
# 13:         CGA E14.5 414794             88           79
# 14:         CGT E14.5 489134             87           78
# 15:         CGG E14.5 436625             90           79
# 16:         CGC E14.5 337312             92           81
# 17:         CGA E15.5 414809             88           79
# 18:         CGT E15.5 487451             86           78
# 19:         CGG E15.5 437485             90           79
# 20:         CGC E15.5 337224             91           80
# 21:         CGA E16.5 405271             88           79
# 22:         CGT E16.5 476146             86           78
# 23:         CGG E16.5 423609             90           78
# 24:         CGC E16.5 327193             92           80
# 25:         CGA    P0 436008             87           78
# 26:         CGT    P0 514761             86           77
# 27:         CGG    P0 464827             89           77
# 28:         CGC    P0 363881             91           80

data_cgy[mouse_tissue == "Forebrain", .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE))), by = .(age)][order(age)]
#      age       N pct_met_median pct_met_mean
# 1: E10.5 1828846             90           81
# 2: E11.5 1685444             89           80
# 3: E12.5 3297230             89           80
# 4: E14.5 1677865             89           79
# 5: E15.5 1676969             88           79
# 6: E16.5 1632219             89           79
# 7:    P0 1779477             88           78


fb <- data_cgy[mouse_tissue == "Forebrain", .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE))), by = .(context_cgy, age, rep)][order(-pct_met_mean)]

dcast(fb, age ~ context_cgy + rep, value.car = "pct_met_mean")
#      age CGA_1 CGA_2 CGC_1 CGC_2 CGG_1 CGG_2 CGT_1 CGT_2
# 1: E10.5    81    82    80    80    80    80    81    81
# 2: E11.5    81    81    81    81    80    80    80    80
# 3: E12.5    80    80    81    81    80    80    79    79
# 4: E14.5    80    79    81    80    79    79    79    78
# 5: E15.5    79    79    80    80    78    79    78    78
# 6: E16.5    79    79    80    80    78    78    78    77
# 7:    P0    78    78    80    80    77    78    77    77

fwrite(data_cgy[, .(.N, pct_met_median = as.double(median(pct_met, na.rm=TRUE)), pct_met_mean = round(mean(pct_met, na.rm=TRUE), 2)), by = .(mouse_tissue, context_cgy, age, rep)][order(-pct_met_mean)], file = "~/20190809_mouse_tissue_chr19_CGN.txt", sep = "\t", row.names=FALSE, quote=FALSE)
```
