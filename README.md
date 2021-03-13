# Colorectal-Cancer Early onset 
---
* [Introduction](#Introduction)
* [Resource](#Resource)
* [Pipeline](#Pipeline)

<a name="Introduction"/>

# Introduction

We analyzed somatic mutation and clincial data from the GENIE project (Genomics Evidence Neoplasia Information Exchange, Release 9.0-public
To cite this release: doi.org/10.7303/syn24179657), with the gole of to investigate "Racial differences in somatic mutations among patients with early-onset colorectal cancer". This study was led by Drs. Andreana N. Holowatyj and Xingyi Guo, as well as contributed by colleagues at Vanderbilt, Meharry Medical College, Harvard Medical School, and Fundación Jiménez Díaz University Hospital.  The dbGap approaved project (PI: Xingyi Guo): ID 24541 - titled	"Characterizing driver somatic mutations and genes of early onset cancer types in a pan-cancer analysis".

Andreana N. Holowatyj, PhD, MS1,3*; Wanqing Wen, MD1; Timothy Gibbs, MS4; 
Kay M. Washington, MD, PhD2,3; Cathy Eng, MD1,3; Paulette D. Chandler, MD, MPH5; Jose Perea, MD, PhD6; Wei Zheng, MD, PhD1,3; Xingyi Guo, PhD1,3*

Author Affiliations:
1Department of Medicine, 2Department of Pathology, Microbiology and Immunology, Vanderbilt University Medical Center, Nashville, TN
3Vanderbilt-Ingram Cancer Center, Nashville, TN
4Meharry Medical College, Nashville, TN
5Division of Preventive Medicine, Brigham and Women’s Hospital, Harvard Medical School, Boston, MA
6Department of Surgery, Fundación Jiménez Díaz University Hospital, Madrid, Spain
* contact authors: Drs. Holowatyj and Guo.


# Resource

### R1. Study cohorts of Colorectal cancer patients, including demographic characteristics: Age at sequencing report, Sex, histological subtype, sequencing center/assays, and detailed tumor sample descrption.  

### R2. Somatic mutation data with functional annotation including: 
 ```
3'UTR
5'UTR
Frame_Shift_Del
Frame_Shift_Ins
In_Frame_Del
In_Frame_Ins
Intron
Missense_Mutation
Nonsense_Mutation
Nonstop_Mutation
Nonstop_Mutation
RNA
Silent
Splice_Region
Splice_Site
Targeted_Region
Translation_Start_Site
Nosilent mutations:
mutation<-c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Del","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Site","Translation_Start_Site")

```


<a name="Pipeline"/>

# Pipeline 
---


## step1: build a matrix of somatic mutation and clinical data 

```
mut.dat<-fread("PATH/GENIE_9/data_mutations_extended.txt",head=TRUE,sep = "\t")
cliP<- read.table("PATH/GENIE_9/data_clinical_patient.txt",head=TRUE,sep="\t")
cliS<- fread("PATH/GENIE_9/data_clinical_sample.txt",sep="\t",head=TRUE)
...
```

## step2: load gene assay and coverage information 
1) Need to adjust tumor mutational burden (TMB) and infer MSI/MSS status based on a cutoff of TBM. 
2) To perform additional QC if assays cover less than 50K bps. 
3) We additional examined the boxplot of TMB by assays to remove outliners.
```
PANEL<-read.table("PATH/GENIE_9/genomic_information_noFalse.txt",head=TRUE,sep="\t")
TMB_ref<-read.table("PATH../CRC../GENIE_PANEL.COV",head=TRUE)
```

more information see main.R 
