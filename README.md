# SplicingVariants_Beta
Developed and Maintained by Kaoru Ito (splicing.variant@gmail.com)

This document is for the paper "Identification of Pathigenic Gene Mutations in LMNA and MYBPC3 That Alter RNA splicing".
In the paper, we chose possibile splicing variants by in-silico prediction and tested them by cell-based splicing assay.
Several perl and R scripts placed here were employed to perform the analysis.

*** In the paper, we focused on splcing variants that create / lose splice site. If you'd like to know about exon-skipping variants, please check SPANR (http://tools.genes.toronto.edu/)***

#Big Picture
1) Calculate the score and choose possible candidates for splicing variants
  Regress_Score.v0.88.R
2) Design constructs for the candidates to perform cell-based splicing assay
  ConstructDesigner.v0.93.R
3) Perform cell-based splicing assay (transfect constructs to cells, RNA extraction from the cells, prepare libraries for Miseq run
4) Count normal / aberrant splciing products and calculate p-values from raw-fastq files (alignment not required)
  Make.inputframe.2.pl
  SpliceConstructSearchGrepV1.5.pl

#Rationale for the cell-based splicing assay (p-values obtained from SpliceConstructSearchGrepV1.5.pl)
 Because p-value of the assay reflects the degree of aberrant splicing which changes continuously, we investigated the relationship between the p-values and clinical diagnosis in variants from clinical databases (Supplementary Fig. 3).  As expected, we found a significant correlation between -log10 p-value of the assay and clinical diagnosis (missense variants included p=8.2e-06, missense variants excluded p=1.2e-06, kendall’s rank correlation test).  Since clinically likely-benign and benign variants never surpassed p<0.001 and just missense variants in likely-pathogenic and pathogenic groups had p>0.001,  we defined the threshold p=0.001 to distinguish pathogenic variants from benign variants. 
 
