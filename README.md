# SplicingVariants_Beta
Developed and Maintained by Kaoru Ito (splicing.variant[at]gmail.com)

This document is for the paper "Mutations in LMNA and MYBPC3 That Alter RNA splicing" in PNAS.
In the paper, we chose possibile splicing variants by in-silico prediction and tested them by cell-based splicing assay.
Several perl and R scripts placed here were employed to perform the analysis.

*** In the paper, we focused on splcing variants that create / lose splice site. If you'd like to know about exon-skipping variants, please check SPANR (http://tools.genes.toronto.edu/)***

# Workflow
1) Calculate scores to assess the possibility for a splicing variant and choose candidates for cell-based splicing assay<br>
  -- Regress_Score.v0.**.R<br>
2) Design constructs for the candidates to perform cell-based splicing assay<br>
  -- ConstructDesigner.v0.**.R<br>
3) Perform cell-based splicing assay (transfect constructs to cells, RNA extraction from the cells, prepare libraries for Miseq run<br>
4) Count normal / aberrant splciing products and calculate p-values from raw-fastq files (alignment not required)<br>
  -- Make.inputframe.*.pl<br>
  -- SpliceConstructSearchGrepV1.*.pl<br>
*means developing version. Please use the latest one.

To run each script, please read README_regress.score.md, README_construct.designer.md and README_spliceconstructsearch.md.

# Rationale for the cell-based splicing assay (p-values obtained from SpliceConstructSearchGrepV1.5.pl)
  Because p-value of the assay reflects the degree of aberrant splicing which changes continuously, we investigated the relationship between the p-values and clinical diagnosis in variants from clinical databases (Supplementary Fig.).  As expected, we found a significant correlation between -log10 p-value of the assay and severity of the clinical diagnosis (missense variants included p=8.2e-06, missense variants excluded p=1.2e-06, kendall’s rank correlation test).  Since clinically likely-benign and benign variants never surpassed p<0.001 and just missense variants (which can be deleterious without affectingh splicing) in likely-pathogenic and pathogenic groups had p>0.001,  we defined the threshold p=0.001 to distinguish pathogenic variants from benign variants in terms of a splicing variant. 
 ![Supplementary Figure](https://github.com/SplicingVariant/SplicingVariants_Beta/blob/master/Supplementary%20Figure.JPG)

