# SplicingVariants_Beta
Developed and Maintained by Kaoru Ito (splicing.variant@gmail.com)

This document is for the paper "Identification of Pathigenic Gene Mutations in LMNA and MYBPC3 That Alter RNA splicing".
In the paper, we chose possibile splicing variants by in-silico prediction and tested them by cell-based splicing assay.
Several perl and R scripts placed here were employed to perform the analysis.

*** In the paper, we focused on splcing variants that create / lose splice site. If you'd like to know about exon-skipping variants, please check SPANR (http://tools.genes.toronto.edu/)***

#Big Picture
1) Calculate the score and choose possible candidates for splicing variants
2) Design constructs for the candidates to perform cell-based splicing assay
3) Perform cell-based splicing assay (transfect constructs to cells, RNA extraction from the cells, prepare libraries for Miseq run
4) Count normal / aberrant splciing products and calculate p-values from raw-fastq files (alignment not required)
