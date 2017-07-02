# SplicingVariants_Beta
Developed and Maintained by Kaoru Ito (splicing.variant[at]gmail.com)

This document is for the paper **"Identification of Pathogenic Mutations in *LMNA* and *MYBPC3* That Alter RNA splicing" (Kaoru Ito and Parth Patel et al.) in PNAS 2017**.<br>
In the paper, we chose possibile splicing variants by *in-silico* prediction and tested them by cell-based splicing assay.
Several perl and R scripts placed here were employed to perform the analysis.

*** In the paper, we focused on splcing variants that create / lose splice site. If you'd like to know about exon-skipping variants, please check SPANR (http://tools.genes.toronto.edu/)***

# What is "Splicing Variant"(Splice-Altering Variant)?
**In the intepretation of genetic variants, synonymous mutation is considered benign. Also missense mutaiton is on the fence in terms of its deleteriousness, whereas stop-gain/loss, framn-shift and splice-site (GT/AG-broken) mutations are considered damaging.
These intepretations are NOT ALWAYS TRUE.**<br>
 ![Supplementary Figure](https://github.com/SplicingVariant/SplicingVariants_Beta/blob/master/WhichMutationIsPathogenic.png)
 
Messenger RNA splicing, where intron regions flanked by a splice donor site and a splice acceptor site are cut out, occurs before protein syhthesis by mRNA translation.<br>
In the mRNA splicing procedure, it's known that a nucleotide alteration that doens't change amino acide sequence can activate a cryptic splice site, which results in creating an aberrant splice donor/acceptor site in the middle of an exon / intron, followed by the exon truncation / intron retention in the protein.
Additionally, not only GT/AG broken but also nearby nucleotide changes (donor:-3~+6bp, acceptor:-20~+3 from the splice juction) can disrupt the function of a canonical splice site.<br>
 ![Supplementary Figure](https://github.com/SplicingVariant/SplicingVariants_Beta/blob/master/SplicingVariant.png)
 
Several papers reported that such variants cause severe mendelian disroders, such as progeria syndrome and dilated cardiomyopathy 
However, because a method for detecting such splicing-altering variants in a high-throughput pipeline had not been developed, these variants were overlooled in usual NGS settings.<br>
Here we present a high-throughput-friendly method to detecte candidates for a splicing-altering variants.
Additionally we developed a cell-based splicing assay utilizing NGS technology for mulitiplexing analysis, where construct design and assessment of splicing alteration are automated.
  
# Workflow
1) Calculate scores to assess the possibility for a splicing variant and choose candidates for cell-based splicing assay<br>
  -- Regress_Score.v0.##.R<br>
2) Design constructs for the candidates to perform cell-based splicing assay<br>
  -- ConstructDesigner.v0.##.R<br>
3) Perform cell-based splicing assay (transfect constructs to cells, RNA extraction from the cells, prepare libraries for Miseq run<br>
4) Count normal / aberrant splciing products and calculate p-values from raw-fastq files (alignment not required)<br>
  -- Make.inputframe.#.pl<br>
  -- SpliceConstructSearchGrepV1.##.pl<br>
\* # means developing version. Please use the latest one.

To run each script, please read README\_regress.score.md, README\_construct.designer.md and README\_spliceconstructsearch.md.

# Rationale for the cell-based splicing assay (p-values obtained from SpliceConstructSearchGrepV1.5.pl)
  Because p-value of the assay reflects the degree of aberrant splicing which changes continuously, we investigated the relationship between the p-values and clinical diagnosis in variants from clinical databases (Supplementary Fig.).  As expected, we found a significant correlation between -log10 p-value of the assay and severity of the clinical diagnosis (missense variants included p=8.2e-06, missense variants excluded p=1.2e-06, kendall’s rank correlation test).  Since clinically likely-benign and benign variants never surpassed p<0.001 and just missense variants (which can be deleterious without affectingh splicing) in likely-pathogenic and pathogenic groups had p>0.001,  we defined the threshold p=0.001 to distinguish pathogenic variants from benign variants in terms of a splicing variant. 
 ![Supplementary Figure](https://github.com/SplicingVariant/SplicingVariants_Beta/blob/master/Supplementary%20Figure.JPG)

