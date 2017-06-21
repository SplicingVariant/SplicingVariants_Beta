README_regress.score.md
===========
Written by Kaoru Ito 

<br />

This document mentions about how to use Regress_Score.v0.97.R.  
Regress_Score.v0.97.R is a script to calculate MaxEnt scores to assess the possibility for a splicing variant.  Please note that the algorithm employed in this script considers just splice donor / acceptor site gain / loss.   


Although the older version utilized Biomart server provided by ENSEMBL, the latest version can run stand-alone, which means no internet access required. But you need to prepare several files to do so. (All files placed on the SplicingVariants_Beta folder) 

##What are required to run the script?
-----------
1. **R 3.3.1 or later (could run on older versions of R, but no warranty)**  
https://www.r-project.org/  
2. **Biocondutor package:BSgenome.Hsapiens.UCSC.hg19 when the genome coordinate is hg19 (GRCh37)**  
https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html  
3. **Bioconductor package:BSgenome.Hsapiens.UCSC.hg38 when the genome coordinate is hg38 (GRCh38)**  
https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html  
4. **SS files**  
Pease make a folder named "ssfiles" and place files from ssfiles.zip (provided here)  
These files are required to calculate scores.  
5. **usehg19 files when the genome coordinate is hg19 (GRCh37)**   
Pease make a folder named "usehg19" and place files from usehg19.tar.gz (provided here)
These files contain hg19 gene and transcript information.  
6. **usehg38 files when when the genome coordinate is hg38 (GRCh38)**  
Pease make a folder named "usehg38" and place files from usehg38.tar.gz (provided here)
These files contain hg38 gene and transcript information.

##File format that contains variant information
----------
1. **mutfile**   
A unique format for this script. Use -mutfile to provide the mutfile. The file format is as follows.  

```js:
chr1:156105901C>T ENST00000368300  
chr1:156105564G>A ENST00000368300
```

One variant per line.  
chr1:156105564G>A (Variant name) and ENST00000368300 (ENSEMBL transcriptID) should be separated by a "space".

If you don't provide ENSEMBL transcriptID, this script automatically finds all transcripts that cover a given variant. So you can provide variant information like   

```js:
chr1:156105901C>T  
chr1:156105564G>A
```

In this case, the procedure to search transcripts is sometimes time-consuming. To avoid the timeloss, you can indicate one transcriptID for all variants using -esbTranscriptID, or just take a canonical transcript (longest protein-coding transcript for a gene) using -takeCanonicalTranscript.

Also you can provide variant information using cDNA position like

```js:
c.640-10A>G ENST00000368300
c.1146C>T ENST00000368300
c.513G>A ENST00000368300
```

Again, c.1146C>T and ENST00000368300 should be separated by a "SPACE". 
 
As long as you present the transcriptID using -esbTranscriptID, the mut file can be without transcriptID information like below.

```js:
c.1146C>T
c.513G>A
```

2. **vcf file**  
A vcf file can be provided using -vcffile  
This script searches all transcripts included in the INFO section.   
In case of no transcript information provided  in the vcf file, this script automatically searches all transcripts that cover a given variant. To avoid time-consuming for the transcript search, you can make use of  -esbTranscriptID or -takeCanonicalTranscript.
  
##Basic Usage and Examples
-------------
Please make sure that all required software were installed and that ssfile, usehg19 or usehg38 folders were prepared.

Just for a mutfile analysis, please type a command on command line interface (bash on UNIX, terminal on Mac or Cygwin on Windows)

```js:
R --slave --vanilla --args -mutfile [mutfile] -output [raw output file] -summarizeresult [basename for summarized result] -sjdbout [output filename for sjdb out] -skipRegressScore -skipSRE -refdirectory []folder that contains ssfiles e.g. /Volumes/ssfiles ] -usehg19 [folder that contains usehg19files e.g. /Volumes/usehg19]  -skipcDNApos < Regress_Score.v0.97.R
```

The above command generates 3 result files.  

1.  **[raw output file]**  
 This file contains all calculated scores. Since a score calculation bin (9bp-long for splice donor site, 23bp-long for splice acceptor site) slides around a mutation, the mutation has several score results.
2. **[basename for summarized result].SigLines.txt**  
 This file contains most significant lines for a variant.  
i.e. 1 line for aberrant splice donor site gain, 1 line for aberrant splice acceptor site gain, 1 line for canonical spice donor site loss (if exists), 1line for canonical splice acceptor site loss (if exists).  
Additionally, "Kind","Possibility" and "Decision" columns provided for decision-making.
Using these parameters, you can select candidates for the cell-based splicing assay.

3. **[output filename for sjdb out]**  
 This file contains splice junction data for STAR (RNAseq aligner) sjdb file.
 When you calculate scores and confirm aberrant splcing using RNA-seq, please provide this sjdb file to STAR when alignment.


##Flag index for Regress_Score.v.0.97.R
----------------------
**-mutfile** (chrpos format mut file with transcript ID) e.g.chr1:123G>A ENST000001234 in v0.5 just 1 line accept  
**-seq** (a file contains sequence data)  exon should be UPPERCASE, intron should be LOWERCASE, you can indicate mutation like (A|B), where A means refseq and B mean alt seq
                                       Also you can indicate nucleutide of interest like ATGC[ATGC]ATGCa              
**-vcffile** (vcf file) indicates a variant information file as a input 
     *If TRANSCIPT_ID=ENSTXXXXXXXXX, pick the ensembl ID for transcript ID  
**-output** (filename to be saved as a result)  
**-marthost** (indicates biomaRt server)  
**-esbTranscriptID** (transcriptID)  when you check cDNA variant list, you should provide transcriptID  
**-nomart**        if you indicate this switch, connecting biomaRt is being skipped.  
**-canonicals**  indicate this script to process canonical 5ss and 3ss around the given mutation, if the mutation is not on the exon, this function does not work.  
**-entiregene** indicates this script wull calculate entire gene score   
**-ss5** indicates this script to calculate ss5 scores  
**-ss3** indicates this script to calculate ss3 scores    
     if omitted the both above automatically calculate both
**-wide** (integer) indicates to make scan region added (integer) around the mutation  
**-lseq** (filename) indicagtes transcript info generated by GetSeq(TranscriptID,1) if you employ this, you don't need to connect biomart server , which make run time short.  
                 whme you use it, I recommned you to use -esbTranscriptID to indicate the name of transcript.  
**-skipRegressScore** indicates skipping regress score calculation to make time short. The score will be zero  
   The regress score is an integrated score of MaxEnt and splicing regulatory elements, which is better than just MaxEnt score or just splicing regulatory elements.  
**-skipSRE** indicates skipping checking # of splicing regulatory elements. The # of elements will be zero.  
**-localseq** (foldername) is a siwtch to reduce the access to biomaRt server. In concrete,  
                      1) if required seq is on the indicated local folder, employ it  
                      3) if not, retrieve the seq from biomaRt server and save the seq on the indicated local server.  
**-refdirectory** (directory where maxent and other files are located) :when these files are no loacted on the current directory , please indicate their locations  
**-useGTF19** (directory that contains Homo_Sapiens.GRCh37.75,gtf for ensemble hg19 gtf file downloaded extarct from http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz)  
**-usehg19** (directory that contains hg19.transcripts.txt and hg19.exons.txt) those 2 files are generated from Homo_Sapiens.GRCh37.75,gtf by Break_Ensembl_GTF.pl  
**-useGRCh38**  (directory that contains GRCh38.transcripts.txt and GRCh38.exons.txt) those 2 files are generated from Homo_Sapiens.GRCh38.87,gtf by Break_EnsemblGRCh38_GTF.pl  
**-summarizeresult** (prefix of output file for summarizeresult) one variant 9 scores for 5SS and 23scores for 3SS. Then just pick siginificant lines file *.SigLines.txt will be generated   
**-summarizeresultin** (raw result file of Regress_Score) when you want just to summarize the result, please indicate the result file, which mean skipping calculations.  
**-SummarizeByVariant** :when you indicate -summarizresult, this switch make this script generate files. *.SigLines.txt and *.SigLinesByVariant.txt  
**-TakeOnePerVariant** indicates that the script dumps no-significant lines and integrates significant lines to make one line for one variant  
**-skip1stMiningStep** indicates skipping 1st procedure when you indicate -summarizresult   
**-skip2ndMiningStep** indicates skipping 1st procedure when you indicate -summarizresult  
**-sjdbout** (sjdbfilename for STAR) this switch indicates to make sjdbFile for STAR RNAseq   aligner --sjdbFileChrStartEnd. also result file will have Predicted_SpliceJunction column  
     *when use this switch, you also need to indicate -summarizeresult. That's because this requires the column <Kind>, which is attached by -summarizeresult  
**-sjdbrestrictin** (columnname, value) gives an inclusion creteria when -sjdb being indicated. just lines with the column = value will be included  
**-sjdbrestrictout** (columnname, value) gives an exclusion creteria when -sjdb being indicated. just lines with the column = value will be excluded  
**-sjdbin** (result file of Regress_Score) if you want to skip score calculation and give the result file directory, use this switch.  
**-takeCanonicalTranscript** indicates when multiple transcripts found, take a canonical transcript 
 #the definition of Canonical Transcript means   
  http://asia.ensembl.org/Help/Glossary?id=346
  For human, the canonical transcript for a gene is set according to the following hierarchy: 1. Longest CCDS translation with no stop codons. 2. If no (1), choose the longest Ensembl/Havana merged translation with no stop codons. 3. If no (2), choose the longest translation with no stop codons. 4. If no translation, choose the longest non-protein-coding transcript.  
  *CAUTION this switch currently works when you employ -usehg19  
**-skipcDNApos** indicates skipping providing variant name in cDNA positon format (e.g. c.123A>G)when the original variant name is given in chromosome position format (e.g. chr1:12345A>G)   



