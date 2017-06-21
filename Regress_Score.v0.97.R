Helpdocument<-"
#this script is to score given seq using the newly invented integration method (please see the last part of .docx) 
#
#Usage: R --slave --vanilla --args  -mutfile (chrpos format mut file with transcript ID) -output (output file name) < Regress_Score.v0.5.R
#
#ex) R --slave --vanilla --args --mutfile temp.test.chrpos.txt -output temp.test.chrpos.output.txt <Regress_Score.v0.5.R
#ex) R --slave --vanilla --args -seq testseq.txt -output temp.test.chrpos.output.v0.6.txt -debug <Regress_Score.v0.7.R
#ex) R --slave --vanilla --args --mutfile LMNAmut_for_gBlcok.txt -output LMNAmut_for_gBlcok.txt.regres_score.v0.7.txt -esbTranscriptID ENST00000368300 <Regress_Score.v0.7.R
#ex) R --slave --vanilla --args -canonicalSS --mutfile LMNAmut_for_gBlcok.txt -output LMNAmut_for_gBlcok.txt.regres_score.v0.8.txt -esbTranscriptID ENST00000368300 <Regress_Score.v0.8.R
#ex) R --slave --vanilla --args -canonicalSS -seq LMNAc768GAseq.txt -output LMNAc768GA.regres_score.v0.8.txt -debug <Regress_Score.v0.8.R
#ex) R --slave --vanilla --args --mutfile JamesTTNvariant.txt -output JamesTTNvariant.regres_score.v0.83.txt -entiregene -ss3 -debug < Regress_Score.v0.83.R
#ex) R --slave --vanilla --args --mutfile JamesTTNvariant.txt -output JamesTTNvariant.3000.regres_score.v0.84.txt -wide 3000 -ss3 -debug < Regress_Score.v0.84.R
#ex) R --slave --vanilla --args --mutfile JamesTTNvariant.txt -output JamesTTNvariant.wide1000.regres_score.v0.84.txt -wide 1000 -ss3 -debug < Regress_Score.v0.84.R
#ex) R --slave --vanilla --args --mutfile JamesTTNvariant.txt -output JamesTTNvariant.regres_score.v0.84.txt -entiregene -ss3 -debug -wide 1000 < Regress_Score.v0.84.R
#ex) R --slave --vanilla --args -mutfile test.input1.txt -output test.output.v0.87.txt < Regress_Score.v0.87.R
#ex) R --slave --vanilla --args -mutfile ./Exac/test.txt -output ./Exac/Exac.r0.1-0.3.TTN.SNV.txt.Regress_Score.v0.87.txt -esbTranscriptID ENST00000342992 -lseq ./Exac/ENST00000342992.txt -ss5 -skipRegressScore -skipSRE< Regress_Score.v0.88.R
#ex) R --slave --vanilla --args -mutfile temp.test.temp2.txt -output temp.test.temp2.output.v0.89 -localseq C:/Users/unknown/OneDrive/Copy/MAE-working/5ss-hotspot/test.lseqstore < Regress_Score.v0.89
#ex) R --slave --vanilla --args -mutfile temp.test.temp2.txt -output temp.test.temp2.output.v0.89 -refdirectory C:/Users/unknown/OneDrive/Copy/MAE-working/5ss-hotspot/ssfiles< Regress_Score.v0.90.R
#ex) R --slave --vanilla --args -mutfile temp.test.temp2.txt -output temp.test.temp2.output.v0.91 -useGTF19 E:/WorkDir/GTF < Regress_Score.v0.91.R
#ex) R --slave --vanilla --args -mutfile temp.test.temp2.txt -output temp.test.temp2.output.v0.91 -useGTF19 c:/Users/unknown/Dropbox/riken/Database -localseq ./temp91< Regress_Score.v0.91.R
#ex) R --slave --vanilla --args -mutfile test.reverse.txt -output test.reverse.output.v0.91 -useGTF19 c:/Users/unknown/Dropbox/riken/Database -localseq ./temp91< Regress_Score.v0.91.R
#ex) R --slave --vanilla --args -mutfile test.error.txt -output test.error.output.v0.91 -useGTF19 c:/Users/unknown/Dropbox/riken/Database -localseq ./temp91< Regress_Score.v0.91.R
#ex) R --slave --vanilla --args -mutfile test.error2.txt -output test.error2.output.v0.91 -useGTF19 E:/WorkDir/GTF < Regress_Score.v0.91.R
#ex) R --slave --vanilla --args -mutfile temp.test.temp2.txt -output temp.test.temp2.output.v0.92.txt -usehg19 c:/Users/unknown/Dropbox/riken/Database < Regress_Score.v0.92.R
#
#ex)  R --slave --vanilla --args -summarizeresultin alireza.temp.txt.regres_score.v0.88.txt -summarizeresult alireza.temp.txt.regress_score.v0.88 < TRegress_Score.v0.93.R
#ex)  R --slave --vanilla --args -summarizeresultin alireza.temp.txt.regres_score.v0.88.txt -summarizeresult alireza.temp.txt.regress_score.v0.88 -SummarizeByVariant < Regress_Score.v0.93beta.R
#ex)  R --slave --vanilla --args -summarizeresultin split.variant.list_0000.regress_score.v0.88.txt -summarizeresult temp.test -TakeOnePerVariant < Regress_Score.v0.93beta.R
#
#ex) R --slave --vanilla --args -sjdbin (result file of Regress score) -sjdbout (outout filename for sjdb file) -output (usual result outpt) -sjdbrestrictin (columnname,value) -sjdbrestrictout ( columnname,value)
#
#ex) R --slave --vanilla --args -mutfile C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/NEW60s.txt -output C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/NEW60s.output.v0.93.txt -summarizeresult C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/NEW60s.output.v0.93 -sjdbout C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/NEW60s.output.v0.93.sjdbtab -skipRegressScore -skipSRE -refdirectory C:/Users/unknown/OneDrive/Copy/MAE-working/5ss-hotspot/ssfiles -usehg19 c:/Users/unknown/Dropbox/riken/Database -takeCanonicalTranscript< C:/Users/unknown/OneDrive/Copy/MAE-working/5ss-hotspot/Regress_Score.v0.93zeta.R
#ex) R --slave --vanilla --args -mutfile E:/WorkDir/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/NEW60s.txt -output temp.test.txt -summarizeresult temp.test -sjdbout temp.test.sjdbtab -skipRegressScore -skipSRE -refdirectory E:/WorkDir/OneDrive/Copy/MAE-working/5ss-hotspot/ssfiles -usehg19 E:/WorkDir/Dropbox/riken/Database -takeCanonicalTranscript -skipcDNApos< E:/WorkDir/OneDrive/Copy/MAE-working/5ss-hotspot/Regress_Score.v0.93-+.R
#ex) R --slave --vanilla --args -mutfile C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/NEW60s.txt -output C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/temp.output.v0.94.txt -summarizeresult C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/temp.output.v0.94 -sjdbout C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/temp.output.v0.94.sjdbtab -skipRegressScore -skipSRE -refdirectory C:/Users/unknown/OneDrive/Copy/MAE-working/5ss-hotspot/ssfiles -usehg19 c:/Users/unknown/Dropbox/riken/Database -takeCanonicalTranscript -skipcDNApos< C:/Users/unknown/OneDrive/Copy/MAE-working/5ss-hotspot/Regress_Score.v0.94.R
#ex) R --slave --vanilla --args -mutfile C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/NEW60s.txt -output C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/temp.output.v0.94beta.txt -summarizeresult C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/temp.output.v0.94beta -sjdbout C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/temp.output.v0.94beta.sjdbtab -skipRegressScore -skipSRE -refdirectory C:/Users/unknown/OneDrive/Copy/MAE-working/5ss-hotspot/ssfiles -usehg19 c:/Users/unknown/Dropbox/riken/Database -takeCanonicalTranscript -skipcDNApos< C:/Users/unknown/OneDrive/Copy/MAE-working/5ss-hotspot/Regress_Score.v0.94beta.R
#ex) R --slave --vanilla --args -vcffile C:/Users/unknown/Dropbox/riken/workdir/test.test1.short.vcf -output  C:/Users/unknown/Dropbox/riken/workdir/test.test1.short.vcf.output.v0.95.txt -summarizeresult  C:/Users/unknown/Dropbox/riken/workdir/test.test1.short.vcf.output.v0.95 -sjdbout  C:/Users/unknown/Dropbox/riken/workdir/test.test1.short.vcf.output.v0.94beta.sjdbtab -skipRegressScore -skipSRE -refdirectory C:/Users/unknown/OneDrive/Copy/MAE-working/5ss-hotspot/ssfiles -usehg19 c:/Users/unknown/Dropbox/riken/Database -takeCanonicalTranscript -skipcDNApos< C:/Users/unknown/OneDrive/Copy/MAE-working/5ss-hotspot/Regress_Score.v0.95.R
#ex) R --slave --vanilla --args -mutfile C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/NEW60s.txt -output C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/temp.output.v0.96-.txt -summarizeresult C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/temp.output.v0.96- -sjdbout C:/Users/unknown/OneDrive/Copy/Minigeneassay/no-vector/CHD-denovo_matome/temp.output.v0.96-.sjdbtab -skipRegressScore -skipSRE -refdirectory C:/Users/unknown/OneDrive/Copy/MAE-working/5ss-hotspot/ssfiles -usehg19 c:/Users/unknown/Dropbox/riken/Database -takeCanonicalTranscript -skipcDNApos< C:/Users/unknown/OneDrive/Copy/MAE-working/5ss-hotspot/Regress_Score.v0.96-.R
#
#Options:
#-seq (a file contains sequence data)  exon should be UPPERCASE, intron should be LOWERCASE, you can indicate mutation like (A|B), where A means refseq and B mean alt seq
#                                       Also you can indicate nucleutide of interest like ATGC[ATGC]ATGCa              
#-vcffile (vcf file) indicates a variant information file as a input 
#     **If TRANSCIPT_ID=ENSTXXXXXXXXX, pick the ensembl ID for transcript ID
#-mutfile (chrpos format mut file with transcript ID) e.g.chr1:123G>A ENST000001234 in v0.5 just 1 line accept
#-output (filename to be saved as a result)
#-marthost (indicates biomaRt server)
#-esbTranscriptID (transcriptID)  when you check cDNA variant list, you should provide transcriptID
#-nomart        if you indicate this switch, connecting biomaRt is being skipped.
#-canonicalSS  indicate this script to process canonical 5ss and 3ss around the given mutation, if the mutation is not on the exon, this function does not work.
#-entiregene indicates this script wull calculate entire gene score 
#-ss5 indicates this script to calculate ss5 scores
#-ss3 indicates this script to calculate ss3 scores
#     if omitted the both above automatically calculate both
#-wide (integer) indicates to make scan region added (integer) around the mutation
#-lseq (filename) indicagtes transcript info generated by GetSeq(TranscriptID,1) if you employ this, you don't need to connect biomart server , which make run time short.
#                 whme you use it, I recommned you to use -esbTranscriptID to indicate the name of transcript.
#-skipRegressScore indicates skipping regress score calculation to make time short. The score will be zero
#-skipSRE indicates skipping checking # of splicing regulatory elements. The # of elements will be zero
#-localseq (foldername) is a siwtch to reduce the access to biomaRt server. In concrete,
#                      1) if required seq is on the indicated local folder, employ it
#                      3) if not, retrieve the seq from biomaRt server and save the seq on the indicated local server.
#-refdirectory (directory where maxent and other files are located) :when these files are no loacted on the current directory , please indicate their locations
#-useGTF19 (directory that contains Homo_Sapiens.GRCh37.75,gtf for ensemble hg19 gtf file downloaded extarct from http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz)
#-usehg19 (directory that contains hg19.transcripts.txt and hg19.exons.txt) those 2 files are generated from Homo_Sapiens.GRCh37.75,gtf by Break_Ensembl_GTF.pl
#-useGRCh38  (directory that contains GRCh38.transcripts.txt and GRCh38.exons.txt) those 2 files are generated from Homo_Sapiens.GRCh38.87,gtf by Break_EnsemblGRCh38_GTF.pl
#
#-summarizeresult (prefix of output file for summarizeresult) one variant 9 scores for 5SS and 23scores for 3SS. Then just pick siginificant lines file *.SigLines.txt will be generated 
#-summarizeresultin (raw result file of Regress_Score) when you want just to summarize the result, please indicate the result file, which mean skipping calculations.
#-SummarizeByVariant :when you indicate -summarizresult, this switch make this script generate files. *.SigLines.txt and *.SigLinesByVariant.txt
#-TakeOnePerVariant indicates that the script dumps no-significant lines and integrates significant lines to make one line for one variant
#-skip1stMiningStep indicates skipping 1st procedure when you indicate -summarizresult 
#-skip2ndMiningStep indicates skipping 1st procedure when you indicate -summarizresult
#
#-sjdbout (sjdbfilename for STAR) this switch indicates to make sjdbFile for STAR RNAseq aligner --sjdbFileChrStartEnd. also result file will have Predicted_SpliceJunction column
#      ****when use this switch, you also need to indicate -summarizeresult. That's because this requires the column <Kind>, which is attached by -summarizeresult
##-sjdbrestrictin (columnname, value) gives an inclusion creteria when -sjdb being indicated. just lines with the column = value will be included
#-sjdbrestrictout (columnname, value) gives an exclusion creteria when -sjdb being indicated. just lines with the column = value will be excluded
#-sjdbin (result file of Regress_Score) if you want to skip score calculation and give the result file directory, use this switch.
#
#-takeCanonicalTranscript indicates when multiple transcripts found, take a canonical transcript 
# #the definition of Canonical Transcript means 
#  http://asia.ensembl.org/Help/Glossary?id=346
#  For human, the canonical transcript for a gene is set according to the following hierarchy: 1. Longest CCDS translation with no stop codons. 2. If no (1), choose the longest Ensembl/Havana merged translation with no stop codons. 3. If no (2), choose the longest translation with no stop codons. 4. If no translation, choose the longest non-protein-coding transcript.
#  ##CAUTION this switch currently works when you employ -usehg19
#-skipcDNApos indicates skipping providing variant name in cDNA positon format (e.g. c.123A>G)when the original variant name is given in chromosome position format (e.g. chr1:12345A>G) 
#
#most of element search sub routines came from DataCollecterPartA.v0.9.R and they had been modified to fit this current situation.
#most of constructing ref seq and searching alt seq code came from 5ss_cal_ver1.19.R
#
#History
#v0.5 launch ver.
#v0.6 can process submitted sequence.   can process nearby authentic SS.
#v0.7 can process both chr pos and cDNA format
#v0.8 can indicate to process canonical splice site using -canonicalSS
#v0.81 can indicate to process entire gene using -entiregene
#v0.82 can indicate either 5ss or 3ss calculation
#v0.83 bug-fix for reverse strand gene (fix the process going up accoding to phypos-> according to cDNApos)
#v0.84 -wide check added
#v0.85 when transcripID NOT given in chrpos mode, automatically find all transcritps for it and try all of them
#         (change sub-routine findtranscript and chrposmut)
#     skip variants NOT on the transcript e.g. downstream varint or UTRs...
#     bug-fix
#v0.86 check increase or decrease # of ESE,ESS, Grun and ISE 250bp around the mut for comprehensive analysis. just columns added in the result. So you can employ old-scripts for v0.84 or 0.85
#v0.87 bug-fix  v0.86 counts # of splicing elements on the putative exon-intron generated by putative splicing variants. In this ver, just count # of elements around the variant
#       #ESE,ESS derived from exon where the mutation is.      #Grun,ISE derived from intron where the mutation is.  e.g.)if the mutation is on exon, #Grun, #ISE must be 0.
# Accept ENSEMBL transcriptID 9e.g.ENST0000001234 and RefseqmRNAID e.g.NM_001195597  (if NM_001195597.2 given, automatically take just NM_001195597 to avoid error)
#     when 1 NM ID corresponds to multiple ENST ID, just take the 1st one (->GetSeq)
#v0.88 stand alone mode added, which means when you process one transcript and you have seq file for it, just indicate that and read lseq info from file and not connect to biomaRt
#     it can speed up process and run more nodes
#v0.89 -localseq (foldername) switch added, by which sequence retreivig from biomaRT will be omitted if the local folder has the seq file. if not, reireive from biomaRt
#v0.90 -refdirectory (foldername) switch added. the indicated folder should contain ssfiles
#v0.91 -useGTF19 (directory that contains Homo_Sapiens.GRCh37.75,gtf) swithch added. when use this switch, hg19 unsembleGTF file is employed. then sequence was retrived from ucsc human hg19. That means biotmaRT not required anymore
#       STABLE fix forward & reverse
#       STABLE one gene for several geneid e.g.HSD17B8  if you make use of this switch, no protein transcirpt will be ignored.
#       STABLE for one exon gene (skipped)
#       STABLE for patch ID as long as GetSeq function contains it.
#v0.92 -usehg19 (directory that contains hg19.transcripts.txt and hg19.exons.txt) those 2 files are extracted from hg19 GTF file. That means no biomaRt required.
#v0.93beta -summarizeresult and related switches added from TakeSigLines3.R
#   beta means not tested if the script ran without error
#v0.93theta -sjdbout implemented from sjdbMake2.R
#v0.93zeta -takeCanonicalTranscript added but works just with -usehg19
#v0.93 you can indicate different transcriptIDs for cDNA-mut frame 
#v0.94 phypos and cDNApos provided automatically
#     provide information for #exon/#last exon  #CAUTION currently, works just when -usehg19
#     Decision for Significant variant came to consider the position of the variant i.e. a variant on the last exon can not create a splice donor since the corresponding acceptor does not exist. or vice versa
#v0.95 -useGRCh38 added
#        provide  in-frame or frame-shift change info in -sjdbout command
#      -vcffile added
#v0.96 for RNAseqcompare, Predicted_NormalJunction added  #CAUTION in case of LOSS variant, you will find nearby junction for comparison (normal junction found in Predicted_SpliceJunction)
#v0.97 bug-fix for vcf read
#In the future version, Using STAR RNAseq aligners sjdb file, p-value will be attached for the splicing variant
#
"


#Subs
inputparameters<-function(inputlist,inputparas) {  #input parameters from commnand line
  if (length(inputlist) != length(inputparas)) stop("The length NOT match in the input parameters\n")
  inputparasrev<-rep(-9,length(inputparas))  #check value if the value is already input or not
  for (i in 1:length(commandArgs(trailingOnly=TRUE))) { 
    for (j in 1:length(inputlist)) {
      if (commandArgs(trailingOnly=TRUE)[i] == inputlist[[j]] & inputparasrev[j]==-9) {
        inputparasrev[j]<-1
        if (inputparas[j]==0) {  #input sw 1
          inputlist[[j]]<-1 
        } else { #input a series of value or character
          inputlist[[j]]<-commandArgs(trailingOnly=TRUE)[i+1]
        }
      }
    }    
  }
  
  for (i in 1:length(inputparasrev)) {  #when no input, then assign -9 for null parameter (-9 is a null value)
    if (inputparasrev[i]==-9) inputlist[[i]]<-as.numeric(-9)
  }
  return (inputlist)
}

positionSeq<-function(Seq,start,end, oristart,oriend,strand) {  # this script is to figure out where actual start and end position in the given sequence and return corresponding region indicated by start-end   #given that seq spans oristart-oriend, calculate where start and end locates on the seq.
  if (start>0 & end>0) {
    if (strand==1) { # forward strand
      acstart<-start-oristart+1
      acend<-end-oristart+1
    } else { #reverse strand
      acstart<-oriend-end+1
      acend<-oriend-start+1
    }
    return(substr(Seq,acstart,acend))
  } else {
    return("")
  }
}

#sub routine for 5ss and 3ss maxent
makemaxentscores3<-function(){
  list<-c('me2x3acc1','me2x3acc2','me2x3acc3','me2x3acc4','me2x3acc5','me2x3acc6','me2x3acc7','me2x3acc8','me2x3acc9')
  metables<-list()
  num<-0
  initial<-1
  for (i in list) {
    vec<-c()
    SCOREF<-readLines(con=directory(i))
    for (j in 1:length(SCOREF)) {
      vec<-c(vec,as.numeric(SCOREF[j]))
    }
    new<-list(val=I(vec))
    names(new)<-num
    if (initial==1) {metables<-new;initial<-0;} else metables<-c(metables,new)  # row correspond to n and col correspond to num
    num<-num+1
  }
  return (metables)
}

makemaxentscores5<-function(){
  me2x5<-read.table(directory("me2x5"),header=F,stringsAsFactors=F)  #if you want to make the running time short, you should make this variable in the first of your script
  seq<-read.table(directory("splice5sequences"),header=F,stringsAsFactors=F)
  scorematrix<-cbind(seq,me2x5)
  names(scorematrix)<-c("Seq","Score")
  return(scorematrix)
}

hashseq<-function(seq){
  seq<-toupper(seq)
  seq<-gsub("A","0",seq)
  seq<-gsub("C","1",seq)
  seq<-gsub("G","2",seq)
  seq<-gsub("T","3",seq)
  seqa<-unlist(strsplit(seq,""))
  sum<-1
  len<-nchar(seq)
  four<-c(1,4,16,64,256,1024,4096,16384)
  i<-0
  while (i<len) {
    sum<-sum+as.numeric(seqa[i+1])*four[len-i]
    i<-i+1
  }
  return(sum)
}

getrest<-function(seq) {
  return(paste(substr(seq,1,18),substr(seq,21,23),sep=""))
}

scoreconsensus3sub<-function(seq) {
  seqa<-unlist(strsplit(seq,""))
  bgd<-list(A=0.27,C=0.23,G=0.23,T=0.27)
  cons1<-list(A=0.9903,C=0.0032,G=0.0034,T=0.0030)
  cons2<-list(A=0.0027,C=0.0037,G=0.9905,T=0.0030)
  addscore<-as.numeric(cons1[seqa[19]])*as.numeric(cons2[seqa[20]])/(as.numeric(bgd[seqa[19]])*as.numeric(bgd[seqa[20]]))
  return(addscore)
}

scoreconsensus5sub<-function(inputseq) {
  sc34<-data.frame(ATGC=c("A","C","G","T"),
                   bgd=c(0.27,0.23,0.23,0.27),
                   cons1=c(0.004,0.0032,0.9896,0.0032),
                   cons2=c(0.0034,0.0039,0.0042,0.9884))
  chr3<-substr(inputseq,4,4)
  chr4<-substr(inputseq,5,5)
  score<-sc34[sc34$ATGC==chr3,]$cons1*sc34[sc34$ATGC==chr4,]$cons2/(sc34[sc34$ATGC==chr3,]$bgd*sc34[sc34$ATGC==chr4,]$bgd); 
  return(score)
}

maxentscore3sub<-function(seq,metables) {
  sc0<-as.numeric(unlist(metables[1])[hashseq(substr(seq,1,7))])
  sc1<-as.numeric(unlist(metables[2])[hashseq(substr(seq,8,14))])
  sc2<-as.numeric(unlist(metables[3])[hashseq(substr(seq,15,21))])
  sc3<-as.numeric(unlist(metables[4])[hashseq(substr(seq,5,11))])
  sc4<-as.numeric(unlist(metables[5])[hashseq(substr(seq,12,18))])
  sc5<-as.numeric(unlist(metables[6])[hashseq(substr(seq,5,7))])
  sc6<-as.numeric(unlist(metables[7])[hashseq(substr(seq,8,11))])
  sc7<-as.numeric(unlist(metables[8])[hashseq(substr(seq,12,14))]) 
  sc8<-as.numeric(unlist(metables[9])[hashseq(substr(seq,15,18))]) 
  
  finalscore<-sc0*sc1*sc2*sc3*sc4/(sc5*sc6*sc7*sc8)
  return(finalscore)
}

#main function of maxent
maxent3score<-function(seq) {
  return(log2(scoreconsensus3sub(seq)*maxentscore3sub(getrest(seq),metables)))
}

maxent5score<-function (seq) {
  getrest<-paste(substr(seq,1,3),substr(seq,6,9),sep="")
  return(log2(scoreconsensus5sub(seq)*as.numeric(scorematrix[scorematrix$Seq==getrest,]$Score)))
}

#checking accessory introns e.g. G run, ISE
intronscan<-function(seq,beforeafter,skip) {   #seq -intronseq ,beforeafter  the intron exists before a exon of interest (-1) or after a exon of interest (1) ,skip 20(5 end of exon) or 6 (3 end of exon)
  output<-list(Grun=0,ISE=0)  #reply scores of intron elements
  gflag<-0
  scanstart<-1+(beforeafter==1)*skip
  scanend<-nchar(seq)-(beforeafter==-1)*skip
  for (i in scanstart:scanend) {
    basepos<-ifelse(beforeafter==1,i,nchar(seq)-i+1)
    if (substr(seq,i,i)=="G") { #checking G run
      gflag<-gflag+1
      if (i==nchar(seq) & gflag>2) {
        output$n_Grun<-output$n_Grun+1
        pos<-basepos-gflag*(beforeafter==1)  #pos means the distance from the exon
        output$Grun<-output$Grun+PosScore(pos,ifelse(beforeafter==1,Eweight$aGrun,Eweight$bGrun))
      }
    } else if (gflag<3) {
      gflag<-0
    } else {
      output$n_Grun<-output$n_Grun+1
      pos<-basepos-gflag*(beforeafter==1)  #pos means the distance from the exon
      output$Grun<-output$Grun+PosScore(pos,ifelse(beforeafter==1,Eweight$aGrun,Eweight$bGrun))
      gflag<-0
    }
    tempISE9<-ifelse ((i+9-1)<=scanend,substr(seq,i,i+9-1),"dummy") #checking ISE9 seq
    if (length(ISE9seq[ISE9seq==tempISE9])>0) {
      output$n_ISE<-output$n_ISE+1
      pos<-basepos-8*(beforeafter==-1) #pos means the distance from the exon
      output$ISE<-output$ISE+PosScore(pos,ifelse(beforeafter==1,Eweight$aISE,Eweight$bISE))
     }
    tempISE10<-ifelse ((i+10-1)<=scanend,substr(seq,i,i+10-1),"dummy") #checking ISE10 seq
    if (length(ISE10seq[ISE10seq==tempISE10])>0) {
      output$n_ISE<-output$n_ISE+1
      pos<-basepos-9*(beforeafter==-1) #pos means the distance from the exon
      output$ISE<-output$ISE+PosScore(pos,ifelse(beforeafter==1,Eweight$aISE,Eweight$bISE))
    }    
  }
  return(output)
}

#exonscan 3ss-seq & 5ss-seq, ESE and ESS on the exon, G run on both introns surrounding the exon
exonscan<-function(seq,from) {  #Seq a sequence from the edge of exon,  from 5 or 3
  maxlen<-nchar(seq)
  exonlen<-maxlen
  output<-list(ESE=0,ESS=0)  #reply scores of exon elements
  for (j in 1:(maxlen-min(ESElength,ESSlength)))  {#checking ESE &ESS seqs on the exon
    
    tempESE<-ifelse ((j+ESElength-1)<=maxlen,substr(seq,j,j+ESElength-1),"dummy")
    tempESS<-ifelse ((j+ESSlength-1)<=maxlen, substr(seq,j,j+ESSlength-1),"dummy")
    #ESpos<-ifelse(Strand==1,start+j-1,end-j+1)  #ES position is given by phypos
    ESpos<-j  #ES position is given by the distance from 5'start of the exon  start pos is 1
    if (length(ESEseq[ESEseq==tempESE])>0) {
      ESposR<-exonlen-ESpos-ESElength+1 #ESposR means a distance from 3' end of the exon
      pos<-ifelse(from==5,ESpos,ESposR)
      output$ESE<-output$ESE+PosScore(pos,ifelse(from==5,Eweight$ESE5,Eweight$ESE3))
    }
    if (length(ESSseq[ESSseq==tempESS])>0) {
      ESposR<-exonlen-ESpos-ESSlength+1 #ESposR means a distance from 3' end of the exon
      pos<-ifelse(from==5,ESpos,ESposR)
      output$ESS<-output$ESS+PosScore(pos,ifelse(from==5,Eweight$ESS5,Eweight$ESS3))
    }
  }
  return(output)
}

PosScore<-function(pos,ScoreVec) {  #calculate a score for each element, ScoreVec is "Mean,SD" (should be given as charcter) for gaussian model.  multi Mean,Sds are acceptable e.g. mean,sd,mean,sd,....
  meansd<-as.numeric(unlist(strsplit(ScoreVec,",")))
  score<-sum(dnorm(pos,mean=meansd[seq(1,length(meansd),by=2)],sd=meansd[seq(2,length(meansd),by=2)]))
  return(score)
}

SSscan<-function(seq,ExPosInTheContext,SS)  {#seq a sequence including 5ss or 3ss, ExPosInTheContext exon position(Start,End) in the given seq  SS should be 5(5Splice site) or 3 (3Spice site) 
  startend<-as.numeric(unlist(strsplit(ExPosInTheContext,",")))
  seq<-toupper(seq)
  output<-list(Maxent=0,SSscore=0,SSseq="")
  if (SS==5) { #5ss processing
    ss5seq<-substr(seq,startend[2]-2,startend[2]+6)
    if (nchar(ss5seq)==9) output$Maxent<-maxent5score(ss5seq) else output$Maxent<--99
    #cat("Submit exon:",substr(seq,startend[1],startend[2]), " Submit intron:",substr(seq,startend[2]+1,nchar(seq)),sep="")
    output$SSscore<--99
    if (IN$skipRegressScore==-9) {
      Escore<-exonscan(substr(seq,startend[1],startend[2]),3)
      Iscore<-intronscan(substr(seq,startend[2]+1,nchar(seq)),1,6)
      logOdds<-RegCoef5$Intercept+RegCoef5$MAXENT5*output$Maxent+RegCoef5$ESE3*Escore$ESE+RegCoef5$ESS3*Escore$ESS+RegCoef5$aGrun*Iscore$Grun+RegCoef5$aISE*Iscore$ISE
      if (output$Maxent!=-99) output$SSscore<-1/(1+exp(-logOdds))
    }
    output$SSseq<-paste(substr(ss5seq,1,3),tolower(substr(ss5seq,4,9)),sep="")
    #cat("SS5seq:",substr(ss5seq,1,3),tolower(substr(ss5seq,4,9)),"\n",sep="")
  } else { #3ss processing
    ss3seq<-substr(seq,startend[1]-20,startend[1]+2)
    #cat("SS3seq:",ss3seq,"\n")
    if (nchar(ss3seq)==23) output$Maxent<-maxent3score(ss3seq) else output$Maxent<--99
    #cat("Submit exon:",substr(seq,startend[1],startend[2]), " Submit intron:",substr(seq,1,startend[1]-1),sep="")
    output$SSscore<--99
    if (IN$skipRegressScore==-9) {
      Escore<-exonscan(substr(seq,startend[1],startend[2]),5)
      Iscore<-intronscan(substr(seq,1,startend[1]-1),-1,20)
      logOdds<-RegCoef3$Intercept+RegCoef3$MAXENT3*output$Maxent+RegCoef3$ESE5*Escore$ESE+RegCoef3$ESS5*Escore$ESS+RegCoef3$bGrun*Iscore$Grun+RegCoef3$bISE*Iscore$ISE
      if (output$Maxent!=-99) output$SSscore<-1/(1+exp(-logOdds)) 
    }
    output$SSseq<-paste(tolower(substr(ss3seq,1,20)),substr(ss3seq,21,23),sep="")
    #cat("SS3seq:",tolower(substr(ss3seq,1,20)),substr(ss3seq,21,23),"\n",sep="")
  } 
  return(output)
}

SEEscan<-function(seq) {
  output<-list(nESE=0,nESS=0)
  maxlen<-nchar(seq)
  for (j in 1:(maxlen-min(ESElength,ESSlength)))  {#checking ESE &ESS seqs on the exon
    tempESE<-ifelse ((j+ESElength-1)<=maxlen,substr(seq,j,j+ESElength-1),"dummy")
    tempESS<-ifelse ((j+ESSlength-1)<=maxlen, substr(seq,j,j+ESSlength-1),"dummy")
    #ESpos<-ifelse(Strand==1,start+j-1,end-j+1)  #ES position is given by phypos
    ESpos<-j  #ES position is given by the distance from 5'start of the exon  start pos is 1
    if (length(ESEseq[ESEseq==tempESE])>0) output$nESE<-output$nESE+1
    if (length(ESSseq[ESSseq==tempESS])>0) output$nESS<-output$nESS+1
  }
  return(output)
}

SEIscan<-function(seq) {
  output<-list(nGrun=0,nISE=0)
  gflag<-0
  for (i in 1:nchar(seq)) {
    if (substr(seq,i,i)=="G") { #checking G run
      gflag<-gflag+1
      if (i==nchar(seq) & gflag>2) {
        output$nGrun<-output$nGrun+1
      }
    } else if (gflag<3) {
      gflag<-0
    } else {
      output$nGrun<-output$nGrun+1
      gflag<-0
    }
    tempISE9<-ifelse ((i+9-1)<=nchar(seq),substr(seq,i,i+9-1),"dummy") #checking ISE9 seq
    if (length(ISE9seq[ISE9seq==tempISE9])>0) {
      output$nISE<-output$nISE+1
    }
    tempISE10<-ifelse ((i+10-1)<=nchar(seq),substr(seq,i,i+10-1),"dummy") #checking ISE10 seq
    if (length(ISE10seq[ISE10seq==tempISE10])>0) {
      output$nISE<-output$nISE+1
    }    
  }
  return(output)
}

SEscan<-function(newseq,cDNApos) {
  exonseq<-"";
  intronseq<-"";
  #name intron as -1 , -2 ,-3...
  if (cDNApos[1]<0 | IN$skipSRE==1) {
    return(list(nESE=-99,nESS=-99,nGrun=-99,nISE=-99))
  } else {
    #name intron to identify each intron separately
    exons<-unique(newseq$exon)
    exons<-sort(exons[exons>0])
    for (i in exons) {
      minpos<-min(which(newseq$exon==i))
      maxpos<-max(which(newseq$exon==i))
      if (i<max(exons)) {
        nextmin<-min(which(newseq$exon==i+1))
        newseq$exon[(maxpos+1):(nextmin-1)]=-i
      }
      if (i==1) {
        if (minpos>1) newseq$exon[1:(minpos-1)]=0
      } 
    }
  takeregion<-unique(newseq$exon[cDNApos[1]:cDNApos[2]])
  takeexon<-takeregion[takeregion>0]
  takeintron<-takeregion[takeregion<1]
  SEE<-list(nESE=0,nESS=0)
  if (length(takeexon)>0) {
    seq<-""
    for (i in takeexon) {
      seq<-paste(seq,paste(newseq$Seq[newseq$exon==i],collapse=""),"XXXXXXXXXXXX")
    }
    SEE<-SEEscan(seq)
  }
  SEI<-list(nGrun=0,nISE=0)
  if (length(takeintron)>0) {
    seq<-""
    for (i in takeintron) {
      seq<-paste(seq,paste(newseq$Seq[newseq$exon==i],collapse=""),"XXXXXXXXXXXX")
    }
    SEI<-SEIscan(seq)
  }
  return(list(nESE=SEE$nESE,nESS=SEE$nESS,nGrun=SEI$nGrun,nISE=SEI$nISE))
  }
}

SS5scan<-function(newseq,Start,End,cDNApos) {
  #cat("Seq: ",tolower(paste(newseq$Seq[1:(Start-1)],collapse="")),paste(newseq$Seq[Start:End],collapse=""),tolower(paste(newseq$Seq[(End+1):nrow(newseq)],collapse="")), " Start:",Start," End:",End,"\n")
  result5<-data.frame()
  if (Start>End) {temp<-c(Start,End);Start<-temp[2];End<-temp[1]}
  for (j in Start:End) {
    exonend<-j+2    #putative exon end generated by the mutation
    if (exonend>=3 & exonend<=nrow(newseq)) { 
    exonrank<-newseq$exon[exonend] # #exon where the above exon end exists
    if (exonrank==-9 ) exonrank<-max(newseq$exon[1:exonend])
    if (exonrank>0) {
      exonstart<-min(which(newseq$exon==exonrank))   #  exon start pos in the context of seq
    } else {
      exonstart<-1
    }
    a<-which(newseq$exon==(exonrank+1));if (length(a)>0) intronend<-min(a)-1 else  intronend<-nrow(newseq)
    ###landscape is   seq(exonstart....exonend....intronend)
    exonstart_posincontext<-1
    if ((exonend-exonstart+1)<250) {
      takestart<-exonstart;exonend_posincontext<-exonend-exonstart+1
      #cat("5within\n")
    } else {
      takestart<-exonend-250;exonend_posincontext<-250+1
    }
    if ((intronend-exonend)<250) {
      takeend<-intronend
      #cat("3within\n")
    } else {
      takeend<-exonend+250
    }
    seq<-paste(newseq$Seq[takestart:takeend],collapse ="")
    #cat("j;",j,"takestart",takestart,"takeend ", takeend,"\n")
    #cat("Seq: ",tolower(substr(seq,1,exonstart_posincontext-1)),substr(seq,exonstart_posincontext,exonend_posincontext),tolower(substr(seq,exonend_posincontext,nchar(seq))),"\n",sep="")
    output<-SSscan(seq,paste(exonstart_posincontext,exonend_posincontext,sep=","),5)
    cDNAboundary<-paste("c.",newseq$cDNApos[exonend],"-",newseq$cDNApos[exonend+1],sep="")
    phyposboundary<-paste("chr:",newseq$chr[exonend],":",newseq$phypos[exonend],"-",newseq$phypos[exonend+1],sep="")
    output2<-SEscan(newseq,cDNApos)
    #cat("len newseqphypos",length(newseq$phypos[j])," nrow output:",unlist(output), " length cdnaboundary:",length(cDNAboundary)," length phyposboundary:",length(phyposboundary),"\n")
    new<-data.frame(phypos=newseq$phypos[j],MaxEntScoreSS5=output$Maxent,SS5score=output$SSscore,SS5cDNAboundary=cDNAboundary,SS5phyposboundary=phyposboundary,SS5seq=output$SSseq,ESE5=output2$nESE,ESS5=output2$nESS,Grun5=output2$nGrun,ISE5=output2$nISE)
    result5<-rbind(result5,new)
    }
  }
  return(result5)
}

SS3scan<-function(newseq,Start,End,cDNApos) {
  #cat("Seq: ",tolower(paste(newseq$Seq[1:(Start-1)],collapse="")),paste(newseq$Seq[Start:End],collapse=""),tolower(paste(newseq$Seq[(End+1):nrow(newseq)],collapse="")), " Start:",Start," End:",End,"\n")
  result3<-data.frame()
  if (Start>End) {temp<-c(Start,End);Start<-temp[2];End<-temp[1]}
  for (j in Start:End) {
    exonstart<-j+20    #putative exon end generated by the mutation
    if (exonstart>=21 & exonstart<=nrow(newseq)) { 
      #cat("exonstart:",exonstart,"\n")
    exonrank<-newseq$exon[exonstart] # #exon where the above exon end exists
    if (exonrank==-9) {
      temp<-sort(unique(newseq$exon[exonstart:nrow(newseq)]))
      if (length(temp)>1) exonrank<-temp[2] 
    }  
    if (exonrank>0) {
      exonend<-max(which(newseq$exon==exonrank))   #  exon end pos in the context of seq
    } else {
      exonend<-nrow(newseq)
    }
    a<-which(newseq$exon==(exonrank-1));if (length(a)>0) intronstart<-max(a)+1 else  intronstart<-1
    ###landscape is   seq(intronstart...exonstart...exonend)
    
    if ((exonstart-intronstart)<250) {
      takestart<-intronstart;exonstart_posincontext<-exonstart-intronstart+1
    } else {
      takestart<-exonstart-250;exonstart_posincontext<-250+1
    }
    if ((exonend-exonstart)<250) {
      takeend<-exonend;
    } else {
      takeend<-exonstart+250-1
    }
    seq<-paste(newseq$Seq[takestart:takeend],collapse ="")
    exonend_posincontext<-nchar(seq)
    #cat("Seq: ",tolower(substr(seq,1,exonstart_posincontext-1)),substr(seq,exonstart_posincontext,exonend_posincontext),tolower(substr(seq,exonend_posincontext,nchar(seq))),"\n",sep="")
    output<-SSscan(seq,paste(exonstart_posincontext,exonend_posincontext,sep=","),3)
    cDNAboundary<-paste("c.",newseq$cDNApos[exonstart-1],"-",newseq$cDNApos[exonstart],sep="")
    phyposboundary<-paste("chr:",newseq$chr[exonstart],":",newseq$phypos[exonstart-1],"-",newseq$phypos[exonstart],sep="")
    output2<-SEscan(newseq,cDNApos)
    new<-data.frame(phypos=newseq$phypos[j],MaxEntScoreSS3=output$Maxent,SS3score=output$SSscore,SS3cDNAboundary=cDNAboundary,SS3phyposboundary=phyposboundary,SS3seq=output$SSseq,ESE3=output2$nESE,ESS3=output2$nESS,Grun3=output2$nGrun,ISE3=output2$nISE)
    result3<-rbind(result3,new)
    }
  }
  return(result3)
}

GetSeq<-function(transciptID,cDNAsw) { #ID -transcriptID, Seq -each allele, phypos -physical position,cDNApos-cDNA position, exon-#exon
  if (IN$lseq==-9) {
    transcriptkind<-""
    if (length(grep("ENST",transciptID))>0) transcriptkind<-"ensembl_transcript_id"
    if (length(grep("NM_",transciptID))>0) {
      transcriptkind<-"refseq_mrna"
      if (length(grep("\\.",transciptID))>0) {
        temp<-unlist(strsplit(transciptID,"\\."))
        transciptID<-temp[1]
      }
    }
    if (IN$localseq!=-9) {
      files<-list.files(path=IN$localseq,pattern=transciptID)
      if (transciptID%in%files) {
        cat("Found Seq data for transcipt:",transciptID,"on the folder:",IN$localseq,"\n")
        lseq<-read.table(paste0(IN$localseq,"/",transciptID),sep="\t",header=T,stringsAsFactors=F)
        return(lseq)
      }
    }
    if (IN$useGTF19!=-9) {#v0.91  transcriptkind should be ensembl_transcript_id
      trID<-as.character(transciptID)
      gene_id<-as.character(getGenePositions(unifyJuncs(getSpliceTable(extractTranscript(ens,trID))))$gene_id[1])
      gm<-geneModel(ens, gene_id)
      tr<-getTranscript(gm,trID)
      trstartend<-attr(tr,"stcodon")
      if (trstartend[1]<trstartend[2]) {      #stcodon[1] indicates the 1st pos of start codon and stdodon[2] indicates the 1st codon of stop codon. Then to cover the entire stop codon, stcodon[2] should be +2
        trstartend[2]<-trstartend[2]+2 
      } else {  #in case of reverse trancript
        trstartend[1]<-trstartend[1]+2;trstartend<-c(trstartend[2],trstartend[1])
      }
      exons<-getExonData(tr)
      exons$seqid<-sapply(exons$seqid,function(x){ #in case of the chrmosome is patch id e.g.HG1079_PATCH
        if (x%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT")) return(x)
        transframe<-data.frame(patchid=c("HG104_HG975_PATCH","HG1079_PATCH","HG1082_HG167_PATCH","HG1211_PATCH","HG122_PATCH","HG1257_PATCH","HG1304_PATCH","HG1308_PATCH","HG1350_HG959_PATCH","HG1433_PATCH","HG1436_HG1432_PATCH"
                                         ,"HG1459_PATCH","HG1497_PATCH","HG1501_PATCH","HG185_PATCH","HG19_PATCH","HG243_PATCH","HG280_PATCH","HG306_PATCH","HG329_PATCH","HG344_PATCH","HG388_HG400_PATCH","HG414_PATCH","HG745_PATCH"
                                         ,"HG873_PATCH","HG979_PATCH","HG991_PATCH","HSCHR16_2_CTG3_1","HSCHR17_3_CTG4","HSCHR17_4_CTG4","HSCHR19LRC_COX1_CTG1","HSCHR19LRC_LRC_I_CTG1","HSCHR19LRC_LRC_S_CTG1","HSCHR19LRC_PGF1_CTG1"
                                         ,"HSCHR19LRC_PGF2_CTG1","HSCHR2_2_CTG12","HSCHR22_2_CTG1","HSCHR4_2_CTG9","HSCHR6_MHC_COX","HSCHR6_MHC_DBB","HSCHR6_MHC_SSTO","HSCHR19LRC_LRC_J_CTG","HG871_PATCH","HG1595_PATCH"
                                         ,"HG1592_PATCH","HG748_PATCH","HG237_PATCH","HG186_PATCH","HG1091_PATCH","HG174_HG254_PATCH","HG357_PATCH","HG115_PATCH","HG1699_PATCH","HG1502_PATCH","HG311_PATCH"),
                               chr=c("8","19","5","10","11","7","6","7","19","X","X","X","X","9","17","8","8","3","11","22","12","11","11","17","11","10","3","16","17","17","19","19","19","19","19","2","22","4","6","6","6"
                                     ,"19","10","12","14","17","21","3","3","4","6","7","8","9","10")
                                )
        newid<-as.character(transframe$chr[transframe$patchid==x])
        if (length(newid)==0) {
          match<-regexpr("CHR(\\d+|X|Y|MT)",x)
          if (match>0) newid<-substr(x,match+3,match+attr(match,"match.length")-1)
        }
        cat("PATCH ID:", x, "was replaced by ",newid,"\n")
        return (newid)
      })
      exonsinfo<-data.frame(ensembl_transcript_id=rep(trID,nrow(exons)),ensembl_exon_id=exons$exon_number,rank=exons$exon_number,chromosome_name=sapply(exons$seqid,function(x){paste0("chr",x)}),exon_chrom_start=exons$start,exon_chrom_end=exons$end)
      startend<-c(exons$start,exons$end)
      transcriptinfo<-data.frame(hgnc_symbol=attributes(geneModel(ens, gene_id))$gene_name[1],ensembl_transcript_id=transciptID,chromosome_name=exons$seqid[1],transcript_start=min(startend),transcript_end=max(startend),strand=ifelse(exons$strand[1]=="+",1,-1))
      #cat(paste0("chr",exons$seqid[1]),"\t",min(startend),max(startend),transcriptinfo$strand,"\n")##debug
      lseqtemp <-data.frame(transcript_exon_intron=as.character(toString(getSeq(genome,names=paste0("chr",exons$seqid[1]),start=min(startend),end=max(startend),strand=exons$strand[1]))),ensembl_transcript_id=trID)
      lseqtemp$transcript_exon_intron<-as.character(lseqtemp$transcript_exon_intron)
      lseqtemp$ensembl_transcript_id<-as.character(lseqtemp$ensembl_transcript_id)
    } else if (IN$usehg19!=-9) {#v0.92
      trID<-as.character(transciptID)
      preexons<-hg19exons[hg19exons$transcript_id==trID,]
      #if (nrow(preexons)<2) {  #since makedesion in v0.94 or later dump this kind of transcript, do not dump it here
      #  cat("Transcript:",trID,"has",nrow(preexons),"exon, which means no possibility for an aberrant splice site candidate.\n")
      #  return(data.frame())
      #}
      exons<-data.frame()
      currentn<-1
      for (exonnum in 1:nrow(preexons)) {
        if (nrow(preexons[preexons$exon_number==currentn,])>0) {
          exons<-rbind(exons,preexons[preexons$exon_number==currentn,])
          currentn<-currentn+1
        }
      }
      trstartend<-c(hg19transcripts$cdstart[hg19transcripts$transcript_id==trID][1],hg19transcripts$cdend[hg19transcripts$transcript_id==trID][1])
      if (trstartend[1]==-9) {
        trstartend<-c(hg19transcripts$start[hg19transcripts$transcript_id==trID][1],hg19transcripts$end[hg19transcripts$transcript_id==trID][1]) #in case of non-protein coding
        cat("Non-protein coding transcript:",trID, "found. Then cDNA position in the transcript means nDNA position.\n")
      }
      exonsinfo<-data.frame(ensembl_transcript_id=rep(trID,nrow(exons)),ensembl_exon_id=exons$exon_id,rank=exons$exon_number,chromosome_name=rep(paste0("chr",hg19transcripts$chr[hg19transcripts$transcript_id==trID][1]),nrow(exons)),exon_chrom_start=exons$start,exon_chrom_end=exons$end)
      transcriptinfo<-data.frame(hgnc_symbol=hg19transcripts$gene_name[hg19transcripts$transcript_id==trID][1],ensembl_transcript_id=trID,chromosome_name=hg19transcripts$chr[hg19transcripts$transcript_id==trID][1],transcript_start=hg19transcripts$start[hg19transcripts$transcript_id==trID][1],transcript_end=hg19transcripts$end[hg19transcripts$transcript_id==trID][1],strand=ifelse(hg19transcripts$strand[hg19transcripts$transcript_id==trID][1]=="+",1,-1))
      lseqtemp <-data.frame(transcript_exon_intron=as.character(toString(getSeq(genome,names=exonsinfo$chr[1],start=hg19transcripts$start[hg19transcripts$transcript_id==trID][1],end=hg19transcripts$end[hg19transcripts$transcript_id==trID][1],strand=hg19transcripts$strand[hg19transcripts$transcript_id==trID][1]))),ensembl_transcript_id=trID)
      lseqtemp$transcript_exon_intron<-as.character(lseqtemp$transcript_exon_intron)
      lseqtemp$ensembl_transcript_id<-as.character(lseqtemp$ensembl_transcript_id)
      
    } else {
      transcriptinfo<-getBM(attributes=c("hgnc_symbol","ensembl_transcript_id","chromosome_name","transcript_start","transcript_end","strand"),
                        transcriptkind, values=transciptID, mart=mart)
      if (nrow(transcriptinfo)==0 ) return ( data.frame() )
      if(nrow(transcriptinfo)>1 ) {
        cat("Multiple transcript found ",transcriptinfo$ensembl_transcript_id,"from given ID:",transciptID, "\tThen take 1st one...\n")
        transcriptinfo<-transcriptinfo[1,]
        transcriptkind<-"ensembl_transcript_id"
        transciptID<-transcriptinfo$ensembl_transcript_id
      }
      exonsinfo<-getBM(attributes=c("ensembl_transcript_id","ensembl_exon_id","rank","chromosome_name","exon_chrom_start","exon_chrom_end"),
                     transcriptkind, values=transciptID, mart=mart)
      utr5info<-getBM(attributes=c("5_utr_start","5_utr_end"),transcriptkind, values=transciptID, mart=mart)
      utr3info<-getBM(attributes=c("3_utr_start","3_utr_end"),transcriptkind, values=transciptID, mart=mart)
      utr5info<-utr5info[!is.na(utr5info[,1]),]
      utr3info<-utr3info[!is.na(utr3info[,1]),]
      utrinfo<-data.frame(start=c(utr5info$"5_utr_start",utr3info$"3_utr_start"),end=c(utr5info$"5_utr_end",utr3info$"3_utr_end"))
      lseqtemp     <- getSequence(id = transciptID, type=transcriptkind, mart =mart, seqType = "transcript_exon_intron") 
    }
  if (nrow(lseqtemp)>0) {
    lseqlength<-nchar(lseqtemp$transcript_exon_intron)
  } else {
    lseqlength<-0
    cat("TranscriptID:", transciptID,"was NOT found in the sequence database.\n")
  }
  if (transcriptinfo$strand==-1) {  #for reverse strand genes
    cat("Reverse strand gene found.\n")
    reverseflag<<-1 #universal variable to represent the gene is reversestrand genes
    temp<-transcriptinfo$transcript_start
    temp2<-transcriptinfo$transcript_end
    transcriptinfo$transcript_start<-temp2
    transcriptinfo$transcript_end<-temp
  } else {
    reverseflag<<-0
  }
  #cat("Lseqlength:",lseqlength,"phypos start,end:",transcriptinfo$transcript_start,",",transcriptinfo$transcript_end," transcript length by start end:",length(transcriptinfo$transcript_start:transcriptinfo$transcript_end),"\n")
  if (lseqlength==length(transcriptinfo$transcript_start:transcriptinfo$transcript_end)) {# check the consistency between annotated length (by start and end position) and actual sequence}
    lseq<-data.frame(GENE=rep(transcriptinfo$hgnc_symbol,lseqlength),ID=rep(transciptID,lseqlength),chr=rep(transcriptinfo$chromosome_name,lseqlength),phypos=transcriptinfo$transcript_start:transcriptinfo$transcript_end,cDNApos=rep(-9,lseqlength),exon=rep(-9,lseqlength),Seq=I(unlist(strsplit(lseqtemp$transcript_exon_intron,""))))
    for (i in 1:nrow(exonsinfo)) {  #mark exon number
      start<-which(lseq$phypos==exonsinfo[i,]$exon_chrom_start)
      end<-which(lseq$phypos==exonsinfo[i,]$exon_chrom_end)
      lseq[start:end,]$exon<-i
    }
    if (cDNAsw) {
      cdna.frame<-data.frame()
      for (i in 1:nrow(exonsinfo)) {  #integrate exon info and utr info then define cDNA region
        start<-exonsinfo[i,]$exon_chrom_start
        end<-exonsinfo[i,]$exon_chrom_end
        all<-0
        if (IN$useGTF19!=-9 | IN$usehg19!=-9) {
          #print(exonsinfo)#debug
          #cat("trstartend:",trstartend,"\n")#debug
          if (start<trstartend[1] & trstartend[1]<=end) start<-trstartend[1]
          if (start<=trstartend[2] & trstartend[2]<end) end<-trstartend[2]
          if (end< trstartend[1] | trstartend[2]<start) all<-1 #in case of the entire exon not included in translation
        } else {
        if (nrow(utrinfo)>0) for (j in 1:nrow(utrinfo)) {
          if (start<utrinfo[j,]$end & utrinfo[j,]$end<end) start<-utrinfo[j,]$end+1
          if (start<utrinfo[j,]$start & utrinfo[j,]$start<end) end<-utrinfo[j,]$start-1
          if (utrinfo[j,]$start<=start & end<=utrinfo[j,]$end) all<-1
        }
        }  
        if (all==0) {
          new<-data.frame(start=start,end=end,exon=i)
          #cat("New region start:",start," end:",end,"exon:",i)
          cdna.frame<-rbind(cdna.frame,new)
        }  
      } 
      if (transcriptinfo$strand==-1) {  #for reverse strand genes
        cdna.frame<-data.frame(start=cdna.frame$end,end=cdna.frame$start,exon=cdna.frame$exon)
      }
      cDNApos<-1
      #write.table(lseq,file="debug.lseq.txt",col.names=T,row.names=F,quote=F,sep="\t")#debug
      #write.table(cdna.frame,file="debug.cdna.frame.txt",col.names=T,row.names=F,quote=F,sep="\t")#debug
      for (i in 1:nrow(cdna.frame)) { #actual cDNA mark
        start<-which(lseq$phypos==cdna.frame[i,]$start)
        end<-which(lseq$phypos==cdna.frame[i,]$end)
        lseq[start:end,]$cDNApos<-cDNApos:(cDNApos+abs(end-start))
        cDNApos<-cDNApos+abs(end-start)+1
        new<-data.frame(ID=transciptID,cDNApos=paste("c.",cDNApos-1,"_",cDNApos,sep=""),exon=paste("5ss_exon",cdna.frame[i,]$exon,"_",cdna.frame[i,]$exon+1))
        exon.frame<<-rbind(exon.frame,new)  #record original 5ss between exons
      } 
    }
  } else {  #if inconsistency between  between annotated length (by start and end position) and actual sequence was found 
    cat("The length based on the sequence:",lseqlength,"\t The start and end position:",transcriptinfo$transcript_start,",",transcriptinfo$transcript_end,"\t The length based on the start and end positions:",length(transcriptinfo$transcript_start:transcriptinfo$transcript_end),"\nInconsistency was found. Then an empty sequence returned to avoid further error.\n")
    lseq<-data.frame()
  }
  } else {
    lseq<-masterlseq
  }
  if (IN$localseq!=-9) {
    files<-list.files(path=IN$localseq,pattern=transciptID)
    if (!(transciptID%in%files)) {
      cat("NOT found Seq data for transcipt:",transciptID,"on the folder:",IN$localseq, "So write it on the foleder",paste0(IN$localseq,"/",transciptID),"\n")
      write.table(lseq,file=paste0(IN$localseq,"/",transciptID),sep="\t",quote=F,col.names=T,row.names=F)
    }
  }
  return(lseq)
}

findtranscript <- function (chr,pos,mutname) {
  #targetval<-paste(mut.frame[i,]$chr,":",mut.frame[i,]$pos,":",mut.frame[i,]$pos,sep="")
  targetval<-paste(chr,":",pos,":",pos,sep="")
  cat ("Searching Transcripts for",mutname, "\tSubmitting string:",targetval,"\t" )
  if (IN$useGTF19!=-9) { #v0.91
    geneids<-gt$gene_id[gt$seqid==as.character(chr) & gt$start<=as.numeric(pos) & as.numeric(pos)<=gt$end]
    cat("geneid",geneids,"\n")
    transcripts<-c()
    
    for (geneid in geneids) {
      gm<-try(geneModel(ens, geneid))  #in case of the gene has just one exon, geneModel says that error! at least 2 exons needed. try means to avoid script stopping because of this
      if (class(gm)=="try-error") next
      #print (gm)
      #ll<-extractByGeneId(ens, geneid)#oldver
      #print (ll)#oldver
      pretranscripts <-attr(gm,"transcripts")#names(tableTranscript.id(ll))#oldver
      #cat("Candidate Transcripts:",pretranscripts,"\n")

      for (trID in pretranscripts) {
        tr<-getTranscript(gm,trID)
        startend<-attr(tr,"stcodon")
        cat("Candidate Transcript:",trID,"startend:",startend," ")
        start<-min(startend);end<-max(startend)
        if (!is.na(start) & !is.na(end)) if (start<=pos & pos<=end) {transcripts<-c(transcripts,trID);cat("Match.")}
        cat("\n")
      }
    }
    transcripts<-data.frame(ensembl_transcript_id=transcripts)
    transcripts$ensembl_transcript_id<-as.character(transcripts$ensembl_transcript_id)
  } else if (IN$usehg19!=-9){ #v0.92
    transcripts<-hg19transcripts[hg19transcripts$chr==as.character(chr) & hg19transcripts$start<=as.numeric(pos) &  as.numeric(pos)<=hg19transcripts$end,]
    if (nrow(transcripts)>0) {
      if (IN$takeCanonicalTranscript!=-9) { #take canonical transcript 1,protein-coding longest-> 2.protein non-coding longest #v0.93zeta
        cat("Try to find Canonical Transcript:")
        temptranscripts<-transcripts[transcripts$cdstart>-9,] #choose just protein coding
        if (nrow(temptranscripts)>0) { #in case of protein coding found
            cat("Protein Coding Transcript Found ")
            transcripts<-temptranscripts #just take protein coding
        } else {
            cat("Non-Protein Coding Transcript Found ")
        }
        transcripts$length<-apply(transcripts,1,function(x){
          if (as.numeric(x[7])>-9) return(as.numeric(x[8])-as.numeric(x[7])) else return(as.numeric(x[6])-as.numeric(x[5]))
        })
        transcripts<-transcripts[transcripts$length==max(transcripts$length),]
        transcripts<-transcripts[1,]
      }
      transcripts<-data.frame(ensembl_transcript_id=transcripts$transcript_id)
      transcripts$ensembl_transcript_id<-as.character(transcripts$ensembl_transcript_id)
    } else {
      transcripts<-data.frame()
    }
    
  } else { #biomaRt
    transcripts<-getBM(attributes=c("ensembl_transcript_id","transcript_start","transcript_end"),
                     "chromosomal_region", values=targetval, mart=mart)
    transcripts<-transcripts[transcripts$transcript_start <= pos & pos <= transcripts$transcript_end,]
  }
  if (nrow(transcripts)>0) {
    cat("Found Transcripts:",transcripts$ensembl_transcript_id,"\n")
    return(transcripts$ensembl_transcript_id)
  } 
  cat("No trascripts for",mutname,"\n")
  return(c())
}

chrposmut<-function(muttemp){  #build a matrix for mutation e.g. chr1:11233A>G (TranscriptID)
  mut.frame<-data.frame()
  for (i in 1:length(muttemp)) {
    muttemp[i]<-gsub("^(\\s+)|(\\s+)$","",muttemp[i])
    pretemp<-unlist(strsplit(muttemp[i]," "))
    if (length(pretemp)>1) transID<-pretemp[2] else transID<-NA
    temp<-unlist(strsplit(pretemp[1],":"))
    chr<-temp[1]
    chr<-sub("chr","",chr)
    pos<-gsub("(A|T|G|C|>)","",temp[2])
    alleles <-sub("\\d+","",temp[2])
    temp3<-unlist(strsplit(alleles,">"))
    ref<-temp3[1]
    alt<-temp3[2]
    #skipflag<-0  #In ver0.9 chr1:112233_A>AAG or indels are ignored  #In ver1.12 indels in vcf file came to be accepted
    cat ("Reading Line #",i," :",muttemp[i])
    #if (nchar(ref)>1 | nchar(alt)>1) {skipflag<-1;cat("\tSkipped\n");}  
    #if (skipflag==0) {
    if (nchar(ref)==1 & nchar(alt)==1 & length(grep("^(A|T|G|C|.)+$",ref))>0 & length(grep("^(A|T|G|C|.)+$",alt))>0) { #in the current ver, just SNV accepted
      cat("\tAccepted\n")
      if (is.na(transID)) {
        if (IN$esbTranscriptID!=-9) { #when transID NOT given in Mut file and esbTranscriptID DOES given in -esbTranscriptID switch
          transID<-IN$esbTranscriptID
          cat("TranscriptID given in -esbTranscriptID:",IN$esbTranscriptID,"applied\n")
        } else {
          transID<-findtranscript(chr,pos,muttemp[i])
        }
      } 
      tlen<-length(transID)
      new.frame<-data.frame()
      if (tlen>0) new.frame<-data.frame(kind=I(rep("chrpos",tlen)),chr=I(rep(chr,tlen)),pos=I(rep(as.integer(pos),tlen)),ref=I(rep(ref,tlen)),alt=I(rep(alt,tlen)),name=I(rep(paste("chr",chr,":",pos,ref,">",alt,sep=""),tlen)),transID=I(transID))
      if (nrow(mut.frame)==0) mut.frame<-new.frame else mut.frame<-rbind(mut.frame,new.frame)
      #cat("mut.frame chr:",chr,"pos:",pos,"ref:",ref,"alt:",alt,"\n")
    } else {
      cat("\tSkipped\n")
    }
  }
  if (nrow(mut.frame)>0) mut.frame<-mut.frame[order(mut.frame$transID),]  #because it takes time when transID differs from the previous one
  return(mut.frame)  
}

cDNAmut<-function(muttemp){ #build a matrix for mutation e.g. c.123A>G
  mut.frame<-data.frame()
  for (i in 1:length(muttemp)) {
    pos<-NA
    ref<-NA
    alt<-NA
    kind<-NA
    skipflag<-1  #In ver0.6 c.-14-51G>A c41+2T>A c53-11_53-7delCTT or ins are ignored
    pretemp<-unlist(strsplit(muttemp[i]," "))#v0.93 changed
    if (length(pretemp)>1) transID<-pretemp[2] else transID<-IN$esbTranscriptID#v0.93 changed
    temp<-regexpr("((\\-)*\\d+((\\+|\\-)\\d+)*)(_(\\-)*\\d+((\\+|\\-)\\d+)*)*",muttemp[i])
    pos<-substr(muttemp[i],as.numeric(temp),as.numeric(temp)+attr(temp,"match.length")-1)  #extract pos information like (-)123(+|-123)(_(-)345(+|-345))
    temp<-regexpr("(A|T|G|C)>(A|T|G|C)",muttemp[i])
    tempins<-regexpr("ins(A|T|G|C)+",muttemp[i])
    tempdel<-regexpr("del",muttemp[i])
    tempdelins<-regexpr("delins(A|T|G|C)+",muttemp[i])
    tempdup<-regexpr("dup",muttemp[i])
    tempdelins2<-regexpr("del(A|T|G|C)+(\\s+)*ins(A|T|G|C)+",muttemp[i])
    if (as.numeric(temp)!=-1) {   # G>A format
      ref<-substr(muttemp[i],as.numeric(temp),as.numeric(temp)) 
      alt<-substr(muttemp[i],as.numeric(temp)+2,as.numeric(temp)+2) 
      skipflag<-0
      kind<-"coding"
    } else if (as.numeric(tempdelins)!=-1 & attr(tempdelins,"match.length")>6) { #delinsATGC format
      ref<-"delins"
      alt<-substr(muttemp[i],as.numeric(tempdelins)+6,as.numeric(tempdelins)+attr(tempdelins,"match.length")-1)
      skipflag<-0
      kind<-"coding_delins"
    } else if (as.numeric(tempdelins2)!=-1 ) { #delATGCinsATGC format
      ref<-"del,ins"
      del<-NA
      ins<-NA
      temp<-regexpr("del(A|T|G|C)+",muttemp[i])  # take del seq
      if (as.numeric(temp)!=-1 & attr(temp,"match.length")>3 ) del<-substr(muttemp[i],as.numeric(temp)+3,as.numeric(temp)+attr(temp,"match.length")-1)
      temp<-regexpr("ins(A|T|G|C)+",muttemp[i])  #if dup seq given, take it
      if (as.numeric(temp)!=-1 & attr(temp,"match.length")>3 ) ins<-substr(muttemp[i],as.numeric(temp)+3,as.numeric(temp)+attr(temp,"match.length")-1)
      if (!is.na(del) & !is.na(ins)) {
        alt<-paste(del,ins,sep=",")
        skipflag<-0
        kind<-"coding_del_ins"
      }
    } else if (as.numeric(tempins)!=-1 & attr(tempins,"match.length")>3) {  #insATGC format
      ref<-"ins"  
      alt<-substr(muttemp[i],as.numeric(tempins)+3,as.numeric(tempins)+attr(tempins,"match.length")-1)
      skipflag<-0
      kind<-"coding_ins"
    } else if (as.numeric(tempdel)!=-1) {  #del(ATGC) format
      ref<-"del"
      alt<-""
      temp<-regexpr("del(A|T|G|C)+",muttemp[i])  #if del seq given, take it
      if (as.numeric(temp)!=-1 & attr(temp,"match.length")>3 ) alt<-substr(muttemp[i],as.numeric(temp)+3,as.numeric(temp)+attr(temp,"match.length")-1)
      skipflag<-0
      kind<-"coding_del"
    } else if (as.numeric(tempdup)!=-1) {  #dup(ATGC) format
      ref<-"dup"
      alt<-""
      temp<-regexpr("dup(A|T|G|C)+",muttemp[i])  #if dup seq given, take it
      if (as.numeric(temp)!=-1 & attr(temp,"match.length")>3 ) alt<-substr(muttemp[i],as.numeric(temp)+3,as.numeric(temp)+attr(temp,"match.length")-1)
      temp<-regexpr("dup\\d+",muttemp[i])  #if #dup given, take it
      if (as.numeric(temp)!=-1 & attr(temp,"match.length")>3 ) alt<-substr(muttemp[i],as.numeric(temp)+3,as.numeric(temp)+attr(temp,"match.length")-1)
      skipflag<-0
      kind<-"coding_dup"
    } 
    if (skipflag==0) {
      #      new.frame<-data.frame(kind=I("chrpos"),chr=I(chr),pos=as.integer(pos),ref=I(ref),alt=I(alt),name=I(paste("chr",chr,":",pos,ref,">",alt,sep="")),transID=I(transID))
      
      
      new.frame<-data.frame(kind=I(kind),chr=I(NA),pos=I(pos),ref=I(ref),alt=I(alt),name=I(paste("c.",pos,ref,ifelse(length(grep(">",muttemp[i]))>0,">",""),alt,sep="")),transID=I(transID))#v0.93 changed
      if (nrow(mut.frame)==0) mut.frame<-new.frame else mut.frame<-rbind(mut.frame,new.frame)
    } 
  }
  return(mut.frame)  
}

cDNAposdetermine<-function(pos) {
  pos<-as.character(pos)
  if (length(grep("_",pos))>0) c<-as.character(unlist(strsplit(pos,"_"))) else c<-as.character(pos)
  if (length(c)>1) {
    if (length(grep("\\D",c[1]))==0 & length(grep("\\D",c[2]))==0) {
      c1temp<-as.numeric(c[1])
      c2temp<-as.numeric(c[2])
      if ((c1temp>c2temp) & (c2temp  == (c1temp %% 10^nchar(c[2]) +1 ) ) )  {  #if given c.123_4insAA then complement second position
        cat ("Digits could be omitted in second position of ",pos, ". Therefore it was complemented automatically.\n")
        c[2]<-as.character(c1temp+1) 
      }
    }
    newpos<-c(cDNApossub(c[1]),cDNApossub(c[2])) 
  }
  else newpos<-c(cDNApossub(c[1]),cDNApossub(c[1]))
  
  #cat("pos1:",newpos[1]," pos2:", newpos[2],"\n")
  return (newpos)
}

cDNApossub<-function(pos) {
  pos<-as.character(pos)
  minusflag<-0 #if you have c.-12, then 1
  jointflag<-0 #if you have 12+12 ,then 1 if you have 12-12 then 2   
  if (length(grep("^-",pos))>0) {minusflag<-1; pos <-sub("-","",pos);}
  if (length(grep("\\-",pos))>0) {c<-unlist(strsplit(pos,"\\-"));jointflag<-2;} else if (length(grep("\\+",pos))>0) {c<-unlist(strsplit(pos,"\\+"));jointflag<-1;} else c<-pos
  if (minusflag) {
    temp<-which(lseq$cDNApos==1)
    newpos<-temp-as.numeric(c[1])
  } else newpos<-which(lseq$cDNApos==as.numeric(c[1]))
  if (jointflag>0) {
    if (jointflag==1) newpos<-newpos+as.numeric(c[2]) else newpos<-newpos-as.numeric(c[2])
  }
  return(newpos)
}

filltheblank<-function (refseq,altseq,word="#",forward=FALSE){  #this function returns filled sentence with WORD if the length between refseq and altseq is different . e.g.  ref=AA  alt=AATTTG then  return AA####  (forward=FALSE)  ####AA (forward=TRUE)     
  returnword<-""
  diff<-0
  if (nchar(refseq) > nchar(altseq)) {
    returnword<-altseq
    diff<-nchar(refseq) - nchar(altseq)
  } else {
    returnword<-refseq
    diff<-nchar(altseq) - nchar(refseq)
  }
  addword<-""
  if (diff>0) for (i in 1:diff) {addword<-paste(addword,word,sep="");}
  returnword<- ifelse (forward,paste(addword, returnword,sep=""),paste(returnword,addword,sep="") )
  return(returnword)
}

reversenuc<-function(text) {  #reverse nucleotide and order of sequence   ex)input ATGC  -> output GCAT
  cseq<-c()
  returntext<-""
  for (i in 1:nchar(text)) {
    c<-substr(text,i,i)
    if (c=="t") newc<-"a"
    if (c=="a") newc<-"t"
    if (c=="c") newc<-"g"
    if (c=="g") newc<-"c"
    if (c==".") newc<-"."
    if (c=="T") newc<-"A"
    if (c=="A") newc<-"T"
    if (c=="C") newc<-"G"
    if (c=="G") newc<-"C"
    if (c==".") newc<-"."
    if (c=="(") newc<-")"
    if (c==")") newc<-"("
    if (c=="[") newc<-"]"
    if (c=="]") newc<-"["
    cseq<-c(cseq,newc)
  }
  for (i in nchar(text):1) {
    returntext<-paste(returntext,cseq[i],sep="")
  }
  return(returntext)
}

assignseq<-function(addstr,lseq,chr,startpos) {
  justbeforeintron<-0
  pos<-nrow(lseq)
  if (pos==0) {
    cDNApos<-1
    exon<-0
    justbeforeintron<-1
  } else {
    cDNApos<-max(lseq$cDNApos)+1; if (cDNApos==(-9+1)) cDNApos<-1
    exon<-max(lseq$exon); if(exon==-9) exon<-0
    if (lseq$cDNApos[pos]==-9) justbeforeintron<-1
  }
  pos<-pos+1
  phypos<-pos+startpos-1
  if (IN$seqinverse==1) phypos<-startpos+length(lseqtemp2)-pos
  if (length(grep("a|t|g|c|n",addstr))>0) { #intron part
    addstr<-toupper(addstr)
    new<-data.frame(ID="Givenseq",Seq=I(addstr),chr=chr,phypos=I(phypos),cDNApos=-9,exon=-9)
  } else if (length(grep("A|T|G|C|N",addstr))>0) { #exon part
    if (justbeforeintron) exon<-exon+1
    new<-data.frame(ID="Givenseq",Seq=I(addstr),chr=chr,phypos=I(phypos),cDNApos=cDNApos,exon=exon)
  }
  return(new)
}

minmaxfunc<-function(phypos) {
    min<-min(phypos)
    max<-max(phypos)
  if (reverseflag) {
    minmax<-list(min=max,max=min)
  } else {
    minmax<-list(min=min,max=max)
  }
 return (minmax)
}

directory<-function(file) {
  return(paste(IN$refdirectory,file,sep="/"))
}

takeresult5<-function(dat,mode) {  #mode 0-> new5ss   1->broken canonical5ss
  #      new<-data.frame(phypos=newseq$phypos[j],MaxEntScoreSS5=output$Maxent,SS5score=output$SSscore,SS5cDNAboundary=cDNAboundary,SS5phyposboundary=phyposboundary,SS5seq=output$SSseq)
  #write.table(refresult,file="temp.refresult.txt",quote=F,sep="\t",col.names=T,row.names=F)
  #write.table(altresult,file="temp.altresult.txt",quote=F,sep="\t",col.names=T,row.names=F)
  #ref5<-dat$Ref_MaxEntScoreSS5
  #alt5<-dat$Alt_MaxEntScoreSS5
  #ref3<-dat$Ref_MaxEntScoreSS3
  #alt3<-dat$Alt_MaxEntScoreSS3
  ss5output<-0
  dat<-dat[dat$Alt_MaxEntScoreSS5>-999,] #delete null line
  if (mode==0) {
    if (length(grep("c.\\d+--9",as.character(dat$Ref_SS5cDNAboundary)))>0) dat<-dat[-grep("c.\\d+--9",as.character(dat$Ref_SS5cDNAboundary)),]#exclude canonical SS 
    pos<-which(dat$Alt_MaxEntScoreSS5==max(dat$Alt_MaxEntScoreSS5)) [1]
    #if (length(grep("c.\\d+--9",as.character(dat$Ref_SS5cDNAboundary[pos])))==0) { #in case of creating 5ss, 1st pass
    cat("Peak Alt MaxEnt Score Pos:",pos," Alt score:",dat$Alt_MaxEntScoreSS5[pos], "Ref score:",dat$Ref_MaxEntScoreSS5[pos],"\n")
    if ((dat$Ref_MaxEntScoreSS5[pos]-dat$Alt_MaxEntScoreSS5[pos])<0) {#in case of creating 5ss, determine
      ss5output<-1
    } else {
      cat("No Splice Site Creating Peak found in MaxEnt 5SS Score.\n")
      ss5output<-1 # all result will be reported in spite of negative results
    } 
    #}
  } else if (mode==1) {
    if (length(grep("c.\\d+--9",dat$Ref_SS5cDNAboundary))>0) { #in case of canonical exon-intron boundary involved, 
      pos<-grep("c.\\d+--9",dat$Ref_SS5cDNAboundary)[1]
      #write.table(refresult,file="temp.debug.txt",row.names=F,col.names=T,sep="\t",quote=F)
      cat("Exon-Intron Boundary MaxEnt Score Pos:",pos," Alt score:",dat$Alt_MaxEntScoreSS5[pos], "Ref score:",dat$Ref_MaxEntScoreSS5[pos],"\n")
      if ((dat$Ref_MaxEntScoreSS5[pos]-dat$Alt_MaxEntScoreSS5[pos])>0) {#in case of broken 5ss, determine
        ss5output<-1
      } else {
        cat("The Score in the Exon-Intron boundary does not meet the canonical 5SS loss criteria.\n")
      }
    } else {
      cat("No Exon-Intron boundary found in the region MaxEnt 5SS Scanned.\n")
    }
  }
  if (ss5output) {
    return(dat[pos,])
  } else {
    return(data.frame())
  }
}

takeresult3<-function(dat,mode) {  #mode 0-> new5ss   1->broken canonical5ss
  #      new<-data.frame(phypos=newseq$phypos[j],MaxEntScoreSS5=output$Maxent,SS5score=output$SSscore,SS5cDNAboundary=cDNAboundary,SS5phyposboundary=phyposboundary,SS5seq=output$SSseq)
  #write.table(refresult,file="temp.refresult.txt",quote=F,sep="\t",col.names=T,row.names=F)
  #write.table(altresult,file="temp.altresult.txt",quote=F,sep="\t",col.names=T,row.names=F)
  #data.frame(phypos=newseq$phypos[j],MaxEntScoreSS3=output$Maxent,SS3score=output$SSscore,SS3cDNAboundary=cDNAboundary,SS3phyposboundary=phyposboundary,SS3seq=output$SSseq)
  #ref5<-dat$Ref_MaxEntScoreSS5
  #alt5<-dat$Alt_MaxEntScoreSS5
  #ref3<-dat$Ref_MaxEntScoreSS3
  #alt3<-dat$Alt_MaxEntScoreSS3
  dat<-dat[dat$Alt_MaxEntScoreSS3>-999,] #delete null line
  ss3output<-0
  if (mode==0) {
    #if (nrow(altresult[altresult$MaxEntScoreSS5>4,])>0) { #in case of creating 5ss, 1st pass
    if (length(grep("c.-9-\\d+",as.character(dat$Ref_SS3cDNAboundary)))>0) dat<-dat[-grep("c.-9-\\d+",as.character(dat$Ref_SS3cDNAboundary)),]#exclude canonical SS 
    pos<-which(dat$Alt_MaxEntScoreSS3==max(dat$Alt_MaxEntScoreSS3)) [1]
    #if (length(grep("c.-9-\\d+",as.character(Ref_SS3cDNAboundary[pos])))==0) { 
    cat("Peak Alt MaxEnt Score Pos in the Seq:",pos," Alt score:",dat$Alt_MaxEntScoreSS3[pos], "Ref score:",dat$Ref_MaxEntScoreSS3[pos],"\n")
    if ((dat$Ref_MaxEntScoreSS3[pos]-dat$Alt_MaxEntScoreSS3[pos])<0) {#in case of creating 3ss, determine
      ss3output<-1
    } else {
      cat("No Splice Site Creating Peak found in MaxEnt 3SS Score.\n")
      ss3output<-1 # all result will be reported in spite of negative results
    } 
    #} else {
    #  cat("Peak Position is on an Intron-Exon Boundary:",as.character(Ref_SS3cDNAboundary[pos]),"\n")
    #}
  } else if (mode==1) {
    if (length(grep("c.-9-\\d+",as.character(dat$Ref_SS3cDNAboundary)))>0) { #in case of canonical exon-intron boundary involved, 
      pos<-grep("c.-9-\\d+",as.character(dat$Ref_SS3cDNAboundary))[1]
      #write.table(refresult,file="temp.debug.txt",row.names=F,col.names=T,sep="\t",quote=F)
      cat("Intron-Exon Boundary:",as.character(dat$Ref_SS3cDNAboundary[pos]) ,"MaxEnt Score Pos in the Seq:",pos," Alt score:",dat$Alt_MaxEntScoreSS3[pos], "Ref score:",dat$Ref_MaxEntScoreSS3[pos],"\n")
      if ((dat$Ref_MaxEntScoreSS3[pos]-dat$Alt_MaxEntScoreSS3[pos])>0) {#in case of broken 3ss, determine
        ss3output<-1
      } else {
        cat("The Score in the Intron-Exon boundary does not meet the canonical 3SS loss criteria.\n")
        ss3output<-1 # all result will be reported in spite of negative results
      }
    } else {
      cat("No Intron-Exon boundary found in the region MaxEnt 3SS Scanned.\n")
    }
  }
  if (ss3output) {
    return(data.frame(dat[pos,]))
  } else {
    return(data.frame())
  }
}

makedecision<-function(dat) { #this function attaches Probability and Decision columns to result frame
  ref5<-dat$Ref_MaxEntScoreSS5
  alt5<-dat$Alt_MaxEntScoreSS5
  ref3<-dat$Ref_MaxEntScoreSS3
  alt3<-dat$Alt_MaxEntScoreSS3
  
  #cat("nrow",nrow(dat),"names",names(dat),"ref5",ref5,"alt5",alt5,"ref3",ref3,"alt3",alt3,"\n")
  dat$Possibility<-"No"
  dat$Decision<-"No"
  if (!("Exon"%in%names(dat))) {
    Var_Name<-NA
    if (length(grep("chr",dat$Var_name))>0) Var_Name<-dat$Var_name
    if ("VarName_ChrPos"%in%names(dat))  Var_Name<-dat$VarName_ChrPos
    temp<-exonnumber(Var_Name,dat$ID)
    dat$Exon<-as.character(temp$exon)
    dat$Max.Exon<-as.character(temp$maxexon)
    dat$Protein.Coding<-as.character(temp$proteincoding)
  }
  

  
  if (dat$Kind=="5SS_Gain" & (alt5-ref5)>0) {
    if (alt5>0) dat$Possibility<-"Low"
    if (alt5>IN$ss5gainlow) dat$Possibility<-"Mid"
    if (alt5>IN$ss5gainhigh) dat$Possibility<-"High"
    if (alt5>IN$ss5gaintake) dat$Decision<-"YES"
    if (!is.na(dat$Exon)) {
      if (as.character(dat$Exon)==as.character(dat$Max.Exon)) {dat$Possibility<-paste0(dat$Possibility,"_On_the_Last_Exon");dat$Decision<-"No";}
      if (as.character(dat$Exon)=="UTR5") {dat$Possibility<-paste0(dat$Possibility,"_On_the_UTR5");dat$Decision<-"No";}
      if (as.character(dat$Exon)=="UTR3") {dat$Possibility<-paste0(dat$Possibility,"_On_the_UTR3");dat$Decision<-"No";}
    }
  } else if (dat$Kind=="5SS_Loss"){
    if ((ref5-alt5)>0) dat$Possibility<-"Low"
    if ((ref5-alt5)>IN$ss5losslow) dat$Possibility<-"Mid"
    if ((ref5-alt5)>IN$ss5losshigh) dat$Possibility<-"High"
    if ((ref5-alt5)>IN$ss5losstake) dat$Decision<-"YES"
    #if (dat$Ref_SS3cDNAboundary=="c.-9-1") {dat$Possibility<-paste0(dat$Possibility,"_For_the_Acceptor_of_1st_Exon");dat$Decision<-"No";}
  } else if (dat$Kind=="3SS_Gain" & (alt3-ref3)>0){
    if (alt3>0) dat$Possibility<-"Low"
    if (alt3>IN$ss3gainlow) dat$Possibility<-"Mid"
    if (alt3>IN$ss3gainhigh) dat$Possibility<-"High"
    if (alt3>IN$ss3gaintake) dat$Decision<-"YES"
    if (!is.na(dat$Exon)) {
      if (as.character(dat$Exon)=="1") {dat$Possibility<-paste0(dat$Possibility,"_On_the_1st_Exon");dat$Decision<-"No";}
      if (as.character(dat$Exon)=="UTR5") {dat$Possibility<-paste0(dat$Possibility,"_On_the_UTR5");dat$Decision<-"No";}
      if (as.character(dat$Exon)=="UTR3") {dat$Possibility<-paste0(dat$Possibility,"_On_the_UTR3");dat$Decision<-"No";}
    }
    
  } else if (dat$Kind=="3SS_Loss"){
    if ((ref3-alt3)>0) dat$Possibility<-"Low"
    if ((ref3-alt3)>IN$ss3losslow) dat$Possibility<-"Mid"
    if ((ref3-alt3)>IN$ss3losshigh) dat$Possibility<-"High"
    if ((ref3-alt3)>IN$ss3losstake) dat$Decision<-"YES"
    if (dat$Ref_SS3cDNAboundary=="c.-9-1") {dat$Possibility<-paste0(dat$Possibility,"_For_the_Acceptor_of_1st_Exon");dat$Decision<-"No";}
  }
  
  if (!is.na(dat$Exon)) {
    if (as.numeric(dat$Max.Exon)==1) {dat$Possibility<-paste0(dat$Possibility,"_Transcript_With_Just_One_Exon");dat$Decision<-"No";}
  }
  dat<-dat[,-which(names(dat)=="SID")]
  return(dat)
}


TakeSigLines<-function(resultfilename,outputprefix) {
  dat<-read.table(resultfilename,sep="\t",header=T, stringsAsFactors = F)
  #dat$Kind<-NA
  #dat$Probability<-NA
  #dat$Decison<-NA
  if (IN$skip1stMiningStep==-9) {
    dat$SID<-apply(dat,1,function(x){paste0(x["Var_name"],"_",x["ID"])})
    SIDs<-unique(dat$SID)
    
    #big loop start
    cat("\nSummaizing Process Starts...\n")
    result.frame<-data.frame()
    for (ID in SIDs) {
      cat("Processing unique ID:",ID,"\n")
      unit<-dat[dat$SID==ID,] # take one group according to ID (Varname+transcriptID)
      result<-takeresult5(unit,1);   if (nrow(result)>0) {result$Kind<-"5SS_Loss"; result.frame<-rbind(result.frame,makedecision(result))}
      result<-takeresult3(unit,1);   if (nrow(result)>0) {result$Kind<-"3SS_Loss"; result.frame<-rbind(result.frame,makedecision(result))}
      result<-takeresult5(unit,0);   if (nrow(result)>0) {result$Kind<-"5SS_Gain"; result.frame<-rbind(result.frame,makedecision(result))}
      result<-takeresult3(unit,0);   if (nrow(result)>0) {result$Kind<-"3SS_Gain"; result.frame<-rbind(result.frame,makedecision(result))}
    }
    
    #finish
    cat("Writing a result to:",paste0(outputprefix,".SigLines.txt"),"\n")
    write.table(result.frame,file=paste0(outputprefix,".SigLines.txt"),sep="\t",col.names=T,row.names=F,quote=F)
    if (IN$sjdbout!=-9) {
      cat("Checking Splicing gap for ",paste0(outputprefix,".SigLines.txt"),"\n")
      sjdb(paste0(outputprefix,".SigLines.txt"),paste0(outputprefix,".SigLines.txt"))
      result.frame<-read.table(paste0(outputprefix,".SigLines.txt"),sep="\t",header=T,stringsAsFactors = F)
    }
    
  } else {
    cat("1st procedure was skipped.\n")
    result.frame<-dat
  }
  
  #integrate results with diffrent transcript ID for a variant
  if (IN$SummarizeByVariant==1 | IN$TakeOnePerVariant==1) {
    if (IN$skip2ndMiningStep==-9) {
      sum.result.frame<-data.frame()
      dat<-result.frame
      dat$SID<-apply(dat,1,function(x){paste0(x["Var_name"],"_",x["Kind"],"_",x["Ref_SS5phyposboundary"])})
      SIDs<-unique(dat$SID)
      
      cat("Summarizing lines by variant...\n")
      #big loop start
      for (ID in SIDs) {
        cat("Processing unique ID:",ID,"\n")
        unit<-dat[dat$SID==ID,] # take one group according to ID (Varname+transcriptID)
        
        #priotize the result of Loss  <-since the 1st procedure focus on priotizing loss factors, no more need for these lines.
        #if (nrow(unit)>1 & length(grep("Loss",unit$Kind))>0) {
        #  tempunit<-data.frame();use<-c()
        #  if (length(which(unit$Kind=="5SS_Loss"))>0) {pos<-which(unit$Kind=="5SS_Loss");for (j in pos){tempunit<-rbind(tempunit,unit[pos,]);use<-c(use,pos)}}
        #  if (length(which(unit$Kind=="3SS_Loss"))>0) {pos<-which(unit$Kind=="3SS_Loss");for (j in pos){tempunit<-rbind(tempunit,unit[pos,]);use<-c(use,pos)}}
        #  for (j in 1:nrow(unit)) {
        #    if (!(j%in%use)) tempunit<-rbind(tempunit,unit[j,])
        #  }
        #  unit<-tempunit
        #}
        
        newunit<-unit[1,] 
        if (nrow(unit)>1) for (n_row in 2:nrow(unit)) {
          for (n_col in 1:ncol(unit)) {
            if (is.na(newunit[1,n_col])) newunit[1,n_col]<-"NA"
            if (is.na(unit[n_row,n_col])) unit[n_row,n_col]<-"NA"
            if (newunit[1,n_col]!=unit[n_row,n_col]) newunit[1,n_col]<-paste0(newunit[1,n_col],",",unit[n_row,n_col])
          }
        }
        sum.result.frame<-rbind(sum.result.frame,newunit)
      }
      
      sum.result.frame<-sum.result.frame[,-which(names(sum.result.frame)=="SID")]
      
      if (IN$SummarizeByVariant==1) {
        cat("Writing a summarized result by variant to:",paste0(outputprefix,".SigLinesByVariant.txt"),"\n")
        write.table(sum.result.frame,file=paste0(outputprefix,".SigLinesByVariant.txt"),sep="\t",col.names=T,row.names=F,quote=F)
      }
    } else {
      cat("Summarize By Variant procedure was skipped.\n")
      sum.result.frame<-dat
    }
    if (IN$TakeOnePerVariant==1) {
      dat<-sum.result.frame
      sum.result.frame<-data.frame()
      SIDs<-unique(dat$Var_name)
      cat("Making lines where one line is assigned to one variant...\n")
      #big loop start
      for (ID in SIDs) {
        cat("Processing Var_name:",ID,"\n")
        preunit<-dat[dat$Var_name==ID,] # take one group according to ID (Varname+transcriptID)
        unit<-preunit[preunit$Possibility!="No",]
        if (nrow(unit)==0) {unit<-preunit} #in case of all no significant, just take all; if you do not want to take this, please delete this line.
        if (nrow(unit)>0) {
          newunit<-unit[1,]
          if (nrow(unit)>1) {
            for (n_row in 2:nrow(unit)) {
              for (n_col in 1:ncol(unit)) {
                if (is.na(newunit[1,n_col])) newunit[1,n_col]<-"NA"
                if (is.na(unit[n_row,n_col])) unit[n_row,n_col]<-"NA"
                if (newunit[1,n_col]!=unit[n_row,n_col]) newunit[1,n_col]<-paste0(newunit[1,n_col],",",unit[n_row,n_col])
              }
            }
          }
          sum.result.frame<-rbind(sum.result.frame,newunit)
        } 
      }
      
      cat("Writing a summarized result with 1 line for 1 variant to:",paste0(outputprefix,".OneLineByVariant.txt"),"\n")
      write.table(sum.result.frame,file=paste0(outputprefix,".OneLineByVariant.txt"),sep="\t",col.names=T,row.names=F,quote=F)
    }
  }
  
}

sjdb<-function(resultfilename,outputfilename) {
  lines<-read.table(resultfilename,sep="\t",header=T,stringsAsFactors = F)
  lines$Predicted_SpliceJunction<-NA
  lines$Predicted_NormalJunction<-NA #v0.96
  lines$Predicted_Effect<-NA
  sjdbFile<-c() #for sjdbFile
  
  if (!("Kind"%in%names(lines))) {
    cat("The column <Kind> not found in the file:",resultfilename,". Then splice junction information is not able to be attached.\n")
    return()
  }
  
  #mainloop
  for (i in 1:nrow(lines)) {
    line<-lines[i,]
    cat("Processing #line",i,"Var:",line$Var_name,"\t")
    IDs<-unlist(strsplit(line$ID,","))
    Gap<-unlist(strsplit(ifelse (length(grep("5SS",line$Kind))>0,unlist(strsplit(line$Ref_SS5phyposboundary,",")),unlist(strsplit(line$Ref_SS3phyposboundary,","))),":|-"))
    #cat("Gap",Gap)
    
    if (IN$sjdbrestrictin!=-9) {
      spl<-unlist(strsplit(IN$sjdbrestrictin,","))
      colpos<-which(names(line)==spl[1])
      if (length(colpos)>0) {
        if (line[1,colpos]!=spl[2]) {cat("Skip line by restrictin value:",line[1,colpos],"required:",spl[2],"\n");next;}
      } else {
        cat("No Column found for restrict in\t")
      }
    }
    if (IN$sjdbrestrictout!=-9) {
      spl<-unlist(strsplit(IN$sjdbrestrictout,","))
      colpos<-which(names(line)==spl[1])
      if (length(colpos)>0) {
        if (line[1,colpos]==spl[2]) {cat("Skip line by restrictout value:",line[1,colpos],"required:",spl[2],"\n");next;}
      } else {
        cat("No Column found for restrict out\t")
      }
    }
    for (nID in 1:length(IDs)) {
      
      ID<-IDs[nID]
      #lseq<-GetSeq(ID,1)
      #strand<-ifelse(lseq$phypos[1]<lseq$phypos[2],"+","-")
      if (length(grep("5SS",line$Kind))>0) {  #bug-fix part
        intpos<-as.numeric(Gap[length(Gap)]) #for 5SS 1st-intron is in later position regardress of the strand +/-
      } else {
        intpos<-as.numeric(Gap[length(Gap)-1]) #for 3SS last-intron is formert position regardless of the strand +/-
      }
      #cat("intpos",intpos)
      #pos<-which(lseq$phypos==intpos) #intron start pos
      
      #cat("pos",pos)
      #pos2<-pos #intron end pos
      #search splice acceptor
      intpos<-as.numeric(intpos)
      canonicalintpos<--9 #canonical splice donor or acceptor start site
      intpos2<--9 #if not found , set -9
      if (IN$usehg19!=-9) {
        trID<-as.character(ID)
        preexons<-hg19exons[hg19exons$transcript_id==trID,]
        if (nrow(preexons)<2) {
          cat("Transcript:",trID,"has",nrow(preexons),"exon, which means no possibility for an aberrant splice site candidate.\n")
          next
        }
        exons<-data.frame()
        currentn<-1
        for (exonnum in 1:nrow(preexons)) {
          if (nrow(preexons[preexons$exon_number==currentn,])>0) {
            exons<-rbind(exons,preexons[preexons$exon_number==currentn,])
            currentn<-currentn+1
          }
        }
        exoninfo<-data.frame(ensembl_exon_id=exons$exon_id,rank=exons$exon_number,exon_chrom_start=exons$start,exon_chrom_end=exons$end,strand=ifelse(hg19transcripts$strand[hg19transcripts$transcript_id==trID][1]=="+",1,-1))
      } else { #use biomaRt
        exoninfo<-getBM(attributes=c("ensembl_exon_id","rank","exon_chrom_start","exon_chrom_end","strand"),"ensembl_transcript_id", ID, mart)
      }
      #if (trID=="ENST00000225538") write.table(exoninfo,file="temp.exons.ENST00000225538.info.txt",sep="\t",col.names=T,row.names=F,quote=F)#debug
      strand<-ifelse(exoninfo$strand[1]>0,"+","-")
      if (exoninfo$strand[1]>0) { #forward strand
        if (length(grep("5SS",line$Kind))>0) {  
          edge<-sort(exoninfo$exon_chrom_start[exoninfo$exon_chrom_start>intpos])
          if (length(edge)>0) {intpos2<-edge[1]-1} else {cat("Corresponding exon edge not found because of the Variant located on the edge\n");next;}
          canonicaledge<-sort(exoninfo$exon_chrom_end[exoninfo$exon_chrom_end<intpos2])
          canonicalintpos<-canonicaledge[length(canonicaledge)]+1
        } else {
          edge<-sort(exoninfo$exon_chrom_end[exoninfo$exon_chrom_end<intpos])
          if (length(edge)>0) {intpos2<-edge[length(edge)]+1} else {cat("Corresponding exon edge not found because of the Variant located on the edge\n");next;}
          canonicaledge<-sort(exoninfo$exon_chrom_start[exoninfo$exon_chrom_start>intpos2])
          canonicalintpos<-canonicaledge[1]-1
        }
        
      } else { #reverse strand
        if (length(grep("5SS",line$Kind))>0) {  
          edge<-sort(exoninfo$exon_chrom_end[exoninfo$exon_chrom_end<intpos])
          if (length(edge)>0) {intpos2<-edge[length(edge)]+1} else {cat("Corresponding exon edge not found because of the Variant located on the edge\n");next;}
          canonicaledge<-sort(exoninfo$exon_chrom_start[exoninfo$exon_chrom_start>intpos2])
          canonicalintpos<-canonicaledge[1]-1
        } else {
          edge<-sort(exoninfo$exon_chrom_start[exoninfo$exon_chrom_start>intpos])
          if (length(edge)>0) {intpos2<-edge[1]-1}  else {cat("Corresponding exon edge not found because of the Variant located on the edge\n");next;}
          canonicaledge<-sort(exoninfo$exon_chrom_end[exoninfo$exon_chrom_end<intpos2])
          canonicalintpos<-canonicaledge[length(canonicaledge)]+1
          
        }
        
      }
      if (intpos2>0) {# if not found, intpos2 is -9. then pass the output
        if (intpos<intpos2) sj<-paste0(intpos,"_",intpos2) else sj<-paste0(intpos2,"_",intpos)
        newgap<-paste0("chr",line$chr,"_",sj,"_",strand)
        if (length(grep(newgap,lines$Predicted_SpliceJunction[i]))==0) {
          lines$Predicted_SpliceJunction[i]<-paste0(lines$Predicted_SpliceJunction[i],",",newgap)
          sjdbFile<-c(sjdbFile,gsub("_","\t",newgap))
        }
        #cat (paste0("chr",line$chr,"\t",gsub("_","\t",sj),"\t",strand),"\t") #for debug
        if (canonicalintpos>0) {
          canonicalgap<-abs(canonicalintpos-intpos2)+1
          aberrantgap<-abs(intpos-intpos2)+1
          diff<-aberrantgap-canonicalgap
          if (diff==0) {
            lines$Predicted_Effect[i]<-paste0("Canonical_Splice_Site:",diff,"bp")
            exonstartend<-sort(c(exoninfo$exon_chrom_start,exoninfo$exon_chrom_end))
            normalgap<-NA
            if (length(exonstartend[exonstartend>max(canonicalintpos,intpos2)])>2) {
              remain<-exonstartend[exonstartend>max(canonicalintpos,intpos2)]
              normalgap<-paste0("chr",line$chr,"_",(remain[2]+1),"_",(remain[3]-1),"_",strand)
            } else if (length(exonstartend[exonstartend<min(canonicalintpos,intpos2)])>2) {
              remain<-exonstartend[exonstartend<min(canonicalintpos,intpos2)]
              max<-length(remain)
              normalgap<-paste0("chr",line$chr,"_",(remain[max-2]+1),"_",(remain[max-1]-1),"_",strand)
            }
            #when canonical splice site, nearby splice junction provided
            if (!is.na(normalgap)){
              if (length(grep(normalgap,lines$Predicted_NormalJunction[i]))==0) lines$Predicted_NormalJunction[i]<-paste0(lines$Predicted_NormalJunction[i],",",normalgap)
              sjdbFile<-c(sjdbFile,gsub("_","\t",normalgap))
            }
          } else {
            prefix<-ifelse(abs(diff)%%3==0,"in-frame","frame-shift")
            postfix<-ifelse(diff>0,"_extension:","_deletion:")
            lines$Predicted_Effect[i]<-paste0(prefix,postfix,diff,"bp")
            if (canonicalintpos<intpos2) normalsj<-paste0(canonicalintpos,"_",intpos2) else normalsj<-paste0(intpos2,"_",canonicalintpos)
            normalgap<-paste0("chr",line$chr,"_",normalsj,"_",strand)
            if (length(grep(normalgap,lines$Predicted_NormalJunction[i]))==0) lines$Predicted_NormalJunction[i]<-paste0(lines$Predicted_NormalJunction[i],",",normalgap)
          }
        }
      }
    }
    lines$Predicted_SpliceJunction[i]<-sub("NA,","",lines$Predicted_SpliceJunction[i])
    lines$Predicted_NormalJunction[i]<-sub("NA,","",lines$Predicted_NormalJunction[i])
    cat("\n")
  }#mainloop end
  cat("Writing a result with sjGap information to:",IN$output,"\n")
  write.table(lines,file=outputfilename,sep="\t",col.names=T,row.names=F,quote=F)
  cat("Writing a sjdbFile.tab to:",IN$sjdbout,"\n")
  sjdbFile<-unique(sjdbFile)
  writeLines(sjdbFile,con=IN$sjdbout,sep="\n")
}


#v0.94subs
#cdnatophypos(mut.frame[i,],lseq)
cdnatophypos<-function(mut.frame,lseq) { #in the current ver, just SNV can be converted into phypos
  str<-unlist(strsplit(mut.frame$pos,"_"))[1] #in the current ver, second position is ignored. e.g. c.1158-28_1161del,  1161del will be ignored.
  if (length(grep("+|-",str))>0) {
    str2<-unlist(strsplit(str,ifelse(length(grep("+",str))>0,"+","-")))
    temp<-regexpr("\\d+",str2[1]) 
    basepos<-as.numeric(substr(str2[1],as.numeric(temp),as.numeric(temp)+attr(temp,"match.length")-1))
    temp<-regexpr("\\d+",str2[2]) 
    addnum<-as.numeric(substr(str2[2],as.numeric(temp),as.numeric(temp)+attr(temp,"match.length")-1))
  } else {
    temp<-regexpr("\\d+",str) 
    basepos<-as.numeric(substr(str,as.numeric(temp),as.numeric(temp)+attr(temp,"match.length")-1))
    addnum<-0
  }
  pos<- lseq$phypos[which(lseq$cDNApos==basepos)]
  ref<-mut.frame$ref
  alt<-mut.frame$alt
  if (reverseflag==1) {
    ref<-reversenuc(ref)
    alt<-reversenuc(alt)
    if (length(grep("+",str))>0) {
      pos=pos-addnum; # - reverse strand   + forward strand
    } else if (length(grep("-",str))>0) {
      pos=pos+addnum; # + reverse strand   - forward strand
    }
  } else {
    if (length(grep("+",str))>0) {
      pos=pos+addnum; # - reverse strand   + forward strand
    } else if (length(grep("-",str))>0) {
      pos=pos-addnum; # + reverse strand   - forward strand
    }
  }
  phypos<-paste0("chr",lseq$chr[1],":",pos,ref,">",alt)
  return(phypos)
}

#phypostocdna(mut.frame[i,],lseq)
phypostocdna<-function(mut.frame,lseq) {
  basepos<-as.numeric(mut.frame$pos)
  ref<-mut.frame$ref
  alt<-mut.frame$alt
  cdnapos<-NA
  #cat ("chrpos:",chrpos,"basepos:",basepos,"\n")
  if (lseq$cDNApos[which(lseq$phypos==basepos)]!=-9) {
    cdnapos<-lseq$cDNApos[which(lseq$phypos==basepos)]
  } else {
    step<-0
    flag<-0
    while (flag==0) {
      step<-step+1
      if (length(which(lseq$phypos==(basepos-step)))>0) if (lseq$cDNApos[which(lseq$phypos==(basepos-step))]!=-9) flag<-1 
      if (length(which(lseq$phypos==(basepos+step)))>0) if (lseq$cDNApos[which(lseq$phypos==(basepos+step))]!=-9) flag<-2 
    }
    add1<-"+";add2<-"-"
    if (reverseflag==1) {add1<-"-";add2<-"+";}
    if (flag==1) cdnapos<-paste0(lseq$cDNApos[which(lseq$phypos==(basepos-step))],add1,step) #because this is for reverse strand gene
    if (flag==2) cdnapos<-paste0(lseq$cDNApos[which(lseq$phypos==(basepos+step))],add2,step) #because this is for reverse strand gene
  }
  if (reverseflag==1) {ref<-reversenuc(ref);alt<-reversenuc(alt);}
  cdnapos<-paste0("c.",cdnapos,ref,">",alt)
  return(cdnapos)
}

exonnumber<-function(phypos,trID) { #in the current ver, work when -usehg19 mode
  str<-unlist(strsplit(phypos,":"))[2]  #phypos should be chr1:123456A>G or something
  temp<-regexpr("\\d+",str) 
  basepos<-as.numeric(substr(str,as.numeric(temp),as.numeric(temp)+attr(temp,"match.length")-1)) #just take physical position
  onexonnum<-NA
  maxexonnum<-NA
  protein.coding<-NA
  if (IN$usehg19!=-9) {
    protein.coding<-"Yes"
    preexons<-hg19exons[hg19exons$transcript_id==trID,]
    exons<-data.frame()
    currentn<-1
    for (exonnum in 1:nrow(preexons)) {
      if (nrow(preexons[preexons$exon_number==currentn,])>0) {
        addexon<-preexons[preexons$exon_number==currentn,]
        addexon<-addexon[1,]
        exons<-rbind(exons,addexon)
        currentn<-currentn+1
      }
    }
    for (exonnum in 1:nrow(exons)) {
      if (exons$start[exonnum]<=basepos & basepos<=exons$end[exonnum]) { 
        onexonnum<-exons$exon_number[exonnum] 
      } else if (exonnum<nrow(exons)) {
        if (hg19transcripts$strand[hg19transcripts$transcript_id==trID][1]=="+" & exons$end[exonnum]<basepos & basepos<exons$start[exonnum+1]) onexonnum<-paste0("I",exons$exon_number[exonnum],"-",exons$exon_number[exonnum]+1)
        if (hg19transcripts$strand[hg19transcripts$transcript_id==trID][1]=="-" & exons$end[exonnum+1]<basepos & basepos<exons$start[exonnum]) onexonnum<-paste0("I",exons$exon_number[exonnum],"-",exons$exon_number[exonnum]+1)
      }
    }
    if (hg19transcripts$strand[hg19transcripts$transcript_id==trID][1]=="+") {
      if (basepos<exons$start[1]) onexonnum<-"UTR5"
      if (basepos>exons$end[nrow(exons)]) onexonnum<-"UTR3"
    } else {
      if (basepos<exons$start[nrow(exons)]) onexonnum<-"UTR3"
      if (basepos>exons$end[1]) onexonnum<-"UTR5"
    }
    maxexonnum<-as.character(nrow(exons))
    if (hg19transcripts$cdstart[hg19transcripts$transcript_id==trID][1]==-9) {
      protein.coding<-"No"
    }
  }
  return(data.frame(exon=onexonnum,maxexon=maxexonnum,proteincoding=protein.coding))
}

readvcf<-function(vcffile) { # if vcf file given 
  cat("Reading vcf file...\n")
  muttemp<-readLines(con=IN$vcffile)  #interpret mutfile
  head<-c()
  newline<-c()
  for (i in 1:length(muttemp)) {
    if (length(grep("^##",muttemp[i]))>0) next;
    if (length(grep("^#",muttemp[i]))>0) {   #take a header information 
      head<-unlist(strsplit(muttemp[i],"\\s+"))
      head[1]<-sub("#","",head[1])
      chr<-which(head=="CHROM")
      pos<-which(head=="POS")
      ref<-which(head=="REF")
      alt<-which(head=="ALT")
      filt<-which(head=="FILTER")
      info<-which(head=="INFO")
      next
    }
    spl<-unlist(strsplit(muttemp[i],"\\s+"))
    
    #if (spl[filt]=="PASS") {  #just take variants with filter pass
      temp<-regexpr("TRANSCRIPT_ID=ENST\\d+;",spl[info])
      transID<-c()
      while (regexpr("TRANSCRIPT_ID=ENST\\d+;",spl[info])>0 ) {
        temp<-regexpr("TRANSCRIPT_ID=ENST\\d+;",spl[info])
        tempID<-substr(spl[info],as.numeric(temp),as.numeric(temp)+attr(temp,"match.length")-1)
        spl[info]<-sub(tempID,"",spl[info])
        tempID<-gsub("TRANSCRIPT_ID=|;","",tempID)
        transID<-c(transID,tempID)
      }
      transID<-unique(transID)
      if (length(grep("^chr",spl[chr]))>0) spl[chr]<-sub("chr","",spl[chr])
      alts<-unlist(strsplit(spl[alt],","))
      if (length(transID)>0) {
        for (trID in transID) {
          for (altallele in alts) {
            newline<-c(newline,paste0("chr",spl[chr],":",spl[pos],spl[ref],">",altallele," ",trID))
          }
        }
      } else {
        for (altallele in alts) {
          newline<-c(newline,paste0("chr",spl[chr],":",spl[pos],spl[ref],">",altallele))
        }
      }
    #}
  }
  return(newline)
}



#main
cat("This script is to calculate the integrated score of 5ss or 3ss\n")

#Input parametes -mutfile (chrpos format mut file with transcript ID) -output 
inputlist<-list(helpflag="-help",debug="-debug",marthost="-marthost",mutfile="-mutfile",output="-output",seq="-seq",seqinverse="-seqinverse",nomart="-nomart",
                esbTranscriptID="-esbTranscriptID",canonicalSS="-canonicalSS",entiregene="-entiregene",ss5="-ss5",ss3="-ss3",wide="-wide",lseq="-lseq",
                skipRegressScore="-skipRegressScore",skipSRE="-skipSRE",localseq="-localseq",refdirectory="-refdirectory",useGTF19="-useGTF19",usehg19="-usehg19"
                ,ss5gainlow="-ss5gainlow",ss5gainhigh="-ss5gainhigh",ss3gainlow="-ss3gainlow",ss3gainhigh="-ss3gainhigh"
                ,ss5losslow="-ss5losslow",ss5losshigh="-ss5losshigh",ss3losslow="-ss3losslow",ss3losshigh="-ss3losshigh"
                ,ss5gaintake="-ss5gaintake",ss3gaintake="-ss3gaintake",ss5losstake="-ss5losstake",ss3losstake="-ss3losstake",SummarizeByVariant="-SummarizeByVariant"
                ,TakeOnePerVariant="-TakeOnePerVariant",skip1stMiningStep="-skip1stMiningStep",skip2ndMiningStep="-skip2ndMiningStep"
                ,summarizeresult="-summarizeresult",summarizeresultin="-summarizeresultin",
                sjdbout="-sjdbout",sjdbrestrictin="-sjdbrestrictin",sjdbrestrictout="-sjdbrestrictout",sjdbin="-sjdbin",takeCanonicalTranscript="-takeCanonicalTranscript"
                ,skipcDNApos="-skipcDNApos",usehg38="-usehg38",vcffile="-vcffile")
inputparas<-c(0,0,1,1,1,1,0,0,
              1,0,0,0,0,1,1,
              0,0,1,1,1,1
              ,1,1,1,1
              ,1,1,1,1
              ,1,1,1,1,0,0,0,0
              ,1,1,
              1,1,1,1,0
              ,0,1,1)  #the length of parameters to be input in the command line 0-switch yes or no 1- one variable 2-two variables
IN<-inputparameters(inputlist,inputparas)
if (IN$helpflag==1) stop("\n",Helpdocument)
if (IN$output==-9) IN$output<-"Regress_Score.v0.93.R.result.txt"
if (IN$ss5==-9 & IN$ss3==-9) {IN$ss5<-1;IN$ss3<-1}
if (IN$refdirectory==-9) IN$refdirectory<-"./ssfiles"
if (IN$summarizeresult==-9) IN$summarizeresult=IN$summarizeresultin

#belows are theresholds for -summarizeresult
if (IN$ss5gainlow==-9) IN$ss5gainlow<-4.1   #low threshold for probability      (0) low possible        (low) possible    (high) probable
if (IN$ss5gainhigh==-9) IN$ss5gainhigh<-7.1 #high threshold for probability
if (IN$ss3gainlow==-9) IN$ss3gainlow<-7.3  #low threshold for probability
if (IN$ss3gainhigh==-9) IN$ss3gainhigh<-8.7 #high threshold for probability

if (IN$ss5losslow==-9) IN$ss5losslow<-0   #low threshold for probability      (0) low possible        (low) possible    (high) probable
if (IN$ss5losshigh==-9) IN$ss5losshigh<-1.5 #high threshold for probability
if (IN$ss3losslow==-9) IN$ss3losslow<-0  #low threshold for probability
if (IN$ss3losshigh==-9) IN$ss3losshigh<-1.5 #high threshold for probability

if (IN$ss5gaintake==-9) IN$ss5gaintake<-4.1 #threshold for taking the variant  (Decision)
if (IN$ss5losstake==-9) IN$ss5losstake<-1.5
if (IN$ss3gaintake==-9) IN$ss3gaintake<-7.3
if (IN$ss3losstake==-9) IN$ss3losstake<-1.5
##

cat("Parameters:\n")
for (i in names(IN)) {
  cat(i,":",unlist(IN[which(names(IN)==i)]),"\n")
}

#initiation for using biomaRT
#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
#biocLite("BiocUpgrade")
#to use developer ver. download source file, then type install.packages("biomaRt_2.23.0.tar.gz", repos = NULL, type="source")

masterlseq<-data.frame()
if (IN$useGTF19!=-9) {#v0.91 newly added
  IN$nomart<-1
  library(refGenome)
  ens<-ensemblGenome()
  cat("Using hg19 coordinated GTF file from ENSEMBL.\nReading GTF file from the directory",IN$useGTF19,"...\n")
  basedir(ens)<-IN$useGTF19
  read.gtf(ens,"Homo_sapiens.GRCh37.75.gtf")
  gt<-getGeneTable(ens)
  cat("Also BSgenome.Hsapiens.UCSC.hg19 is employed for hg19 sequence.\n")
  library(BSgenome.Hsapiens.UCSC.hg19)
  genome<-BSgenome.Hsapiens.UCSC.hg19
}
if (IN$usehg19!=-9) {#v0.91 newly added
  IN$nomart<-1
  cat("Using hg19 coordinated informatoin from ENSEMBL.\nReading files from the directory",IN$usehg19,"...\n")
  hg19transcripts<-read.table(paste0(IN$usehg19,"/hg19.transcripts.txt"),sep="\t",header=T,stringsAsFactors = F)
  hg19exons<-read.table(paste0(IN$usehg19,"/hg19.exons.txt"),sep="\t",header=T,stringsAsFactors = F)
  cat("Also BSgenome.Hsapiens.UCSC.hg19 is employed for hg19 sequence.\n")
  library(BSgenome.Hsapiens.UCSC.hg19)
  genome<-BSgenome.Hsapiens.UCSC.hg19
}
if (IN$usehg38!=-9) {#v0.91 newly added
  IN$nomart<-1
  cat("Using hg38 coordinated informatoin from ENSEMBL.\nReading files from the directory",IN$usehg38,"...\n")
  hg19transcripts<-read.table(paste0(IN$usehg38,"/hg38.transcripts.txt"),sep="\t",header=T,stringsAsFactors = F)
  hg19exons<-read.table(paste0(IN$usehg38,"/hg38.exons.txt"),sep="\t",header=T,stringsAsFactors = F)
  cat("Also BSgenome.Hsapiens.UCSC.hg38 is employed for hg38 sequence.\n")
  library(BSgenome.Hsapiens.UCSC.hg38)
  genome<-BSgenome.Hsapiens.UCSC.hg38
  IN$usehg19<-"DummyOn" #mimic hg19 mode
}
if (IN$lseq!=-9) {
  IN$nomart<-1
  masterlseq<-read.table(IN$lseq,header=T,sep="\t",stringsAsFactors=F)
  if (IN$esbTranscriptID==-9) IN$esbTranscriptID<-IN$lseq
}
if (IN$nomart==-9) {
  library(biomaRt) #if you need downloading biomaRt, please type as follows
  cat("Connecting biomaRt server...\n")
  marthostname<-ifelse(IN$marthost==-9,"feb2014.archive.ensembl.org",IN$marthost) #feb2014.archive.ensembl.org   is hg19 latest ver host
  cat("Host name is",marthostname,"\n")
  mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host=marthostname)  # biomart for gene annotation
  #host="may2009.archive.ensembl.org",  hg18
  #host="Sep2013.archive.ensembl.org", hg19
} else {
  cat("Skip connecting biomaRt server.\n")
}

if (IN$summarizeresultin!=-9) {
  TakeSigLines(IN$summarizeresultin,IN$summarizeresult)
  cat("Done.")
  quit("no")
}

if (IN$sjdbin!=-9) {
  sjdb(IN$sjdbin,IN$output)
  cat("Done.")
  quit("no")
}

#preparing 3ss precalculated matrix
cat("Preparing maxentscan3ss parameters...\n")
metables<-makemaxentscores3()

#preparing 5ss precalculated matrix
cat("Preparing maxentscan5ss parameters...\n")
scorematrix<-makemaxentscores5()


#Read ESE and ESS consensus sequence
cat("Reading ESE and ESS sequences...\n")
ESEseq<-readLines(con=directory("ESE.txt"))
ESElength<-nchar(ESEseq[1])

ESSseq<-readLines(con=directory("fas-hex3.txt"))
ESSlength<-nchar(ESSseq[1])
cat("Reading ISE sequences...\n")
ISEseq<-readLines(con=directory("FAS-ISE.txt"))
ISE9seq<-ISEseq[nchar(ISEseq)==9]
ISE9length<-9
ISE10seq<-ISEseq[nchar(ISEseq)==10]
ISE9length<-10

#set weight of each element and parameters in regression model
Eweight<-list(ESE5="47.84491,60.92733",ESE3="48.83209,60.70861",
              bGrun="44.53022,83.98231",aGrun="0.00000,58.65735",
              bISE="81.18222,14.99889",aISE="22,15,82,15",
              ESS5="42.17946,65.45010,203.29782,27.52051",ESS3="-63.3802067,15.25197,0.7984965,31.62696,86.4669410,59.74309,202.9456879,22.73068")
RegCoef5<-list(Intercept=-4.745344,MAXENT5=0.492312,ESE3=17.136712,ESS3=-20.178786,aGrun=80.972702,aISE=11.225050)
RegCoef3<-list(Intercept=-4.320291,MAXENT3=0.407358,ESE5=19.912353,ESS5=-31.136150,bGrun=61.247241,bISE=25.560157)

reverseflag<-0
muttemp<-c()
lseq<-data.frame()
newcDNApos<-c()
if (IN$seq!=-9) { #exon should be UPPERCASE, intron should be LOWERCASE, you can indicate mutation like (A|B), where A means refseq and B mean alt seq
                  #                                       Also you can indicate nucleutide of interest like ATGC[ATGC]ATGCa   
  cat("loading sequence information from:",IN$seq,"\n")
  data<-readLines(con=IN$seq)
  chr<-"-9"
  startpos<-1
  #clean up lines
  lseqtemp<-""
  if (length(data)>1) { 
    for (i in 1:length(data)) {
     if (length(grep(">",data[i]))) {
       temp<-regexpr("chr\\w+:\\d+",data[i]) #instead of variant position, actual 5ss gap taken
       if (as.numeric(temp)>0) {
          chrpos<-unlist(strsplit(substr(data[i],as.numeric(temp),as.numeric(temp)+attr(temp,"match.length")-1),":"))
          chr<-sub("chr","",chrpos[1])
          startpos<-as.numeric(chrpos[2])
       } 
     } else {
       lseqtemp<-paste(lseqtemp,data[i],sep="")
     }
    }
  } else {
   lseqtemp<-data[1]
  }
  ##clearnup the strings
  lseqtemp2<-c()
  if (IN$seqinverse==1) {start<-nchar(lseqtemp);end<-1} else {start<-1;end<-nchar(lseqtemp)} #if inverse switch on, then take characters inversely
  for (i in start:end) {  #construct string including necessary information ( and discard unnecessary characters)
    if (length(grep("a|t|g|c|n|A|T|G|C|N|(|)|[|]|\\|",substr(lseqtemp,i,i)))>0) {
      addstr<-substr(lseqtemp,i,i)
      if (IN$seqinverse==1) addstr<-reversenuc(addstr)
      lseqtemp2<-c(lseqtemp2,addstr)
    }
  }
  #construct sequence frame
  refposstart<-0
  refposend<-0
  concernstart<-0
  concernend<-0
  altseq<-""
  i<-1
  while (i<=length(lseqtemp2)) {
      addstr<-lseqtemp2[i]
      pos<-nrow(lseq)+1
      if (addstr=="(") {  #indicate mut like (ref|alt)
        refposstart<-pos
        i<-i+1
        tempchr<-lseqtemp2[i]
        loopflag<-1
        altcollect<-0
        while (loopflag) {
          if (tempchr==")") {
            loopflag<-0
          } else if (altcollect) {
            altseq<-paste(altseq,tempchr,sep="")
          } else if (tempchr=="|") {
            refposend<-nrow(lseq)
            altcollect<-1
          } else  {
            lseq<-rbind(lseq,assignseq(tempchr,lseq,chr,startpos))
          }
          if (loopflag) i<-i+1
          tempchr<-lseqtemp2[i]
        }
      } else if (addstr=="[") {
        concernstart<-pos
      } else if (addstr=="]") {
        concernend<-pos-1
      } else if (length(grep("a|t|g|c|n|A|T|G|C|N",addstr))>0) { 
        lseq<-rbind(lseq,assignseq(addstr,lseq,chr,startpos))
      }
      i<-i+1
  }
  if ((concernstart*concernend)!=0) {
    newcDNApos<-c(concernstart,concernend)
  } else if ((refposstart*refposend) !=0) {
    phypos<-refposstart+startpos-1
    if (IN$seqinverse==1) phypos<-startpos+length(lseqtemp2)-refposstart
    muttemp<-paste("chr",chr,":",phypos,paste(lseq$Seq[refposstart:refposend],collapse=""),">",altseq,sep="")
  }
}

if (IN$mutfile!=-9) {
  cat("Reading a mutation file:",IN$mutfile,"\n")
  muttemp<-readLines(con=IN$mutfile)
} else if (IN$vcffile!=-9) {
  cat("Reading a vcf file:",IN$vcffile,"\n")
  muttemp<-readvcf(IN$vcffile)
}

mut.frame<-data.frame()
if (length(muttemp)>0) {
  if (length(grep("c\\.",muttemp))>0) mut.frame<-cDNAmut(muttemp)
  if (length(grep("chr",muttemp))>0) mut.frame<-chrposmut(muttemp) #new.frame<-data.frame(kind=I("chrpos"),chr=I(chr),pos=as.integer(pos),ref=I(ref),alt=I(alt),name=I(paste("chr",chr,":",pos,ref,">",alt,sep="")),transID=I(transID))
}

exon.frame<-data.frame(ID=I("DUMMY"),cDNApos=I("DUMMY"),exon=I("DUMMY")) #record original 5ss between exons
result.frame<-data.frame()

if (IN$debug==1) {
  write.table(lseq,file=paste("debug.lseq.",basename(IN$output),sep=""),col.names=T,row.names=F,quote=F)
  write.table(mut.frame,file=paste("debug.mut.table.",basename(IN$output),sep=""),col.names=T,row.names=F,quote=F)
}

#mut.frame loop
preTransID<-"dummy"
if (nrow(mut.frame)>0) { 
  mut.frame$VarName_cDNAPos<-NA #v0.94
  mut.frame$VarName_ChrPos<-NA #v0.94
  mut.frame$Exon<-NA
  mut.frame$Max.Exon<-NA
  mut.frame$Protein.Coding<-NA
  for (i in 1:nrow(mut.frame)) {  #main loop for mut.frame
    cat("Processing",mut.frame$name[i],"   ",mut.frame$transID[i],"...\n")
  getseq<-0;if (!is.na(mut.frame$transID[i])) if (mut.frame$transID[i]!=preTransID) getseq<-1
  if (nrow(lseq)==0 | getseq) {lseq<-GetSeq(mut.frame$transID[i],1);preTransID<-mut.frame$transID[i]}
  if (length(grep("coding",mut.frame[i,]$kind))>0) {
    newcDNApos<-cDNAposdetermine(mut.frame[i,]$pos) #which(lseq$cDNApos==mut.frame[i,]$pos)
    mut.frame$VarName_cDNAPos[i]<-mut.frame$name[i]
    mut.frame$VarName_ChrPos[i]<-cdnatophypos(mut.frame[i,],lseq)
  }
if (mut.frame[i,]$kind=="chrpos" ) {
  newcDNApos<-which(lseq$phypos==mut.frame[i,]$pos);newcDNApos<-c(newcDNApos,newcDNApos)
  if (IN$skipcDNApos==-9) mut.frame$VarName_cDNAPos[i]<-phypostocdna(mut.frame[i,],lseq)
  mut.frame$VarName_ChrPos[i]<-mut.frame$name[i]
}
if (reverseflag==1 & mut.frame[i,]$kind=="chrpos") {
  cat("ChrPos information in a reverseStarnd gene. Therefore, change the allele info for cDNA.\n")
  mut.frame[i,]$ref<-reversenuc(mut.frame[i,]$ref)
  mut.frame[i,]$alt<-reversenuc(mut.frame[i,]$alt)
}
if (!is.na(mut.frame$VarName_ChrPos[i])) {
  temp<-exonnumber(as.character(mut.frame$VarName_ChrPos[i]),mut.frame$transID[i])
  mut.frame$Exon[i]<-as.character(temp$exon)
  mut.frame$Max.Exon[i]<-as.character(temp$maxexon)
  mut.frame$Protein.Coding[i]<-as.character(temp$proteincoding)
}
if (length(newcDNApos)>0) { #if cDNApos assigned correctly
newseq<-lseq
if (mut.frame[i,]$kind=="chrpos" & length(grep("^(A|T|G|C)$",mut.frame[i,]$ref))==1  & length(grep("^(A|T|G|C)$",mut.frame[i,]$alt))==1 | mut.frame[i,]$kind=="coding" ) {
  newseq[newcDNApos[1],]$Seq<-mut.frame[i,]$alt;indelflag<-0;
} else if (mut.frame[i,]$kind=="coding_ins") {
  insert<-data.frame()
  for (j in 1:nchar(mut.frame[i,]$alt)) {
    temp<-lseq[1,]
    temp$Seq<-substr(mut.frame[i,]$alt,j,j)
    temp$phypos<-paste(lseq[newcDNApos[1],]$phypos,"+ins",j,sep="")
    temp$cDNApos<-paste(lseq[newcDNApos[1],]$cDNApos,"+ins",j,sep="")
    temp$exon<-lseq[newcDNApos[1],]$exon
    insert<-rbind(insert,temp)
  }
  newseq<-rbind(newseq[1:newcDNApos[1],],insert,newseq[newcDNApos[2]:nrow(newseq),])
  newcDNApos[1]<-newcDNApos[1]+1
  newcDNApos[2]<-newcDNApos[2]+nrow(insert)-1
  
} else if (mut.frame[i,]$kind=="coding_del") {
  delseq<-""
  for (j in newcDNApos[1]:newcDNApos[2]) {delseq<-paste(delseq,newseq[j,]$Seq,sep="");}
  #cat("Deleted Seq is :",toupper(delseq),"mut.alt:",mut.frame[i,]$alt,"\n")####
  if (nchar(mut.frame[i,]$alt)==0 | is.na(mut.frame[i,]$alt)) cat("Deleted Seq is :",toupper(delseq),"\n") else if (toupper(delseq)!=toupper(mut.frame[i,]$alt)) cat("Warning. Given del seq:",toupper(mut.frame[i,]$alt),"is different from actual cDNA seq:",toupper(delseq),"\n")
  newseq<-rbind(newseq[1:(newcDNApos[1]-1),],newseq[(newcDNApos[2]+1):nrow(newseq),])
  #newcDNApos[1]<-newcDNApos[1]
  newcDNApos[2]<-newcDNApos[1]-1
} else if (mut.frame[i,]$kind=="coding_delins") {
  insert<-data.frame()
  for (j in 1:nchar(mut.frame[i,]$alt)) {
    temp<-lseq[1,]
    temp$Seq<-substr(mut.frame[i,]$alt,j,j)
    temp$phypos<-paste(lseq[newcDNApos[1],]$phypos,"+ins",i,sep="")
    temp$cDNApos<-paste(lseq[newcDNApos[1],]$cDNApos,"+ins",j,sep="")
    temp$exon<-lseq[newcDNApos[1],]$exon
    insert<-rbind(insert,temp)
  }
  newseq<-rbind(newseq[1:(newcDNApos[1]-1),],insert,newseq[(newcDNApos[2]+1):nrow(newseq),])
  newcDNApos[2]<-newcDNApos[1]+nrow(insert)-1
} else if (mut.frame[i,]$kind=="coding_del_ins") {
  delins<-unlist(strsplit(mut.frame[i,]$alt,","))
  delseq<-"";
  for (j in newcDNApos[1]:newcDNApos[2]) {delseq<-paste(delseq,newseq[j,]$Seq,sep="");}
  if (toupper(delseq)!=toupper(delins[1])) cat("Warning. Given del seq:",toupper(mut.frame[i,]$alt),"is different from actual cDNA seq:",toupper(delseq),"\n")
  insert<-data.frame()
  for (j in 1:nchar(mut.frame[i,]$alt)) {
    temp<-lseq[1,]
    temp$Seq<-substr(mut.frame[i,]$alt,j,j)
    temp$phypos<-paste(lseq[newcDNApos[1],]$phypos,"+ins",j,sep="")
    temp$cDNApos<-paste(lseq[newcDNApos[1],]$cDNApos,"+ins",j,sep="")
    temp$exon<-lseq[newcDNApos[1],]$exon
    insert<-rbind(insert,temp)
  }
  newseq<-rbind(newseq[1:(newcDNApos[1]-1),],insert,newseq[(newcDNApos[2]+1):nrow(newseq),])
  newcDNApos[2]<-newcDNApos[1]+nrow(insert)-1
} else if (mut.frame[i,]$kind=="coding_dup") {
  insertseq<-"";
  for (j in newcDNApos[1]:newcDNApos[2]) {
    insertseq<-paste(insertseq,toupper(lseq[j,]$Seq),sep="")
  }
  if (length(grep("\\d+",mut.frame[i,]$alt))>0) {
    insertseq2<-"";
    for (j in 1:(as.numeric(mut.frame[i,]$alt)-1)) {
      insertseq2<-paste(insertseq2,insertseq,sep="")
    }
  } else if (nchar(mut.frame[i,]$alt)==0) {
    insertseq2<-insertseq
  } else {
    insertseq2<-toupper(mut.frame[i,]$alt)
    #cat ("insertseq",insertseq,"insertseq2",insertseq2,"i",i,"ref",mut.frame[i,]$ref,"alt",class(mut.frame[i,]$alt),"\n")
    #warnings()
    if (insertseq!=substr(insertseq2,1,nchar(insertseq))) cat ("Warning. Given dup seq:",insertseq2,"is different from actual cDNA seq:",insertseq,"\n")
  }
  insert<-data.frame()
  for (j in 1:nchar(insertseq2)) {
    temp<-lseq[1,]
    temp$Seq<-substr(insertseq2,j,j)
    temp$phypos<-paste(lseq[newcDNApos[1],]$phypos,"+dup",j,sep="")
    temp$cDNApos<-paste(lseq[newcDNApos[2],]$cDNApos,"+dup",j,sep="")
    temp$exon<-lseq[newcDNApos[2],]$exon
    insert<-rbind(insert,temp)
  }
  newseq<-rbind(newseq[1:(newcDNApos[2]),],insert,newseq[(newcDNApos[2]+1):nrow(newseq),])
  newcDNApos[1]<-newcDNApos[2]+1
  newcDNApos[2]<-newcDNApos[2]+nrow(insert)
} else if (mut.frame[i,]$kind=="chrpos") { #in case of indel of vcf file
  delseq<-"";
  ncommon<-0;
  
  if (mut.frame[i,]$ref=="." ) mut.frame[i,]$ref<-""
  if (mut.frame[i,]$alt=="." ) mut.frame[i,]$alt<-""
  ref <-mut.frame[i,]$ref
  alt <-mut.frame[i,]$alt
  add<-ifelse(nchar(mut.frame[i,]$alt)>0,1,0)  #for length(alt)==0 
  if (reverseflag==1) {
    temprefalt<-filltheblank(ref,alt,word="#",forward=TRUE)
    if (nchar(ref)<nchar(alt)) ref<-temprefalt else if (nchar(ref)>nchar(alt)) alt<-temprefalt
    conflag<-1;
    newref<-""  # new ref seq where common charcters omitted
    newalt<-"" # new alt seq where common characters omitted
    for (j in nchar(ref):1) {
      if ( substr(ref,j,j)==substr(alt,j,j) & conflag==1 ) {ncommon<-ncommon+1} else {conflag<-0;newref<-paste(ifelse(substr(ref,j,j)!="#",substr(ref,j,j),""),newref,sep="");newalt<-paste(ifelse(substr(alt,j,j)!="#",substr(alt,j,j),""),newalt,sep="")}
    }
    cat("newref:",newref," newalt:",newalt,"\n")
    refstart<-newcDNApos[1]-nchar(mut.frame[i,]$ref) #just left position of substitution
    refend<-newcDNApos[1]-ncommon+1 #just right position of substitution
  } else {
    temprefalt<-filltheblank(ref,alt,word="#",forward=FALSE)
    if (nchar(ref)<nchar(alt)) ref<-temprefalt else if (nchar(ref)>nchar(alt)) alt<-temprefalt
    conflag<-1;
    newref<-""  # new ref seq where common charcters omitted
    newalt<-"" # new alt seq where common characters omitted
    for (j in 1:nchar(ref)) {
      if ( substr(ref,j,j)==substr(alt,j,j) & conflag==1 ) {ncommon<-ncommon+1} else {conflag<-0;newref<-paste(newref,ifelse(substr(ref,j,j)!="#",substr(ref,j,j),""),sep="");newalt<-paste(newalt,ifelse(substr(alt,j,j)!="#",substr(alt,j,j),""),sep="")}
    }
    cat("newref:",newref," newalt:",newalt,"\n")
    refstart<-newcDNApos[1]+ncommon-1 #just left position of substitution
    refend<-newcDNApos[1]+nchar(mut.frame[i,]$ref) #just right position of substitution
  }
  newcDNApos[1]<-refstart+1
  newcDNApos[2]<-refstart+nchar(newalt)
  
  for (j in (refstart+1-ifelse(reverseflag==0,ncommon,0)):(refend-1+ifelse(reverseflag==1,ncommon,0))) {delseq<-paste(delseq,newseq[j,]$Seq,sep="");}
  if (length(grep("^(A|T|G|C)+$",mut.frame[i,]$ref))==1 & mut.frame[i,]$ref == delseq | nchar(mut.frame[i,]$ref)==0) {# at 1st, check if ref seq is consistent with given ref information
    cat("Refseq:",mut.frame[i,]$ref ,"  Altseq:",mut.frame[i,]$alt,"\n")
    insert<-data.frame() #make insert seq  from alt seq
    subcount<-1 #substitution count
    insseq<-""
    if (nchar(newalt)>0) for (j in 1:nchar(newalt)) {
      temp<-lseq[1,]
      temp$Seq<-substr(newalt,j,j)
      temp$phypos<-paste(lseq[refstart,]$phypos,ifelse(reverseflag==1,"-","+"),"ins",j,sep="")
      tempcDNApos<-ifelse(lseq[refstart,]$cDNApos==-9,-9,lseq[refstart,]$cDNApos)
      temp$cDNApos<-paste(tempcDNApos,"+ins",j,sep="")
      temp$exon<-lseq[refstart,]$exon
      insert<-rbind(insert,temp)
    }
    newseq<-rbind(newseq[1:refstart,],insert,newseq[refend:nrow(newseq),])
    
  } else {
    cat("Warning. Given ref seq in vcf:",toupper(mut.frame[i,]$ref),"is different from actual cDNA seq:",toupper(delseq),"\nAnalysis continued, though.\n")
  }
}

#lseq<-data.frame(GENE=rep(transcriptinfo$hgnc_symbol,lseqlength),ID=rep(transciptID,lseqlength),chr=rep(transcriptinfo$chromosome_name,lseqlength),phypos=transcriptinfo$transcript_start:transcriptinfo$transcript_end,cDNApos=rep(-9,lseqlength),exon=rep(-9,lseqlength),Seq=I(unlist(strsplit(lseqtemp$transcript_exon_intron,""))))
#SSscan<-function(seq,ExPosInTheContext,SS)  {#seq a sequence including 5ss or 3ss, ExPosInTheContext exon position(Start,End) in the given seq  SS should be 5(5Splice site) or 3 (3Spice site) 
if (IN$debug==1) write.table(newseq,file=paste("debug.newseq.",basename(IN$output),sep=""),col.names=T,row.names=F,quote=F)

if (IN$wide!=-9) {
  cat("Scan region ", IN$wide, "bp attached around the indicated mutation area.\n")
  newcDNApos[1]<-newcDNApos[1]-as.numeric(IN$wide)
  newcDNApos[2]<-newcDNApos[2]+as.numeric(IN$wide)
  if (newcDNApos[1]<1) newcDNApos[1]<-1
  if (newcDNApos[2]>nrow(newseq)) newcDNApos[2]<-nrow(newseq)
}

#scan for 5ss
result5<-data.frame()
if (IN$ss5==1) {
  cat("Checking the potential of 5\'splice site...\n")
  Start<-newcDNApos[1]-9+1
  if (Start<1) Start<-1  #indicate the region where the mutation affects
  End <-newcDNApos[2]
  if (End > (nrow(newseq)-9+1)) End<-nrow(newseq)-9+1
  #cat("SS5 search pos1:",newcDNApos[1], " pos2:",newcDNApos[2], "Start:",Start," End", End, "\n")
  result5<-SS5scan(newseq,Start,End,newcDNApos)
  if (IN$debug==1) write.table(result5,file=paste("debug.result5.",basename(IN$output),sep=""),col.names=T,row.names=F,quote=F)
}


#scan for 3ss
result3<-data.frame()
if (IN$ss3==1) {
  cat("Checking the potential of 3\'splice site...\n")
  Start<-newcDNApos[1]-23+1
  if (Start<1) Start<-1  #indicate the region where the mutation affects
  End <-newcDNApos[2]
  if (End > (nrow(newseq)-23+1)) end<-nrow(newseq)-23+1
  #cat("SS3 search pos1:",newcDNApos[1], " pos2:",newcDNApos[2], "Start:",Start," End", End, "\n")
  result3<-SS3scan(newseq,Start,End,newcDNApos)
  if (IN$debug==1) write.table(result3,file=paste("debug.result3.",basename(IN$output),sep=""),col.names=T,row.names=F,quote=F)
}



if (IN$ss3==-9) {
  minmax<-minmaxfunc(result5$phypos)
} else if (IN$ss5==-9) {
  minmax<-minmaxfunc(result3$phypos)
} else {
  minmax<-minmaxfunc(c(result3$phypos,result5$phypos))
}
min<-minmax$min
max<-minmax$max

#-entiregene
if (IN$entiregene==1) {
  cat("Checking entire sequence of the gene. Then it could take time...\n")
  min<-lseq$phypos[1]
  max<-lseq$phypos[nrow(lseq)]
}

#calc scores in ref seq
cat("Calculating Scores in the referrence sequence...\n")

lseqmin<-which(lseq$phypos==min)
lseqmax<-which(lseq$phypos==max)
if (IN$ss5==1) {
  result5ref<-SS5scan(lseq,lseqmin,lseqmax,newcDNApos)
  if (IN$debug==1) write.table(result5ref,file=paste("debug.result5ref.",basename(IN$output),sep=""),col.names=T,row.names=F,quote=F)
}
if (IN$ss3==1) {
  result3ref<-SS3scan(lseq,lseqmin,lseqmax,newcDNApos)
  if (IN$debug==1) write.table(result3ref,file=paste("debug.result3ref.",basename(IN$output),sep=""),col.names=T,row.names=F,quote=F)
}

#integrate result
cat("Summarizing the results...\n")

for (j in min:max) {
  if (IN$ss5==1) {
  ss5pos<-which(result5ref$phypos==j)
  if (length(ss5pos)>0) {
    ss5lineref<-result5ref[ss5pos,-1]  #remove nrow column
  } else {
    ss5lineref<-data.frame(MaxEntScoreSS5=-999,SS5score=-999,SS5cDNAboundary="-",SS5phyposboundary="-",SS5seq="-",ESE5=-999,ESS5=-999,Grun5=-999,ISE5=-999)
  }
  names(ss5lineref)<-paste("Ref_",names(ss5lineref),sep="")
  ss5pos<-which(result5$phypos==j)
  if (length(ss5pos)>0) {
    ss5line<-result5[ss5pos,-1]  #remove nrow column
  } else {
    ss5line<-data.frame(MaxEntScoreSS5=-999,SS5score=-999,SS5cDNAboundary="-",SS5phyposboundary="-",SS5seq="-",ESE5=-999,ESS5=-999,Grun5=-999,ISE5=-999)
  }
  names(ss5line)<-paste("Alt_",names(ss5line),sep="")
  if (ss5lineref$Ref_MaxEntScoreSS5!=-999 & ss5line$Alt_MaxEntScoreSS5!=-999 ) diff5<-data.frame(Diff_MaxEntScoreSS5=(ss5line$Alt_MaxEntScoreSS5-ss5lineref$Ref_MaxEntScoreSS5),Diff_SS5score=(ss5line$Alt_SS5score-ss5lineref$Ref_SS5score),
                                                                                                 Diff_ESE5=(ss5line$Alt_ESE5-ss5lineref$Ref_ESE5),Diff_ESS5=(ss5line$Alt_ESS5-ss5lineref$Ref_ESS5),Diff_Grun5=(ss5line$Alt_Grun5-ss5lineref$Ref_Grun5),Diff_ISE5=(ss5line$Alt_ISE5-ss5lineref$Ref_ISE5))
  else diff5<-data.frame(Diff_MaxEntScoreSS5=-999,Diff_SS5score=-999,
                         Diff_ESE5=-999,Diff_ESS5=-999,Diff_Grun5=-999,Diff_ISE5=-999)
  }

  if (IN$ss3==1) {
  ss3pos<-which(result3ref$phypos==j)
  if (length(ss3pos)>0) {
    ss3lineref<-result3ref[ss3pos,-1]#remove nrow column
  } else {
    ss3lineref<-data.frame(MaxEntScoreSS3=-999,SS3score=-999,SS3cDNAboundary="-",SS3phyposboundary="-",SS3seq="-",ESE3=-999,ESS3=-999,Grun3=-999,ISE3=-999)
  }
  names(ss3lineref)<-paste("Ref_",names(ss3lineref),sep="")
  ss3pos<-which(result3$phypos==j)
  if (length(ss3pos)>0) {
    ss3line<-result3[ss3pos,-1]#remove nrow column
  } else {
    ss3line<-data.frame(MaxEntScoreSS3=-999,SS3score=-999,SS3cDNAboundary="-",SS3phyposboundary="-",SS3seq="-",ESE3=-999,ESS3=-999,Grun3=-999,ISE3=-999)
  }
  names(ss3line)<-paste("Alt_",names(ss3line),sep="")
  #cat("#row newseq:",nrow(newseq[j,]), " #row ss5line:",nrow(ss5line),"#row ss3line:",nrow(ss3line),"\n")
  if (ss3lineref$Ref_MaxEntScoreSS3!=-999 & ss3line$Alt_MaxEntScoreSS3!=-999 ) diff3<-data.frame(Diff_MaxEntScoreSS3=(ss3line$Alt_MaxEntScoreSS3-ss3lineref$Ref_MaxEntScoreSS3),Diff_SS3score=(ss3line$Alt_SS3score-ss3lineref$Ref_SS3score),
                                                                                                 Diff_ESE3=(ss3line$Alt_ESE3-ss3lineref$Ref_ESE3),Diff_ESS3=(ss3line$Alt_ESS3-ss3lineref$Ref_ESS3),Diff_Grun3=(ss3line$Alt_Grun3-ss3lineref$Ref_Grun3),Diff_ISE3=(ss3line$Alt_ISE3-ss3lineref$Ref_ISE3))
  else diff3<-data.frame(Diff_MaxEntScoreSS3=-999,Diff_SS3score=-999,
                         Diff_ESE3=-999,Diff_ESS3=-999,Diff_Grun3=-999,Diff_ISE3=-999)
  }
  predata<-data.frame(Var_name=mut.frame$name[i],VarName_cDNAPos=mut.frame$VarName_cDNAPos[i],VarName_ChrPos=mut.frame$VarName_ChrPos[i],Exon=mut.frame$Exon[i],Max.Exon=as.character(mut.frame$Max.Exon[i]),Protein.Coding=as.character(mut.frame$Protein.Coding[i]))
  if (IN$ss5==1 & IN$ss3==1) {
    new<-cbind(predata,newseq[which(newseq$phypos==j),],ss5lineref,ss5line,diff5,ss3lineref,ss3line,diff3)
  } else if (IN$ss5==1 & IN$ss3==-9) {
    new<-cbind(predata,newseq[which(newseq$phypos==j),],ss5lineref,ss5line,diff5)
  } else {
    new<-cbind(predata,newseq[which(newseq$phypos==j),],ss3lineref,ss3line,diff3)
  }
  result.frame<-rbind(result.frame,new)
}

#if -canonicalSS on
if (IN$canonicalSS==1 & newseq$exon[newcDNApos[1]]!=-9 & newseq$exon[newcDNApos[2]]!=-9) {  #the variant should be on an exon
  cat("Calculating Score in the canonical splice sites nearby...\n")
  #searchin for 5 prime splice site
  nexon<-newseq$exon[newcDNApos[1]]
  if (IN$ss5==1) {
  pos5ss<-0
  for (j in newcDNApos[2]:nrow(newseq)) {
    if (newseq$exon[j]==-9) {pos5ss<-j-1; break}
  }
  if (pos5ss>2) {
    result5ss_in_alt<-SS5scan(newseq,pos5ss-2,pos5ss-2,c(-9,-9))
    refpos<-which(lseq$phypos==newseq$phypos[pos5ss])
    result5ss_in_ref<-SS5scan(lseq,refpos-2,refpos-2,c(-9,-9))
    ss5lineref<-result5ss_in_ref[,-1]
    ss5line<-result5ss_in_alt[,-1]
    names(ss5lineref)<-paste("Ref_",names(ss5lineref),sep="")
    names(ss5line)<-paste("Alt_",names(ss5line),sep="")
    if (ss5lineref$Ref_MaxEntScoreSS5!=-999 & ss5line$Alt_MaxEntScoreSS5!=-999 ) diff5<-data.frame(Diff_MaxEntScoreSS5=(ss5line$Alt_MaxEntScoreSS5-ss5lineref$Ref_MaxEntScoreSS5),Diff_SS5score=(ss5line$Alt_SS5score-ss5lineref$Ref_SS5score),
                                                                                                   Diff_ESE5=(ss5line$Alt_ESE5-ss5lineref$Ref_ESE5),Diff_ESS5=(ss5line$Alt_ESS5-ss5lineref$Ref_ESS5),Diff_Grun5=(ss5line$Alt_Grun5-ss5lineref$Ref_Grun5),Diff_ISE5=(ss5line$Alt_ISE5-ss5lineref$Ref_ISE5))
    else diff5<-data.frame(Diff_MaxEntScoreSS5=-999,Diff_SS5score=-999,
                           Diff_ESE5=-999,Diff_ESS5=-999,Diff_Grun5=-999,Diff_ISE5=-999)
    if (IN$ss3==1) {
      new<-cbind(data.frame(Var_name=mut.frame$name[i]),newseq[pos5ss-2,],ss5lineref,ss5line,diff5,
               data.frame(Ref_MaxEntScoreSS3=-999,Ref_SS3score=-999,Ref_SS3cDNAboundary="-",Ref_SS3phyposboundary="-",Ref_SS3seq="-",Alt_MaxEntScoreSS3=-999,Alt_SS3score=-999,Alt_SS3cDNAboundary="-",Alt_SS3phyposboundary="-",Alt_SS3seq="-",Diff_MaxEntScoreSS3=-999,Diff_SS3score=-999,Diff_ESE3=-999,Diff_ESS3=-999,Diff_Grun3=-999,Diff_ISE3=-999))
    } else {
      new<-cbind(data.frame(Var_name=mut.frame$name[i]),newseq[pos5ss-2,],ss5lineref,ss5line,diff5)
    }
      result.frame<-rbind(result.frame,new)
  }
  }
  
  if (IN$ss3==1) {
  pos3ss<-0
  for (j in newcDNApos[1]:1) {
    if (newseq$exon[j]==-9) {pos3ss<-j+1; break}
  }
  if (pos3ss>20)  {
    result3ss_in_alt<-SS3scan(newseq,pos3ss-20,pos3ss-20,c(-9,-9))
    refpos<-which(lseq$phypos==newseq$phypos[pos3ss])
    result3ss_in_ref<-SS3scan(lseq,refpos-20,refpos-20,c(-9,-9))
    ss3lineref<-result3ss_in_ref[,-1]
    ss3line<-result3ss_in_alt[,-1]
    names(ss3lineref)<-paste("Ref_",names(ss3lineref),sep="")
    names(ss3line)<-paste("Alt_",names(ss3line),sep="")
    if (ss3lineref$Ref_MaxEntScoreSS3!=-999 & ss3line$Alt_MaxEntScoreSS3!=-999 ) diff3<-data.frame(Diff_MaxEntScoreSS3=(ss3line$Alt_MaxEntScoreSS3-ss3lineref$Ref_MaxEntScoreSS3),Diff_SS3score=(ss3line$Alt_SS3score-ss3lineref$Ref_SS3score),
                                                                                                   Diff_ESE3=(ss3line$Alt_ESE3-ss3lineref$Ref_ESE3),Diff_ESS3=(ss3line$Alt_ESS3-ss3lineref$Ref_ESS3),Diff_Grun3=(ss3line$Alt_Grun3-ss3lineref$Ref_Grun3),Diff_ISE3=(ss3line$Alt_ISE3-ss3lineref$Ref_ISE3))
    else diff3<-data.frame(Diff_MaxEntScoreSS3=-999,Diff_SS3score=-999,Diff_ESE3=-999,Diff_ESS3=-999,Diff_Grun3=-999,Diff_ISE3=-999)
    if (IN$ss5==1) {
    new<-cbind(data.frame(Var_name=mut.frame$name[i]),newseq[pos3ss-20,],data.frame(Ref_MaxEntScoreSS5=-999,Ref_SS5score=-999,Ref_SS5cDNAboundary="-",Ref_SS5phyposboundary="-",Ref_SS5seq="-",Alt_MaxEntScoreSS5=-999,Alt_SS5score=-999,Alt_SS5cDNAboundary="-",Alt_SS5phyposboundary="-",Alt_SS5seq="-",Diff_MaxEntScoreSS5=-999,Diff_SS5score=-999,Diff_ESE5=-999,Diff_ESS5=-999,Diff_Grun5=-999,Diff_ISE5=-999),
               ss3lineref,ss3line,diff3)
    } else {
      new<-cbind(data.frame(Var_name=mut.frame$name[i]),newseq[pos3ss-20,],ss3lineref,ss3line,diff3)
    }
    result.frame<-rbind(result.frame,new)
  }
  }
}
} else {  #if cDNApos not assigned
  cat("Variant:",mut.frame$name[i], " is NOT on the indicated transcript:",mut.frame$transID[i], "propably because it is on the downstream or UTR regions.\n")
}

}  #mut.frame BIG i loop end
} else {  #just treat ref seq
  cat("Calculating Scores just in the referrence sequence...\n")
  if (IN$ss5==1) {
  cat("Checking the potential of 5\'splice site...\n")
  Start<-newcDNApos[1]-9+1
  if (Start<1) Start<-1  #indicate the region where the mutation affects
  End <-newcDNApos[2]
  if (End > (nrow(lseq)-9+1)) End<-nrow(lseq)-9+1
  result5ref<-SS5scan(lseq,Start,End,newcDNApos)
  }
  if (IN$ss3==1) {
  cat("Checking the potential of 3\'splice site...\n")
  Start<-newcDNApos[1]-23+1
  if (Start<1) start<-1  #indicate the region where the mutation affects
  End <-newcDNApos[2]
  if (End > (nrow(lseq)-23+1)) end<-nrow(lseq)-23+1
  result3ref<-SS3scan(lseq,Start,End,newcDNApos)
  }
  
  #integrate result
  cat("Summarizing the results...\n")
  if (IN$ss3==-9) {
    minmax<-minmaxfunc(result5ref$phypos)
  } else if (IN$ss5==-9) {
    minmax<-minmaxfunc(result3ref$phypos)
  } else {
    minmax<-minmaxfunc(c(result3ref$phypos,result5ref$phypos))
  }
  min<-minmax$min
  max<-minmax$max
  
  for (j in min:max) {
    if (IN$ss5==1) {
    ss5pos<-which(result5ref$phypos==j)
    if (length(ss5pos)>0) {
      ss5lineref<-result5ref[ss5pos,-1]  #remove nrow column
    } else {
      ss5lineref<-data.frame(MaxEntScoreSS5=-999,SS5score=-999,SS5cDNAboundary="-",SS5phyposboundary="-",SS5seq="-",ESE5="-",ESS5="-",Grun5="-",ISE5="-")
    }
    names(ss5lineref)<-paste("Ref_",names(ss5lineref),sep="")
    }
    if (IN$ss3==1) {
    ss3pos<-which(result3ref$phypos==j)
    if (length(ss3pos)>0) {
      ss3lineref<-result3ref[ss3pos,-1]#remove nrow column
    } else {
      ss3lineref<-data.frame(MaxEntScoreSS3=-999,SS3score=-999,SS3cDNAboundary="-",SS3phyposboundary="-",SS3seq="-",ESE3="-",ESS3="-",Grun3="-",ISE3="-")
    }
    names(ss3lineref)<-paste("Ref_",names(ss3lineref),sep="")
    }
    if (IN$ss5==1 & IN$ss3==1) {
      new<-cbind(lseq[which(lseq$phypos==j),],ss5lineref,ss3lineref)
    } else if (IN$ss5==1 & IN$ss3==-9) {
      new<-cbind(lseq[which(lseq$phypos==j),],ss5lineref)
    } else {
      new<-cbind(lseq[which(lseq$phypos==j),],ss3lineref)
    }
    
    result.frame<-rbind(result.frame,new)
  }
}


#save result
cat("Writing splice site estimation result to:",IN$output,"\n")
write.table(result.frame,file=IN$output,row.names=F,col.names=T,sep="\t",quote=F)

#attache splice junction info if required  ##sjdb requires Kind column which is attached in TakeLines function
#if (IN$sjdbout!=-9) {
#  sjdb(IN$output)
#}#############sjdb will be done in summarizeresult *.SigLines.txt

#summarize result if required
if (IN$summarizeresult!=-9) {
  TakeSigLines(IN$output,IN$summarizeresult)
}

cat("Done.\n")

warnings()











