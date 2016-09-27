Helpdocument<-"
#this script is to make a construct for in-vitro splicing assay
#
#Usage: R --slave --vanilla --args  -mutfile (chrpos format mut file with transcript ID) -output (output file name) < Regress_Score.v0.5.R
#
#ex) R --slave --vanilla --args -refdirectory ../../MAE-working/5ss-hotspot/ssfiles -mutfile test.variant.list.txt -ss5 -output temp.result.txt -esbTranscriptID ENST00000368300 <ConstructDesigner.v.0.5.R
#ex) R --slave --vanilla --args -refdirectory ../../MAE-working/5ss-hotspot/ssfiles -mutfile test.variant.list.txt -ss5 -output temp.result.v0.6.txt -esbTranscriptID ENST00000368300 -readbarcodes 8bpBARCODE.txt -seq5 CTGGTTTAGTGAACCGTCAG -seq3 ATCTAGATAACTGATCATAATCAGCCATACCACATTTGTAGAGGTTTTACTTGCTTTAAAAAACCTCCCACAXXXXXXXXTGAACCTGAAACATAAAATGAATGCAATTGTTGTTGTTAACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTGCATTCTAGTTGTGGTTTGTCCAAACTCATCAATGTATCTTA <ConstructDesigner.v.0.6.R
#ex) R --slave --vanilla --args -refdirectory ../../MAE-working/5ss-hotspot/ssfiles -mutfile test.variant.list.txt -ss5 -output temp.result.v0.6.txt -esbTranscriptID ENST00000368300 -readbarcodes 8bpBARCODE.txt -seq5 CTGGTTTAGTGAACCGTCAG -seq3 ATCTAGATAACTGATCATAATCAGCCATACCACATTTGTAGAGGTTTTACTTGCTTTAAAAAACCTCCCACAXXXXXXXXTGAACCTGAAACATAAAATGAATGCAATTGTTGTTGTTAACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTGCATTCTAGTTGTGGTTTGTCCAAACTCATCAATGTATCTTA -randombarcode -constructlength 20,50,50,45,55,50  <ConstructDesigner.v.0.6.R
#ex) R --slave --vanilla --args -refdirectory ../../MAE-working/5ss-hotspot/ssfiles -mutfile test.variant.list.txt -ss5 -output temp.result.v0.6.txt -esbTranscriptID ENST00000368300 -readbarcodes 8bpBARCODE.txt -seq5 CTGGTTTAGTGAACCGTCAG -seq3 ATCTAGATAACTGATCATAATCAGCCATACCACATTTGTAGAGGTTTTACTTGCTTTAAAAAACCTCCCACAXXXXXXXXTGAACCTGAAACATAAAATGAATGCAATTGTTGTTGTTAACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTGCATTCTAGTTGTGGTTTGTCCAAACTCATCAATGTATCTTA -randombarcode -constructlength 40,30,30,45,55,40  <ConstructDesigner.v.0.6.R
#ex) R --slave --vanilla --args -refdirectory ../../MAE-working/5ss-hotspot/ssfiles -mutfile LMNA-N500candidates_ENST00000368300.txt -ss5 -output LMNA-N500candidates_ENST00000368300.New5SS.result.v0.7.txt -esbTranscriptID ENST00000368300 -barcode AA,AT,AG,AC,TA,TT,TG,TC,GA,GT,GG,GC,CA,CT,CG,CC -seq5 TTACGCCAAGTTATTTAGGTGACA -seq3 XXATCTAGATAACTGATCATAATCAGCCATACCACATTTGT -randombarcode -constructlength 85,85,85,85,55,40 -id 6,2 <ConstructDesigner.v.0.7.R
#ex) R --slave --vanilla --args -refdirectory ../../MAE-working/5ss-hotspot/ssfiles -mutfile LMNA-N500candidates_Exac.txt -ss5 -output LMNA-N500candidates_Exac.NEW5SS.result.v0.7.txt -barcode AA,AT,AG,AC,TA,TT,TG,TC,GA,GT,GG,GC,CA,CT,CG,CC -seq5 TTACGCCAAGTTATTTAGGTGACA -seq3 XXATCTAGATAACTGATCATAATCAGCCATACCACATTTGT -randombarcode -constructlength 85,85,85,85,55,40 -id 6,2<ConstructDesigner.v.0.7.R
#ex) R --slave --vanilla --args -refdirectory /Users/meizhu/Copy/MAE-working/5ss-hotspot/ssfiles -mutfile Var_LMM_2ndbatch.txt -ss5 -output Var_LMM_2ndbatch.result.v0.8.txt -esbTranscriptID ENST00000368300 -barcode AA,AT,AG,AC,TA,TT,TG,TC,GA,GT,GG,GC,CA,CT,CG,CC -seq5 TTACGCCAAGTTATTTAGGTGACA -seq3 XXATCTAGATAACTGATCATAATCAGCCATACCACATTTGT -randombarcode -constructlength 85,85,85,85,55,40 -id 6,2 <../ConstructDesigner.v.0.8.R
#ex) R --slave --vanilla --args -refdirectory /Users/meizhu/Copy/MAE-working/5ss-hotspot/ssfiles -mutfile Var_LMM1664_2ndbatch.txt -ss5 -output Var_LMM1664_2ndbatch.result.v0.8.txt -esbTranscriptID ENST00000368300 -barcode AA,AT,AG,AC,TA,TT,TG,TC,GA,GT,GG,GC,CA,CT,CG,CC -seq5 TTACGCCAAGTTATTTAGGTGACA -seq3 XXATCTAGATAACTGATCATAATCAGCCATACCACATTTGT -randombarcode -constructlength 55,35,0,128,98,119 -id 6,2 <../ConstructDesigner.v.0.8.R
#ex) R --slave --vanilla --args -refdirectory /Users/meizhu/Copy/MAE-working/5ss-hotspot/ssfiles -mutfile Var_Exac_2ndbatch.txt -ss5 -output Var_Exac_2ndbatch.result.v0.8.txt -esbTranscriptID ENST00000368300 -barcode AA,AT,AG,AC,TA,TT,TG,TC,GA,GT,GG,GC,CA,CT,CG,CC -seq5 TTACGCCAAGTTATTTAGGTGACA -seq3 XXATCTAGATAACTGATCATAATCAGCCATACCACATTTGT -randombarcode -constructlength 85,85,85,85,55,40 -id 6,2 <../ConstructDesigner.v.0.8.R
#ex) R --slave --vanilla --args -refdirectory /Users/meizhu/Copy/MAE-working/5ss-hotspot/ssfiles -summarizeConstructs Var_LMM_2ndbatch.result.v0.8.txt,Var_Exac_2ndbatch.result.v0.8.txt,Var_LMM1664_2ndbatch.result.v0.8.txt -output LMNA_Var_2ndbatch.result.v0.8.txt -barcode AA,AT,AG,AC,TA,TT,TG,TC,GA,GT,GG,GC,CA,CT,CG,CC -seq5 TTACGCCAAGTTATTTAGGTGACA -seq3 XXATCTAGATAACTGATCATAATCAGCCATACCACATTTGT -id 6,2 <../ConstructDesigner.v.0.8.R
#ex) R --slave --vanilla --args -refdirectory /Users/meizhu/Copy/MAE-working/5ss-hotspot/ssfiles -mutfile LMM_Candidates.txt -ss5 -output Var_LMM_MYBPC3.result.v0.81.txt -esbTranscriptID ENST00000545968 -barcode AA,AT,AG,AC,TA,TT,TG,TC,GA,GT,GG,GC,CA,CT,CG,CC -seq5 TTACGCCAAGTTATTTAGGTGACA -seq3 XXATCTAGATAACTGATCATAATCAGCCATACCACATTTGT -randombarcode -constructlength 85,85,85,85,55,40 -id 6,2 -limitintron 200 -limit2ndexon 120 <../ConstructDesigner.v.0.81.R
#ex) R --slave --vanilla --args -refdirectory /Users/meizhu/Copy/MAE-working/5ss-hotspot/ssfiles -mutfile Exac_Candidates.txt -ss5 -output Var_Exac_MYBPC3.result.v0.81.txt -barcode AA,AT,AG,AC,TA,TT,TG,TC,GA,GT,GG,GC,CA,CT,CG,CC -seq5 TTACGCCAAGTTATTTAGGTGACA -seq3 XXATCTAGATAACTGATCATAATCAGCCATACCACATTTGT -randombarcode -constructlength 85,85,85,85,55,40 -id 6,2 -limitintron 200 -limit2ndexon 120 <../ConstructDesigner.v.0.81.R
#ex) R --slave --vanilla --args -refdirectory /Users/meizhu/Copy/MAE-working/5ss-hotspot/ssfiles -summarizeConstructs Var_LMM_MYBPC3.result.v0.81.txt,Var_Exac_MYBPC3.result.v0.81.txt -output LMNA_Var_MYBPC3.result.v0.81.txt -barcode AA,AT,AG,AC,TA,TT,TG,TC,GA,GT,GG,GC,CA,CT,CG,CC -seq5 TTACGCCAAGTTATTTAGGTGACA -seq3 XXATCTAGATAACTGATCATAATCAGCCATACCACATTTGT -id 6,2 <../ConstructDesigner.v.0.81.R
#ex) R --slave --vanilla --args -input LMNA_Var_MYBPC3.result.v0.81.txt -output LMNA_Var_MYBPC3.result.v0.81.deldup.txt <../ConstructDesigner.DelDup.R
#ex) R --slave --vanilla --args -refdirectory /Users/meizhu/Copy/MAE-working/5ss-hotspot/ssfiles -summarizeConstructs LMNA_Var_MYBPC3.result.v0.81.deldup.txt -output LMNA_Var_MYBPC3.result.v0.81.deldup+.txt -barcode AA,AT,AG,AC,TA,TT,TG,TC,GA,GT,GG,GC,CA,CT,CG,CC -seq5 TTACGCCAAGTTATTTAGGTGACA -seq3 XXATCTAGATAACTGATCATAATCAGCCATACCACATTTGT -id 6,2 <../ConstructDesigner.v.0.81.R
#ex) R --slave --vanilla --args -refdirectory /Users/meizhu/Copy/MAE-working/5ss-hotspot/ssfiles -mutfile LMNA-MYBPC3-mut.txt -ss3 -output LMNA_MYBPC3.ss3.result.v0.9.txt -barcode AA,AT,AG,AC,TA,TT,TG,TC,GA,GT,GG,GC,CA,CT,CG,CC -seq5 TTACGCCAAGTTATTTAGGTGACA -seq3 XXATCTAGATAACTGATCATAATCAGCCATACCACATTTGT -randombarcode -constructlength 40,55,85,85,85,85 -id 6,2 -limitintron 200 -limit1stexon 120 <../ConstructDesigner.v.0.9.R
#ex) R --slave --vanilla --args -refdirectory /Users/meizhu/Copy/MAE-working/5ss-hotspot/ssfiles -mutfile LMNA-MYBPC3-mut.txt -ss3 -output LMNA_MYBPC3.ss3.result.v0.9.txt -barcode AA,AT,AG,AC,TA,TT,TG,TC,GA,GT,GG,GC,CA,CT,CG,CC -seq5 TTACGCCAAGTTATTTAGGTGACA -seq3 XXATCTAGATAACTGATCATAATCAGCCATACCACATTTGT -randombarcode -constructlength 40,55,85,85,85,85 -id 6,2 -limitintron 200 -limit1stexon 120 -summarizeConstructs current -gblocksubmit gBlockOrder_LMNA_MYBPC3.ss3.result+.v0.91.txt <../ConstructDesigner.v.0.92.R
#Options:
#-refdirectory (directory where maxent and other files are located) :when these files are no loacted on the current directory , please indicate their locations
#-mutfile (chrpos format mut file with transcript ID) e.g.chr1:123G>A ENST000001234 in v0.5 just 1 line accept
#-output (filename to be saved as a result)
#-marthost (indicates biomaRt server)
#
#-seq5 (indicates 5 prime sequence): the default seq is the 20-bp tail seq of CMV promotor
#-seq3 (indicates seq attached to 3 end): the default seq is SV40polyA (240bp)  
### in seq5 or seq3, if you want to attach barcode, please indicate the location where you want to insert barcode using X (capital X)
#-barcode (seq):indicates to attach the barcode you indicated.  you can give greater than 2 barcodes using comman-delimiter (e.g. ATG,GTG,CAC,GGG)
#-randombarcode : barcodes are randomly selected  ##if not barcodes are employed according to the order
#          * default setting is no barcode
#-readbarcodes (file) : indicates a file that contains barcodes. (barcodes should be separated by \n  in the file)
#-constructlength 
#                 in SS5
#                   1-#bp of exon before the mutation               2-#bp of exon after the mutation,
#                   3-#bp of exon before canonical splice site      4-#bp of intron adjacent to 3 end of 1st exon 
#                   5-#bp of intron adjacent to 5 end of 2nd exon)  6-#bp of 2nd exon (from 5 end)
#                   defualt setting is 85,85,85,85,55,40  in -ss5 mode
#                 in SS3
#                   1-#bp in 1st exon                               2-#bp in intron1
#                     if variant on exon  
#                       3-#bp in intron2                                4-#bp in Exon2 after canonical acceptor site
#                       5-#bp in exon2 before mutation                  6-#bp in exon after mut
#                     if variant on intron
#                       3-#bp-intron2 before mut                        4-#bp in intron2 after mut
#                       5-#bp-intron2 before canonical acceptor site    6-#bp in exon2
#           default setting is 40,55,85,85,85,85 in -ss3 mode
#-ss5 indicates this script to make a construct for 5 prime splice site
#-ss3 indicates this script to make a construct for 3 prime splice site
#-id (integer1,integer2) indicates to generate Specific Identifier for the construct using the barcode and the peripheral seqs. In concrete, (integer1 bp before barcode)+(barcode)+(integer2 bp after barcode) will be an identifier.
#                 default setting is 0,0
#-seqgrep (integer1,integer2) indicates # of additionl seqs (integer1 for left side and integer2 for right side) around the mutation for seq grep, which is seqs when you want to check reference seq, mutant seq, normal splice and aberrant splice
#                 default values are 8,8
#-summarizeConstructs (current or filenames delimited with comma)   if you give the word current, this script summarizes the result it generates current, if you give filenames (when you provide more than one file, they should be delimited by comma), given files that are generated by this script will be summarized.
#               In concrete, In constructs with a same transcriptID, First_Exon, Normal_SpliceGap, we need to prepare just one reference construct for them. So Grouping constructs then giving barcode in a sophisticated way
#               *** when you just summarize generated results, please give -id, -ss5 and -ss3 and specify barcode by -barcode or -readbarcodes
#-limitintron (integer) when you want to limit the length of intron, indicate the maximam length of intron  in both SS5 and SS3 mode
#-limit2ndexon (integer) when you want to limit the length of 2ndexon, indicate the maximam length of 2ndexon  in SS5 mode
#-limit1stexon (interger) when you want to limit the length of 1stexon, indicate the maximam length of 1stexon  in SS3 mode
#-gblocksubmit (filename) when you want to make a file for gBlock order, please indicate the filename for it. the file constains construct file name and actual seqs. Redundant constructs will be removed because ONE VARIANT COULD HAVE 2LINES WHEN IT WORKS AS CREATE AND BROKEN
#-scoreresult (filename) when you have a result of Regress_Scorev0.85 or later, you can omit score calculation and make working time short.
#                       if you use this switch, you do not need to use -mutfile 
#basic architecture is from Regress_Score.v0.8.R
#most of element search sub routines came from DataCollecterPartA.v0.9.R and they had been modified to fit this current situation.
#most of constructing ref seq and searching alt seq code came from 5ss_cal_ver1.19.R
#
#output file format is
#Var_Name Construct_Seq Seq_NormalSplice Seq_AberrantSplice Length_NormalSplice Length_AberrantSplice Seq_for_Primer3,Ref_MaxEnt,Ref_Integrated,ALt_MaxEnt,Alt_Integrated,Diff_MaxENt,Diff_Integrated  
#In construct seq EXONs are UPPER LETTER, introns are Lower letter, 
# attached 5 seq ,3seq and barcode are the opposite letter of the adjacent element.(if neigherhood is exon, these element will be lower case)(if barcode is embedded in lowercase 3 attached seq, barcode is uppercase)
#In NormalSplice and AberrantSplice seqs, the region that will be spliced out is surrrounded by () 
#In Seq for primer3, the seqs inside of the 5 attached seq and 3 attached seq will be surrounded [ ]. (if barcode indicated, the barcode will be included in the content of [ ])
#Length_NormalSplice indicates 1st exon-2ndexon length if normal splice occurs
#Length_AberrantSplice indicates 1st exon-2ndexon length if aberrant splice occurs
#
#History
#v0.5 launch ver.
#v0.6 barcode implemented (no barcode implement in v0.5)  just -ss5 mode implemented   (-ss3 not implemented)
#v0.61 ID mode added when the length of barcode is short (-id)
#v0.7 
#bug-fix in broken-5SS variants
#bug-fix in variants on a intron (near splice donor site)
#bug-fix when the variant is on the last exon (can not make a construct because of no 2nd-exon)
#give intuitive construct summary ( e.g. mutant point (A/T)  aberrant splice donor site{  aberrant splice acceptor site} normal splice donor site[ normal splice acceptor site]  )
#give gap position and mutant allele position
#when given primers, calculate gap position, mut allele position and predicted product lengths
#v0.8
#In constructs with a same transcriptID, First_Exon, Normal_SpliceGap, we need to prepare just one reference construct for them. So Grouping constructs then giving barcode in a sophisticated way using -summarizeConstructs switch
#bug fix in v0.7 -- mutant on exon & broken SS -->SHORT exon2nd
#V0.81 
#when intron length is >200bp, splicing assay sometimes failed in N500batch2. So add a switch to limit the length of intron <=x bp   -limitintron (integer)
#I guess longer 2nd exon could also prevent us from reading construct correctly because RT primer is just R-primer. So I added a switch to limit then length of 2nd exon -limit2ndexon (integer)
#v0.9
#can generate constructs for 3SS variation  -ss3
#add and refine takeseq subroutine
#refine intronlimit subroutine (when combined componets, the latter length become zero. but when down sizing, zero could be a problem in case of not combined. then fixed the problem)
#in the future -ss5 and -ss3 in constructmake will be able to be combined beucase they share many components
#
#bigfix in changeseq function. In concrete, if surrounding seq around ID is not specific, ID change works wrong. 
#                              When barcode is not unique, take former or later one according to the location of barcode (5ss seq or 3ss seq)
#In -summarizeConstructs function, add the function to delete duplicates (ConstructDesigner.DelDup.R modified for both 5SS and 3SS use)
v0.91
#detect duplicated IDs and inform it in -summarizeConstructs .
v0.92
#g-block submit mode implemented by -gblocksubmit (outfilename)
#v0.93
#if you have the result of regress score v0.85 or later, this script understand the result and omit the score calculation using -scoreresult () switch,  which makes work time short
#UNDER CONSTRUCTION
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
  output<-list(Maxent=0,SSscore=0)
  if (SS==5) { #5ss processing
    ss5seq<-substr(seq,startend[2]-2,startend[2]+6)
    if (nchar(ss5seq)==9) output$Maxent<-maxent5score(ss5seq) else output$Maxent<--99
    #cat("Submit exon:",substr(seq,startend[1],startend[2]), " Submit intron:",substr(seq,startend[2]+1,nchar(seq)),sep="")
    Escore<-exonscan(substr(seq,startend[1],startend[2]),3)
    Iscore<-intronscan(substr(seq,startend[2]+1,nchar(seq)),1,6)
    logOdds<-RegCoef5$Intercept+RegCoef5$MAXENT5*output$Maxent+RegCoef5$ESE3*Escore$ESE+RegCoef5$ESS3*Escore$ESS+RegCoef5$aGrun*Iscore$Grun+RegCoef5$aISE*Iscore$ISE
    if (output$Maxent!=-99) output$SSscore<-1/(1+exp(-logOdds)) else output$SSscore<--99
    output$SSseq<-paste(substr(ss5seq,1,3),tolower(substr(ss5seq,4,9)),sep="")
    #cat("SS5seq:",substr(ss5seq,1,3),tolower(substr(ss5seq,4,9)),"\n",sep="")
  } else { #3ss processing
    ss3seq<-substr(seq,startend[1]-20,startend[1]+2)
    #cat("SS3seq:",ss3seq,"\n")
    if ((nchar(ss3seq)==23) & length(grep("N",ss3seq))==0) output$Maxent<-maxent3score(ss3seq) else output$Maxent<--99
    #cat("Submit exon:",substr(seq,startend[1],startend[2]), " Submit intron:",substr(seq,1,startend[1]-1),sep="")
    Escore<-exonscan(substr(seq,startend[1],startend[2]),5)
    Iscore<-intronscan(substr(seq,1,startend[1]-1),-1,20)
    logOdds<-RegCoef3$Intercept+RegCoef3$MAXENT3*output$Maxent+RegCoef3$ESE5*Escore$ESE+RegCoef3$ESS5*Escore$ESS+RegCoef3$bGrun*Iscore$Grun+RegCoef3$bISE*Iscore$ISE
    if (output$Maxent!=-99) output$SSscore<-1/(1+exp(-logOdds)) else output$SSscore<--99
    output$SSseq<-paste(tolower(substr(ss3seq,1,20)),substr(ss3seq,21,23),sep="")
    #cat("SS3seq:",tolower(substr(ss3seq,1,20)),substr(ss3seq,21,23),"\n",sep="")
  } 
  return(output)
}

SS5scan<-function(newseq,Start,End) {
  #cat("Seq: ",tolower(paste(newseq$Seq[1:(Start-1)],collapse="")),paste(newseq$Seq[Start:End],collapse=""),tolower(paste(newseq$Seq[(End+1):nrow(newseq)],collapse="")), " Start:",Start," End:",End,"\n")
  result5<-data.frame()
  if (Start>End) {temp<-c(Start,End);Start<-temp[2];End<-temp[1]}
  for (j in Start:End) {
    exonend<-j+2    #putative exon end generated by the mutation
    if (exonend>3 & exonend<=nrow(newseq)) { 
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
      #cat("len newseqphypos",length(newseq$phypos[j])," nrow output:",unlist(output), " length cdnaboundary:",length(cDNAboundary)," length phyposboundary:",length(phyposboundary),"\n")
      new<-data.frame(phypos=newseq$phypos[j],MaxEntScoreSS5=output$Maxent,SS5score=output$SSscore,SS5cDNAboundary=cDNAboundary,SS5phyposboundary=phyposboundary,SS5seq=output$SSseq)
      result5<-rbind(result5,new)
    }
  }
  return(result5)
}

SS3scan<-function(newseq,Start,End) {
  #cat("Seq: ",tolower(paste(newseq$Seq[1:(Start-1)],collapse="")),paste(newseq$Seq[Start:End],collapse=""),tolower(paste(newseq$Seq[(End+1):nrow(newseq)],collapse="")), " Start:",Start," End:",End,"\n")
  result3<-data.frame()
  if (Start>End) {temp<-c(Start,End);Start<-temp[2];End<-temp[1]}
  for (j in Start:End) {
    exonstart<-j+20    #putative exon end generated by the mutation
    if (exonstart>21 & exonstart<=nrow(newseq)) { 
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
      new<-data.frame(phypos=newseq$phypos[j],MaxEntScoreSS3=output$Maxent,SS3score=output$SSscore,SS3cDNAboundary=cDNAboundary,SS3phyposboundary=phyposboundary,SS3seq=output$SSseq)
      result3<-rbind(result3,new)
    }
  }
  return(result3)
}

GetSeq<-function(transciptID,cDNAsw) { #ID -transcriptID, Seq -each allele, phypos -physical position,cDNApos-cDNA position, exon-#exon
  transcriptinfo<-getBM(attributes=c("hgnc_symbol","ensembl_transcript_id","chromosome_name","transcript_start","transcript_end","strand"),
                        "ensembl_transcript_id", values=transciptID, mart=mart)
  if (nrow(transcriptinfo)==0 ) return ( data.frame() )
  exonsinfo<-getBM(attributes=c("ensembl_transcript_id","ensembl_exon_id","rank","chromosome_name","exon_chrom_start","exon_chrom_end"),
                   "ensembl_transcript_id", values=transciptID, mart=mart)
  utr5info<-getBM(attributes=c("5_utr_start","5_utr_end"),"ensembl_transcript_id", values=transciptID, mart=mart)
  utr3info<-getBM(attributes=c("3_utr_start","3_utr_end"),"ensembl_transcript_id", values=transciptID, mart=mart)
  utr5info<-utr5info[!is.na(utr5info[,1]),]
  utr3info<-utr3info[!is.na(utr3info[,1]),]
  utrinfo<-data.frame(start=c(utr5info$"5_utr_start",utr3info$"3_utr_start"),end=c(utr5info$"5_utr_end",utr3info$"3_utr_end"))
  lseqtemp     <- getSequence(id = transciptID, type="ensembl_transcript_id", mart =mart, seqType = "transcript_exon_intron") 
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
        if (nrow(utrinfo)>0) for (j in 1:nrow(utrinfo)) {
          if (start<utrinfo[j,]$end & utrinfo[j,]$end<end) start<-utrinfo[j,]$end+1
          if (start<utrinfo[j,]$start & utrinfo[j,]$start<end) end<-utrinfo[j,]$start-1
          if (utrinfo[j,]$start<=start & end<=utrinfo[j,]$end) all<-1
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
  return(lseq)
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
    if (length(grep("^(A|T|G|C|.)+$",ref))>0 & length(grep("^(A|T|G|C|.)+$",alt))>0) {
      cat("\tAccepted\n")
      new.frame<-data.frame(kind=I("chrpos"),chr=I(chr),pos=as.integer(pos),ref=I(ref),alt=I(alt),name=I(paste("chr",chr,":",pos,ref,">",alt,sep="")),transID=I(transID))
      if (nrow(mut.frame)==0) mut.frame<-new.frame else mut.frame<-rbind(mut.frame,new.frame)
      #cat("mut.frame chr:",chr,"pos:",pos,"ref:",ref,"alt:",alt,"\n")
    } else {
      cat("\tSkipped\n")
    }
  }
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
      
      
      new.frame<-data.frame(kind=I(kind),chr=I(NA),pos=I(pos),ref=I(ref),alt=I(alt),name=I(paste("c.",pos,ref,ifelse(length(grep(">",muttemp[i]))>0,">",""),alt,sep="")),transID=I(IN$esbTranscriptID))
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

makenewseq<-function(lseq,mut,newcDNApos) {
newseq<-lseq
if (mut$kind=="chrpos" & length(grep("^(A|T|G|C)$",mut$ref))==1  & length(grep("^(A|T|G|C)$",mut$alt))==1 | mut$kind=="coding" ) {
  newseq[newcDNApos[1],]$Seq<-mut$alt;indelflag<-0;
} else if (mut$kind=="coding_ins") {
  insert<-data.frame()
  for (j in 1:nchar(mut$alt)) {
    temp<-lseq[1,]
    temp$Seq<-substr(mut$alt,j,j)
    temp$phypos<-paste(lseq[newcDNApos[1],]$phypos,"+ins",j,sep="")
    temp$cDNApos<-paste(lseq[newcDNApos[1],]$cDNApos,"+ins",j,sep="")
    temp$exon<-lseq[newcDNApos[1],]$exon
    insert<-rbind(insert,temp)
  }
  newseq<-rbind(newseq[1:newcDNApos[1],],insert,newseq[newcDNApos[2]:nrow(newseq),])
  newcDNApos[1]<-newcDNApos[1]+1
  newcDNApos[2]<-newcDNApos[2]+nrow(insert)-1
  
} else if (mut$kind=="coding_del") {
  delseq<-""
  for (j in newcDNApos[1]:newcDNApos[2]) {delseq<-paste(delseq,newseq[j,]$Seq,sep="");}
  #cat("Deleted Seq is :",toupper(delseq),"mut.alt:",mut$alt,"\n")####
  if (nchar(mut$alt)==0 | is.na(mut$alt)) cat("Deleted Seq is :",toupper(delseq),"\n") else if (toupper(delseq)!=toupper(mut$alt)) cat("Warning. Given del seq:",toupper(mut$alt),"is different from actual cDNA seq:",toupper(delseq),"\n")
  newseq<-rbind(newseq[1:(newcDNApos[1]-1),],newseq[(newcDNApos[2]+1):nrow(newseq),])
  #newcDNApos[1]<-newcDNApos[1]
  newcDNApos[2]<-newcDNApos[1]-1
} else if (mut$kind=="coding_delins") {
  insert<-data.frame()
  for (j in 1:nchar(mut$alt)) {
    temp<-lseq[1,]
    temp$Seq<-substr(mut$alt,j,j)
    temp$phypos<-paste(lseq[newcDNApos[1],]$phypos,"+ins",i,sep="")
    temp$cDNApos<-paste(lseq[newcDNApos[1],]$cDNApos,"+ins",j,sep="")
    temp$exon<-lseq[newcDNApos[1],]$exon
    insert<-rbind(insert,temp)
  }
  newseq<-rbind(newseq[1:(newcDNApos[1]-1),],insert,newseq[(newcDNApos[2]+1):nrow(newseq),])
  newcDNApos[2]<-newcDNApos[1]+nrow(insert)-1
} else if (mut$kind=="coding_del_ins") {
  delins<-unlist(strsplit(mut$alt,","))
  delseq<-"";
  for (j in newcDNApos[1]:newcDNApos[2]) {delseq<-paste(delseq,newseq[j,]$Seq,sep="");}
  if (toupper(delseq)!=toupper(delins[1])) cat("Warning. Given del seq:",toupper(mut$alt),"is different from actual cDNA seq:",toupper(delseq),"\n")
  insert<-data.frame()
  for (j in 1:nchar(mut$alt)) {
    temp<-lseq[1,]
    temp$Seq<-substr(mut$alt,j,j)
    temp$phypos<-paste(lseq[newcDNApos[1],]$phypos,"+ins",j,sep="")
    temp$cDNApos<-paste(lseq[newcDNApos[1],]$cDNApos,"+ins",j,sep="")
    temp$exon<-lseq[newcDNApos[1],]$exon
    insert<-rbind(insert,temp)
  }
  newseq<-rbind(newseq[1:(newcDNApos[1]-1),],insert,newseq[(newcDNApos[2]+1):nrow(newseq),])
  newcDNApos[2]<-newcDNApos[1]+nrow(insert)-1
} else if (mut$kind=="coding_dup") {
  insertseq<-"";
  for (j in newcDNApos[1]:newcDNApos[2]) {
    insertseq<-paste(insertseq,toupper(lseq[j,]$Seq),sep="")
  }
  if (length(grep("\\d+",mut$alt))>0) {
    insertseq2<-"";
    for (j in 1:(as.numeric(mut$alt)-1)) {
      insertseq2<-paste(insertseq2,insertseq,sep="")
    }
  } else if (nchar(mut$alt)==0) {
    insertseq2<-insertseq
  } else {
    insertseq2<-toupper(mut$alt)
    #cat ("insertseq",insertseq,"insertseq2",insertseq2,"i",i,"ref",mut$ref,"alt",class(mut$alt),"\n")
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
} else if (mut$kind=="chrpos") { #in case of indel of vcf file
  delseq<-"";
  ncommon<-0;
  
  if (mut$ref=="." ) mut$ref<-""
  if (mut$alt=="." ) mut$alt<-""
  ref <-mut$ref
  alt <-mut$alt
  add<-ifelse(nchar(mut$alt)>0,1,0)  #for length(alt)==0 
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
    refstart<-newcDNApos[1]-nchar(mut$ref) #just left position of substitution
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
    refend<-newcDNApos[1]+nchar(mut$ref) #just right position of substitution
  }
  newcDNApos[1]<-refstart+1
  newcDNApos[2]<-refstart+nchar(newalt)
  
  for (j in (refstart+1-ifelse(reverseflag==0,ncommon,0)):(refend-1+ifelse(reverseflag==1,ncommon,0))) {delseq<-paste(delseq,newseq[j,]$Seq,sep="");}
  if (length(grep("^(A|T|G|C)+$",mut$ref))==1 & mut$ref == delseq | nchar(mut$ref)==0) {# at 1st, check if ref seq is consistent with given ref information
    cat("Refseq:",mut$ref ,"  Altseq:",mut$alt,"\n")
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
    cat("Warning. Given ref seq in vcf:",toupper(mut$ref),"is different from actual cDNA seq:",toupper(delseq),"\nAnalysis continued, though.\n")
  }
}
return(list(newseq=newseq,pos=newcDNApos))
}

makeseqframe<-function(sentence) { #make the sentence to lseq-type dataframe
  preexon<-0
  exon<--9
  re<-data.frame()
  cDNApointer<-1
  for (i in 1:nchar(sentence)) {
    chr<-substr(sentence,i,i)
    exon<--9
    cDNApos<--9
    if (chr %in% c("A","T","G","C")) {
      exon<-preexon+1
      cDNApos<-cDNApointer
      cDNApointer<-cDNApointer+1
    } else if (i>1) if (re$exon[i-1]>0) {
      preexon<-1
    }
    
    temp<-data.frame(ID="Givenseq",chr=1,phypos=i,cDNApos=cDNApos,exon=exon,Seq=chr)
    re<-rbind(re,temp)
  }
  return(re)
}

takeresult5<-function(refresult,altresult,mode) {  #mode 0-> new5ss   1->broken canonical5ss
  #      new<-data.frame(phypos=newseq$phypos[j],MaxEntScoreSS5=output$Maxent,SS5score=output$SSscore,SS5cDNAboundary=cDNAboundary,SS5phyposboundary=phyposboundary,SS5seq=output$SSseq)
  #write.table(refresult,file="temp.refresult.txt",quote=F,sep="\t",col.names=T,row.names=F)
  #write.table(altresult,file="temp.altresult.txt",quote=F,sep="\t",col.names=T,row.names=F)
  ss5output<-0
  if (mode==0) {
    #if (nrow(altresult[altresult$MaxEntScoreSS5>4,])>0) { #in case of creating 5ss, 1st pass
      pos<-which(altresult$MaxEntScoreSS5==max(altresult$MaxEntScoreSS5)) 
      cat("Peak Alt MaxEnt Score Pos:",pos," Alt score:",altresult$MaxEntScoreSS5[pos], "Ref score:",refresult$MaxEntScoreSS5[pos],"\n")
      if ((refresult$MaxEntScoreSS5[pos]-altresult$MaxEntScoreSS5[pos])<0) {#in case of creating 5ss, determine
        ss5output<-1
      } else {
        cat("No Splice Site Creating Peak found in MaxEnt Score.\n")
      } 
    #}
  } else if (mode==1) {
    if (length(grep("c.\\d+--9",refresult$SS5cDNAboundary))>0) { #in case of canonical exon-intron boundary involved, 
      pos<-grep("c.\\d+--9",refresult$SS5cDNAboundary)
      #write.table(refresult,file="temp.debug.txt",row.names=F,col.names=T,sep="\t",quote=F)
      cat("Exon-Intron Boundary MaxEnt Score Pos:",pos," Alt score:",altresult$MaxEntScoreSS5[pos], "Ref score:",refresult$MaxEntScoreSS5[pos],"\n")
      if ((refresult$MaxEntScoreSS5[pos]-altresult$MaxEntScoreSS5[pos])>0) {#in case of broken 5ss, determine
        ss5output<-1
      } else {
        cat("The Score in the Exon-Intron boundary does not meet the canonical splice site broken criteria.\n")
      }
    } else {
      cat("No Exon-Intron boundary found in the region MaxEnt Scanned.\n")
    }
  }
  if (ss5output) {
    return(data.frame(pos=I(refresult$phypos[pos]),Ref_MaxEntScoreSS5=refresult$MaxEntScoreSS5[pos],
                      Ref_SS5score=refresult$SS5score[pos],Alt_MaxEntScoreSS5=altresult$MaxEntScoreSS5[pos],
                      Alt_SS5score=altresult$SS5score[pos],Diff_MaxEntScoreSS5=(altresult$MaxEntScoreSS5[pos]-refresult$MaxEntScoreSS5[pos]),
                      Diff_SS5score=(altresult$SS5score[pos]-refresult$SS5score[pos]),Boundary=refresult$SS5phyposboundary[pos]  ))
  } else {
    return(data.frame(pos=NA,Ref_MaxEntScoreSS5=NA,Ref_SS5score=NA,Alt_MaxEntScoreSS5=NA,Alt_SS5score=NA,Diff_MaxEntScoreSS5=NA,Diff_SS5score=NA,Boundary=NA))
  }
}

takeresult3<-function(refresult,altresult,mode) {  #mode 0-> new5ss   1->broken canonical5ss
  #      new<-data.frame(phypos=newseq$phypos[j],MaxEntScoreSS5=output$Maxent,SS5score=output$SSscore,SS5cDNAboundary=cDNAboundary,SS5phyposboundary=phyposboundary,SS5seq=output$SSseq)
  #write.table(refresult,file="temp.refresult.txt",quote=F,sep="\t",col.names=T,row.names=F)
  #write.table(altresult,file="temp.altresult.txt",quote=F,sep="\t",col.names=T,row.names=F)
  #data.frame(phypos=newseq$phypos[j],MaxEntScoreSS3=output$Maxent,SS3score=output$SSscore,SS3cDNAboundary=cDNAboundary,SS3phyposboundary=phyposboundary,SS3seq=output$SSseq)
  ss3output<-0
  if (mode==0) {
    #if (nrow(altresult[altresult$MaxEntScoreSS5>4,])>0) { #in case of creating 5ss, 1st pass
    pos<-which(altresult$MaxEntScoreSS3==max(altresult$MaxEntScoreSS3)) 
    if (length(grep("c.-9-\\d+",as.character(refresult$SS3cDNAboundary[pos])))==0) { 
        cat("Peak Alt MaxEnt Score Pos in the Seq:",pos," Alt score:",altresult$MaxEntScoreSS3[pos], "Ref score:",refresult$MaxEntScoreSS3[pos],"\n")
        if ((refresult$MaxEntScoreSS3[pos]-altresult$MaxEntScoreSS3[pos])<0) {#in case of creating 3ss, determine
          ss3output<-1
        } else {
          cat("No Splice Site Creating Peak found in MaxEnt Score.\n")
        } 
    } else {
      cat("Peak Position is on an Intron-Exon Boundary:",as.character(refresult$SS3cDNAboundary[pos]),"\n")
    }
  } else if (mode==1) {
    if (length(grep("c.-9-\\d+",as.character(refresult$SS3cDNAboundary)))>0) { #in case of canonical exon-intron boundary involved, 
      pos<-grep("c.-9-\\d+",as.character(refresult$SS3cDNAboundary))
      #write.table(refresult,file="temp.debug.txt",row.names=F,col.names=T,sep="\t",quote=F)
      cat("Intron-Exon Boundary:",as.character(refresult$SS3cDNAboundary[pos]) ,"MaxEnt Score Pos in the Seq:",pos," Alt score:",altresult$MaxEntScoreSS3[pos], "Ref score:",refresult$MaxEntScoreSS3[pos],"\n")
      if ((refresult$MaxEntScoreSS3[pos]-altresult$MaxEntScoreSS3[pos])>0) {#in case of broken 3ss, determine
        ss3output<-1
      } else {
        cat("The Score in the Intron-Exon boundary does not meet the canonical splice site broken criteria.\n")
      }
    } else {
      cat("No Intron-Exon boundary found in the region MaxEnt Scanned.\n")
    }
  }
  if (ss3output) {
    return(data.frame(pos=I(refresult$phypos[pos]),Ref_MaxEntScoreSS3=refresult$MaxEntScoreSS3[pos],
                      Ref_SS3score=refresult$SS3score[pos],Alt_MaxEntScoreSS3=altresult$MaxEntScoreSS3[pos],
                      Alt_SS3score=altresult$SS3score[pos],Diff_MaxEntScoreSS3=(altresult$MaxEntScoreSS3[pos]-refresult$MaxEntScoreSS3[pos]),
                      Diff_SS3score=(altresult$SS3score[pos]-refresult$SS3score[pos]),Boundary=refresult$SS3phyposboundary[pos]  ))
  } else {
    return(data.frame(pos=NA,Ref_MaxEntScoreSS3=NA,Ref_SS3score=NA,Alt_MaxEntScoreSS3=NA,Alt_SS3score=NA,Diff_MaxEntScoreSS3=NA,Diff_SS3score=NA,Boundary=NA))
  }
}




takeseq<-function(seq,pos,len1,len2, mode) {
  exonrev<-seq$exon[pos] #actual exon 
  exon<-exonrev  #if mut is on an intron, search nearby exon
  temppos<-pos
  #cat ("given exon:",exon,"\t")
  if (exon>0 & IN$ss3==1) exon<-exon-1 #because mut is on 2nd exon in case of 3ss  ##here, exon should be set to 1st exon 
  if (exon<0) {  #if mut on intron , search 1st exon in case of both 5ss and 3ss
    while (exon<0) {
      temppos<-temppos-1
      exon<-seq$exon[temppos]
    }
  } 
  preexon<-exon
  chars<-""
  reseq<-""
  #cat ("Searched exon:",exon,"\tmode",mode,"Len1 Len2:",len1,len2,"\t")
  #write.table(seq,file="temp.takeseq.txt",quote=F,sep="\t",)
  if (mode==1) {#take seq in the 2nd exon
    if (IN$limit2ndexon!=-9) if (len1>IN$limit2ndexon) len1<-IN$limit2ndexon
    startpos<-min(which(seq$exon==(exon+1)))
    endpos<-max(which(seq$exon==(exon+1)))
    exon2seq<-paste(seq$Seq[startpos:endpos],collapse="")
    if (nchar(exon2seq)<=len1) {
      reseq<-exon2seq
      aclen[6]<<-nchar(reseq)
    } else {
      reseq<-substr(exon2seq,1,len1)
    }
    aclen[6]<<-nchar(reseq)
  } else if (mode==2) { #take intron  
    startpos<-max(which(seq$exon==exon))  #pos of the last nucleotide in the 1st exon
    endpos<-min(which(seq$exon==(exon+1)))  #pos of the 1st nucleotide in the 2nd exon
    intronseq<-paste(seq$Seq[(startpos+1):(endpos-1)],collapse="")      #take all intron         
    if (nchar(intronseq)<=(len1+len2)) { #if intron is shorter than indicated, take all
      reseq<-intronseq
      if (IN$ss5==1) {aclen[4]<<-nchar(intronseq);aclen[5]<<-0}
      if (IN$ss3==1) {aclen[2]<<-nchar(intronseq);aclen[3]<<-0}  #just in case of mut is on exon 
    } else { #take first len1 and the last len2
      reseq<-paste(substr(intronseq,1,len1),substr(intronseq,nchar(intronseq)-len2+1,nchar(intronseq)),sep="")
      if (IN$ss5==1) {aclen[4]<<-len1;aclen[5]<<-len2}
      if (IN$ss3==1) {aclen[2]<<-len1;aclen[3]<<-len2} #just in case of mut is on exon 
    } 
  } else if (mode==3) { #take seq after mutation  #include the mutation!
    exonseqall<-""
    #take all exon after mutation
    for (i in pos:nrow(seq)) {
      if (seq$exon[i]==-9) break
      exonseqall<-paste(exonseqall,seq$Seq[i],sep="")
    }
    if (nchar(exonseqall)<=(len1+len2)) { #if exon after mutation is shorter than indicated
      reseq<-exonseqall
      aclen[2]<<-nchar(exonseqall)
      aclen[3]<<-0
    } else { #take first len1 and the last len2
      reseq<-paste(substr(exonseqall,1,len1),substr(exonseqall,nchar(exonseqall)-len2+1,nchar(exonseqall)),sep="")
      aclen[2]<<-len1
      aclen[3]<<-len2
    }
    
  } else if (mode==4) { #take seq before mutation #does not include the mutation
    for (i in 1:len1) {
      if ((pos-i)==0) break
      nexon<-seq$exon[pos-i]
      if (nexon!=exon) break
      reseq<-paste(reseq, seq$Seq[pos-i],sep="")
      aclen[1]<<-nchar(reseq)
    }
    temp<-unlist(strsplit(reseq,""))
    reseq<-paste(temp[length(temp):1],collapse="")
  } else if (mode==5) { #take intron before mut ** include the mutation
    for (i in (temppos+1):pos) {
      reseq<-paste(reseq, seq$Seq[i],sep="")
    }
    aclen[3]<<-nchar(reseq)
  } else if (mode==6) { #take intron after mut , intron attached to exon2  in case that mut is on intron
    endpos<-min(which(seq$exon==(exon+1)))-1
    if (endpos==pos) {  #in case of 0bp long in intron after mut (not incl mut)
      aclen[4]<<-0;aclen[5]<<-0
      return("")
    }
    intronseq<-paste(seq$Seq[(pos+1):endpos],collapse="")
    #cat("pos:",pos,"endpos:",endpos,"Intronseq:",intronseq,"\n")
    #write.table(seq,file="temp.intronseq.seq.txt",sep="\t",quote=F,col.names=F,row.names=F)
    if (length(grep("NA",intronseq))>0) stop("NA intron\n")
    if (nchar(intronseq)<=(len1+len2)) { #if intron is shorter than indicated, take all
      reseq<-intronseq
      aclen[4]<<-nchar(intronseq)
      aclen[5]<<-0
    } else { #take first len1 and the last len2
      reseq<-paste(substr(intronseq,1,len1),substr(intronseq,nchar(intronseq)-len2+1,nchar(intronseq)),sep="")
      aclen[4]<<-len1
      aclen[5]<<-len2
    } 
  } else if (mode==7) {  #take 1st exon just before mut
    #cat("temppos:",temppos, "len1:",len1,"\n")
    start<-(temppos-len1+1)
    if (start<1) start<-1
      for (i in start:temppos) {
        if (seq$exon[i]==exon) {
            reseq<-paste0(reseq,seq$Seq[i])
        }
     }
     aclen[1]<<-nchar(reseq)
  } else if (mode==8) { #take 1st exon for 3SS
    if (IN$limit1stexon!=-9) if (len1>IN$limit1stexon) len1<-IN$limit1stexon
    exon1stseq<-paste(seq[seq$exon==exon,]$Seq,collapse="")#here take entire 1st exon
    #cat("exon1stseq",exon1stseq,"len1:",len1,"\n")
    reseq<-substr(exon1stseq,nchar(exon1stseq)-len1+1,nchar(exon1stseq))  
    aclen[1]<<-nchar(reseq)
  } else if (mode==9) {  #take 2nd exon before mut (incl mut) and 2nd exon after canonical spice site for 3SS # mut should be on 2nd exon!
    startpos<-min(which(seq$exon==(exon+1)))
    exon2nd<-paste(seq$Seq[startpos:pos],collapse="")
    if (nchar(exon2nd)<=(len1+len2)) {  #if within the limit length, take all
      reseq<-exon2nd
      aclen[4]<<-nchar(reseq);aclen[5]<<-0
    } else {   #if greater than  the limit length, take the left side and the right side respectively 
      reseq<-paste0(substr(exon2nd,1,len1),substr(exon2nd,nchar(exon2nd)-len2+1,nchar(exon2nd)))
      aclen[4]<<-len1;aclen[5]<<-len2
    }
  } else if (mode==10) {  #take exon2nd after mut (not incl ,mut) for 3ss
    endpos<-max(which(seq$exon==(exon+1)))
    #cat ("pos;",pos,"endpos",endpos,"\n")
    exon2nd<-paste(seq$Seq[(pos+1):endpos],collapse="") #take all exon after mut
    
    if (nchar(exon2nd)<len1) {  #if actual exon is shorter than indicated
      reseq<-exon2nd
    } else {
      reseq<-substr(exon2nd,1,len1)
    }
    aclen[6]<<-nchar(reseq)
  } else if (mode==11) { #take 1st intron+intron before mut (incl mut) for 3ss in case that mut in on intron
    statpos<-max(which(seq$exon==exon))+1
    intron<-paste(seq$Seq[statpos:pos],collapse="") #take all intron before mut (incl mut)
    if (nchar(intron)<=(len1+len2)) {  #if actual intron before mut is shorter than indicated
      reseq<-intron
      aclen[2]<<-nchar(reseq);aclen[3]<<-0
    } else {
      reseq<-paste0(substr(intron,1,len1),substr(intron,nchar(intron)-len2+1,nchar(intron)))
      aclen[2]<<-len1;aclen[3]<<-len2
    }
  }
  #cat("Mode:",mode,"refseq:",reseq,"\n")
  return(reseq)
}


  
#refconst<-constructmake(lseq,newcDNApos[1])
constructmake<-function(seq,pos,mut,conlen) { #make construcrt according to the indication
  if (IN$ss5==1) {
    exon<-seq$exon[pos]
    temppos<-pos
    #cat("pos:",pos,"exon:",seq$exon[pos],"\n")
    if (exon<0) {
      while (exon<0) {
        temppos<-temppos-1
        if (temppos==0) {
          cat("The given mutation is located before 1st exon. Therefore 3SS construct can NOT be designed because of 1st exon not being able to be assigned.\n")
          return (data.frame(Var_Name=NA )  )
        }
        exon<-seq$exon[temppos]
      }
    } 
    if (exon==max(seq$exon)) {
      cat("The given mutation is located on the last exon. Therefore 5SS construct can NOT be designed because of 2nd exon not being able to be assigned.\n")
      return (data.frame(Var_Name=NA )  )
    }
    
    #cat("Construct Length:",conlen,"\n")
    limit<-sum(conlen)
    aclen<<-c(0,0,0,0,0,0) # length for mutation on exon  #mut on exon new splice site weighted    exon1before,exon1after(+mut),intron,exon2nd
    temp<-as.integer( (conlen[1]+conlen[3]+conlen[4])/2)
    tempfinal<-conlen[1]+conlen[3]+conlen[4]-temp
    conlencsds<-c(temp,conlen[2],0,tempfinal,conlen[5],conlen[6]) #length for mutataion on exon and weight canonical splice donor site
    aclencsds<<-c(0,0,0,0,0,0)                 #1stExonbeforeMut,1stExonAfterMut,0(blank),intron attached to 1st exon, intron attached to 2nd exon, 2nd exon
    #get exon seq which the mutation resides in
    #if exon2nd < conlen[5] _> the remaing length will move to intron1st
    exon2nd<-takeseq(seq,pos,conlen[6],0,1)
    remain<-conlen[6]-nchar(exon2nd)
    conlen[6]<-aclen[6]
    aclencsds[6]<<-aclen[6]
    exon2ndcsds<-exon2nd
    conlencsds[6]<-aclencsds[6]
    if (seq$exon[pos]>0) { #when mut is on an exon
      if (remain>0) {
        quarter<-as.integer(remain/4)
        quarterfinal<-remain-quarter*3
        triplet<-as.integer(remain/3)
        tripletfinal<-remain-triplet*2
        #cat("#bp of Remain:",remain,"\t Quarter:",quarter,"\tQuarterFinal:",quarterfinal,"\n")
        conlen[4]<-conlen[4]+quarter;      conlen[3]<-conlen[3]+quarter;      conlen[2]<-conlen[2]+quarter;    conlen[1]<-conlen[1]+quarterfinal
        conlencsds[1]<-conlencsds[1]+triplet; conlencsds[2]<-conlencsds[2]+triplet; conlencsds[4]<-conlencsds[4]+tripletfinal;
      }
      #if intron1 + intron2 < conlen[3]+conlen[4]  -> the reamining length will move to exon1 region
      
      # canonical splice site weighted part start
      aclenstore<-aclen #store aclen because csds make use of the same subroutine of aclen and change the parameter
      exon1aftercsds<-takeseq(seq,pos,conlencsds[2],0,3) #include mutation #exon after mut and exon before canonical splice site. assume that mut is located near canonical splice site.
      aclencsds[2]<<-aclen[2];       remaincsds<-conlencsds[2]-aclencsds[2]
      if (remaincsds>0) {
        double<-as.integer(remaincsds/2)
        doublefinal<-remaincsds-double
        conlencsds[1]<-conlencsds[1]+double; conlencsds[4]<-conlencsds[4]+doublefinal
      }
      temp<-limitIntron(c(conlencsds[4],conlencsds[5]));conlencsds[4]<-temp[1];conlencsds[5]<-temp[2]
      introncsds<-takeseq(seq,pos,conlencsds[4],conlencsds[5],2)
      aclencsds[4]<<-aclen[4];aclencsds[5]<<-aclen[5]
      remaincsds<-conlencsds[4]-aclencsds[4]+conlencsds[5]-aclencsds[5]
      conlencsds[4]<-aclencsds[4];conlencsds[5]<-aclencsds[5]
      if (remaincsds>0) {
        conlencsds[1]<-conlencsds[1]+remaincsds
      }
      exon1beforecsds<-takeseq(seq,pos,conlencsds[1],0,4)
      aclencsds[1]<<-aclen[1];conlencsds[1]<-aclencsds[1]
      remaincsds<-limit-nchar(exon1beforecsds)-nchar(exon1aftercsds)-nchar(introncsds)-nchar(exon2ndcsds)
      if (remaincsds>0) {
        conlencsds[4]<-conlencsds[4]+remaincsds
        if (IN$limitintron!=-9) if ((conlencsds[4]+conlencsds[5])>IN$limitintron) {conlencsds[4]<-IN$limitintron-conlencsds[5]}
        introncsds<-takeseq(seq,pos,conlencsds[4],conlencsds[5],2)
        aclencsds[4]<<-aclen[4];aclencsds[5]<<-aclen[5];conlencsds[4]<-aclencsds[4];conlencsds[5]<-aclencsds[5]
        remaincsds<-limit-nchar(exon1beforecsds)-nchar(exon1aftercsds)-nchar(introncsds)-nchar(exon2ndcsds)
        if (remaincsds>0) {
          double<-as.integer(remaincsds/2)
          doublefinal<-remaincsds-double
          conlencsds[5]<-conlencsds[5]+double; conlencsds[6]<-conlencsds[6]+doublefinal
          if (IN$limitintron!=-9) if ((conlencsds[4]+conlencsds[5])>IN$limitintron) {conlencsds[5]<-IN$limitintron-conlencsds[4]}
          introncsds<-takeseq(seq,pos,conlencsds[4],conlencsds[5],2)
          exon2ndcsds<-takeseq(seq,pos,conlencsds[6],0,1)
          aclencsds[4]<<-aclen[4];aclencsds[5]<<-aclen[5];aclencsds[6]<<-aclen[6];conlencsds[4]<-aclencsds[4];conlencsds[5]<-aclencsds[5];conlencsds[6]<-aclencsds[6]
          remaincsds<-limit-nchar(exon1beforecsds)-nchar(exon1aftercsds)-nchar(introncsds)-nchar(exon2ndcsds)
          if (remaincsds>0) {
            conlencsds[6]<-conlencsds[6]+remaincsds
            exon2ndcsds<-takeseq(seq,pos,conlencsds[6],0,1)
            aclencsds[6]<<-aclen[6];conlencsds[6]<-aclencsds[6]
         }
        }
      }
      aclen<<-aclenstore  # canonical splice site weighted part end 
      temp<-limitIntron(c(conlen[4],conlen[5]));conlen[4]<-temp[1];conlen[5]<-temp[2]
      intron<-takeseq(seq,pos,conlen[4],conlen[5],2)  #intron 5 part and 3 part. if actual intron is shorter than indicated , merge and output
      remain<-conlen[4]+conlen[5]-nchar(intron)
      conlen[4]<-aclen[4];conlen[5]<-aclen[5]
      if (remain>0) {
        quarter<-as.integer(remain/3)
        quarterfinal<-remain-quarter*2
        #cat("#bp of Remain:",remain,"\t Quarter:",quarter,"\tQuarterFinal:",quarterfinal,"\n")
        conlen[3]<-conlen[3]+quarter;      conlen[2]<-conlen[2]+quarter
        conlen[1]<-conlen[1]+quarterfinal
      }
      #if exon after ends within conlen[2] limit, the remaining length will move to exon before.
      exon1after<-takeseq(seq,pos,conlen[2],conlen[3],3)  #include mutation #exon after mut and exon before canonical splice site. if actual exon after mut and exon before cano splice site is shorter than indicated, then merged and output
      conlen[2]<-aclen[2];conlen[3]<-aclen[3]
      conlen[1]<-limit-nchar(exon1after)-nchar(intron)-nchar(exon2nd)
      exon1before<-takeseq(seq,pos,conlen[1],0,4)
      conlen[1]<-aclen[1]
      remain<-limit-nchar(exon1before)-nchar(exon1after)-nchar(intron)-nchar(exon2nd)
      if (remain>0) {   #remainings redistribution
        quarter<-as.integer(remain/3)
        quarterfinal<-remain-quarter*2
        #cat("#bp of Remain:",remain,"\t Quarter:",quarter,"\tQuarterFinal:",quarterfinal,"\n")
        conlen[4]<-conlen[4]+quarter;      conlen[3]<-conlen[3]+quarter;
        conlen[2]<-conlen[2]+quarterfinal
        exon1after<-takeseq(seq,pos,conlen[2],conlen[3],3)
        conlen[2]<-aclen[2];conlen[2]<-aclen[2]
        if (IN$limitintron!=-9) if ((conlen[4]+conlen[5])>IN$limitintron) {conlen[4]<-IN$limitintron-conlen[5]}
        intron<-takeseq(seq,pos,conlen[4],conlen[5],2)
        conlen[4]<-aclen[4];conlen[5]<-aclen[5]
        remain<-limit-nchar(exon1before)-nchar(exon1after)-nchar(intron)-nchar(exon2nd)
        
        if (remain>0) {
          quarter<-as.integer(remain/2)
          quarterfinal<-remain-quarter
          #cat("#bp of Remain:",remain,"\t Quarter:",quarter,"\tQuarterFinal:",quarterfinal,"\n")
          conlen[6]<-conlen[6]+quarter;
          conlen[5]<-conlen[5]+quarterfinal
          if (IN$limitintron!=-9) if ((conlen[4]+conlen[5])>IN$limitintron) {conlen[5]<-IN$limitintron-conlen[4]}
          intron<-takeseq(seq,pos,conlen[4],conlen[5],2)
          exon2nd<-takeseq(seq,pos,conlen[6],0,1)
          conlen[4]<-aclen[4];conlen[5]<-aclen[5];conlen[6]<-aclen[6]
          remain<-limit-nchar(exon1before)-nchar(exon1after)-nchar(intron)-nchar(exon2nd)
          if (remain>0) {
            conlen[6]<-conlen[6]+remain
            #cat("#bp of Remain:",remain," conlen6",conlen[6]," aclen6",aclen[6],"exon1before",exon1before,"exon1after",exon1after,"intron","\n")
            exon2nd<-takeseq(seq,pos,conlen[6],0,1)
            conlen[6]<-aclen[6]
          }
        }
      }
    } else { #when mut is on an intron
      #cat("I am here\n")
      conlen<-c(tempfinal,0,conlen[2],temp,conlen[5],conlen[6]) #length for mutataion on intron and weight canonical splice donor site
      #aclen<<-c(0,0,0,0,0,0)                 #1stExon,0(dummy),intron between 1st exon and mut,intron after mut, intron attached to 2nd exon, 2nd exon
      if (remain>0) {
        quarter<-as.integer(remain/3)
        quarterfinal<-remain-quarter*2
        #cat("#bp of Remain:",remain,"\t Quarter:",quarter,"\tQuarterFinal:",quarterfinal,"\n")
        conlen[4]<-conlen[4]+quarter;conlen[3]<-conlen[3]+quarter;
        conlen[1]<-conlen[1]+quarterfinal
      }
      #if intron1 + intron2 < conlen[3]+conlen[4]  -> the reamining length will move to exon1 region
      intronbeforemut<-takeseq(seq,pos,conlen[3],0,5)  #intron 5 part before mut and include mut
      remain<-conlen[3]-aclen[3]
      conlen[3]<-aclen[3]
      if (remain>0) {
        quarter<-as.integer(remain/2)
        quarterfinal<-remain-quarter
        #cat("#bp of Remain:",remain,"\t Quarter:",quarter,"\tQuarterFinal:",quarterfinal,"\n")
        conlen[4]<-conlen[4]+quarter
        conlen[1]<-conlen[1]+quarterfinal
      }
      intronaftermut<-takeseq(seq,pos,conlen[4],conlen[5],6) #intron5 part after mut (NOT include mut) and intron3 part
      
      temp<-limitIntron(c(conlen[3],conlen[4],conlen[5]));conlen[3]<-temp[1];conlen[4]<-temp[2];conlen[5]<-temp[3]
      intronbeforemut<-takeseq(seq,pos,conlen[3],0,5) 
      intronaftermut<-takeseq(seq,pos,conlen[4],conlen[5],6)
      conlen[3]<-aclen[3];conlen[4]<-aclen[4];conlen[5]<-aclen[5]
      
      #cat("I am here\n")
      remain<-limit-nchar(intronbeforemut)-nchar(intronaftermut)-nchar(exon2nd)
      #conlen[4]<-aclen[4];conlen[5]<-aclen[5]
      if (remain>0) {
        conlen[1]<-conlen[1]+remain
      }
      exon1<-takeseq(seq,pos,conlen[1],0,7)  #just 1st exon before mut
      conlen[1]<-aclen[1]
      #cat("here1 conlen:",conlen,"\naclen:",aclen,"\n")
      remain<-limit-nchar(exon1)-nchar(intronbeforemut)-nchar(intronaftermut)-nchar(exon2nd)
      if (remain>0) {   #remainings redistribution
        conlen[4]<-conlen[4]+remain
        if (IN$limitintron!=-9) if ((conlen[3]+conlen[4]+conlen[5])>IN$limitintron) {conlen[4]<-IN$limitintron-conlen[3]-conlen[5]}
        intronaftermut<-takeseq(seq,pos,conlen[4],conlen[5],6) 
        conlen[4]<-aclen[4];conlen[5]<-aclen[5]
        #cat("here2 conlen:",conlen,"\naclen:",aclen,"\n")
        remain<-limit-nchar(exon1)-nchar(intronbeforemut)-nchar(intronaftermut)-nchar(exon2nd)
        if (remain>0) {
          quarter<-as.integer(remain/2)
          quarterfinal<-remain-quarter
          #cat("#bp of Remain:",remain,"\t Quarter:",quarter,"\tQuarterFinal:",quarterfinal,"\n")
          conlen[6]<-conlen[6]+quarter;
          conlen[5]<-conlen[5]+quarterfinal
          if (IN$limitintron!=-9) if ((conlen[3]+conlen[4]+conlen[5])>IN$limitintron) {conlen[5]<-IN$limitintron-conlen[3]-conlen[4]}
          intronaftermut<-takeseq(seq,pos,conlen[4],conlen[5],6) 
          exon2nd<-takeseq(seq,pos,conlen[6],0,1)
          conlen[4]<-aclen[4];conlen[5]<-aclen[5];conlen[6]<-aclen[6]
          #cat("here3 conlen:",conlen,"\naclen:",aclen,"\n")
          remain<-limit-nchar(exon1)-nchar(intronbeforemut)-nchar(intronaftermut)-nchar(exon2nd)
          if (remain>0) {
            conlen[6]<-conlen[6]+remain
            #cat("#bp of Remain:",remain," conlen6",conlen[6]," aclen6",aclen[6],"exon1before",exon1before,"exon1after",exon1after,"intron","\n")
            exon2nd<-takeseq(seq,pos,conlen[6],0,1)
            conlen[6]<-aclen[6]
            #cat("here4 conlen:",conlen,"\naclen:",aclen,"\n")
          }
        }
      }
      #cat("finale conlen:",conlen,"\naclen:",aclen,"\n")
    }
    #assemble entire construct
    #mut on exon new splice site weighted    exon1before,exon1after(+mut),intron,exon2nd
    #mut on exon canonical splice donor site weighted exon1beforecsds,exon1aftercsds(+mut),introncsds,exon2ndcsds
    #mut on intron, exon1,intronbeforemut(+mut), intronaftermut,exon2nd
    #if ((seq$exon[pos]>0) & ((nchar(exon2nd)*nchar(intron)*nchar(exon1after)*nchar(exon1before)==0) | (nchar(exon2ndcsds)*nchar(introncsds)*nchar(exon1aftercsds)*nchar(exon1beforecsds)==0) )) {
    #  cat("Zero-length component found.\nNew Splice Site Weighte Construct   Exon before the mutation:",nchar(exon1before),"Exon after the mutation:",nchar(exon1after),"Intron:",nchar(intron),"2nd Exon:",nchar(exon2nd),
    #  "\n Canonical Splice Donor Site Weighted Construct   Exon before the mutation:",nchar(exon1beforecsds),"Exon after the mutation:",nchar(exon1aftercsds),"Intron:",nchar(introncsds),"2nd Exon:",nchar(exon2ndcsds),"\n")
    #  return (data.frame(Var_Name=NA )  )  #if any of elements lacks, null returned
    #} else if ((seq$exon[pos]<0) & (nchar(exon2nd)*nchar(intronaftermut)*nchar(intronbeforemut)*nchar(exon1)==0)) {
    #  cat("Zero-length component found.\n1st Exon:",nchar(exon1),"Intron before the mutation:",nchar(intronbeforemut),"Intron after the mutation:",nchar(intronaftermut),"2nd Exon:",nchar(exon2nd),"\n")
    #  return (data.frame(Var_Name=NA )  )  #if any of elements lacks, null returned
    #
  
    #barcode embed
    #barcodesYES<-0  #if employ barcode or not
    #barcodes<-c()  #barcodes to be employed
    #barcodesPOS<-1   #pointer in the barcodes
    #temp<-regexpr("\\d+",mut.change[i,]$Var_5SplicePoint) #instead of variant position, actual 5ss gap taken
    #x<-as.numeric(substr(mut.change[i,]$Var_5SplicePoint,as.numeric(temp),as.numeric(temp)+attr(temp,"match.length")-1))
    #
    refseq5<-IN$seq5;    refseq3<-IN$seq3
    refcurrentbarcode<-NA
    altseq5<-IN$seq5;    altseq3<-IN$seq3
    altcurrentbarcode<-NA
    
    if (barcodesYES) {
      barcodesPOSadd<-0;
      barcodesPOS2<-barcodesPOS+1
      if (barcodesPOS2>length(barcodes)) barcodesPOS2<-1
      temp3<-regexpr("(X)+",IN$seq3)
      temp3s<-as.numeric(temp3);temp3e<-temp3s+attr(temp3,"match.length")-1
      if (temp3s>0) {
        refseq3<-paste(substr(IN$seq3,1,temp3s-1),barcodes[barcodesPOS],substr(IN$seq3,temp3e+1,nchar(IN$seq3)),sep="")
        altseq3<-paste(substr(IN$seq3,1,temp3s-1),barcodes[barcodesPOS2],substr(IN$seq3,temp3e+1,nchar(IN$seq3)),sep="")
        barcodesPOSadd<-2
        refcurrentbarcode<-barcodes[barcodesPOS]
        altcurrentbarcode<-barcodes[barcodesPOS2]
      }
      
      temp5<-regexpr("(X)+",IN$seq5)
      temp5s<-as.numeric(temp5);temp5e<-temp5s+attr(temp5,"match.length")-1
      if (temp5s>0) {
        refseq5<-paste(substr(IN$seq5,1,temp5s-1),barcodes[barcodesPOS],substr(IN$seq5,temp5e+1,nchar(IN$seq5)),sep="")
        altseq5<-paste(substr(IN$seq5,1,temp5s-1),barcodes[barcodesPOS2],substr(IN$seq5,temp5e+1,nchar(IN$seq5)),sep="")
        barcodesPOSadd<-2
        refcurrentbarcode<-barcodes[barcodesPOS]
        altcurrentbarcode<-barcodes[barcodesPOS2]
      }
      barcodesPOS<<-barcodesPOS+barcodesPOSadd
      if (barcodesPOS>length(barcodes)) {
        barcodesPOS<<-barcodesPOS-length(barcodes)
        cat("barcode cycle returned to the initial\n")
      }
    }
    
    if (seq$exon[pos]>0) {
      #cat("here1\n")
      returnseq<-paste(tolower(refseq5),toupper(exon1before),toupper(exon1after),tolower(intron),toupper(exon2nd),tolower(refseq3),sep="")
      returnmutseq<-paste(tolower(altseq5),toupper(exon1before),toupper(exon1after),tolower(intron),toupper(exon2nd),tolower(altseq3),sep="")
      mutpos<-nchar(refseq5)+nchar(exon1before)+1
      
      returnseqcsds<-paste(tolower(refseq5),toupper(exon1beforecsds),toupper(exon1aftercsds),tolower(introncsds),toupper(exon2ndcsds),tolower(refseq3),sep="")
      returnmutseqcsds<-paste(tolower(altseq5),toupper(exon1beforecsds),toupper(exon1aftercsds),tolower(introncsds),toupper(exon2ndcsds),tolower(altseq3),sep="")
      mutposcsds<-nchar(refseq5)+nchar(exon1beforecsds)+1
    } else {
      #cat("here2\n")
      returnseq<-paste(tolower(refseq5),toupper(exon1),tolower(intronbeforemut),tolower(intronaftermut),toupper(exon2nd),tolower(refseq3),sep="")
      returnmutseq<-paste(tolower(altseq5),toupper(exon1),tolower(intronbeforemut),tolower(intronaftermut),toupper(exon2nd),tolower(altseq3),sep="")
      mutpos<-nchar(refseq5)+nchar(exon1)+nchar(intronbeforemut)
    }
    #check score in ref seq
    cat("Checking the potential of 5\'splice site  in the referrence construct...\n")
    Start<-mutpos-9+1
    if (Start<1) Start<-1  #indicate the region where the mutation affects
    End <-mutpos
    seqtemp<-makeseqframe(returnseq)
    #cat("start;",Start,"end;",End, "\n")
    #write.table(seq,file="temp.seq.txt",sep="\t",quote=F,col.names=T,row.names=F)
    refresult5temp<-SS5scan(seqtemp,Start,End)
    #make mut allele  ###in this ver. just SNV permitted
    if (toupper(mut$ref)==toupper(substr(returnmutseq,mutpos,mutpos))) {
      if (seq$exon[pos]>0) {
        returnmutseq<-paste(substr(returnmutseq,1,mutpos-1),toupper(mut$alt),substr(returnmutseq,mutpos+1,nchar(returnmutseq)),sep="")
      } else {
        returnmutseq<-paste(substr(returnmutseq,1,mutpos-1),tolower(mut$alt),substr(returnmutseq,mutpos+1,nchar(returnmutseq)),sep="")
      }
      
    } else {
      cat("Ref allele not found in the actual sequence.\n Ref in the mutation file:",mut$ref,"Ref in the actual sequence:",substr(returnmutseq,mutpos,mutpos),"\n")
      cat ("Sequence in Reference Construct:",substr(returnseq,1,mutpos-1),"[",substr(returnseq,mutpos,mutpos),"]",substr(returnseq,mutpos+1,nchar(returnseq)),"\n",sep="")
      cat ("Sequence in Mutant Construct before MutGenesis:",substr(returnmutseq,1,mutpos-1),"[",substr(returnmutseq,mutpos,mutpos),"]",substr(returnmutseq,mutpos+1,nchar(returnmutseq)),"\n",sep="")
      return (    data.frame(Var_Name=NA )   )
    }
    seqtemp<-makeseqframe(returnmutseq)
    altresult5temp<-SS5scan(seqtemp,Start,End)

    #write.table(result5,file="temp.result5.txt",sep="\t",quote=F,col.names=T,row.names=F)
    
    #cat("ssgap:",ssgap,"\n")
    #cat("pos:",pos,"exon:",seq$exon[pos],"\n")
    if (seq$exon[pos]>0 ) {  #in case that mutation on an exon
      #cat("here1\n")
      result5<-takeresult5(refresult5temp,altresult5temp,0)
      ssgap<-as.numeric(result5$pos)+2
            
      Start<-mutposcsds-9+1;if (Start<1) Start<-1 
      End <-mutposcsds
      seqtemp<-makeseqframe(returnseqcsds)
      refresult5csdstemp<-SS5scan(seqtemp,Start,End)
      returnmutseqcsds<-paste(substr(returnmutseqcsds,1,mutposcsds-1),mut$alt,substr(returnmutseqcsds,mutposcsds+1,nchar(returnmutseqcsds)),sep="")
      seqtemp<-makeseqframe(returnmutseqcsds)
      altresult5csdstemp<-SS5scan(seqtemp,Start,End)
      result5csds<-takeresult5(refresult5csdstemp,altresult5csdstemp,1)
      returnframe<-data.frame()
      if (!is.na(result5$pos)) {      
        PredictedNormalSplice<-paste(tolower(refseq5),toupper(exon1before),toupper(exon1after),"[",tolower(intron),"]",toupper(exon2nd),tolower(refseq3),sep="")
        exon1aftermut<-paste(mut$alt,substr(exon1after,2,nchar(exon1after)),sep="")
        PredictedNormalSpliceMut<-paste(tolower(altseq5),toupper(exon1before),toupper(exon1aftermut),"[",tolower(intron),"]",toupper(exon2nd),tolower(altseq3),sep="")
        PredictedNormalLength<-nchar(returnseq)-nchar(intron)
        PredictedNormalLengthMut<-nchar(returnmutseq)-nchar(intron)
        #temp<-paste(exon2nd,altseq3,sep="")
        PredictedAberrantSplice<-paste(substr(returnseq,1,ssgap),"[",substr(returnseq,ssgap+1,nchar(paste(refseq5,exon1before,exon1after,intron,sep=""))),"]",substr(returnseq,nchar(paste(refseq5,exon1before,exon1after,intron,sep=""))+1,nchar(returnseq)),sep="")
        PredictedAberrantSpliceMut<-paste(substr(returnmutseq,1,ssgap),"[",substr(returnmutseq,ssgap+1,nchar(paste(altseq5,exon1before,exon1aftermut,intron,sep=""))),"]",substr(returnmutseq,nchar(paste(altseq5,exon1before,exon1aftermut,intron,sep=""))+1,nchar(returnmutseq)),sep="")
        PredictedAberrantLength<-PredictedNormalLength-(nchar(paste(refseq5,exon1before,exon1after))-ssgap)+2
        PredictedAberrantLengthMut<-PredictedNormalLengthMut-(nchar(paste(altseq5,exon1before,exon1aftermut))-ssgap)+2
        #(data.frame(pos=refresult$phypos[pos],Ref_MaxEntScoreSS5=refresult$MaxEntScoreSS5[pos],Ref_SS5score=refresult$SS5score[pos],Alt_MaxEntScoreSS5=altresult$MaxEntScoreSS5[pos],Alt_SS5score=altresult$SS5score[pos],Diff_MaxEntScoreSS5=(altresult$MaxEntScoreSS5[pos]-refresult$MaxEntScoreSS5[pos]),Diff_SS5score=(altresult$SS5score[pos]-refresult$SS5score[pos]) ))
        rid<-substr(paste(exon2nd,refseq3,sep=""),nchar(exon2nd)+1-nid[1],nchar(exon2nd)+nchar(barcodes[1])+nid[2])
        aid<-substr(paste(exon2nd,altseq3,sep=""),nchar(exon2nd)+1-nid[1],nchar(exon2nd)+nchar(barcodes[1])+nid[2])
        PredictedNormalSplicedOutPos<-paste0(nchar(paste0(refseq5,exon1before,exon1after))+1,"_",nchar(paste0(refseq5,exon1before,exon1after,intron)))
        PredictedAberantSplicedOutPos<-paste0(ssgap+1,"_",nchar(paste0(altseq5,exon1before,exon1aftermut,intron)))
        content<-paste("Attached5Seq:",nchar(refseq5),";1stExon_BeforeTheMut:",aclen[1],";1stExon_On_AfterTheMut:",aclen[2],";1stExon_Before_Canonical5SS:",aclen[3],";Intron_5:",aclen[4],";Intron_3:",aclen[5],";2ndExon:",aclen[6],";Attached3Seq:",nchar(refseq3),";Remain:",limit-sum(aclen),sep="")
        temp<-paste0(substr(returnseq,1,ssgap),"{",substr(returnseq,ssgap+1,nchar(paste0(refseq5,exon1before,exon1after))),"[",tolower(intron),"]}",toupper(exon2nd),tolower(IN$seq3))
        temp2<-paste0(substr(temp,1,mutpos-1+1),"(",mut$ref,"/",mut$alt,")",substr(temp,mutpos+1+1,nchar(temp)))  #give (mut/alt) info *+1 means a shift to rightside because of inserted {
        summaryseq<-paste0(tolower(IN$seq5),substr(temp2,nchar(IN$seq5)+1,nchar(temp2)))  #reset 5seq
        normalseqgrep<-substr(returnseq,mutpos-seqgrep[1],mutpos+seqgrep[2])
        mutantseqgrep<-substr(returnmutseq,mutpos-seqgrep[1],mutpos+seqgrep[2])
        nosplicegap_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1before,exon1after))-seqgrep[1]+1,nchar(paste0(refseq5,exon1before,exon1after))+seqgrep[2]),",",
                                substr(returnseq,nchar(paste0(refseq5,exon1before,exon1after,intron))-seqgrep[1]+1,nchar(paste0(refseq5,exon1before,exon1after,intron))+seqgrep[2]) )
        nosplicegap_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1before,exon1after))-seqgrep[1]+1,nchar(paste0(altseq5,exon1before,exon1after))+seqgrep[2]),",",
                                      substr(returnmutseq,nchar(paste0(altseq5,exon1before,exon1after,intron))-seqgrep[1]+1,nchar(paste0(altseq5,exon1before,exon1after,intron))+seqgrep[2]) )
        normalsplicegrep_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1before,exon1after))-seqgrep[1]+1,nchar(paste0(refseq5,exon1before,exon1after))),
                                         substr(returnseq,nchar(paste0(refseq5,exon1before,exon1after,intron))+1,nchar(paste0(refseq5,exon1before,exon1after,intron))+seqgrep[2]))
        normalsplicegrep_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1before,exon1after))-seqgrep[1]+1,nchar(paste0(altseq5,exon1before,exon1after))),
                                         substr(returnmutseq,nchar(paste0(altseq5,exon1before,exon1after,intron))+1,nchar(paste0(altseq5,exon1before,exon1after,intron))+seqgrep[2]))
        aberrantsplicegrep_normalseq<-paste0(substr(returnseq,ssgap-seqgrep[1]+1,ssgap),
                                           substr(returnseq,nchar(paste0(refseq5,exon1before,exon1after,intron))+1,nchar(paste0(refseq5,exon1before,exon1after,intron))+seqgrep[2]))
        aberrantsplicegrep_mutantseq<-paste0(substr(returnmutseq,ssgap-seqgrep[1]+1,ssgap),
                                           substr(returnmutseq,nchar(paste0(altseq5,exon1before,exon1after,intron))+1,nchar(paste0(altseq5,exon1before,exon1after,intron))+seqgrep[2]))
        returnframe<-data.frame(Var_Name=mut$name,Kind="Creating_Aberrant5SS",First_Exon=exon,Mutation_location=mutpos,Normal_SpliceGap=PredictedNormalSplicedOutPos,Aberrant_SpliceGap=PredictedAberantSplicedOutPos,Summary_Seq=summaryseq,
                              Ref_Construct=returnseq,Ref_Barcode=refcurrentbarcode,Ref_ID=rid,Ref_Content=content,Alt_Construct=returnmutseq,Alt_Barcode=altcurrentbarcode,Alt_ID=aid,Alt_Content=content,
                              Ref_MaxEntScoreSS5=result5$Ref_MaxEntScoreSS5,Ref_SS5score=result5$Ref_SS5score,Alt_MaxEntScoreSS5=result5$Alt_MaxEntScoreSS5,Alt_SS5score=result5$Alt_SS5score,Diff_MaxEntScoreSS5=result5$Diff_MaxEntScoreSS5,Diff_SS5score=result5$Diff_SS5score,
                              grep_NormalSeq=normalseqgrep,grep_MutantSeq=mutantseqgrep,grep_NoSplice_In_NormalSeq=nosplicegap_normalseq,NoSplice_In_MutantSeq=nosplicegap_mutantseq,
                              grep_NormalSplice_In_NormalSeq=normalsplicegrep_normalseq,grep_NormalSplice_IN_MutantSeq=normalsplicegrep_mutantseq,grep_AberrantSplice_In_NormalSeq=aberrantsplicegrep_normalseq,grep_AberrantSplice_IN_MutantSeq=aberrantsplicegrep_mutantseq,
                              Ref_Norm_Splice=PredictedNormalSplice,Ref_Norm_Splice_Length=PredictedNormalLength,Ref_Aberrant_Splice=PredictedAberrantSplice,Ref_Aberrant_Splice_Length=PredictedAberrantLength,
                            Mut_Norm_Splice=PredictedNormalSpliceMut,Mut_Norm_Splice_Length=PredictedNormalLengthMut,Mut_Aberrant_Splice=PredictedAberrantSpliceMut,Mut_Aberrant_Splice_Length=PredictedAberrantLengthMut)
        cat("Creating_Aberrant5SS",content,"\n")
      }
      if (!is.na(result5csds$pos)) {
        #mut on exon canonical splice donor site weighted exon1beforecsds,exon1aftercsds(+mut),introncsds,exon2ndcsds
        PredictedNormalSplice<-paste(tolower(refseq5),toupper(exon1beforecsds),toupper(exon1aftercsds),"[",tolower(introncsds),"]",toupper(exon2ndcsds),tolower(refseq3),sep="")
        exon1aftermutcsds<-paste(mut$alt,substr(exon1aftercsds,2,nchar(exon1aftercsds)),sep="")
        PredictedNormalSpliceMut<-paste(tolower(altseq5),toupper(exon1beforecsds),toupper(exon1aftermutcsds),"[",tolower(introncsds),"]",toupper(exon2ndcsds),tolower(altseq3),sep="")
        PredictedNormalLength<-nchar(returnseqcsds)-nchar(introncsds)
        PredictedNormalLengthMut<-nchar(returnmutseqcsds)-nchar(introncsds)
        #temp<-paste(exon2nd,altseq3,sep="")
        #In broken canonical 5SS, normal splice would not occur. So just give seq without splicing
        PredictedAberrantSplice<-returnseqcsds
        PredictedAberrantSpliceMut<-returnmutseqcsds
        PredictedAberrantLength<-nchar(returnseqcsds)
        PredictedAberrantLengthMut<-nchar(returnmutseqcsds)
        #(data.frame(pos=refresult$phypos[pos],Ref_MaxEntScoreSS5=refresult$MaxEntScoreSS5[pos],Ref_SS5score=refresult$SS5score[pos],Alt_MaxEntScoreSS5=altresult$MaxEntScoreSS5[pos],Alt_SS5score=altresult$SS5score[pos],Diff_MaxEntScoreSS5=(altresult$MaxEntScoreSS5[pos]-refresult$MaxEntScoreSS5[pos]),Diff_SS5score=(altresult$SS5score[pos]-refresult$SS5score[pos]) ))
        rid<-substr(paste(exon2ndcsds,refseq3,sep=""),nchar(exon2ndcsds)+1-nid[1],nchar(exon2ndcsds)+nchar(barcodes[1])+nid[2])
        aid<-substr(paste(exon2ndcsds,altseq3,sep=""),nchar(exon2ndcsds)+1-nid[1],nchar(exon2ndcsds)+nchar(barcodes[1])+nid[2])
        PredictedNormalSplicedOutPos<-paste0(nchar(paste0(refseq5,exon1beforecsds,exon1aftercsds))+1,"_",nchar(paste0(refseq5,exon1beforecsds,exon1aftercsds,introncsds)))
        PredictedAberantSplicedOutPos<-"No_Expected_Aberrant_Splice"
      
        content<-paste("Attached5Seq:",nchar(refseq5),";1stExon_BeforeTheMut:",aclencsds[1],";1stExon_On_AfterTheMut_Before_Canonical5SS:",aclencsds[2],";Intron_5:",aclencsds[4],";Intron_3:",aclencsds[5],";2ndExon:",aclencsds[6],";Attached3Seq:",nchar(refseq3),";Remain:",limit-sum(aclencsds),sep="")
        temp<-paste(tolower(IN$seq5),toupper(exon1beforecsds),toupper(exon1aftercsds),"[",tolower(introncsds),"]",toupper(exon2ndcsds),tolower(IN$seq3),sep="")
        summaryseq<-paste0(substr(temp,1,mutposcsds-1),"(",mut$ref,"/",mut$alt,")",substr(temp,mutposcsds+1,nchar(temp)))  #give (mut/alt) info
        normalseqgrep<-substr(returnseqcsds,mutposcsds-seqgrep[1],mutposcsds+seqgrep[2])
        mutantseqgrep<-substr(returnmutseqcsds,mutposcsds-seqgrep[1],mutposcsds+seqgrep[2])
        nosplicegap_normalseq<-paste0(substr(returnseqcsds,nchar(paste0(refseq5,exon1beforecsds,exon1aftercsds))-seqgrep[1]+1,nchar(paste0(refseq5,exon1beforecsds,exon1aftercsds))+seqgrep[2]),",",
                              substr(returnseqcsds,nchar(paste0(refseq5,exon1beforecsds,exon1aftercsds,introncsds))-seqgrep[1]+1,nchar(paste0(refseq5,exon1beforecsds,exon1aftercsds,introncsds))+seqgrep[2]) )
        nosplicegap_mutantseq<-paste0(substr(returnmutseqcsds,nchar(paste0(altseq5,exon1beforecsds,exon1aftercsds))-seqgrep[1]+1,nchar(paste0(altseq5,exon1beforecsds,exon1aftercsds))+seqgrep[2]),",",
                              substr(returnmutseqcsds,nchar(paste0(altseq5,exon1beforecsds,exon1aftercsds,introncsds))-seqgrep[1]+1,nchar(paste0(altseq5,exon1beforecsds,exon1aftercsds,introncsds))+seqgrep[2]) )
        normalsplicegrep_normalseq<-paste0(substr(returnseqcsds,nchar(paste0(refseq5,exon1beforecsds,exon1aftercsds))-seqgrep[1]+1,nchar(paste0(refseq5,exon1beforecsds,exon1aftercsds))),
                                   substr(returnseqcsds,nchar(paste0(refseq5,exon1beforecsds,exon1aftercsds,introncsds))+1,nchar(paste0(refseq5,exon1beforecsds,exon1aftercsds,introncsds))+seqgrep[2]))
        normalsplicegrep_mutantseq<-paste0(substr(returnmutseqcsds,nchar(paste0(altseq5,exon1beforecsds,exon1aftercsds))-seqgrep[1]+1,nchar(paste0(altseq5,exon1beforecsds,exon1aftercsds))),
                                   substr(returnmutseqcsds,nchar(paste0(altseq5,exon1beforecsds,exon1aftercsds,introncsds))+1,nchar(paste0(altseq5,exon1beforecsds,exon1aftercsds,introncsds))+seqgrep[2]))
        aberrantsplicegrep_normalseq<-"No_Expected_Aberrant_Splice"
        aberrantsplicegrep_mutantseq<-"No_Expected_Aberrant_Splice"
        newframe<-data.frame(Var_Name=mut$name,Kind="Broken_Canonical5SS",First_Exon=exon,Mutation_location=mutposcsds,Normal_SpliceGap=PredictedNormalSplicedOutPos,Aberrant_SpliceGap=PredictedAberantSplicedOutPos,Summary_Seq=summaryseq,
                             Ref_Construct=returnseqcsds,Ref_Barcode=refcurrentbarcode,Ref_ID=rid,Ref_Content=content,Alt_Construct=returnmutseqcsds,Alt_Barcode=altcurrentbarcode,Alt_ID=aid,Alt_Content=content,
                             Ref_MaxEntScoreSS5=result5csds$Ref_MaxEntScoreSS5,Ref_SS5score=result5csds$Ref_SS5score,Alt_MaxEntScoreSS5=result5csds$Alt_MaxEntScoreSS5,Alt_SS5score=result5csds$Alt_SS5score,Diff_MaxEntScoreSS5=result5csds$Diff_MaxEntScoreSS5,Diff_SS5score=result5csds$Diff_SS5score,
                             grep_NormalSeq=normalseqgrep,grep_MutantSeq=mutantseqgrep,grep_NoSplice_In_NormalSeq=nosplicegap_normalseq,NoSplice_In_MutantSeq=nosplicegap_mutantseq,
                             grep_NormalSplice_In_NormalSeq=normalsplicegrep_normalseq,grep_NormalSplice_IN_MutantSeq=normalsplicegrep_mutantseq,grep_AberrantSplice_In_NormalSeq=aberrantsplicegrep_normalseq,grep_AberrantSplice_IN_MutantSeq=aberrantsplicegrep_mutantseq,
                             Ref_Norm_Splice=PredictedNormalSplice,Ref_Norm_Splice_Length=PredictedNormalLength,Ref_Aberrant_Splice=PredictedAberrantSplice,Ref_Aberrant_Splice_Length=PredictedAberrantLength,
                                Mut_Norm_Splice=PredictedNormalSpliceMut,Mut_Norm_Splice_Length=PredictedNormalLengthMut,Mut_Aberrant_Splice=PredictedAberrantSpliceMut,Mut_Aberrant_Splice_Length=PredictedAberrantLengthMut)
        returnframe<-rbind(returnframe,newframe)
        cat("Broken_Canonical5SS",content,"\n")
      }
      if (nrow(returnframe)==0) return  (data.frame(Var_Name=NA ) ) else return(returnframe)
    } else {  #in case that mutation on an intron  **should consider both new 5SS and broken 5SS
      #mut on intron, exon1,intronbeforemut(+mut), intronaftermut,exon2nd
      #cat("here2\n")
      result5<-takeresult5(refresult5temp,altresult5temp,0)  #check for new 5SS
      ssgap<-as.numeric(result5$pos)+2
      result5csds<-takeresult5(refresult5temp,altresult5temp,1) #check for broken 5SS
      returnframe<-data.frame()
      PredictedNormalSplice<-paste(tolower(refseq5),toupper(exon1),"[",tolower(intronbeforemut),tolower(intronaftermut),"]",toupper(exon2nd),tolower(refseq3),sep="")
      intronbeforemutWithMut<-paste(substr(intronbeforemut,1,nchar(intronbeforemut)-1),mut$alt,sep="")
      PredictedNormalSpliceMut<-paste(tolower(altseq5),toupper(exon1),"[",tolower(intronbeforemutWithMut),tolower(intronaftermut),"]",toupper(exon2nd),tolower(altseq3),sep="")
      PredictedNormalLength<-nchar(returnseq)-nchar(intronbeforemut)-nchar(intronaftermut)
      PredictedNormalLengthMut<-nchar(returnmutseq)-nchar(intronbeforemutWithMut)-nchar(intronaftermut)
      rid<-substr(paste(exon2nd,refseq3,sep=""),nchar(exon2nd)+1-nid[1],nchar(exon2nd)+nchar(barcodes[1])+nid[2])
      aid<-substr(paste(exon2nd,altseq3,sep=""),nchar(exon2nd)+1-nid[1],nchar(exon2nd)+nchar(barcodes[1])+nid[2])
      PredictedNormalSplicedOutPos<-paste0(nchar(paste0(refseq5,exon1))+1,"_",nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut)))
      content<-paste("Attached5Seq:",nchar(refseq5),";1stExon:",aclen[1],";Intron_On_Before_Mutation:",aclen[3],";Intron_After_Mutation:",aclen[4],";Intron_3:",aclen[5],";2ndExon:",aclen[6],";Attached3Seq:",nchar(refseq3),";Remain:",limit-sum(aclen),sep="")
      if (!is.na(result5$pos) & ssgap<=nchar(paste0(tolower(refseq5),toupper(exon1))) ) {  #new gap should be located before canonical Splice Donor Site
        #temp<-paste(exon2nd,altseq3,sep="") 
        PredictedAberrantSplice<-paste(substr(returnseq,1,ssgap),"[",substr(returnseq,ssgap+1,nchar(paste(refseq5,exon1,intronbeforemut,intronaftermut,sep=""))),"]",substr(returnseq,nchar(paste(refseq5,exon1,intronbeforemut,intronaftermut,sep=""))+1,nchar(returnseq)),sep="")
        PredictedAberrantSpliceMut<-paste(substr(returnmutseq,1,ssgap),"[",substr(returnmutseq,ssgap+1,nchar(paste(altseq5,exon1,intronbeforemutWithMut,intronaftermut,sep=""))),"]",substr(returnmutseq,nchar(paste(altseq5,exon1,intronbeforemutWithMut,intronaftermut,sep=""))+1,nchar(returnmutseq)),sep="")
        PredictedAberrantLength<-PredictedNormalLength-(nchar(paste(refseq5,exon1))-ssgap)+2
        PredictedAberrantLengthMut<-PredictedNormalLengthMut-(nchar(paste(altseq5,exon1))-ssgap)+2
        #(data.frame(pos=refresult$phypos[pos],Ref_MaxEntScoreSS5=refresult$MaxEntScoreSS5[pos],Ref_SS5score=refresult$SS5score[pos],Alt_MaxEntScoreSS5=altresult$MaxEntScoreSS5[pos],Alt_SS5score=altresult$SS5score[pos],Diff_MaxEntScoreSS5=(altresult$MaxEntScoreSS5[pos]-refresult$MaxEntScoreSS5[pos]),Diff_SS5score=(altresult$SS5score[pos]-refresult$SS5score[pos]) ))
        PredictedAberantSplicedOutPos<-paste0(ssgap+1,"_",nchar(paste0(altseq5,exon1,intronbeforemutWithMut,intronaftermut)))
        temp<-paste0(substr(returnseq,1,ssgap),"{",substr(returnseq,ssgap+1,nchar(paste0(refseq5,exon1))),"[",tolower(intronbeforemut),tolower(intronaftermut),"]}",toupper(exon2nd),tolower(IN$seq3))
        temp2<-paste0(substr(temp,1,mutpos-1+2),"(",mut$ref,"/",mut$alt,")",substr(temp,mutpos+1+2,nchar(temp)))  #give (mut/alt) info *+2 for inserted { & [
        summaryseq<-paste0(tolower(IN$seq5),substr(temp2,nchar(IN$seq5)+1,nchar(temp2)))  #reset 5seq
        normalseqgrep<-substr(returnseq,mutpos-seqgrep[1],mutpos+seqgrep[2])
        mutantseqgrep<-substr(returnmutseq,mutpos-seqgrep[1],mutpos+seqgrep[2])
        nosplicegap_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1))-seqgrep[1]+1,nchar(paste0(refseq5,exon1))+seqgrep[2]),",",
                              substr(returnseq,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))-seqgrep[1]+1,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]) )
        nosplicegap_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1))-seqgrep[1]+1,nchar(paste0(altseq5,exon1))+seqgrep[2]),",",
                              substr(returnmutseq,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))-seqgrep[1]+1,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]) )
        normalsplicegrep_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1))-seqgrep[1]+1,nchar(paste0(refseq5,exon1))),
                                   substr(returnseq,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+1,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]))
        normalsplicegrep_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1))-seqgrep[1]+1,nchar(paste0(altseq5,exon1))),
                                   substr(returnmutseq,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))+1,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]))
        aberrantsplicegrep_normalseq<-paste0(substr(returnseq,ssgap-seqgrep[1]+1,ssgap),
                                     substr(returnseq,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+1,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]))
        aberrantsplicegrep_mutantseq<-paste0(substr(returnmutseq,ssgap-seqgrep[1]+1,ssgap),
                                     substr(returnmutseq,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))+1,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]))
        returnframe<-data.frame(Var_Name=mut$name,Kind="Creating_Aberrant5SS",First_Exon=exon,Mutation_location=mutpos,Normal_SpliceGap=PredictedNormalSplicedOutPos,Aberrant_SpliceGap=PredictedAberantSplicedOutPos,Summary_Seq=summaryseq,
                                Ref_Construct=returnseq,Ref_Barcode=refcurrentbarcode,Ref_ID=rid,Ref_Content=content,Alt_Construct=returnmutseq,Alt_Barcode=altcurrentbarcode,Alt_ID=aid,Alt_Content=content,
                                Ref_MaxEntScoreSS5=result5$Ref_MaxEntScoreSS5,Ref_SS5score=result5$Ref_SS5score,Alt_MaxEntScoreSS5=result5$Alt_MaxEntScoreSS5,Alt_SS5score=result5$Alt_SS5score,Diff_MaxEntScoreSS5=result5$Diff_MaxEntScoreSS5,Diff_SS5score=result5$Diff_SS5score,
                                grep_NormalSeq=normalseqgrep,grep_MutantSeq=mutantseqgrep,grep_NoSplice_In_NormalSeq=nosplicegap_normalseq,NoSplice_In_MutantSeq=nosplicegap_mutantseq,
                                grep_NormalSplice_In_NormalSeq=normalsplicegrep_normalseq,grep_NormalSplice_IN_MutantSeq=normalsplicegrep_mutantseq,grep_AberrantSplice_In_NormalSeq=aberrantsplicegrep_normalseq,grep_AberrantSplice_IN_MutantSeq=aberrantsplicegrep_mutantseq,
                                Ref_Norm_Splice=PredictedNormalSplice,Ref_Norm_Splice_Length=PredictedNormalLength,Ref_Aberrant_Splice=PredictedAberrantSplice,Ref_Aberrant_Splice_Length=PredictedAberrantLength,
                        Mut_Norm_Splice=PredictedNormalSpliceMut,Mut_Norm_Splice_Length=PredictedNormalLengthMut,Mut_Aberrant_Splice=PredictedAberrantSpliceMut,Mut_Aberrant_Splice_Length=PredictedAberrantLengthMut)
        cat("Creating_Aberrant5SS",content,"\n")
      }
      if (!is.na(result5csds$pos)) {
        #mut on exon canonical splice donor site weighted exon1beforecsds,exon1aftercsds(+mut),introncsds,exon2ndcsds
        #mut on intron, exon1,intronbeforemut(+mut), intronaftermut,exon2nd
        #temp<-paste(exon2nd,altseq3,sep="")
        #In broken canonical 5SS, normal splice would not occur. So just give seq without splicing
        PredictedAberrantSplice<-returnseq
        PredictedAberrantSpliceMut<-returnmutseq
        PredictedAberrantLength<-nchar(returnseq)
        PredictedAberrantLengthMut<-nchar(returnmutseq)
        #(data.frame(pos=refresult$phypos[pos],Ref_MaxEntScoreSS5=refresult$MaxEntScoreSS5[pos],Ref_SS5score=refresult$SS5score[pos],Alt_MaxEntScoreSS5=altresult$MaxEntScoreSS5[pos],Alt_SS5score=altresult$SS5score[pos],Diff_MaxEntScoreSS5=(altresult$MaxEntScoreSS5[pos]-refresult$MaxEntScoreSS5[pos]),Diff_SS5score=(altresult$SS5score[pos]-refresult$SS5score[pos]) ))
        PredictedAberantSplicedOutPos<-"No_Expected_Aberrant_Splice"
        temp<-paste(tolower(IN$seq5),toupper(exon1),"[",tolower(intronbeforemut),tolower(intronaftermut),"]",toupper(exon2nd),tolower(IN$seq3),sep="")
        summaryseq<-paste0(substr(temp,1,mutpos-1+1),"(",mut$ref,"/",mut$alt,")",substr(temp,mutpos+1+1,nchar(temp)))  #give (mut/alt) info *+1 for inserted [
        normalseqgrep<-substr(returnseq,mutpos-seqgrep[1],mutpos+seqgrep[2])
        mutantseqgrep<-substr(returnmutseq,mutpos-seqgrep[1],mutpos+seqgrep[2])
        nosplicegap_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1))-seqgrep[1]+1,nchar(paste0(refseq5,exon1))+seqgrep[2]),",",
                                      substr(returnseq,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))-seqgrep[1]+1,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]) )
        nosplicegap_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1))-seqgrep[1]+1,nchar(paste0(altseq5,exon1))+seqgrep[2]),",",
                                      substr(returnmutseq,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))-seqgrep[1]+1,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]) )
        normalsplicegrep_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1))-seqgrep[1]+1,nchar(paste0(refseq5,exon1))),
                                           substr(returnseq,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+1,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]))
        normalsplicegrep_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1))-seqgrep[1]+1,nchar(paste0(altseq5,exon1))),
                                           substr(returnmutseq,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))+1,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]))
        aberrantsplicegrep_normalseq<-"No_Expected_Aberrant_Splice"
        aberrantsplicegrep_mutantseq<-"No_Expected_Aberrant_Splice"
        newframe<-data.frame(Var_Name=mut$name,Kind="Broken_Canonical5SS",First_Exon=exon,Mutation_location=mutpos,Normal_SpliceGap=PredictedNormalSplicedOutPos,Aberrant_SpliceGap=PredictedAberantSplicedOutPos,Summary_Seq=summaryseq,
                             Ref_Construct=returnseq,Ref_Barcode=refcurrentbarcode,Ref_ID=rid,Ref_Content=content,Alt_Construct=returnmutseq,Alt_Barcode=altcurrentbarcode,Alt_ID=aid,Alt_Content=content,
                             Ref_MaxEntScoreSS5=result5csds$Ref_MaxEntScoreSS5,Ref_SS5score=result5csds$Ref_SS5score,Alt_MaxEntScoreSS5=result5csds$Alt_MaxEntScoreSS5,Alt_SS5score=result5csds$Alt_SS5score,Diff_MaxEntScoreSS5=result5csds$Diff_MaxEntScoreSS5,Diff_SS5score=result5csds$Diff_SS5score,
                             grep_NormalSeq=normalseqgrep,grep_MutantSeq=mutantseqgrep,grep_NoSplice_In_NormalSeq=nosplicegap_normalseq,NoSplice_In_MutantSeq=nosplicegap_mutantseq,
                             grep_NormalSplice_In_NormalSeq=normalsplicegrep_normalseq,grep_NormalSplice_IN_MutantSeq=normalsplicegrep_mutantseq,grep_AberrantSplice_In_NormalSeq=aberrantsplicegrep_normalseq,grep_AberrantSplice_IN_MutantSeq=aberrantsplicegrep_mutantseq,
                             Ref_Norm_Splice=PredictedNormalSplice,Ref_Norm_Splice_Length=PredictedNormalLength,Ref_Aberrant_Splice=PredictedAberrantSplice,Ref_Aberrant_Splice_Length=PredictedAberrantLength,
                     Mut_Norm_Splice=PredictedNormalSpliceMut,Mut_Norm_Splice_Length=PredictedNormalLengthMut,Mut_Aberrant_Splice=PredictedAberrantSpliceMut,Mut_Aberrant_Splice_Length=PredictedAberrantLengthMut)
        returnframe<-rbind(returnframe,newframe)
        cat("Broken_Canonical5SS",content,"\n")
      }
      if (nrow(returnframe)==0) returnframe<-data.frame(Var_Name=NA )
      return(returnframe)
    }  
  }  #ss5==1 end
  
  if (IN$ss3==1) {
    exon<-seq$exon[pos]
    if (exon==1) {
      cat("The given mutation is located on 1st exon. Therefore 3SS construct can NOT be designed because of 1st exon not being able to be assigned.\n")
      return (data.frame(Var_Name=NA )  )
    }
    if (exon>1) {
      exon<-exon-1  #in 3ss, the variant is located on 2nd exon, so # of 1st exon calculated
    }
    temppos<-pos
    #cat("pos:",pos,"exon:",seq$exon[pos],"\n")
    if (exon<0) {
      while (exon<0) {
        temppos<-temppos-1  #even in case of 3ss, find 1st exon like 5ss
        if (temppos==0) {
          cat("The given mutation is located before 1st exon. Therefore 3SS construct can NOT be designed because of 1st exon not being able to be assigned.\n")
          return (data.frame(Var_Name=NA )  )
        }
        exon<-seq$exon[temppos]
      }
    } 
    if (exon==max(seq$exon)) {
      cat("The given mutation is located on the last exon. Therefore 3SS construct can NOT be designed because of 2nd exon not being able to be assigned.\n")
      return (data.frame(Var_Name=NA )  )
    }
    
    #cat("Construct Length:",conlen,"\n")
    limit<-sum(conlen)
    aclen<<-c(0,0,0,0,0,0) # length for mutation on exon  #mut on exon new splice site weighted    exon1/intron1|intron2/exon2afterAcceptorSite/exon2beforeMut/Exon2AfterMut
    #conlencsds<-c(conlen[1],conlen[2],conlen[5],conlen[6],conlen[4],conlen[3]) #length for mutataion on intron #exon1/intron1/Intron2beforeMut/Intron2AfterMut/Intron2beforeAS/exon2
    #aclencsds<<-c(0,0,0,0,0,0)                 #exon1/intron1/Intron2beforeMut/Intron2AfterMut/Intron2beforeAS/exon2
    #get exon seq which the mutation resides in
    #if exon1st < conlen[1] _> the remaing length will move to intron2nd
    exon1<-takeseq(seq,pos,conlen[1],0,8)
    remain<-conlen[1]-nchar(exon1)
    conlen[1]<-aclen[1]
    #aclencsds[1]<<-aclen[1]
    #exon1csds<-exon1
    #conlencsds[1]<-aclencsds[1]  
    if (seq$exon[pos]>0) { #when mut is on an exon
      if (remain>0) {
        quarter<-as.integer(remain/4)
        quarterfinal<-remain-quarter*3
        #cat("#bp of Remain:",remain,"\t Quarter:",quarter,"\tQuarterFinal:",quarterfinal,"\n")
        conlen[3]<-conlen[3]+quarter;      conlen[4]<-conlen[4]+quarter;      conlen[5]<-conlen[5]+quarter;    conlen[6]<-conlen[6]+quarterfinal
        ###conlencsds[3]<-conlencsds[3]+quarter;      conlencsds[4]<-conlencsds[4]+quarter;      conlencsds[5]<-conlencsds[5]+quarter;    conlencsds[6]<-conlencsds[6]+quarterfinal
      }
       
      # canonical splice site weighted part start
      #aclenstore<-aclen #store aclen because csds make use of the same subroutine of aclen and change the parameter
      #exon1aftercsds<-takeseq(seq,pos,conlencsds[2],0,3) #include mutation #exon after mut and exon before canonical splice site. assume that mut is located near canonical splice site.
      #aclencsds[2]<<-aclen[2];       remaincsds<-conlencsds[2]-aclencsds[2]
      #if (remaincsds>0) {
      #  double<-as.integer(remaincsds/2)
      #  doublefinal<-remaincsds-double
      #  conlencsds[1]<-conlencsds[1]+double; conlencsds[4]<-conlencsds[4]+doublefinal
      #}
      ##########really!?
      #temp<-limitIntron(c(conlencsds[4],conlencsds[5]));conlencsds[4]<-temp[1];conlencsds[5]<-temp[2]
      #introncsds<-takeseq(seq,pos,conlencsds[4],conlencsds[5],2)
      #aclencsds[4]<<-aclen[4];aclencsds[5]<<-aclen[5]
      #remaincsds<-conlencsds[4]-aclencsds[4]+conlencsds[5]-aclencsds[5]
      #conlencsds[4]<-aclencsds[4];conlencsds[5]<-aclencsds[5]
      #if (remaincsds>0) {
      #  conlencsds[1]<-conlencsds[1]+remaincsds
      #}
      #exon1beforecsds<-takeseq(seq,pos,conlencsds[1],0,4)
      #aclencsds[1]<<-aclen[1];conlencsds[1]<-aclencsds[1]
      #remaincsds<-limit-nchar(exon1beforecsds)-nchar(exon1aftercsds)-nchar(introncsds)-nchar(exon2ndcsds)
      #if (remaincsds>0) {
      #  conlencsds[4]<-conlencsds[4]+remaincsds
      #  if (IN$limitintron!=-9) if ((conlencsds[4]+conlencsds[5])>IN$limitintron) {conlencsds[4]<-IN$limitintron-conlencsds[5]}
      #  introncsds<-takeseq(seq,pos,conlencsds[4],conlencsds[5],2)
      #  aclencsds[4]<<-aclen[4];aclencsds[5]<<-aclen[5];conlencsds[4]<-aclencsds[4];conlencsds[5]<-aclencsds[5]
      #  remaincsds<-limit-nchar(exon1beforecsds)-nchar(exon1aftercsds)-nchar(introncsds)-nchar(exon2ndcsds)
      #  if (remaincsds>0) {
      #    double<-as.integer(remaincsds/2)
      #    doublefinal<-remaincsds-double
      #    conlencsds[5]<-conlencsds[5]+double; conlencsds[6]<-conlencsds[6]+doublefinal
      #    if (IN$limitintron!=-9) if ((conlencsds[4]+conlencsds[5])>IN$limitintron) {conlencsds[5]<-IN$limitintron-conlencsds[4]}
      #    introncsds<-takeseq(seq,pos,conlencsds[4],conlencsds[5],2)
      #    exon2ndcsds<-takeseq(seq,pos,conlencsds[6],0,1)
      #    aclencsds[4]<<-aclen[4];aclencsds[5]<<-aclen[5];aclencsds[6]<<-aclen[6];conlencsds[4]<-aclencsds[4];conlencsds[5]<-aclencsds[5];conlencsds[6]<-aclencsds[6]
      #    remaincsds<-limit-nchar(exon1beforecsds)-nchar(exon1aftercsds)-nchar(introncsds)-nchar(exon2ndcsds)
      #    if (remaincsds>0) {
      #      conlencsds[6]<-conlencsds[6]+remaincsds
      #      exon2ndcsds<-takeseq(seq,pos,conlencsds[6],0,1)
      #      aclencsds[6]<<-aclen[6];conlencsds[6]<-aclencsds[6]
      #    }
      #  }
      #}
      #aclen<<-aclenstore  # canonical splice site weighted part end 
      
      temp<-limitIntron(c(conlen[2],conlen[3]));conlen[2]<-temp[1];conlen[3]<-temp[2]
      intron<-takeseq(seq,pos,conlen[2],conlen[3],2)  #intron 5 part and 3 part. if actual intron is shorter than indicated , merge and output
      remain<-conlen[2]+conlen[3]-nchar(intron)
      conlen[2]<-aclen[2];conlen[3]<-aclen[3]
      if (remain>0) {
        quarter<-as.integer(remain/3)
        quarterfinal<-remain-quarter*2
        #cat("#bp of Remain:",remain,"\t Quarter:",quarter,"\tQuarterFinal:",quarterfinal,"\n")
        conlen[4]<-conlen[4]+quarter;      conlen[5]<-conlen[5]+quarter
        conlen[6]<-conlen[6]+quarterfinal
      }
      #if exon beforemut within conlen[4] limit, the remaining length will move to exon after.
      exon2before<-takeseq(seq,pos,conlen[4],conlen[5],9)  #include mutation #exon before mut and exon after canonical splice acceptor site. if actual exon before mut and exon after cano splice acceptor site is shorter than indicated, then merged and output
      
      conlen[4]<-aclen[4];conlen[5]<-aclen[5]
      conlen[6]<-limit-nchar(exon2before)-nchar(intron)-nchar(exon1)
      exon2after<-takeseq(seq,pos,conlen[6],0,10)
      conlen[6]<-aclen[6]
  
      remain<-limit-nchar(exon2before)-nchar(exon2after)-nchar(intron)-nchar(exon1)
      if (remain>0) {   #remainings redistribution
        quarter<-as.integer(remain/3)
        quarterfinal<-remain-quarter*2
        #cat("#bp of Remain:",remain,"\t Quarter:",quarter,"\tQuarterFinal:",quarterfinal,"\n")
        conlen[4]<-conlen[4]+quarter;      conlen[5]<-conlen[5]+quarter;
        conlen[3]<-conlen[3]+quarterfinal
        exon2before<-takeseq(seq,pos,conlen[4],conlen[5],9) 
        conlen[4]<-aclen[4];conlen[5]<-aclen[5]
        if (IN$limitintron!=-9) if ((conlen[2]+conlen[3])>IN$limitintron) {conlen[3]<-IN$limitintron-conlen[2]}
        intron<-takeseq(seq,pos,conlen[2],conlen[3],2)
        conlen[2]<-aclen[2];conlen[3]<-aclen[3]
        remain<-limit-nchar(exon2before)-nchar(exon2after)-nchar(intron)-nchar(exon1)
        
        if (remain>0) {
          quarter<-as.integer(remain/2)
          quarterfinal<-remain-quarter
          #cat("#bp of Remain:",remain,"\t Quarter:",quarter,"\tQuarterFinal:",quarterfinal,"\n")
          conlen[2]<-conlen[2]+quarter;
          conlen[1]<-conlen[1]+quarterfinal
          if (IN$limitintron!=-9) if ((conlen[2]+conlen[3])>IN$limitintron) {conlen[2]<-IN$limitintron-conlen[3]}
          intron<-takeseq(seq,pos,conlen[2],conlen[3],2)
          exon1<-takeseq(seq,pos,conlen[1],0,8)
          conlen[1]<-aclen[1];conlen[2]<-aclen[2];conlen[3]<-aclen[3];
          remain<-limit-nchar(exon2before)-nchar(exon2after)-nchar(intron)-nchar(exon1)
          if (remain>0) {
            conlen[1]<-conlen[1]+remain
            #cat("#bp of Remain:",remain," conlen6",conlen[6]," aclen6",aclen[6],"exon1before",exon1before,"exon1after",exon1after,"intron","\n")
            exon1<-takeseq(seq,pos,conlen[1],0,8)
            conlen[1]<-aclen[1]
          }
        }
      }
    } else { #when mut is on an intron
      #cat("I am here\n")
      conlen<-c(conlen[1],conlen[2],conlen[5],conlen[6],conlen[4],conlen[3]) #length for mutataion on intron 
      #exon1/intron1/Intron2beforeMut/Intron2AfterMut/Intron2beforeAS/exon2
      if (remain>0) {
        quarter<-as.integer(remain/4)
        quarterfinal<-remain-quarter*3
        #cat("#bp of Remain:",remain,"\t Quarter:",quarter,"\tQuarterFinal:",quarterfinal,"\n")
        conlen[2]<-conlen[2]+quarter;conlen[3]<-conlen[3]+quarter;conlen[4]<-conlen[4]+quarter;
        conlen[1]<-conlen[1]+quarterfinal
      }
      intronbeforemut<-takeseq(seq,pos,conlen[2],conlen[3],11)  #intron 5 part before mut and include mut
      remain<-conlen[2]+conlen[3]-aclen[2]-aclen[3]
      conlen[2]<-aclen[2];conlen[3]<-aclen[3]
 
      if (remain>0) {
        quarter<-as.integer(remain/3)
        quarterfinal<-remain-quarter*2
        #cat("#bp of Remain:",remain,"\t Quarter:",quarter,"\tQuarterFinal:",quarterfinal,"\n")
        conlen[4]<-conlen[4]+quarter;conlen[5]<-conlen[5]+quarter
        conlen[6]<-conlen[6]+quarterfinal
      }
      intronaftermut<-takeseq(seq,pos,conlen[4],conlen[5],6) #intron5 part after mut (NOT include mut) and intron3 part
      
      temp<-limitIntron(c(conlen[2],conlen[3],conlen[4] ,conlen[5]));conlen[2]<-temp[1];conlen[3]<-temp[2];conlen[4]<-temp[3];conlen[5]<-temp[4]
      intronbeforemut<-takeseq(seq,pos,conlen[2],conlen[3],11) 
      intronaftermut<-takeseq(seq,pos,conlen[4],conlen[5],6)
      conlen[2]<-aclen[2];conlen[3]<-aclen[3];conlen[4]<-aclen[4];conlen[5]<-aclen[5]
      remain<-limit-nchar(intronbeforemut)-nchar(intronaftermut)-nchar(exon1)-conlen[6]
      #cat("limit;",limit,"ib;",nchar(intronbeforemut),"iaf",nchar(intronaftermut),"exon1;",nchar(exon1),"remain:",remain,"\t")
      if (remain>0) {
        conlen[6]<-conlen[6]+remain
      }
      #cat("I m here!\n")
      exon2nd<-takeseq(seq,pos,conlen[6],0,1) 
      conlen[6]<-aclen[6]
      #cat("I m here too!\n")
      #cat("here1 conlen:",conlen,"\naclen:",aclen,"\n")
      #cat("I m here!\n")
      remain<-limit-nchar(exon1)-nchar(intronbeforemut)-nchar(intronaftermut)-nchar(exon2nd)
      if (remain>0) {   #remainings redistribution
        quarter<-as.integer(remain/3)
        quarterfinal<-remain-quarter*2
        conlen[3]<-conlen[3]+quarter
        conlen[4]<-conlen[4]+quarter
        conlen[5]<-conlen[5]+quarterfinal
        if (IN$limitintron!=-9) if ((conlen[2]+conlen[3]+conlen[4]+conlen[5]) > IN$limitintron) {conlen[4]<-IN$limitintron-conlen[2]-conlen[3]-conlen[5]}
        intronbeforemut<-takeseq(seq,pos,conlen[2],conlen[3],11) 
        intronaftermut<-takeseq(seq,pos,conlen[4],conlen[5],6) 
        conlen[2]<-aclen[2];conlen[3]<-aclen[3];conlen[4]<-aclen[4];conlen[5]<-aclen[5]
        #cat("here2 conlen:",conlen,"\naclen:",aclen,"\n")
        remain<-limit-nchar(exon1)-nchar(intronbeforemut)-nchar(intronaftermut)-nchar(exon2nd)
        if (remain>0) {
          quarter<-as.integer(remain/3)
          quarterfinal<-remain-quarter*2
          #cat("#bp of Remain:",remain,"\t Quarter:",quarter,"\tQuarterFinal:",quarterfinal,"\n")
          conlen[1]<-conlen[1]+quarter;conlen[3]<-conlen[3]+quarter;
          conlen[2]<-conlen[2]+quarterfinal
          if (IN$limitintron!=-9) if ((conlen[2]+conlen[3]+conlen[4]+conlen[5])>IN$limitintron) {
            while ((conlen[2]+conlen[3]+conlen[4]+conlen[5])>IN$limitintron) { 
              conlen[2]<-conlen[2]-1;conlen[3]<-conlen[3]-1;
            }
          }
          intronbeforemut<-takeseq(seq,pos,conlen[2],conlen[3],11) 
          exon1<-takeseq(seq,pos,conlen[1],0,7) 
          conlen[1]<-aclen[1];conlen[2]<-aclen[2];conlen[3]<-aclen[3]
          #cat("here3 conlen:",conlen,"\naclen:",aclen,"\n")
          remain<-limit-nchar(exon1)-nchar(intronbeforemut)-nchar(intronaftermut)-nchar(exon2nd)
          if (remain>0) {
            conlen[1]<-conlen[1]+remain
            #cat("#bp of Remain:",remain," conlen6",conlen[6]," aclen6",aclen[6],"exon1before",exon1before,"exon1after",exon1after,"intron","\n")
            exon1<-takeseq(seq,pos,conlen[1],0,7) 
            conlen[1]<-aclen[1]
            #cat("here4 conlen:",conlen,"\naclen:",aclen,"\n")
          }
        }
      }
      #cat("finale conlen:",conlen,"\naclen:",aclen,"\n")
    }
    #assemble entire construct
    #mut on exon new splice site weighted    exon1before,exon1after(+mut),intron,exon2nd
    #mut on exon canonical splice donor site weighted exon1beforecsds,exon1aftercsds(+mut),introncsds,exon2ndcsds
    #mut on intron, exon1,intronbeforemut(+mut), intronaftermut,exon2nd
    #if ((seq$exon[pos]>0) & ((nchar(exon2nd)*nchar(intron)*nchar(exon1after)*nchar(exon1before)==0) | (nchar(exon2ndcsds)*nchar(introncsds)*nchar(exon1aftercsds)*nchar(exon1beforecsds)==0) )) {
    #  cat("Zero-length component found.\nNew Splice Site Weighte Construct   Exon before the mutation:",nchar(exon1before),"Exon after the mutation:",nchar(exon1after),"Intron:",nchar(intron),"2nd Exon:",nchar(exon2nd),
    #  "\n Canonical Splice Donor Site Weighted Construct   Exon before the mutation:",nchar(exon1beforecsds),"Exon after the mutation:",nchar(exon1aftercsds),"Intron:",nchar(introncsds),"2nd Exon:",nchar(exon2ndcsds),"\n")
    #  return (data.frame(Var_Name=NA )  )  #if any of elements lacks, null returned
    #} else if ((seq$exon[pos]<0) & (nchar(exon2nd)*nchar(intronaftermut)*nchar(intronbeforemut)*nchar(exon1)==0)) {
    #  cat("Zero-length component found.\n1st Exon:",nchar(exon1),"Intron before the mutation:",nchar(intronbeforemut),"Intron after the mutation:",nchar(intronaftermut),"2nd Exon:",nchar(exon2nd),"\n")
    #  return (data.frame(Var_Name=NA )  )  #if any of elements lacks, null returned
    #
    
    #barcode embed
    #barcodesYES<-0  #if employ barcode or not
    #barcodes<-c()  #barcodes to be employed
    #barcodesPOS<-1   #pointer in the barcodes
    #temp<-regexpr("\\d+",mut.change[i,]$Var_5SplicePoint) #instead of variant position, actual 5ss gap taken
    #x<-as.numeric(substr(mut.change[i,]$Var_5SplicePoint,as.numeric(temp),as.numeric(temp)+attr(temp,"match.length")-1))
    #

    refseq5<-IN$seq5;    refseq3<-IN$seq3
    refcurrentbarcode<-NA
    altseq5<-IN$seq5;    altseq3<-IN$seq3
    altcurrentbarcode<-NA
    
    if (barcodesYES) {
      barcodesPOSadd<-0;
      barcodesPOS2<-barcodesPOS+1
      if (barcodesPOS2>length(barcodes)) barcodesPOS2<-1
      temp3<-regexpr("(X)+",IN$seq3)
      temp3s<-as.numeric(temp3);temp3e<-temp3s+attr(temp3,"match.length")-1
      if (temp3s>0) {
        refseq3<-paste(substr(IN$seq3,1,temp3s-1),barcodes[barcodesPOS],substr(IN$seq3,temp3e+1,nchar(IN$seq3)),sep="")
        altseq3<-paste(substr(IN$seq3,1,temp3s-1),barcodes[barcodesPOS2],substr(IN$seq3,temp3e+1,nchar(IN$seq3)),sep="")
        barcodesPOSadd<-2
        refcurrentbarcode<-barcodes[barcodesPOS]
        altcurrentbarcode<-barcodes[barcodesPOS2]
      }
      
      temp5<-regexpr("(X)+",IN$seq5)
      temp5s<-as.numeric(temp5);temp5e<-temp5s+attr(temp5,"match.length")-1
      if (temp5s>0) {
        refseq5<-paste(substr(IN$seq5,1,temp5s-1),barcodes[barcodesPOS],substr(IN$seq5,temp5e+1,nchar(IN$seq5)),sep="")
        altseq5<-paste(substr(IN$seq5,1,temp5s-1),barcodes[barcodesPOS2],substr(IN$seq5,temp5e+1,nchar(IN$seq5)),sep="")
        barcodesPOSadd<-2
        refcurrentbarcode<-barcodes[barcodesPOS]
        altcurrentbarcode<-barcodes[barcodesPOS2]
      }
      barcodesPOS<<-barcodesPOS+barcodesPOSadd
      if (barcodesPOS>length(barcodes)) {
        barcodesPOS<<-barcodesPOS-length(barcodes)
        cat("barcode cycle returned to the initial\n")
      }
    }
    
    if (seq$exon[pos]>0) {
      #cat("here1\n")
      returnseq<-paste(tolower(refseq5),toupper(exon1),tolower(intron),toupper(exon2before),toupper(exon2after),tolower(refseq3),sep="")
      returnmutseq<-paste(tolower(altseq5),toupper(exon1),tolower(intron),toupper(exon2before),toupper(exon2after),tolower(altseq3),sep="")
      mutpos<-nchar(refseq5)+nchar(exon1)+nchar(intron)+nchar(exon2before)
      
      #returnseqcsds<-paste(tolower(refseq5),toupper(exon1beforecsds),toupper(exon1aftercsds),tolower(introncsds),toupper(exon2ndcsds),tolower(refseq3),sep="")
      #returnmutseqcsds<-paste(tolower(altseq5),toupper(exon1beforecsds),toupper(exon1aftercsds),tolower(introncsds),toupper(exon2ndcsds),tolower(altseq3),sep="")
      #mutposcsds<-nchar(refseq5)+nchar(exon1beforecsds)+1
    } else {
      #cat("here2\n")
      returnseq<-paste(tolower(refseq5),toupper(exon1),tolower(intronbeforemut),tolower(intronaftermut),toupper(exon2nd),tolower(refseq3),sep="")
      returnmutseq<-paste(tolower(altseq5),toupper(exon1),tolower(intronbeforemut),tolower(intronaftermut),toupper(exon2nd),tolower(altseq3),sep="")
      mutpos<-nchar(refseq5)+nchar(exon1)+nchar(intronbeforemut)
    }
    #check score in ref seq
    cat("Checking the potential of 3\'splice site  in the referrence construct...\n")
    Start<-mutpos-23+1
    if (Start<1) Start<-1  #indicate the region where the mutation affects
    End <-mutpos
    seqtemp<-makeseqframe(returnseq)
    #cat("start;",Start,"end;",End, "\n")
    #write.table(seq,file="temp.seq.txt",sep="\t",quote=F,col.names=T,row.names=F)
    refresult3temp<-SS3scan(seqtemp,Start,End)

    #make mut allele  ###in this ver. just SNV permitted
    if (toupper(mut$ref)==toupper(substr(returnmutseq,mutpos,mutpos))) {
      if (seq$exon[pos]>0) {
        returnmutseq<-paste(substr(returnmutseq,1,mutpos-1),toupper(mut$alt),substr(returnmutseq,mutpos+1,nchar(returnmutseq)),sep="")
      } else {
        returnmutseq<-paste(substr(returnmutseq,1,mutpos-1),tolower(mut$alt),substr(returnmutseq,mutpos+1,nchar(returnmutseq)),sep="")
      }
      
    } else {
      cat("Ref allele not found in the actual sequence.\n Ref in the mutation file:",mut$ref,"Ref in the actual sequence:",substr(returnmutseq,mutpos,mutpos),"\n")
      cat ("Sequence in Reference Construct:",substr(returnseq,1,mutpos-1),"[",substr(returnseq,mutpos,mutpos),"]",substr(returnseq,mutpos+1,nchar(returnseq)),"\n",sep="")
      cat ("Sequence in Mutant Construct before MutGenesis:",substr(returnmutseq,1,mutpos-1),"[",substr(returnmutseq,mutpos,mutpos),"]",substr(returnmutseq,mutpos+1,nchar(returnmutseq)),"\n",sep="")
      return (    data.frame(Var_Name=NA )   )
    }
    seqtemp<-makeseqframe(returnmutseq)
    altresult3temp<-SS3scan(seqtemp,Start,End)

    #write.table(result5,file="temp.result5.txt",sep="\t",quote=F,col.names=T,row.names=F)
    
    #cat("ssgap:",ssgap,"\n")
    #cat("pos:",pos,"exon:",seq$exon[pos],"\n")
    if (seq$exon[pos]>0 ) {  #in case that mutation on an exon
      #cat("here1\n")
      result3<-takeresult3(refresult3temp,altresult3temp,0)  #create 3ss
      result3csds<-takeresult3(refresult3temp,altresult3temp,1)  #broken 3ss
      ssgap<-as.numeric(result3$pos)+19

      #Start<-mutposcsds-9+1;if (Start<1) Start<-1 
      #End <-mutposcsds
      #seqtemp<-makeseqframe(returnseqcsds)
      #refresult5csdstemp<-SS5scan(seqtemp,Start,End)
      #returnmutseqcsds<-paste(substr(returnmutseqcsds,1,mutposcsds-1),mut$alt,substr(returnmutseqcsds,mutposcsds+1,nchar(returnmutseqcsds)),sep="")
      #seqtemp<-makeseqframe(returnmutseqcsds)
      #altresult5csdstemp<-SS5scan(seqtemp,Start,End)
      #result5csds<-takeresult5(refresult5csdstemp,altresult5csdstemp,1)
      returnframe<-data.frame()
      if (!is.na(result3$pos)) {      
        PredictedNormalSplice<-paste(tolower(refseq5),toupper(exon1),"[",tolower(intron),"]",toupper(exon2before),toupper(exon2after),tolower(refseq3),sep="")
        exon2beforemut<-paste(substr(exon2before,1,nchar(exon2before)-1),mut$alt,sep="")
        PredictedNormalSpliceMut<-paste(tolower(altseq5),toupper(exon1),"[",tolower(intron),"]",toupper(exon2beforemut),toupper(exon2after),tolower(altseq3),sep="")
        PredictedNormalLength<-nchar(returnseq)-nchar(intron)
        PredictedNormalLengthMut<-nchar(returnmutseq)-nchar(intron)
        #temp<-paste(exon2nd,altseq3,sep="")
        PredictedAberrantSplice<-paste0(substr(returnseq,1,nchar(paste0(refseq5,exon1))),"[",substr(returnseq,nchar(paste0(refseq5,exon1))+1,ssgap),"]",substr(returnseq,ssgap+1,nchar(returnseq)) )
        PredictedAberrantSpliceMut<-paste0(substr(returnmutseq,1,nchar(paste0(altseq5,exon1))),"[",substr(returnmutseq,nchar(paste0(altseq5,exon1))+1,ssgap),"]",substr(returnmutseq,ssgap+1,nchar(returnmutseq)) )
       
        PredictedAberrantLength<-nchar(returnseq)-ssgap+nchar(paste0(refseq5,exon1))  #PredictedNormalLength-(nchar(paste(refseq5,exon1before,exon1after))-ssgap)+2
        PredictedAberrantLengthMut<-nchar(returnmutseq)-ssgap+nchar(paste0(altseq5,exon1)) #PredictedNormalLengthMut-(nchar(paste(altseq5,exon1before,exon1aftermut))-ssgap)+2
        #(data.frame(pos=refresult$phypos[pos],Ref_MaxEntScoreSS5=refresult$MaxEntScoreSS5[pos],Ref_SS5score=refresult$SS5score[pos],Alt_MaxEntScoreSS5=altresult$MaxEntScoreSS5[pos],Alt_SS5score=altresult$SS5score[pos],Diff_MaxEntScoreSS5=(altresult$MaxEntScoreSS5[pos]-refresult$MaxEntScoreSS5[pos]),Diff_SS5score=(altresult$SS5score[pos]-refresult$SS5score[pos]) ))
        rid<-substr(paste(exon2after,refseq3,sep=""),nchar(exon2after)+1-nid[1],nchar(exon2after)+nchar(barcodes[1])+nid[2]) #ID
        aid<-substr(paste(exon2after,altseq3,sep=""),nchar(exon2after)+1-nid[1],nchar(exon2after)+nchar(barcodes[1])+nid[2]) #ID
        PredictedNormalSplicedOutPos<-paste0(nchar(paste0(refseq5,exon1))+1,"_",nchar(paste0(refseq5,exon1,intron)))
        PredictedAberantSplicedOutPos<-paste0(nchar(paste0(refseq5,exon1))+1,"_",ssgap)
        content<-paste("Attached5Seq:",nchar(refseq5),";1stExon:",aclen[1],";Intron_5:",aclen[2],";Intron_3:",aclen[3],";2ndExon_After_SpliceAcceptorSite:",aclen[4],";2ndExon_Before_Mutation(Incl Mut):",aclen[5],";2ndExon_After_Mutation:",aclen[6],";Attached3Seq:",nchar(refseq3),";Remain:",limit-sum(aclen),sep="")
        if (ssgap<nchar(paste0(refseq5,exon1,intron))) {  #in case that aberrant acceptor is located before canonical acceptor site   #######check the script in 5SS
          temp<-paste0(substr(returnseq,1,nchar(paste0(refseq5,exon1)) ),"{[",substr(returnseq,nchar(paste0(refseq5,exon1))+1,ssgap),"}",substr(returnseq,ssgap+1,nchar(paste0(refseq5,exon1,intron))),"]",substr(returnseq,nchar(paste0(refseq5,exon1,intron))+1,nchar(paste0(refseq5,exon1,intron,exon2before,exon2after))),tolower(IN$seq3))
        } else { #in case that aberrant acceptor is located after canonical acceptor site
          temp<-paste0(substr(returnseq,1,nchar(paste0(refseq5,exon1)) ),"{[",substr(returnseq,nchar(paste0(refseq5,exon1))+1,nchar(paste0(refseq5,exon1,intron))),"]",substr(returnseq,nchar(paste0(refseq5,exon1,intron))+1,ssgap),"}",substr(returnseq,ssgap+1,nchar(paste0(refseq5,exon1,intron,exon2before,exon2after))),tolower(IN$seq3))
        }
        add<-2 #shift for {[
        if (mutpos>ssgap) add<-add+1  #shift for }
        if (mutpos>nchar(paste0(refseq5,exon1,intron))) add<-add+1 #shift for ]
        temp2<-paste0(substr(temp,1,mutpos-1+add),"(",mut$ref,"/",mut$alt,")",substr(temp,mutpos+1+add,nchar(temp)))  #give (mut/alt) info *+add means a shift to rightside because of inserted {[ and ] or }
        summaryseq<-paste0(tolower(IN$seq5),substr(temp2,nchar(IN$seq5)+1,nchar(temp2)))  #reset 5seq
        normalseqgrep<-substr(returnseq,mutpos-seqgrep[1],mutpos+seqgrep[2])
        mutantseqgrep<-substr(returnmutseq,mutpos-seqgrep[1],mutpos+seqgrep[2])
        nosplicegap_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1))-seqgrep[1]+1,nchar(paste0(refseq5,exon1))+seqgrep[2]),",",
                                      substr(returnseq,nchar(paste0(refseq5,exon1,intron))-seqgrep[1]+1,nchar(paste0(refseq5,exon1,intron))+seqgrep[2]) )
        nosplicegap_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1))-seqgrep[1]+1,nchar(paste0(altseq5,exon1))+seqgrep[2]),",",
                                      substr(returnmutseq,nchar(paste0(altseq5,exon1,intron))-seqgrep[1]+1,nchar(paste0(altseq5,exon1,intron))+seqgrep[2]) )
        normalsplicegrep_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1))-seqgrep[1]+1,nchar(paste0(refseq5,exon1))),
                                           substr(returnseq,nchar(paste0(refseq5,exon1,intron))+1,nchar(paste0(refseq5,exon1,intron))+seqgrep[2]))
        normalsplicegrep_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1))-seqgrep[1]+1,nchar(paste0(altseq5,exon1))),
                                           substr(returnmutseq,nchar(paste0(altseq5,exon1,intron))+1,nchar(paste0(altseq5,exon1,intron))+seqgrep[2]))
        aberrantsplicegrep_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1))-seqgrep[1]+1,nchar(paste0(refseq5,exon1))),
                                             substr(returnseq,ssgap+1,ssgap+seqgrep[2]) )
        aberrantsplicegrep_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1))-seqgrep[1]+1,nchar(paste0(altseq5,exon1))),
                                             substr(returnmutseq,ssgap+1,ssgap+seqgrep[2]) )
        #cat("Var_Name",length(mut$name),"Kind","Creating_Aberrant3SS","First_Exon",length(exon),"Mutation_location",length(mutpos),"Normal_SpliceGap",length(PredictedNormalSplicedOutPos),"Aberrant_SpliceGap",length(PredictedAberantSplicedOutPos),"Summary_Seq",length(summaryseq),
        #    "Ref_Construct",length(returnseq),"Ref_Barcode",length(refcurrentbarcode),"Ref_ID",length(rid),"Ref_Content",length(content),"Alt_Construct",length(returnmutseq),"Alt_Barcode",length(altcurrentbarcode),"Alt_ID",length(aid),"Alt_Content",length(content),
        #    "Ref_MaxEntScoreSS3",length(result3$Ref_MaxEntScoreSS3),"Ref_SS3score",length(result3$Ref_SS3score),"Alt_MaxEntScoreSS3",length(result3$Alt_MaxEntScoreSS3),"Alt_SS3score",length(result3$Alt_SS5score),"Diff_MaxEntScoreSS3",length(result3$Diff_MaxEntScoreSS3),"Diff_SS3score",length(result3$Diff_SS3score),
        #    "grep_NormalSeq",length(normalseqgrep),"grep_MutantSeq",length(mutantseqgrep),"grep_NoSplice_In_NormalSeq",length(nosplicegap_normalseq),"NoSplice_In_MutantSeq",length(nosplicegap_mutantseq),
        #    "grep_NormalSplice_In_NormalSeq",length(normalsplicegrep_normalseq),"grep_NormalSplice_IN_MutantSeq",length(normalsplicegrep_mutantseq),"grep_AberrantSplice_In_NormalSeq",length(aberrantsplicegrep_normalseq),"grep_AberrantSplice_IN_MutantSeq",length(aberrantsplicegrep_mutantseq),
        #    "Ref_Norm_Splice",length(PredictedNormalSplice),"Ref_Norm_Splice_Length",length(PredictedNormalLength),"Ref_Aberrant_Splice",length(PredictedAberrantSplice),"Ref_Aberrant_Splice_Length",length(PredictedAberrantLength),
        #    "Mut_Norm_Splice",length(PredictedNormalSpliceMut),"Mut_Norm_Splice_Length",length(PredictedNormalLengthMut),"Mut_Aberrant_Splice",length(PredictedAberrantSpliceMut),"Mut_Aberrant_Splice_Length",length(PredictedAberrantLengthMut),"\n")
        
        returnframe<-data.frame(Var_Name=mut$name,Kind="Creating_Aberrant3SS",First_Exon=exon,Mutation_location=mutpos,Normal_SpliceGap=PredictedNormalSplicedOutPos,Aberrant_SpliceGap=PredictedAberantSplicedOutPos,Summary_Seq=summaryseq,
                                Ref_Construct=returnseq,Ref_Barcode=refcurrentbarcode,Ref_ID=rid,Ref_Content=content,Alt_Construct=returnmutseq,Alt_Barcode=altcurrentbarcode,Alt_ID=aid,Alt_Content=content,
                                Ref_MaxEntScoreSS3=result3$Ref_MaxEntScoreSS3,Ref_SS3score=result3$Ref_SS3score,Alt_MaxEntScoreSS3=result3$Alt_MaxEntScoreSS3,Alt_SS3score=result3$Alt_SS3score,Diff_MaxEntScoreSS3=result3$Diff_MaxEntScoreSS3,Diff_SS3score=result3$Diff_SS3score,
                                grep_NormalSeq=normalseqgrep,grep_MutantSeq=mutantseqgrep,grep_NoSplice_In_NormalSeq=nosplicegap_normalseq,NoSplice_In_MutantSeq=nosplicegap_mutantseq,
                                grep_NormalSplice_In_NormalSeq=normalsplicegrep_normalseq,grep_NormalSplice_IN_MutantSeq=normalsplicegrep_mutantseq,grep_AberrantSplice_In_NormalSeq=aberrantsplicegrep_normalseq,grep_AberrantSplice_IN_MutantSeq=aberrantsplicegrep_mutantseq,
                                Ref_Norm_Splice=PredictedNormalSplice,Ref_Norm_Splice_Length=PredictedNormalLength,Ref_Aberrant_Splice=PredictedAberrantSplice,Ref_Aberrant_Splice_Length=PredictedAberrantLength,
                                Mut_Norm_Splice=PredictedNormalSpliceMut,Mut_Norm_Splice_Length=PredictedNormalLengthMut,Mut_Aberrant_Splice=PredictedAberrantSpliceMut,Mut_Aberrant_Splice_Length=PredictedAberrantLengthMut)
        cat("Creating_Aberrant3SS",content,"\n")
      }
      if (!is.na(result3csds$pos)) {  #broken 3SS part
        PredictedNormalSplice<-paste(tolower(refseq5),toupper(exon1),"[",tolower(intron),"]",toupper(exon2before),toupper(exon2after),tolower(refseq3),sep="")
        exon2beforemut<-paste(substr(exon2before,1,nchar(exon2before)-1),mut$alt,sep="")
        PredictedNormalSpliceMut<-paste(tolower(altseq5),toupper(exon1),"[",tolower(intron),"]",toupper(exon2beforemut),toupper(exon2after),tolower(altseq3),sep="")
        PredictedNormalLength<-nchar(returnseq)-nchar(intron)
        PredictedNormalLengthMut<-nchar(returnmutseq)-nchar(intron)
        #temp<-paste(exon2nd,altseq3,sep="")
        #In broken canonical 5SS, normal splice would not occur. So just give seq without splicing
        PredictedAberrantSplice<-returnseq
        PredictedAberrantSpliceMut<-returnmutseq
        PredictedAberrantLength<-nchar(returnseq)
        PredictedAberrantLengthMut<-nchar(returnmutseq)
        #(data.frame(pos=refresult$phypos[pos],Ref_MaxEntScoreSS5=refresult$MaxEntScoreSS5[pos],Ref_SS5score=refresult$SS5score[pos],Alt_MaxEntScoreSS5=altresult$MaxEntScoreSS5[pos],Alt_SS5score=altresult$SS5score[pos],Diff_MaxEntScoreSS5=(altresult$MaxEntScoreSS5[pos]-refresult$MaxEntScoreSS5[pos]),Diff_SS5score=(altresult$SS5score[pos]-refresult$SS5score[pos]) ))
        rid<-substr(paste(exon2after,refseq3,sep=""),nchar(exon2after)+1-nid[1],nchar(exon2after)+nchar(barcodes[1])+nid[2]) #ID
        aid<-substr(paste(exon2after,altseq3,sep=""),nchar(exon2after)+1-nid[1],nchar(exon2after)+nchar(barcodes[1])+nid[2]) #ID
        PredictedNormalSplicedOutPos<-paste0(nchar(paste0(refseq5,exon1))+1,"_",nchar(paste0(refseq5,exon1,intron)))
        PredictedAberantSplicedOutPos<-"No_Expected_Aberrant_Splice"
        content<-paste("Attached5Seq:",nchar(refseq5),";1stExon:",aclen[1],";Intron_5:",aclen[2],";Intron_3:",aclen[3],";2ndExon_After_SpliceAcceptorSite:",aclen[4],";2ndExon_Before_Mutation(Incl Mut):",aclen[5],";2ndExon_After_Mutation:",aclen[6],";Attached3Seq:",nchar(refseq3),";Remain:",limit-sum(aclen),sep="")
        temp<-paste0(substr(returnseq,1,nchar(paste0(refseq5,exon1)) ),"[",substr(returnseq,nchar(paste0(refseq5,exon1))+1,nchar(paste0(refseq5,exon1,intron))),"]",substr(returnseq,nchar(paste0(refseq5,exon1,intron))+1,nchar(paste0(refseq5,exon1,intron,exon2before,exon2after))),tolower(IN$seq3))
        add<-1 #shift for [
        if (mutpos>nchar(paste0(refseq5,exon1,intron))) add<-add+1 #shift for ]
        temp2<-paste0(substr(temp,1,mutpos-1+add),"(",mut$ref,"/",mut$alt,")",substr(temp,mutpos+1+add,nchar(temp)))  #give (mut/alt) info *+add means a shift to rightside because of inserted {[ and ] or }
        summaryseq<-paste0(tolower(IN$seq5),substr(temp2,nchar(IN$seq5)+1,nchar(temp2)))  #reset 5seq
        normalseqgrep<-substr(returnseq,mutpos-seqgrep[1],mutpos+seqgrep[2])
        mutantseqgrep<-substr(returnmutseq,mutpos-seqgrep[1],mutpos+seqgrep[2])
        nosplicegap_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1))-seqgrep[1]+1,nchar(paste0(refseq5,exon1))+seqgrep[2]),",",
                                      substr(returnseq,nchar(paste0(refseq5,exon1,intron))-seqgrep[1]+1,nchar(paste0(refseq5,exon1,intron))+seqgrep[2]) )
        nosplicegap_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1))-seqgrep[1]+1,nchar(paste0(altseq5,exon1))+seqgrep[2]),",",
                                      substr(returnmutseq,nchar(paste0(altseq5,exon1,intron))-seqgrep[1]+1,nchar(paste0(altseq5,exon1,intron))+seqgrep[2]) )
        normalsplicegrep_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1))-seqgrep[1]+1,nchar(paste0(refseq5,exon1))),
                                           substr(returnseq,nchar(paste0(refseq5,exon1,intron))+1,nchar(paste0(refseq5,exon1,intron))+seqgrep[2]))
        normalsplicegrep_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1))-seqgrep[1]+1,nchar(paste0(altseq5,exon1))),
                                           substr(returnmutseq,nchar(paste0(altseq5,exon1,intron))+1,nchar(paste0(altseq5,exon1,intron))+seqgrep[2]))
        aberrantsplicegrep_normalseq<-"No_Expected_Aberrant_Splice"
        aberrantsplicegrep_mutantseq<-"No_Expected_Aberrant_Splice"
        newframe<-data.frame(Var_Name=mut$name,Kind="Broken_Canonical3SS",First_Exon=exon,Mutation_location=mutpos,Normal_SpliceGap=PredictedNormalSplicedOutPos,Aberrant_SpliceGap=PredictedAberantSplicedOutPos,Summary_Seq=summaryseq,
                             Ref_Construct=returnseq,Ref_Barcode=refcurrentbarcode,Ref_ID=rid,Ref_Content=content,Alt_Construct=returnmutseq,Alt_Barcode=altcurrentbarcode,Alt_ID=aid,Alt_Content=content,
                             Ref_MaxEntScoreSS3=result3csds$Ref_MaxEntScoreSS3,Ref_SS3score=result3csds$Ref_SS3score,Alt_MaxEntScoreSS3=result3csds$Alt_MaxEntScoreSS3,Alt_SS3score=result3csds$Alt_SS3score,Diff_MaxEntScoreSS3=result3csds$Diff_MaxEntScoreSS3,Diff_SS3score=result3csds$Diff_SS3score,
                             grep_NormalSeq=normalseqgrep,grep_MutantSeq=mutantseqgrep,grep_NoSplice_In_NormalSeq=nosplicegap_normalseq,NoSplice_In_MutantSeq=nosplicegap_mutantseq,
                             grep_NormalSplice_In_NormalSeq=normalsplicegrep_normalseq,grep_NormalSplice_IN_MutantSeq=normalsplicegrep_mutantseq,grep_AberrantSplice_In_NormalSeq=aberrantsplicegrep_normalseq,grep_AberrantSplice_IN_MutantSeq=aberrantsplicegrep_mutantseq,
                             Ref_Norm_Splice=PredictedNormalSplice,Ref_Norm_Splice_Length=PredictedNormalLength,Ref_Aberrant_Splice=PredictedAberrantSplice,Ref_Aberrant_Splice_Length=PredictedAberrantLength,
                             Mut_Norm_Splice=PredictedNormalSpliceMut,Mut_Norm_Splice_Length=PredictedNormalLengthMut,Mut_Aberrant_Splice=PredictedAberrantSpliceMut,Mut_Aberrant_Splice_Length=PredictedAberrantLengthMut)
        returnframe<-rbind(returnframe,newframe)
        cat("Broken_Canonical3SS",content,"\n")
      }
      if (nrow(returnframe)==0) return  (data.frame(Var_Name=NA ) ) else return(returnframe)
    } else {  #in case that mutation on an intron  **should consider both new 5SS and broken 5SS
      #mut on intron, exon1,intronbeforemut(+mut), intronaftermut,exon2nd
      #cat("here2\n")
      result3<-takeresult3(refresult3temp,altresult3temp,0)  #check for new 5SS
      ssgap<-as.numeric(result3$pos)+19
      result3csds<-takeresult3(refresult3temp,altresult3temp,1) #check for broken 5SS
      returnframe<-data.frame()
      PredictedNormalSplice<-paste(tolower(refseq5),toupper(exon1),"[",tolower(intronbeforemut),tolower(intronaftermut),"]",toupper(exon2nd),tolower(refseq3),sep="")
      intronbeforemutWithMut<-paste(substr(intronbeforemut,1,nchar(intronbeforemut)-1),mut$alt,sep="")
      PredictedNormalSpliceMut<-paste(tolower(altseq5),toupper(exon1),"[",tolower(intronbeforemutWithMut),tolower(intronaftermut),"]",toupper(exon2nd),tolower(altseq3),sep="")
      PredictedNormalLength<-nchar(returnseq)-nchar(intronbeforemut)-nchar(intronaftermut)
      PredictedNormalLengthMut<-nchar(returnmutseq)-nchar(intronbeforemutWithMut)-nchar(intronaftermut)
      rid<-substr(paste(exon2nd,refseq3,sep=""),nchar(exon2nd)+1-nid[1],nchar(exon2nd)+nchar(barcodes[1])+nid[2])
      aid<-substr(paste(exon2nd,altseq3,sep=""),nchar(exon2nd)+1-nid[1],nchar(exon2nd)+nchar(barcodes[1])+nid[2])
      PredictedNormalSplicedOutPos<-paste0(nchar(paste0(refseq5,exon1))+1,"_",nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut)))
      content<-paste("Attached5Seq:",nchar(refseq5),";1stExon:",aclen[1],";Intron_5:",aclen[2],";Intron_On_Before_Mutation:",aclen[3],";Intron_After_Mutation:",aclen[4],";Intron_Before_Canonical_AccepterSite:",aclen[5],";2ndExon:",aclen[6],";Attached3Seq:",nchar(refseq3),";Remain:",limit-sum(aclen),sep="")
      
      if (!is.na(result3$pos)) { 
        #temp<-paste(exon2nd,altseq3,sep="") 
        PredictedAberrantSplice<-paste0(substr(returnseq,1,nchar(paste0(refseq5,exon1))),"[",substr(returnseq,nchar(paste0(refseq5,exon1))+1,ssgap),"]",substr(returnseq,ssgap+1,nchar(returnseq)) )
        PredictedAberrantSpliceMut<-paste0(substr(returnmutseq,1,nchar(paste0(altseq5,exon1))),"[",substr(returnmutseq,nchar(paste0(altseq5,exon1))+1,ssgap),"]",substr(returnmutseq,ssgap+1,nchar(returnmutseq)) )
        PredictedAberrantLength<-nchar(returnseq)-ssgap+nchar(paste0(refseq5,exon1))  #PredictedNormalLength-(nchar(paste(refseq5,exon1before,exon1after))-ssgap)+2
        PredictedAberrantLengthMut<-nchar(returnmutseq)-ssgap+nchar(paste0(altseq5,exon1)) #PredictedNormalLengthMut-(nchar(paste(altseq5,exon1before,exon1aftermut))-ssgap)+2
        #(data.frame(pos=refresult$phypos[pos],Ref_MaxEntScoreSS5=refresult$MaxEntScoreSS5[pos],Ref_SS5score=refresult$SS5score[pos],Alt_MaxEntScoreSS5=altresult$MaxEntScoreSS5[pos],Alt_SS5score=altresult$SS5score[pos],Diff_MaxEntScoreSS5=(altresult$MaxEntScoreSS5[pos]-refresult$MaxEntScoreSS5[pos]),Diff_SS5score=(altresult$SS5score[pos]-refresult$SS5score[pos]) ))
        PredictedAberantSplicedOutPos<-paste0(nchar(paste0(refseq5,exon1))+1,"_",ssgap)
        
        if (ssgap<nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))) {  #in case that aberrant acceptor is located before canonical acceptor site   #######check the script in 5SS
          temp<-paste0(substr(returnseq,1,nchar(paste0(refseq5,exon1)) ),"{[",substr(returnseq,nchar(paste0(refseq5,exon1))+1,ssgap),"}",substr(returnseq,ssgap+1,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))),"]",substr(returnseq,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+1,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut,exon2nd))),tolower(IN$seq3))
        } else { #in case that aberrant acceptor is located after canonical acceptor site
          temp<-paste0(substr(returnseq,1,nchar(paste0(refseq5,exon1)) ),"{[",substr(returnseq,nchar(paste0(refseq5,exon1))+1,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))),"]",substr(returnseq,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+1,ssgap),"}",substr(returnseq,ssgap+1,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut,exon2nd))),tolower(IN$seq3))
        }
        add<-2 #shift for {[
        if (mutpos>ssgap) add<-add+1  #shift for }
        if (mutpos>nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))) add<-add+1 #shift for ]
        temp2<-paste0(substr(temp,1,mutpos-1+add),"(",mut$ref,"/",mut$alt,")",substr(temp,mutpos+1+add,nchar(temp)))  #give (mut/alt) info *+add means a shift to rightside because of inserted {[ and ] or }
        summaryseq<-paste0(tolower(IN$seq5),substr(temp2,nchar(IN$seq5)+1,nchar(temp2)))  #reset 5seq
        ###kokomadekita!!!
        normalseqgrep<-substr(returnseq,mutpos-seqgrep[1],mutpos+seqgrep[2])
        mutantseqgrep<-substr(returnmutseq,mutpos-seqgrep[1],mutpos+seqgrep[2])
        nosplicegap_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1))-seqgrep[1]+1,nchar(paste0(refseq5,exon1))+seqgrep[2]),",",
                                      substr(returnseq,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))-seqgrep[1]+1,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]) )
        nosplicegap_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1))-seqgrep[1]+1,nchar(paste0(altseq5,exon1))+seqgrep[2]),",",
                                      substr(returnmutseq,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))-seqgrep[1]+1,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]) )
        normalsplicegrep_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1))-seqgrep[1]+1,nchar(paste0(refseq5,exon1))),
                                           substr(returnseq,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+1,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]))
        normalsplicegrep_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1))-seqgrep[1]+1,nchar(paste0(altseq5,exon1))),
                                           substr(returnmutseq,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))+1,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]))
        aberrantsplicegrep_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1))-seqgrep[1]+1,nchar(paste0(refseq5,exon1))),
                                             substr(returnseq,ssgap+1,ssgap+seqgrep[2]) )
        aberrantsplicegrep_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1))-seqgrep[1]+1,nchar(paste0(altseq5,exon1))),
                                             substr(returnmutseq,ssgap+1,ssgap+seqgrep[2]) )
        returnframe<-data.frame(Var_Name=mut$name,Kind="Creating_Aberrant3SS",First_Exon=exon,Mutation_location=mutpos,Normal_SpliceGap=PredictedNormalSplicedOutPos,Aberrant_SpliceGap=PredictedAberantSplicedOutPos,Summary_Seq=summaryseq,
                                Ref_Construct=returnseq,Ref_Barcode=refcurrentbarcode,Ref_ID=rid,Ref_Content=content,Alt_Construct=returnmutseq,Alt_Barcode=altcurrentbarcode,Alt_ID=aid,Alt_Content=content,
                                Ref_MaxEntScoreSS3=result3$Ref_MaxEntScoreSS3,Ref_SS3score=result3$Ref_SS3score,Alt_MaxEntScoreSS3=result3$Alt_MaxEntScoreSS3,Alt_SS3score=result3$Alt_SS3score,Diff_MaxEntScoreSS3=result3$Diff_MaxEntScoreSS3,Diff_SS3score=result3$Diff_SS3score,
                                grep_NormalSeq=normalseqgrep,grep_MutantSeq=mutantseqgrep,grep_NoSplice_In_NormalSeq=nosplicegap_normalseq,NoSplice_In_MutantSeq=nosplicegap_mutantseq,
                                grep_NormalSplice_In_NormalSeq=normalsplicegrep_normalseq,grep_NormalSplice_IN_MutantSeq=normalsplicegrep_mutantseq,grep_AberrantSplice_In_NormalSeq=aberrantsplicegrep_normalseq,grep_AberrantSplice_IN_MutantSeq=aberrantsplicegrep_mutantseq,
                                Ref_Norm_Splice=PredictedNormalSplice,Ref_Norm_Splice_Length=PredictedNormalLength,Ref_Aberrant_Splice=PredictedAberrantSplice,Ref_Aberrant_Splice_Length=PredictedAberrantLength,
                                Mut_Norm_Splice=PredictedNormalSpliceMut,Mut_Norm_Splice_Length=PredictedNormalLengthMut,Mut_Aberrant_Splice=PredictedAberrantSpliceMut,Mut_Aberrant_Splice_Length=PredictedAberrantLengthMut)
        cat("Creating_Aberrant3SS",content,"\n")
      }
      if (!is.na(result3csds$pos)) {
        #mut on exon canonical splice donor site weighted exon1beforecsds,exon1aftercsds(+mut),introncsds,exon2ndcsds
        #mut on intron, exon1,intronbeforemut(+mut), intronaftermut,exon2nd
        #temp<-paste(exon2nd,altseq3,sep="")
        #In broken canonical 5SS, normal splice would not occur. So just give seq without splicing
        PredictedAberrantSplice<-returnseq
        PredictedAberrantSpliceMut<-returnmutseq
        PredictedAberrantLength<-nchar(returnseq)
        PredictedAberrantLengthMut<-nchar(returnmutseq)
        #(data.frame(pos=refresult$phypos[pos],Ref_MaxEntScoreSS5=refresult$MaxEntScoreSS5[pos],Ref_SS5score=refresult$SS5score[pos],Alt_MaxEntScoreSS5=altresult$MaxEntScoreSS5[pos],Alt_SS5score=altresult$SS5score[pos],Diff_MaxEntScoreSS5=(altresult$MaxEntScoreSS5[pos]-refresult$MaxEntScoreSS5[pos]),Diff_SS5score=(altresult$SS5score[pos]-refresult$SS5score[pos]) ))
        PredictedAberantSplicedOutPos<-"No_Expected_Aberrant_Splice"
        temp<-paste(tolower(IN$seq5),toupper(exon1),"[",tolower(intronbeforemut),tolower(intronaftermut),"]",toupper(exon2nd),tolower(IN$seq3),sep="")
        add<-1 #shift for [
        if (mutpos>nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))) add<-add+1 #shift for ]  ####check 5ss in order to confirm the same procedure conducted in 5ss, maybe not done
        summaryseq<-paste0(substr(temp,1,mutpos-1+add),"(",mut$ref,"/",mut$alt,")",substr(temp,mutpos+1+add,nchar(temp)))  #give (mut/alt) info *+add for inserted [ or ]
        normalseqgrep<-substr(returnseq,mutpos-seqgrep[1],mutpos+seqgrep[2])
        mutantseqgrep<-substr(returnmutseq,mutpos-seqgrep[1],mutpos+seqgrep[2])
        nosplicegap_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1))-seqgrep[1]+1,nchar(paste0(refseq5,exon1))+seqgrep[2]),",",
                                      substr(returnseq,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))-seqgrep[1]+1,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]) )
        nosplicegap_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1))-seqgrep[1]+1,nchar(paste0(altseq5,exon1))+seqgrep[2]),",",
                                      substr(returnmutseq,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))-seqgrep[1]+1,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]) )
        normalsplicegrep_normalseq<-paste0(substr(returnseq,nchar(paste0(refseq5,exon1))-seqgrep[1]+1,nchar(paste0(refseq5,exon1))),
                                           substr(returnseq,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+1,nchar(paste0(refseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]))
        normalsplicegrep_mutantseq<-paste0(substr(returnmutseq,nchar(paste0(altseq5,exon1))-seqgrep[1]+1,nchar(paste0(altseq5,exon1))),
                                           substr(returnmutseq,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))+1,nchar(paste0(altseq5,exon1,intronbeforemut,intronaftermut))+seqgrep[2]))
        aberrantsplicegrep_normalseq<-"No_Expected_Aberrant_Splice"
        aberrantsplicegrep_mutantseq<-"No_Expected_Aberrant_Splice"
        newframe<-data.frame(Var_Name=mut$name,Kind="Broken_Canonical3SS",First_Exon=exon,Mutation_location=mutpos,Normal_SpliceGap=PredictedNormalSplicedOutPos,Aberrant_SpliceGap=PredictedAberantSplicedOutPos,Summary_Seq=summaryseq,
                             Ref_Construct=returnseq,Ref_Barcode=refcurrentbarcode,Ref_ID=rid,Ref_Content=content,Alt_Construct=returnmutseq,Alt_Barcode=altcurrentbarcode,Alt_ID=aid,Alt_Content=content,
                             Ref_MaxEntScoreSS3=result3csds$Ref_MaxEntScoreSS3,Ref_SS3score=result3csds$Ref_SS3score,Alt_MaxEntScoreSS3=result3csds$Alt_MaxEntScoreSS3,Alt_SS3score=result3csds$Alt_SS3score,Diff_MaxEntScoreSS3=result3csds$Diff_MaxEntScoreSS3,Diff_SS3score=result3csds$Diff_SS3score,
                             grep_NormalSeq=normalseqgrep,grep_MutantSeq=mutantseqgrep,grep_NoSplice_In_NormalSeq=nosplicegap_normalseq,NoSplice_In_MutantSeq=nosplicegap_mutantseq,
                             grep_NormalSplice_In_NormalSeq=normalsplicegrep_normalseq,grep_NormalSplice_IN_MutantSeq=normalsplicegrep_mutantseq,grep_AberrantSplice_In_NormalSeq=aberrantsplicegrep_normalseq,grep_AberrantSplice_IN_MutantSeq=aberrantsplicegrep_mutantseq,
                             Ref_Norm_Splice=PredictedNormalSplice,Ref_Norm_Splice_Length=PredictedNormalLength,Ref_Aberrant_Splice=PredictedAberrantSplice,Ref_Aberrant_Splice_Length=PredictedAberrantLength,
                             Mut_Norm_Splice=PredictedNormalSpliceMut,Mut_Norm_Splice_Length=PredictedNormalLengthMut,Mut_Aberrant_Splice=PredictedAberrantSpliceMut,Mut_Aberrant_Splice_Length=PredictedAberrantLengthMut)
        returnframe<-rbind(returnframe,newframe)
        cat("Broken_Canonical3SS",content,"\n")
      }
      if (nrow(returnframe)==0) returnframe<-data.frame(Var_Name=NA )
      return(returnframe)
    }  
  }  #ss3==1 end
}

takebarcodes<-function(barcodegrid,previous,n) {
  returnbarcode<-takeminbarcodes(previous)
  if (n>1) for (i in 2:n) {
    returnbarcode<-c(returnbarcode,takeminbarcodes(returnbarcode))
  }
  return(returnbarcode)
}

takeminbarcodes<-function(previous) {
  if (length(previous)>0) {
    scores<-sapply(previous,function(x){
      temp<-barcodegrid[barcodegrid$Var1==x,]
      temp<-temp[order(as.character(temp$Var2)),]
      return(temp$distance)
    })
    scores<-apply(scores,1,function(x){sum(x)})
    tempbarcodes<-barcodes[order(barcodes)]
    #cat("previous:",previous,"\n")
    #cat(tempbarcodes,",",scores)
    seed<-tempbarcodes[which.min(scores)]
  } else {
    seed<-barcodes[1]
  }
  #cat(seed,"\n")
  return(seed)
}

changeSeq<-function(seq,grep1,grep2,barcodelen,barcode) {
  barcode<-tolower(barcode)
  grep1<-tolower(grep1)
  grep2<-tolower(grep2)
  temp<-regexpr(paste0(grep1,paste(rep(".",barcodelen),collapse=""),grep2),seq)
  if (as.integer(temp)<0) stop("Seq around Barcode:",paste0(grep1,paste(rep(".",barcodelen),collapse=""),grep2)," NOT found in the Seq:",seq,"\n")
  temps<-as.numeric(temp);tempe<-temps+attr(temp,"match.length")-1
  newseq<-paste0(substr(seq,1,temps-1),paste(rep("Z",(nchar(grep1)+nchar(barcode)+nchar(grep2))),collapse=""), substr(seq,tempe+1,nchar(seq)))
  tempdupcheck<-regexpr(paste0(grep1,paste(rep(".",barcodelen),collapse=""),grep2),newseq)
  if (as.integer(tempdupcheck)>0) {
    cat("Seq around Barcode:",paste0(grep1,paste(rep(".",barcodelen),collapse=""),grep2)," found twice in the Seq:",seq,"\n")
    temp3<-regexpr("(X)+",IN$seq3)
    if (as.integer(temp3)>0) {
      cat("Since Barcode is located in the 3ss seq, take the later one.\n")
      temps<-as.numeric(tempdupcheck);tempe<-temps+attr(tempdupcheck,"match.length")-1
    }
    temp5<-regexpr("(X)+",IN$seq5)
    if (as.integer(temp5)>0) {
      cat("Since Barcode is located in the 5SS seq, take the first one.\n")
    }
  }
  
  
  return( paste0(substr(seq,1,temps-1+nchar(grep1)), barcode, substr(seq,tempe-nchar(grep2)+1,nchar(seq))) )  
}

changeID<-function(ID,barcode) {
  barcode<-toupper(barcode)
  return( paste0(substr(ID,1,nid[1]),barcode,substr(ID,nchar(ID)-nid[2]+1,nchar(ID)))  ) 
}

limitIntron<-function(vec) {
  if (IN$limitintron!=-9) {
    if(sum(vec)>IN$limitintron) {
      cat("Intron is longer than the limitintron:",IN$limitintron,"/ The sum of the Intron length:",sum(vec),"\n")
      nvec<-which(vec>30)  #just reduce components with length>30 because ,max length of uniside of consensus seq is 20bp (3ss 20bp+3bp)
      if (length(nvec)>0 & length(nvec)<length(vec) ) {
        #cat("vec:",vec,"nvec:",nvec,"\n")
        alpha<-(IN$limitintron-vec[-nvec])/sum(vec[nvec])
        vec[nvec]<-as.integer(vec[nvec]*alpha)
      } else if (length(nvec)==length(vec) ){
        alpha<-IN$limitintron/sum(vec)
        vec<-as.integer(vec*alpha)
      }
      for (i in 2:length(vec)) {  #give the length if the latter length become 0 when united contents.   
        if (vec[i]==0) {
          half<-as.integer(vec[i-1]/2)
          resthalf<-vec[i-1]-half
          vec[i-1]<-resthalf
          vec[i]<-half
        } 
      }
    }
  } 
  return(vec)
}

dupIDcheck<-function(resultframe) {
  dupids<-data.frame()
  group<-unique(as.character(resultframe$SplitID))
  for (i in 1:length(group)) {
    refpos<-which(as.character(resultframe$SplitID)==group[i])
    refframe<-data.frame(Var_Name=I(c(paste(as.character(resultframe$Var_Name[refpos]),collapse=","),as.character(resultframe$Var_Name[refpos]))),
                         Kind=I(c("Reference_Construct",as.character(resultframe$Kind[refpos]))),
                         ID=I(c(as.character(resultframe$Ref_ID[refpos[1]]),as.character(resultframe$Alt_ID[refpos]))))
    if (i<length(group)) for (j in (i+1):length(group)) {
      compos<-which(as.character(resultframe$SplitID)==group[j])
      compframe<-data.frame(Var_Name=I(c(paste(as.character(resultframe$Var_Name[compos]),collapse=","),as.character(resultframe$Var_Name[compos]))),
                            Kind=I(c("Reference_Construct",as.character(resultframe$Kind[compos]))),
                            ID=I(c(as.character(resultframe$Ref_ID[compos[1]]),as.character(resultframe$Alt_ID[compos]))))
      for (k in 1:nrow(refframe)) {
        hit<-which(compframe$ID==refframe$ID[k])
        if (length(hit)>0) {
          cat("SplitID:",group[i],"Var_Name:",refframe$Var_Name[k],"Kind:",refframe$Kind[k],"ID:",refframe$ID[k]," is equal to SplitID:",group[j],"Var_Name:",compframe$Var_Name[hit],"Kind:",compframe$Kind[hit],"ID:",compframe$ID[hit],"\n")
          temp<-data.frame(ref=group[i],alt=group[j])
          dupids<-rbind(dupids,temp)
        }
      }
    }
  }
  return(dupids)
}


#main
cat("This script is to calculate the integrated score of 5ss or 3ss\n")
inputlist<-list(helpflag="-help",debug="-debug",marthost="-marthost",mutfile="-mutfile",output="-output",seq5="-seq5",seq3="-seq3",barcode="-barcode",randombarcode="-randombarcode",
                readbarcodes="-readbarcodes",esbTranscriptID="-esbTranscriptID",constructlength="-constructlength",ss5="-ss5",ss3="-ss3",refdirectory="-refdirectory",id="-id",seqgrep="-seqgrep",
                summarizeConstructs="-summarizeConstructs",limitintron="-limitintron",limit2ndexon="-limit2ndexon",limit1stexon="-limit1stexon",
                gblocksubmit="-gblocksubmit",scoreresult="-scoreresult")
inputparas<-c(0,0,1,1,1,1,1,1,0,
              1,1,1,0,0,1,1,1,
              1,1,1,1,1,1)  #the length of parameters to be input in the command line 0-switch yes or no 1- one variable 2-two variables
IN<-inputparameters(inputlist,inputparas)
if (IN$helpflag==1) stop("\n",Helpdocument)
if (IN$output==-9) IN$output<-"ConstructDesigner.v0.9.R.result.txt"
if (IN$ss5==-9 & IN$ss3==-9 & IN$summarizeConstructs==-9) stop("Please indicate -ss5 or -ss3\n")
if (IN$ss5>0 & IN$ss3>0) stop("Please indicate one of -ss5 or -ss3\n")
if (IN$mutfile==-9 & IN$summarizeConstructs==-9) stop("Please indicate a mutation file by -mutfile\n")
if (IN$refdirectory==-9) IN$refdirectory<-"./ssfiles"
if (IN$constructlength==-9) {
  if (IN$ss5==1) IN$constructlength<-"85,85,85,85,55,40"
  if (IN$ss3==1) IN$constructlength<-"40,55,85,85,85,85"
}
if (IN$seq5==-9) IN$seq5<-""
if (IN$seq3==-9) IN$seq3<-""
if (IN$id==-9) IN$id<-"0,0"
if (IN$seqgrep==-9) IN$seqgrep<-"8,8"
IN$limitintron<-as.numeric(IN$limitintron)
IN$limit2ndexon<-as.numeric(IN$limit2ndexon)
IN$limit1stexon<-as.numeric(IN$limit1stexon)
cat("Parameters:\n")
for (i in names(IN)) {
  cat(i,":",unlist(IN[which(names(IN)==i)]),"\n")
}

nid<-as.integer(unlist(strsplit(IN$id,",")))
seqgrep<-as.integer(unlist(strsplit(IN$seqgrep,",")))
#barcode construction
barcodesYES<-0  #if employ barcode or not
barcodes<-c()  #barcodes to be employed
barcodesPOS<-1   #pointer in the barcodes
if (IN$barcode!=-9) {
  barcodes<-unlist(strsplit(IN$barcode,","))
  if (length(barcodes)==0) stop("Give barcodes in the -barcode switch.\n")
  barcodesYES<-1
}
if (IN$readbarcodes!=-9) {
  barcodes<-readLines(con=IN$readbarcodes)
  if (length(barcodes)==0) stop("Indicate a file that contains barcodes in the -readbarcodesswitch.\n")
  barcodesYES<-1
}
if (barcodesYES & IN$randombarcode!=-9) {
  cat("The order of barcodes are being shuffled for randamization...\n")
  barcodes<-sample(barcodes,length(barcodes),replace=F)
}
if (barcodesYES) {
  temp3<-length(grep("X",IN$seq3))
  temp5<-length(grep("X",IN$seq5))
  if ((temp3+temp5)==0) stop("Please indicate where barcodes should be placed in the 5seq or the 3seq.\n")
  cat("The barcodes we employ are as follows\n")
  for (i in barcodes) {
    cat(i,"\t")
  }
  cat("\n")
}



if (IN$summarizeConstructs=="current" | IN$summarizeConstructs==-9) {

conlen<-as.numeric(unlist(strsplit(IN$constructlength,",")))
#cat("Construct Length:",conlen,"\n")
#initiation for using biomaRT
library(biomaRt) #if you need downloading biomaRt, please type as follows
cat("Connecting biomaRt server...\n")
marthostname<-ifelse(IN$marthost==-9,"feb2014.archive.ensembl.org",IN$marthost) #feb2014.archive.ensembl.org   is hg19 latest ver host
cat("Host name is",marthostname,"\n")
mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host=marthostname)  # biomart for gene annotation
#host="may2009.archive.ensembl.org",  hg18
#host="Sep2013.archive.ensembl.org", hg19

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

#making mutation frame
muttemp<-readLines(con=IN$mutfile)
mut.frame<-data.frame()
if (length(muttemp)>0) {
  if (length(grep("c\\.",muttemp))>0) mut.frame<-cDNAmut(muttemp)
  if (length(grep("chr",muttemp))>0) mut.frame<-chrposmut(muttemp) #new.frame<-data.frame(kind=I("chrpos"),chr=I(chr),pos=as.integer(pos),ref=I(ref),alt=I(alt),name=I(paste("chr",chr,":",pos,ref,">",alt,sep="")),transID=I(transID))
} else {
  stop("The length of mutation file is zero.\n")
}
if (nrow(mut.frame)==0)  stop("Given mutation file is invalid.\n")
reverseflag<-0
lseq<-data.frame()
newcDNApos<-c()
exon.frame<-data.frame()
aclen<-c(0,0,0,0,0,0)
resultframe<-data.frame()
#mut.frame loop
preTransID<-"dummy"
for (i in 1:nrow(mut.frame)) {  #main loop for mut.frame
    cat("Processing",mut.frame$name[i],"...\n")
    #get sequence and position of the mutation in the seq
    getseq<-0;if (!is.na(mut.frame$transID[i])) if (mut.frame$transID[i]!=preTransID) getseq<-1
    if (nrow(lseq)==0 | getseq) {lseq<-GetSeq(mut.frame$transID[i],1);preTransID<-mut.frame$transID[i]}
    if (length(grep("coding",mut.frame[i,]$kind))>0) newcDNApos<-cDNAposdetermine(mut.frame[i,]$pos) #which(lseq$cDNApos==mut.frame[i,]$pos)
    if (mut.frame[i,]$kind=="chrpos" ) {newcDNApos<-which(lseq$phypos==mut.frame[i,]$pos);newcDNApos<-c(newcDNApos,newcDNApos);}
    if (reverseflag==1 & mut.frame[i,]$kind=="chrpos") {
      cat("ChrPos information in a reverseStarnd gene. Therefore, change the allele info for cDNA.\n")
      mut.frame[i,]$ref<-reversenuc(mut.frame[i,]$ref)
      mut.frame[i,]$alt<-reversenuc(mut.frame[i,]$alt)
    }
    #generate ref construct
    #write.table(lseq,file="temp.lseq.txt",sep="\t",col.names=T,row.names=F,quote=F)
    const<-constructmake(lseq,newcDNApos[1],mut.frame[i,],conlen)
    const<-cbind(data.frame(TranscriptID=rep(mut.frame$transID[i],nrow(const))),const)
    #write.table(const,file="temp.const.txt",sep="\t",col.names=T,row.names=F,quote=F)
    for (numcon in 1:nrow(const)) {
      if (!is.na(const$Var_Name[numcon])) {
        resultframe<-rbind(resultframe,const[numcon,])
      }
    }
}  #the end of for (i in 1:nrow(mut.frame)) {  #main loop for mut.frame 


} else { #if (IN$summarizeConstructs== current | IN$summarizeConstructs==-9) end
  files<-unlist(strsplit(IN$summarizeConstructs,","))
  resultframe<-data.frame()
  for (file in files) {
    cat("Reading the result file:",file," for summarizing\n")
    new<-read.table(file,header=T,stringsAsFactors=F,sep="\t")
    resultframe<-rbind(resultframe,new)
  }
}

cat("Writing a result to ",IN$output,"\n")
write.table(resultframe,file=IN$output,sep="\t",quote=F,col.names=T,row.names=F)

if (IN$summarizeConstructs!=-9) {
  cat("\nSummarizing Constructs\nFirst, deleting duplicated lines...\n")
  resultframe<-read.table(IN$output,sep="\t",header=T,stringsAsFactors = F)
  delrow<-c()
  for (i in 1:nrow(resultframe)) {
    temp<-which(resultframe$TranscriptID==resultframe$TranscriptID[i] 
                & resultframe$Kind==resultframe$Kind[i] 
                & resultframe$First_Exon==resultframe$First_Exon[i] 
                & resultframe$Mutation_location==resultframe$Mutation_location[i] 
                & resultframe$Summary_Seq==resultframe$Summary_Seq[i])
    temp<-temp[temp>i]
    if (length(temp)>0) {
      cat("Line#:",i,"Variant:",resultframe$Var_Name[i],resultframe$Kind[i], "is the same as Line#:",temp,"Variant:",resultframe$Var_Name[temp],resultframe$Kind[temp],"\n")
      delrow<-c(delrow,temp) 
    }
  }
  if (length(delrow)>0) resultframe<-resultframe[-delrow,]
  
  cat("Summarizing and rearranging the results...\n")
  barcodegrep1<-"" #beforebarcode
  barcodegrep2<-"" #afterbarcode
  barcodelen<-0
  #if (barcodesYES) {   
  #prepare seqs around barcode to replace barcode
  barcodelen<-nchar(barcodes[1])
  temp3<-regexpr("(X)+",IN$seq3)
  temp3s<-as.numeric(temp3);temp3e<-temp3s+attr(temp3,"match.length")-1
  if (temp3s>0) {
    barcodegrep1<-substr(IN$seq3,1,temp3s-1)
    barcodegrep2<-substr(IN$seq3,temp3e+1,nchar(IN$seq3))
  }
  
  temp5<-regexpr("(X)+",IN$seq5)
  temp5s<-as.numeric(temp5);temp5e<-temp5s+attr(temp5,"match.length")-1
  if (temp5s>0) {
    barcodegrep1<-substr(IN$seq5,1,temp5s-1)
    barcodegrep2<-substr(IN$seq5,temp5e+1,nchar(IN$seq5))
  }
  #prepare bacode matrix to calculate the distance
  barcodegrid<-expand.grid(barcodes,barcodes)
  barcodegrid$distance<-0
  barcodegrid$distance<-apply(barcodegrid,1,function(x) {
    c1<-x[1];c2<-x[2]
    penalty<-0
    for (i in 1:nchar(c1)) {
      if (substr(c1,i,i)==substr(c2,i,i)) penalty<-penalty+1
    }
    if (penalty==nchar(c1)) penalty<-nchar(c1)*1000 #if all the same, give the most penalty
    return(penalty)
  })  #a<-barcodegrid[barcodegrid$Var1=="AA",]   a<-a[order(a$distance),] #take distance==0 because of no overlap
  #}
  
  resultframe<-cbind(resultframe,split=paste0(resultframe$TranscriptID,"_",resultframe$First_Exon,"_",resultframe$Normal_SpliceGap),split2=paste0(substr(resultframe$Ref_ID,1,nid[1]),substr(resultframe$Ref_ID,nchar(resultframe$Ref_ID)-nid[2]+1,nchar(resultframe$Ref_ID)))) #exclude resultframe$Normal_SpliceGap
  splitframe<-split(resultframe,resultframe$split2) #split by ID
  resultframe<-data.frame()
  n_elements<-0
  usedbarcodes<-c()
  for (i in 1:length(splitframe)){
    temp<-splitframe[[i]] #pick out lines grouped by transcriptID, firstexona dn splicegap
    #prepare one barcode for ref and the same n of barcodes for mut construct
    n_ref<-unique(as.character(temp$split))
    cat("Split group by ID without barcode:",as.character(temp$split2[1]),"   # of reference constructs:",length(n_ref) ,"   # of mutant constructs:",nrow(temp),"\n")
    cbarcodes<-takebarcodes(barcodegrid,previous=usedbarcodes,n=(nrow(temp)+length(n_ref)))  #pick barcodes
    n_elements<-n_elements+nrow(temp)+length(n_ref)
    for (j in 1:length(n_ref)) {
      n<-which(temp$split==n_ref[j])
      n1<-n[1]
      Ref_Barcode<-cbarcodes[j]
      Ref_Construct<-changeSeq(temp$Ref_Construct[n1],barcodegrep1,barcodegrep2,barcodelen,Ref_Barcode)
      Ref_Norm_Splice<-changeSeq(temp$Ref_Norm_Splice[n1],barcodegrep1,barcodegrep2,barcodelen,Ref_Barcode)
      Ref_Aberrant_Splice<-changeSeq(temp$Ref_Norm_Splice[n1],barcodegrep1,barcodegrep2,barcodelen,Ref_Barcode)
      Ref_ID<-changeID(temp$Ref_ID[n1],Ref_Barcode)
      
      temp$Ref_Barcode[n]<-Ref_Barcode
      temp$Ref_Construct[n]<-Ref_Construct
      temp$Ref_Norm_Splice[n]<-Ref_Norm_Splice
      temp$Ref_Aberrant_Splice[n]<-Ref_Aberrant_Splice
      temp$Ref_ID[n]<-Ref_ID
    }
    for (j in 1:nrow(temp)) {
      temp$Alt_Barcode[j]<-cbarcodes[j+length(n_ref)]
      temp$Alt_Construct[j]<-changeSeq(temp$Alt_Construct[j],barcodegrep1,barcodegrep2,barcodelen,temp$Alt_Barcode[j])
      temp$Mut_Norm_Splice[j]<-changeSeq(temp$Mut_Norm_Splice[j],barcodegrep1,barcodegrep2,barcodelen,temp$Alt_Barcode[j])
      temp$Mut_Aberrant_Splice[j]<-changeSeq(temp$Mut_Norm_Splice[j],barcodegrep1,barcodegrep2,barcodelen,temp$Alt_Barcode[j])
      temp$Alt_ID[j]<-changeID(temp$Alt_ID[j],temp$Alt_Barcode[j])
    }
    resultframe<-rbind(resultframe,temp)
    usedbarcodes<-c(usedbarcodes,cbarcodes)
  }
  resultframe<-cbind(data.frame(SplitID=resultframe$split),resultframe)  
  resultframe<-resultframe[,-which(names(resultframe)%in%c("split","split2"))]
  resultframe<-resultframe[order(resultframe$Var_Name),]
  resultframe<-resultframe[order(resultframe$Normal_SpliceGap),]
  resultframe<-resultframe[order(as.numeric(resultframe$First_Exon)),]
  resultframe<-resultframe[order(resultframe$TranscriptID),]
  cat("#unique ref_ID:",length(unique(resultframe$Ref_ID)),"   #unique mut_ID:",length(unique(resultframe$Alt_ID)), "   #given barcodes:",n_elements,"\n")
  #checking uniqueness of barcodes
  if (n_elements==length(unique(c(as.character(resultframe$Ref_ID),as.character(resultframe$Alt_ID))))) cat("Redundancy check passed.\n") else {
    cat("Cauation!!! Redundancy found in constructIDs.\n")
    dupIDs<-dupIDcheck(resultframe)
    
  }
  cat("Writing a summarized result to ",IN$output,"\n")
  write.table(resultframe,file=IN$output,sep="\t",quote=F,col.names=T,row.names=F)
  
  
}

if (IN$gblocksubmit!=-9) {  #making a file for gBlock Order
  file<-read.table(IN$output,sep="\t",header=T,stringsAsFactors = F)
  cat("\nMaking a file for gBlock Order, where you can find construct names and seqs...\n")
  newframe<-data.frame()
  for (i in 1:nrow(file)) {
    name<-gsub(":|\\.|>|\\+|\\-","_",file$Var_Name[i])
    new<-data.frame(Name=I(paste0("REF_",name)),ID=I(file$Ref_ID[i]),seq=I(file$Ref_Construct[i]))
    newframe<-rbind(newframe,new)
    new<-data.frame(Name=I(paste0("ALT_",name)),ID=I(file$Alt_ID[i]),seq=I(file$Alt_Construct[i]))
    newframe<-rbind(newframe,new)
  }
  
  new2frame<-data.frame()
  skips<-c()
  for (i in 1:nrow(newframe)) {
    temp<-which(skips==i)
    if (length(temp)>0) next
    dupcheck<-which(newframe$ID==newframe$ID[i])
    dupcheck<-dupcheck[dupcheck>i]
    if (length(dupcheck)>1) {
      skips<-c(skips,dupcheck)
      names<-c()
      for (j in dupcheck) {
        names<-c(names,sub("(REF_|ALT_)","",newframe$Name[j]))
      }
      newframe$Name[i]<-paste0(newframe$Name[i],"_",paste(names,collapse="_"))
    } 
    new<-data.frame(Name=newframe$Name[i],Seq=newframe$seq[i])
    new2frame<-rbind(new2frame,new)
  }
  
  cat("Deleting duplicated lines...\n")
  delrow<-c()
  for (i in 1:nrow(new2frame)) {
    temp<-which(new2frame$Name==new2frame$Name[i] )
    temp<-temp[temp>i]
    if (length(temp)>0) {
      cat("Line#:",i,"Variant:",new2frame$Name[i], "is the same as Line#:",temp,"Variant:",new2frame$Name[temp],"\n")
      delrow<-c(delrow,temp) 
    }
  }
  if (length(delrow)>0) new2frame<-new2frame[-delrow,]
  
  cat("Writing a gBlock Order File to:",IN$gblocksubmit,"\n")
  write.table(new2frame,file=IN$gblocksubmit,sep="\t",col.names=T,row.names=F,quote=F)
  cat("Done.\n")
}
warnings()

