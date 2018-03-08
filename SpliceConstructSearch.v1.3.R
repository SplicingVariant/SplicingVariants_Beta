Helpdocument<-"
#this scirpt is to find barcode from fastq, then find corresponding seq in the another pair-end read, finally merge them and align it to the ref seq.
#
#Usage: R --slave --vanilla --args -inputframe (filename of inputframe) -fastqF (R1 fastq file) -fastqR (R2 fastq file) -output (output file name) <SpliceConstructSearch.v0.5.R
# 
#ex) R --slave --vanilla --args -inputframe Seq200inputframe.txt -fastqF /groups/seidman/solexa/bpf/150313_M00455_0336_ADY04_FC01415/Kaoru/fastq/150313-FC01415-L1-NI-$I--Kaoru--LMNA--Hs--R1.fq -fastqR /groups/seidman/solexa/bpf/150313_M00455_0336_ADY04_FC01415/Kaoru/fastq/150313-FC01415-L1-NI-$I--Kaoru--LMNA--Hs--R2.fq -output 150313-FC01415-L1-NI-$I.v0.7.output.txt -centerpenalty -displayalign -takethreshold -80 -nogapbonus < SpliceConstructSearch.v0.7.R
#ex) R --slave --vanilla --args -inputframe Seq200inputframe.txt -fastqF /groups/seidman/solexa/bpf/150313_M00455_0336_ADY04_FC01415/Kaoru/fastq/150313-FC01415-L1-NI-$I--Kaoru--LMNA--Hs--R1.fq -fastqR /groups/seidman/solexa/bpf/150313_M00455_0336_ADY04_FC01415/Kaoru/fastq/150313-FC01415-L1-NI-$I--Kaoru--LMNA--Hs--R2.fq -output 150313-FC01415-L1-NI-$I.v0.7.output.revcomp.txt -centerpenalty -displayalign -takethreshold -80 -nogapbonus -reversecomp < SpliceConstructSearch.v0.7.R
#ex) R --slave --vanilla --args -inputframe Seq200inputframe.txt -fastqF /groups/seidman/solexa/bpf/150313_M00455_0336_ADY04_FC01415/Kaoru/fastq/150313-FC01415-L1-NI-$I--Kaoru--LMNA--Hs--R1.fq -fastqR /groups/seidman/solexa/bpf/150313_M00455_0336_ADY04_FC01415/Kaoru/fastq/150313-FC01415-L1-NI-$I--Kaoru--LMNA--Hs--R2.fq -output 150313-FC01415-L1-NI-$I.v0.9.output.revcomp.txt -centerpenalty -displayalign -valiablethreshold 0.8 < SpliceConstructSearch.v0.9.R
#ex) R --slave --vanilla --args -inputframe GCGGCGAA_Seq500inputframe.txt -fastqF testR1.fq -fastqR testR2.fq -output GCGGCGAA_Seq500inputframe.txt.LIB015021_GEN00036527_CTTGTA.fastq.v1.1.output.txt -anonymousspliceassign -reversecomp< SpliceConstructSearch.v1.3.R
#ex) R --slave --vanilla --args -inputframe GCGGCGAA_Seq500inputframe.txt -fastqF testR1.fq -fastqR testR2.fq -output GCGGCGAA_Seq500inputframe.txt.LIB015021_GEN00036527_CTTGTA.fastq.v1.1.output.txt -anonymousspliceassign -directReads < SpliceConstructSearch.v1.3.R
#
##the format of inputframe should be as follows
#Barcode\tSeq\n
#(actual barcode)\t(actual refseq)\n
#and it should be tab-delimted
#
#-swithches
#-help shows help and quit.
#-debug indicates to save internal variables when jobs are excuted
#-inputframe (filename) indicates the input data frame
#-fastqF (filemname) indicates the fastq R1 file
#-fastqR (filename) indicates the fastq R2 file
#-output (filename) indicates the output file name
#-mergethreshold (integer) indicates how many charaters should be overlapped between fastqF and fastqR fragemnt when they are merged
#                 default value:10
#-centerpenalty indicates to give penalty score when found gap in the middle of the aligned sequneces
#-takethreshold (integer) indicates if the alignment score is greater than this value, assign the fragment to the ref score
#                 defalt vakue is :-60  Optimal value should be varied according to the lengths of ref seqs and given seqs.
#                 roughly calculated as follows  -10*2-2*(nchar(min-ref seq)-readlength)+1*(readlength)*(%exprected match )-3*(readlength)*(1-%exprected match )  ##-10 gapopening penality*2both side ;-2 gap extension penalty; +1 match score
#-valiablethreshold (ratio of expected match 0-1 e.g. 0.8 ) indicates threshold is calculated by each ref seq sets like -10*2-2*(nchar(min-ref seq)-readlength)+1*(readlength)*(%exprected match )-3*(readlength)*(1-%exprected match ) 
#-displayalign indicates to show alignment result and score to modulate the -takethreshold 
#-nogapbonus indicagtes to give bonus point if the subject seq does not show gap ( show little gap)
#-permitgapsize (integer) indicates how long gap should be permitted when conducting -centerpenalty and -nogapbonus
#                 default value is 3
#-reversecomp indicates that reverse comp barcode are also employed in seq search (in addition to given barcode) 
#-unknowncollect indicates to save reads assigned to unknown  filename is outputfilename.BARCODEunknownseq.txt
#
#-anonymousspliceassign indicates that
#                           1, employ just reference seq in the 1st column
#                           2, assign merged reads to the above refseq
#                           3, find 2 blocks which should be longer than indicated -blocklength (integer)
#                           4, check intrnal gap point    like aaa---ttt  gap4-6 
#                           5, finally correct data in 4. and groups reads according to internal gap
#                           6, output representive seq and the frequency, and the acutual seqs (the former and the later saved separately) 
#-blocklength (integer) indicates the length of 1 block in -anonymousspliceassign switch . default value is 10
#-blockaccuracy (0-1) indicates hwo accurate the block are compared to corresponding ref seq. default value is 0.95
#-unknownanonymousspliceassign indicates that the script conduct -anonymousspliceassign for unknown-assigned reads
#-inputreadsanonymousspliceassign (filename) indicates that the script conduct -anonymousspliceassign for read in (filename)
#                                     #NOTE THAT 1st seq of (filename) WILL BE EMPLOYED AS REF SEQ
#History
#v0.5 launch ver.
#v0.6 change matching algorithm  break in the center -> big minus
#v0.7 bug-fixed , add swithes     
#v0.8 reverse-complementary mode implemented
#v0.9 consider the length difference between ref seqs (fix the tendency that longer gets more penalty in alignment)
#     re-assign takethroshold
#     when tied score found, then Freq_unknown will increase.
#v0.91 in the request of Jon, unknown seq collecter implemented
#v1.0  anonymous splice assginment mode implement
#      change reverse complement mode  (old) just deal with reverse comp (new) deal with forward AND REVERSE COMP
#v1.1 when barcode search and merging procedure, do not discard paired seq when merging failed
#v1.2 because retreiveing seq procedure from fastq files takes much time when big fastq given, perl wrapper SpliceCOnstructSearchWrapperV2.pl made and specific switch for it has been added.
#     -directReads indicates reads already aligned reads from barcode_fastqR1 and barcode_fastqR2 files
#     for parallel computing using directReads, please use SpliceConstructSearchWrapperV2.pl
# ex)perl SpliceConstructSearchWrapperV2.pl GCGGCGAA_Seq500inputframe.txt testR1.fq testR2.fq R --slave --vanilla --args -inputframe <inputframe> -fastqF testR1.fq -fastqR testR2.fq -output <inputframe>.LIB015021_GEN00036527_CTTGTA.fastq.v1.1.output.txt -anonymousspliceassign -directReads < SpliceConstructSearch.v1.2.R
#     note that -directReads mode requires you to indicate -fastqF and -fastqR, because a filename which contains directReads should be barcode_fastqF pr barcode_fastqR, which means these files should be located at the current directory
#v1.3 bug-fixed in v1.1update, which mean v1.1 and v1.2 usual mode (except -directReads) doesn not work.
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

penaltycalc<-function(subject){
  temp<-unlist(strsplit(subject,""))
  len<-length(temp)
  start<-0
  accum<-0
  minus<-data.frame() #find blank
  for (i in 1:len) {
    if (temp[i]=="-") {
      if(start==0) start<-i
      if (i==len) {
        new<-data.frame(start=start,end=i)
        minus<-rbind(minus,new)
      }
    } else {
      if(start>0) {
        new<-data.frame(start=start,end=i-1)
        minus<-rbind(minus,new)
        start<-0
      }
    }
  }
  #calculate score
  score<-0
  
  if (nrow(minus)>0 ) {
    for (i in 1:nrow(minus)) {
      start<-minus$start[i]
      end<-minus$end[i]
      if ((end-start+1)>IN$permitgapsize) {
        loc<-(start+end)/2
        score<-score+2000*dnorm(loc,mean=len/2,sd=len/10)
      }
    }
  }
  return(score)
}


searchpairedreads<-function(i,barcode){
cat("Searching Sequences including #",i,"Barcode:",barcode,"\n")
if (IN$directReads==-9) {
hit1<-grep(barcode,fastqF)
hit2<-grep(barcode,fastqR)
cat("#hit sequences is",length(hit1), "in R1 and",length(hit2)," in R2\nThen Searching Sequences in the Corresponding pair-end seqs\n")
#find corresponding pair-end read
if (length(hit1)>0) {
  hit1rev<-sapply(hit1,function(x) {
    temp<-unlist(strsplit(fastqF[x-1],":"))
    ID1<-temp[6]
    ID2<-unlist(strsplit(temp[7]," "))[1]
    ID<-paste(ID1,ID2,sep=":")
    temp2<-grep(ID,fastqR)
    rev<-sapply(temp2,function(x){
      temp<-unlist(strsplit(fastqR[x],":"))
      ID1<-temp[6]
      ID2<-unlist(strsplit(temp[7]," "))[1]
      IDrev<-paste(ID1,ID2,sep=":")
      return(ifelse(ID==IDrev,x+1,NA))
    })
    rev<-rev[!is.na(rev)]
    return(ifelse(length(rev)==0,NA,rev[1]))
  })
  cat("Found the corresponding pair-end reads #", length(hit1rev[!is.na(hit1rev)]),"for R1 #",length(hit1),"\n")
  #merge F & R
  cat("Merging Sequences from R1\n")
  mergeseq1<-apply(data.frame(hitF=hit1,hitR=hit1rev),1,function(x) {
    seq1<-fastqF[x[1]]
    if (!is.na(x[2])) {# if paired seq exist ,then merge
      seq2<-reversenuc(fastqR[x[2]])
      temp<-regexpr(barcode,seq1) #instead of variant position, actual 5ss gap taken
      start1<-as.numeric(temp)
      end1<-as.numeric(temp)+attr(temp,"match.length")-1
      local<-pairwiseAlignment(DNAString(seq1), DNAString(seq2), type = "local", substitutionMatrix = mat,gapOpening = -5, gapExtension = -2)
      common<-as.character(attr(local,"pattern"))
      #cat("R1 seq:",seq1," Paird R2 Seq:",seq2," Common Seq:",common,"\n")
      if (nchar(common)>=IN$mergethreshold) { #if common string is longer than indicated
        temp<-regexpr(common,seq1) 
        start1<-as.numeric(temp)
        end1<-as.numeric(temp)+attr(temp,"match.length")-1
        temp<-regexpr(common,seq2) 
        start2<-as.numeric(temp)
        end2<-as.numeric(temp)+attr(temp,"match.length")-1
        if (min(start2,end2)>0) {  #if actually found common seq in seq2     i.e. negative value means not found in the seq
          length1<-nchar(seq1)
          length2<-nchar(seq2)
          if (start1<start2) { #take seq2 first then take seq1
            endseq<-end1+1
            if (endseq>length1) endseq<-length1
            seq1<-paste(substr(seq2,1,end2),substr(seq1,endseq,length1),sep="")
            seq2<-NA
          } else { #take seq1 first then take seq2
            endseq<-end2+1
            if (endseq>length2) endseq<-length2
            seq1<-paste(substr(seq1,1,end1),substr(seq2,endseq,length2),sep="")
            seq2<-NA
          }
          #cat("Merged Seq:",seq1,"\n")
        } else {
          #seq1<-c(seq1,seq2) #cat("No merge\n") #in case of no merge, reply seq1 and seq2 seprately
        }
      } else {
        #seq1<-c(seq1,seq2) #cat("No merge\n") #in case of no merge, reply seq1 and seq2 seprately
      } 
    } 
    return(list(seq1=seq1,seq2=seq2))
  })
  mergeseq1<-unlist(mergeseq1)
  mergeseq1<-mergeseq1[!is.na(mergeseq1)]
}
if (length(hit2)>0) {
  hit2rev<-sapply(hit2,function(x) {
    temp<-unlist(strsplit(fastqR[x-1],":"))
    ID1<-temp[6]
    ID2<-unlist(strsplit(temp[7]," "))[1]
    ID<-paste(ID1,ID2,sep=":")
    temp2<-grep(ID,fastqF)
    rev<-sapply(temp2,function(x){
      temp<-unlist(strsplit(fastqF[x],":"))
      ID1<-temp[6]
      ID2<-unlist(strsplit(temp[7]," "))[1]
      IDrev<-paste(ID1,ID2,sep=":")
      return(ifelse(ID==IDrev,x+1,NA))
    })
    rev<-rev[!is.na(rev)]
    return(ifelse(length(rev)==0,NA,rev[1]))
  })
  cat("Found the corresponding pair-end reads #", length(hit2rev[!is.na(hit2rev)]),"for R2 #",length(hit2),"\n")
  #merge F & R
  cat("Merging Sequences from R2\n")
  mergeseq2<-apply(data.frame(hitF=hit2,hitR=hit2rev),1,function(x) {
    seq1<-fastqR[x[1]]
    if (!is.na(x[2])) {# if paired seq exist ,then merge
      seq2<-reversenuc(fastqF[x[2]])
      temp<-regexpr(barcode,seq1) #instead of variant position, actual 5ss gap taken
      start1<-as.numeric(temp)
      end1<-as.numeric(temp)+attr(temp,"match.length")-1
      local<-pairwiseAlignment(DNAString(seq1), DNAString(seq2), type = "local", substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)
      common<-as.character(attr(local,"pattern"))
      #cat("R1 seq:",seq1," Paird R2 Seq:",seq2," Common Seq:",common,"\n")
      if (nchar(common)>=IN$mergethreshold) {
        temp<-regexpr(common,seq1) 
        start1<-as.numeric(temp)
        end1<-as.numeric(temp)+attr(temp,"match.length")-1
        temp<-regexpr(common,seq2) 
        start2<-as.numeric(temp)
        end2<-as.numeric(temp)+attr(temp,"match.length")-1
        if (min(start2,end2)>0) {
          length1<-nchar(seq1)
          length2<-nchar(seq2)
          if (start1<start2) { #take seq2 first then take seq1
            endseq<-end1+1
            if (endseq>length1) endseq<-length1
            seq1<-paste(substr(seq2,1,end2),substr(seq1,endseq,length1),sep="")
            seq2<-NA
          } else { #take seq1 first then take seq2
            endseq<-end2+1
            if (endseq>length2) endseq<-length2
            seq1<-paste(substr(seq1,1,end1),substr(seq2,endseq,length2),sep="")
            seq2<-NA
          }
          #cat("Merged Seq:",seq1,"\n")
        } else {
          #cat("No merge\n")
        }
      } 
    } 
    return(list(seq1=seq1,seq2=seq2))
  })
  mergeseq2<-unlist(mergeseq2)
  mergeseq2<-mergeseq2[!is.na(mergeseq2)]
}
if (length(hit1)>0 & length(hit2)>0) mergeseq<-c(mergeseq1,mergeseq2) 
else if (length(hit1)>0 & length(hit2)==0) mergeseq<-mergeseq1 
else if (length(hit1)==0 & length(hit2)>0) mergeseq<-mergeseq2
else if (length(hit1)==0 & length(hit2)==0) mergeseq<-c()
return(mergeseq)
} else { #if you employ SpliceConstructSearchWrapperV2 and indicate -directReads swithch
  fileF<-paste(barcode,"_",IN$fastqF,sep="")
  fileR<-paste(barcode,"_",IN$fastqR,sep="")
  cat("Loading Foward reads from:", fileF,"Reverse reads from:",fileR,"\n")
  seq1<-readLines(con=fileF)
  seq2<-readLines(con=fileR)
  cat("#Lines in F:",length(seq1),"#Lines in R:",length(seq2),"\n")
  cat("Merging Sequences...\n")
  #cu<-0
  mergeseq<-apply(cbind(seq1,seq2),1,function(x) {
    #cat("cu:",cu,"\t")
    seq1<-x[1]
    if (x[2]!="NA") {# if paired seq exist ,then merge
      seq2<-reversenuc(x[2])
      temp<-regexpr(barcode,seq1) #instead of variant position, actual 5ss gap taken
      start1<-as.numeric(temp)
      end1<-as.numeric(temp)+attr(temp,"match.length")-1
      local<-pairwiseAlignment(DNAString(seq1), DNAString(seq2), type = "local", substitutionMatrix = mat,gapOpening = -5, gapExtension = -2)
      common<-as.character(attr(local,"pattern"))
      #cat("R1 seq:",seq1," Paird R2 Seq:",seq2," Common Seq:",common,"\n")
      if (nchar(common)>=IN$mergethreshold) { #if common string is longer than indicated
        temp<-regexpr(common,seq1) 
        start1<-as.numeric(temp)
        end1<-as.numeric(temp)+attr(temp,"match.length")-1
        temp<-regexpr(common,seq2) 
        start2<-as.numeric(temp)
        end2<-as.numeric(temp)+attr(temp,"match.length")-1
        if (min(start2,end2)>0) {  #if actually found common seq in seq2     i.e. negative value means not found in the seq
          #cat("success\n")
          length1<-nchar(seq1)
          length2<-nchar(seq2)
          if (start1<start2) { #take seq2 first then take seq1
            endseq<-end1+1
            if (endseq>length1) endseq<-length1
            seq1<-paste(substr(seq2,1,end2),substr(seq1,endseq,length1),sep="")
            seq2<-NA
          } else { #take seq1 first then take seq2
            endseq<-end2+1
            if (endseq>length2) endseq<-length2
            seq1<-paste(substr(seq1,1,end1),substr(seq2,endseq,length2),sep="")
            seq2<-NA
          }
          #cat("Merged Seq:",seq1,"\n")
        } else {
          #cat("fail\n")
          #seq1<-c(seq1,seq2) #cat("No merge\n") #in case of no merge, reply seq1 and seq2 seprately
        }
      } else {
        #cat("fail\n")
        #seq1<-c(seq1,seq2) #cat("No merge\n") #in case of no merge, reply seq1 and seq2 seprately
      } 
    } 
    #cat("mergedSeq:",seq1,"\n")
    return(list(seq1=seq1,seq2=seq2))
  })
  mergeseq<-unlist(mergeseq)
  mergeseq<-mergeseq[!is.na(mergeseq)]
  #writeLines(mergeseq,con="temp.txt",sep="\n")
  return(mergeseq)
}
}



findgap<-function(refseq,seq) {
  refseq<-gsub(" ","",refseq)
  seq<-gsub(" ","",seq)
  #cat("Seq:",seq,"\n")
  globalAlign <-pairwiseAlignment(DNAString(refseq), DNAString(seq), substitutionMatrix = mat,gapOpening = -5, gapExtension = -2)
  #cat("GLOBAL ALIGN\n")
  if (IN$displayalign) print(globalAlign)
  pattern<-as.character(attr(globalAlign,"pattern"))
  patternstart<-as.numeric(attr(attr(attr(globalAlign,"pattern"),"range"),"start"))-1
  subject<-as.character(attr(globalAlign,"subject"))
  result<-list(gap=NA,seq=subject)
  block<-data.frame()
  cont<-0
  start<-0
  discon<-0
  for (i in 1:nchar(subject)) {
    if (substr(subject,i,i)!="-") {
      discon<-0
      if (cont==0) start<-i 
      cont<-cont+1
      if (i==nchar(subject) & cont>=IN$blocklength & discon<2) {
        new<-data.frame(start=start+patternstart,end=i-1+patternstart,length=i-start)
        block<-rbind(block,new)
      }
    } else {
      if (cont>=IN$blocklength & discon==1) {
        new<-data.frame(start=start+patternstart,end=i-2+patternstart,length=i-1-start)
        block<-rbind(block,new)
      } 
      discon<-discon+1
      if (discon>1) {cont<-0;start<-0}
    }
    #cat(i,"cont:",cont, "start:",start,"discon:",discon,"\n")
  } 
  if (nrow(block)==0) return(result)
  passblock<-data.frame()
  for (i in 1:nrow(block)) {
    ss1<-substr(refseq,block[i,]$start,block[i,]$end)
    ss2<-substr(subject,block[i,]$start-patternstart,block[i,]$end-patternstart)
    s1<-DNAString(ss1)
    s2<-DNAString(gsub("-","",ss2))
    max<-1*(block[i,]$length+1)*IN$blockaccuracy-3*(1-IN$blockaccuracy)
    if (IN$displayalign) cat("Block Check  RefBlock:",ss1, " CompareBlock:",ss2,"\nThreshold:",max)
    globalAlign <-pairwiseAlignment(s1,s2, substitutionMatrix = mat,gapOpening = 0, gapExtension = -3)
    if (as.numeric(attr(globalAlign,"score")) > max) passblock<-rbind(passblock,data.frame(start=block[i,]$start,end=block[i,]$end,length=block[i,]$length,id=i))
    if (IN$displayalign) cat("  actual score:",as.numeric(attr(globalAlign,"score")),"\n")
  }
  if (nrow(passblock)<2) return(result)
  passblock<-passblock[order(passblock$length,decreasing=T),]
  if (abs(passblock[1,]$id-passblock[2,]$id) > 1 ) return(result)  #2 blocks should be adjecent to each other
  if (passblock[1,]$end < passblock[2,]$start) {
    result<-list(gap=paste((passblock[1,]$end+1),"_",(passblock[2,]$start-1),sep=""),seq=subject)
  } else {
    result<-list(gap=paste((passblock[2,]$end+1),"_",(passblock[1,]$start-1),sep=""),seq=subject)
  }
  return(result)
}


#main 
#Input parametes -mutfile (chrpos format mut file with transcript ID) -output 

inputlist<-list(helpflag="-help",debug="-debug",inputframe="-inputframe",fastqF="-fastqF",fastqR="-fastqR",output="-output",
    mergethreshold="-mergethreshold",centerpenalty="-centerpenalty",takethreshold="-takethreshold",displayalign="-displayalign",
    nogapbonus="-nogapbonus",permitgapsize="-permitgapsize",reversecomp="-reversecomp",valiablethreshold="-valiablethreshold",
    unknowncollect="-unknowncollect",anonymousspliceassign="-anonymousspliceassign",blocklength="-blocklength",blockaccuracy="-blockaccuracy",
    unknownanonymousspliceassign="-unknownanonymousspliceassign",inputreadsanonymousspliceassign="-inputreadsanonymousspliceassign",
    directReads="-directReads")
inputparas<-c(0,0,1,1,1,1,
              1,0,1,0,
              0,1,0,1,
              0,0,1,1,
              0,1,0)  #the length of parameters to be input in the command line 0-switch yes or no 1- one variable 2-two variables
                        
                        

IN<-inputparameters(inputlist,inputparas)
if (IN$helpflag==1) stop("\n",Helpdocument)
if (IN$output==-9) IN$output<-"SpliceConstructSearch.v0.5.R.result.txt"
if (IN$mergethreshold==-9) IN$mergethreshold<-10
if (IN$takethreshold==-9) IN$takethreshold<--100
if (IN$permitgapsize==-9) IN$permitgapsize<-3
if (IN$blocklength==-9) IN$blocklength<-10
if (IN$blockaccuracy==-9) IN$blockaccuracy<-0.95
if (IN$directReads!=-9) {
    cat("Read assigned reads from each file corresponding to the barcode\n")
    IN$reversecomp<--9  #because readF file includes forward strand file and readR file includes reverse strand filem reversecomp switch not required.
}

cat("Parameters:\nNegative Nine is null value,\n")
for (i in names(IN)) {
  cat(i,":",unlist(IN[which(names(IN)==i)]),"\n")
}
if (IN$displayalign==-9) IN$displayalign<-0
if (IN$centerpenalty==-9) IN$centerpenalty<-0
if (IN$nogapbonus==-9) IN$nogapbonus<-0
if (IN$valiablethreshold==-9) IN$valiablethreshold<-0

IN$takethreshold<-as.numeric(IN$takethreshold)
IN$permitgapsize<-as.numeric(IN$permitgapsize)
IN$valiablethreshold<-as.numeric(IN$valiablethreshold)
IN$blocklength<-as.numeric(IN$blocklength)
IN$blockaccuracy<-as.numeric(IN$blockaccuracy)

if (IN$inputreadsanonymousspliceassign==-9) {
  cat("Reading input frame from:",IN$inputframe,"\n")
  inputframe<-read.table(IN$inputframe,sep="\t",header=T,stringsAsFactors = F)
  if ( IN$directReads==-9) {
    cat("Reading fastq R1 file from:",IN$fastqF,"\n")
    fastqF<-readLines(con=IN$fastqF)
    cat("Reading fastq R2 file from:",IN$fastqR,"\n")
    fastqR<-readLines(con=IN$fastqR)
  }
}
resultframe<-data.frame()
anonymousframe<-data.frame()

library(Biostrings)
mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)


#MAIN loop
if (IN$inputreadsanonymousspliceassign==-9) {
for (i in 1:nrow(inputframe)) {
  #search barcode and read strings
  barcode<-inputframe[i,1]
  mergeseq<-searchpairedreads(i,barcode)
  if (IN$reversecomp==1) {
    cat("Searching Reverse Complement Sequence of Barcode...\n")
    mergeseqrev<-searchpairedreads(i,reversenuc(barcode))
    if (length(mergeseqrev)>0) {
     mergeseqrev<-sapply(mergeseqrev,reversenuc)
     mergeseq<-c(mergeseq,mergeseqrev)
    }  
  }
  anoref<-inputframe[i,2] #in anonymous search mode, just refseq in 2nd column matters 
  new<-inputframe[i,]
  if (length(mergeseq)>0) {
    if(IN$debug==1) writeLines(mergeseq,con=paste("debug.each.renew.mergeseq.",barcode,".",IN$output,sep=""),sep="\n")
    
    #-anonymousspliceassign indicates that
    #-blocklength (integer) indicates the length of 1 block in -anonymousspliceassign switch
    #-unknownanonymousspliceassign indicates that the script conduct -anonymousspliceassign for unknown-assigned reads
    #-inputreadsanonymousspliceassign (filename) indicates that the script conduct -anonymousspliceassign for read in (filename)
    
    
    if (IN$anonymousspliceassign==-9) {#matched assign mode
    cat("Searching most matched refseq...\n")
    #blast ref seq vs mergeseq
    match<-sapply(mergeseq,function(seq){
      
      reflen<-sapply(2:ncol(inputframe),function(j){ 
        return(nchar(gsub(" ","",inputframe[i,j])))
      })
      minlen<-min(reflen)
      adjustscore<-sapply(2:ncol(inputframe),function(j) {
        return((reflen[j-1]-minlen)*2) # 2 is equal to gap extension penalty
      })
      score<-c()
      if (IN$valiablethreshold) {
        IN$takethreshold<-1*min(minlen,nchar(gsub(" ","",seq)))*IN$valiablethreshold-3*min(minlen,nchar(gsub(" ","",seq)))*(1-IN$valiablethreshold)-10*2-2*abs(minlen-nchar(gsub(" ","",seq)))
        if (IN$displayalign) cat("Variable threshold appiled:",IN$takethreshold,"\n")
      }
      for (j in 2:ncol(inputframe)) {
        refseq<-DNAString(gsub(" ","",inputframe[i,j]))
        seq<-DNAString(gsub(" ","",seq))
        globalAlign <-pairwiseAlignment(refseq, seq, substitutionMatrix = mat,gapOpening = -10, gapExtension = -2)
        addscore<-as.numeric(attr(globalAlign,"score"))
        addscore<-addscore+adjustscore[j-1]
        if (IN$displayalign) {
          cat("Align with Refseq",j-1,"\n")
          print (globalAlign)
          cat("Adjust score by length:",addscore,"\n")
        }
        if (IN$centerpenalty) {
          subject<-as.character(attr(globalAlign,"subject"))
          penalty<-penaltycalc(subject)
          if (IN$displayalign) cat("Center gap penalty: -",penalty,"\n",sep="")
          addscore<-addscore-penalty
          if (IN$nogapbonus) {
            addscore<-addscore-addscore*(addscore<0)*(penalty==0)
            if (penalty==0) cat("Nogap bonus added\n")
          }
        }
        if (IN$displayalign) cat("Final adjust score:",addscore,"\n")
        score<-c(score,addscore)
      }
      pos<-which(score==max(score))
      if (max(score)<=IN$takethreshold | length(pos)>1) pos<-ncol(inputframe)
      return(pos)
    })
    match<-as.numeric(unlist(match))
    if (IN$displayalign) cat("match:",match,"\n")
    if (IN$unknowncollect==1 | IN$unknownanonymousspliceassign==1) {
      unknownseq<-mergeseq[which(match==ncol(inputframe))]
      if (IN$unknowncollect==1) writeLines(unknownseq,con=paste(IN$output,".",barcode,"unknownseq.txt",sep=""),sep="\n")
      if (IN$unknownanonymousspliceassign==1) {
        cat("Searching anonymous gap of unknown reads to refseq...\n")
        gap<-sapply(mergeseq,function(seq){
          return(findgap(anoref,seq))
        })
        gap<-as.data.frame(t(gap))  #transform to data.frame  gap(gap=gap (ex 1_2) seq=aligned seq)
        gap<-data.frame(barcode=barcode,gap=unlist(gap$gap),seq=unlist(gap$seq))
        anonymousframe<-rbind(anonymousframe,data.frame(barcode=barcode,gap="Reference",seq=anoref),gap)
      }
    }
    
    freqtable<-data.frame(table(match))
    if (IN$debug==1) write.table(match,file="debug.match.table.txt",quote=F,sep="\t",col.names=T,row.names=F)
    for (j in 2:(ncol(inputframe)+1)) {
      id<-j-1
      pos<-which(freqtable$match==id)
      if (length(pos)>0) freq<-freqtable[pos,2] else freq<-0
      new<-cbind(new,data.frame(freq=freq))
      if (j<=ncol(inputframe)) names(new)[ncol(inputframe)+id]<-paste("Freq_",names(inputframe)[j],sep="")
      else names(new)[ncol(inputframe)+id]<-"Freq_Unknown"
    }
    } else { #anonymous assign mode
      cat("Searching anonymous gap to refseq...\n")
      gap<-sapply(mergeseq,function(seq){
        return(findgap(anoref,seq))
      })
      gap<-as.data.frame(t(gap))  #transform to data.frame  gap(gap=gap (ex 1_2) seq=aligned seq)
      gap<-data.frame(barcode=barcode,gap=unlist(gap$gap),seq=unlist(gap$seq))
      anonymousframe<-rbind(anonymousframe,data.frame(barcode=barcode,gap="Reference",seq=anoref),gap)
    }
  } else {
    if (IN$anonymousspliceassign==-9) {#matched assign mode
    for (j in 2:(ncol(inputframe)+1)) {
      id<-j-1
      new<-cbind(new,data.frame(freq=0))
      if (j<=ncol(inputframe)) names(new)[ncol(inputframe)+id]<-paste("Freq_",names(inputframe)[j],sep="")
      else names(new)[ncol(inputframe)+id]<-"Freq_Unknown"
    }
    } else { #anonymous assign mode
      gap<-data.frame(barcode=barcode,gap=NA,seq=NA)
      anonymousframe<-rbind(anonymousframe,data.frame(barcode=barcode,gap="Reference",seq=anoref),gap)
    }
  }
  
  resultframe<-rbind(resultframe,new)
  if (IN$debug==1) write.table(resultframe,file=paste("debug.each.renew.result.",IN$output,sep=""),sep="\t",col.names=T,row.names=F,quote=F)
}
} else { #inputreadsanonymousspliceassign mode
  seq<-readLines(con=IN$inputreadsanonymousspliceassign)
  cat("Searching anonymous gap from given reads to refseq...\n")
  anoref<-seq[1]
  mergeseq<-seq[-1]
  gap<-sapply(mergeseq,function(seq){
    return(findgap(anoref,seq))
  })
  gap<-as.data.frame(t(gap))  #transform to data.frame  gap(gap=gap (ex 1_2) seq=aligned seq)
  gap<-data.frame(barcode=IN$inputreadsanonymousspliceassign,gap=unlist(gap$gap),seq=unlist(gap$seq))
  anonymousframe<-rbind(anonymousframe,data.frame(barcode=IN$inputreadsanonymousspliceassign,gap="Reference",seq=anoref),gap)
  if (IN$debug==1) write.table(anonymousframe,file=paste("rawdata_anonymousalign_",IN$output,sep=""),sep="\t",col.names=T,row.names=F,quote=F)
}

#result

if (nrow(resultframe)>0) {
  cat("Writing result to:",IN$output,"\n")
  write.table(resultframe,file=IN$output,sep="\t",col.names=T,row.names=F,quote=F)
}
if (nrow(anonymousframe)>0) {
  #summarize the result
  barcodes<-unique(as.character(anonymousframe$barcode))
  if (IN$debug==1) cat("Barcodes:",barcodes,"\n")
  sumanonymous<-data.frame()
  for (barcode in barcodes) {
    seqs<-anonymousframe[anonymousframe$barcode==barcode,]
    refseq<-as.character(seqs[which(as.character(seqs$gap)=="Reference"),]$seq)[1]
    if (length(which(as.character(seqs$gap)=="Reference")>0)) seqs<-seqs[-which(as.character(seqs$gap)=="Reference"),]
    nafreq<-nrow(seqs[seqs$gap==NA,])
    list<-unlist(table(seqs$gap))
    if (length(list[list==0])>0) list<-list[-which(list==0)]
    list<-sort(list,decreasing=T)
    new<-data.frame(barcode=barcode,gap="Reference",seq=refseq,freq=0)
    sumanonymous<-rbind(sumanonymous,new)
    if (length(list)>0) {
      for (j in 1:length(list)) {
        gap<-names(list)[j]
        ngap<-as.numeric(unlist(strsplit(gap,"_")))
        
        new<-data.frame(barcode=barcode,gap=gap,freq=as.numeric(list[j]),seq=paste(substr(refseq,1,ngap[1]-1),paste(rep("-",ngap[2]-ngap[1]+1),collapse=""),substr(refseq,ngap[2]+1,nchar(refseq)),sep=""))
        sumanonymous<-rbind(sumanonymous,new)
      }
    }
    new<-data.frame(barcode=barcode,gap=NA,seq=NA,freq=nafreq)
    sumanonymous<-rbind(sumanonymous,new)
  }
  cat("Writing anonymous alignment summary to: anonymousalign_",IN$output,"\n",sep="")
  write.table(sumanonymous,file=paste("anonymousalign_",IN$output,sep=""),sep="\t",col.names=T,row.names=F,quote=F)
  cat("Writing anonymous alignment raw data to: rawdata_anonymousalign_",IN$output,"\n",sep="")
  write.table(anonymousframe,file=paste("rawdata_anonymousalign_",IN$output,sep=""),sep="\t",col.names=T,row.names=F,quote=F)
}
warnings()
