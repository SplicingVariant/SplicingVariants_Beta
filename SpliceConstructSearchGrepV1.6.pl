#SpliceConstructSearchGrepV1.5.pl
#this script is to find consensus seq using grep command and count the lines.
#
#Caution!!!!  Before run, type [module load  seq/flash/1.2.11] and load flash
#Caution!!!!   Beofre run, load R module if you are trying to run this script on orchestra server.
#
#Usage
#In v0.5,command -inputframe () -output () -basefolder ()
#In v0.6 or 0.7, command -R1 (fastqR1) -R2 (fastqR2) -inputframe() -flashout (header filename of FLASH output) -output () -basefolder()
#ex)perl SpliceConstructSearchGrepV0.5.pl -inputframe N500-1stbatch_for_Align_revised+Grep.txt -output N500-1stbatcg_for_Align_revised+Grep.v0.5.result.txt -basefolder /groups/seidman/solexa/bpf/150514_M00451_0056_AF41Y_FC01501/Kaoru/fastq
#ex)perl SpliceConstructSearchGrepV0.6B.pl -inputframe N500-1stbatch_for_Align_revised+Grep.txt -output temp.txt -R1 LIB015021_GEN00036527_GATCAG_L001_R1.fastq -R2 LIB015021_GEN00036527_GATCAG_L001_R2.fastq -flashout temp.flash.txt
#ex)perl SpliceConstructSearchGrepV0.6B.pl -inputframe N500-1stbatch_for_Align_revised+Grep.txt -output temp.txt -R1 LIB015021_GEN00036527_GATCAG_L001_R1.fastq -R2 LIB015021_GEN00036527_GATCAG_L001_R2.fastq -flashin temp.flash.txt   ***you can skip merging procedure if you have a flash file
#ex)perl SpliceConstructSearchGrepV0.7.pl -inputframe N500-1stbatch_for_Align_revised+Grep.txt -output temp.txt -R1 /n/data1/hms/genetics/seidman/solexa/bpf/150514_M00451_0056_AF41Y_FC01501/Kaoru/fastq/LIB015515_GEN00038389_CTTGTA_L001_R1.fastq -R2 /n/data1/hms/genetics/seidman/solexa/bpf/150514_M00451_0056_AF41Y_FC01501/Kaoru/fastq/LIB015515_GEN00038389_CTTGTA_L001_R2.fastq -flashout temp.flash.txt
#ex)perl SpliceConstructSearchGrepV0.8.pl -inputframe N500-1stbatch_for_Align_revised+Grep.txt -output temp.v0.8.txt -R1 /n/data1/hms/genetics/seidman/solexa/bpf/150514_M00451_0056_AF41Y_FC01501/Kaoru/fastq/LIB015515_GEN00038389_CTTGTA_L001_R1.fastq -R2 /n/data1/hms/genetics/seidman/solexa/bpf/150514_M00451_0056_AF41Y_FC01501/Kaoru/fastq/LIB015515_GEN00038389_CTTGTA_L001_R2.fastq -flashin temp.flash.txt -contamicheck 4,5
#ex)perl SpliceConstructSearchGrepV0.9.pl -inputframe N500-1stbatch_for_Align_revised+Grep.txt -output temp.v0.9.txt -R1 /n/data1/hms/genetics/seidman/solexa/bpf/150514_M00451_0056_AF41Y_FC01501/Kaoru/fastq/LIB015515_GEN00038389_CTTGTA_L001_R1.fastq -R2 /n/data1/hms/genetics/seidman/solexa/bpf/150514_M00451_0056_AF41Y_FC01501/Kaoru/fastq/LIB015515_GEN00038389_CTTGTA_L001_R2.fastq -flashout temp.flash.v0.9.txt -contamicheck 4,5
#ex)perl ../../SpliceConstructSearchGrepV1.0.pl -inputframe ../dm09_MYBPC3_inputframe.txt -output ../dm09_MYBPC3_inputframe.txt.Align_Grep.v1.0.txt -R1 150921_FC01738_L1_4_DM09--MYBPC3-1--Hs--R1.fastq -R2 150921_FC01738_L1_4_DM09--MYBPC3-1--Hs--R2.fastq -flashout temp.150921_FC01738_L1_4_DM09--MYBPC3-1--Hs.flash.v1.0.txt -contamicheck 4
#ex)perl ../../SpliceConstructSearchGrepV1.1.pl -inputframe ../dm09_MYBPC3_inputframe.txt -output ../dm09_MYBPC3_inputframe.txt.Align_Grep.v1.1.txt -R1 150921_FC01738_L1_4_DM09--MYBPC3-1--Hs--R1.fastq -R2 150921_FC01738_L1_4_DM09--MYBPC3-1--Hs--R2.fastq -flashin temp.150921_FC01738_L1_4_DM09--MYBPC3-1--Hs.flash.v1.0.txt -contamicheck 4 -minseq ../LMNA_Var_MYBPC3.result.v0.81.final.txt,3,8,8
#ex)perl ../../SpliceConstructSearchGrepV1.11.pl -inputframe ../dm09_MYBPC3_inputframe.txt -output ../dm09_MYBPC3_inputframe.txt.Align_Grep.v1.11.txt -R1 150921_FC01738_L1_4_DM09--MYBPC3-1--Hs--R1.fastq -R2 150921_FC01738_L1_4_DM09--MYBPC3-1--Hs--R2.fastq -flashin temp.150921_FC01738_L1_4_DM09--MYBPC3-1--Hs.flash.v1.0.txt -contamicheck 4 -minseq ../LMNA_Var_MYBPC3.result.v0.81.final.txt,1,8,8,8
#ex)perl ../SpliceConstructSearchGrepV1.2.pl -addstats MYBPC3_inputframe.txt.Align_Grep.v1.1+.txt
#ex) perl ../../SpliceConstructSearchGrepV1.3.pl -inputframe dm10_LMNA_MYBPC3.ss3.1sttry_inputframe.txt -output dm10_LMNA_MYBPC3.ss3.1sttry_inputframe.Align_Grep.v1.3.txt -R1 160111-FC01937-L1-DM10--Kaoru--LMNA-MYBPC3-3SS-2--Hs--trimmed-pair1.fastq -R2 160111-FC01937-L1-DM10--Kaoru--LMNA-MYBPC3-3SS-2--Hs--trimmed-pair2.fastq -flashin temp.160111-FC01937-L1-DM10--Kaoru--LMNA-MYBPC3-3SS-2--Hs--trimmed.flash.v1.0.txt -contamicheck 4
#ex) perl ../../SpliceConstructSearchGrepV1.3.pl -inputframe dm09_LMNA_MYBPC3.ss3.1sttry_inputframe.txt -output dm09_LMNA_MYBPC3.ss3.1sttry_inputframe.Align_Grep.v1.3minseq.txt -R1 160111-FC01937-L1-DM09--Kaoru--LMNA-MYBPC3-3SS-1--Hs--trimmed-pair1.fastq -R2 160111-FC01937-L1-DM09--Kaoru--LMNA-MYBPC3-3SS-1--Hs--trimmed-pair2.fastq -flashin temp.160111-FC01937-L1-DM09--Kaoru--LMNA-MYBPC3-3SS-1--Hs--trimmed.flash.v1.0.txt -contamicheck 4 -minseq ../LMNA_MYBPC3.ss3.result+.v0.91.1sttry.txt,1,8,8,8
#ex) perl ../../SpliceConstructSearchGrepV1.5.pl -inputframe dm09_LMNA_MYBPC3.ss3.1sttry_inputframe.txt -output dm09_LMNA_MYBPC3.ss3.1sttry_inputframe.Align_Grep.v1.3minseq.txt -R1 160111-FC01937-L1-DM09--Kaoru--LMNA-MYBPC3-3SS-1--Hs--trimmed-pair1.fastq -R2 160111-FC01937-L1-DM09--Kaoru--LMNA-MYBPC3-3SS-1--Hs--trimmed-pair2.fastq -flashin temp.160111-FC01937-L1-DM09--Kaoru--LMNA-MYBPC3-3SS-1--Hs--trimmed.flash.v1.0.txt -contamicheck 4 -addstats -countreads -getannotation LMNA_MYBPC3.ss3.result+.v0.91.1sttry.txt
#ex)/Users/meizhu/Copy/Minigeneassay/no-vector/15-LMNA-MYBPC3-3SS perl ../SpliceConstructSearchGrepV1.5.pl -output LMNA-MYBPC3-3SS-1sttry_inputframe.Align_Grep.v1.3.stats.txt -getannotation LMNA_MYBPC3.ss3.result+.v0.91.1sttry.txt
#Switches
#-inputframe (filename) 
# the file format should be as follows (In both v0.5 and v0.6 the same)
#Construct_Name	ID(SeqToIdentifyTheConstruct)	File1	File2	SeqToIdentifyTheMUTANT	SeqToGrep1 SeqToGre2
#File1 -> fastqR1   File2-> fastqR2  RefSeqs_Around_MutPoint-> seqs you want to grep ....
# * this script greps given seq and the reverse compliment
# * In v0.6 or later, the columns File1 and File2 will be ignored.
#-output (filename) contains the result
#-basefolder (folder) if the R1 and R2 files are not on the current folder, you can indicate the folder on which fastqR1 and R2 files are
#-flashout (filename) filename for R1 & R2 fastq merging result
#-flashin (filename) if you have made a flashout file and want to count it, use this switch (skip merging phase and save time)
#-contamicheck (columns of seqs for checking in the inputframe file, delimited by comma when you indicate several) if you want to check contamination ( which means involvement of other construct that you don not expect), please indicate # of columns of seqs where other construct seqs found. note that the first column is 0 as a perl way.
#-minseq (ConstructDesingerFileV0.8orlater [,#attach extra bp to trimmed seq, #maxlimit of cut for one side, #mer of 5 side seq,#mer of 3side seq]) indicates that grep file will be trimmed as long as it keeps the uniqueness. To get refseq, you need to indicate ConstructDesingerFile.
#																				Also you can indicate ( or omit) #mer of 5 side seq,#mer of 3side seq for the script to understand how to trim grepseq. The default setting is 4,4 (grep is 8bp long , which consists of 4bp 5 prime seq and 4bp 3 prime seq.
#-addstats  when you have a sheet generated by SpliceConstructSearchGrepV1.11.pl or before and want to add statistics, please indicate the file. Then <FILENAME>.stats.<EXT> will be made.
#															when you count the reads without addstats switch , addstats automatically conducted and fileoutput and fileoutput.stats.txt will be generated. 
#									**caution if you indicate just -addstats, this script does NOT count reads
#				you need to indicate the file name for stat calculation using -output switch (file should be made by SpliceConstructSearchGrepV1.11.pl or before)
#-countreads indicates to count the reads. if not, this script does NOT count the reads
#					In v1.3 you don not need to indicate this , because it always count the reads except -addstats
#-getannotation (SpliceConstructfile1,SpliceConstructfile2,... ) when you want to add information 
#						original file which is being added information should be indicated by -output switch
#						you can select several splice construct files by comma-delimited file name indication
#						**Caution new file name is the same as original output file name
#-getannotationlaxid    indicates when -getannoation conducted, construct identifier will be just Varname. ( if not, it wiil be VarNamr+ConstructID)
#history
#v0.5 launch ver.
#v0.6 using FLASH, merge FASTQ R1&R2 then grep
#v0.6B makes use of the another way to pick out seq lines
#v0.7 bug fix
#v0.8 comtamination check function added
#v0.9 make unmatched reads merging process faster, add % calc feature
#v1.0 for SS broken, <Anonymous>CTTGAGGCAG in AberrantSplicing column implemented.
#  if <Anonymous> found in the grep column , grep using strings except <Anonymous> then grep -v using the #-1, -2 columns, then count #lines.
#	So please place NO-splice and Normal-Splice stirngs before <Anonymous>Aberrant splicing 
#v1.1 for short read, this script automatically reduce grep seq as long as it maintains the uniqueness in the seq, in order to improve the signal.
#	please implement this function using the switch -minseq 
#v1.11  refine the -minseq, can attach extra mers to unique seq
#v1.2 renew statistics to improve the interpretation
#      add command to attach statistics that already calculated
#v1.3 showing grep status
#		bug-fix for trim seq
#     -addstats does NOIT work on Orchestra, but works on Mac because of the status of module installation. 
#v1.5 some modification for addstats and countreads
#	add -annotation swith
#v1.6 --getannotation problem
#     If you fix ID-duplicate when order constructs,  ID in construct make file and ID in gblockorderfile will differ from each other.
#	   Then --getannotation indicates construct make file, some constructs info will be NA because of the above reason. (because -getannotation employ VarName+ConstructID)
#     to above the above situation, -getannotationlaxid added. when you employ this, construct identifier will be not VarName+constructID but just Varname. 
#     -getannotation problem 2
#       old var construct designer does not change + or - -> _ , while newer ver changes + or - -> _
#       Therefore 2kinds of construct name can exist.
#       So I made -getannotation deal with both.
#future direction
#  parallel run   change temp file name
#  can indicate flash script location

print "V1.6.\n";
use strict;
use warnings;
use File::Basename;
use lib qw(/home/ki47/perl/libs/share/perl5);  #for orchestra., Stattistics::R modele is there; if you are not sure how to install the module, see my wiki.
use Statistics::R;


#getting parameters
my $inputframe;
my $basefolder="-9";
my $output="temp.spliceconstruct.search.output.txt";
my $R1;
my $R2;
my $flashout="temp.flash.txt";
my $flashin="-9";
my $contamicheck="-9";
my $minseq="-9";
my $addstats="-9";
my $countreads="-9";
my $getannotation="-9";
my $getannotationlaxid=-9;
for my $I (0..$#ARGV) {
 $inputframe=$ARGV[$I+1] if ($ARGV[$I] eq "-inputframe");
 $basefolder=$ARGV[$I+1] if ($ARGV[$I] eq "-basefolder");
 $output=$ARGV[$I+1] if ($ARGV[$I] eq "-output");
 $R1=$ARGV[$I+1] if ($ARGV[$I] eq "-R1");
 $R2=$ARGV[$I+1] if ($ARGV[$I] eq "-R2");
 $flashout=$ARGV[$I+1] if ($ARGV[$I] eq "-flashout");
 $flashin=$ARGV[$I+1] if ($ARGV[$I] eq "-flashin");
 $contamicheck=$ARGV[$I+1] if ($ARGV[$I] eq "-contamicheck");
 $minseq=$ARGV[$I+1] if ($ARGV[$I] eq "-minseq");
 $addstats=$output if ($ARGV[$I] eq "-addstats");
 $countreads=$ARGV[$I+1] if ($ARGV[$I] eq "-countreads");
 $getannotation=$ARGV[$I+1] if ($ARGV[$I] eq "-getannotation");
 $getannotationlaxid=1 if ($ARGV[$I] eq "-getannotationlaxid");
}

my $nleft=8; #mer of 3 side
my $nright=8;  #mer of 5 side
my $maxcut=8;
my $attachextra=0;
my %refseq;  #refseq for searching in minseq function  key is varname ; content is corresponding seq
my %spseq1; 
my %spseq2;
my @othergreps;

if ($countreads ne "-9") {
if ($minseq ne "-9") {
	my $pos_Var_Name=10000;
	my $pos_Ref_Construct=10000;
	my $pos_Alt_Construct=10000;
	my $pos_Ref_Norm_Splice=10000;
	my $pos_Ref_Aberrant_Splice=10000;
	my $pos_Mut_Norm_Splice=10000;
	my $pos_Mut_Aberrant_Splice=10000;
	
	my $constfile;  #filename of Construct Designer File from which refseq is being retreaved
	my @SPL=split(/,/,$minseq);
	$constfile=$SPL[0];
	if ($#SPL>1) {
		$attachextra=$SPL[1];$maxcut=$SPL[2];$nleft=$SPL[3];$nright=$SPL[4];
	}
	print "Minimizing Grep Switch On.\nConstruct Designer File:$constfile\t#max cut:$maxcut\t#mer of 5\'side of Seq:$nleft\t#mer of 3\'side of Seq:$nleft\nReading $constfile and checking entire sequence of each construct...\n";
	#extraacting ref seq for each construct
	open (FILE, $constfile) or die "$!";
	my $head=<FILE>;
	chomp($head);
	$head=~s/\r+//g;
	my @head=split(/\t/,$head);
	for my $I (0..$#head) {
		$pos_Var_Name=$I if ($head[$I] eq "Var_Name");
		$pos_Ref_Construct=$I if ($head[$I] eq "Ref_Construct");
		$pos_Alt_Construct=$I if ($head[$I] eq "Alt_Construct");
		$pos_Ref_Norm_Splice=$I if ($head[$I] eq "Ref_Norm_Splice");
		$pos_Ref_Aberrant_Splice=$I if ($head[$I] eq "Ref_Aberrant_Splice");
		$pos_Mut_Norm_Splice=$I if ($head[$I] eq "Mut_Norm_Splice");
		$pos_Mut_Aberrant_Splice=$I if ($head[$I] eq "Mut_Aberrant_Splice");
	}
	die "Header info in the construct designer file is incorrect\n" if (($pos_Var_Name+$pos_Ref_Construct+$pos_Alt_Construct+$pos_Ref_Norm_Splice+$pos_Ref_Aberrant_Splice+$pos_Mut_Norm_Splice+$pos_Mut_Aberrant_Splice)>9999);
	while (my $I=<FILE>) {
		chomp($I);
		$I=~s/\r+//g;
		my @SPL2=split(/\t/,$I);
		my $newname=$SPL2[$pos_Var_Name];
		$newname=~s/(\.|:|>)/_/g;
		$refseq{"REF_".$newname}=$SPL2[$pos_Ref_Construct];
		$spseq1{"REF_".$newname}=&atgc($SPL2[$pos_Ref_Norm_Splice]);
		$spseq2{"REF_".$newname}=&atgc($SPL2[$pos_Ref_Aberrant_Splice]);
		$refseq{"ALT_".$newname}=$SPL2[$pos_Alt_Construct];
		$spseq1{"ALT_".$newname}=&atgc($SPL2[$pos_Mut_Norm_Splice]);
		$spseq2{"ALT_".$newname}=&atgc($SPL2[$pos_Mut_Aberrant_Splice]);
	}
}


if ($flashin eq "-9") {
	#merge R1 and R2 using flash
	print "Merging R1 and R2 using flash...\n";
	my $command="/n/data1/hms/genetics/seidman/ki47/flash/FLASH-1.2.11/flash $R1 $R2 2>&1 | tee flash.log";
	system($command);
#   - out.extendedFrags.fastq      The merged reads.
#   - out.notCombined_1.fastq      Read 1 of mate pairs that were not merged.
#   - out.notCombined_2.fastq      Read 2 of mate pairs that were not merged.
#   - out.hist                     Numeric histogram of merged read lengths.
#   - out.histogram                Visual histogram of merged read lengths.

	print "Combining unmatched R1 and R2 lines...\nPreprocessing...\n";
	system('grep -A1 "^@" out.notCombined_1.fastq | grep -v "^@" | grep -v "^--" >temp.spliceconstructsearch.file0_1.txt');
	system('grep -A1 "^@" out.notCombined_2.fastq | grep -v "^@" | grep -v "^--" >temp.spliceconstructsearch.file0_2.txt');
	open (FILE1, "./temp.spliceconstructsearch.file0_1.txt") or die "$!";
	open (FILE2, "./temp.spliceconstructsearch.file0_2.txt") or die "$!";
	open (FILEOUT, "> temp.spliceconstructsearch.file1.txt");
	print "Merging Start\n";
	my $counter=0;
	while (my $read1=<FILE1>) {
		my $read2=<FILE2>;
		#my $rev=&revcom($read2);
		chomp($read1);
		chomp($read2);
		#chomp($rev);
		print FILEOUT $read1."NNN".$read2."\n";
		$counter++;
		warn "$counter\t" if ($counter % 10000 == 0);
	}
	close(FILE1);close(FILE2);close(FILEOUT);

	system('grep -A1 "^@" out.extendedFrags.fastq | grep -v "^@" | grep -v "^--" >temp.spliceconstructsearch.file2.txt');
	system("cat temp.spliceconstructsearch.file1.txt temp.spliceconstructsearch.file2.txt > $flashout");
	system("rm temp.spliceconstructsearch.file?.txt");
	system("rm out.extendedFrags.fastq");
	system("rm out.notCombined_?.fastq");
	system("rm temp.spliceconstructsearch.file0_?.txt");
	system("rm out.hist");system("rm out.histogram");
} else {
	$flashout=$flashin;
}

print "\nCounting started...\n";
open (FILE, $inputframe) or die "$!";
open (FILEOUT, "> $output");
#head output
my @FILE=<FILE>;
close(FILE);
my $head=$FILE[0];
chomp($head);
$head=~s/\r+//g;
print FILEOUT $head;
my @hSPL=split(/\t/,$head);
my $element=$#hSPL-4;
my $r="#Reads_";
my @nCol;
@nCol=split(/,/,$contamicheck) if ($contamicheck ne "-9");
print FILEOUT "\t$r"."_Total";
for my $J (0..$element) {
	my $flag=1;
	if ($contamicheck ne "-9") {for my $K (@nCol) {$flag=0 if ($K eq ($J+4));}}  #skip contamicheck column because these will be checked later 
	print FILEOUT "\t$r"."_$hSPL[4+$J]" if ($flag);
}
if ($contamicheck ne "-9") {
	for my $cI (1..$#FILE) {
		my $I=$FILE[$cI];
		my @SPL=split(/\t/,$I);
		my $rr="$r$SPL[0]";
		for my $J (@nCol) {
			print FILEOUT "\t$rr"."_$hSPL[$J]";
		}
	}
}



print FILEOUT "\n";
#while loop
for my $cI (1..$#FILE) {
#while (my $I=<FILE>) {
	my $I=$FILE[$cI];
	chomp($I);
	$I=~s/\r+//g;
	print FILEOUT $I;
	my @SPL=split(/\t/,$I);
	print "Processing Sample: $SPL[0]\n";
	##Construct_Name	ID	File1	File2	RefSeqs_Around_MutPoint	MutSeqs_Around_MutPoint	SeqsWhenNormaltSplicingHappens	SeqsWhenAberrantSplicingHappens
		#my $grepfile=$SPL[2+$j];
		#$grepfile="$basefolder/".basename($grepfile) if ($basefolder ne "-9") ;
		#egrep "ATCGGCTAAGGAGCCCAGAG|ACATGGAGATCCACGCCTACCGCAAGCTCTTGGAGGGCG" out.notCombined_1.fastq  | head
	my $ID=$SPL[1];
	my $revID=&revcom($ID);
	my $command="egrep \"$ID|$revID\" $flashout > temp.spliceconstructsearch.file.$ID.txt";
	system($command);
	$command="wc temp.spliceconstructsearch.file.$ID.txt";
	my $stout=`$command`;my @sSPL=split(/\s+/,$stout);
	my $count=$sSPL[1];
	print FILEOUT "\t$count";
	@othergreps=();
	for my $k (0..$element) {push(@othergreps,uc($SPL[4+$k]));}
	for my $k (0..$element) {
			my $flag=1;
			if ($contamicheck ne "-9") {for my $l (@nCol) {$flag=0 if (($k+4)==$l);}}  #skip contamicheck column because these will be checked later 
			if ($flag) {
				my $voidgrep=0;
				my $void;
				my $voidrev;
				my $void2;
				my $voidrev2;
				my $myseq;
				if ($minseq ne "-9") {
					$myseq=&trimseq($SPL[4+$k],$SPL[0]) ;
				} else {
					$myseq=$SPL[4+$k];
				}
				if ($myseq=~/<Anonymous>/) {
					$myseq=~s/<Anonymous>//;
					$voidgrep=1;
					if ($minseq eq "-9") {$void=&atgc($SPL[3+$k]);} else {$void=&atgc(&trimseq($SPL[3+$k],$SPL[0]));}
					$voidrev=&revcom($void);
					if ($minseq eq "-9") {$void2=&atgc($SPL[2+$k]);} else {$void2=&atgc(&trimseq($SPL[2+$k],$SPL[0]));}
					$voidrev2=&revcom($void2);
				} 
				my $greped1=&atgc($myseq);
				my $greped2=&revcom($greped1);
				my @sSPL;
				if ($voidgrep) {
					print "Sequence Search with Anonymous - Grep: $greped1 - $greped2 Vgrep:$void - $voidrev, $void2 - $voidrev2\n";
					my $command="egrep \"$greped1|$greped2\" temp.spliceconstructsearch.file.$ID.txt | egrep -v \"$void|$voidrev|$void2|$voidrev2\" | wc";	my $stout=`$command`;@sSPL=split(/\s+/,$stout);
				} else {
					print "Sequence Search - Grep: $greped1 - $greped2\n";
					my $command="egrep \"$greped1|$greped2\" temp.spliceconstructsearch.file.$ID.txt | wc";	my $stout=`$command`;@sSPL=split(/\s+/,$stout);
				}
				my $fcount=$sSPL[1];
				my $temp=$fcount/$count*100;
				my $round=sprintf("%.1f",$temp);
				print FILEOUT "\t$round% ($fcount)";
			}
	}
	if ($contamicheck ne "-9") {
		for my $cK (1..$#FILE) {
			my $K=$FILE[$cK];
			my @SPL=split(/\t/,$K);
			for my $J (@nCol) {
				my $greped1=&atgc($SPL[$J]);
				my $greped2=&revcom($greped1);
				print "Contamination Check - Grep: $greped1 - $greped2\n";
				my $command="egrep \"$greped1|$greped2\" temp.spliceconstructsearch.file.$ID.txt | wc";	my $stout=`$command`;my @sSPL=split(/\s+/,$stout);
				my $fcount=$sSPL[1];
				my $temp=$fcount/$count*100;
				my $round=sprintf("%.1f",$temp);
				print FILEOUT "\t$round% ($fcount)";
				#print FILEOUT "\t$rr"."_$hSPL[$J]";
			}
		}
	}	
	print FILEOUT "\n";
}#for $cI end
#close(FILE);
close(FILEOUT);
print "Counting finished.\n";
}


#get annotation part
if ($getannotation ne "-9") {
	my %Var_Name;
	my %Kind;
	my %TranscriptID;
	my %Ref_MaxEntScore;
	my %Alt_MaxEntScore;
	my %Diff_MaxEntScore;
	my %First_Exon;
	
	my $pos_Var_Name=9999;
	my $pos_Kind=9999;
	my $pos_TranscriptID=9999;
	my $pos_Ref_MaxEntScore=9999;
	my $pos_Alt_MaxEntScore=9999;
	my $pos_Diff_MaxEntScore=9999;
	my $pos_First_Exon=9999;
	my $pos_refID=9999;
	my $pos_altID=9999;
	
	my @constfile=split(/,/,$getannotation); #as long as tab-delimited, pleural files can be indicated
	print "Construct Designer Files for Annotation: @constfile\n";
	for my $constfile (@constfile) {
		open (FILE, $constfile) or die "$!";
		my $head=<FILE>;
		chomp($head);
		$head=~s/\r+//g;
		my @head=split(/\t/,$head);
		for my $I (0..$#head) {
			$pos_Var_Name=$I if ($head[$I] eq "Var_Name");
			$pos_Kind=$I if ($head[$I] eq "Kind");
			$pos_TranscriptID=$I if ($head[$I] eq "TranscriptID");
			$pos_Ref_MaxEntScore=$I if ($head[$I] =~/Ref_MaxEntScore/);
			$pos_Alt_MaxEntScore=$I if ($head[$I] =~/Alt_MaxEntScore/);
			$pos_Diff_MaxEntScore=$I if ($head[$I] =~/Diff_MaxEntScore/);
			$pos_First_Exon=$I if ($head[$I] eq "First_Exon");
			$pos_refID=$I if ($head[$I] eq "Ref_ID");
			$pos_altID=$I if ($head[$I] eq "Alt_ID");
		}
		die "Header info in the construct designer file is incorrect\n" if (($pos_Var_Name+$pos_Kind+$pos_TranscriptID+$pos_Ref_MaxEntScore+$pos_Alt_MaxEntScore+$pos_Diff_MaxEntScore+$pos_First_Exon+$pos_refID+$pos_altID)>9999);
		while (my $I=<FILE>) {
			chomp($I);
			$I=~s/\r+//g;
			my @SPL2=split(/\t/,$I);
				my $newnameref;
				my $newnamealt;
			if ($getannotationlaxid==1) {
				$newnameref="$SPL2[$pos_Var_Name]";
				$newnamealt="$SPL2[$pos_Var_Name]";
			} else {
				$newnameref="$SPL2[$pos_Var_Name]$SPL2[$pos_refID]";
				$newnamealt="$SPL2[$pos_Var_Name]$SPL2[$pos_altID]";
			}
			$newnameref=~s/(\.|:|>)/_/g;  #new name will be ID for annotation
			$newnamealt=~s/(\.|:|>)/_/g; 
			$Var_Name{$newnameref}=$SPL2[$pos_Var_Name];$Var_Name{$newnamealt}=$SPL2[$pos_Var_Name];
			$Kind{$newnameref}=$SPL2[$pos_Kind];$Kind{$newnamealt}=$SPL2[$pos_Kind];
			$TranscriptID{$newnameref}=$SPL2[$pos_TranscriptID];$TranscriptID{$newnamealt}=$SPL2[$pos_TranscriptID];
			$Ref_MaxEntScore{$newnameref}=$SPL2[$pos_Ref_MaxEntScore];$Ref_MaxEntScore{$newnamealt}=$SPL2[$pos_Ref_MaxEntScore];
			$Alt_MaxEntScore{$newnameref}=$SPL2[$pos_Alt_MaxEntScore];$Alt_MaxEntScore{$newnamealt}=$SPL2[$pos_Alt_MaxEntScore];
			$Diff_MaxEntScore{$newnameref}=$SPL2[$pos_Diff_MaxEntScore];$Diff_MaxEntScore{$newnamealt}=$SPL2[$pos_Diff_MaxEntScore];
			$First_Exon{$newnameref}=$SPL2[$pos_First_Exon];$First_Exon{$newnamealt}=$SPL2[$pos_First_Exon];
			#if + or - involved in Varname
			$newnameref=~s/(\+|\-)/_/g;  #new name will be ID for annotation
			$newnamealt=~s/(\+|\-)/_/g; 
			$Var_Name{$newnameref}=$SPL2[$pos_Var_Name];$Var_Name{$newnamealt}=$SPL2[$pos_Var_Name];
			$Kind{$newnameref}=$SPL2[$pos_Kind];$Kind{$newnamealt}=$SPL2[$pos_Kind];
			$TranscriptID{$newnameref}=$SPL2[$pos_TranscriptID];$TranscriptID{$newnamealt}=$SPL2[$pos_TranscriptID];
			$Ref_MaxEntScore{$newnameref}=$SPL2[$pos_Ref_MaxEntScore];$Ref_MaxEntScore{$newnamealt}=$SPL2[$pos_Ref_MaxEntScore];
			$Alt_MaxEntScore{$newnameref}=$SPL2[$pos_Alt_MaxEntScore];$Alt_MaxEntScore{$newnamealt}=$SPL2[$pos_Alt_MaxEntScore];
			$Diff_MaxEntScore{$newnameref}=$SPL2[$pos_Diff_MaxEntScore];$Diff_MaxEntScore{$newnamealt}=$SPL2[$pos_Diff_MaxEntScore];
			$First_Exon{$newnameref}=$SPL2[$pos_First_Exon];$First_Exon{$newnamealt}=$SPL2[$pos_First_Exon];

		}
	}
	open (FILE, $output) or die "$!";
	my $newfile=$output;
	$newfile=~s/(\S+)\.(\S+)$/$1\.annot\.$2/; #rename for file with stats
	open (FILEOUT, "> $newfile");
	my $head=<FILE>;
	chomp($head);
	$head=~s/\r+//g;
	print FILEOUT "$head\tVar_Name\tREF_ALT\tKind\tTranscriptID\tFirst_Exon\tRef_MaxEntScore\tAlt_MaxEntScore\tDiff_MaxEntScore\n";
	my @head=split(/\t/,$head);
	my $p_name=9999;
	my $p_ID=9999;
	my $p_abesp=9999;
	for my $I (0..$#head) {
		$p_name=$I if ($head[$I] eq "Construct_Name");
		$p_ID=$I if ($head[$I] eq "ID");
		$p_abesp=$I if ($head[$I] eq "AberrantSplicing");
	}
	while (my $I=<FILE>) {
		chomp($I);
		$I=~s/\r+//g;
		print FILEOUT "$I\t";
		my @SPL=split(/\t/,$I);
		my $cname=$SPL[$p_name];
		$cname=~s/(REF_|ALT_)//g;
		if ($getannotationlaxid==-9) {
			$cname="$cname$SPL[$p_ID]";
		}
		my $Kind="SS_Create";
		my $REFALT="REF";
		$Kind="SS_Broken" if ($SPL[$p_abesp]=~/<Anonymous>/);
		$REFALT="ALT" if ($SPL[$p_name]=~/ALT_/);
		if (exists($Var_Name{$cname})) {
			print FILEOUT "$Var_Name{$cname}\t$REFALT\t$Kind\t$TranscriptID{$cname}\t$First_Exon{$cname}\t$Ref_MaxEntScore{$cname}\t$Alt_MaxEntScore{$cname}\t$Diff_MaxEntScore{$cname}\n";
		} else {
			print FILEOUT "NA\t$REFALT\t$Kind\tNA\tNA\tNA\tNA\tNA\n";
		}
	}
	close(FILEOUT);
	close(FILE);
	system("mv $newfile $output");
} #get annotation part end




#$addstats part
if ($addstats ne "-9") {
print "Calculating Statistics...\n";
open (FILE, $addstats) or die "$!";
my $newname=$addstats;
$newname=~s/(\S+)\.(\S+)$/$1\.stats\.$2/; #rename for file with stats
open (FILEOUT, "> $newname") or die "$!";
my $head=<FILE>;
chomp($head);
$head=~s/\r+//g;
my @head=split(/\t/,$head);
my $p_name=9999;
my $p_Total=9999;
my $p_NoSplicing=9999;
my $p_NormaltSplicing=9999;
my $p_Aberrant=9999;
my $p_AberrantSplicing=9999;
for my $I (0..$#head) {
	$p_name=$I if ($head[$I] eq "Construct_Name");
	$p_Total=$I if ($head[$I] eq "#Reads__Total");
	$p_NoSplicing=$I if ($head[$I] eq "#Reads__NoSplicing");
	$p_NormaltSplicing=$I if ($head[$I] eq "#Reads__NormaltSplicing");
	$p_Aberrant=$I if ($head[$I] eq "#Reads__AberrantSplicing");
	$p_AberrantSplicing=$I if ($head[$I] eq "AberrantSplicing");
}
#print "$p_Total\t$p_NoSplicing\t$p_NormaltSplicing\t$p_Aberrant\t$p_AberrantSplicing\n";
die "Irregular header information.\n" if (($p_Total+$p_NoSplicing+$p_NormaltSplicing+$p_Aberrant+$p_AberrantSplicing)>=9999);
print FILEOUT "$head\tKind\tN_AnalyzableReads\tN_NoSplice\tN_NormalSplice\tN_AberrantSplice\tPercent_NoSplice\tPercent_NormalSplice\tPercent_AberrantSplice\tRatio_Aberrant_Normal\tRatio_NOandAberrant_Normal\tp_value\tPathogenic_Score\n";

my $pre_per_nosplice;
my $pre_per_normalsplice;
my $pre_per_aberrantsplice;

while (my $I=<FILE>) {
	chomp($I);
	$I=~s/\r+//g;
	print FILEOUT "$I\t";
	my @SPL=split(/\t/,$I);
	my $kind="SS_Create";
	$kind="SS_Broken" if ($SPL[$p_AberrantSplicing]=~/<Anonymous>/);
	my $total=$SPL[$p_Total];
	my $nosplice=&extractnum($SPL[$p_NoSplicing]);
	my $normalsplice=&extractnum($SPL[$p_NormaltSplicing]);
	my $aberrantsplice=&extractnum($SPL[$p_Aberrant]);
	my $sum=$nosplice+$normalsplice+$aberrantsplice; #analyzable reads
	$sum=0.001 if ($sum==0); # to avoid division by zero
	my $per_nosplice=sprintf("%.1f", $nosplice/$sum*100);
	my $per_normalsplice=sprintf("%.1f", $normalsplice/$sum*100);
	my $per_aberrantsplice=sprintf("%.1f", $aberrantsplice/$sum*100);
	my $ratio;
	my $ratio2;
	if ($normalsplice>0) {
		$ratio=$aberrantsplice/$normalsplice;
		$ratio2=($aberrantsplice+$nosplice)/$normalsplice; 
	} else {$ratio="NA";$ratio2="NA";}
	my $p_val=1;
	my $pathogenic_score=0;
	if ($SPL[$p_name]=~/ALT_/) {
		$p_val= &fishertest(int($per_aberrantsplice),int($per_normalsplice),int($pre_per_aberrantsplice),int($pre_per_normalsplice)) if ($kind eq "SS_Create");
		$p_val= &fishertest(int($per_normalsplice),int($per_aberrantsplice+$per_nosplice),int($pre_per_normalsplice),int($pre_per_nosplice+$pre_per_aberrantsplice)) if ($kind eq "SS_Broken");
		if ((($per_aberrantsplice+$per_normalsplice)>0) and (($pre_per_aberrantsplice+$pre_per_normalsplice)>0)) {
			$pathogenic_score=sprintf("%.1f",($per_aberrantsplice/($per_aberrantsplice+$per_normalsplice)-$pre_per_aberrantsplice/($pre_per_aberrantsplice+$pre_per_normalsplice))*100/2) if ($kind eq "SS_Create");
		}
		if ((($per_aberrantsplice+$per_normalsplice+$per_nosplice)>0) and (($pre_per_aberrantsplice+$pre_per_normalsplice+$pre_per_nosplice)>0)) {
			$pathogenic_score=sprintf("%.1f",($pre_per_normalsplice/($per_aberrantsplice+$per_normalsplice+$per_nosplice)-$per_normalsplice/($pre_per_aberrantsplice+$pre_per_normalsplice+$pre_per_nosplice))*100/2) if ($kind eq "SS_Broken");
		}
	}
	
	print FILEOUT "$kind\t$sum\t$nosplice\t$normalsplice\t$aberrantsplice\t$per_nosplice\t$per_normalsplice\t$per_aberrantsplice\t$ratio\t$ratio2\t$p_val\t$pathogenic_score\n";
	$pre_per_nosplice=$per_nosplice;
	$pre_per_normalsplice=$per_normalsplice;
	$pre_per_aberrantsplice=$per_aberrantsplice;

}

close(FILEOUT);
} #addstats end


print "Done.\n";

sub atgc {
	my $input=shift;
	$input=uc $input;
	#print "Input seq:$input\n";
	my $return="";;
	for my $I (0..length($input)) {
 		my $ch=substr($input,$I,1);
 		#print "subtract character:$ch\n";
 		if ($ch=~/(A|T|G|C)/) {
 			$return="$return$ch";
 		}
	}
	return($return);

}
sub revcom {
	my $input=shift;
	$input=uc $input;
	#print "Input seq:$input\n";
	my $return="";;
	for my $I (0..length($input)) {
 		my $ch=substr($input,length($input)-$I,1);
 		#print "subtract character:$ch\n";
 		if ($ch=~/(A|T|G|C)/) {
 			my $add;
 			$add="A" if ($ch eq "T");
 			$add="T" if ($ch eq "A");
 			$add="G" if ($ch eq "C");
 			$add="C" if ($ch eq "G");
 			$return="$return$add";
 		}
	}
	return($return);
}

sub trimseq {  #&trimseq($SPL[4+$k],$SPL[0])
#my $nleft=4; #mer of 3 side
#my $nright=4;  #mer of 5 side
	my $seq=uc(shift);
	my $revseq=$seq;
	my $varname=shift;
	my $refseq=uc($refseq{$varname});
	my $spseq1=uc($spseq1{$varname});
	my $spseq2=uc($spseq2{$varname});
	print "Given Seq:$seq\tCut:";
	my $unique=1;
	my $rcut=0;my $lcut=0;
	if (&checkuniqe($refseq,$spseq1,$spseq2,$seq)==0) {
		print "Given seq is NOT unique\n";
		return ($seq);
	};
	
    if ($seq=~/^(<ANONYMOUS>)/) {
    	$seq=~s/^(<ANONYMOUS>)//;
    	
    	while (($unique==1) and ($rcut<(length($seq)-1)) and ($rcut<=$maxcut) ) {
    		my $tempseq=substr($seq,0,length($seq)-1);
    		$unique=&checkuniqe($refseq,$spseq1,$spseq2,$tempseq);
    		if ($unique==1) {
    			$seq=$tempseq;
	    		$rcut++;
    			print "r$rcut ";
    		}
    	}
		$seq="<Anonymous>$seq";    	
    } elsif ($seq=~/(<ANONYMOUS>)$/) {
    	$seq=~s/(<ANONYMOUS>)$//;
    	while (($unique==1) and ($lcut<(length($seq)-1)) and ($lcut<=$maxcut) ) {
    		my $tempseq=substr($seq,1,length($seq)-1);
    		$unique=&checkuniqe($refseq,$spseq1,$spseq2,$tempseq);
    		if ($unique==1) {
    			$seq=$tempseq;
    			$lcut++;
    			print "l$lcut ";
    		}
    	}
		$seq="$seq<Anonymous>";    	
    } else {
		while (($unique==1) and ($rcut<($nright-1)) and ($lcut<($nleft-1)) and (($rcut+$lcut)<(length($seq)-2)) and ($rcut<=$maxcut)  and ($lcut<=$maxcut) ) {
			my $tempseq=substr($seq,1,length($seq)-2);
    		$unique=&checkuniqe($refseq,$spseq1,$spseq2,$tempseq);
    		if ($unique==1) {
    			$seq=$tempseq; 
    			$rcut++;	$lcut++;
    			print "l$lcut r$rcut "; 		
    		} else {
    			$tempseq=substr($seq,1,length($seq)-1);
    			$unique=&checkuniqe($refseq,$spseq1,$spseq2,$tempseq);
    			if ($unique==1) {
    				$seq=$tempseq; 
    				$lcut++;		
    				print "l$lcut ";
    			} else {
    				$tempseq=substr($seq,0,length($seq)-1);
    				$unique=&checkuniqe($refseq,$spseq1,$spseq2,$tempseq);
    				if ($unique==1) {
    					$seq=$tempseq; 
    					$rcut++;	
    					print "r$rcut ";	
    				}
    			}
    		}
		}
    	
    } 
    if ($attachextra>0) {
    	print "\tAttach $attachextra ExtraSeq to UniqueSeq.";
    	$lcut=$lcut-$attachextra;
    	$rcut=$rcut-$attachextra;
    	$lcut=0 if ($lcut<0);
    	$rcut=0 if ($rcut<0);
    }	
    $seq=substr($revseq,$lcut,length($revseq)-$rcut-$lcut);
	print "\tSeq for Search:$seq\n";
	return ($seq);
}

sub checkuniqe{ #=&checkuniqe($refseq,$spseq1,$spseq2,$tempseq);  
 my $refseq=uc(shift);#refseq1  no splicing seq
 my $spseq1=uc(shift);#refseq2  normal splicing seq
 my $spseq2=uc(shift);#refseq3  aberrant splicing seq
 my $tempseq=uc(shift); #short seq that is tested if it is unique in the above 3 refseqs
 my $flagref=1; #uniquflag for refseq1
 my $flagsp1=1;#uniquflag for refseq2
 my $flagsp2=1;#uniquflag for refseq2
 my $flaggrep=1;
 my $grephit=0;
 #print " othergroups:",join(",",@othergreps)," ";
 for my $I (@othergreps) { #check the redundancy with other grep characters 
 	$grephit++ if ($I=~/$tempseq/);
 }
 $flaggrep=0 if ($grephit>1);
 
 if ($refseq=~/$tempseq/) {
    $refseq=~s/$tempseq//;
    if ($refseq=~/$tempseq/) {
    	$flagref=0; #tempseq found twice in the refseq
    } 
 }
 if ($spseq1=~/$tempseq/) {
    $spseq1=~s/$tempseq//;
    if ($spseq1=~/$tempseq/) {
    	$flagsp1=0; #tempseq found twice in the refseq
    } 
 }
 if ($spseq2=~/$tempseq/) {
    $spseq2=~s/$tempseq//;
    if ($spseq2=~/$tempseq/) {
    	$flagsp2=0; #tempseq found twice in the refseq
    } 
 }

 return($flagref*$flagsp1*$flagsp2*$flaggrep);
 
}

sub extractnum {
	my $a=shift;
	$a=~/\S+% \((\d+)\)/;
	return($1);
}


sub fishertest{
	my $caseyes=shift;
	my $caseno=shift;
	my $controlyes=shift;
	my $controlno=shift;
	my $R = Statistics::R->new();
	#my $caseno = $ncase - $caseyes;
	#my $controlno = $ncontrol - $controlyes;
	$R->set('a1', $caseyes);
	$R->set('a2', $caseno);
	$R->set('a3', $controlyes);
	$R->set('a4', $controlno);
	$R->run(q`b <- matrix(c(a1,a2,a3,a4),nrow=2)`);
	$R->run(q`c <- fisher.test(b)`);
	$R->run(q`d <- c$p.value`);
	$R->run(q`e <- c$conf.int[1]`);
	$R->run(q`f <- c$conf.int[2]`);
	$R->run(q`g <- as.numeric(c$estimate)`);
	my $pval = $R->get('d');
	my $L95 = $R->get('e');
	my $U95= $R->get('f');
	my $OR=$R->get('g');
	$R->stop();
	return($pval);
}
