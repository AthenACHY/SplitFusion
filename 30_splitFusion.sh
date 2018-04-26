#!/bin/bash

## from BAM SA hits to structure variation (fusion, cDNA or large indel)
## input:
##	- bam/id_folder (current folder)
##	- bam/id_folder/consolidated.sam 
##
## output:
##	- strV_read_ID
##	- sa.cir0 (pre-filtered)
##	- sa.cir.panel (filtered)
##		## at least 25 bp exculsive between hit1 and hit2
##    		## at least 5 (or x$) split reads
##		## invovles targeting site (+/- 1kb)
##		## filtered within fusion < minBreakDistance


# enter your working directory
. ../run.info.sh
subii='SUBID_'
mkdir -p $subii
cd  $subii

    ## determine parallel job number from the cpuBWA field of samplesheet file
    head -n 1 ../../iiFreq.txt > _head.1
    gawk '{for (i=1; i<=NF; i++) {if ($i ~ /cpuBWA/) {print $i,i}}}' _head.1 | sed 's/.* /cpuField=/' > _t.sh
    . _t.sh
    cpuBWA=$(grep $subii ../../iiFreq.txt | cut -f $cpuField | sed 's: .*::')

#==========================
if [ ! -f breakpoint ]; then
#==========================

## get reads with SA
$REPPATH/samtools/samtools view ../../bam/$subii.consolidated.bam | grep 'SA:' > _sa.sam


    $REPPATH/samtools/samtools view --threads $cpuBWA -T $DEPATH/Homo_sapiens_assembly19.fasta -bS _sa.sam > _sa.bam
    bedtools bamtobed -cigar -i _sa.bam > _sa.bed0
	## bedtools uses 0-base, change to 1-base:
	awk '{start = $2+1; $2=start; print}' _sa.bed0 | tr ' ' '\t' > _sa.bed
	sort -k4,4b _sa.bed > _sa.bed.s

    ## read length
    awk '{n=length($10); print $1,n}' _sa.sam > _sa.len
	sort -k1,1b -k2,2nr _sa.len > _sa.len.s0
	sort -k1,1b -u _sa.len.s0 > _sa.len.s

    join -1 1 -2 4 _sa.len.s _sa.bed.s > _sa.len2

    ## MQ 
    echo | awk -v minMQ=$minMQ '{if ($6 >= minMQ) print}' _sa.len2 | tr ' ' '\t' > _sa.bed.mq

## corrected read length for small size clipping which is either
		## 1)  likely residual adaptor
		## 2)  shorter splice fragments (< minMappLength) that could not be mapped by BWA MEM

		# transform CIGAR strings
	sed 's:\([0-9]\+\)\([SH]\)\([0-9A-Z]\+\)$:\1\t\2\t\3:' _sa.bed.mq > _sa.SMH01
	sed 's:\([0-9]\+\)\([SH]\)$:\t\1\t\2:' _sa.SMH01 > _sa.SMH02

		# two splits (10 fileds only) to _tmpOK
		# complicated ones to _tmp1 for further processing
	awk '{if ($12 ~ /[HS]/) {print > "_tmp1"} else {print $1,$2,$3,$4,$5,$6,$7,$8$9$10 > "_tmpOK"} }' _sa.SMH02 

		# mid split
	minMapLength=20
	echo | awk -v minMapLength=$minMapLength '{if ($8 >= minMapLength && $11 >= minMapLength) {print > "split.mid"} else {print > "_tmp2"} }' _tmp1

		# remove headi/tail HS and correct effective query read length
    	awk '{if ($8 >= $11) {$2 = $2-$11; newcigar=$8$9$10} else {$2 = $2-$8; newcigar=$10$11$12} print $0,newcigar}' _tmp2 | cut -d ' ' -f 1-7,13 >  _tmpOK2

		# Corrected SMH, with left size and SMH type separated for later processing. 
	cat _tmpOK _tmpOK2 | tr ' ' '\t' | sed 's:\([0-9]\+\)\([SMH]\)\([0-9A-Z]\+\)$:\1\t\2\t\3:' > _sa.SMHc


	##==========================================
	## calculate query start, end
    	awk '{if ($9 == "M") print}' _sa.SMHc | tr ' ' '\t' | cut -f 1-9 > _M
    	awk '{if ($9 != "M") print}' _sa.SMHc | tr ' ' '\t' | cut -f 1-9 > _HS
	
	## print chr.pos 1,2,3,4 in read order, re-calculate read query start and end for '-'
		## add read part (i.e left/right at field $12, breadkpoint chr ($13) and pos ($14)

	    ## if left of query is mapped
		    ## Left: $7+ and M; $7- and HS
		    ## Right: $7- and M; $7+ and HS
	    awk '{if ($7 =="+") {qstart=1; qend=$5-$4+1; print $0,qstart,qend,"left",$3,$5}}' _M | tr ' ' '\t' > _left
	    awk '{if ($7 =="-") {qstart=1; qend=$5-$4+1; mstart=$5; mend=$4; $4=mstart; $5=mend; print $0,qstart,qend,"left",$3,$5}}' _HS | tr ' ' '\t' >> _left

	    ## if right of query is mapped
	    awk '{if ($7 =="+") {qstart=$8+1; qend=$8+$5-$4; print $0,qstart,qend,"right",$3,$4}}' _HS | tr ' ' '\t' > _right
	    awk '{if ($7 =="-") {qstart=$2-$5+$4; qend=$2; mstart=$5; mend=$4; $4=mstart; $5=mend; print $0,qstart,qend,"right",$3,$4}}' _M | tr ' ' '\t' >> _right

	## adjust breakpoint position to group breakpoints fuzzy window (5 bases)
	cat _left _right > _left_right0
	sort -k13,13b -k14,14n _left_right0 > _left_right

			    awk '{if ($13==pre13) {
				    d=$14-pre14; if (d<5){
							    posAdj=preAdj
						    } else {posAdj=$14}	
				    }else{posAdj=$14}
				    ; pre13=$13; pre14=$14; preAdj=posAdj; print $0,posAdj
				    }' _left_right > _left_right.adj
			    awk '{print $0,$14}' _left_right > _left_right.adj

	    ## join left and right to form a query read
		## for pair-end read 2, reverse left <-> right
	    grep 'left'  _left_right.adj > _lefts
	    grep 'right' _left_right.adj > _rights

		sort -k1,1b _lefts > _leftsrt
		sort -k1,1b _rights > _rightsrt

		## $1: readID
		## $2 - $15: left
		## $16 - $29: right
	    join _leftsrt _rightsrt > _left_right.j

	    ##====== breakpoint 
	    ## add at $30, $31, $32
	    awk '{pLeft=$13"_"$15; pRight=$27"_"$29; if (pLeft < pRight){pp=pLeft"-"pRight} else {pp=pRight"-"pLeft}; print $0,pLeft,pRight,pp}' _left_right.j > breakpoint0
		sort -k1,1b breakpoint0 > breakpoint
#==========================
fi
#==========================

		## reads with middle split
		cut -f1 split.mid | sort -u | awk '{print $0,"_mid"}' > _mid.id
		
		join -a1 breakpoint _mid.id > breakpoint2
		grep '_mid$' breakpoint2 | sed 's:_mid$::' > breakpoint.w.mid
		grep -v '_mid$' breakpoint2 > breakpoint.wo.mid

		## 1. at least minMapLength 
		## 2. at least minExclusive bp exculsive between hit1 and hit2
		## 3. no larger than maxQueryGap
		## 4. less than max overlapping length

		## reads without middle split
		    ## left:  $10-------$11
		    ## right:       $24------$25
		 echo | awk -v minMapLength=$minMapLength -v minExclusive=$minExclusive -v maxQueryGap=$maxQueryGap -v maxOverlap=$maxOverlap \
		    '{ gap = $24-$11-1; 
		       overlap = $11-$24+1;
			    if (  ($11-$10 >= minMapLength && $25-$24 >= minMapLength) \
				    && ($24-$10 >= minExclusive && $25-$11 >= minExclusive && gap <= maxQueryGap && overlap <= maxOverlap) \
			) {print $0,overlap}}' \
			     breakpoint.wo.mid > _sa.fu0

		## reads with middle split, turn off maxQueryGap by let = 100
		    ## left:  $10-------$11
		    ## right:		       $24------$25
		 echo | awk -v minMapLength=$minMapLength -v minExclusive=$minExclusive -v maxQueryGap=100 \
		    '{ gap = $24-$11-1; 
		       overlap = $11-$24+1;
			    if (  ($11-$10 >= minMapLength && $25-$24 >= minMapLength) \
				    && ($24-$10 >= minExclusive && $25-$11 >= minExclusive && gap <= maxQueryGap && overlap <= maxOverlap) \
			) {print $0,overlap}}' \
			     breakpoint.w.mid >> _sa.fu0

			## keep fusionID, gap, overlap
			cut -d ' ' -f 32,33 _sa.fu0 > fusionID.overlap


############################
##====== A CRITICAL FILTER:
## Default in run.info.sh: StrVarMinUniqMap=3
    ## or uncomment below and define the value:
    # StrVarMinUniqMap=new.value

	## stagger sites on left $4
	# sort by point-point, left-breakpoint, left-stagger.point
	sort -k32,32b  -k30,30b -k4,4n _sa.fu0 > _sa.fu0s
	sort -k32,32b  -k30,30b -k4,4n -u  _sa.fu0s > _fu_1u	

		## if stagger sites >= minUniqMap
                echo | awk -v minUniqMap=$StrVarMinUniqMap '{if (($32==pre32 && $30==pre30) && ($4 != pre4))
                            {cnt=cnt+1} else {cnt=1};
                            if (cnt == minUniqMap) {print $32,"stgLeft"};
                            pre32=$32; pre30=$30; pre4=$4;
				}' _fu_1u | sort -u > fu.breakpointID.stgLeft
		join -2 32 fu.breakpointID.stgLeft _sa.fu0s > fu.info.stgLeft	
	
	## stagger sites on right $19
	sort -k32,32b -k31,31b -k19,19n -u  _sa.fu0s > _fu_2u	

		## if stagger sites >= minUniqMap
                echo | awk -v minUniqMap=$StrVarMinUniqMap '{if (($32==pre32 && $31==pre31) && ($19 != pre19))
                            {cnt=cnt+1} else {cnt=1};
                            if (cnt == minUniqMap) {print $32,"stgRight"};
                            pre32=$32; pre31=$31; pre19=$19;
				}' _fu_2u | sort -u > fu.breakpointID.stgRight
		join -2 32 fu.breakpointID.stgRight _sa.fu0s > fu.info.stgRight	


	##=== combine fusionID
	cat fu.breakpointID.stgLeft fu.breakpointID.stgRight | cut -d ' '  -f1 | sort -u > fu.breakpointID
	cat fu.info.stgLeft fu.info.stgRight | sort -k1,1b  > fu.info.stgLR0
	join fu.breakpointID fu.info.stgLR0 > fu.info.stgLR

	
##=================================================
##=== Annotate breakpoint gene, exon, cDNA position

		## if stagger sites >= minUniqMap, print breakpoint ($3, $4), and overlap-removed point (orp). Keep both points for later judgement
		### to correctly annotate exon and cDNA position (some alignment ends in intron could be overlap with, and should belong to, next split alignment)
		### , add 6 bases to orp, which needs to be accounted for in later frame determination
		awk '{
			if ($34 >=0){
				    if ($9=="+") {orpL = $7-$34-6};
				    if ($9=="-") {orpL = $7+$34+6};
				    if ($23=="+") {orpR = $20+$34+6};
				    if ($23=="-") {orpR = $20-$34-6};
				}else{
				    if ($9=="+") {orpL = $7-6};
				    if ($9=="-") {orpL = $7+6};
				    if ($23=="+") {orpR = $20+6};
				    if ($23=="-") {orpR = $20-6};
				}
			print $5"_"orpL,$5,orpL,$7,$9,$19"_"orpR,$19,orpR,$20,$23,$34,$1,$2,$3
				}' fu.info.stgLR > fu.breakpoint.for.anno0 

		awk '{print $2,$3,$3,"A","A"}' fu.breakpoint.for.anno0 > _tmp1
		awk '{print $7,$8,$8,"A","A"}' fu.breakpoint.for.anno0 > _tmp2
		cat _tmp1 _tmp2 | sort -u > fu.breakpoint.for.anno

		perl $REPPATH/annovar/table_annovar.pl fu.breakpoint.for.anno $REPPATH/annovar/humandb/ -buildver hg19 -out fu.anno -remove -protocol refGene -operation g -nastring NA
		Rscript $SplitFusionPath/scripts/chr.pos.anno.extraction.R fu.anno ## generate .ext0
		sort -k1,1b fu.anno.ext0 > fu.anno.ext
		
		## readID orp
		cut -d ' ' -f 1,2,4,5,11,14 fu.breakpoint.for.anno0 | sort -u | sort -k1,1b > fu.orpL
		join fu.orpL fu.anno.ext > fu.anno.left

		cut -d ' ' -f 6,7,9,10,11,14 fu.breakpoint.for.anno0 | sort -u | sort -k1,1b > fu.orpR
		join fu.orpR fu.anno.ext > fu.anno.right

		    ## anno middle split
		    awk '{mid = $4 + 10; print $3,mid,mid,"A","A",$1,$4,$5,$6}' split.mid > _mid.for.anno0
		    tr ' ' '\t' < _mid.for.anno0 | cut -f1-5 | sort -u > mid.for.anno

		    perl $REPPATH/annovar/table_annovar.pl mid.for.anno $REPPATH/annovar/humandb/ -buildver hg19 -out mid.anno -remove -protocol refGene -operation g -nastring NA
		    Rscript $SplitFusionPath/scripts/chr.pos.anno.extraction.R mid.anno ## generate .ext0

		    sed 's: :_:' _mid.for.anno0 | sort -k1,1b > _mid.for.anno1
		    sort -k1,1b mid.anno.ext0 > _mid.anno.ext
		    join _mid.for.anno1 _mid.anno.ext | cut -d ' ' -f5- > mid.anno2

## In-frame status determination
		Rscript $SplitFusionPath/scripts/breakpoint.anno.R
		
##=== End of StrVarMinUniqMap
	


#################################	
## Prepare data for circos plot...to add later
#################################	

rm _*

rm ../_job_.$subii
