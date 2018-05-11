###2028_04_26 SplitFusion Python 2.7###
import os
import sys
import subprocess
from subprocess import Popen, PIPE
import HTSeq
import argparse
parser = argparse.ArgumentParser(description='SplitFusion v1')
import itertools
from itertools import groupby


###1 get all SA reads from bam ###
### sort SA reads with read ID ###
wkdir="/home/athena/Splitfusion/SplitFusion_ZLZ/example-run-folder"
in_bam="A01-P701.consolidated.bam"

if "bam" in in_bam:
    in_samSA=in_bam.replace("bam", "sam")
    Grep_SA="samtools view " + in_bam  + " | grep \'SA:\' > " + wkdir + "/" + in_samSA
    Grep_=Popen(Grep_SA, shell=True)
    Grep_.communicate()
    Sort_SA="sort -k1 " + in_samSA + " > " +  wkdir + "/" + in_samSA + ".sorted" 
    Sort_=Popen(Sort_SA, shell=True)
    Sort_.communicate()   
elif "sam" in in_bam:
    Sort_SA="sort -k1 " + in_bam + " > " +  wkdir + "/" + in_samSA + ".sorted" 
    Sort_=Popen(Sort_SA, shell=True)
    Sort_.communicate()   
else:
    print "please input bam file"
    sys.exit
    

### 2 sort out left and right partition of the read ###        
def query_orientation_check(bundle, gap_max, overlap_max, min_exclusive):
    """check overlapping bp and no. of gaps between matches of alignment"""
    """correct for h-clipped base that has 0 cooridinate, messing up the orientations"""
    """ignore gap at ends"""
    """filter reads if they are not exclusive >= 25 bp"""
    Matches=HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    overlap=0
    gap=0
    no_align=len(bundle)
    match_cood=[]
    for i in range(no_align):
        ind=i
        H_length=0
        seq_length=0
        for c in bundle[i].cigar:
            if c.type=="M":            
                Matches[HTSeq.GenomicInterval("q", H_length+c.query_from, H_length+c.query_to, c.ref_iv.strand) ] += str(ind) 
                match_cood.append([[c.ref_iv.chrom, c.ref_iv.start], H_length+ c.query_from, [c.ref_iv.chrom, c.ref_iv.end], H_length+ c.query_to])                
            if c.type=="H":            
                H_length+=c.size
            seq_length+=c.size
    cigar_pattern=zip([len(s) for s in list(Matches[HTSeq.GenomicInterval("q", 0, seq_length, "+")])], [len(s) for s in reversed(list(Matches[HTSeq.GenomicInterval("q", 0, seq_length, "-")]))]) 
    cigar_pattern_counts=[[_, sum(1 for _ in group)] for _, group in groupby(cigar_pattern)]  
    exclusive=([i for i in cigar_pattern_counts if i[0]==(0, 1) or i[0]==(1, 0)])
    left_exclusive=exclusive[0][1]
    right_exclusive=exclusive[-1][1]
    overlap += sum([i[1] for i in cigar_pattern_counts if sum(i[0])>no_align-1])
    gap += sum([i[1] for i in cigar_pattern_counts[1:-1] if sum(i[0])==0])
    match_cood=sorted(match_cood, key=lambda x: (x[1], x[0][1]))
    ### deal with overlap: left_split - overlap, right_split + overlap ###
    if no_align ==2:    
        match_cood[0][2][1]-=overlap
        match_cood[-1][0][1]+=overlap
    match_cood[0][0]=tuple(match_cood[0][0])
    match_cood[0][2]=tuple(match_cood[0][2])
    match_cood[-1][0]=tuple(match_cood[-1][0])
    match_cood[-1][2]=tuple(match_cood[-1][2])
    if gap < gap_max and overlap < overlap_max and left_exclusive>=min_exclusive and right_exclusive>=min_exclusive and len(exclusive) >1:
        return (True, match_cood)
    else:
        return (False, None) 


 

###1 only get SA bundles >= mapq ###
in_sam=HTSeq.SAM_Reader(wkdir + "/" + in_samSA + ".sorted")
gap_max=2
overlap_max=10
mapq=20
M_length=25
min_exclusive=25

SA_bundles={}
### only save alignment of a read if > mapq###
for a in HTSeq.bundle_multiple_alignments(in_sam):
    j=[i for i in a if i.aQual >= mapq]
    l=[]
    for aln in a:
        match_base=0
        for cigar in aln.cigar:
            if cigar.type=="M":
                match_base+=cigar.size
        if match_base>=M_length:
                l.append(cigar)
    if len(j) > 1 and len(l) >1:
        query_pass, match_cood=query_orientation_check(j, gap_max, overlap_max, min_exclusive)
        if query_pass:
            SA_bundles[j[0].read.name]=[j, match_cood]
            

###2.5 validate with the answer ###
b=open("/home/athena/Splitfusion/SplitFusion/example-run-folder/fusion-output/A01-P701/breakpoint", "r")
breakpoints=[]
for row in b:
    row=row.replace("\n", "")
    row=row.split(" ")
    breakpoints.append(row)

fu_Overlap="/home/athena/Splitfusion/SplitFusion_ZLZ/example-run-folder/fusionID.overlap"
fusion_ID=[]
b=open(fu_Overlap, "r")
for row in b:
    row=row.replace("\n", "")
    row=row.split(" ")
    fusion_ID.append(row[0])

fusion_1=[i[0] for i in breakpoints if i[-1] in fusion_ID] 

SA_names=SA_bundles.keys()
left_behind=set(fusion_1) - set(SA_names)
### passed reads sample [:5] ###
### ['CL100059760L2C008R091_453362', 'CL100059760L2C003R079_93091', 'CL100059760L2C006R063_212134', 'CL100059760L2C011R064_252699', 'CL100059760L2C008R037_409337']

### failed false negative compared to fusionID_overlap ###
### 'CL100059760L2C007R018_517271' --- failed the exclusive filter but was in the fusionIDoverlap file... ###
### 'CL100059760L2C015R076_146009' 
### 'CL100059760L2C012R028_513555'
### 'CL100059760L2C016R049_90964'
for i in list(left_behind)[100:500]:
    bundle=[a for a in HTSeq.bundle_multiple_alignments(in_sam) if a[0].read.name==i][0]
    j=[i for i in bundle if i.aQual >= mapq]
    q_pass, split_coor=query_orientation_check(j, gap_max, overlap_max, min_exclusive)
    if q_pass:
        print i    


####3. define, sort and cluster split points ###
### report readID, left split, right split ###
left_breakpoints=HTSeq.GenomicArrayOfSets( "auto", stranded=False )
right_breakpoints=HTSeq.GenomicArrayOfSets( "auto", stranded=False )
for i in SA_bundles.keys():
    left_breakpoints[HTSeq.GenomicInterval(SA_bundles[i][1][0][2][0], SA_bundles[i][1][0][2][1], SA_bundles[i][1][0][2][1]+1, strand='.')]+=i
    right_breakpoints[HTSeq.GenomicInterval(SA_bundles[i][1][-1][0][0], SA_bundles[i][1][-1][0][1], SA_bundles[i][1][-1][0][1]+1,strand='.')]+=i

###group and filter splits###
###3.1 group_splits ###
def group_splits(breakpoints):
    """group_split if within 5bp apart"""
    """only keep group having unique mapped reads >= 3"""
    groups=[]
    split_group=[]
    previous=None
    for i in breakpoints.steps():
        if len(i[1])>0:
            if previous==None:
                split_group.append(i)
            elif previous!=None:
                if previous[0].chrom==i[0].chrom and abs(previous[0].start - i[0].start) <= 5:
                    split_group.append(i)
                else:
                    groups.append(split_group)
                    split_group=[i]
            previous=i
    return groups

left_splits_clustered=group_splits(left_breakpoints)
right_splits_clustered=group_splits(right_breakpoints)

###3.2 filter for stagerring ends ###
def filter_stagerring_ends(splits_clustered, SA_bundles, orientation):
    """remove grouped reads if no. staggering < 3"""
    filtered_list=[]
    for i in splits_clustered:
        reads=[]
        bk_list=[]
        for j in i:
            reads.extend(list(j[1]))
        left_stagger=[]
        right_stagger=[]
        for r in reads:
            left_stagger.append(SA_bundles[r][1][0][0])
            right_stagger.append(SA_bundles[r][1][-1][2])
        if orientation == "left":
### so left_stagger > 3 and right_stagger ==1###
            stag_support=len(set(left_stagger))
            if stag_support >=3 and len(set(right_stagger))==1:
                for j in i:
                   bk_list.append((j[0], stag_support))
        if orientation == "right":
### so right_stagger > 3 and left_stagger ==1###
            stag_support=len(set(right_stagger))
            if stag_support >=3 and len(set(left_stagger))==1:
                for j in i:
                    bk_list.append((j[0], stag_support))
        if len(bk_list)>0:
            filtered_list.append(bk_list)
    return filtered_list


left_splits_clustered0=filter_stagerring_ends(left_splits_clustered, SA_bundles, "left")
right_splits_clustered0=filter_stagerring_ends(right_splits_clustered, SA_bundles, "right")

###4 identify the other end of breakpoints from splits_clustered ###
def collect_reads_from_splits(a_split_list, breakpoints):     
    boundary=[]
    for GIO in a_split_list:
        r=[]   
        GI=GIO[0]
        score=GIO[1]
        for j in breakpoints[GI].steps():
            r.extend(list(j[1]))
        break_boundary=retrieve_break_from_reads(SA_bundles, r)
        boundary.append([GIO[1], break_boundary])
    return boundary

def retrieve_break_from_reads(SA_bundles, reads):
    """retrieve both left and right boundaries"""
    """allow left/right split to link to different fusion"""
    break_cood=[]
    for i in reads:
        break_cood.append("-".join(("_".join(str(w) for w in SA_bundles[i][1][0][2]), "_".join(str(w) for w in SA_bundles[i][1][1][0]))))
    break_cood.sort()
    return [_ for _, group in groupby(break_cood)]

def define_break_points_boundary(splits_clustered, breakpoints, SA_bundles, orientation):
    """extract genomic interval from splits_clustered"""
    """extract reads ID from breakpoints array"""
    """retrieve breakpoint info from SA_bundles"""
    """group breakpoints - assign groups"""
    """sort breakpoints by support"""
    break_list=[]    
    for i in splits_clustered:
        reads=collect_reads_from_splits(i, breakpoints)
        split_support=sum([u[0] for u in reads])
        i=sorted(i, key=lambda w: (w[0].chrom, w[0].start))
        split_bin="-".join(("_".join((str(i[0][0].chrom), str(i[0][0].start))), "_".join((str(i[-1][0].chrom), str(i[-1][0].start)))))
        break_list.append([split_bin, reads, split_support])
    return break_list

for i in left_splits_clustered0:
        reads=collect_reads_from_splits(i, breakpoints)
        split_support=sum([u[0] for u in reads])
        i=sorted(i, key=lambda w: (w[0].chrom, w[0].start))
        split_bin="-".join(("_".join((str(i[0][0].chrom), str(i[0][0].start))), "_".join((str(i[-1][0].chrom), str(i[-1][0].start)))))
        break_list.append([split_bin, reads, split_support])

left_breaklist=define_break_points_boundary(left_splits_clustered0, left_breakpoints, SA_bundles, "left")
right_breaklist=define_break_points_boundary(right_splits_clustered0, right_breakpoints, SA_bundles, "right")

###5 input annotation of exons###       
### gene features from UCSC table browser, Genes and gene prediction  , UCSC genes, known Genes, output as all fields from selected table ###
### building gnee model is computationally expensive, so loop through breaklist to retrieve genes spanning breakpoints, then build gene model for each gene, then retrieve relevant codon-partition for the breaks ###
### the feat.txt file appears to be 0 base###

GFF_file='/home/athena/Refseq/HG19/UCSCGenes_HG19_feat.txt'

def read_in_feature_file(GFF_file):
    """read in tab delimited file"""
    GFF_dict={}
    transcript=HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    with open(GFF_file, "r") as infile:
        for row in infile:
            row=row.replace("\n", "")
            if "#" in row:
                header=row.replace("#", "")
                header=header.split("\t")
            else:
                row=row.split("\t")
                row_dict={k:v for k, v in zip(header,row)}
                if row_dict["proteinID"]!="":
                    row_dict["chrom"]=row_dict["chrom"].replace("chr", "")
                    row_dict["txStart"]=int(row_dict["txStart"])
                    row_dict["txEnd"]=int(row_dict["txEnd"])    
                    GFF_dict[row_dict["name"]]=row_dict
                    transcript[HTSeq.GenomicInterval( row_dict["chrom"], row_dict["txStart"], row_dict["txEnd"], "." )]+=row_dict["name"]
    return GFF_dict, transcript

GFF_dict, transcript=read_in_feature_file(GFF_file)

#### Refseq version###
GFF_file='/home/athena/Refseq/HG19/NCBIRef_seq_HG19_fest.txt'
def read_in_feature_file(GFF_file):
    """read in tab delimited file"""
    GFF_dict={}
    transcript=HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    with open(GFF_file, "r") as infile:
        for row in infile:
            row=row.replace("\n", "")
            if "#" in row:
                header=row.replace("#", "")
                header=header.split("\t")
            else:
                row=row.split("\t")
                row_dict={k:v for k, v in zip(header,row)}
                if row_dict["cdsStart"]!= row_dict["cdsEnd"]:
                    row_dict["chrom"]=row_dict["chrom"].replace("chr", "")
                    row_dict["txStart"]=int(row_dict["txStart"])
                    row_dict["txEnd"]=int(row_dict["txEnd"])    
                    GFF_dict[row_dict["name"]]=row_dict
                    transcript[HTSeq.GenomicInterval( row_dict["chrom"], row_dict["txStart"], row_dict["txEnd"], "." )]+=row_dict["name"]
    return GFF_dict, transcript

GFF_dict, transcript=read_in_feature_file(GFF_file)

####5.0.1 build model for gene transcript ###
def build_gene_model(g, GFF_dict):
    """return gene model of a gene"""
    """define with codon_no, and codon partition"""
    gene_model=HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    exon_no=int(GFF_dict[g]['exonCount'])
    exon_start=[int(j) for j in GFF_dict[g]['exonStarts'].split(",")[:exon_no]]
    exon_end=[int(j) for j in GFF_dict[g]['exonEnds'].split(",")[:exon_no]]
    if GFF_dict[g]['strand']=="-":
        start_codon=int(GFF_dict[g]['cdsEnd'])
        stop_codon=int(GFF_dict[g]['cdsStart'])
        exon_start=list(reversed(exon_start))
        exon_end=list(reversed(exon_end))
    else:
        start_codon=int(GFF_dict[g]['cdsStart'])
        stop_codon=int(GFF_dict[g]['cdsEnd'])
    start_exon=[[s, e] for s, e in zip(exon_start, exon_end) if s<=start_codon and e>start_codon][0]
    end_exon=[[s, e] for s, e in zip(exon_start, exon_end) if s<stop_codon and e>=stop_codon][0]
    in_between_codon=[[s, e] for s, e in zip(exon_start, exon_end)]
    Start_index=in_between_codon.index(start_exon)
    End_index=in_between_codon.index(end_exon)
    if GFF_dict[g]['strand']=="-":
        start_exon=(start_exon[0], start_codon)
        end_exon=(stop_codon, end_exon[1])
    else:
        start_exon=[start_codon, start_exon[1]]
        end_exon=[end_exon[0], stop_codon]
    exons_cood=[start_exon]
    exons_cood.extend(in_between_codon[Start_index+1:End_index])
    exons_cood.append(end_exon) 
    cDNA_part=0 
    exon_no=0  
    codon_n=1
    codon_partition=0
    if GFF_dict[g]['strand']=="-":
        for i in exons_cood:
            exon_no+=1
            for location in list(reversed(range(i[0], i[1]))):
                codon_partition+=1
                cDNA_part+=1
                if codon_partition==4:
                    codon_n+=1
                    codon_partition=1
                in_name=str(exon_no)+"_"+str(codon_n)+"_"+str(cDNA_part)+"_"+str(codon_partition)
                gene_model[HTSeq.GenomicInterval(GFF_dict[g]['chrom'], location, location+1)]+=in_name
    else:
        for i in exons_cood:
            exon_no+=1
            for location in range(i[0], i[1]):
                codon_partition+=1
                cDNA_part+=1
                if codon_partition==4:
                    codon_n+=1
                    codon_partition=1
                in_name=str(exon_no)+"_"+str(codon_n)+"_"+str(cDNA_part)+"_"+str(codon_partition)
                gene_model[HTSeq.GenomicInterval(GFF_dict[g]['chrom'], location, location+1)]+=in_name
    return gene_model

###5.1 Annotate breakpoint###                
def extract_break_genes(bk_point, transcript):
    """find genes where breaks happen."""
    """build gene model"""
    """if break in exon, report exon number and codon position"""
    """if break in intron, report closest intron"""
    bk_point=bk_point.split("-")
    break5=bk_point[0].split("_")
    break3=bk_point[-1].split("_")
    gene5=list(transcript[HTSeq.GenomicPosition(break5[0], int(break5[1]), ".")])
    gene3=list(transcript[HTSeq.GenomicPosition(break3[0], int(break3[1]), ".")])
    return HTSeq.GenomicPosition(break5[0], int(break5[1]), "."), HTSeq.GenomicPosition(break3[0], int(break3[1]), "."), gene5, gene3

def determine_break_location(break_5, g_model):
    """identify where break position contain"""
    """ if in intron: set([]), identify intron no. """
    bk_location=g_model[break_5]
    if bk_location==set([]):
        try:
            upstream5=list([i for i in g_model.steps() if i[0].end < break_5.start][-1][1])[0]
            upstream3=list([i for i in g_model.steps() if i[0].start > break_5.start][0][1])[0]    
            bk_location=upstream5.split("_")[0] + "_" + "intron" + "_" + upstream3.split("_")[0] 
        except:
            bk_location="UTR"      
    else:
        bk_location=list(bk_location)[0] 
    return bk_location

def get_split_info_from_GFF(GFF_dict, transcript, breaklist):
    """export gene-exon cood for each split point"""
    """test if the split points are in-frame"""
    breaks_per_bin={}
    genes_model={}
    for i in breaklist:
        breaks_per_bin[i[0]]={}   
        breaks=i[1]
        for b in breaks:
            for bk_point in b[1]:
                ###print bk_point
                break_5, break_3, gene5, gene3=extract_break_genes(bk_point, transcript)
                #### have to sort out codes here ####
                break_5_locations=[]
                break_3_locations=[]
                for g in gene5:
                    ###print g          
                    try:
                ### if g already in genes_model then do the annotation###
                        g_model=genes_model[g]
                        break_5_locations.append([g, determine_break_location(break_5, g_model)])
                    except:
                ### if not in genes_model then build gene model and do annotation
                        g_model=build_gene_model(g, GFF_dict)
                        genes_model[g]=g_model
                        determine_break_location(break_5, g_model)
                        break_5_locations.append([g, determine_break_location(break_5, g_model)])
                for g in gene3:
                    ###print g 
                    try:
                ### if g already in genes_model then do the annotation###
                        g_model=genes_model[g]      
                        break_3_locations.append([g, determine_break_location(break_3, g_model)])
                    except:
                ### if not in genes_model then build gene model and do annotation
                        g_model=build_gene_model(g, GFF_dict)
                        genes_model[g]=g_model
                        determine_break_location(break_5, g_model)
                        break_3_locations.append([g, determine_break_location(break_3, g_model)])             
                try:
                    breaks_per_bin[i[0]][bk_point]=[break_5_locations, break_3_locations]
                except:
                    breaks_per_bin[i[0]]={bk_point:[break_5_locations, break_3_locations]}
    return breaks_per_bin
                

left_breaks_per_bin=get_split_info_from_GFF(GFF_dict, transcript, left_breaklist)
right_breaks_per_bin=get_split_info_from_GFF(GFF_dict, transcript, right_breaklist)

#Test code to ensure the exon_no and cds inference is correct
# g='NM_004521.2'
# g_model=build_gene_model(g, GFF_dict)
# g_model[HTSeq.GenomicPosition(GFF_dict[g]['chrom'], 32317354)]
# test_left=left_breaklist[:1]
#[i for i in left_breaks_per_bin.keys() if "10_" in i[0]]
# left_breaks_per_bin=get_split_info_from_GFF(GFF_dict, transcript, test_left)
# left_breaks_per_bin['15_88669605-15_88669610']
# left_breaks_per_bin['10_43607554-10_43607556']

# right_breaks_per_bin['10_43608912-10_43608920']
# right_breaks_per_bin['1_27951600-1_27951600']

### 5 merge left and right_split, filter if present in both side ###
def organise_split_junction(right_breaklist, right_breaks_per_bin,left_breaklist, left_breaks_per_bin, orientation, GFF_dict):
    """ if left and right junctions were common, remove the bin"""
    """input the break-list and breaks_per_bin"""
    """for each break junction, output the set of left/right split annotations"""
    """column: #bin left/right split unimolsupport subsplit leftjunctiongene leftjunctionlocation rightjuctiongene right junctionlocation"""
    
    for i in break_list:
        break_key=i[0]
        break_support=i[-1]
        c=break_per_bin[break_key]
        for bk_points in c.keys():
            subsplit=bk_point
            left_jun=c[bk_points][0]
            right_jun=c[bk_points][-1]
            junctions=def combine_junctions(left_jun, right_jun, GFF_dict)




