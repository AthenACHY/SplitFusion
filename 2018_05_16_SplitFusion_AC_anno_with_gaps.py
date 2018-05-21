###2028_04_26 SplitFusion Python 2.7###

###0import modules ###


###1 get all SA reads from bam ###
### sort SA reads with read ID ###

### 2 sort out left and right partition of the read ###  
def calculate_overlaps_between_alignment(partial_patterns, i):
    """assume overlapping only at the end of each alignment"""
    """left_fragment has left_overlap at the 3' end of the fragment """
    """right_fragment has right_overlap at the 5'end of the fragment"""
    left_overlap=0
    right_overlap=0
    right_fragment_ind=None
    left_fragment_ind=None
    for a in partial_patterns:
        if a.index(i)==0:
            right_fragment_ind=int(list(a[-1])[0])
            right_overlap=len([s for s in a if s==i])
        elif list(reversed(a)).index(i)==0:
            left_fragment_ind=int(list(a[0])[0])
            left_overlap=len([s for s in a if s==i])
    return left_fragment_ind, left_overlap, right_fragment_ind, right_overlap

def calculate_exclusive(cigar_pattern, match_cood, cigar_pattern_counts, no_align):
    left_exclusive=[i for i in cigar_pattern_counts if len(i[0])==1][0][1]
    right_exclusive=[i for i in cigar_pattern_counts if len(i[0])==1][-1][1] 
    return left_exclusive, right_exclusive 

###report query partition correctly for each alignment ###

def create_cigar_object(bundle):
    Matches=HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    no_align=len(bundle)
    match_cood=[]
    seq_lengths=[]
    for i in range(no_align):
        ind=i
        H_length=0
        seq_length=0
        for c in bundle[i].cigar:
            if c.type=="M":            
                Matches[HTSeq.GenomicInterval("q", H_length+c.query_from, H_length+c.query_to, c.ref_iv.strand) ] += str(ind)
                match_cood.append([[c.ref_iv.chrom, c.ref_iv.start_d], H_length+ c.query_from, [c.ref_iv.chrom, c.ref_iv.end_d], H_length+ c.query_to, c.ref_iv.strand])                
            if c.type=="H":            
                H_length+=c.size
            seq_length+=c.size
        seq_lengths.append(seq_length)
    return no_align, Matches, match_cood, max(seq_lengths)


def correct_match_cood(match_cood, cigar_pattern, no_align):
    """correct for query partition 0-base strand +"""
    cigar_pattern_ind=zip(range(seq_length), cigar_pattern)
    for n in range(no_align):
        q=str(n)
        part_n=[w for w in cigar_pattern_ind if q in w[1]]
        match_cood[n][1]=part_n[0][0]
        match_cood[n][3]=part_n[-1][0]
    return match_cood


def caculate_cigar_features(Matches, seq_length, match_cood, no_align):  
    cigar_pattern=[p.union(s) for p, s in zip(list(Matches[HTSeq.GenomicInterval("q", 0, seq_length, "+")]), reversed(list(Matches[HTSeq.GenomicInterval("q", 0, seq_length, "-")])))]   
    cigar_pattern_counts=[[_, sum(1 for _ in group)] for _, group in groupby(cigar_pattern)] 
    overlap = sum([i[1] for i in cigar_pattern_counts if len(i[0])>1])
    gap = sum([i[1] for i in cigar_pattern_counts[1:-1] if len(i[0])==0])
    match_cood0=correct_match_cood(match_cood, cigar_pattern, no_align)
    return match_cood0, cigar_pattern, cigar_pattern_counts, overlap, gap


def query_orientation_check(bundle, gap_max, overlap_max, min_exclusive):
    """check overlapping bp and no. of gaps between matches of alignment"""
    """correct for h-clipped base that has 0 cooridinate, messing up the orientations"""
    """ignore gap at ends"""
    """filter reads if they are not exclusive >= 25 bp"""
    no_align, Matches, match_cood, seq_length=create_cigar_object(bundle)
    match_cood, cigar_pattern, cigar_pattern_counts, overlap, gap=caculate_cigar_features(Matches, seq_length, match_cood, no_align)
###print cigar_pattern_counts
###unite + and - strand###
    if gap < gap_max and overlap < overlap_max and no_align >=2:
        left_exclusive, right_exclusive=calculate_exclusive(cigar_pattern, match_cood, cigar_pattern_counts, no_align)
        if left_exclusive>=min_exclusive and right_exclusive>=min_exclusive:
#### add in overlap info ###        
### correct for query partition for - strand ###
            for m in match_cood:
                m[0]=tuple(m[0])
                m[2]=tuple(m[2])
            match_cood_ind=[int(list(ind[0])[0]) for ind in cigar_pattern_counts if len(ind[0])==1]
            match_cood=[match_cood[ind] for ind in match_cood_ind] 
            match_cood=[match_cood[0], match_cood[-1]]
    else:
        ###print "fail gap and overlap filter"
        return (False, None)         
    if gap < gap_max and overlap < overlap_max and left_exclusive>=min_exclusive and right_exclusive>=min_exclusive:
        return (True, match_cood)
    else:
        return (False, None) 

def build_SA_bundles(in_sam):
    SA_bundles={}
    for a in HTSeq.bundle_multiple_alignments(in_sam):
        j=[]
        for aln in a:
            match_base=0
            for cigar in aln.cigar:
                if cigar.type=="M":
                    match_base+=cigar.size
            if match_base>=M_length and aln.aQual >=mapq:
                j.append(aln)
        if len(j) > 1:
            query_pass, match_cood=query_orientation_check(j, gap_max, overlap_max, min_exclusive)            
            if query_pass:
                SA_bundles[j[0].read.name]=[j, match_cood]
    return SA_bundles

### test code ###
###bundle=[a for a in HTSeq.bundle_multiple_alignments(in_sam) if a[0].read.name=="ZHENGLABSEQ002_1:11107:19590:10337"][0]
###ALK_keys= [i for i in SA_bundles.keys() if len(SA_bundles[i][1])>2 and SA_bundles[i][0][0].iv.chrom == "2"]
###SA_bundles['ZHENGLABSEQ002_1:11111:16408:20346']

####3. define, sort and cluster split points ###
### report readID, left split, right split ###
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



###3.2 filter for stagerring ends ###
def filter_stagerring_ends(splits_clustered, SA_bundles, orientation):
    """remove grouped reads if no. staggering < 3"""
    """relax filter to include breakpoints that had fused to different part of chromosome"""
    filtered_list=[]
    for i in splits_clustered:
        reads=[]
        bk_list=[]
        for j in i:
            reads.extend(list(j[1]))
        left_stagger=[]
        right_stagger=[]
        left_split=[]
        right_split=[]
        for r in reads:
            left_stagger.append(SA_bundles[r][1][0][0])
            right_stagger.append(SA_bundles[r][1][-1][2])
            left_split.append(SA_bundles[r][1][0][2])
            right_split.append(SA_bundles[r][1][-1][0])
        if orientation == "left":
### so left_stagger > 3 pass first filter layer ###
            stag_support=len(set(left_stagger))
### then if right_split == 1 and right_stagger ==1 ### 
            if stag_support >=3:
                if len(set(right_split))==1 and len(set(right_stagger))==1:
                    for j in i:
                        bk_list.append((j[0], stag_support))
                elif len(set(right_split))>1:
### multiple breakpoints so relax the filter for right stagger ###
                    for j in i:
                        bk_list.append((j[0], stag_support))
        if orientation == "right":
### so right_stagger > 3 and left_stagger ==1###
            stag_support=len(set(right_stagger))
            if stag_support >=3:
                if len(set(left_split))==1 and len(set(left_stagger))==1:
                    for j in i:
                        bk_list.append((j[0], stag_support))
                elif len(set(left_split))>1:
### multiple breakpoints so relax the filter for left stagger ###
                    for j in i:
                        bk_list.append((j[0], stag_support))
        if len(bk_list)>0:
            filtered_list.append(bk_list)
    return filtered_list




#### Refseq version###
###GFF_file='/home/athena/Refseq/HG19/NCBIRef_seq_HG19_fest.txt'
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



####4.0.1 build model for gene transcript ###
def build_gene_model(g, GFF_dict):
    """return gene model of a gene"""
    """define with codon_no, and codon partition"""
    gene_model=HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    exon_no=int(GFF_dict[g]['exonCount'])
    exon_start=[int(j) for j in GFF_dict[g]['exonStarts'].split(",")[:exon_no]]
    exon_end=[int(j) for j in GFF_dict[g]['exonEnds'].split(",")[:exon_no]]
    ###print g
    if GFF_dict[g]['strand']=="-":
        start_codon=int(GFF_dict[g]['cdsEnd'])
        stop_codon=int(GFF_dict[g]['cdsStart'])
        exon_start=list(reversed(exon_start))
        exon_end=list(reversed(exon_end))
        start_exon=[[s, e] for s, e in zip(exon_start, exon_end) if s<start_codon and e>=start_codon][0]
        end_exon=[[s, e] for s, e in zip(exon_start, exon_end) if s<=stop_codon and e>stop_codon][0]
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
                cDNA_part+=1
                if codon_partition==3:
                    codon_n+=1
                    codon_partition=0
                in_name=str(exon_no)+"_"+str(codon_n)+"_"+str(cDNA_part)+"_"+str(codon_partition)
                codon_partition+=1
                gene_model[HTSeq.GenomicInterval(GFF_dict[g]['chrom'], location, location+1)]+=in_name
    else:
        for i in exons_cood:
            exon_no+=1
            for location in range(i[0], i[1]):
                cDNA_part+=1
                if codon_partition==3:
                    codon_n+=1
                    codon_partition=0
                in_name=str(exon_no)+"_"+str(codon_n)+"_"+str(cDNA_part)+"_"+str(codon_partition)
                codon_partition+=1
                gene_model[HTSeq.GenomicInterval(GFF_dict[g]['chrom'], location, location+1)]+=in_name
    return gene_model



###4 identify the both end of breakpoints from splits_clustered ###
def collect_reads_from_splits(a_split_list, breakpoints):     
    r=set([])  
    for GIO in a_split_list:
        GI=GIO[0]
        score=GIO[1]
        for j in breakpoints[GI].steps():
            r.update(j[1])
    return list(r)


def group_splits_reads(reads, SA_bundles):
    """summarize splits and leave an example per split for annotation """
    """assume same left-right_split, same overlaps and gaps number in read"""
    left_right_splits={}
    for r in reads:
        match_cood=SA_bundles[r][1]
        boundary="-".join(["_".join([str(x) for x in SA_bundles[r][1][0][2]]), "_".join([str(x) for x in SA_bundles[r][1][-1][0]])])
        left_right_splits[boundary]=match_cood
    return left_right_splits

def combine_left_right_clusteres(left_splits_clustered0, right_splits_clustered, SA_bundles, left_breakpoints, right_breakpoints):
    """remove redundancy of splits invovled in multiple fusions"""
    left_right_splits={}
    for i in left_splits_clustered0:
        reads=collect_reads_from_splits(i, left_breakpoints)
        left_right_splits.update(group_splits_reads(reads, SA_bundles))
    for r in right_splits_clustered0:
        reads=collect_reads_from_splits(r, right_breakpoints)
        left_right_splits.update(group_splits_reads(reads, SA_bundles))
    return left_right_splits


###4.1 Annotate breakpoint### 
def extract_break_genes(breaks, transcript):
    """find genes where breaks happen."""
    """build gene model"""
    """if break in exon, report exon number and codon position"""
    """if break in intron, report closest intron"""
    break5=breaks[0][2]
    break3=breaks[1][0]
    gene5=list(transcript[HTSeq.GenomicPosition(break5[0], int(break5[1]), ".")])
    gene3=list(transcript[HTSeq.GenomicPosition(break3[0], int(break3[1]), ".")])
    return HTSeq.GenomicPosition(break5[0], int(break5[1]), "."), HTSeq.GenomicPosition(break3[0], int(break3[1]), "."), gene5, gene3

 
def determine_intron_5(break_5, g, genes_model, breaks):
    ### break_5 very likely at intron location, find the upstream and the downstream exons ###
        upstream5=[i for i in genes_model[g[0]].steps() if i[0].end < break_5.start]
        if upstream5==[]:
            bk5_anno="5_UTR"
            dis_from_5anno=0
        downstream5=[i for i in genes_model[g[0]].steps() if i[0].start > break_5.start]
        if downstream5 ==[]:
            bk5_anno="3_UTR"
            dis_from_5anno=0
        if upstream5!=[] and downstream5 !=[]:
            upstream5=upstream5[-1]
            downstream5=downstream5[0]
            downstream5_dis=downstream5[0].start - break_5.start          
            upstream5_dis=break_5.start - upstream5[0].start
        ### pick annotation closetest to the break### 
        ### in relative to the orientatino of the query sequence ###
            if breaks[0][4]=="+":
            ###dis from anno, anno upstream of bk == +, anno downstream of bk == - ###
                if upstream5_dis < downstream3_dis:
                    bk5_anno=list(upstream5[1])[0]
                    dis_from_5anno=upstream5_dis 
                else:  
                    bk5_anno=list(downstream5[1])[0]        
                    dis_from_5anno=downstream5_dis * -1
            else:
            ###dis from anno, anno upstream of bk == -, anno downstream of bk == + ###
                if upstream5_dis < downstream3_dis:
                    bk5_anno=list(upstream5[1])[0]
                    dis_from_5anno=upstream5_dis * -1
                else:  
                    bk5_anno=list(downstream5[1])[0]        
                    dis_from_5anno=downstream5_dis 
        return bk5_anno, dis_from_5anno, upstream5, downstream5

def determine_intron_3(break_3, g, genes_model, breaks):
### break_3 is very likely at intron location, find the upstream and the downstream exons ###
        upstream3=[i for i in genes_model[g[1]].steps() if i[0].end < break_3.start]
        if upstream3==[]:
            bk3_anno="5_UTR"
            dis_from_3anno=0
        downstream3=[i for i in genes_model[g[1]].steps() if i[0].start > break_3.start]
        if downstream3 ==[]:
            bk3_anno="3_UTR"
            dis_from_3anno=0
        if upstream3!=[] and downstream3 !=[]:
            upstream3=upstream3[-1]
            downstream3=downstream3[0]
            downstream3_dis=downstream3[0].start - break_3.start          
            upstream3_dis=break_3.start - upstream3[0].start
        ### pick annotation closetest to the break###
        ### in relative to the orientatino of the query sequence ###
            if breaks[1][4]=="+":
            ###dis from anno, anno upstream of bk == -, anno downstream of bk == +###
                if upstream3_dis < downstream3_dis:
                    bk3_anno=list(upstream3[1])[0]
                    dis_from_3anno=upstream3_dis * -1
                else:  
                    bk3_anno=list(downstream3[1])[0]        
                    dis_from_3anno=downstream3_dis
            else:
                if upstream5_dis < downstream3_dis:
                    bk3_anno=list(upstream3[1])[0]
                    dis_from_3anno=upstream3_dis 
                else:  
                    bk3_anno=list(downstream3[1])[0]        
                    dis_from_3anno=downstream3_dis * -1
### codon_gaps in how to calculate frame_shift ### 
        return bk3_anno, dis_from_3anno, upstream5, downstream3


def name_intron(bk5_ref_upstream5, bk5_ref_downstream3, bk5_gene_orientation):
    if bk5_gene_orientation=="+":
        intron_name="intron_" + list(bk5_ref_upstream5[1])[0].split("_")[0]
    else:
        intron_name="intron_" + list(bk5_ref_downstream3[1])[0].split("_")[0]
    return intron_name


def determine_frame(bk5_anno, dis_from_5anno, bk3_anno, dis_from_3anno, bk5_ref_upstream5, bk5_ref_downstream3, bk3_ref_upstream5, bk3_ref_downstream3, breaks, bk5_gene_orientation, bk3_gene_orientation):
    """put the annotations and query together to determine the mutation"""
###bp to consider + == gap, - == overlap ###
    bps=breaks[1][1]-breaks[0][3]
    bk_5_length=breaks[0][3]-breaks[0][1]
    bk_3_length=breaks[1][3]-breaks[1][1] 
    bk5_frame_status=None
    bk3_frame_status=None
### report the bk as intronic and abondon frame-shift calculation if dis > alignment length or gaps/overlap ###
    if dis_from_5anno > bk_5_length:
### annotation from break_5 is out of the read, so in theory break_5 is in intron###
        bk5_frame_status="intron"
    elif dis_from_5anno <0 and abs(dis_from_5anno) > abs(bps):
        bk5_frame_status="intron"
    if dis_from_3anno <0 and abs(dis_from_3anno) > abs(bps):
        bk3_frame_status="intron"
    elif dis_from_3anno > bk_3_length:
        bk3_frame_status="intron"
    if  bk5_frame_status=="intron":
        bk5_frame_status=name_intron(bk5_ref_upstream5, bk5_ref_downstream3, bk5_gene_orientation)
    if  bk3_frame_status=="intron":
        bk3_frame_status=name_intron(bk3_ref_upstream5, bk3_ref_downstream3, bk3_gene_orientation)
### calculate frame status if both bkpoints are not intronic ###
    if  "UTR" not in bk5_anno and "UTR" not in bk3_anno:
        if bk5_frame_status==None:
            bk5_frame_status="exon" + bk5_anno.split("_")[0]
        if bk3_frame_status==None:
            bk3_frame_status="exon" + bk3_anno.split("_")[0]
    if  "exon" in bk5_frame_status and "exon" in bk3_frame_status:
#### bk very close to exons of in exons and require frame determination ###
        bk5_cDNA_part=int(bk5_anno.split("_")[-1])
        bk3_cDNA_part=int(bk3_anno.split("_")[-1])      
        gaps=dis_from_5anno +  bps  +  dis_from_3anno
### what is in frame???
        if (bk5_cDNA_part + bk3_cDNA_part+ gaps)%3 == 0:
                frame_status="inframe"
        else:
                frame_status="out_frame"
    else:
        frame_status="NAN"
    return bk5_frame_status, bk3_frame_status, frame_status



def determine_break_location(break_5, break_3, gene5, gene3, genes_model, GFF_dict, breaks):
    """identify where break position contain"""
    """return break5-break3 annotation"""
    """update gens_model"""
### what if there is a list of gene? ###
    splits_combinations=[]
    for g in list(itertools.product(gene5, gene3)):
        print g
        try:
            bk5_location=genes_model[g[0]][break_5]
        except:
            genes_model[g[0]]=build_gene_model(g[0], GFF_dict)
        try:
            bk3_location=genes_model[g[1]][break_3]
        except:
            genes_model[g[1]]=build_gene_model(g[1], GFF_dict)
        bk5_location=genes_model[g[0]][break_5]
        bk3_location=genes_model[g[1]][break_3]
### genes Name ###
        bk5_gene=GFF_dict[g[0]]['name2']
        bk3_gene=GFF_dict[g[1]]['name2']
        bk5_gene_orientation=GFF_dict[g[0]]['strand']
        bk3_gene_orientation=GFF_dict[g[1]]['strand']
        if bk5_location==set([]):
            bk5_anno, dis_from_5anno, bk5_ref_upstream5, bk5_ref_downstream3=determine_intron_5(break_5, g, genes_model, breaks)
        else:
            bk5_anno=list(bk5_location)[0]
            dis_from_5anno=0
            bk5_ref_upstream5=None
            bk5_ref_downstream3=None
        if bk3_location==set([]):
            bk3_anno, dis_from_3anno, bk3_ref_upstream5, bk3_ref_downstream3=determine_intron_3(break_3, g, genes_model, breaks)
        else:
            bk3_anno=list(bk3_location)[0]
            dis_from_3anno=0
            bk3_ref_upstream5=None
            bk3_ref_downstream3=None
### claculate frames for fusion###
        bk5_frame_status, bk3_frame_status, frame_status=determine_frame(bk5_anno, dis_from_5anno, bk3_anno, dis_from_3anno, bk5_ref_upstream5, bk5_ref_downstream3, bk3_ref_upstream5, bk3_ref_downstream3, breaks, bk5_gene_orientation, bk3_gene_orientation)
        splits_combinations.append([bk5_gene, bk5_anno, bk3_gene, bk3_anno, frame_status])
    return splits_combinations, genes_model


def extract_boundary_from_reads(GFF_dict, transcript, left_right_splits, genes_model, target_genes):
    """export gene-exon cood for each split point"""
    """test if the split points are in-frame"""
    breaks_output=[]
    genes_model={}
    for i in left_right_splits.keys():  
        breaks=left_right_splits[i]
        break_5, break_3, gene5, gene3=extract_break_genes(breaks, transcript)
### Only annotate target_genes if available###
###
###
### annotate break5 and break3 together and down size gene_model if come from the same gene ###
        splits_combinations, gene_model0=determine_break_location(break_5, break_3, gene5, gene3, genes_model, GFF_dict, breaks)
        genes_model.update(gene_model0)        
        break_output.extend(splits_combinations)
    return break_output



.append("\t".join([i, bk5_gene, bk5_anno, bk3_gene, bk3_anno, frame_status]))

            

### execute code if called ###
if __name__ == "__main__":
    import os, sys, subprocess, HTSeq, argparse, itertools
    from subprocess import Popen, PIPE
    from itertools import groupby 
    parser = argparse.ArgumentParser(description='SplitFusion v1')   
    parser = argparse.ArgumentParser(description='SplitFusion: identifying chromosomal rearrangement via target-deep sequencing RNA/DNA data')
    parser.add_argument('--wkdir', help='dir for intermediate and final output file(s).', action="store", dest="wkdir")
    parser.add_argument('-i', help='full path for the input sam/bam file for fusion detection', action="store", dest="in_file")
    parser.add_argument('--ref', help='full path for the reference genomic feature file', action="store", dest="GFF_file")
    parser.add_argument('-o', help='prefix of the output file, default=Split_fusion.out', action="store", dest="out_txt", default='Split_fusion.out')
    parser.add_argument('--gapmax', help='Maximum number of gaps allowed in the alignments of split read', type=int, action="store", dest="gap_max", default=2)
    parser.add_argument('--overlapmax', help='Maximum number of overlaps allowed in the alignments of split read', type=int, action="store", dest="overlap_max", default=10)
    parser.add_argument('--mapq', help='Minimum mapping quality for an alignment used for analysis', type=int, action="store", dest="mapq", default=20)
    parser.add_argument('--minexclusive', help='Minimum exclusive bases allowed among alignments in a split-read', type=int, action="store", dest="min_exclusive", default=25)
    parser.add_argument('--mlength', help='Minimum number of bases matches in an alignment to be used in the analysis', type=int, action="store", dest="M_length", default=25)
    args = parser.parse_args()
    wkdir=args.wkdir
    infile=args.in_file
    in_samSA=[i for i in infile.split("/") if "bam" in i or "sam" in i][0]
    in_samSA=in_samSA.replace("bam", "sam")
    if "bam" in infile:
        Grep_SA="samtools view " + infile  + " | grep \'SA:\' > " + wkdir + "/" + in_samSA
        Grep_=Popen(Grep_SA, shell=True)
        Grep_.communicate()
        Sort_SA="sort -k1 " + wkdir + "/" +  in_samSA + " > " +  wkdir + "/" + in_samSA + ".sorted" 
        Sort_=Popen(Sort_SA, shell=True)
        Sort_.communicate()   
    elif "sam" in in_bam:
        Sort_SA="grep \'SA:\'" + infile + "| sort -k1 > " +  wkdir + "/" + in_samSA + ".sorted" 
        Sort_=Popen(Sort_SA, shell=True)
        Sort_.communicate()   
    else:
        print "please input bam file"
        sys.exit   
    ###Set parameters for filtering split-reads###
    in_sam=HTSeq.SAM_Reader(wkdir + "/" + in_samSA + ".sorted")
    gap_max=args.gap_max
    overlap_max=args.overlap_max
    mapq=args.mapq
    M_length=args.M_length
    min_exclusive=args.M_length
    ### only save alignment of a read if > mapq###
    SA_bundles=build_SA_bundles(in_sam)
    left_breakpoints=HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    right_breakpoints=HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    for i in SA_bundles.keys():
        left_breakpoints[HTSeq.GenomicInterval(SA_bundles[i][1][0][2][0], SA_bundles[i][1][0][2][1], SA_bundles[i][1][0][2][1]+1, strand='.')]+=i
        right_breakpoints[HTSeq.GenomicInterval(SA_bundles[i][1][-1][0][0], SA_bundles[i][1][-1][0][1], SA_bundles[i][1][-1][0][1]+1,strand='.')]+=i
    left_splits_clustered=group_splits(left_breakpoints)
    right_splits_clustered=group_splits(right_breakpoints)
    left_splits_clustered0=filter_stagerring_ends(left_splits_clustered, SA_bundles, "left")
    right_splits_clustered0=filter_stagerring_ends(right_splits_clustered, SA_bundles, "right")
    left_right_splits=combine_left_right_clusteres(left_splits_clustered0, right_splits_clustered, SA_bundles, left_breakpoints, right_breakpoints)

    GFF_dict, transcript=read_in_feature_file(args.GFF_file)
    genes_model={}
###
    

    
    

    bk_output= extract_boundary_from_reads(GFF_dict, transcript, left_right_splits, genes_model, target_genes)

    outfile=wkdir + "/" + args.out_txt
    o=open(outfile, "w")
    o.write(left_output)
    o.write(right_output)
    o.close()


### test_code ###
### python 2018_05_15_SplitFusion_AC_MA.py \
###-i /home/athena/Splitfusion/SplitFusion_ZLZ/example-run-folder/A01-P701.consolidated.bam \
###-o A01-P701.split_Fusion.txt \
###--ref /home/athena/Refseq/HG19/NCBIRef_seq_HG19_fest.txt \
###--wkdir /home/athena/Splitfusion/SplitFusion_ZLZ/example-run-folder

###python 2018_05_15_SplitFusion_AC_MA.py \
###-i /home/athena/Splitfusion/SplitFusion_ZLZ/example-run-folder/A02-P702.consolidated.bam \
###-o 2018_05_15_A02-P702.split_Fusion.txt \
###--ref /home/athena/Refseq/HG19/NCBIRef_seq_HG19_fest.txt \
###--wkdir /home/athena/Splitfusion/SplitFusion_ZLZ/example-run-folder



