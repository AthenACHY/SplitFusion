import HTSeq, itertools
from itertools import groupby
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
    """remove split_location_that fused to more than 20 different places -- potential artefacts???"""
    left_right_splits={}
    for r in reads:
        match_cood=SA_bundles[r][1]
        boundary="-".join(["_".join([str(x) for x in SA_bundles[r][1][0][2]]), "_".join([str(x) for x in SA_bundles[r][1][-1][0]])])
        left_right_splits[boundary]=match_cood
    if len(left_right_splits)<=20:
        return left_right_splits
    else:
        return {}

def combine_left_right_clusteres(left_splits_clustered0, right_splits_clustered0, SA_bundles, left_breakpoints, right_breakpoints):
    """remove redundancy of splits invovled in multiple fusions"""
    left_right_splits={}
    for i in left_splits_clustered0:
        reads=collect_reads_from_splits(i, left_breakpoints)
        left_right_splits.update(group_splits_reads(reads, SA_bundles))
    for r in right_splits_clustered0:
        reads=collect_reads_from_splits(r, right_breakpoints)
        left_right_splits.update(group_splits_reads(reads, SA_bundles))
    return left_right_splits

###3+4 wrapper###
def group_and_filter_splits(left_breakpoints, right_breakpoints, SA_bundles):
    """Filter splits according to uniqmapping at staggering ends"""
    """merge left and right splits together for annotations"""
    left_splits_clustered=group_splits(left_breakpoints)
    right_splits_clustered=group_splits(right_breakpoints)
    left_splits_clustered0=filter_stagerring_ends(left_splits_clustered, SA_bundles, "left")
    right_splits_clustered0=filter_stagerring_ends(right_splits_clustered, SA_bundles, "right")
    left_right_splits=combine_left_right_clusteres(left_splits_clustered0, right_splits_clustered0, SA_bundles, left_breakpoints, right_breakpoints)
    return left_splits_clustered0, right_splits_clustered0, left_right_splits


