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


def correct_match_cood(match_cood, cigar_pattern, no_align,seq_length):
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
    match_cood0=correct_match_cood(match_cood, cigar_pattern, no_align, seq_length)
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

def build_alignment_bundles(in_sam):
    """import in the name-sorted sam file as SA_bundles"""
    """create chormosomal_feature_arrays of the left/right breakpoints"""
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
    left_breakpoints=HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    right_breakpoints=HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    for i in SA_bundles.keys():
        left_breakpoints[HTSeq.GenomicInterval(SA_bundles[i][1][0][2][0], SA_bundles[i][1][0][2][1], SA_bundles[i][1][0][2][1]+1, strand='.')]+=i
        right_breakpoints[HTSeq.GenomicInterval(SA_bundles[i][1][-1][0][0], SA_bundles[i][1][-1][0][1], SA_bundles[i][1][-1][0][1]+1,strand='.')]+=i
    return SA_bundles, left_breakpoints, right_breakpoints

