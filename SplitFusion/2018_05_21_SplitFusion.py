if __name__ == "__main__":
    import os, sys, subprocess, HTSeq, argparse, itertools
    from subprocess import Popen, PIPE
    from itertools import groupby 
    from Split_Fusion import Multiple_alignment, Group_split, Fusion_annotation
    parser = argparse.ArgumentParser(description='SplitFusion v1')   
    parser = argparse.ArgumentParser(description='SplitFusion: identifying chromosomal rearrangement via target-deep sequencing RNA/DNA data')
    parser.add_argument('--wkdir', help='dir for intermediate and final output file(s).', action="store", dest="wkdir")
    parser.add_argument('-i', help='full path for the input sam/bam file for fusion detection', action="store", dest="in_file")
    parser.add_argument('--targets', help='full path for targets_genes', action="store", dest="target_genes", default=None)
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
    SA_bundles, left_breakpoints, right_breakpoints=Multiple_alignment.build_alignment_bundles(in_sam)
`   left_splits_clustered0, right_splits_clustered0, left_right_splits=Group_split.group_and_filter_splits(left_breakpoints, right_breakpoints, SA_bundles)
    GFF_dict, transcript=Fusion_annotation.read_in_feature_file(args.GFF_file)
    genes_model={}
    bk_output=Fusion_annotation.extract_boundary_from_reads(GFF_dict, transcript, left_right_splits, genes_model, args.target_genes)
### output all posible combinations of the transcript models for each breakpoints###
    bk_output=list(set(bk_output)
    bk_output="\n".join(bk_output)
    outfile=wkdir + "/" + args.out_txt
    o=open(outfile, "w")
    o.write(bk_output)
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



