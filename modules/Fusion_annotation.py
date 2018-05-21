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
        pre_exon=len([[s, e] for s, e in zip(exon_start, exon_end) if e>=start_codon])
        end_exon=[[s, e] for s, e in zip(exon_start, exon_end) if s<=stop_codon and e>stop_codon][0]
    else:
        start_codon=int(GFF_dict[g]['cdsStart'])
        stop_codon=int(GFF_dict[g]['cdsEnd'])
        start_exon=[[s, e] for s, e in zip(exon_start, exon_end) if s<=start_codon and e>start_codon][0]
        pre_exon=len([[s, e] for s, e in zip(exon_start, exon_end) if s<=start_codon])
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
    exon_no=pre_exon-1  
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
                if upstream5_dis < downstream5_dis:
                    bk5_anno=list(upstream5[1])[0]
                    dis_from_5anno=upstream5_dis 
                else:  
                    bk5_anno=list(downstream5[1])[0]        
                    dis_from_5anno=downstream5_dis * -1
            else:
            ###dis from anno, anno upstream of bk == -, anno downstream of bk == + ###
                if upstream5_dis < downstream5_dis:
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
                if upstream3_dis < downstream3_dis:
                    bk3_anno=list(upstream3[1])[0]
                    dis_from_3anno=upstream3_dis 
                else:  
                    bk3_anno=list(downstream3[1])[0]        
                    dis_from_3anno=downstream3_dis * -1
### codon_gaps in how to calculate frame_shift ### 
        return bk3_anno, dis_from_3anno, upstream3, downstream3


def name_intron(bk5_ref_upstream5, bk5_ref_downstream3, bk5_gene_orientation):
    if bk5_gene_orientation=="+":
        intron_name="intron_" + list(bk5_ref_upstream5[1])[0].split("_")[0]
    else:
        intron_name="intron_" + list(bk5_ref_downstream3[1])[0].split("_")[0]
    return intron_name




###4.1 Annotate breakpoint### 
def determine_frame(bk5_anno, dis_from_5anno, bk3_anno, dis_from_3anno, bk5_ref_upstream5, bk5_ref_downstream3, bk3_ref_upstream5, bk3_ref_downstream3, breaks, bk5_gene_orientation, bk3_gene_orientation):
    """put the annotations and query together to determine the mutation"""
###bp to consider + == gap, - == overlap ###
###    print bk5_anno, dis_from_5anno, bk3_anno, dis_from_3anno, breaks, bk5_gene_orientation, bk3_gene_orientation
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
    else:
        bk5_frame_status="UTR"
        bk3_frame_status="UTR"
    if  "exon" in bk5_frame_status and "exon" in bk3_frame_status:
#### bk very close to exons of in exons and require frame determination ###
        bk5_cDNA_part=int(bk5_anno.split("_")[-2]) 
        bk3_codon_phrase=int(bk3_anno.split("_")[-1])           
        gaps=dis_from_5anno +  bps  +  dis_from_3anno
### what is in frame???
        if (bk5_cDNA_part + gaps)%3 == bk3_codon_phrase:
                frame_status="inframe"
        else:
                frame_status="out_frame"
    else:
        frame_status="NAN"
    return bk5_frame_status, bk3_frame_status, frame_status



def determine_break_location(break_5, break_3, gene5, gene3, genes_model, GFF_dict, breaks, gene_list, a_split):
    """identify where break position contain"""
    """return break5-break3 annotation"""
    """update gens_model"""
### what if there is a list of gene? ###
    splits_combinations=[]
    for g in gene_list:
###        print g
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
        splits_combinations.append("\t".join([a_split, bk5_gene, bk5_anno, bk3_gene, bk3_anno, frame_status])+"\n")
    return splits_combinations, genes_model





def extract_genes_targets(target_genes):
    targets=[]
    with open(target_genes, "r") as t:
        for row in t:
            targets.append(row.replace("\n", ""))
    return targets

def extract_boundary_from_reads(GFF_dict, transcript, left_right_splits, genes_model, target_genes):
    """export gene-exon cood for each split point"""
    """test if the split points are in-frame"""
    if target_genes !=None:
        targets=extract_genes_targets(target_genes)
    breaks_output=[]
    genes_model={}
    for a_split in left_right_splits.keys():  
        breaks=left_right_splits[a_split]
        break_5, break_3, gene5, gene3=extract_break_genes(breaks, transcript)
        gene_list=list(itertools.product(gene5, gene3))
        if target_genes !=None:
            gene_list=[g for g in gene_list if GFF_dict[g[0]]["name2"] in targets or GFF_dict[g[1]]["name2"] in targets]
            gene_list
        if len(gene_list) >0:
            splits_combinations, gene_model0=determine_break_location(break_5, break_3, gene5, gene3, genes_model, GFF_dict, breaks, gene_list, a_split)
            genes_model.update(gene_model0)        
            breaks_output.extend(splits_combinations)
    return breaks_output



