### protospacer template ###
### read in coordinate ###
### cut out regions from reference genome ###
### align to 100bp and 200bp encompassing region to find a protospacer like sequence ###

import string
import operator
import pyfaidx
import re
import swalign
import sys

match = 2
mismatch = -1
scoring = swalign.NucleotideScoringMatrix(match, mismatch)
sw = swalign.LocalAlignment(scoring, gap_penalty=-100, gap_extension_penalty=-100, prefer_gap_runs=True)  

### load fasta and phrase reference ###
reference_genome="/home/a/Documents/Refseq/GRCh_38.fa"
genome = pyfaidx.Fasta(reference_genome)

#intab="/home/a/Documents/Edited_seq/20191120_brekpoints_information/Lib010.genomic.breakpoint.candidates.txt"
#outtab="/home/a/Documents/Edited_seq/20191120_brekpoints_information/Lib010.genomic.breakpoint.candidates_spacers.txt"
#protospecer_template="GACCCCCTCCACCCCGCCTC"
def read_in_breakpoint_AB(intab):
	breakpoints=[]
	f=open(intab, "r")
	for row in f:
		row=row.replace("\n", "")
		row=row.split("\t")
		try:
			row[2]=int(row[2])
			row[4]=int(row[4])
			breakpoints.append(row)
		except:
			continue
	return breakpoints

def assign_breakpoint_region(breakpoint, margin):
	chrA=breakpoint[1]
	start_A=breakpoint[2]-margin
	end_A=breakpoint[2]+margin
	if start_A<0:
		start_A=0
	chrB=breakpoint[3]
	start_B=breakpoint[4]-margin
	end_B=breakpoint[4]+margin
	if start_B<0:
		start_B=0
	return chrA, start_A, end_A, chrB, start_B, end_B

def hamming_distance2(chaine1, chaine2):
    return len(list(filter(lambda x : ord(x[0])^ord(x[1]), zip(chaine1, chaine2))))


def refine_alignment(ref_seq_rc, protospecer_template):
	'''compare protospacer_template with all possible seq with same size within the ref_seq_rc'''
	'''change numbers in case the protospacer length changes'''
	protospacer_len=len(protospecer_template)
	end=len(ref_seq_rc)-len(protospecer_template)
	ref_seq_rc_list=[ref_seq_rc[x:x+protospacer_len] for x in range(1, end, 1)]
	hamming_d=[hamming_distance2(i, protospecer_template) for i in ref_seq_rc_list]
	min_mismatches=min(hamming_d)
	min_mismatches_no=hamming_d.count(min_mismatches)
	align_list=zip(hamming_d, ref_seq_rc_list)
	align_list.sort()
	if min_mismatches_no >1:
		# prioritize sites with NGG #
		align_list=[j for j in align_list if j[0]==min_mismatches]
		off_targets_pam_location=[ref_seq_rc_list.index(j[1])+protospacer_len+1 for j in align_list ]
		off_targets_pams_3rd=[ref_seq_rc[j+2:j+3] for j in off_targets_pam_location]
		if off_targets_pams_3rd.count("G") ==1:
			best_protospacer=align_list[off_targets_pams_3rd.index("G")]
		else:
			off_targets_pams_2nd=[ref_seq_rc[j+1:j+2] for j in off_targets_pam_location]
			if off_targets_pams_2nd.count("G") ==1:	
				best_protospacer=align_list[off_targets_pams_2nd.index("G")]
			elif off_targets_pams_2nd.count("A") ==1 and off_targets_pams_2nd.count("G") ==0:		
				best_protospacer=align_list[off_targets_pams_2nd.index("A")]
			else:
				best_protospacer=align_list[0] 	
	else:
		best_protospacer=align_list[0]
	off_target_seq_start=ref_seq_rc_list.index(best_protospacer[1])+1 ### switch to 1-base genome ###
	off_target_spacer_seq=ref_seq_rc[off_target_seq_start:off_target_seq_start+protospacer_len+3]
	return best_protospacer[0], off_target_spacer_seq, off_target_seq_start


def find_protospacer(sw, ref_seq, ref_seq_rc, protospecer_template):
	aln=sw.align(protospecer_template, ref_seq)
	aln_rc=sw.align(protospecer_template, ref_seq_rc)
	if aln.matches > aln_rc.matches:
		orientation="+"
		mismatches, off_target_seq, off_target_seq_start=refine_alignment(ref_seq, protospecer_template)
	else:
		mismatches, off_target_seq, off_target_seq_start=refine_alignment(ref_seq_rc, protospecer_template)
		orientation="-"
	pam_seq=off_target_seq[-3:]
	return mismatches, off_target_seq, pam_seq, off_target_seq_start, orientation

def generate_ref_seq(chrA, start_A, end_A):
	ref_seq=str(genome[chrA][start_A:end_A].seq)
	ref_seq=ref_seq.upper()
	ref_seq_rc=str(genome[chrA][start_A:end_A].reverse.complement.seq)
	ref_seq_rc=ref_seq_rc.upper()
	return ref_seq, ref_seq_rc 

def find_seq_homology(sw, genome, breakpoint, protospecer_template, margin):
	chrA, start_A, end_A, chrB, start_B, end_B=assign_breakpoint_region(breakpoint, margin)
	try:
		ref_seq, ref_seq_rc=generate_ref_seq(chrA, start_A, end_A)
		mismatches, off_target_seq, pam_seq, off_target_seq_start, orientation=find_protospacer(sw, ref_seq, ref_seq_rc, protospecer_template)
		if orientation =="+":
			protospacer_alignment_start=start_A+off_target_seq_start
		else:
			protospacer_alignment_start=end_A - off_target_seq_start - len(protospecer_template)
		off_target_A="\t".join([str(mismatches), off_target_seq, pam_seq, orientation, str(protospacer_alignment_start)])
		ref_seq, ref_seq_rc=generate_ref_seq(chrB, start_B, end_B)
		mismatches, off_target_seq, pam_seq, off_target_seq_start, orientation=find_protospacer(sw, ref_seq, ref_seq_rc, protospecer_template)
		if orientation =="+":
			protospacer_alignment_start=start_B+off_target_seq_start
		else:
			protospacer_alignment_start=end_B - off_target_seq_start - len(protospecer_template)
		off_target_B="\t".join([str(mismatches), off_target_seq, pam_seq, orientation, str(protospacer_alignment_start)])
		return off_target_A + "\t" + off_target_B
	except:
		return "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"

def write_out_analysis(outtab, breakpoints, off_target_info_100, off_target_info_200):
	header="NumReads\tchrA\tpositionA\tchrB\tpositionB\tStructureVariation\tIntraChrDistance\tmismatches_A_100\toff_target_seq_A_100\tpam_seq_A_100\torientation\tprotospacer_startA_100\tmismatches_B_100\toff_target_seq_B_100\tpam_seq_B_100\torientation\tprotospacer_startB_100\tmismatches_A_200\toff_target_seq_A_200\tpam_seq_A_200\torientation\tprotospacer_startA_200\tmismatches_B_200\toff_target_seq_B_200\tpam_seq_B_200\torientation\tprotospacer_startB_200\n"
	output=zip(["\t".join([str(j) for j in b]) for b in breakpoints], off_target_info_100, off_target_info_200)
	output=["\t".join(t) for t in output]
	output="\n".join(output)
	o=open(outtab, "w")
	o.write(header)
	o.write(output)
	o.close()




if __name__ == '__main__':
	intab=sys.argv[1]
	print intab
	outtab=sys.argv[2]
	protospecer_template=sys.argv[3]
	breakpoints=read_in_breakpoint_AB(intab)
	off_target_info_100=[]
	off_target_info_200=[]
	for b in breakpoints:
		off_target_info_100.append(find_seq_homology(sw, genome, b, protospecer_template, 50))
	for b in breakpoints:
		off_target_info_200.append(find_seq_homology(sw, genome, b, protospecer_template, 100))
	write_out_analysis(outtab, breakpoints, off_target_info_100, off_target_info_200)

### command ###
# python2 /home/a/Documents/Edited_seq/2019_11_20_VEGFA/2019_11_28_find_spacer.py /home/a/Documents/Edited_seq/20191120_brekpoints_information/Lib010.genomic.breakpoint.candidates.txt /home/a/Documents/Edited_seq/20191120_brekpoints_information/Lib010.genomic.breakpoint.candidates_spacer.txt GACCCCCTCCACCCCGCCTC

# python2 /home/a/Documents/Edited_seq/2019_11_20_VEGFA/2019_11_28_find_spacer.py /home/a/Documents/Edited_seq/20191120_brekpoints_information/Lib011.genomic.breakpoint.candidates.txt /home/a/Documents/Edited_seq/20191120_brekpoints_information/Lib011.genomic.breakpoint.candidates_spacer.txt GACCCCCTCCACCCCGCCTC
# python2 /home/a/Documents/Edited_seq/2019_11_20_VEGFA/2019_11_28_find_spacer.py /home/a/Documents/Edited_seq/20191120_brekpoints_information/Lib012.genomic.breakpoint.candidates.txt /home/a/Documents/Edited_seq/20191120_brekpoints_information/Lib012.genomic.breakpoint.candidates_spacer.txt GACCCCCTCCACCCCGCCTC
# python2 /home/a/Documents/Edited_seq/2019_11_20_VEGFA/2019_11_28_find_spacer.py /home/a/Documents/Edited_seq/20191120_brekpoints_information/Lib013.genomic.breakpoint.candidates.txt /home/a/Documents/Edited_seq/20191120_brekpoints_information/Lib013.genomic.breakpoint.candidates_spacer.txt GACCCCCTCCACCCCGCCTC
# python2 /home/a/Documents/Edited_seq/2019_11_20_VEGFA/2019_11_28_find_spacer.py /home/a/Documents/Edited_seq/20191120_brekpoints_information/Lib014.genomic.breakpoint.candidates.txt /home/a/Documents/Edited_seq/20191120_brekpoints_information/Lib014.genomic.breakpoint.candidates_spacer.txt GACCCCCTCCACCCCGCCTC
