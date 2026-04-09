
import sys
import os


import pysam
import bisect



def refine_matched_coords(matching_ref_coords, start_pos):

	cur_coord = start_pos -1
	for i in range(len(matching_ref_coords)):
		if matching_ref_coords[i] == None:
			matching_ref_coords[i] = cur_coord
		else:
			cur_coord = matching_ref_coords[i] 

def consensus_localization(consensus_bam, chrom, coord1, coord2):
	samfile = pysam.AlignmentFile(consensus_bam)
	sam_iter = samfile.fetch(chrom, coord1, coord2)
	read = next(sam_iter)
	#print(read)

	matching_ref_coords = read.get_reference_positions(full_length=True)
	refine_matched_coords(matching_ref_coords, read.pos)
	#print(read.reference_name)
	read_pos1 = bisect.bisect_left(matching_ref_coords, coord1)
	read_pos2 = bisect.bisect_left(matching_ref_coords, coord2)
	return(read.query_sequence, read_pos1, read_pos2)


def viewconsensus_in_region(consensus_bam, chrom, begin, end):
	fullseq, read_pos1, read_pos2 = consensus_localization(consensus_bam, chrom, begin, end)
	name = "_".join(map(str, ["consensus", begin, end]))

	print(">"+name+"\n"+fullseq[read_pos1:read_pos2])



if __name__ == '__main__':
	
	operation = sys.argv[1]
	if operation == "CONSENSUS_REGION":
		consensus_bam = sys.argv[2]
		chrom = sys.argv[3]
		begin = int(sys.argv[4])
		end = int(sys.argv[5])
		viewconsensus_in_region(consensus_bam, chrom, begin, end)
