#!/bin/bash
set -e
set -x
set -u 

PN=$1 # NA21309
HAP_1_FA_GZ=$2 # /odinn/tmp/dorukb/scratch/MSA/review.NG/replication/AllofUs/sample_assembly_map/HPRC_test/NA21309.paternal.f1_assembly_v2.fa.gz
HAP_2_FA_GZ=$3 # /odinn/tmp/dorukb/scratch/MSA/review.NG/replication/AllofUs/sample_assembly_map/HPRC_test/NA21309.maternal.f1_assembly_v2.fa.gz
REFERENCE=$4 # /nfs/odinn/tmp/dorukb/scratch/MSA/review.NG/replication/AllofUs/references/chr3.fa
STOP_AT="$5"

OUT_M1_COUNT_FN="${PN}.motif_counts.csv" # outputdir/NA21309.motif_counts.csv
OUT_LENGTH_FN="${PN}.lengths.csv" # outputdir/NA21309.lengths.csv


# progs
SAMTOOLS=samtools
MINIMAP=minimap2


# python and R progs
MSA_VIEW_PY=/usr/local/bin/view_region.py

# vars
CHROM_REGION="chr3"
BEGIN_REGION=169770169
END_REGION=169771552
M1="ATATATAAAAATTTATATTTATATATC"
M1_REVCOMP="GATATATAAATATAAATTTTTATATAT"


check_stop() {
  local name="$1"
  if [[ "$STOP_AT" == "$name" ]]; then
    echo "Stopping at step: $name"
    exit 0
  fi
}


# Get both alleles.
# Allele 1
python ${MSA_VIEW_PY} CONCAT_FASTA_LINES ${HAP_1_FA_GZ} | grep --no-group-separator  -B 1 -P "${M1}|${M1_REVCOMP}" > ${HAP_1_FA_GZ}.tomap.fa
python ${MSA_VIEW_PY} CONCAT_FASTA_LINES ${HAP_2_FA_GZ} | grep --no-group-separator  -B 1 -P "${M1}|${M1_REVCOMP}" > ${HAP_2_FA_GZ}.tomap.fa

check_stop "unzip"

if [[ -s "${HAP_1_FA_GZ}.tomap.fa" && -s "${HAP_2_FA_GZ}.tomap.fa" ]]; then

	${MINIMAP} -a -x asm5 --cs -t1 -z 3000,1500 ${REFERENCE} ${HAP_1_FA_GZ}.tomap.fa > ${HAP_1_FA_GZ}.tomap.unsorted.bam 
	${SAMTOOLS} sort --threads 1 ${HAP_1_FA_GZ}.tomap.unsorted.bam >  ${HAP_1_FA_GZ}.tomap.bam
	${SAMTOOLS} index ${HAP_1_FA_GZ}.tomap.bam
	python ${MSA_VIEW_PY} CONSENSUS_REGION ${HAP_1_FA_GZ}.tomap.bam  ${CHROM_REGION} ${BEGIN_REGION} ${END_REGION} > ${HAP_1_FA_GZ}.tomap.bam.allele.fa

	check_stop "allele1"

	# Allele 2
	${MINIMAP} -a -x asm5 --cs -t1 -z 3000,1500 ${REFERENCE} ${HAP_2_FA_GZ}.tomap.fa | ${SAMTOOLS} sort --threads 1 >  ${HAP_2_FA_GZ}.tomap.bam
	${SAMTOOLS} index ${HAP_2_FA_GZ}.tomap.bam
	python ${MSA_VIEW_PY} CONSENSUS_REGION ${HAP_2_FA_GZ}.tomap.bam  ${CHROM_REGION} ${BEGIN_REGION} ${END_REGION} > ${HAP_2_FA_GZ}.tomap.bam.allele.fa


	# Now get info on the alleles.
	HAP1_M1_COUNT=$(cat ${HAP_1_FA_GZ}.tomap.bam.allele.fa | awk '$0 !~ ">" ' |  awk -F"${M1}" '{print  (NF-1)  }' ) 
	HAP1_LENGTH=$(cat ${HAP_1_FA_GZ}.tomap.bam.allele.fa | awk '$0 !~ ">" ' |  awk  '{print  length($1) }' )

	HAP2_M1_COUNT=$(cat ${HAP_2_FA_GZ}.tomap.bam.allele.fa | awk '$0 !~ ">" ' |  awk -F"${M1}" '{print  (NF-1)  }' ) 
	HAP2_LENGTH=$(cat ${HAP_2_FA_GZ}.tomap.bam.allele.fa | awk '$0 !~ ">" ' |  awk  '{print  length($1) }' )

	# Final outputs.
	echo -e "${PN}\t${HAP1_M1_COUNT}\t${HAP2_M1_COUNT}" > ${OUT_M1_COUNT_FN}
	echo -e "${PN}\t${HAP1_LENGTH}\t${HAP2_LENGTH}" > ${OUT_LENGTH_FN}
else
	echo -e "${PN}\tNA\tNA" > ${OUT_M1_COUNT_FN}
	echo -e "${PN}\tNA\tNA" > ${OUT_LENGTH_FN}
fi;
