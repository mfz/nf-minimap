
# progs
SAMTOOLS=/nfs/fs1/bioinfo/apps-x86_64/samtools/1.21.10/bin/samtools
MINIMAP=/nfs/odinn/users/dorukb/programs/minimap2-2.28_x64-linux/minimap2


# python and R progs
MSA_VIEW_PY=/nfs/odinn/users/dorukb/manuscripts/tandem_repeats_cancer.NG/review_20260210/scripts/florian/view_region.py


# files and folders
TERC_REGION_REFERENCE_FA=/nfs/odinn/tmp/dorukb/scratch/MSA/review.NG/replication/AllofUs/references/chr3.fa

.ONESHELL:
SHELL := /bin/bash


%.fa.3q26.2-TR.PacBio_asm_allele.info.done: %.fa
	mkdir -p $$(dirname "$@")
	#set -e
	#set -euo pipefail
	@M1="ATATATAAAAATTTATATTTATATATC"
	@M1_REVCOMP="GATATATAAATATAAATTTTTATATAT"
	cat $(word 1,$^) |  awk '/^>/{if (seq) print seq; print; seq=""; next} {seq=seq $$0} END{if (seq) print seq}' | grep -B 1 -P "$$M1|$$M1_REVCOMP" > $(word 1,$^).tomap.fa
	$(MINIMAP) -a -x asm5 --cs -t1 -z 3000,1500 $(TERC_REGION_REFERENCE_FA) $(word 1,$^).tomap.fa  | $(SAMTOOLS) sort --threads 1 >  $(word 1,$^).tomap.bam
	$(SAMTOOLS) index $(word 1,$^).tomap.bam
	

	@CHROM_REGION="chr3"
	@BEGIN_REGION=169770169
	@END_REGION=169771552
	python $(MSA_VIEW_PY) CONSENSUS_REGION $(word 1,$^).tomap.bam  $$CHROM_REGION $$BEGIN_REGION $$END_REGION > $(word 1,$^).tomap.bam.allele.fa
	cat $(word 1,$^).tomap.bam.allele.fa  | awk '$$0 !~ ">" ' |  awk -v P=$(PN) -v H=$(HAPNAME) -F"$$M1" '{print P "\t" H "\t" (NF-1)  }' >  $(word 1,$^).tomap.bam.allele.motif_count.csv
	cat $(word 1,$^).tomap.bam.allele.fa  | awk '$$0 !~ ">" ' |  awk -v P=$(PN) -v H=$(HAPNAME)  '{print P "\t" H "\t" length($$1) }' >  $(word 1,$^).tomap.bam.allele.length.csv

	paste -d "\t"  $(word 1,$^).tomap.bam.allele.motif_count.csv  $(word 1,$^).tomap.bam.allele.length.csv | awk '{print $$1 "\t" $$2 "\t" $$3 "\t" $$6}' > $(word 1,$^).tomap.bam.allele.info.csv 
	echo "done." >  $@