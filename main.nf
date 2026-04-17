#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


def buildSamplesChannel() {
    if (params.manifest) {
        return Channel
            .fromPath(params.manifest, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                tuple(row.researchId as String, row.assembly_hap1_fa, row.assembly_hap2_fa)
            }
    }

}

process PACBIO_ASM_ALLELE_INFO {
    tag "${sample_id}"

    cpus 1
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(hap1), path(hap2)
    path(reference_fa)

    output:
    tuple val(sample_id),
        path("*.tomap.*"),
        path("*.motif_counts.csv"),
        path("*.lengths.csv")
    
    script:
    """
    PN="${sample_id}"
    HAP_1_FA_GZ="${hap1}" 
    HAP_2_FA_GZ="${hap2}" 
    REFERENCE="${reference_fa}" 
    OUT_M1_COUNT_FN="${sample_id}.motif_counts.csv" 
    OUT_LENGTH_FN="${sample_id}.lengths.csv" 


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


    # Get both alleles.
    # Allele 1
    zcat \${HAP_1_FA_GZ} | awk '/^>/{if (seq) print seq; print; seq=""; next} {seq=seq \$0} END {if (seq) print seq}' | grep -B 1 -P "\${M1}|\${M1_REVCOMP}" > \${HAP_1_FA_GZ}.tomap.fa
    zcat \${HAP_2_FA_GZ} | awk '/^>/{if (seq) print seq; print; seq=""; next} {seq=seq \$0} END {if (seq) print seq}' | grep -B 1 -P "\${M1}|\${M1_REVCOMP}" > \${HAP_2_FA_GZ}.tomap.fa

    if [[ -s "\${HAP_1_FA_GZ}.tomap.fa" && -s "\${HAP_2_FA_GZ}.tomap.fa" ]]; then

	    \${MINIMAP} -a -x asm5 --cs -t1 -z 3000,1500 \${REFERENCE} \${HAP_1_FA_GZ}.tomap.fa | \${SAMTOOLS} sort --threads 1 >  \${HAP_1_FA_GZ}.tomap.bam
	    \${SAMTOOLS} index \${HAP_1_FA_GZ}.tomap.bam
	    python \${MSA_VIEW_PY} CONSENSUS_REGION \${HAP_1_FA_GZ}.tomap.bam  \${CHROM_REGION} \${BEGIN_REGION} \${END_REGION} > \${HAP_1_FA_GZ}.tomap.bam.allele.fa

	    # Allele 2
	    \${MINIMAP} -a -x asm5 --cs -t1 -z 3000,1500 \${REFERENCE} \${HAP_2_FA_GZ}.tomap.fa | \${SAMTOOLS} sort --threads 1 >  \${HAP_2_FA_GZ}.tomap.bam
	    \${SAMTOOLS} index \${HAP_2_FA_GZ}.tomap.bam
	    python \${MSA_VIEW_PY} CONSENSUS_REGION \${HAP_2_FA_GZ}.tomap.bam  \${CHROM_REGION} \${BEGIN_REGION} \${END_REGION} > \${HAP_2_FA_GZ}.tomap.bam.allele.fa


	    # Now get info on the alleles.
	    HAP1_M1_COUNT=\$(cat \${HAP_1_FA_GZ}.tomap.bam.allele.fa | awk '$0 !~ ">" ' |  awk -F"\${M1}" '{print  (NF-1)  }' | bc ) 
	    HAP1_LENGTH=\$(cat \${HAP_1_FA_GZ}.tomap.bam.allele.fa | awk '$0 !~ ">" ' |  awk  '{print  length($1) }' | bc )

	    HAP2_M1_COUNT=\$(cat \${HAP_2_FA_GZ}.tomap.bam.allele.fa | awk '$0 !~ ">" ' |  awk -F"\${M1}" '{print  (NF-1)  }' | bc ) 
	    HAP2_LENGTH=\$(cat \${HAP_2_FA_GZ}.tomap.bam.allele.fa | awk '$0 !~ ">" ' |  awk  '{print  length($1) }' | bc )

	    # Final outputs.
	    echo -e "\${PN}\t\${HAP1_M1_COUNT}\t\${HAP2_M1_COUNT}" > \${OUT_M1_COUNT_FN}
	    echo -e "\${PN}\t\${HAP1_LENGTH}\t\${HAP2_LENGTH}" > \${OUT_LENGTH_FN}
    else
	    echo -e "\${PN}\tNA\tNA" > \${OUT_M1_COUNT_FN}
	    echo -e "\${PN}\tNA\tNA" > \${OUT_LENGTH_FN}
    fi;
    """
}

workflow {
    samples_ch = buildSamplesChannel()
    reference_ch = Channel.value(file(params.reference_fa, checkIfExists: true))

    PACBIO_ASM_ALLELE_INFO(samples_ch, reference_ch)
}
