#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


def buildSamplesChannel() {
    if (params.manifest) {
        return Channel
            .fromPath(params.manifest, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                tuple(row.research_id as String, row.assembly_hap1_fa, row.assembly_hap2_fa)
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
    path(reference_fa_fai)

    output:
    tuple val(sample_id),
        path("*.tomap.*"),
        path("*.motif_counts.csv"),
        path("*.lengths.csv")
    
    script:
    """
    process.sh ${sample_id} ${hap1} ${hap2} ${reference_fa} 
    """
}

workflow {
    samples_ch = buildSamplesChannel()
    
    ref= params.reference_fa
    ref_fai = "${params.reference_fa}.fai"

    PACBIO_ASM_ALLELE_INFO(samples_ch, ref, ref_fai, "unzip")
}
