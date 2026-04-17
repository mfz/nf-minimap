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

    output:
    tuple val(sample_id),
        path("*.tomap.*"),
        path("*.motif_counts.csv"),
        path("*.lengths.csv")
    
    script:
    """
    process.sh ${hap1} ${hap2} ${reference_fa} ${sample_id}
    """
}

workflow {
    samples_ch = buildSamplesChannel()
    reference_ch = Channel.value(file(params.reference_fa, checkIfExists: true))

    PACBIO_ASM_ALLELE_INFO(samples_ch, reference_ch)
}
