#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


def buildSamplesChannel() {
    if (params.manifest) {
        return Channel
            .fromPath(params.manifest, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def researchId = row.research_id
                def hapnames = ['hap1', 'hap2']
                def fastas = [
                    file(row.assembly_hap1_fa, checkIfExists: true),
                    file(row.assembly_hap2_fa, checkIfExists: true),
                ]
                tuple(researchId as String, researchId as String, hapnames, fastas)
            }
    }

}

process PACBIO_ASM_ALLELE_INFO {
    tag "${sample_id}"

    cpus 1
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), val(pn), val(hapnames), path(fastas)
    path reference_fa

    output:
    tuple val(sample_id),
        path("*.tomap.fa"),
        path("*.tomap.bam"),
        path("*.tomap.bam.bai"),
        path("*.tomap.bam.allele.fa"),
        path("*.tomap.bam.allele.motif_count.csv"),
        path("*.tomap.bam.allele.length.csv"),
        path("*.tomap.bam.allele.info.csv"),
        path("*.3q26.2-TR.PacBio_asm_allele.info.done")

    script:
    def motif = params.motif
    def motifRevcomp = params.motif_revcomp
    def chrom = params.chrom_region
    def begin = params.begin_region
    def end = params.end_region
    def hapNameArray = hapnames.collect { "\"${it}\"" }.join(' ')
    def fastaNameArray = fastas.collect { "\"${it.name}\"" }.join(' ')
    """
    set -euo pipefail

    M1="${motif}"
    M1_REVCOMP="${motifRevcomp}"
    HAP_NAMES=( ${hapNameArray} )
    FASTA_FILES=( ${fastaNameArray} )

    write_empty_result() {
        local prefix="\$1"
        local hap_name="\$2"
        local bam="\$prefix.tomap.bam"
        local fai="${reference_fa}.fai"

        if [[ ! -f "\$fai" ]]; then
            samtools faidx "${reference_fa}"
        fi

        awk 'BEGIN {OFS="\\t"} {print "@SQ", "SN:" \$1, "LN:" \$2}' "\$fai" \\
          | samtools view -b -o "\$bam" -

        samtools index "\$bam"
        : > "\$prefix.tomap.bam.allele.fa"
        printf "%s\\t%s\\t0\\n" "${pn}" "\$hap_name" > "\$prefix.tomap.bam.allele.motif_count.csv"
        printf "%s\\t%s\\t0\\n" "${pn}" "\$hap_name" > "\$prefix.tomap.bam.allele.length.csv"
        printf "%s\\t%s\\t0\\t0\\n" "${pn}" "\$hap_name" > "\$prefix.tomap.bam.allele.info.csv"
        echo "done." > "\$prefix.3q26.2-TR.PacBio_asm_allele.info.done"
    }

    for idx in "\${!FASTA_FILES[@]}"; do
        HAP_NAME="\${HAP_NAMES[\$idx]}"
        FASTA="\${FASTA_FILES[\$idx]}"
        PREFIX="${sample_id}.\${HAP_NAME}"

        if [[ "\$FASTA" == *.gz ]]; then
            gzip -cd "\$FASTA" | awk -v m1="\$M1" -v m2="\$M1_REVCOMP" '
                /^>/ {
                    if (header && seq ~ (m1 "|" m2)) {
                        print header
                        print seq
                    }
                    header=\$0
                    seq=""
                    next
                }
                {
                    seq=seq \$0
                }
                END {
                    if (header && seq ~ (m1 "|" m2)) {
                        print header
                        print seq
                    }
                }
            ' > "\$PREFIX.tomap.fa"
        else
            awk -v m1="\$M1" -v m2="\$M1_REVCOMP" '
                /^>/ {
                    if (header && seq ~ (m1 "|" m2)) {
                        print header
                        print seq
                    }
                    header=\$0
                    seq=""
                    next
                }
                {
                    seq=seq \$0
                }
                END {
                    if (header && seq ~ (m1 "|" m2)) {
                        print header
                        print seq
                    }
                }
            ' "\$FASTA" > "\$PREFIX.tomap.fa"
        fi

        if [[ ! -s "\$PREFIX.tomap.fa" ]]; then
            write_empty_result "\$PREFIX" "\$HAP_NAME"
            continue
        fi

        minimap2 -a -x asm5 --cs -t ${task.cpus} -z 3000,1500 "${reference_fa}" "\$PREFIX.tomap.fa" \\
          > "\$PREFIX.tomap.sam" 2> "\$PREFIX.tomap.minimap2.log"

        if [[ ! -s "\$PREFIX.tomap.sam" ]] || ! samtools view -H "\$PREFIX.tomap.sam" >/dev/null 2>&1; then
            echo "minimap2 did not produce a valid SAM header for \$PREFIX" >&2
            if [[ -s "\$PREFIX.tomap.minimap2.log" ]]; then
                cat "\$PREFIX.tomap.minimap2.log" >&2
            fi
            write_empty_result "\$PREFIX" "\$HAP_NAME"
            rm -f "\$PREFIX.tomap.sam"
            continue
        fi

        samtools sort --threads ${task.cpus} -o "\$PREFIX.tomap.bam" "\$PREFIX.tomap.sam"
        rm -f "\$PREFIX.tomap.sam"

        samtools index "\$PREFIX.tomap.bam"

        if [[ "\$(samtools view -c "\$PREFIX.tomap.bam" "${chrom}:${begin}-${end}")" -eq 0 ]]; then
            write_empty_result "\$PREFIX" "\$HAP_NAME"
            continue
        fi

        python "${params.msa_view_py}" CONSENSUS_REGION "\$PREFIX.tomap.bam" "${chrom}" "${begin}" "${end}" \\
          > "\$PREFIX.tomap.bam.allele.fa"

        awk '\$0 !~ ">"' "\$PREFIX.tomap.bam.allele.fa" \\
          | awk -v P="${pn}" -v H="\$HAP_NAME" -F"\${M1}" '{print P "\\t" H "\\t" (NF-1)}' \\
          > "\$PREFIX.tomap.bam.allele.motif_count.csv"

        awk '\$0 !~ ">"' "\$PREFIX.tomap.bam.allele.fa" \\
          | awk -v P="${pn}" -v H="\$HAP_NAME" '{print P "\\t" H "\\t" length(\$1)}' \\
          > "\$PREFIX.tomap.bam.allele.length.csv"

        paste -d "\\t" "\$PREFIX.tomap.bam.allele.motif_count.csv" "\$PREFIX.tomap.bam.allele.length.csv" \\
          | awk '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$6}' > "\$PREFIX.tomap.bam.allele.info.csv"

        echo "done." > "\$PREFIX.3q26.2-TR.PacBio_asm_allele.info.done"
    done
    """

}

workflow {
    samples_ch = buildSamplesChannel()
    reference_ch = Channel.value(file(params.reference_fa, checkIfExists: true))

    PACBIO_ASM_ALLELE_INFO(samples_ch, reference_ch)
}
