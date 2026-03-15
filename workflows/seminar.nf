/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { FASTQC }           from '../modules/nf-core/fastqc/main'
include { TRIMGALORE }       from '../modules/nf-core/trimgalore/main'
include { STAR_ALIGN }       from '../modules/nf-core/star/align/main'
include { SALMON_QUANT }     from '../modules/nf-core/salmon/quant/main'
include { QUALIMAP_RNASEQ }  from '../modules/nf-core/qualimap/rnaseq/main'
include { DUPRADAR }         from '../modules/nf-core/dupradar/main'
include { MULTIQC }          from '../modules/nf-core/multiqc/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SEMINAR {

    take:
    ch_samplesheet

    main:
    ch_gtf_tuple      = Channel.value([[id: 'gtf'], file(params.gtf)])
    ch_gtf_path       = Channel.value(file(params.gtf))
    ch_star_index     = Channel.value([[id: 'star_index'], file(params.star_index)])
    ch_salmon_index   = Channel.value(file(params.salmon_index))
    ch_transcriptome  = Channel.value(file(params.transcriptome))
    ch_alignment_mode = Channel.value(false)
    ch_lib_type       = Channel.value(false)

    ch_multiqc_files = Channel.empty()

    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { it[1] })

    TRIMGALORE(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.zip.map { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.map { it[1] })

    STAR_ALIGN(TRIMGALORE.out.reads, ch_star_index, ch_gtf_tuple, false)
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.map { it[1] })

    QUALIMAP_RNASEQ(STAR_ALIGN.out.bam, ch_gtf_tuple)
    ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_RNASEQ.out.results.map { it[1] })

    DUPRADAR(STAR_ALIGN.out.bam, ch_gtf_tuple)
    ch_multiqc_files = ch_multiqc_files.mix(DUPRADAR.out.multiqc.map { it[1] })

    SALMON_QUANT(TRIMGALORE.out.reads, ch_salmon_index, ch_gtf_path, ch_transcriptome, ch_alignment_mode, ch_lib_type)
    ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.results.map { it[1] })

    ch_multiqc_files = ch_multiqc_files.collect()
    ch_multiqc_config       = Channel.value([])
    ch_extra_multiqc_config = Channel.value([])
    ch_multiqc_logo         = Channel.value([])
    ch_replace_names        = Channel.value([])
    ch_sample_names         = Channel.value([])

    MULTIQC(
        ch_multiqc_files,
        ch_multiqc_config,
        ch_extra_multiqc_config,
        ch_multiqc_logo,
        ch_replace_names,
        ch_sample_names
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
