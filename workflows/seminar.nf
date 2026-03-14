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

workflow SEMINAR {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_gtf = Channel.value(file(params.gtf))
    ch_star_index = Channel.value(file(params.star_index))
    ch_salmon_index = Channel.value(file(params.salmon_index))
    ch_transcriptome = Channel.value(file(params.transcriptome))
    ch_multiqc_files = Channel.empty()


     
    FASTQC(ch_samplesheet)
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip)
    
    TRIMGALORE(ch_samplesheet)
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log)    
    
    STAR_ALIGN(TRIMGALORE.out.reads,ch_star_index,ch_gtf,false)
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log)    

    QUALIMAP_RNASEQ(STAR_ALIGN.out.bam, ch_gtf)
    ch_versions = ch_versions.mix(QUALIMAP_RNASEQ.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_RNASEQ.out.report)

    DUPRADAR(STAR_ALIGN.out.bam, ch_gtf)
    ch_versions = ch_versions.mix(DUPRADAR.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(DRUPRADAR.out.plot)

    SALMON_QUANT(TRIMGALORE.out.reads, ch_salmon_index, ch_gtf, ch_transcriptome, false, false)
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.log)

    MULTIQC(ch_multiqc_files)
    ch_versions = ch_versions.mix(MULTIQC.out.versions) 

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'seminar_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
