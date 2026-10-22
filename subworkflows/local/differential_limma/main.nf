//
// Differential gene expression with limma/voom
//
// Uses the nf-core modules CUSTOM_MATRIXFILTER (low-count filtering) and
// LIMMA_DIFFERENTIAL (TMM normalisation, voom, lmFit, eBayes, topTable). Runs once per
// row of the contrasts sheet that is shared with rMATS.
//

include { CUSTOM_MATRIXFILTER } from '../../../modules/nf-core/custom/matrixfilter/main'
include { LIMMA_DIFFERENTIAL  } from '../../../modules/nf-core/limma/differential/main'

workflow DIFFERENTIAL_LIMMA {
    take:
    ch_counts        // channel: [ val(meta), path(gene_counts.tsv) ], merged gene x sample matrix from rnaseq
    ch_samples       // channel: [ val(sample_id), val(condition) ], in samplesheet order
    ch_contrastsheet // channel: [ contrast:, treatment:, control: ]

    main:

    //
    // Sample sheet for the R templates: one row per sample (technical replicates have
    // already been merged by rnaseq, so duplicate sample ids are dropped)
    //
    ch_de_samplesheet = ch_samples
        .unique { sample_id, _condition -> sample_id }
        .map { sample_id, condition -> "${sample_id},${condition}" }
        .collectFile(name: 'limma_samplesheet.csv', seed: 'sample,condition', newLine: true, sort: false)
        .map { csv -> [ [ id: 'samplesheet' ], csv ] }
        .first()

    //
    // MODULE: CUSTOM_MATRIXFILTER
    //
    // Drops the gene_name annotation column, orders the samples as in the sample sheet and
    // removes lowly expressed genes (see ext.args in conf/modules/limma.config)
    //
    CUSTOM_MATRIXFILTER (
        ch_counts.first(),
        ch_de_samplesheet
    )

    //
    // MODULE: LIMMA_DIFFERENTIAL
    //
    // contrast_variable is always the samplesheet 'condition' column; the model formula
    // and contrast string are left unset so the template builds ~ condition itself
    //
    ch_limma_contrasts = ch_contrastsheet
        .map { row -> [ [ id: row.contrast, treatment: row.treatment, control: row.control ], 'condition', row.control, row.treatment, '', '' ] }

    LIMMA_DIFFERENTIAL (
        ch_limma_contrasts,
        ch_de_samplesheet
            .combine(CUSTOM_MATRIXFILTER.out.filtered.map { _meta, tsv -> tsv })
            .first()
    )

    emit:
    filtered_counts   = CUSTOM_MATRIXFILTER.out.filtered          // channel: [ val(meta), path(tsv) ]
    results           = LIMMA_DIFFERENTIAL.out.results            // channel: [ val(meta), path(*.limma.results.tsv) ]
    md_plot           = LIMMA_DIFFERENTIAL.out.md_plot            // channel: [ val(meta), path(png) ]
    normalised_counts = LIMMA_DIFFERENTIAL.out.normalised_counts  // channel: [ val(meta), path(tsv) ]
    model             = LIMMA_DIFFERENTIAL.out.model              // channel: [ val(meta), path(txt) ]
}
