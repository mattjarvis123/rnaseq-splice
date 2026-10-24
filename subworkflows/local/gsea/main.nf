//
// Gene set enrichment of the limma/voom and rMATS results (fgsea). rnaseq-splice addition.
//
//   FGSEA_DE       : preranked GSEA of every tested gene, ranked by the limma statistic
//   FGSEA_SPLICING : over-representation of genes with significant rMATS events
//                    (higher- and lower-inclusion lists), background = genes with a tested event
//
// Either input channel may be empty (e.g. --limma without --rmats); that half is then skipped.
//
include { FGSEA_DE       } from '../../../modules/local/fgsea/de/main'
include { FGSEA_SPLICING } from '../../../modules/local/fgsea/splicing/main'

workflow GSEA {

    take:
    ch_limma_results // channel: [ val(meta), path(<contrast>.limma.results.tsv) ]
    ch_gene_map      // channel: path(gene counts with gene_id, gene_name), value channel
    ch_mats          // channel: [ val(meta), [ path(*.MATS.{JC,JCEC}.txt) ] ]
    ch_gmt           // channel: [ path(gmt), ... ], value channel
    rmats_counts     // string: 'JC' (junction counts) or 'JCEC'

    main:
    FGSEA_DE (
        ch_limma_results,
        ch_gene_map,
        ch_gmt
    )

    // Keep only the requested count type, e.g. SE.MATS.JC.txt, MXE.MATS.JC.txt, ...
    ch_mats_counts = ch_mats.map { meta, files ->
        [ meta, files.findAll { f -> f.name.endsWith(".MATS.${rmats_counts}.txt") } ]
    }
    FGSEA_SPLICING (
        ch_mats_counts,
        ch_gmt
    )

    emit:
    de_results       = FGSEA_DE.out.results        // channel: [ val(meta), path(*.fgsea_de.tsv) ]
    de_ranks         = FGSEA_DE.out.ranks          // channel: [ val(meta), path(*.fgsea_ranks.tsv) ]
    splicing_results = FGSEA_SPLICING.out.results // channel: [ val(meta), path(*.fgsea_splicing.tsv) ]
    splicing_genes   = FGSEA_SPLICING.out.genes   // channel: [ val(meta), path(*.fgsea_splicing_genes.tsv) ]
}
