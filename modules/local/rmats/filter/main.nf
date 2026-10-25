//
// Read-coverage filter for the rMATS results of one contrast. rnaseq-splice addition.
//
// Keeps events with at least `--min_reads` junction reads (IJC + SJC) in every replicate of
// both conditions and recomputes the Benjamini-Hochberg FDR over the retained events
// (rMATS's own FDR is kept in the column FDR_unfiltered). The unfiltered rMATS tables are
// published unchanged; the filtered ones, with the same file names, go to <contrast>/filtered/.
// Script: bin/rmats_filter.R (base R, so it runs in the rMATS container).
//
process RMATS_FILTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/4697fb2978d58ad7209a2d7b267bc02cf56cb8c256d787767b482c97346330d6/data'
        : 'community.wave.seqera.io/library/r-pairadise_rmats:4a32921359285283'}"

    input:
    tuple val(meta), path(mats, stageAs: 'mats/*')   // the <TYPE>.MATS.{JC,JCEC}.txt files of one contrast

    output:
    tuple val(meta), path("${prefix}/filtered/*.MATS.{JC,JCEC}.txt"), emit: mats
    tuple val(meta), path("${prefix}/filtered/filter_summary.tsv")   , emit: summary
    // No `eval` version output: Nextflow numbers eval outputs across the whole pipeline, so a new
    // one renumbers those of other processes and invalidates their -resume cache. The script is
    // base R in the rMATS container, whose versions RMATS_POST already reports.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def files = (mats instanceof List ? mats : [mats]).join(',')
    """
    rmats_filter.R \\
        --mats ${files} \\
        --outdir ${prefix}/filtered \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/filtered
    for event in SE MXE A3SS A5SS RI; do
        touch ${prefix}/filtered/\${event}.MATS.JC.txt
        touch ${prefix}/filtered/\${event}.MATS.JCEC.txt
    done
    touch ${prefix}/filtered/filter_summary.tsv
    """
}
