//
// Preranked gene set enrichment analysis (fgsea) of one limma/voom contrast.
// Script: bin/fgsea_de.R (also runs on its own, outside the pipeline). rnaseq-splice addition.
//
process FGSEA_DE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-fgsea:1.32.2--r44he5774e6_0' :
        'biocontainers/bioconductor-fgsea:1.32.2--r44he5774e6_0' }"

    input:
    tuple val(meta), path(results)   // <contrast>.limma.results.tsv
    path gene_map                    // TSV with gene_id and gene_name (merged gene counts)
    path gmt                         // one or more GMT files of gene-symbol gene sets

    output:
    tuple val(meta), path("*.fgsea_de.tsv")   , emit: results
    tuple val(meta), path("*.fgsea_ranks.tsv"), emit: ranks
    tuple val("${task.process}"), val('fgsea'), eval('Rscript -e "cat(as.character(packageVersion(\'fgsea\')))"'), topic: versions, emit: versions_fgsea

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gmts   = (gmt instanceof List ? gmt : [gmt]).join(',')
    """
    fgsea_de.R \\
        --results ${results} \\
        --gene_map ${gene_map} \\
        --gmt ${gmts} \\
        --contrast ${prefix} \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fgsea_de.tsv ${prefix}.fgsea_ranks.tsv
    """
}
