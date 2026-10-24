//
// Gene set over-representation analysis (fgsea::fora) of the genes with significant rMATS
// events of one contrast, separately for higher- and lower-inclusion events.
// Script: bin/fgsea_splicing.R (also runs on its own, outside the pipeline). rnaseq-splice addition.
//
process FGSEA_SPLICING {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-fgsea:1.32.2--r44he5774e6_0' :
        'biocontainers/bioconductor-fgsea:1.32.2--r44he5774e6_0' }"

    input:
    tuple val(meta), path(mats, stageAs: 'mats/*')   // the five <TYPE>.MATS.JC.txt files of one contrast
    path gmt                                          // one or more GMT files of gene-symbol gene sets

    output:
    tuple val(meta), path("*.fgsea_splicing.tsv")      , emit: results
    tuple val(meta), path("*.fgsea_splicing_genes.tsv"), emit: genes
    tuple val("${task.process}"), val('fgsea'), eval('Rscript -e "cat(as.character(packageVersion(\'fgsea\')))"'), topic: versions, emit: versions_fgsea

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def files  = (mats instanceof List ? mats : [mats]).join(',')
    def gmts   = (gmt instanceof List ? gmt : [gmt]).join(',')
    """
    fgsea_splicing.R \\
        --mats ${files} \\
        --gmt ${gmts} \\
        --contrast ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fgsea_splicing.tsv ${prefix}.fgsea_splicing_genes.tsv
    """
}
