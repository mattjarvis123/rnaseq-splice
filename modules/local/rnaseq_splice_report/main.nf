//
// rnaseq-splice results report (R Markdown / flexdashboard HTML). rnaseq-splice addition.
//
// Final step of the pipeline: lays the limma/voom, rMATS, fgsea and MultiQC outputs out in the
// same folder structure as the published results and renders assets/report/render_report.R on
// them. Writes <prefix>.html and <prefix>_tables/ (full TSV tables and the KEGG map images the
// report links to); keep the two together.
//
// No `eval` version output on purpose: Nextflow numbers eval outputs across the whole pipeline,
// so adding one would renumber those of other processes and invalidate their -resume cache.
//
process RNASEQ_SPLICE_REPORT {
    tag "${prefix}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_pandoc_r-rmarkdown_r-flexdashboard_pruned:a00b89db443fbdef' :
        'community.wave.seqera.io/library/r-base_pandoc_r-rmarkdown_r-flexdashboard_pruned:23c75384d7452c07' }"

    input:
    path report_src , stageAs: 'report_src'         // assets/report: render_report.R, *.Rmd, *.js, *.tsv
    path samplesheet, stageAs: 'samplesheet.csv'
    path contrasts  , stageAs: 'contrasts.csv'
    path gene_counts, stageAs: 'inputs/gene_counts.tsv'   // gene_id, gene_name, ... (names only)
    path limma      , stageAs: 'inputs/limma/*'     // <contrast>.limma.results.tsv, <contrast>.normalised_counts.tsv
    val  rmats_ids                                  // contrast id of each folder in rmats_dirs, same order
    path rmats_dirs , stageAs: 'inputs/rmats/*'     // one rMATS output folder per contrast
    path gsea       , stageAs: 'inputs/gsea/*'      // <contrast>.fgsea_*.tsv
    path multiqc    , stageAs: 'inputs/multiqc/*'   // MultiQC *_data folder
    path versions   , stageAs: 'inputs/software_versions.yml'
    val  run_params                                 // map of pipeline parameters quoted in the report
    val  aligner

    output:
    path "*.html"   , emit: html
    path "*_tables" , emit: tables, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: 'rnaseq_splice_report'
    def ids    = rmats_ids instanceof List ? rmats_ids : [rmats_ids]
    def dirs   = rmats_dirs instanceof List ? rmats_dirs : (rmats_dirs ? [rmats_dirs] : [])
    def link_rmats = [ids, dirs].transpose().collect { id, d ->
        "ln -s \"\$PWD/${d}\" \"results/${aligner}/rmats/${id}\""   // d is the staged path, inputs/rmats/<folder>
    }.join('\n    ')
    // the parameters the report quotes, as the params JSON the pipeline writes to pipeline_info
    def json = '{\n' + run_params.collect { k, v ->
        def val = v == null ? 'null' : (v instanceof Number || v instanceof Boolean) ? v.toString() : '"' + v.toString().replace('\\', '\\\\').replace('"', '\\"') + '"'
        "  \"${k}\": ${val}"
    }.join(',\n') + '\n}'
    """
    # Results folder in the layout of the published results, built from links to the inputs
    mkdir -p results/${aligner}/limma results/${aligner}/rmats results/${aligner}/gsea results/multiqc results/pipeline_info
    [ -e inputs/gene_counts.tsv ] && ln -s "\$PWD/inputs/gene_counts.tsv" results/${aligner}/salmon.merged.gene_counts.tsv
    for f in inputs/limma/*.limma.results.tsv inputs/limma/*.normalised_counts.tsv; do
        [ -e "\$f" ] || continue
        b=\$(basename "\$f"); ct=\${b%.limma.results.tsv}; ct=\${ct%.normalised_counts.tsv}
        mkdir -p results/${aligner}/limma/\$ct && ln -s "\$PWD/\$f" results/${aligner}/limma/\$ct/\$b
    done
    ${link_rmats}
    for f in inputs/gsea/*.fgsea_*.tsv; do
        [ -e "\$f" ] || continue
        b=\$(basename "\$f"); ct=\${b%%.fgsea_*}
        mkdir -p results/${aligner}/gsea/\$ct && ln -s "\$PWD/\$f" results/${aligner}/gsea/\$ct/\$b
    done
    for d in inputs/multiqc/*; do [ -e "\$d" ] && ln -s "\$PWD/\$d" results/multiqc/; done
    [ -e inputs/software_versions.yml ] && ln -s "\$PWD/inputs/software_versions.yml" results/pipeline_info/nf_core_rnaseq_software_mqc_versions.yml
    cat <<'END_OF_PARAMS' > results/pipeline_info/params_report.json
${json}
END_OF_PARAMS

    # The report's R code reads the same environment variables as an interactive render
    export HOME=\$PWD
    Rscript report_src/render_report.R \\
        --results_dir results \\
        --samplesheet samplesheet.csv \\
        --contrasts contrasts.csv \\
        --aligner ${aligner} \\
        --output ${prefix}.html \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: 'rnaseq_splice_report'
    """
    touch ${prefix}.html
    mkdir ${prefix}_tables
    """
}
