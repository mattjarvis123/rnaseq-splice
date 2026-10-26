#!/usr/bin/env Rscript
# Render the rnaseq-splice results report.
#
# Usage:
#   Rscript render_report.R --results_dir results --samplesheet samplesheet.csv \
#       --contrasts contrasts.csv [--output rnaseq_splice_report.html] [--aligner star_salmon] \
#       [--params_file params.yaml] \
#       [--fdr_de 0.01] [--lfc_de 0] [--fdr_splicing 0.01] [--dpsi_splicing 0] [--fdr_gsea 0.05] \
#       [--rmats_counts JC] [--heatmap_top_n 50] [--event_heatmap_per_type 15] [--event_heatmap_scale z|psi] \
#       [--lookup_events_per_gene 3] [--kegg_maps true] [--kegg_top_n 8] [--kegg_cache <dir>] \
#       [--table_max_rows 5000] [--plot_max_points 20000] [--plot_ns_points 2000] \
#       [--title "rnaseq-splice report"] [--project ""] [--results_label <results folder name shown in the report; default: basename of --results_dir>]
#
# Significance cut-offs (the pipeline itself never filters; only the report applies them). Each one is
# taken from the first of these that sets it:
#   1. this script's option            --fdr_de  --lfc_de  --fdr_splicing  --dpsi_splicing  --fdr_gsea
#   2. --params_file (the pipeline's    report_de_fdr  report_de_lfc  report_splicing_fdr
#      params.yaml or a params JSON)    report_splicing_dpsi  gsea_fdr
#   3. the run's own record, <results_dir>/pipeline_info/params_*.json (latest), same names as 2
#   4. defaults: FDR 0.01, no log2FC cut-off, FDR 0.01, no dPSI cut-off, gene set FDR 0.05
#
# Full tables of significant genes/events are written as TSV to <output name>_tables/ next to the
# report and linked from it; keep that folder together with the HTML when sharing.
#
# Needs: rmarkdown, flexdashboard, knitr, dplyr, tidyr, ggplot2, plotly, DT, htmltools,
#        UpSetR, eulerr, gridExtra, yaml, jsonlite (and pandoc).
#
# KEGG maps: for significant gene sets from MSigDB's KEGG collection (c2.cp.kegg / kegg_legacy) in the
# GSEA results, KEGG pathway maps are coloured by log2FC (expression) or lead-event dPSI (splicing).
# KGML/PNG files are read from --kegg_cache (default: R_user_dir("rnaseq-splice", "cache")/kegg) and
# downloaded from rest.kegg.jp when missing. KEGG's API is free for academic use only.
#
# Gene set enrichment pages are added when <results_dir>/<aligner>/gsea/<contrast>/ holds the
# pipeline's --gsea output (*.fgsea_de.tsv, *.fgsea_splicing.tsv).

# The report contains non-ASCII characters (ΔPSI, −); containers often default to the C locale
invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "C.UTF-8")))

args <- commandArgs(trailingOnly = TRUE)
opts <- list(
    results_dir = "results", samplesheet = "samplesheet.csv", contrasts = "contrasts.csv",
    output = "rnaseq_splice_report.html", aligner = "star_salmon", params_file = "",
    fdr_de = NA, lfc_de = NA, fdr_splicing = NA, dpsi_splicing = NA, fdr_gsea = NA,
    rmats_counts = "JC", heatmap_top_n = "50", event_heatmap_per_type = "15", event_heatmap_scale = "z", lookup_events_per_gene = "3",
    kegg_maps = "true", kegg_top_n = "8", kegg_cache = "",
    table_max_rows = "5000", plot_max_points = "20000", plot_ns_points = "2000",
    title = "rnaseq-splice report", project = "", results_label = ""
)
i <- 1
while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (!key %in% names(opts)) stop("Unknown option: ", args[i])
    opts[[key]] <- args[i + 1]
    i <- i + 2
}

# Render from a scratch copy so knitr's intermediate files stay out of the results folder
script_dir <- local({
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    if (length(f)) dirname(normalizePath(f)) else getwd()
})
work <- tempfile("rnaseq_splice_report_")
dir.create(work)
for (f in list.files(script_dir, pattern = "\\.(Rmd|js|tsv)$|^kegg_overlay\\.R$", full.names = TRUE)) {
    file.copy(f, work)
}

# Resolve paths now: render() evaluates its arguments after changing directory
# normalizePath() leaves a not-yet-existing file's path relative, so make the folder and normalise that
dir.create(dirname(opts$output), showWarnings = FALSE, recursive = TRUE)
output      <- file.path(normalizePath(dirname(opts$output)), basename(opts$output))
results_dir <- normalizePath(opts$results_dir, mustWork = TRUE)
samplesheet <- normalizePath(opts$samplesheet, mustWork = TRUE)
contrasts   <- normalizePath(opts$contrasts, mustWork = TRUE)
dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
tables_link <- paste0(sub("\\.html$", "", basename(output)), "_tables")

# ---- significance cut-offs: script option > --params_file > run's params JSON > default -------------
read_params <- function(f) {
    if (!nzchar(f)) return(list())
    if (!file.exists(f)) stop("--params_file not found: ", f)
    if (grepl("\\.json$", f, ignore.case = TRUE)) jsonlite::fromJSON(f) else yaml::read_yaml(f)
}
file_params <- read_params(opts$params_file)
run_json    <- sort(list.files(file.path(results_dir, "pipeline_info"), pattern = "^params_.*\\.json$", full.names = TRUE))
run_params  <- if (length(run_json)) jsonlite::fromJSON(run_json[length(run_json)]) else list()
cutoffs <- list(  # report option = c(pipeline parameter, default)
    fdr_de        = c("report_de_fdr",        0.01),
    lfc_de        = c("report_de_lfc",        0),
    fdr_splicing  = c("report_splicing_fdr",  0.01),
    dpsi_splicing = c("report_splicing_dpsi", 0),
    fdr_gsea      = c("gsea_fdr",             0.05)
)
for (k in names(cutoffs)) {
    pipeline_name <- cutoffs[[k]][1]
    if (!is.na(opts[[k]])) {
        value <- opts[[k]]; source <- "command line"
    } else if (!is.null(file_params[[pipeline_name]])) {
        value <- file_params[[pipeline_name]]; source <- basename(opts$params_file)
    } else if (!is.null(run_params[[pipeline_name]])) {
        value <- run_params[[pipeline_name]]; source <- file.path("pipeline_info", basename(run_json[length(run_json)]))
    } else {
        value <- cutoffs[[k]][2]; source <- "default"
    }
    opts[[k]] <- as.numeric(value)
    message(sprintf("  %-14s = %-6s (%s: %s)", k, opts[[k]], source, pipeline_name))
}
tables_dir  <- file.path(dirname(output), tables_link)
rmarkdown::render(
    file.path(work, "rnaseq_splice_report.Rmd"),
    output_file   = basename(output),
    output_dir    = dirname(output),
    intermediates_dir = work,
    knit_root_dir = getwd(),
    params = list(
        results_dir   = results_dir,
        results_label = opts$results_label,
        aligner       = opts$aligner,
        samplesheet   = samplesheet,
        contrasts     = contrasts,
        fdr_de        = as.numeric(opts$fdr_de),
        fdr_splicing  = as.numeric(opts$fdr_splicing),
        lfc_de        = as.numeric(opts$lfc_de),
        dpsi_splicing = as.numeric(opts$dpsi_splicing),
        table_max_rows  = as.integer(opts$table_max_rows),
        plot_max_points = as.integer(opts$plot_max_points),
        plot_ns_points  = as.integer(opts$plot_ns_points),
        tables_dir    = tables_dir,
        tables_link   = tables_link,
        fdr_gsea      = as.numeric(opts$fdr_gsea),
        rmats_counts  = opts$rmats_counts,
        heatmap_top_n = as.integer(opts$heatmap_top_n),
        event_heatmap_per_type = as.integer(opts$event_heatmap_per_type),
        event_heatmap_scale = opts$event_heatmap_scale,
        lookup_events_per_gene = as.integer(opts$lookup_events_per_gene),
        kegg_maps     = tolower(opts$kegg_maps) %in% c("true", "t", "1", "yes"),
        kegg_top_n    = as.integer(opts$kegg_top_n),
        kegg_cache    = if (nzchar(opts$kegg_cache)) normalizePath(opts$kegg_cache, mustWork = FALSE) else "",
        title         = opts$title,
        project       = opts$project
    ),
    envir = new.env(),
    quiet = TRUE
)
message("Report written to ", output)
