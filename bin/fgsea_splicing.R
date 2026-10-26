#!/usr/bin/env Rscript
#
# fgsea_splicing.R - gene set over-representation analysis of rMATS results (rnaseq-splice)
#
# Builds two gene lists from the rMATS events of one contrast and tests each against one or more
# GMT files with fgsea::fora (hypergeometric test, BH-adjusted):
#   Higher inclusion : genes with >= 1 significant event with dPSI > 0 (more inclusion in treatment)
#   Lower inclusion  : genes with >= 1 significant event with dPSI < 0
# A gene with significant events in both directions is counted in both lists. The background
# (universe) is every gene with at least one tested rMATS event, so the test asks which gene sets
# are over-represented among spliced genes relative to genes where splicing could be measured.
#
# Usage:
#   fgsea_splicing.R --mats SE.MATS.JC.txt,MXE.MATS.JC.txt,A3SS.MATS.JC.txt,A5SS.MATS.JC.txt,RI.MATS.JC.txt \
#       --gmt h.all.symbols.gmt,c2.cp.reactome.symbols.gmt --contrast <contrast> \
#       [--event_fdr 0.01] [--event_dpsi 0] [--min_size 10] [--max_size 500] [--fdr 0.05] [--outdir .]
#
#   --mats        comma-separated rMATS *.MATS.JC.txt (or JCEC) files of one contrast
#   --event_fdr   rMATS FDR below which an event is significant
#   --event_dpsi  |dPSI| above which an event is significant (0 = no effect-size cut-off)
#
# Output (in --outdir):
#   <contrast>.fgsea_splicing.tsv        one row per gene set and gene list, sorted by p-value
#   <contrast>.fgsea_splicing_genes.tsv  the gene lists and universe that were tested

suppressPackageStartupMessages({
    library(fgsea)
    library(data.table)
})

# ---- arguments -----------------------------------------------------------------------------
opts <- list(mats = NULL, gmt = NULL, contrast = NULL, event_fdr = "0.01", event_dpsi = "0",
             min_size = "10", max_size = "500", fdr = "0.05", outdir = ".")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2 != 0) stop("Arguments must be given as --name value pairs")
for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[i])
    if (!key %in% names(opts)) stop("Unknown option: ", args[i])
    opts[[key]] <- args[i + 1]
}
for (req in c("mats", "gmt", "contrast")) {
    if (is.null(opts[[req]])) stop("Missing required option --", req)
}
event_fdr  <- as.numeric(opts$event_fdr)
event_dpsi <- as.numeric(opts$event_dpsi)
fdr        <- as.numeric(opts$fdr)
dir.create(opts$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- gene sets -----------------------------------------------------------------------------
read_gmts <- function(files) {
    files <- trimws(strsplit(files, ",")[[1]])
    sets <- list(); source <- character()
    for (f in files) {
        p <- gmtPathways(f)
        new <- setdiff(names(p), names(sets))          # first file wins on duplicate names
        sets <- c(sets, p[new])
        source <- c(source, setNames(rep(sub("\\.gmt$", "", basename(f)), length(new)), new))
    }
    list(sets = lapply(sets, function(g) unique(g[nzchar(g)])), source = source)
}
gs <- read_gmts(opts$gmt)
message(sprintf("Read %d gene sets from %s", length(gs$sets), opts$gmt))

# ---- rMATS events -> gene lists ------------------------------------------------------------
mats_files <- trimws(strsplit(opts$mats, ",")[[1]])
events <- rbindlist(lapply(mats_files, function(f) {
    # sep/dec fixed: with 2 replicates per group the count columns ("18,5") look like decimal-comma
    # numbers, and fread's auto-detection would then switch to dec = "," and read FDR/dPSI as text
    d <- fread(f, sep = "\t", dec = ".", select = c("geneSymbol", "FDR", "IncLevelDifference"), data.table = TRUE)
    d[, event_type := sub("\\.MATS\\..*$", "", basename(f))]
    d
}))
events[, gene := gsub('"', "", geneSymbol)]
events <- events[!is.na(gene) & nzchar(gene) & gene != "NA" & !is.na(FDR)]
events[, sig := FDR < event_fdr & abs(IncLevelDifference) > event_dpsi]

universe <- unique(events$gene)
gene_lists <- list(
    `Higher inclusion` = unique(events[sig == TRUE & IncLevelDifference > 0, gene]),
    `Lower inclusion`  = unique(events[sig == TRUE & IncLevelDifference < 0, gene])
)
message(sprintf("Universe: %d genes with a tested event; %d with higher-inclusion and %d with lower-inclusion events",
                length(universe), length(gene_lists[[1]]), length(gene_lists[[2]])))

fwrite(rbind(data.frame(list = "universe", gene = universe),
             do.call(rbind, lapply(names(gene_lists), function(n) data.frame(list = n, gene = gene_lists[[n]])))),
       file.path(opts$outdir, paste0(opts$contrast, ".fgsea_splicing_genes.tsv")), sep = "\t")

# ---- fora ----------------------------------------------------------------------------------
out <- rbindlist(lapply(names(gene_lists), function(n) {
    genes <- gene_lists[[n]]
    if (!length(genes)) return(NULL)
    fo <- fora(pathways = gs$sets, genes = genes, universe = universe,
               minSize = as.integer(opts$min_size), maxSize = as.integer(opts$max_size))
    if (!nrow(fo)) return(NULL)
    data.table(
        contrast        = opts$contrast,
        gene_list       = n,
        pathway         = fo$pathway,
        collection      = unname(gs$source[fo$pathway]),
        size            = fo$size,                  # gene set size within the universe
        overlap         = fo$overlap,
        list_size       = length(genes),
        universe_size   = length(universe),
        fold_enrichment = (fo$overlap / length(genes)) / (fo$size / length(universe)),
        pval            = fo$pval,
        padj            = fo$padj,
        significant     = !is.na(fo$padj) & fo$padj < fdr,
        overlap_genes   = vapply(fo$overlapGenes, paste, character(1), collapse = ";")
    )
}))
if (!nrow(out)) {
    out <- data.table(contrast = character(), gene_list = character(), pathway = character(),
                      collection = character(), size = integer(), overlap = integer(), list_size = integer(),
                      universe_size = integer(), fold_enrichment = numeric(), pval = numeric(),
                      padj = numeric(), significant = logical(), overlap_genes = character())
}
setorder(out, pval)
fwrite(out, file.path(opts$outdir, paste0(opts$contrast, ".fgsea_splicing.tsv")), sep = "\t")

message(sprintf("%s: %d gene sets enriched among higher-inclusion genes, %d among lower-inclusion genes (FDR < %s)",
                opts$contrast, sum(out$significant & out$gene_list == "Higher inclusion"),
                sum(out$significant & out$gene_list == "Lower inclusion"), fdr))
