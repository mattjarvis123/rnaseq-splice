#!/usr/bin/env Rscript
#
# fgsea_de.R - preranked gene set enrichment analysis of a limma/voom contrast (rnaseq-splice)
#
# Ranks every tested gene by its limma statistic (moderated t by default) and runs fgsea
# (multilevel) against one or more GMT files of gene-symbol gene sets. A gene set is called
# Up (NES > 0) or Down (NES < 0) when its BH-adjusted p-value is below --fdr.
#
# Usage:
#   fgsea_de.R --results <contrast>.limma.results.tsv --gene_map salmon.merged.gene_counts.tsv \
#       --gmt h.all.symbols.gmt,c2.cp.reactome.symbols.gmt --contrast <contrast> \
#       [--rank_by t] [--min_size 10] [--max_size 500] [--fdr 0.05] [--seed 42] [--threads 1] \
#       [--outdir .]
#
#   --results   limma results table (gene_id, logFC, AveExpr, t, P.Value, adj.P.Val, B)
#   --gene_map  any TSV with gene_id and gene_name columns (e.g. the merged Salmon gene counts);
#               limma reports Ensembl ids, the GMT files use gene symbols
#   --gmt       comma-separated GMT files (gene symbols), e.g. MSigDB *.symbols.gmt
#   --rank_by   limma column used to rank genes: t (default) or logFC
#
# Output (in --outdir):
#   <contrast>.fgsea_de.tsv     one row per tested gene set, sorted by p-value
#   <contrast>.fgsea_ranks.tsv  the ranked gene list that was tested (gene symbol, statistic)
#
# Genes without a gene symbol are dropped; when several Ensembl ids share a symbol the one with
# the largest |statistic| is kept.

suppressPackageStartupMessages({
    library(fgsea)
    library(data.table)
})

# ---- arguments -----------------------------------------------------------------------------
opts <- list(results = NULL, gene_map = NULL, gmt = NULL, contrast = NULL, rank_by = "t",
             min_size = "10", max_size = "500", fdr = "0.05", seed = "42", threads = "1", outdir = ".")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2 != 0) stop("Arguments must be given as --name value pairs")
for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[i])
    if (!key %in% names(opts)) stop("Unknown option: ", args[i])
    opts[[key]] <- args[i + 1]
}
for (req in c("results", "gene_map", "gmt", "contrast")) {
    if (is.null(opts[[req]])) stop("Missing required option --", req)
}
min_size <- as.integer(opts$min_size)
max_size <- as.integer(opts$max_size)
fdr      <- as.numeric(opts$fdr)
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

# ---- ranked gene list ----------------------------------------------------------------------
res <- fread(opts$results, data.table = FALSE)
if (!opts$rank_by %in% colnames(res)) stop("--rank_by column '", opts$rank_by, "' not in ", opts$results)
map <- unique(fread(opts$gene_map, select = c("gene_id", "gene_name"), data.table = FALSE))
map <- map[!duplicated(map$gene_id), ]

res$gene_name <- map$gene_name[match(res$gene_id, map$gene_id)]
res$stat      <- res[[opts$rank_by]]
res <- res[!is.na(res$gene_name) & nzchar(res$gene_name) & res$gene_name != res$gene_id & !is.na(res$stat), ]
res <- res[order(-abs(res$stat)), ]
res <- res[!duplicated(res$gene_name), ]               # keep the strongest id per symbol
ranks <- sort(setNames(res$stat, res$gene_name), decreasing = TRUE)
message(sprintf("Ranked %d genes by limma %s", length(ranks), opts$rank_by))

fwrite(data.frame(gene = names(ranks), stat = unname(ranks)),
       file.path(opts$outdir, paste0(opts$contrast, ".fgsea_ranks.tsv")), sep = "\t")

# ---- fgsea ---------------------------------------------------------------------------------
set.seed(as.integer(opts$seed))
fg <- fgsea(pathways = gs$sets, stats = ranks, minSize = min_size, maxSize = max_size,
            eps = 0, nproc = as.integer(opts$threads))

out <- data.frame(
    contrast         = opts$contrast,
    pathway          = fg$pathway,
    collection       = unname(gs$source[fg$pathway]),
    size             = fg$size,
    ES               = fg$ES,
    NES              = fg$NES,
    pval             = fg$pval,
    padj             = fg$padj,
    log2err          = fg$log2err,
    direction        = ifelse(!is.na(fg$padj) & fg$padj < fdr, ifelse(fg$NES > 0, "Up", "Down"), "NS"),
    leading_edge_n   = lengths(fg$leadingEdge),
    leading_edge     = vapply(fg$leadingEdge, paste, character(1), collapse = ";"),
    stringsAsFactors = FALSE
)
out <- out[order(out$pval, -abs(out$NES)), ]
fwrite(out, file.path(opts$outdir, paste0(opts$contrast, ".fgsea_de.tsv")), sep = "\t")

message(sprintf("%s: %d gene sets tested, %d Up and %d Down at FDR < %s",
                opts$contrast, nrow(out), sum(out$direction == "Up"), sum(out$direction == "Down"), fdr))
