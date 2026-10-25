#!/usr/bin/env Rscript
#
# rmats_filter.R - read-coverage filter for rMATS results (rnaseq-splice)
#
# Keeps an rMATS event only if every replicate of both conditions has at least --min_reads
# junction reads supporting it (inclusion + skipping counts, IJC + SJC). Events below that
# cannot give a reliable PSI, but rMATS still tests them and they inflate the list of
# significant events with large, noisy dPSI values.
#
# rMATS computes FDR with Benjamini-Hochberg over all events of a type, so after filtering the
# FDR is recomputed (BH over the retained events of the file). This is independent filtering:
# coverage does not depend on which condition the reads came from, so dropping low-coverage
# events before the correction is valid and gives the retained events more power. rMATS's own
# FDR is kept in the column FDR_unfiltered.
#
# Two columns are appended to each table; every rMATS column is kept unchanged otherwise:
#   min_replicate_reads  smallest IJC + SJC of any replicate (the value the filter tests)
#   FDR_unfiltered       rMATS's FDR, computed over all events before filtering
#
# Usage:
#   rmats_filter.R --mats SE.MATS.JC.txt,SE.MATS.JCEC.txt,... [--min_reads 10] [--outdir filtered]
#
#   --mats       comma-separated rMATS <TYPE>.MATS.{JC,JCEC}.txt files
#   --min_reads  minimum IJC + SJC in every replicate (0 keeps every event and leaves FDR as is)
#   --outdir     where the filtered tables (same file names) and filter_summary.tsv are written
#
# Base R only, so it runs in the rMATS container.

opts <- list(mats = NULL, min_reads = "10", outdir = "filtered")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2 != 0) stop("Arguments must be given as --name value pairs")
for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[i])
    if (!key %in% names(opts)) stop("Unknown option: ", args[i])
    opts[[key]] <- args[i + 1]
}
if (is.null(opts$mats)) stop("Missing required option --mats")
min_reads <- as.numeric(opts$min_reads)
if (is.na(min_reads) || min_reads < 0) stop("--min_reads must be a number >= 0")
dir.create(opts$outdir, showWarnings = FALSE, recursive = TRUE)

# Smallest per-replicate IJC + SJC of each event. rMATS writes one comma-separated count per
# replicate; a condition with no samples (single-condition runs) has an empty field and is ignored.
group_min <- function(ijc, sjc) {
    i <- strsplit(as.character(ijc), ",", fixed = TRUE)
    s <- strsplit(as.character(sjc), ",", fixed = TRUE)
    mapply(function(a, b) {
        a <- suppressWarnings(as.numeric(a)); b <- suppressWarnings(as.numeric(b))
        if (!length(a) && !length(b)) return(Inf)              # empty condition
        if (length(a) != length(b) || anyNA(a) || anyNA(b)) return(NA_real_)
        min(a + b)
    }, i, s, USE.NAMES = FALSE)
}

files <- trimws(strsplit(opts$mats, ",", fixed = TRUE)[[1]])
summary_rows <- list()
for (f in files) {
    # quote = "" keeps rMATS's quoted GeneID/geneSymbol strings exactly as written
    d <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE, quote = "",
                    colClasses = "character", na.strings = character())
    out_file <- file.path(opts$outdir, basename(f))
    n_all <- nrow(d)
    if (n_all == 0) {
        hdr <- readLines(f, n = 1)
        writeLines(paste(hdr, "min_replicate_reads", "FDR_unfiltered", sep = "\t"), out_file)
        summary_rows[[f]] <- data.frame(file = basename(f), events = 0, kept = 0, removed = 0,
                                        significant_fdr05_before = 0, significant_fdr05_after = 0)
        next
    }
    cov <- pmin(group_min(d$IJC_SAMPLE_1, d$SJC_SAMPLE_1), group_min(d$IJC_SAMPLE_2, d$SJC_SAMPLE_2))
    cov[is.infinite(cov)] <- NA
    keep <- if (min_reads > 0) !is.na(cov) & cov >= min_reads else rep(TRUE, n_all)

    p   <- suppressWarnings(as.numeric(d$PValue))
    fdr <- suppressWarnings(as.numeric(d$FDR))
    d$min_replicate_reads <- ifelse(is.na(cov), "NA", format(cov, scientific = FALSE, trim = TRUE))
    d$FDR_unfiltered <- d$FDR
    d <- d[keep, , drop = FALSE]

    # --statoff (single-condition) runs have no p-values; their FDR column is left as rMATS wrote it
    p_kept <- p[keep]
    new_fdr <- fdr[keep]
    if (min_reads > 0 && any(!is.na(p_kept))) {
        new_fdr <- p.adjust(p_kept, method = "BH")
        d$FDR <- ifelse(is.na(new_fdr), "NA", as.character(new_fdr))
    }
    write.table(d, out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

    summary_rows[[f]] <- data.frame(
        file = basename(f), events = n_all, kept = sum(keep), removed = n_all - sum(keep),
        significant_fdr05_before = sum(!is.na(fdr) & fdr < 0.05),
        significant_fdr05_after  = sum(!is.na(new_fdr) & new_fdr < 0.05))
}

s <- do.call(rbind, summary_rows)
s$min_reads <- min_reads
write.table(s, file.path(opts$outdir, "filter_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (r in seq_len(nrow(s))) {
    message(sprintf("%-20s %7d events, %7d kept (IJC + SJC >= %g in every replicate); FDR < 0.05: %d before, %d after",
                    s$file[r], s$events[r], s$kept[r], min_reads, s$significant_fdr05_before[r], s$significant_fdr05_after[r]))
}
