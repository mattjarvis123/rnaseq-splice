# kegg_overlay.R - colour KEGG pathway maps with log2FC or dPSI values (rnaseq-splice report)
#
# Uses the pathview Bioconductor package: Luo W, Brouwer C. Pathview: an R/Bioconductor package for
# pathway-based data integration and visualization. Bioinformatics. 2013;29(14):1830-1831.
# doi:10.1093/bioinformatics/btt285. KEGG maps: Kanehisa M, Goto S. Nucleic Acids Res. 2000;28(1):27-30.
#
# Reworked from Matt Jarvis's pathview.R (R_Functions/pathview.R, pathPlot()): KEGG's KGML and
# PNG files are fetched from the KEGG REST API, gene nodes are mapped and coloured with pathview's
# node.map() / node.color() / render.kegg.node(), and the colour key and a title are drawn on top.
# Changes for the report:
#   - KGML/PNG files are kept in a cache folder and downloaded only when missing, so reruns
#     (and runs without internet, once the cache is filled) do not contact KEGG;
#   - output is a PNG at the map's native size (for embedding in the HTML report);
#   - genes come in as Ensembl ids and are mapped to Entrez ids with org.Hs.eg.db.
#
# Note: KEGG's REST API and map images are free for academic use; commercial use needs a KEGG
# licence (https://www.pathway.jp).

suppressPackageStartupMessages({
    library(pathview)
    library(png)
})

# Make sure <pathway_id>.xml and <pathway_id>.png are in cache_dir; returns TRUE when both are there
kegg_fetch <- function(pathway_id, cache_dir) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    files <- c(xml = file.path(cache_dir, paste0(pathway_id, ".xml")),
               png = file.path(cache_dir, paste0(pathway_id, ".png")))
    urls  <- c(xml = paste0("https://rest.kegg.jp/get/", pathway_id, "/kgml"),
               png = paste0("https://rest.kegg.jp/get/", pathway_id, "/image"))
    for (k in names(files)) {
        if (file.exists(files[[k]]) && file.size(files[[k]]) > 0) next
        tmp <- paste0(files[[k]], ".part")
        ok <- tryCatch({
            utils::download.file(urls[[k]], tmp, mode = "wb", quiet = TRUE)
            TRUE
        }, error = function(e) FALSE, warning = function(w) FALSE)
        if (ok && file.exists(tmp) && file.size(tmp) > 0) file.rename(tmp, files[[k]]) else unlink(tmp)
        Sys.sleep(0.35)  # KEGG asks for no more than ~3 requests per second
    }
    all(file.exists(files) & file.size(files) > 0)
}

# Ensembl gene ids (version suffix allowed) -> Entrez ids; the first Entrez id per Ensembl id
ensembl_to_entrez <- function(ensembl) {
    ens <- sub("\\.[0-9]+$", "", ensembl)
    m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = unique(ens),
                                                columns = "ENTREZID", keytype = "ENSEMBL"))
    m <- m[!is.na(m$ENTREZID) & !duplicated(m$ENSEMBL), ]
    m$ENTREZID[match(ens, m$ENSEMBL)]
}

# Draw one coloured KEGG map.
#   values      named numeric vector, names = Entrez ids (several genes on one node are summed,
#               as in pathPlot(); use "max.abs" to keep the strongest instead)
#   limit       colour scale runs from -limit (low colour) to +limit (high colour); a function is
#               called with the values of the genes on this map and returns the limit
#   title       text drawn at the top of the map
#   out_png     output file
# Returns c(genes = genes on the map, with_value = those with a value), or NULL when the KEGG
# files are missing and cannot be downloaded.
kegg_map_png <- function(values, pathway_id, cache_dir, out_png, title = pathway_id, limit = 1,
                         node_sum = "max.abs", bins = 10,
                         low = "cornflowerblue", mid = "gray", high = "orange2", na_col = "transparent",
                         key_label = NULL) {
    if (!kegg_fetch(pathway_id, cache_dir)) return(NULL)
    xml_file <- file.path(cache_dir, paste0(pathway_id, ".xml"))
    png_file <- file.path(cache_dir, paste0(pathway_id, ".png"))

    # node naming and position info for each gene in the pathway (pathview::node.info)
    node_data <- node.info(xml_file)
    keep <- node_data$type %in% c("gene", "enzyme", "compound", "ortholog") &
        !is.na(node_data$x + node_data$y + node_data$width + node_data$height)
    node_data <- lapply(node_data, "[", keep)

    map_genes <- unique(unlist(node_data$kegg.names[node_data$type == "gene"]))
    if (is.function(limit)) limit <- limit(values[names(values) %in% map_genes])
    plot_data <- node.map(mol.data = values, node.data = node_data, node.types = "gene",
                          node.sum = node_sum, entrez.gnodes = 1)
    cols <- node.color(plot_data, limit = limit, bins = bins, both.dirs = TRUE, trans.fun = NULL,
                       discrete = FALSE, low = low, mid = mid, high = high, na.col = na_col)

    img    <- readPNG(png_file)
    width  <- ncol(img)
    height <- nrow(img)
    if (length(dim(img)) == 3 && dim(img)[3] == 4) img <- img[, , 1:3]   # drop an alpha channel
    if (length(dim(img)) == 2) img <- array(img, c(dim(img), 3))          # grey-scale image

    # leave room above the map for the title
    top_pad <- 60
    png(out_png, width = width, height = height + top_pad, res = 96)
    op <- par(mar = c(0, 0, 0, 0))
    plot(c(0, width), c(0, height + top_pad), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes = FALSE)
    img <- pathview:::render.kegg.node(plot.data = plot_data, cols.ts = cbind(cols)[, 1], img = img,
                                       same.layer = TRUE, type = "gene")
    rasterImage(img, 0, 0, width, height, interpolate = FALSE)
    col.key(limit = limit, bins = bins, both.dirs = TRUE, discrete = FALSE,
            graph.size = c(width = width, height = height), node.size = c(46, 17),
            key.pos = "topright", cex = 0.8, lwd = 1, low = low, mid = mid, high = high, align = "n")
    if (!is.null(key_label)) text(width - 10, height - 2, key_label, adj = c(1, 1), cex = 0.8)
    text(8, height + top_pad / 2, title, adj = c(0, 0.5), cex = 1, font = 2)
    par(op)
    invisible(dev.off())
    c(genes = length(map_genes), with_value = sum(map_genes %in% names(values)), limit = limit)
}
