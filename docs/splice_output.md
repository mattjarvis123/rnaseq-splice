# rnaseq-splice: Output

This page describes the output that rnaseq-splice adds to nf-core/rnaseq. All other folders (`fastqc/`, `trimgalore/`,
`star_salmon/` alignments and quantification, `multiqc/`, `pipeline_info/` and so on) are described in the
[nf-core/rnaseq output documentation](output.md).

Paths are relative to `--outdir`. `<aligner>` is the aligner folder, `star_salmon` by default, and `<contrast>` is a
contrast name from the contrasts sheet.

- [Differential splicing: rMATS](#differential-splicing-rmats)
- [Differential expression: limma/voom](#differential-expression-limmavoom)
- [Gene set enrichment: fgsea](#gene-set-enrichment-fgsea)
- [Results report](#results-report)

## Differential splicing: rMATS

<details markdown="1">
<summary>Output files</summary>

- `<aligner>/rmats/<contrast>/`
  - `{SE,MXE,A3SS,A5SS,RI}.MATS.JC.txt`: test results per event type, using junction reads only.
  - `{SE,MXE,A3SS,A5SS,RI}.MATS.JCEC.txt`: the same, also counting reads on the exon bodies.
  - `JC.raw.input.*.txt`, `JCEC.raw.input.*.txt`: read counts per event and sample, as given to the statistical model.
  - `fromGTF.*.txt`, `fromGTF.novelJunction.*.txt`, `fromGTF.novelSpliceSite.*.txt`: every event detected, with its
    coordinates, including novel junctions and splice sites.
  - `summary.txt`: number of events and significant events per type.
  - `filtered/`: the `*.MATS.JC.txt` and `*.MATS.JCEC.txt` tables after the read-coverage filter
    (`--rmats_min_reads`), with the FDR recomputed over the kept events, and `filter_summary.tsv` (events kept and
    removed, and significant events before and after filtering, per table).
- `<aligner>/rmats/<contrast>.log`: the rMATS post-step log.
- `<aligner>/rmats/prep/`: per-sample read outcomes from the rMATS prep step.

</details>

Each row of a `*.MATS.*.txt` table is one splicing event. The most useful columns are:

| Column                       | Meaning                                                                                              |
| ---------------------------- | ---------------------------------------------------------------------------------------------------- |
| `GeneID`, `geneSymbol`       | Gene of the event.                                                                                   |
| `chr`, `strand` and exon coordinates | Location of the event (0-based starts). The coordinate columns differ by event type.         |
| `IJC_SAMPLE_1`, `SJC_SAMPLE_1` | Inclusion and skipping junction counts in the treatment samples, comma-separated in samplesheet order. |
| `IJC_SAMPLE_2`, `SJC_SAMPLE_2` | The same for the control samples.                                                                  |
| `IncLevel1`, `IncLevel2`     | Per-sample PSI (percent spliced in, 0–1) in the treatment and control samples.                       |
| `IncLevelDifference`         | ΔPSI, mean treatment PSI minus mean control PSI.                                                     |
| `PValue`, `FDR`              | rMATS test p-value and FDR. In `filtered/`, `FDR` is recomputed over the kept events.                |
| `min_replicate_reads`        | `filtered/` only: smallest inclusion + skipping read count across all replicates.                    |
| `FDR_unfiltered`             | `filtered/` only: rMATS's original FDR.                                                              |

The rMATS documentation describes every column: [rMATS-turbo output](https://github.com/Xinglab/rmats-turbo#output).

## Differential expression: limma/voom

<details markdown="1">
<summary>Output files</summary>

- `<aligner>/limma/<contrast>/`
  - `<contrast>.limma.results.tsv`: results for every tested gene.
  - `<contrast>.normalised_counts.tsv`: voom log2-CPM values (TMM-normalised) per gene and sample.
  - `<contrast>.limma.mean_difference.png`: mean-difference (MA) plot.
  - `<contrast>.limma.model.txt`: the design formula and contrast that were fitted.
  - `<contrast>.MArrayLM.limma.rds`: the fitted limma model, to load in R with `readRDS()`.
  - `<contrast>.R_sessionInfo.log`: R and package versions.
- `<aligner>/limma/`
  - `gene_counts.filtered.tsv`: the gene count matrix after the low-count filter.
  - `gene_counts.tests.tsv`, `gene_counts.thresholds.tsv`: which genes passed the filter, and the thresholds used.

</details>

`<contrast>.limma.results.tsv` has one row per gene:

| Column      | Meaning                                                                   |
| ----------- | ------------------------------------------------------------------------- |
| `gene_id`   | Gene id from the GTF (Ensembl id for GENCODE or Ensembl annotation).      |
| `logFC`     | log2 fold change, treatment minus control.                                |
| `AveExpr`   | Average log2-CPM across all samples.                                      |
| `t`         | Moderated t statistic (used to rank genes for GSEA).                      |
| `P.Value`   | p-value.                                                                  |
| `adj.P.Val` | Benjamini-Hochberg adjusted p-value.                                      |
| `B`         | log-odds that the gene is differentially expressed.                       |

Gene symbols are in the merged Salmon gene counts (`<aligner>/salmon.merged.gene_counts.tsv`, column `gene_name`), and
in every table of the results report.

## Gene set enrichment: fgsea

<details markdown="1">
<summary>Output files</summary>

- `<aligner>/gsea/<contrast>/`
  - `<contrast>.fgsea_de.tsv`: preranked GSEA of the limma results, one row per gene set, sorted by p-value.
  - `<contrast>.fgsea_ranks.tsv`: the ranked gene list that was tested (gene symbol and ranking statistic).
  - `<contrast>.fgsea_splicing.tsv`: over-representation of gene sets among genes with higher-inclusion and
    lower-inclusion events, one row per gene set and gene list, sorted by p-value.
  - `<contrast>.fgsea_splicing_genes.tsv`: the gene lists and background (universe) that were tested.

</details>

The tables contain every tested gene set, with the GMT file it came from in `collection`.

- `fgsea_de.tsv`: `NES` is the normalised enrichment score (positive: higher in the treatment), `pval` and `padj` the
  p-value and BH-adjusted p-value, and `direction` is `Up`, `Down` or `NS` at `padj` < `--gsea_fdr` (0.05).
  `leading_edge` lists the genes driving the enrichment.
- `fgsea_splicing.tsv`: `gene_list` is `Higher inclusion` or `Lower inclusion`; `overlap` is the number of spliced genes
  in the set, `fold_enrichment` the ratio of observed to expected overlap, `significant` is `TRUE` at
  `padj` < `--gsea_fdr`, and `overlap_genes` lists the genes.

## Results report

<details markdown="1">
<summary>Output files</summary>

- `report/`
  - `<report_name>.html`: the results report (default `rnaseq_splice_report.html`).
  - `<report_name>_tables/`: the full tables behind the report, as tab-separated files:
    - `<contrast>_DE_genes.tsv`, `<contrast>_DS_events.tsv`: differentially expressed genes and spliced events at the
      report cut-offs, with gene symbols.
    - `across_contrasts_DE_genes.tsv`, `across_contrasts_DS_events.tsv`: genes and events side by side across contrasts.
    - `de_and_ds_genes.tsv`: genes that are both differentially expressed and differentially spliced.
    - `<contrast>_gsea_expression.tsv`, `<contrast>_gsea_splicing.tsv`, `gsea_expression_significant.tsv`,
      `gsea_splicing_significant.tsv`: gene set results.
    - `kegg/`: KEGG pathway map images (PNG), linked from the report.

</details>

The report is a single HTML file with these pages:

| Page                        | Content                                                                                                     |
| --------------------------- | ----------------------------------------------------------------------------------------------------------- |
| Overview                    | Samples, contrasts, cut-offs used, and the number of DE genes and spliced events per contrast.              |
| QC                          | Read counts, mapping and duplication rates from MultiQC.                                                    |
| Sample structure            | PCA and sample correlation of the voom log-CPM values.                                                      |
| DE summary                  | Up and down genes per contrast, overlaps between contrasts (UpSet and Euler plots).                         |
| DE by contrast              | Volcano and MA plots, a searchable results table and a heatmap of the top genes, for each contrast.          |
| Splicing summary            | Significant events per contrast and event type.                                                             |
| Splicing by contrast        | ΔPSI volcano plot, ΔPSI by event type, an event table and a PSI heatmap of the top events per type.         |
| DE × splicing               | log2FC against ΔPSI per gene, showing genes that are differentially expressed, spliced, both or neither.     |
| Across contrasts            | With two or more contrasts: up/down counts, genes shared between contrasts, and log2FC and ΔPSI heatmaps.   |
| Gene sets                   | GSEA summary across contrasts and, per contrast, the top enriched gene sets for expression and splicing.    |
| KEGG maps                   | KEGG pathway maps of significant KEGG gene sets, coloured by log2FC or ΔPSI.                                |
| Gene lookup                 | Type a gene symbol or Ensembl id to see its expression per sample, log2FC per contrast and its top events.  |
| Methods                     | Methods text with the settings and software versions of the run, ready to adapt for a publication.          |

Pages appear only when their input exists: for example the splicing pages need `--rmats` and the gene set pages need
`--gsea`.
