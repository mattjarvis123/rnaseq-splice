# rnaseq-splice: Usage

This page covers what rnaseq-splice adds to nf-core/rnaseq: differential splicing, differential expression, gene set
enrichment and the results report. Everything else (read QC, trimming, alignment, quantification, reference genomes and
the remaining parameters) works as in nf-core/rnaseq 3.27.0; see the
[nf-core/rnaseq usage documentation](usage.md) and the
[parameter reference](https://nf-co.re/rnaseq/3.27.0/parameters).

- [Samplesheet](#samplesheet)
- [Contrasts sheet](#contrasts-sheet)
- [Differential splicing (rMATS)](#differential-splicing-rmats)
- [Differential expression (limma/voom)](#differential-expression-limmavoom)
- [Gene set enrichment (fgsea)](#gene-set-enrichment-fgsea)
- [Results report](#results-report)
- [Command-line booleans](#command-line-booleans)
- [Running on a cluster or shared storage](#running-on-a-cluster-or-shared-storage)
- [Test profile](#test-profile)

## Samplesheet

The samplesheet is the standard nf-core/rnaseq samplesheet with one extra column, `condition`, which assigns each
sample to an experimental group:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,strandedness,condition
DMSO_REP1,dmso1_R1.fastq.gz,dmso1_R2.fastq.gz,auto,DMSO
DMSO_REP2,dmso2_R1.fastq.gz,dmso2_R2.fastq.gz,auto,DMSO
DMSO_REP3,dmso3_R1.fastq.gz,dmso3_R2.fastq.gz,auto,DMSO
DRUG_REP1,drug1_R1.fastq.gz,drug1_R2.fastq.gz,auto,DRUG_1uM
DRUG_REP2,drug2_R1.fastq.gz,drug2_R2.fastq.gz,auto,DRUG_1uM
DRUG_REP3,drug3_R1.fastq.gz,drug3_R2.fastq.gz,auto,DRUG_1uM
```

- `condition` is required when `--rmats` or `--limma` is on and is ignored otherwise.
- Condition names must be valid R names: start with a letter, then use only letters, digits, `.` and `_`
  (`DRUG_1uM` is fine, `1uM` and `drug-1uM` are not).
- The order of samples within a condition matters only for `--rmats_paired_stats` (see below) and for the per-sample
  PSI columns in the rMATS tables, which follow samplesheet order.
- `strandedness: auto` is recommended. The strandedness inferred by nf-core/rnaseq is passed on to rMATS.

## Contrasts sheet

`--contrasts` lists the comparisons to make. Each row produces both differential expression and differential splicing
results:

```csv title="contrasts.csv"
contrast,treatment,control
DRUG_1uM_vs_DMSO,DRUG_1uM,DMSO
DRUG_0.1uM_vs_DMSO,DRUG_0.1uM,DMSO
```

- `treatment` and `control` must match values of the samplesheet `condition` column.
- log2 fold changes and ΔPSI are always **treatment minus control**: a positive log2FC means higher expression in the
  treatment, a positive ΔPSI means more inclusion in the treatment.
- The report shows contrasts in the order of this file, so list doses or time points in increasing order.

## Differential splicing (rMATS)

Switch on with `--rmats`. rMATS needs genome alignments, so use a genome aligner: the default `star_salmon`,
`star_rsem`, `hisat2` or `bowtie2_salmon`.

**How it runs.** For every sample, `rmats.py --task prep` processes the final genome BAM (after the percent-mapped
filter and Picard duplicate marking). Then for every contrast, `rmats.py --task post` combines the treatment samples
(`--b1`) and control samples (`--b2`) and tests five event types: skipped exon (SE), mutually exclusive exons (MXE),
alternative 3' and 5' splice sites (A3SS, A5SS) and retained intron (RI). The strand setting (`--libType`) comes from
the inferred strandedness, and single- or paired-end mode from the samplesheet.

**Settings:**

- `--rmats_read_len` (default `40`): set this to your read length (for example `100` or `150`). `--variable-read-length`
  is always set, so trimmed reads are still counted.
- `--rmats_splice_diff_cutoff` (default `0.0001`): the ΔPSI threshold in rMATS's null hypothesis (`--cstat`).
- `--rmats_paired_stats`: use the paired model (PAIRADISE) for matched designs, for example the same donors or cell
  passages treated and untreated. The n-th treatment sample is paired with the n-th control sample in samplesheet order,
  so list them in matching order.
- `--rmats_novel_splice_site`: also detect unannotated splice sites, with `--rmats_min_intron_len` (50) and
  `--rmats_max_exon_len` (500).

**Coverage filter.** rMATS tests every event it can find, including events supported by one or two reads, which adds
noise and inflates the number of tests. With `--rmats_min_reads` (default `10`), an event is kept only if every replicate
of both conditions has at least that many junction reads (inclusion plus skipping), and the Benjamini-Hochberg FDR is
then recomputed over the kept events. rMATS's own FDR is kept in an `FDR_unfiltered` column. The filtered tables feed
gene set enrichment and the report; the unfiltered rMATS output is published as well. `--rmats_min_reads 0` turns the
filter off.

The report and gene set enrichment use junction counts only (`*.MATS.JC.txt`). The `JCEC` tables, which also count
reads on the exon bodies, are published for reference.

## Differential expression (limma/voom)

Switch on with `--limma`. It needs gene-level quantification from `star_salmon`, `star_rsem` or a pseudo-aligner.

**Input.** limma uses the merged length-scaled gene counts (`salmon.merged.gene_counts_length_scaled.tsv`, from tximport
`lengthScaledTPM`, or the RSEM or pseudo-aligner equivalent). These are the recommended input when Salmon estimates are
used with limma/voom.

**Steps:**

1. **Low-count filter** (nf-core `custom/matrixfilter`): keep a gene if it has at least `--limma_min_count` (default 10)
   counts in at least as many samples as the smallest condition group.
2. **Model** (nf-core `limma/differential`): design `~ 0 + condition` fitted on all samples together, TMM normalisation,
   voom precision weights, `lmFit`, then `eBayes` (with `robust = TRUE` if `--limma_robust` is set).
3. **Results:** one `topTable` per contrast, with Benjamini-Hochberg adjusted p-values.

**Current limits.** The design has no batch, donor or pairing terms, and extra samplesheet columns are not passed to
limma. Fitting all samples together shares variance information across conditions, which suits designs with similar
variability between groups.

## Gene set enrichment (fgsea)

Switch on with `--gsea --gsea_gmt <files>`. It needs `--limma` and/or `--rmats`. `--gsea_gmt` takes one or more
comma-separated GMT files keyed by **gene symbol**, for example from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb):

```yaml
gsea_gmt: h.all.v2026.2.Hs.symbols.gmt,c2.cp.reactome.v2026.2.Hs.symbols.gmt,c2.cp.kegg_legacy.v2026.2.Hs.symbols.gmt
```

Gene sets with 10–500 genes are tested (`--gsea_min_size`, `--gsea_max_size`).

**Expression.** Preranked GSEA (multilevel fgsea) of every tested gene, ranked by the limma moderated t statistic
(`--gsea_rank_by logFC` ranks by fold change instead). Ensembl ids are mapped to symbols with the `gene_name` column of
the merged counts; when several ids share a symbol, the one with the largest |statistic| is kept. A gene set is **Up**
(NES > 0) or **Down** (NES < 0) when its adjusted p-value is below `--gsea_fdr` (0.05).

**Splicing.** Over-representation analysis (`fgsea::fora`) of genes with at least one significant rMATS event
(FDR < `--gsea_splicing_fdr`, 0.01; |ΔPSI| > `--gsea_splicing_dpsi`, 0), run separately for genes with higher and lower
inclusion in the treatment. The background is every gene with at least one tested event, so the test asks which gene
sets are over-represented among spliced genes relative to genes where splicing could be measured.

With no ΔPSI cut-off, a strong splicing perturbation can make most genes with an event significant, and enrichment then
becomes modest. A cut-off of `--gsea_splicing_dpsi 0.1` to `0.2` gives more specific gene lists in that case.

## Results report

The report is rendered as the last step whenever `--limma` and/or `--rmats` is on (`--skip_report` turns it off). It is a
single HTML file, `<outdir>/report/<report_name>.html`, that opens in any browser without a server, with the full tables
next to it in `<report_name>_tables/`. Keep the two together when you move or share the report: the KEGG map images are
loaded from the tables folder.

**Significance cut-offs.** The pipeline never filters by significance. The report uses its own cut-offs to call
differentially expressed genes and spliced events: `--report_de_fdr` (0.01), `--report_de_lfc` (0),
`--report_splicing_fdr` (0.01) and `--report_splicing_dpsi` (0). The defaults keep small changes, which suits gene set
analyses; raise the effect-size cut-offs for shorter, stronger gene lists.

**KEGG pathway maps.** When the GSEA input includes the MSigDB KEGG legacy collection (`c2.cp.kegg_legacy`), the
report draws KEGG pathway maps for the significant KEGG gene sets, coloured by log2FC (expression) or ΔPSI (splicing),
using [pathview](https://bioconductor.org/packages/pathview/). The map files are downloaded from the KEGG REST API
(`rest.kegg.jp`). Set `--report_kegg_cache <folder>` to keep them between runs, which also lets runs without internet
access reuse them; `assets/report/fetch_kegg_cache.sh <folder>` fills a cache with all KEGG legacy maps in advance. Turn
the maps off with `--report_kegg_maps false`.

> [!NOTE]
> KEGG's API and pathway maps are free for academic use; non-academic use requires a licence from
> [Pathway Solutions](https://www.pathway.jp/) (see the [KEGG terms](https://www.kegg.jp/kegg/legal.html)). Check your
> organisation's status before enabling the maps, and cite KEGG when you use them.

**Other options:** `--report_name`, `--report_title`, `--report_project` (shown on the overview page) and
`--report_table_max_rows` (5000 rows per table in the report; the full tables are always written to the tables folder).

**Size.** The report embeds its data, so it grows with the number of genes, events and contrasts (about 9 MB for the
test profile and 20–30 MB for a human experiment with a few contrasts). Rendering a large report needs several GB of
memory for pandoc. The report re-renders on every run, including `-resume` runs, because it records the software
versions of that run.

## Command-line booleans

Under Nextflow's strict syntax (the default from Nextflow 26.04), command-line values reach the pipeline as strings.
rnaseq-splice converts its own options, so `--rmats`, `--limma true` and `--rmats false` all behave as expected. The
upstream nf-core/rnaseq code does not, so an explicit `--skip_<step> false` on the command line is read as true. Put
explicit `false` values in a parameters file instead:

```yaml title="params.yaml"
skip_markduplicates: false
```

```bash
nextflow run mattjarvis123/rnaseq-splice -profile docker -params-file params.yaml
```

A parameters file is the recommended way to run the pipeline anyway: it records the settings of a run in one place.

## Running on a cluster or shared storage

The pipeline runs anywhere Nextflow does. Where Docker is not available, use `-profile apptainer` (or `singularity`),
and give Nextflow a config describing your machine or scheduler with `-c`. Templates are in
[`conf/examples/`](../conf/examples/):

| File                      | For                                                                     |
| ------------------------- | ----------------------------------------------------------------------- |
| `local_docker.config`     | A laptop or workstation with Docker.                                    |
| `server_apptainer.config` | One shared Linux server without a scheduler; tasks run on the server.   |
| `slurm_apptainer.config`  | A SLURM cluster; every task is submitted as a SLURM job.                |
| `pbs_apptainer.config`    | A PBS Pro or Torque/OpenPBS cluster; every task is submitted as a PBS job. |
| `submit_slurm.sh`, `submit_pbs.sh` | Job scripts that start the pipeline itself on a SLURM or PBS cluster. |

Copy the file you need, replace the paths, queue or partition, account and resource limits, and pass it with `-c`.
If your institution has a profile in [nf-core/configs](https://github.com/nf-core/configs), you can use it instead:
`-profile <institution>` loads it as for any nf-core pipeline.

### SLURM and PBS clusters

Nextflow submits every task as its own job, requesting the CPUs, memory and time of the task's resource label
(`process_low`, `process_medium`, `process_high` and so on, set in `conf/base.config` and capped by
`process.resourceLimits`). The pipeline itself runs as a small, long "head" job:

```bash
sbatch submit_slurm.sh   # SLURM
qsub submit_pbs.sh       # PBS Pro or Torque/OpenPBS
```

or from a login node inside `screen` or `tmux`, if your cluster allows long-running processes there. Points to check:

- **Shared filesystem.** The work folder (`-w`), the container cache and the inputs must be on storage that every
  compute node can see, such as cluster scratch, Lustre, GPFS or NFS.
- **Containers.** Nextflow downloads each container image once, from the node running the head job, into the shared
  cache (`apptainer.cacheDir` or `NXF_APPTAINER_CACHEDIR`); compute nodes then read the images from the cache. The head
  job's node therefore needs internet access, or the cache must be filled beforehand from one that has it.
- **Limits.** Keep `process.resourceLimits` within what your partition or queue allows: a task that asks for more than
  the queue maximum stays pending. A human STAR index needs about 32–38 GB of memory.
- **KEGG maps.** The report downloads KEGG map files when they are not in `--report_kegg_cache`. Compute nodes often
  have no internet access, so fill a shared cache first with `assets/report/fetch_kegg_cache.sh <folder>` from a node
  that does, and pass that folder. Missing maps are listed in the report rather than failing the run.
- **Scheduler friendliness.** `executor.queueSize` and `executor.submitRateLimit` in the templates limit how many jobs
  are queued and how fast they are submitted.

### A single server or network storage

- **Resources.** Set `process.resourceLimits` (largest single task) and `executor.cpus` / `executor.memory` (everything
  running at once) to your share of the machine.
- **Disk space.** The work folder needs roughly 3–5 times the size of the FASTQ files. Delete it after a successful run
  once you no longer need `-resume`.
- **Network storage (NAS, SMB/CIFS).** Nextflow creates symbolic links in every task folder. An SMB/CIFS mount without
  the `mfsymlinks` option cannot hold them, so either mount with `mfsymlinks` or put the work folder (`-w`) on local disk.
  Keep the launch folder, where Nextflow writes `.nextflow/` and its logs, on local disk too: Nextflow's run history
  relies on file locking, which is unreliable over SMB. Input files and `--outdir` can stay on the network share.
- **Container cache.** Set `apptainer.cacheDir` (or `NXF_APPTAINER_CACHEDIR`) to a local folder that is kept between
  runs, so container images are downloaded once.
- **Long runs.** Start the pipeline inside `screen` or `tmux` so it survives a lost connection.

## Test profile

`-profile test_rmats` runs the nf-core/rnaseq test data (yeast, five samples, single- and paired-end) with three
conditions and two contrasts, with rMATS and limma switched on:

```bash
nextflow run mattjarvis123/rnaseq-splice -profile test_rmats,docker --outdir test_results
```

The test checks that the steps connect and produce output; the design (one sample in one condition, mixed read types)
is not meant to give biologically meaningful results. Gene set enrichment is not part of the test profile, because it
needs gene set files.
