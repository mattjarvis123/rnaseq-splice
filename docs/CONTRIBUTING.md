# rnaseq-splice: Contributing

Thanks for your interest in improving rnaseq-splice. Bug reports and feature requests are welcome as
[GitHub issues](https://github.com/mattjarvis123/rnaseq-splice/issues). For changes to the steps inherited from
nf-core/rnaseq (read QC, alignment, quantification), please consider contributing to
[nf-core/rnaseq](https://github.com/nf-core/rnaseq) directly; improvements there reach this pipeline with the next merge.

## Making a change

1. Fork the repository and create a branch from `main`.
2. Make your change, following the conventions below.
3. Run the test profile and the formatting checks:

   ```bash
   nextflow run . -profile test_rmats,docker --outdir test_results
   prek run --all-files   # or: pre-commit run --all-files
   ```

4. Update the documentation (`docs/splice_usage.md`, `docs/splice_output.md`), `CHANGELOG.md`, and `README.md` and
   `CITATIONS.md` for new tools.
5. Open a pull request. The Pipeline test and Lint checks run automatically.

## Conventions

- **nf-core layout and style.** DSL2 modules and subworkflows with `meta` maps, tool options in `ext.args` /
  `ext.prefix` in `conf/modules/*.config`, and one `publishDir` per process.
- **Opt-in features.** New analyses get a parameter that defaults to `false`, so the pipeline still runs as
  nf-core/rnaseq without them. Add every new parameter to `nextflow.config` and `nextflow_schema.json`.
- **Keep fork code in separate files.** Put new code in `modules/local/`, `subworkflows/local/`, `bin/` and
  `conf/modules/`, and keep edits to upstream files (`workflows/rnaseq/main.nf`, `nextflow.config`,
  `nextflow_schema.json`, `assets/schema_input.json`, `modules.json`) to small blocks, so new nf-core/rnaseq releases
  merge cleanly.
- **Containers for every process.** Use BioContainers or Seqera community images, pinned to a version.
- **No `eval()` version outputs in new local modules.** Nextflow numbers `eval` outputs across the whole pipeline, so
  adding one renumbers other processes' scripts and invalidates their `-resume` cache. Report versions another way.
- **Credit code taken from elsewhere** with its source and version in comments, in the README and in `CITATIONS.md`.
- **Formatting.** Prettier (version pinned in `.pre-commit-config.yaml`) with the settings in `.prettierrc.yml`:
  print width 120, 4-space indent for JSON and 2 for Markdown and YAML.

## Merging a new nf-core/rnaseq release

```bash
git remote add upstream https://github.com/nf-core/rnaseq.git   # once
git fetch upstream --tags
git checkout -b merge-rnaseq-3.28.0
git merge 3.28.0
```

Expect conflicts in the upstream files listed above. After merging, run the test profile, check the startup banner
and MultiQC texts still name rnaseq-splice, and note the new upstream version in `README.md`, `CHANGELOG.md`, the
manifest description and the startup message.
