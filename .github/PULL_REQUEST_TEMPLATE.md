<!--
# rnaseq-splice pull request

Thanks for contributing! Please describe the change and why it is needed, and tick the
checklist items that apply (delete the rest). See docs/CONTRIBUTING.md for conventions.
-->

## PR checklist

- [ ] This description explains the change and the reason for it.
- [ ] The test profile passes: `nextflow run . -profile test_rmats,docker --outdir <OUTDIR>`.
- [ ] Formatting passes (`prek run --all-files` or `pre-commit run --all-files`), or the Lint check is green.
- [ ] New options are opt-in and added to `nextflow.config` and `nextflow_schema.json`.
- [ ] `docs/splice_usage.md` and `docs/splice_output.md` are updated.
- [ ] `CHANGELOG.md` is updated.
- [ ] `README.md` and `CITATIONS.md` are updated for new tools.
