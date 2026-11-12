# Building the MEDuSA documentation

The teaching sequence comprises two standard package vignettes and a Quarto technical companion:

| Source | HTML output | Role |
|---|---|---|
| `vignettes/medusa-concepts.Rmd` | `medusa-concepts.html` | Conceptual foundations and literature (updated through PR #66) |
| `vignettes/medusa-workflows.Rmd` | `medusa-workflows.html` | Executable example and study recipes (updated through PR #66) |
| `medusa_deconvolution_review_prs_42_53.qmd` | `medusa_deconvolution_review_prs_42_53.html` | Historical technical review for PRs #42–#55+ |

Run the commands below from the repository root. Install the package dependencies first. Rendering requires R, knitr/rmarkdown, Pandoc (also bundled with Quarto), and the Quarto CLI for the technical companion. The current source must be installed for a Quarto render, which executes R in its own process. Merely loading a development checkout in another R process does not update that installation.

A temporary package library avoids changing a regular installation. It should contain no older build of `fmri.pipeline`; dependencies can remain in the normal R libraries. On systems with a different temporary directory, choose equivalent paths.

```sh
mkdir -p /tmp/medusa-lib local/medusa-docs
R CMD INSTALL --library=/tmp/medusa-lib .
R_LIBS=/tmp/medusa-lib Rscript -e 'rmarkdown::render("vignettes/medusa-concepts.Rmd", output_dir="local/medusa-docs")'
R_LIBS=/tmp/medusa-lib Rscript -e 'rmarkdown::render("vignettes/medusa-workflows.Rmd", output_dir="local/medusa-docs")'
R_LIBS=/tmp/medusa-lib quarto render medusa_deconvolution_review_prs_42_53.qmd --to html --output-dir local/medusa-docs
```

For both the R Markdown and Quarto commands, `output_dir` / `--output-dir` is resolved from the repository root. Keep all three HTML files together in `local/medusa-docs/` for their relative navigation links to work. `local/` is excluded from package builds and git tracking. Open `local/medusa-docs/medusa-concepts.html` to begin.

The two Rmd files have standard vignette engine/index metadata and are included by the existing knitr build machinery. A normal package build with vignettes enabled or pkgdown article build can publish them. The QMD lives at the repository root (and is listed in `.Rbuildignore`) as a standalone technical companion rather than requiring the Quarto CLI on pkgdown CI runners. No rendering command in the documents installs dependencies or compiles native code automatically.

The worked vignette executes the small simulation, HRF estimation, continuous deconvolution, impulse and duration epoch regressions, fixed-HRF comparison, and trial-identity checks. Study-specific NIfTI/file/group recipes are deliberately marked `eval=FALSE`. The technical companion executes its simulation chunks, except that the historical speed figure only plots a legacy hard-coded table. Use `inst/benchmarks/` for actual runtime comparisons and record the source revision, hardware, thread count, and session information.

## Alternative HTML renderer

If Quarto cannot write its user-level theme cache (for example, under a restricted execution environment), the technical companion can also be executed and rendered with rmarkdown. This produces the same R demonstrations in a standard R Markdown HTML layout, with the same output filename and relative links:

```sh
R_LIBS=/tmp/medusa-lib Rscript -e 'rmarkdown::render("medusa_deconvolution_review_prs_42_53.qmd", output_format=rmarkdown::html_document(toc=TRUE, toc_depth=3, number_sections=TRUE, code_folding="hide", fig_width=10.5, fig_height=5.5), output_dir="local/medusa-docs")'
```

This alternative needs no Quarto R package or CLI. It does not change the QMD source or install anything. The checked local previews use this fallback for the companion and the standard renderer for both package vignettes.
