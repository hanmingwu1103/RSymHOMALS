# RSymHOMALS

`RSymHOMALS` is the reproducibility package for the manuscript
`Symbolic HOMALS: Visualization of Multi-Valued Categorical Data through Homogeneity Analysis`.
It provides code for fitting Symbolic HOMALS to multi-valued categorical data,
running the simulation study, reproducing the paper figures, and analyzing the
two real datasets used in Section 6.

The package supports three variants:

- Euclidean Symbolic HOMALS
- Wasserstein Symbolic HOMALS
- Hausdorff Symbolic HOMALS

It also includes a naive MCA benchmark for comparison in the real-data
analyses.

## Main functions

- `fit_symhomals()` fits one Symbolic HOMALS model.
- `analyze_symbolic_dataset()` compares multiple methods on one dataset.
- `analyze_symbolic_csv()` loads a set-valued CSV file and analyzes it.
- `simulate_symbolic_data()` generates synthetic multi-valued categorical data.
- `simulate_symhomals_study()` runs the Section 5 simulation study.
- `plot_symbolic_embedding()` draws the multi-panel object-category maps used in
  the paper.
- `sym_example_path()` locates packaged example files under `inst/extdata`.

## Repository layout

- `R/`: package source code
- `man/`: package manual pages
- `inst/extdata/`: packaged CSV inputs and dataset documentation
- `inst/scripts/`: reproducibility scripts for rebuilding the paper outputs

Repo-only items such as raw working files, generated results, and top-level
driver scripts are excluded from the package build through `.Rbuildignore`.

## Reproducibility workflow

If you are working from the repository source tree, run the scripts from the
package root:

```r
source("run_simulation_study.R")
source("run_section5_figures.R")
source("run_real_data_examples.R")
source("run_rebuild_all.R")
```

If the package is installed, equivalent scripts are available under
`system.file("scripts", package = "RSymHOMALS")`.

## Real-data analyses

The two empirical datasets are packaged in `inst/extdata/`:

- `clinical_comorbidity_5var_set_valued.csv`
- `yeast_gene_birdstyle_set_valued.csv`
- `dataset_documentation.docx`

To rerun the real-data analyses directly from R:

```r
library(RSymHOMALS)

clinical <- analyze_symbolic_csv(
  sym_example_path("clinical_comorbidity_5var_set_valued.csv"),
  include_mca = TRUE,
  max_iter = 70L,
  tol = 1e-5
)

yeast <- analyze_symbolic_csv(
  sym_example_path("yeast_gene_birdstyle_set_valued.csv"),
  include_mca = TRUE,
  max_iter = 70L,
  tol = 1e-5
)
```

The helper script `run_real_data_examples.R` writes the method-comparison
tables, score tables, and figure files used in Section 6.

## Session info

The package was most recently checked under the following R session:

```text
R version 4.5.0 (2025-04-11 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default
  LAPACK version 3.12.1

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

loaded via a namespace (and not attached):
[1] compiler_4.5.0
```

`RSymHOMALS` is dependency-light and currently relies only on base R and the
standard `stats`, `graphics`, `grDevices`, `utils`, and `methods` packages that
ship with R.
