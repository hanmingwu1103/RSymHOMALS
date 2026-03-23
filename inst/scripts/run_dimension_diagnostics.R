source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in source_files) {
  source(file)
}

output_dir <- normalizePath(file.path("inst", "results"), winslash = "/", mustWork = FALSE)
latex_fig_dir <- normalizePath(file.path("..", "LaTeX", "figures"), winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(latex_fig_dir, recursive = TRUE, showWarnings = FALSE)

plot_dimension_curve <- function(diag_df, main_title, file_stub) {
  pdf(file.path(latex_fig_dir, paste0(file_stub, "_ple_curve.pdf")), width = 6.5, height = 4.5)
  par(mar = c(4.2, 4.2, 2.4, 1.2))
  plot(
    diag_df$p,
    diag_df$ple,
    type = "b",
    lwd = 2,
    pch = 19,
    col = "#2b8cbe",
    xlab = "Embedding dimension p",
    ylab = "PLE(p)",
    ylim = c(0, max(diag_df$ple, na.rm = TRUE) * 1.05),
    main = main_title
  )
  abline(v = 2, lty = 2, col = "#636363")
  dev.off()
}

clinical <- load_symbolic_csv(sym_example_path("clinical_comorbidity_5var_set_valued.csv"))
yeast <- load_symbolic_csv(sym_example_path("yeast_gene_birdstyle_set_valued.csv"))

clinical_diag <- analyze_dimension_diagnostics(
  responses = clinical$responses,
  n_categories = clinical$n_categories,
  p_max = 5L,
  method = "euclidean",
  max_iter = 70L,
  tol = 1e-5,
  seed = 123
)

yeast_diag <- analyze_dimension_diagnostics(
  responses = yeast$responses,
  n_categories = yeast$n_categories,
  p_max = 5L,
  method = "euclidean",
  max_iter = 70L,
  tol = 1e-5,
  seed = 123
)

write.csv(clinical_diag$diagnostics, file.path(output_dir, "clinical_dimension_diagnostics.csv"), row.names = FALSE)
write.csv(yeast_diag$diagnostics, file.path(output_dir, "yeast_dimension_diagnostics.csv"), row.names = FALSE)

plot_dimension_curve(clinical_diag$diagnostics, "Clinical data", "clinical")
plot_dimension_curve(yeast_diag$diagnostics, "Yeast data", "yeast")

print(clinical_diag$diagnostics)
print(yeast_diag$diagnostics)
