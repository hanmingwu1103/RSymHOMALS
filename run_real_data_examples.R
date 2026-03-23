source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in source_files) {
  source(file)
}

output_dir <- file.path("inst", "results")
latex_fig_dir <- file.path("..", "LaTeX", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(latex_fig_dir, recursive = TRUE, showWarnings = FALSE)

export_method_panel <- function(analysis, base_name, width = 12, height = 10) {
  fit_order <- c("euclidean", "wasserstein", "hausdorff", "mca")
  panel_fits <- analysis$fits[fit_order]
  panel_fits <- panel_fits[!vapply(panel_fits, is.null, logical(1))]

  pdf(file.path(latex_fig_dir, paste0(base_name, "_method_maps.pdf")), width = width, height = height)
  par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.6, 1.2))
  for (nm in names(panel_fits)) {
    plot_symbolic_embedding(panel_fits[[nm]], main = switch(
      nm,
      euclidean = "Euclidean Symbolic HOMALS",
      wasserstein = "Wasserstein Symbolic HOMALS",
      hausdorff = "Hausdorff Symbolic HOMALS",
      mca = "Naive MCA benchmark"
    ))
  }
  dev.off()

  png(file.path(latex_fig_dir, paste0(base_name, "_method_maps.png")), width = 1800, height = 1500, res = 200)
  par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.6, 1.2))
  for (nm in names(panel_fits)) {
    plot_symbolic_embedding(panel_fits[[nm]], main = switch(
      nm,
      euclidean = "Euclidean Symbolic HOMALS",
      wasserstein = "Wasserstein Symbolic HOMALS",
      hausdorff = "Hausdorff Symbolic HOMALS",
      mca = "Naive MCA benchmark"
    ))
  }
  dev.off()
}

export_coordinate_tables <- function(analysis, base_name) {
  for (nm in names(analysis$fits)) {
    fit <- analysis$fits[[nm]]
    obj <- data.frame(
      object = fit$object_labels,
      dim1 = fit$x[, 1],
      dim2 = fit$x[, 2],
      stringsAsFactors = FALSE
    )
    write.csv(obj, file.path(output_dir, paste0(base_name, "_", nm, "_object_scores.csv")), row.names = FALSE)

    cat_df <- do.call(
      rbind,
      lapply(seq_along(fit$y), function(j) {
        data.frame(
          variable = fit$variable_names[j],
          label = fit$category_labels[[j]],
          dim1 = fit$y[[j]][, 1],
          dim2 = fit$y[[j]][, 2],
          stringsAsFactors = FALSE
        )
      })
    )
    write.csv(cat_df, file.path(output_dir, paste0(base_name, "_", nm, "_category_scores.csv")), row.names = FALSE)
  }
}

clinical_analysis <- analyze_symbolic_csv(
  path = sym_example_path("clinical_comorbidity_5var_set_valued.csv"),
  ndim = 2L,
  max_iter = 70L,
  tol = 1e-5,
  seed = 123,
  include_mca = TRUE
)

yeast_analysis <- analyze_symbolic_csv(
  path = sym_example_path("yeast_gene_birdstyle_set_valued.csv"),
  ndim = 2L,
  max_iter = 70L,
  tol = 1e-5,
  seed = 123,
  include_mca = TRUE
)

write.csv(clinical_analysis$comparison, file.path(output_dir, "clinical_method_comparison.csv"), row.names = FALSE)
write.csv(yeast_analysis$comparison, file.path(output_dir, "yeast_method_comparison.csv"), row.names = FALSE)

write.csv(clinical_analysis$summary$variable_summary, file.path(output_dir, "clinical_variable_summary.csv"), row.names = FALSE)
write.csv(yeast_analysis$summary$variable_summary, file.path(output_dir, "yeast_variable_summary.csv"), row.names = FALSE)

write.csv(clinical_analysis$summary$label_frequency, file.path(output_dir, "clinical_label_frequency.csv"), row.names = FALSE)
write.csv(yeast_analysis$summary$label_frequency, file.path(output_dir, "yeast_label_frequency.csv"), row.names = FALSE)

write.csv(clinical_analysis$summary$richest_objects, file.path(output_dir, "clinical_richest_objects.csv"), row.names = FALSE)
write.csv(yeast_analysis$summary$richest_objects, file.path(output_dir, "yeast_richest_objects.csv"), row.names = FALSE)

saveRDS(clinical_analysis, file.path(output_dir, "clinical_analysis.rds"))
saveRDS(yeast_analysis, file.path(output_dir, "yeast_analysis.rds"))

export_method_panel(clinical_analysis, "clinical")
export_method_panel(yeast_analysis, "yeast")

export_coordinate_tables(clinical_analysis, "clinical")
export_coordinate_tables(yeast_analysis, "yeast")

print(clinical_analysis$comparison)
print(yeast_analysis$comparison)
print(clinical_analysis$summary$variable_summary)
print(yeast_analysis$summary$variable_summary)
