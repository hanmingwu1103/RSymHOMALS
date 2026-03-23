source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in source_files) {
  source(file)
}

output_dir <- file.path("inst", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

clinical <- load_symbolic_csv(sym_example_path("clinical_comorbidity_5var_set_valued.csv"))
yeast <- load_symbolic_csv(sym_example_path("yeast_gene_birdstyle_set_valued.csv"))

clinical_sens <- analyze_weight_sensitivity(
  dataset = clinical,
  schemes = c("uniform", "inverse_frequency"),
  methods = c("euclidean", "wasserstein"),
  ndim = 2L,
  max_iter = 70L,
  tol = 1e-5,
  seed = 123
)

yeast_sens <- analyze_weight_sensitivity(
  dataset = yeast,
  schemes = c("uniform", "inverse_frequency"),
  methods = c("euclidean", "wasserstein"),
  ndim = 2L,
  max_iter = 70L,
  tol = 1e-5,
  seed = 123
)

write.csv(clinical_sens$comparison, file.path(output_dir, "clinical_weight_sensitivity_comparison.csv"), row.names = FALSE)
write.csv(clinical_sens$geometry, file.path(output_dir, "clinical_weight_sensitivity_geometry.csv"), row.names = FALSE)
write.csv(yeast_sens$comparison, file.path(output_dir, "yeast_weight_sensitivity_comparison.csv"), row.names = FALSE)
write.csv(yeast_sens$geometry, file.path(output_dir, "yeast_weight_sensitivity_geometry.csv"), row.names = FALSE)

print(clinical_sens$comparison)
print(clinical_sens$geometry)
print(yeast_sens$comparison)
print(yeast_sens$geometry)
