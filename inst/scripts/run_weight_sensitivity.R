source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in source_files) {
  source(file)
}

output_dir <- file.path("inst", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

category_shift_diagnostics <- function(dataset, sensitivity_result, method) {
  ref_fit <- sensitivity_result$analyses[["uniform"]]$fits[[method]]
  cmp_fit <- sensitivity_result$analyses[["inverse_frequency"]]$fits[[method]]
  align <- sym_procrustes(cmp_fit$x, ref_fit$x)
  x_ref <- scale(ref_fit$x, center = TRUE, scale = FALSE)
  object_corr <- stats::cor(c(x_ref), c(align$x_aligned))

  category_ref <- do.call(rbind, ref_fit$y)
  category_cmp <- do.call(rbind, cmp_fit$y)
  category_ref_c <- scale(category_ref, center = TRUE, scale = FALSE)
  category_cmp_aligned <- align$scale * category_cmp %*% align$rotation
  shift_norm <- sqrt(rowSums((category_cmp_aligned - category_ref_c)^2))

  B_list <- sym_binary_matrices(dataset$responses, dataset$n_categories)
  frequencies <- unlist(lapply(B_list, colSums), use.names = FALSE)
  variable_index <- rep(seq_along(dataset$n_categories), times = dataset$n_categories)
  variable_name <- rep(dataset$variable_names, times = dataset$n_categories)
  label <- unlist(dataset$category_labels, use.names = FALSE)

  data.frame(
    method = method,
    variable = variable_index,
    variable_name = variable_name,
    label = label,
    frequency = frequencies,
    shift_norm = shift_norm,
    object_procrustes_correlation = object_corr,
    stringsAsFactors = FALSE
  )[order(-shift_norm), ]
}

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

clinical_shift_details <- do.call(
  rbind,
  lapply(c("euclidean", "wasserstein"), function(method) {
    category_shift_diagnostics(clinical, clinical_sens, method)
  })
)
yeast_shift_details <- do.call(
  rbind,
  lapply(c("euclidean", "wasserstein"), function(method) {
    category_shift_diagnostics(yeast, yeast_sens, method)
  })
)

write.csv(clinical_shift_details, file.path(output_dir, "clinical_weight_sensitivity_category_shifts.csv"), row.names = FALSE)
write.csv(yeast_shift_details, file.path(output_dir, "yeast_weight_sensitivity_category_shifts.csv"), row.names = FALSE)

print(clinical_sens$comparison)
print(clinical_sens$geometry)
print(yeast_sens$comparison)
print(yeast_sens$geometry)
print(head(clinical_shift_details, 10L))
print(head(yeast_shift_details, 10L))
