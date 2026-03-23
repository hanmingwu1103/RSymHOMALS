source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in source_files) {
  source(file)
}

output_dir <- file.path("inst", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

clinical_data <- load_symbolic_csv(sym_example_path("clinical_comorbidity_5var_set_valued.csv"))
yeast_data <- load_symbolic_csv(sym_example_path("yeast_gene_birdstyle_set_valued.csv"))

clinical_interval <- sym_bootstrap_distance_interval(
  responses = clinical_data$responses,
  n_categories = clinical_data$n_categories,
  variable = 1L,
  label_a = "Congestive_heart_failure",
  label_b = "Renal_disease",
  category_labels = clinical_data$category_labels,
  variable_names = clinical_data$variable_names,
  method = "euclidean",
  ndim = 2L,
  B = 200L,
  max_iter = 70L,
  tol = 1e-5,
  seed = 20260324
)

yeast_interval <- sym_bootstrap_distance_interval(
  responses = yeast_data$responses,
  n_categories = yeast_data$n_categories,
  variable = 5L,
  label_a = "repair",
  label_b = "transport",
  category_labels = yeast_data$category_labels,
  variable_names = yeast_data$variable_names,
  method = "euclidean",
  ndim = 2L,
  B = 200L,
  max_iter = 70L,
  tol = 1e-5,
  seed = 20260325
)

write.csv(as.data.frame(clinical_interval), file.path(output_dir, "clinical_inference_example.csv"), row.names = FALSE)
write.csv(as.data.frame(yeast_interval), file.path(output_dir, "yeast_inference_example.csv"), row.names = FALSE)

n_grid <- c(50L, 100L, 200L, 500L)
n_rep <- 20L
boot_B <- 80L
coverage_rows <- vector("list", length(n_grid) * n_rep)
counter <- 1L

for (n_val in n_grid) {
  for (replicate in seq_len(n_rep)) {
    sim_seed <- 910000L + 1000L * match(n_val, n_grid) + replicate
    sim_data <- simulate_symbolic_data(
      n = n_val,
      m = 3L,
      k = 5L,
      ndim = 2L,
      scenario = "geometry",
      seed = sim_seed
    )
    interval <- sym_bootstrap_distance_interval(
      responses = sim_data$responses,
      n_categories = sim_data$n_categories,
      variable = 1L,
      label_a = 1L,
      label_b = 2L,
      category_labels = lapply(sim_data$n_categories, function(kj) paste0("C", seq_len(kj))),
      variable_names = paste0("Y", seq_along(sim_data$n_categories)),
      method = "euclidean",
      ndim = 2L,
      B = boot_B,
      max_iter = 70L,
      tol = 1e-5,
      seed = sim_seed + 500L
    )
    true_distance <- sqrt(sum((sim_data$true_y[[1]][1, ] - sim_data$true_y[[1]][2, ])^2))
    coverage_rows[[counter]] <- data.frame(
      n = n_val,
      replicate = replicate,
      true_distance = true_distance,
      estimate = interval$estimate,
      se_boot = interval$se_boot,
      wald_lower = interval$wald_lower,
      wald_upper = interval$wald_upper,
      covered = interval$wald_lower <= true_distance && true_distance <= interval$wald_upper,
      interval_width = interval$wald_upper - interval$wald_lower,
      n_boot = interval$n_boot,
      stringsAsFactors = FALSE
    )
    counter <- counter + 1L
  }
}

coverage_raw <- do.call(rbind, coverage_rows)
coverage_summary <- do.call(
  rbind,
  lapply(split(coverage_raw, coverage_raw$n), function(df) {
    data.frame(
      n = unique(df$n),
      n_rep = nrow(df),
      coverage = mean(df$covered),
      mean_width = mean(df$interval_width),
      mean_boot = mean(df$n_boot),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(coverage_raw, file.path(output_dir, "inference_coverage_raw.csv"), row.names = FALSE)
write.csv(coverage_summary, file.path(output_dir, "inference_coverage_summary.csv"), row.names = FALSE)

print(clinical_interval)
print(yeast_interval)
print(coverage_summary)
