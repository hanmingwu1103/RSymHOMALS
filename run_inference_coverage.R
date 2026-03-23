source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in source_files) {
  source(file)
}

output_dir <- file.path("inst", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stack_category_matrix <- function(y_list, ndim) {
  do.call(rbind, lapply(y_list, function(Y) Y[, seq_len(ndim), drop = FALSE]))
}

bootstrap_contrast_bundle <- function(
    responses,
    n_categories,
    true_y,
    method = "euclidean",
    ndim = 2L,
    B = 80L,
    max_iter = 70L,
    tol = 1e-5,
    seed = 123) {
  target_matrix <- stack_category_matrix(true_y, ndim)
  target_centered <- scale(target_matrix, center = TRUE, scale = FALSE)

  ref_fit <- fit_symhomals(
    responses = responses,
    n_categories = n_categories,
    ndim = ndim,
    method = method,
    max_iter = max_iter,
    tol = tol,
    seed = seed,
    verbose = FALSE
  )
  ref_matrix <- stack_category_matrix(ref_fit$y, ndim)
  ref_aligned <- sym_procrustes(ref_matrix, target_matrix)$x_aligned

  point_distance <- sqrt(sum((ref_aligned[1L, ] - ref_aligned[2L, ])^2))
  true_distance <- sqrt(sum((target_centered[1L, ] - target_centered[2L, ])^2))
  point_coordinate <- ref_aligned[1L, 1L]
  true_coordinate <- target_centered[1L, 1L]

  n <- length(responses[[1L]])
  boot_distance <- rep(NA_real_, B)
  boot_coordinate <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    boot_responses <- sym_resample_responses(responses, idx)
    fit_b <- tryCatch(
      fit_symhomals(
        responses = boot_responses,
        n_categories = n_categories,
        ndim = ndim,
        method = method,
        max_iter = max_iter,
        tol = tol,
        seed = seed + b,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_b)) {
      boot_matrix <- stack_category_matrix(fit_b$y, ndim)
      boot_aligned <- sym_procrustes(boot_matrix, target_matrix)$x_aligned
      boot_distance[b] <- sqrt(sum((boot_aligned[1L, ] - boot_aligned[2L, ])^2))
      boot_coordinate[b] <- boot_aligned[1L, 1L]
    }
  }

  boot_distance <- boot_distance[is.finite(boot_distance)]
  boot_coordinate <- boot_coordinate[is.finite(boot_coordinate)]
  z_alpha <- stats::qnorm(0.975)

  list(
    distance = list(
      estimate = point_distance,
      true_value = true_distance,
      se_boot = stats::sd(boot_distance),
      lower = point_distance - z_alpha * stats::sd(boot_distance),
      upper = point_distance + z_alpha * stats::sd(boot_distance),
      n_boot = length(boot_distance)
    ),
    coordinate = list(
      estimate = point_coordinate,
      true_value = true_coordinate,
      se_boot = stats::sd(boot_coordinate),
      lower = point_coordinate - z_alpha * stats::sd(boot_coordinate),
      upper = point_coordinate + z_alpha * stats::sd(boot_coordinate),
      n_boot = length(boot_coordinate)
    )
  )
}

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
coverage_rows <- vector("list", length(n_grid) * n_rep * 2L)
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
    bundle <- bootstrap_contrast_bundle(
      responses = sim_data$responses,
      n_categories = sim_data$n_categories,
      true_y = sim_data$true_y,
      method = "euclidean",
      ndim = 2L,
      B = boot_B,
      max_iter = 70L,
      tol = 1e-5,
      seed = sim_seed + 500L
    )
    coverage_rows[[counter]] <- data.frame(
      n = n_val,
      replicate = replicate,
      contrast = "distance",
      true_value = bundle$distance$true_value,
      estimate = bundle$distance$estimate,
      se_boot = bundle$distance$se_boot,
      wald_lower = bundle$distance$lower,
      wald_upper = bundle$distance$upper,
      covered = bundle$distance$lower <= bundle$distance$true_value &&
        bundle$distance$true_value <= bundle$distance$upper,
      interval_width = bundle$distance$upper - bundle$distance$lower,
      n_boot = bundle$distance$n_boot,
      stringsAsFactors = FALSE
    )
    counter <- counter + 1L
    coverage_rows[[counter]] <- data.frame(
      n = n_val,
      replicate = replicate,
      contrast = "linear_coordinate",
      true_value = bundle$coordinate$true_value,
      estimate = bundle$coordinate$estimate,
      se_boot = bundle$coordinate$se_boot,
      wald_lower = bundle$coordinate$lower,
      wald_upper = bundle$coordinate$upper,
      covered = bundle$coordinate$lower <= bundle$coordinate$true_value &&
        bundle$coordinate$true_value <= bundle$coordinate$upper,
      interval_width = bundle$coordinate$upper - bundle$coordinate$lower,
      n_boot = bundle$coordinate$n_boot,
      stringsAsFactors = FALSE
    )
    counter <- counter + 1L
  }
}

coverage_raw <- do.call(rbind, coverage_rows)
coverage_summary <- do.call(
  rbind,
  lapply(split(coverage_raw, list(coverage_raw$contrast, coverage_raw$n), drop = TRUE), function(df) {
    data.frame(
      contrast = unique(df$contrast),
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
