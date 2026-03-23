sym_category_distance <- function(fit, variable, label_a, label_b) {
  if (is.character(variable)) {
    variable <- match(variable, fit$variable_names)
  }
  if (is.na(variable) || variable < 1L || variable > length(fit$y)) {
    stop("`variable` does not identify a valid variable.")
  }
  labels <- fit$category_labels[[variable]]
  idx_a <- if (is.character(label_a)) match(label_a, labels) else as.integer(label_a)
  idx_b <- if (is.character(label_b)) match(label_b, labels) else as.integer(label_b)
  if (anyNA(c(idx_a, idx_b))) {
    stop("Both category labels must be present in the selected variable.")
  }
  ya <- fit$y[[variable]][idx_a, 1:fit$ndim, drop = TRUE]
  yb <- fit$y[[variable]][idx_b, 1:fit$ndim, drop = TRUE]
  sqrt(sum((ya - yb)^2))
}

sym_resample_responses <- function(responses, indices) {
  lapply(responses, function(variable) variable[indices])
}

sym_bootstrap_distance_interval <- function(
    responses,
    n_categories,
    variable,
    label_a,
    label_b,
    category_labels,
    variable_names = NULL,
    method = "euclidean",
    ndim = 2L,
    B = 200L,
    max_iter = 70L,
    tol = 1e-5,
    seed = 123,
    conf_level = 0.95,
    verbose = FALSE) {
  set.seed(seed)
  reference_fit <- fit_symhomals(
    responses = responses,
    n_categories = n_categories,
    ndim = ndim,
    method = method,
    max_iter = max_iter,
    tol = tol,
    seed = seed,
    verbose = verbose
  )
  reference_fit$category_labels <- category_labels
  reference_fit$variable_names <- variable_names
  point_estimate <- sym_category_distance(reference_fit, variable, label_a, label_b)

  n <- length(responses[[1L]])
  boot_values <- rep(NA_real_, B)
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
      fit_b$category_labels <- category_labels
      fit_b$variable_names <- variable_names
      boot_values[b] <- sym_category_distance(fit_b, variable, label_a, label_b)
    }
  }

  boot_values <- boot_values[is.finite(boot_values)]
  if (!length(boot_values)) {
    stop("Bootstrap interval failed because all bootstrap refits returned non-finite values.")
  }

  alpha <- 1 - conf_level
  z_alpha <- stats::qnorm(1 - alpha / 2)
  se_boot <- stats::sd(boot_values)
  normal_ci <- point_estimate + c(-1, 1) * z_alpha * se_boot
  percentile_ci <- stats::quantile(boot_values, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)

  list(
    estimate = point_estimate,
    se_boot = se_boot,
    wald_lower = normal_ci[1],
    wald_upper = normal_ci[2],
    percentile_lower = percentile_ci[1],
    percentile_upper = percentile_ci[2],
    n_boot = length(boot_values),
    method = method,
    variable = variable,
    label_a = if (is.character(label_a)) label_a else category_labels[[variable]][label_a],
    label_b = if (is.character(label_b)) label_b else category_labels[[variable]][label_b]
  )
}
