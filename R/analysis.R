analyze_symbolic_dataset <- function(
    responses,
    n_categories = NULL,
    response_weights = NULL,
    ndim = 2L,
    methods = c("euclidean", "wasserstein", "hausdorff"),
    max_iter = 80L,
    tol = 1e-6,
    seed = 123,
    verbose = FALSE,
    include_mca = FALSE,
    category_labels = NULL,
    variable_names = NULL,
    object_labels = NULL) {
  fits <- vector("list", length(methods))
  names(fits) <- methods
  comparison <- vector("list", length(methods) + if (include_mca) 1L else 0L)
  for (i in seq_along(methods)) {
    fit <- fit_symhomals(
      responses = responses,
      n_categories = n_categories,
      response_weights = response_weights,
      ndim = ndim,
      method = methods[i],
      max_iter = max_iter,
      tol = tol,
      seed = seed + i - 1L,
      verbose = verbose
    )
    fit$category_labels <- category_labels
    fit$variable_names <- variable_names
    fit$object_labels <- object_labels
    fits[[i]] <- fit
    comparison[[i]] <- data.frame(
      method = methods[i],
      objective = fit$objective,
      euclidean_loss = sym_loss_euclidean(fit$x, fit$y, fit$Z_list),
      wasserstein_loss = sym_loss_wasserstein(fit$x, fit$y, fit$Z_list),
      hausdorff_loss = sym_loss_hausdorff(fit$x, fit$y, fit$B_list),
      runtime_sec = fit$runtime_sec,
      iterations = fit$iterations,
      inertia_2d = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  if (include_mca) {
    mca_fit <- fit_naive_mca(
      responses = responses,
      n_categories = n_categories,
      ndim = ndim,
      category_labels = category_labels,
      variable_names = variable_names,
      object_labels = object_labels
    )
    fits[["mca"]] <- mca_fit
    comparison[[length(comparison)]] <- data.frame(
      method = "mca",
      objective = NA_real_,
      euclidean_loss = NA_real_,
      wasserstein_loss = NA_real_,
      hausdorff_loss = NA_real_,
      runtime_sec = mca_fit$runtime_sec,
      iterations = NA_real_,
      inertia_2d = mca_fit$inertia_2d,
      stringsAsFactors = FALSE
    )
  }
  list(
    fits = fits,
    comparison = do.call(rbind, comparison)
  )
}

analyze_symbolic_csv <- function(
    path,
    response_weights = NULL,
    ndim = 2L,
    methods = c("euclidean", "wasserstein", "hausdorff"),
    max_iter = 80L,
    tol = 1e-6,
    seed = 123,
    verbose = FALSE,
    include_mca = TRUE) {
  dataset <- load_symbolic_csv(path)
  analysis <- analyze_symbolic_dataset(
    responses = dataset$responses,
    n_categories = dataset$n_categories,
    response_weights = response_weights,
    ndim = ndim,
    methods = methods,
    max_iter = max_iter,
    tol = tol,
    seed = seed,
    verbose = verbose,
    include_mca = include_mca,
    category_labels = dataset$category_labels,
    variable_names = dataset$variable_names,
    object_labels = dataset$object_labels
  )
  analysis$dataset <- dataset
  analysis$summary <- summarize_symbolic_dataset(dataset)
  analysis
}

sym_compare_embeddings <- function(reference_fit, candidate_fit) {
  align <- sym_procrustes(candidate_fit$x, reference_fit$x)
  X_ref <- scale(reference_fit$x, center = TRUE, scale = FALSE)
  category_ref <- do.call(rbind, reference_fit$y)
  category_cand <- do.call(rbind, candidate_fit$y)
  category_ref <- scale(category_ref, center = TRUE, scale = FALSE)
  category_aligned <- align$scale * category_cand %*% align$rotation
  data.frame(
    object_rmse = sqrt(mean((align$x_aligned - X_ref)^2)),
    category_rmse = sqrt(mean((category_aligned - category_ref)^2)),
    stringsAsFactors = FALSE
  )
}

sym_barycentric_profile <- function(fit) {
  d_star <- if (!is.null(fit$d_star)) fit$d_star else sym_dstar_weights(fit$Z_list)
  X_num <- Reduce(`+`, Map(function(Z, Y) Z %*% Y, fit$Z_list, fit$y))
  X_bar <- sweep(X_num, 1, pmax(d_star, .Machine$double.eps), "/")
  X_bar_c <- sym_weighted_center(X_bar, d_star)
  Xw <- sweep(X_bar_c, 1, sqrt(d_star), "*")
  sv <- svd(Xw, nu = 0, nv = 0)$d
  if (!length(sv)) {
    return(data.frame(p = integer(0), singular_value = numeric(0), ple = numeric(0)))
  }
  inertia <- sv^2
  data.frame(
    p = seq_along(inertia),
    singular_value = sv,
    ple = cumsum(inertia) / sum(inertia),
    stringsAsFactors = FALSE
  )
}

analyze_dimension_diagnostics <- function(
    responses,
    n_categories = NULL,
    response_weights = NULL,
    p_max = 4L,
    method = "euclidean",
    max_iter = 80L,
    tol = 1e-6,
    seed = 123,
    verbose = FALSE) {
  fit <- fit_symhomals(
    responses = responses,
    n_categories = n_categories,
    response_weights = response_weights,
    ndim = max(1L, p_max),
    method = method,
    max_iter = max_iter,
    tol = tol,
    seed = seed,
    verbose = verbose
  )
  diag_rows <- sym_barycentric_profile(fit)
  diag_rows <- diag_rows[diag_rows$p <= p_max, , drop = FALSE]
  diag_rows$objective <- fit$objective
  diag_rows$runtime_sec <- fit$runtime_sec
  diag_rows$iterations <- fit$iterations
  diag_rows$converged <- fit$converged
  list(
    diagnostics = diag_rows,
    fit = fit,
    method = method
  )
}

analyze_weight_sensitivity <- function(
    dataset,
    schemes = c("uniform", "inverse_frequency"),
    methods = c("euclidean", "wasserstein"),
    ndim = 2L,
    max_iter = 80L,
    tol = 1e-6,
    seed = 123,
    verbose = FALSE,
    include_mca = FALSE) {
  if (is.character(dataset)) {
    dataset <- load_symbolic_csv(dataset)
  }
  analyses <- vector("list", length(schemes))
  names(analyses) <- schemes
  for (i in seq_along(schemes)) {
    analyses[[i]] <- analyze_symbolic_dataset(
      responses = dataset$responses,
      n_categories = dataset$n_categories,
      response_weights = sym_weight_scheme(dataset, schemes[i]),
      ndim = ndim,
      methods = methods,
      max_iter = max_iter,
      tol = tol,
      seed = seed + 10L * (i - 1L),
      verbose = verbose,
      include_mca = include_mca,
      category_labels = dataset$category_labels,
      variable_names = dataset$variable_names,
      object_labels = dataset$object_labels
    )
    analyses[[i]]$comparison$weight_scheme <- schemes[i]
  }
  if (!"uniform" %in% schemes) {
    stop("`schemes` must include 'uniform' so sensitivity can be referenced to the baseline fit.")
  }
  reference <- analyses[["uniform"]]$fits
  geometry_rows <- list()
  idx <- 1L
  for (scheme in schemes[schemes != "uniform"]) {
    for (method in intersect(names(reference), names(analyses[[scheme]]$fits))) {
      geom <- sym_compare_embeddings(reference[[method]], analyses[[scheme]]$fits[[method]])
      geometry_rows[[idx]] <- cbind(
        data.frame(
          method = method,
          reference_scheme = "uniform",
          comparison_scheme = scheme,
          stringsAsFactors = FALSE
        ),
        geom
      )
      idx <- idx + 1L
    }
  }
  list(
    analyses = analyses,
    comparison = do.call(rbind, lapply(analyses, `[[`, "comparison")),
    geometry = if (length(geometry_rows)) do.call(rbind, geometry_rows) else data.frame()
  )
}

analyze_scalability_benchmark <- function(
    designs,
    methods = c("euclidean", "wasserstein", "mca"),
    n_rep = 3L,
    scenario = "sparse",
    ndim = 2L,
    max_iter = 70L,
    tol = 1e-5,
    seed = 123,
    verbose = FALSE) {
  if (!is.data.frame(designs) || !all(c("label", "n", "m", "k") %in% names(designs))) {
    stop("`designs` must be a data frame containing columns `label`, `n`, `m`, and `k`.")
  }
  raw_rows <- vector("list", nrow(designs) * n_rep * length(methods))
  counter <- 1L
  for (d in seq_len(nrow(designs))) {
    if (verbose) {
      message("Benchmark design: ", designs$label[d])
    }
    for (replicate in seq_len(n_rep)) {
      sim_seed <- seed + 1000L * d + replicate
      sim_data <- simulate_symbolic_data(
        n = designs$n[d],
        m = designs$m[d],
        k = designs$k[d],
        ndim = ndim,
        scenario = scenario,
        seed = sim_seed
      )
      B_list <- sym_binary_matrices(sim_data$responses, sim_data$n_categories)
      Z_list <- sym_weighted_symbolic_matrices(sim_data$responses, sim_data$n_categories)
      nnz <- sum(vapply(Z_list, function(Z) sum(Z > 0), numeric(1)))
      total_categories <- sum(sim_data$n_categories)
      mean_set_size <- mean(unlist(lapply(sim_data$responses, lengths), use.names = FALSE))
      density <- nnz / (designs$n[d] * total_categories)
      input_size_mb <- as.numeric(utils::object.size(B_list) + utils::object.size(Z_list)) / 1024^2
      for (method in methods) {
        fit <- if (identical(method, "mca")) {
          fit_naive_mca(
            responses = sim_data$responses,
            n_categories = sim_data$n_categories,
            ndim = ndim
          )
        } else {
          fit_symhomals(
            responses = sim_data$responses,
            n_categories = sim_data$n_categories,
            ndim = ndim,
            method = method,
            max_iter = max_iter,
            tol = tol,
            seed = sim_seed,
            verbose = FALSE
          )
        }
        raw_rows[[counter]] <- data.frame(
          label = designs$label[d],
          n = designs$n[d],
          m = designs$m[d],
          k = designs$k[d],
          replicate = replicate,
          scenario = scenario,
          total_categories = total_categories,
          nnz = nnz,
          mean_set_size = mean_set_size,
          density = density,
          input_size_mb = input_size_mb,
          method = method,
          runtime_sec = fit$runtime_sec,
          iterations = if (is.null(fit$iterations)) NA_real_ else fit$iterations,
          converged = if (is.null(fit$converged)) TRUE else fit$converged,
          stringsAsFactors = FALSE
        )
        counter <- counter + 1L
      }
    }
  }
  raw_df <- do.call(rbind, raw_rows)
  design_stats <- stats::aggregate(
    cbind(total_categories, nnz, mean_set_size, density, input_size_mb) ~ label + n + m + k,
    raw_df,
    mean
  )
  runtime_mean <- stats::aggregate(runtime_sec ~ label + method, raw_df, mean)
  runtime_sd <- stats::aggregate(runtime_sec ~ label + method, raw_df, stats::sd)
  names(runtime_mean)[names(runtime_mean) == "runtime_sec"] <- "runtime_mean"
  names(runtime_sd)[names(runtime_sd) == "runtime_sec"] <- "runtime_sd"

  iter_mean <- stats::aggregate(iterations ~ label + method, raw_df, function(x) mean(x, na.rm = TRUE))
  iter_sd <- stats::aggregate(iterations ~ label + method, raw_df, function(x) stats::sd(x, na.rm = TRUE))
  conv_rate <- stats::aggregate(converged ~ label + method, raw_df, function(x) mean(as.numeric(x)))
  names(iter_mean)[names(iter_mean) == "iterations"] <- "iterations_mean"
  names(iter_sd)[names(iter_sd) == "iterations"] <- "iterations_sd"
  names(conv_rate)[names(conv_rate) == "converged"] <- "convergence_rate"

  summary_df <- Reduce(
    function(x, y) merge(x, y, by = intersect(names(x), names(y)), all = TRUE),
    list(runtime_mean, runtime_sd, iter_mean, iter_sd, conv_rate)
  )
  summary_df <- merge(summary_df, design_stats, by = "label", all.x = TRUE)
  summary_df <- summary_df[order(summary_df$nnz, summary_df$method), ]
  rownames(summary_df) <- NULL

  list(
    raw = raw_df,
    summary = summary_df
  )
}
