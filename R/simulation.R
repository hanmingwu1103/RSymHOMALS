simulate_symbolic_data <- function(
    n = NULL,
    m = NULL,
    k = NULL,
    ndim = 2L,
    scenario = c("geometry", "heterogeneous", "sparse"),
    seed = NULL) {
  scenario <- match.arg(scenario)
  cfg <- sym_scenario_config(scenario)
  if (is.null(n)) n <- cfg$n
  if (is.null(m)) m <- cfg$m
  if (is.null(k)) k <- cfg$k
  if (!is.null(seed)) {
    set.seed(seed)
  }
  k_vec <- if (length(k) == 1L) rep(as.integer(k), m) else as.integer(k)
  centers <- sym_cluster_centers(ndim)
  cluster_id <- sample(seq_len(nrow(centers)), n, replace = TRUE)
  true_x <- centers[cluster_id, , drop = FALSE] +
    matrix(stats::rnorm(n * ndim, sd = cfg$x_sd), nrow = n, ncol = ndim)
  true_x <- scale(true_x, center = TRUE, scale = FALSE)
  true_y <- vector("list", m)
  responses <- vector("list", m)
  set_sizes <- matrix(0L, nrow = n, ncol = m)

  for (j in seq_len(m)) {
    anchors <- centers[sample(seq_len(nrow(centers)), k_vec[j], replace = TRUE), , drop = FALSE]
    true_y[[j]] <- anchors +
      matrix(stats::rnorm(k_vec[j] * ndim, sd = cfg$y_sd), nrow = k_vec[j], ncol = ndim)
    responses[[j]] <- vector("list", n)
    for (i in seq_len(n)) {
      size <- sample(cfg$set_sizes, size = 1L, prob = cfg$probs)
      size <- min(size, k_vec[j])
      d2 <- rowSums((sweep(true_y[[j]], 2, true_x[i, ], "-"))^2)
      weights <- exp(-d2 / cfg$temp)
      weights <- weights * exp(stats::rnorm(k_vec[j], sd = cfg$weight_noise))
      if (scenario == "sparse") {
        weights <- weights * (1 + stats::runif(k_vec[j], min = 0, max = cfg$weight_noise))
      }
      idx <- sym_weighted_sample(weights, size = size)
      responses[[j]][[i]] <- idx
      set_sizes[i, j] <- length(idx)
    }
  }

  list(
    responses = responses,
    true_x = true_x,
    true_y = true_y,
    set_sizes = set_sizes,
    n_categories = k_vec,
    scenario = scenario
  )
}

sym_evaluate_simulation_fit <- function(fit, sim_data) {
  aligned <- sym_procrustes(fit$x, sim_data$true_x)
  size_total <- rowSums(sim_data$set_sizes)
  point_error <- rowSums((aligned$x_aligned - scale(sim_data$true_x, center = TRUE, scale = FALSE))^2)
  cardinality_bias <- suppressWarnings(abs(stats::cor(size_total, point_error)))
  if (!is.finite(cardinality_bias)) {
    cardinality_bias <- 0
  }
  data.frame(
    procrustes_rmse = aligned$rmse,
    neighbor_recovery = sym_neighbor_recovery(aligned$x_aligned, scale(sim_data$true_x, center = TRUE, scale = FALSE), k = 10L),
    cardinality_bias = cardinality_bias,
    euclidean_loss = sym_loss_euclidean(fit$x, fit$y, fit$Z_list),
    wasserstein_loss = sym_loss_wasserstein(fit$x, fit$y, fit$Z_list),
    hausdorff_loss = sym_loss_hausdorff(fit$x, fit$y, fit$B_list),
    runtime_sec = fit$runtime_sec,
    iterations = fit$iterations,
    stringsAsFactors = FALSE
  )
}

simulate_symhomals_study <- function(
    n_rep = 12L,
    scenarios = c("geometry", "heterogeneous", "sparse"),
    methods = c("euclidean", "wasserstein", "hausdorff", "mca"),
    ndim = 2L,
    max_iter = 80L,
    tol = 1e-6,
    seed = 20260322,
    verbose = TRUE) {
  raw_results <- vector("list", length(scenarios) * n_rep * length(methods))
  counter <- 1L
  for (scenario in scenarios) {
    cfg <- sym_scenario_config(scenario)
    if (verbose) {
      message("Running scenario: ", scenario)
    }
    for (replicate in seq_len(n_rep)) {
      sim_seed <- seed + 1000L * match(scenario, scenarios) + replicate
      sim_data <- simulate_symbolic_data(
        n = cfg$n,
        m = cfg$m,
        k = cfg$k,
        ndim = ndim,
        scenario = scenario,
        seed = sim_seed
      )
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
        metrics <- sym_evaluate_simulation_fit(fit, sim_data)
        raw_results[[counter]] <- cbind(
          data.frame(
            scenario = scenario,
            replicate = replicate,
            method = method,
            stringsAsFactors = FALSE
          ),
          metrics
        )
        counter <- counter + 1L
      }
    }
  }
  raw_results <- do.call(rbind, raw_results)
  list(
    raw = raw_results,
    summary = sym_summarize_results(raw_results)
  )
}
