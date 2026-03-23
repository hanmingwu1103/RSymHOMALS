analyze_symbolic_dataset <- function(
    responses,
    n_categories = NULL,
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
