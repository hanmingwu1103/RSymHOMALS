fit_naive_mca <- function(
    responses,
    n_categories = NULL,
    ndim = 2L,
    category_labels = NULL,
    variable_names = NULL,
    object_labels = NULL) {
  info <- sym_validate_responses(responses, n_categories = n_categories)
  B_list <- sym_binary_matrices(responses, info$n_categories)
  indicator <- do.call(cbind, B_list)
  timer <- proc.time()[3]

  total_mass <- sum(indicator)
  P <- indicator / total_mass
  r <- rowSums(P)
  c <- colSums(P)
  Dr_inv_half <- diag(1 / sqrt(r), nrow = length(r))
  Dc_inv_half <- diag(1 / sqrt(c), nrow = length(c))
  S <- Dr_inv_half %*% (P - outer(r, c)) %*% Dc_inv_half
  keep <- min(ndim, nrow(S), ncol(S))
  sv <- svd(S, nu = keep, nv = keep)

  x <- Dr_inv_half %*% sv$u[, seq_len(keep), drop = FALSE] %*%
    diag(sv$d[seq_len(keep)], nrow = keep)
  y_combined <- Dc_inv_half %*% sv$v[, seq_len(keep), drop = FALSE] %*%
    diag(sv$d[seq_len(keep)], nrow = keep)

  split_idx <- rep(seq_along(info$n_categories), info$n_categories)
  y <- lapply(seq_along(info$n_categories), function(j) {
    y_combined[split_idx == j, , drop = FALSE]
  })

  fit <- list(
    x = x,
    y = y,
    method = "mca",
    ndim = ncol(x),
    n_categories = info$n_categories,
    responses = responses,
    runtime_sec = unname(as.numeric(proc.time()[3] - timer)),
    inertia_2d = sum(sv$d[seq_len(keep)]^2),
    singular_values = sv$d[seq_len(keep)],
    category_labels = category_labels,
    variable_names = variable_names,
    object_labels = object_labels
  )
  class(fit) <- "naive_mca_fit"
  fit
}
