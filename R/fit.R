sym_fit_euclidean_core <- function(Z_list, ndim, max_iter, tol, seed, verbose) {
  m <- length(Z_list)
  X <- sym_init_scores(Z_list, ndim = ndim, seed = seed)
  loss_history <- numeric(max_iter)
  converged <- FALSE
  Y_list <- vector("list", m)
  for (iter in seq_len(max_iter)) {
    Y_list <- lapply(Z_list, function(Z) sym_ginv(crossprod(Z)) %*% crossprod(Z, X))
    X_tilde <- Reduce(`+`, Map(function(Z, Y) Z %*% Y, Z_list, Y_list)) / m
    X_new <- sym_center_orthonormalize(X_tilde, ndim = ndim)
    loss_history[iter] <- sym_loss_euclidean(X_new, Y_list, Z_list)
    if (verbose) {
      message(sprintf("[euclidean] iter=%d loss=%.6f", iter, loss_history[iter]))
    }
    if (iter > 1L) {
      rel_change <- abs(loss_history[iter - 1L] - loss_history[iter]) / max(1, loss_history[iter - 1L])
      if (rel_change < tol) {
        converged <- TRUE
        X <- X_new
        break
      }
    }
    X <- X_new
  }
  iter_final <- if (converged) iter else max_iter
  Y_list <- lapply(Z_list, function(Z) sym_ginv(crossprod(Z)) %*% crossprod(Z, X))
  list(
    x = X,
    y = Y_list,
    iterations = iter_final,
    converged = converged,
    loss_history = loss_history[seq_len(iter_final)],
    objective = sym_loss_euclidean(X, Y_list, Z_list)
  )
}

sym_fit_wasserstein_core <- function(Z_list, ndim, max_iter, tol, seed, verbose) {
  m <- length(Z_list)
  X <- sym_init_scores(Z_list, ndim = ndim, seed = seed)
  loss_history <- numeric(max_iter)
  converged <- FALSE
  Y_list <- vector("list", m)
  for (iter in seq_len(max_iter)) {
    Y_list <- lapply(Z_list, function(Z) {
      weights <- colSums(Z)
      Y <- crossprod(Z, X)
      for (r in seq_len(nrow(Y))) {
        if (weights[r] > 0) {
          Y[r, ] <- Y[r, ] / weights[r]
        } else {
          Y[r, ] <- 0
        }
      }
      Y
    })
    X_tilde <- Reduce(`+`, Map(function(Z, Y) Z %*% Y, Z_list, Y_list)) / m
    X_new <- sym_center_orthonormalize(X_tilde, ndim = ndim)
    loss_history[iter] <- sym_loss_wasserstein(X_new, Y_list, Z_list)
    if (verbose) {
      message(sprintf("[wasserstein] iter=%d loss=%.6f", iter, loss_history[iter]))
    }
    if (iter > 1L) {
      rel_change <- abs(loss_history[iter - 1L] - loss_history[iter]) / max(1, loss_history[iter - 1L])
      if (rel_change < tol) {
        converged <- TRUE
        X <- X_new
        break
      }
    }
    X <- X_new
  }
  iter_final <- if (converged) iter else max_iter
  Y_list <- lapply(Z_list, function(Z) {
    weights <- colSums(Z)
    Y <- crossprod(Z, X)
    for (r in seq_len(nrow(Y))) {
      if (weights[r] > 0) {
        Y[r, ] <- Y[r, ] / weights[r]
      }
    }
    Y
  })
  list(
    x = X,
    y = Y_list,
    iterations = iter_final,
    converged = converged,
    loss_history = loss_history[seq_len(iter_final)],
    objective = sym_loss_wasserstein(X, Y_list, Z_list)
  )
}

sym_fit_hausdorff_core <- function(Z_list, B_list, ndim, max_iter, tol, seed, verbose) {
  init_fit <- sym_fit_euclidean_core(
    Z_list = Z_list,
    ndim = ndim,
    max_iter = min(20L, max_iter),
    tol = 1e-4,
    seed = seed,
    verbose = FALSE
  )
  X <- init_fit$x
  Y_prev <- init_fit$y
  m <- length(B_list)
  loss_history <- numeric(max_iter)
  converged <- FALSE
  Y_list <- Y_prev
  for (iter in seq_len(max_iter)) {
    A_list <- sym_active_hausdorff(B_list, X, Y_prev)
    Y_list <- Map(function(A, Y_old) sym_update_means_from_matrix(A, X, Y_old), A_list, Y_prev)
    X_tilde <- Reduce(`+`, Map(function(A, Y) A %*% Y, A_list, Y_list)) / m
    X_new <- sym_center_orthonormalize(X_tilde, ndim = ndim)
    loss_history[iter] <- sym_loss_hausdorff(X_new, Y_list, B_list)
    if (verbose) {
      message(sprintf("[hausdorff] iter=%d loss=%.6f", iter, loss_history[iter]))
    }
    if (iter > 1L) {
      rel_change <- abs(loss_history[iter - 1L] - loss_history[iter]) / max(1, loss_history[iter - 1L])
      if (rel_change < tol) {
        converged <- TRUE
        X <- X_new
        break
      }
    }
    X <- X_new
    Y_prev <- Y_list
  }
  iter_final <- if (converged) iter else max_iter
  list(
    x = X,
    y = Y_list,
    iterations = iter_final,
    converged = converged,
    loss_history = loss_history[seq_len(iter_final)],
    objective = sym_loss_hausdorff(X, Y_list, B_list)
  )
}

fit_symhomals <- function(
    responses,
    n_categories = NULL,
    ndim = 2L,
    method = c("euclidean", "wasserstein", "hausdorff"),
    max_iter = 80L,
    tol = 1e-6,
    seed = 123,
    verbose = FALSE) {
  method <- match.arg(method)
  info <- sym_validate_responses(responses, n_categories = n_categories)
  B_list <- sym_binary_matrices(responses, info$n_categories)
  Z_list <- sym_normalize_binary(B_list)
  timer <- proc.time()[3]
  fit <- switch(
    method,
    euclidean = sym_fit_euclidean_core(Z_list, ndim, max_iter, tol, seed, verbose),
    wasserstein = sym_fit_wasserstein_core(Z_list, ndim, max_iter, tol, seed, verbose),
    hausdorff = sym_fit_hausdorff_core(Z_list, B_list, ndim, max_iter, tol, seed, verbose)
  )
  fit$runtime_sec <- unname(as.numeric(proc.time()[3] - timer))
  fit$method <- method
  fit$ndim <- ndim
  fit$n_categories <- info$n_categories
  fit$B_list <- B_list
  fit$Z_list <- Z_list
  fit$responses <- responses
  class(fit) <- "symhomals_fit"
  fit
}

print.symhomals_fit <- function(x, ...) {
  cat("RSymHOMALS fit\n")
  cat("  method:    ", x$method, "\n", sep = "")
  cat("  dimensions:", x$ndim, "\n")
  cat("  iterations:", x$iterations, "\n")
  cat("  converged: ", x$converged, "\n", sep = "")
  cat("  objective: ", sprintf("%.6f", x$objective), "\n", sep = "")
  cat("  runtime:   ", sprintf("%.3f sec", x$runtime_sec), "\n", sep = "")
  invisible(x)
}
