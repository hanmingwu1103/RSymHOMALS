sym_ginv <- function(A, tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(A)) {
    A <- as.matrix(A)
  }
  if (!length(A)) {
    return(A)
  }
  sv <- svd(A)
  if (!length(sv$d)) {
    return(matrix(0, nrow = ncol(A), ncol = nrow(A)))
  }
  cutoff <- tol * max(sv$d)
  d_inv <- ifelse(sv$d > cutoff, 1 / sv$d, 0)
  sv$v %*% (diag(d_inv, nrow = length(d_inv)) %*% t(sv$u))
}

sym_weighted_sample <- function(weights, size) {
  weights <- as.numeric(weights)
  weights[!is.finite(weights)] <- 0
  weights <- pmax(weights, 0)
  n <- length(weights)
  size <- min(size, n)
  if (size <= 0) {
    return(integer(0))
  }
  if (sum(weights) <= 0) {
    return(sort(sample.int(n, size = size, replace = FALSE)))
  }
  sort(sample.int(n, size = size, replace = FALSE, prob = weights))
}

sym_validate_responses <- function(responses, n_categories = NULL) {
  if (!is.list(responses) || !length(responses)) {
    stop("`responses` must be a non-empty list of variables.")
  }
  n <- length(responses[[1L]])
  if (n == 0L) {
    stop("Each variable must contain at least one object.")
  }
  if (!all(vapply(responses, length, integer(1)) == n)) {
    stop("All variables must contain the same number of objects.")
  }
  if (is.null(n_categories)) {
    n_categories <- vapply(
      responses,
      function(variable) max(vapply(variable, max, integer(1))),
      integer(1)
    )
  }
  if (length(n_categories) == 1L) {
    n_categories <- rep(as.integer(n_categories), length(responses))
  }
  if (length(n_categories) != length(responses)) {
    stop("`n_categories` must have length 1 or length equal to `responses`.")
  }
  for (j in seq_along(responses)) {
    for (i in seq_len(n)) {
      idx <- responses[[j]][[i]]
      if (!length(idx)) {
        stop(sprintf("Object %d of variable %d has an empty response set.", i, j))
      }
      if (any(!is.finite(idx)) || any(idx < 1) || any(idx > n_categories[j])) {
        stop(sprintf("Invalid category index in variable %d, object %d.", j, i))
      }
    }
  }
  list(n = n, m = length(responses), n_categories = as.integer(n_categories))
}

sym_binary_matrices <- function(responses, n_categories) {
  n <- length(responses[[1L]])
  lapply(seq_along(responses), function(j) {
    B <- matrix(0, nrow = n, ncol = n_categories[j])
    for (i in seq_len(n)) {
      B[i, unique(as.integer(responses[[j]][[i]]))] <- 1
    }
    B
  })
}

sym_normalize_binary <- function(B_list) {
  lapply(B_list, function(B) {
    rs <- rowSums(B)
    sweep(B, 1, rs, "/")
  })
}

sym_total_cardinality <- function(B_list) {
  Reduce(`+`, lapply(B_list, rowSums))
}

sym_init_scores <- function(Z_list, ndim, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  Z <- do.call(cbind, Z_list)
  Zc <- scale(Z, center = TRUE, scale = FALSE)
  sv <- svd(Zc, nu = min(ndim, nrow(Zc), ncol(Zc)), nv = 0)
  if (length(sv$d) == 0L) {
    X <- matrix(stats::rnorm(nrow(Z) * ndim), nrow(Z), ndim)
  } else {
    X <- sqrt(nrow(Z)) * sv$u[, seq_len(min(ndim, ncol(sv$u))), drop = FALSE]
  }
  sym_center_orthonormalize(X, ndim = ndim)
}

sym_center_orthonormalize <- function(X, ndim = ncol(X)) {
  X <- as.matrix(X)
  n <- nrow(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  sv <- svd(Xc, nu = min(ndim, nrow(Xc), ncol(Xc)), nv = 0)
  U <- if (length(sv$d)) sv$u else matrix(0, nrow = n, ncol = 0)
  if (ncol(U) < ndim) {
    extra <- matrix(stats::rnorm(n * ndim), nrow = n, ncol = ndim)
    basis <- qr.Q(qr(cbind(U, extra)))
    U <- basis[, seq_len(ndim), drop = FALSE]
  } else {
    U <- U[, seq_len(ndim), drop = FALSE]
  }
  Xnew <- sqrt(n) * U
  colnames(Xnew) <- paste0("Dim", seq_len(ndim))
  Xnew
}

sym_update_means_from_matrix <- function(A, X, Y_prev = NULL) {
  weights <- colSums(A)
  Y <- crossprod(A, X)
  for (r in seq_len(nrow(Y))) {
    if (weights[r] > 0) {
      Y[r, ] <- Y[r, ] / weights[r]
    } else if (!is.null(Y_prev)) {
      Y[r, ] <- Y_prev[r, ]
    } else {
      Y[r, ] <- 0
    }
  }
  Y
}

sym_active_hausdorff <- function(B_list, X, Y_list, alpha = 8) {
  lapply(seq_along(B_list), function(j) {
    B <- B_list[[j]]
    Y <- Y_list[[j]]
    A <- matrix(0, nrow = nrow(B), ncol = ncol(B))
    for (i in seq_len(nrow(B))) {
      idx <- which(B[i, ] > 0)
      diffs <- sweep(Y[idx, , drop = FALSE], 2, X[i, ], "-")
      d2 <- rowSums(diffs * diffs)
      logits <- alpha * (d2 - max(d2))
      w <- exp(logits)
      w <- w / sum(w)
      A[i, idx] <- w
    }
    A
  })
}

sym_loss_euclidean <- function(X, Y_list, Z_list) {
  n <- nrow(X)
  m <- length(Z_list)
  total <- 0
  for (j in seq_len(m)) {
    resid <- X - Z_list[[j]] %*% Y_list[[j]]
    total <- total + sum(rowSums(resid * resid))
  }
  total / (n * m)
}

sym_loss_wasserstein <- function(X, Y_list, Z_list) {
  n <- nrow(X)
  m <- length(Z_list)
  total <- 0
  for (j in seq_len(m)) {
    Z <- Z_list[[j]]
    Y <- Y_list[[j]]
    for (i in seq_len(n)) {
      idx <- which(Z[i, ] > 0)
      diffs <- sweep(Y[idx, , drop = FALSE], 2, X[i, ], "-")
      total <- total + sum(Z[i, idx] * rowSums(diffs * diffs))
    }
  }
  total / (n * m)
}

sym_loss_hausdorff <- function(X, Y_list, B_list) {
  n <- nrow(X)
  m <- length(B_list)
  total <- 0
  for (j in seq_len(m)) {
    B <- B_list[[j]]
    Y <- Y_list[[j]]
    for (i in seq_len(n)) {
      idx <- which(B[i, ] > 0)
      diffs <- sweep(Y[idx, , drop = FALSE], 2, X[i, ], "-")
      total <- total + max(rowSums(diffs * diffs))
    }
  }
  total / (n * m)
}

sym_procrustes <- function(X_hat, X_true) {
  X_hat <- as.matrix(X_hat)
  X_true <- as.matrix(X_true)
  X_hat_c <- scale(X_hat, center = TRUE, scale = FALSE)
  X_true_c <- scale(X_true, center = TRUE, scale = FALSE)
  M <- crossprod(X_hat_c, X_true_c)
  sv <- svd(M)
  R <- sv$u %*% t(sv$v)
  scale_factor <- sum(sv$d) / sum(X_hat_c * X_hat_c)
  X_aligned <- scale_factor * X_hat_c %*% R
  list(
    x_aligned = X_aligned,
    rotation = R,
    scale = scale_factor,
    rmse = sqrt(mean((X_aligned - X_true_c)^2))
  )
}

sym_neighbor_recovery <- function(X_hat, X_true, k = 10L) {
  n <- nrow(X_hat)
  k <- min(k, n - 1L)
  d_hat <- as.matrix(stats::dist(X_hat))
  d_true <- as.matrix(stats::dist(X_true))
  overlap <- numeric(n)
  for (i in seq_len(n)) {
    nn_hat <- order(d_hat[i, ])[2:(k + 1L)]
    nn_true <- order(d_true[i, ])[2:(k + 1L)]
    overlap[i] <- length(intersect(nn_hat, nn_true)) / k
  }
  mean(overlap)
}

sym_cluster_centers <- function(ndim) {
  if (ndim == 1L) {
    return(matrix(c(-1.5, 0, 1.5), ncol = 1))
  }
  base <- matrix(
    c(-1.8, -0.5,
      0.0,  1.6,
      1.8, -0.5),
    ncol = 2,
    byrow = TRUE
  )
  if (ndim == 2L) {
    return(base)
  }
  extra <- matrix(0, nrow = nrow(base), ncol = ndim - 2L)
  cbind(base, extra)
}

sym_scenario_config <- function(name) {
  switch(
    name,
    geometry = list(
      n = 180L, m = 4L, k = 6L, temp = 0.35,
      set_sizes = 1:3, probs = c(0.60, 0.25, 0.15),
      x_sd = 0.30, y_sd = 0.55, weight_noise = 0.08
    ),
    heterogeneous = list(
      n = 180L, m = 4L, k = 8L, temp = 0.65,
      set_sizes = 1:5, probs = c(0.10, 0.15, 0.20, 0.25, 0.30),
      x_sd = 0.35, y_sd = 0.75, weight_noise = 0.25
    ),
    sparse = list(
      n = 260L, m = 8L, k = 16L, temp = 0.90,
      set_sizes = 1:4, probs = c(0.25, 0.35, 0.25, 0.15),
      x_sd = 0.35, y_sd = 0.95, weight_noise = 0.35
    ),
    stop("Unknown scenario: ", name)
  )
}

sym_summarize_results <- function(raw_results) {
  metrics <- setdiff(colnames(raw_results), c("scenario", "replicate", "method"))
  out <- NULL
  for (metric in metrics) {
    mean_df <- stats::aggregate(raw_results[[metric]], raw_results[c("scenario", "method")], mean)
    sd_df <- stats::aggregate(raw_results[[metric]], raw_results[c("scenario", "method")], stats::sd)
    names(mean_df)[3] <- paste0(metric, "_mean")
    names(sd_df)[3] <- paste0(metric, "_sd")
    merged <- merge(mean_df, sd_df, by = c("scenario", "method"))
    out <- if (is.null(out)) merged else merge(out, merged, by = c("scenario", "method"))
  }
  out[order(out$scenario, out$method), ]
}
