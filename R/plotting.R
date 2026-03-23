sym_plot_components <- function(fit) {
  if (inherits(fit, "symhomals_fit") || inherits(fit, "naive_mca_fit")) {
    return(list(x = fit$x[, 1:2, drop = FALSE], y = fit$y))
  }
  stop("Unsupported fit object.")
}

plot_symbolic_embedding <- function(
    fit,
    main = NULL,
    show_category_labels = TRUE,
    show_object_labels = FALSE,
    object_cex = 0.6,
    category_cex = 1.0,
    legend_cex = 0.8,
    ...) {
  comp <- sym_plot_components(fit)
  X <- comp$x
  Y_list <- comp$y
  variable_names <- fit$variable_names %||% paste0("Y", seq_along(Y_list))
  category_labels <- fit$category_labels %||% lapply(Y_list, function(Y) paste0("C", seq_len(nrow(Y))))
  object_labels <- fit$object_labels %||% seq_len(nrow(X))

  palette_cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#a6761d")
  pch_vals <- c(17, 15, 18, 0, 2, 5)
  colors <- rep(palette_cols, length.out = length(Y_list))
  pchs <- rep(pch_vals, length.out = length(Y_list))

  all_y <- do.call(rbind, Y_list)
  xlim <- range(c(X[, 1], all_y[, 1]))
  ylim <- range(c(X[, 2], all_y[, 2]))
  pad_x <- 0.08 * diff(xlim)
  pad_y <- 0.08 * diff(ylim)
  if (is.null(main)) {
    main <- if (inherits(fit, "naive_mca_fit")) "Naive MCA" else sprintf("%s Symbolic HOMALS", tools::toTitleCase(fit$method))
  }

  plot(
    X,
    col = grDevices::adjustcolor("grey35", alpha.f = 0.55),
    pch = 16,
    cex = object_cex,
    xlab = "Dim 1",
    ylab = "Dim 2",
    xlim = xlim + c(-pad_x, pad_x),
    ylim = ylim + c(-pad_y, pad_y),
    main = main,
    ...
  )
  if (show_object_labels) {
    graphics::text(X[, 1], X[, 2], labels = object_labels, cex = 0.55, pos = 3, col = "grey25")
  }
  for (j in seq_along(Y_list)) {
    Y <- Y_list[[j]][, 1:2, drop = FALSE]
    graphics::points(Y, col = colors[j], pch = pchs[j], cex = category_cex)
    if (show_category_labels) {
      graphics::text(
        Y[, 1], Y[, 2],
        labels = paste0("Y", j, ":", category_labels[[j]]),
        cex = 0.7,
        pos = 3,
        col = colors[j]
      )
    }
  }
  graphics::legend("topright", legend = variable_names, col = colors, pch = pchs, cex = legend_cex, bty = "n")
  invisible(NULL)
}

plot_symhomals_map <- function(
    fit,
    variable = 1L,
    main = NULL,
    show_labels = FALSE,
    object_col = "grey40",
    category_col = "firebrick",
    object_pch = 16,
    category_pch = 17,
    ...) {
  if (!inherits(fit, "symhomals_fit")) {
    stop("`fit` must be an object returned by `fit_symhomals()`.")
  }
  if (fit$ndim < 2L) {
    stop("`plot_symhomals_map()` requires at least two fitted dimensions.")
  }
  X <- fit$x[, 1:2, drop = FALSE]
  Y <- fit$y[[variable]][, 1:2, drop = FALSE]
  labels <- fit$category_labels[[variable]]
  rng_x <- range(c(X[, 1], Y[, 1]))
  rng_y <- range(c(X[, 2], Y[, 2]))
  if (is.null(main)) {
    main <- sprintf("Symbolic HOMALS Map (%s)", fit$method)
  }
  plot(X, col = object_col, pch = object_pch, xlab = "Dim 1", ylab = "Dim 2",
       xlim = rng_x, ylim = rng_y, main = main, ...)
  graphics::points(Y, col = category_col, pch = category_pch, cex = 1.2)
  if (show_labels) {
    graphics::text(Y[, 1], Y[, 2], labels = labels, pos = 3, col = category_col)
  }
  invisible(NULL)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
