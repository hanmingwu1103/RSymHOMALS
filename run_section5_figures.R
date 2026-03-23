source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in source_files) {
  source(file)
}

output_dir <- file.path("..", "LaTeX", "figures")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

symbolic_method_levels <- c("euclidean", "wasserstein", "hausdorff")
method_levels <- c("euclidean", "wasserstein", "hausdorff", "mca")
method_labels <- c(
  euclidean = "Euclidean",
  wasserstein = "Wasserstein",
  hausdorff = "Hausdorff",
  mca = "Naive MCA"
)
method_colors <- c(
  euclidean = "#1f4e79",
  wasserstein = "#c45508",
  hausdorff = "#6b1f1f",
  mca = "#4d4d4d"
)

align_fit_to_truth <- function(fit, sim_data, variable = 1L) {
  proc <- sym_procrustes(fit$x, sim_data$true_x)
  list(
    x = proc$x_aligned,
    y = proc$scale * fit$y[[variable]][, 1:2, drop = FALSE] %*% proc$rotation,
    true_x = scale(sim_data$true_x[, 1:2, drop = FALSE], center = TRUE, scale = FALSE),
    true_y = sim_data$true_y[[variable]][, 1:2, drop = FALSE]
  )
}

plot_map_panel <- function(obj_points, cat_points, main, object_col = "grey70",
                           category_col = "#c45508", show_labels = TRUE) {
  xlim <- range(c(obj_points[, 1], cat_points[, 1]))
  ylim <- range(c(obj_points[, 2], cat_points[, 2]))
  pad_x <- 0.08 * diff(xlim)
  pad_y <- 0.08 * diff(ylim)
  plot(obj_points, col = object_col, pch = 16, cex = 0.55,
       xlab = "Dim 1", ylab = "Dim 2",
       xlim = xlim + c(-pad_x, pad_x),
       ylim = ylim + c(-pad_y, pad_y),
       main = main)
  points(cat_points, col = category_col, pch = 17, cex = 1.2)
  if (show_labels) {
    text(cat_points[, 1], cat_points[, 2],
         labels = paste0("C", seq_len(nrow(cat_points))),
         col = category_col, pos = 3, cex = 0.8)
  }
}

draw_violin_boxplot <- function(values, labels, colors, ylab, main,
                                axis_ticks = NULL, axis_tick_labels = NULL) {
  at <- seq_along(values)
  y_all <- unlist(values, use.names = FALSE)
  y_range <- range(y_all, finite = TRUE)
  pad <- 0.08 * diff(y_range)
  if (!is.finite(pad) || pad <= 0) {
    pad <- max(0.1, 0.08 * max(abs(y_range), 1))
  }
  plot(
    NA,
    xlim = c(0.4, length(values) + 0.6),
    ylim = y_range + c(-pad, pad),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = ylab,
    main = main
  )
  if (is.null(axis_ticks)) {
    axis_ticks <- pretty(y_range)
    axis_tick_labels <- axis_ticks
  }
  abline(h = axis_ticks, col = "grey92", lwd = 0.8)

  for (i in seq_along(values)) {
    vals <- values[[i]]
    if (length(unique(vals)) > 1L) {
      dens <- density(vals, na.rm = TRUE)
      width <- 0.34 * dens$y / max(dens$y)
      polygon(
        x = c(at[i] - width, rev(at[i] + width)),
        y = c(dens$x, rev(dens$x)),
        col = grDevices::adjustcolor(colors[i], alpha.f = 0.32),
        border = colors[i],
        lwd = 1
      )
    } else {
      rect(
        xleft = at[i] - 0.10,
        xright = at[i] + 0.10,
        ybottom = vals[1] - 0.02,
        ytop = vals[1] + 0.02,
        col = grDevices::adjustcolor(colors[i], alpha.f = 0.32),
        border = colors[i]
      )
    }
  }

  boxplot(
    values,
    at = at,
    add = TRUE,
    xaxt = "n",
    yaxt = "n",
    boxwex = 0.16,
    outline = FALSE,
    col = grDevices::adjustcolor("white", alpha.f = 0.78),
    border = "grey20",
    whisklty = 1,
    staplelty = 1,
    medlwd = 1.2
  )

  points(
    at,
    vapply(values, stats::median, numeric(1), na.rm = TRUE),
    pch = 21,
    bg = colors,
    col = "grey15",
    cex = 0.8
  )

  axis(1, at = at, labels = labels)
  axis(2, at = axis_ticks, labels = axis_tick_labels, las = 1)
  box()
}

scenario_seed <- 20260401
geometry_example <- simulate_symbolic_data(
  n = 180L,
  m = 4L,
  k = 6L,
  ndim = 2L,
  scenario = "geometry",
  seed = scenario_seed
)

geometry_fits <- lapply(seq_along(symbolic_method_levels), function(i) {
  fit_symhomals(
    responses = geometry_example$responses,
    n_categories = geometry_example$n_categories,
    ndim = 2L,
    method = symbolic_method_levels[i],
    max_iter = 70L,
    tol = 1e-5,
    seed = scenario_seed + i,
    verbose = FALSE
  )
})
names(geometry_fits) <- symbolic_method_levels
aligned_fits <- lapply(geometry_fits, align_fit_to_truth, sim_data = geometry_example)

pdf(file.path(output_dir, "sim_geometry_maps.pdf"), width = 10.5, height = 8.2)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.5, 1.2))
plot_map_panel(
  obj_points = scale(geometry_example$true_x[, 1:2, drop = FALSE], center = TRUE, scale = FALSE),
  cat_points = geometry_example$true_y[[1]][, 1:2, drop = FALSE],
  main = "True latent map"
)
for (method in symbolic_method_levels) {
  aligned <- aligned_fits[[method]]
  plot_map_panel(
    obj_points = aligned$x,
    cat_points = aligned$y,
    main = paste(method_labels[[method]], "Symbolic HOMALS")
  )
}
mtext("Representative geometry-recovery example", outer = TRUE, line = -1.8, cex = 1.1)
dev.off()

raw_path <- file.path("inst", "results", "simulation_raw.csv")
if (!file.exists(raw_path)) {
  stop("simulation_raw.csv not found. Run run_simulation_study.R first.")
}
raw_results <- read.csv(raw_path, stringsAsFactors = FALSE)
raw_results$scenario <- factor(
  raw_results$scenario,
  levels = c("geometry", "heterogeneous", "sparse"),
  labels = c("Geometry", "Heterogeneous", "Sparse")
)
raw_results$method <- factor(
  raw_results$method,
  levels = method_levels,
  labels = unname(method_labels[method_levels])
)

pdf(file.path(output_dir, "sim_performance_boxplots.pdf"), width = 9.8, height = 9.6)
par(mfrow = c(3, 2), mar = c(4.2, 4.2, 2.2, 1.1))
for (scenario_name in levels(raw_results$scenario)) {
  dat <- raw_results[raw_results$scenario == scenario_name, , drop = FALSE]
  dat$runtime_plot <- pmax(dat$runtime_sec, 1e-3)
  rmse_values <- lapply(levels(dat$method), function(method_name) {
    dat$procrustes_rmse[dat$method == method_name]
  })
  draw_violin_boxplot(
    values = rmse_values,
    labels = levels(dat$method),
    colors = unname(method_colors[method_levels]),
    ylab = "Procrustes RMSE",
    main = paste(scenario_name, ": geometry recovery")
  )

  dat$runtime_log10 <- log10(dat$runtime_plot)
  runtime_values <- lapply(levels(dat$method), function(method_name) {
    dat$runtime_log10[dat$method == method_name]
  })
  runtime_ticks <- pretty(range(dat$runtime_log10))
  runtime_labels <- formatC(10^runtime_ticks, format = "fg", digits = 3)
  draw_violin_boxplot(
    values = runtime_values,
    labels = levels(dat$method),
    colors = unname(method_colors[method_levels]),
    ylab = "Runtime (sec)",
    main = paste(scenario_name, ": runtime"),
    axis_ticks = runtime_ticks,
    axis_tick_labels = runtime_labels
  )
}
mtext("Simulation performance across scenarios", outer = TRUE, line = -1.8, cex = 1.1)
dev.off()
