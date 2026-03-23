library(RSymHOMALS)

result_dir <- file.path(getwd(), "results")
figure_dir <- file.path(getwd(), "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

method_levels <- c("euclidean", "wasserstein", "hausdorff")
method_labels <- c(
  euclidean = "Euclidean",
  wasserstein = "Wasserstein",
  hausdorff = "Hausdorff"
)
method_colors <- c(
  euclidean = "#1f4e79",
  wasserstein = "#c45508",
  hausdorff = "#6b1f1f"
)

align_fit_to_truth <- function(fit, sim_data, variable = 1L) {
  proc <- RSymHOMALS:::sym_procrustes(fit$x, sim_data$true_x)
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
  graphics::plot(obj_points, col = object_col, pch = 16, cex = 0.55,
                 xlab = "Dim 1", ylab = "Dim 2",
                 xlim = xlim + c(-pad_x, pad_x),
                 ylim = ylim + c(-pad_y, pad_y),
                 main = main)
  graphics::points(cat_points, col = category_col, pch = 17, cex = 1.2)
  if (show_labels) {
    graphics::text(cat_points[, 1], cat_points[, 2],
                   labels = paste0("C", seq_len(nrow(cat_points))),
                   col = category_col, pos = 3, cex = 0.8)
  }
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

geometry_fits <- lapply(seq_along(method_levels), function(i) {
  fit_symhomals(
    responses = geometry_example$responses,
    n_categories = geometry_example$n_categories,
    ndim = 2L,
    method = method_levels[i],
    max_iter = 70L,
    tol = 1e-5,
    seed = scenario_seed + i,
    verbose = FALSE
  )
})
names(geometry_fits) <- method_levels
aligned_fits <- lapply(geometry_fits, align_fit_to_truth, sim_data = geometry_example)

grDevices::pdf(file.path(figure_dir, "sim_geometry_maps.pdf"), width = 10.5, height = 8.2)
graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.5, 1.2))
plot_map_panel(
  obj_points = scale(geometry_example$true_x[, 1:2, drop = FALSE], center = TRUE, scale = FALSE),
  cat_points = geometry_example$true_y[[1]][, 1:2, drop = FALSE],
  main = "True latent map"
)
for (method in method_levels) {
  aligned <- aligned_fits[[method]]
  plot_map_panel(
    obj_points = aligned$x,
    cat_points = aligned$y,
    main = paste(method_labels[[method]], "Symbolic HOMALS")
  )
}
graphics::mtext("Representative geometry-recovery example", outer = TRUE, line = -1.8, cex = 1.1)
grDevices::dev.off()

raw_path <- file.path(result_dir, "simulation_raw.csv")
if (!file.exists(raw_path)) {
  stop("simulation_raw.csv not found. Run run_simulation_study.R first.")
}
raw_results <- utils::read.csv(raw_path, stringsAsFactors = FALSE)
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

grDevices::pdf(file.path(figure_dir, "sim_performance_boxplots.pdf"), width = 9.8, height = 9.6)
graphics::par(mfrow = c(3, 2), mar = c(4.2, 4.2, 2.2, 1.1))
for (scenario_name in levels(raw_results$scenario)) {
  dat <- raw_results[raw_results$scenario == scenario_name, , drop = FALSE]
  dat$runtime_plot <- pmax(dat$runtime_sec, 1e-3)
  graphics::boxplot(
    procrustes_rmse ~ method,
    data = dat,
    col = unname(method_colors[method_levels]),
    border = "grey30",
    ylab = "Procrustes RMSE",
    main = paste(scenario_name, ": geometry recovery")
  )
  graphics::boxplot(
    runtime_plot ~ method,
    data = dat,
    col = unname(method_colors[method_levels]),
    border = "grey30",
    ylab = "Runtime (sec)",
    log = "y",
    main = paste(scenario_name, ": runtime")
  )
}
graphics::mtext("Simulation performance across scenarios", outer = TRUE, line = -1.8, cex = 1.1)
grDevices::dev.off()
