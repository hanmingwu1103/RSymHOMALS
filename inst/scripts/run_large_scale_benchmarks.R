library(RSymHOMALS)

output_dir <- file.path(getwd(), "results")
figure_dir <- file.path(getwd(), "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

designs <- data.frame(
  label = c("B1", "B2", "B3", "B4"),
  n = c(250L, 500L, 1000L, 2000L),
  m = c(5L, 8L, 10L, 12L),
  k = c(20L, 30L, 40L, 60L),
  stringsAsFactors = FALSE
)

benchmark <- analyze_scalability_benchmark(
  designs = designs,
  methods = c("euclidean", "wasserstein", "mca"),
  n_rep = 3L,
  scenario = "sparse",
  ndim = 2L,
  max_iter = 70L,
  tol = 1e-5,
  seed = 123,
  verbose = TRUE
)

utils::write.csv(benchmark$raw, file.path(output_dir, "scalability_benchmark_raw.csv"), row.names = FALSE)
utils::write.csv(benchmark$summary, file.path(output_dir, "scalability_benchmark_summary.csv"), row.names = FALSE)

plot_scalability <- function(summary_df) {
  colors <- c(euclidean = "#1b9e77", wasserstein = "#d95f02", mca = "#4d4d4d")
  methods <- c("euclidean", "wasserstein", "mca")
  method_labels <- c(euclidean = "Euclidean", wasserstein = "Wasserstein", mca = "Naive MCA")

  grDevices::pdf(file.path(figure_dir, "scalability_benchmarks.pdf"), width = 10.5, height = 4.8)
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.6, 1.2))

  plot(NULL,
       xlim = range(summary_df$nnz),
       ylim = range(summary_df$runtime_mean),
       log = "xy",
       xlab = "Total symbolic memberships (nnz)",
       ylab = "Mean runtime (sec)",
       main = "Runtime scaling")
  for (method in methods) {
    df <- summary_df[summary_df$method == method, ]
    graphics::lines(df$nnz, df$runtime_mean, type = "b", pch = 19, lwd = 2, col = colors[method])
    y0 <- pmax(df$runtime_mean - df$runtime_sd, .Machine$double.eps)
    y1 <- df$runtime_mean + df$runtime_sd
    keep <- is.finite(df$runtime_sd) & df$runtime_sd > 0
    if (any(keep)) {
      graphics::arrows(df$nnz[keep], y0[keep], df$nnz[keep], y1[keep], angle = 90, code = 3, length = 0.04, col = colors[method])
    }
  }
  graphics::legend("topleft", legend = method_labels[methods], col = colors[methods], pch = 19, lwd = 2, bty = "n")

  iter_df <- summary_df[summary_df$method %in% c("euclidean", "wasserstein"), ]
  plot(NULL,
       xlim = range(iter_df$nnz),
       ylim = range(iter_df$iterations_mean + iter_df$iterations_sd, na.rm = TRUE),
       xlab = "Total symbolic memberships (nnz)",
       ylab = "Mean iterations",
       main = "Iteration growth")
  for (method in c("euclidean", "wasserstein")) {
    df <- iter_df[iter_df$method == method, ]
    graphics::lines(df$nnz, df$iterations_mean, type = "b", pch = 19, lwd = 2, col = colors[method])
    y0 <- pmax(df$iterations_mean - df$iterations_sd, 0)
    y1 <- df$iterations_mean + df$iterations_sd
    keep <- is.finite(df$iterations_sd) & df$iterations_sd > 0
    if (any(keep)) {
      graphics::arrows(df$nnz[keep], y0[keep], df$nnz[keep], y1[keep], angle = 90, code = 3, length = 0.04, col = colors[method])
    }
  }
  graphics::legend("topleft", legend = method_labels[c("euclidean", "wasserstein")], col = colors[c("euclidean", "wasserstein")], pch = 19, lwd = 2, bty = "n")
  grDevices::dev.off()
}

plot_scalability(benchmark$summary)

print(benchmark$summary)
