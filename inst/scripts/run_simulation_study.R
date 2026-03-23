library(RSymHOMALS)

output_dir <- file.path(getwd(), "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

study <- simulate_symhomals_study(
  n_rep = 100L,
  scenarios = c("geometry", "heterogeneous", "sparse"),
  methods = c("euclidean", "wasserstein", "hausdorff", "mca"),
  ndim = 2L,
  max_iter = 70L,
  tol = 1e-5,
  seed = 20260322,
  verbose = TRUE
)

write.csv(study$raw, file.path(output_dir, "simulation_raw.csv"), row.names = FALSE)
write.csv(study$summary, file.path(output_dir, "simulation_summary.csv"), row.names = FALSE)
saveRDS(study, file.path(output_dir, "simulation_results.rds"))

print(study$summary)
