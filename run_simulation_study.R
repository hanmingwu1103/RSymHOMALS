source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in source_files) {
  source(file)
}

dir.create(file.path("inst", "results"), recursive = TRUE, showWarnings = FALSE)

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

write.csv(study$raw, file = file.path("inst", "results", "simulation_raw.csv"), row.names = FALSE)
write.csv(study$summary, file = file.path("inst", "results", "simulation_summary.csv"), row.names = FALSE)
saveRDS(study, file = file.path("inst", "results", "simulation_results.rds"))

print(study$summary)
