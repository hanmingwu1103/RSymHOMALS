script_files <- c(
  "run_simulation_study.R",
  "run_section5_figures.R",
  "run_inference_coverage.R",
  "run_dimension_diagnostics.R",
  "run_weight_sensitivity.R",
  "run_large_scale_benchmarks.R",
  "run_real_data_examples.R"
)

for (script_file in script_files) {
  message("Running ", script_file, " ...")
  source(script_file, chdir = TRUE)
}
