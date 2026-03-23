sym_parse_set_string <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x)) {
    return("none")
  }
  x <- trimws(as.character(x))
  if (!nzchar(x)) {
    return("none")
  }
  x <- sub("^\\{", "", x)
  x <- sub("\\}$", "", x)
  x <- trimws(x)
  if (!nzchar(x)) {
    return("none")
  }
  values <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  values <- values[nzchar(values)]
  if (!length(values)) {
    return("none")
  }
  unique(values)
}

sym_example_path <- function(file = NULL) {
  root <- system.file("extdata", package = "RSymHOMALS")
  if (!nzchar(root)) {
    root <- file.path("inst", "extdata")
  }
  if (is.null(file)) {
    return(root)
  }
  file.path(root, file)
}

load_symbolic_csv <- function(path) {
  raw <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  variable_cols <- grep("^Y\\d+_", names(raw), value = TRUE)
  if (!length(variable_cols)) {
    stop("No symbolic variable columns matching '^Y\\\\d+_' were found.")
  }

  object_labels <- if ("gene" %in% names(raw)) {
    raw$gene
  } else if (all(c("patient_id", "admission_id") %in% names(raw))) {
    paste0("P", raw$patient_id, "/A", raw$admission_id)
  } else {
    as.character(seq_len(nrow(raw)))
  }

  variable_names <- gsub("^Y\\d+_", "", variable_cols)
  variable_names <- gsub("_", " ", variable_names)

  parsed_sets <- vector("list", length(variable_cols))
  category_labels <- vector("list", length(variable_cols))
  responses <- vector("list", length(variable_cols))

  for (j in seq_along(variable_cols)) {
    parsed_sets[[j]] <- lapply(raw[[variable_cols[j]]], sym_parse_set_string)
    category_labels[[j]] <- sort(unique(unlist(parsed_sets[[j]], use.names = FALSE)))
    responses[[j]] <- lapply(
      parsed_sets[[j]],
      function(values) match(values, category_labels[[j]])
    )
  }

  list(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    data = raw,
    responses = responses,
    category_labels = category_labels,
    variable_names = variable_names,
    variable_cols = variable_cols,
    object_labels = object_labels,
    n_categories = vapply(category_labels, length, integer(1)),
    dataset_name = tools::file_path_sans_ext(basename(path)),
    parsed_sets = parsed_sets
  )
}

summarize_symbolic_dataset <- function(dataset, top_n = 5L) {
  n <- nrow(dataset$data)
  set_sizes <- lapply(dataset$parsed_sets, lengths)
  variable_summary <- data.frame(
    variable = dataset$variable_names,
    n_categories = dataset$n_categories,
    mean_set_size = vapply(set_sizes, mean, numeric(1)),
    median_set_size = vapply(set_sizes, stats::median, numeric(1)),
    max_set_size = vapply(set_sizes, max, numeric(1)),
    stringsAsFactors = FALSE
  )

  label_frequency <- do.call(
    rbind,
    lapply(seq_along(dataset$parsed_sets), function(j) {
      values <- unlist(dataset$parsed_sets[[j]], use.names = FALSE)
      tab <- sort(table(values), decreasing = TRUE)
      data.frame(
        variable = dataset$variable_names[j],
        label = names(tab),
        count = as.integer(tab),
        proportion = as.numeric(tab) / n,
        stringsAsFactors = FALSE
      )
    })
  )

  total_cardinality <- rowSums(do.call(cbind, set_sizes))
  richest_idx <- order(total_cardinality, decreasing = TRUE)[seq_len(min(top_n, n))]
  richest_objects <- data.frame(
    object = dataset$object_labels[richest_idx],
    total_cardinality = total_cardinality[richest_idx],
    stringsAsFactors = FALSE
  )

  list(
    variable_summary = variable_summary,
    label_frequency = label_frequency,
    richest_objects = richest_objects
  )
}
