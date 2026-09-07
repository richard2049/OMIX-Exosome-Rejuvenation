#!/usr/bin/env Rscript

messagef <- function(fmt, ...) {
  message(sprintf(fmt, ...))
}

script_file_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 0) {
    return(normalizePath("src/scripts/process_OMIX007582_Mammal40.R", mustWork = FALSE))
  }
  normalizePath(sub("^--file=", "", hit[[1]]), mustWork = FALSE)
}

default_paths <- function() {
  script_dir <- dirname(script_file_path())
  repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
  output_dir <- file.path(repo_root, "results", "omix007582_rebuild")

  list(
    repo_root = repo_root,
    idat_dir = file.path(repo_root, "data", "RAW", "data", "OMIX007582_idat"),
    metadata = file.path(repo_root, "data", "RAW", "data", "OMIX007582-02.csv"),
    output_dir = output_dir,
    technical_output = file.path(output_dir, "OMIX007582_beta_matrix_rebuilt_technical_ids.csv"),
    partial_output = file.path(output_dir, "OMIX007582_beta_matrix_rebuilt_partial_order_map.csv"),
    summary_output = file.path(output_dir, "OMIX007582_rebuild_summary.csv"),
    debug_rds = file.path(output_dir, "debug_one_sample_bet.rds")
  )
}

usage_text <- function(defaults) {
  paste(
    "Rebuild OMIX007582 Mammal40 beta values from IDAT files using sesame.",
    "",
    "This is a support workflow, not the canonical pipeline path.",
    "It rebuilds a technical-ID beta matrix from the released IDAT files and writes",
    "derived outputs into a dedicated results directory.",
    "",
    "Important guardrail:",
    "- biological sample mapping remains underdetermined in the public release;",
    "- order-based partial mapping is exploratory and disabled by default;",
    "- raw files under data/RAW are never overwritten by this script.",
    "",
    "Usage:",
    "  Rscript src/scripts/process_OMIX007582_Mammal40.R [options]",
    "",
    "Options:",
    sprintf("  --idat-dir PATH              IDAT directory [default: %s]", defaults$idat_dir),
    sprintf("  --metadata PATH              OMIX007582 metadata CSV [default: %s]", defaults$metadata),
    sprintf("  --output-dir PATH            Output directory [default: %s]", defaults$output_dir),
    "  --prep-candidates LIST       Comma-separated sesame prep candidates.",
    "                               Use '<default>' to try openSesame() without an explicit prep.",
    "                               Default: SHCDPM,<default>",
    "  --prefixes LIST              Comma-separated technical IDAT prefixes to process.",
    "                               Example: 207925070004_R01C01,207925070004_R01C02",
    "  --max-prefixes N             Process only the first N selected prefixes.",
    "                               Useful for smoke tests before a full rebuild.",
    "  --enable-partial-mapping     Also write an exploratory order-based partial map.",
    "  --metadata-id-column NAME    Metadata column to use for exploratory partial mapping.",
    "                               Default: auto (tries OriginalSampleName, OriginalSampleName.1, sample)",
    "  --save-debug-rds [PATH]      Save the first successful sesame object as RDS.",
    "  --skip-sesame-cache          Do not refresh sesameData cache before processing.",
    "  --help                       Show this help text and exit.",
    "",
    "Primary outputs:",
    "- OMIX007582_beta_matrix_rebuilt_technical_ids.csv",
    "- OMIX007582_rebuild_summary.csv",
    "- OMIX007582_beta_matrix_rebuilt_partial_order_map.csv (only when explicitly enabled)",
    sep = "\n"
  )
}

print_usage_and_exit <- function(defaults, status = 0L) {
  cat(usage_text(defaults), "\n")
  quit(save = "no", status = status)
}

parse_prep_candidates <- function(value) {
  parts <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) {
    stop("At least one prep candidate must be provided.")
  }
  parts[parts %in% c("<default>", "default", "DEFAULT")] <- NA_character_
  parts
}

parse_prefix_list <- function(value) {
  parts <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) {
    stop("At least one prefix must be provided for --prefixes.")
  }
  parts
}

parse_positive_integer <- function(value, option_name) {
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L || as.character(parsed) != as.character(value)) {
    stop(option_name, " must be a positive integer.")
  }
  parsed
}

parse_args <- function(args, defaults) {
  opts <- list(
    idat_dir = defaults$idat_dir,
    metadata = defaults$metadata,
    output_dir = defaults$output_dir,
    technical_output = defaults$technical_output,
    partial_output = defaults$partial_output,
    summary_output = defaults$summary_output,
    debug_rds = defaults$debug_rds,
    save_debug_rds = FALSE,
    refresh_sesame_cache = TRUE,
    enable_partial_mapping = FALSE,
    metadata_id_column = "auto",
    prefixes = character(0),
    max_prefixes = NA_integer_,
    prep_candidates = c("SHCDPM", NA_character_)
  )

  i <- 1L
  while (i <= length(args)) {
    raw_arg <- args[[i]]
    if (raw_arg %in% c("-h", "--help")) {
      print_usage_and_exit(defaults, status = 0L)
    }
    if (!startsWith(raw_arg, "--")) {
      stop("Unexpected positional argument: ", raw_arg, "\n\n", usage_text(defaults))
    }

    key <- sub("^--", "", raw_arg)
    value <- NULL
    if (grepl("=", key, fixed = TRUE)) {
      pieces <- strsplit(key, "=", fixed = TRUE)[[1]]
      key <- pieces[[1]]
      value <- paste(pieces[-1], collapse = "=")
    }

    consume_value <- function() {
      if (!is.null(value)) {
        return(value)
      }
      if (i >= length(args)) {
        stop("Missing value for --", key, "\n\n", usage_text(defaults))
      }
      next_arg <- args[[i + 1L]]
      if (startsWith(next_arg, "--")) {
        stop("Missing value for --", key, "\n\n", usage_text(defaults))
      }
      i <<- i + 1L
      next_arg
    }

    if (key == "idat-dir") {
      opts$idat_dir <- consume_value()
    } else if (key == "metadata") {
      opts$metadata <- consume_value()
    } else if (key == "output-dir") {
      opts$output_dir <- consume_value()
    } else if (key == "prep-candidates") {
      opts$prep_candidates <- parse_prep_candidates(consume_value())
    } else if (key == "prefixes") {
      opts$prefixes <- parse_prefix_list(consume_value())
    } else if (key == "max-prefixes") {
      opts$max_prefixes <- parse_positive_integer(consume_value(), "--max-prefixes")
    } else if (key == "enable-partial-mapping") {
      opts$enable_partial_mapping <- TRUE
    } else if (key == "metadata-id-column") {
      opts$metadata_id_column <- consume_value()
    } else if (key == "save-debug-rds") {
      opts$save_debug_rds <- TRUE
      if (!is.null(value)) {
        opts$debug_rds <- value
      } else if (i < length(args) && !startsWith(args[[i + 1L]], "--")) {
        i <- i + 1L
        opts$debug_rds <- args[[i]]
      }
    } else if (key == "skip-sesame-cache") {
      opts$refresh_sesame_cache <- FALSE
    } else {
      stop("Unknown option: --", key, "\n\n", usage_text(defaults))
    }

    i <- i + 1L
  }

  opts$technical_output <- file.path(opts$output_dir, basename(defaults$technical_output))
  opts$partial_output <- file.path(opts$output_dir, basename(defaults$partial_output))
  opts$summary_output <- file.path(opts$output_dir, basename(defaults$summary_output))
  if (isTRUE(opts$save_debug_rds)) {
    if (identical(opts$debug_rds, defaults$debug_rds)) {
      opts$debug_rds <- file.path(opts$output_dir, basename(defaults$debug_rds))
    }
  } else {
    opts$debug_rds <- NULL
  }

  opts
}

select_idat_prefixes <- function(detected_prefixes, requested_prefixes, max_prefixes) {
  selected <- detected_prefixes

  if (length(requested_prefixes) > 0) {
    missing <- setdiff(requested_prefixes, detected_prefixes)
    if (length(missing) > 0) {
      stop("Requested IDAT prefixes were not found: ", paste(missing, collapse = ", "))
    }
    selected <- requested_prefixes
  }

  if (!is.na(max_prefixes)) {
    selected <- head(selected, max_prefixes)
  }

  if (length(selected) == 0) {
    stop("No IDAT prefixes selected for processing.")
  }

  selected
}

require_bioc_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) == 0) {
    return(invisible(NULL))
  }

  stop(
    paste0(
      "Missing required R/Bioconductor packages: ", paste(missing, collapse = ", "), ".\n",
      "Install them explicitly before running this support workflow, for example:\n",
      "  if (!requireNamespace(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\")\n",
      "  BiocManager::install(c(", paste(sprintf("\"%s\"", missing), collapse = ", "), "), ask = FALSE, update = FALSE)"
    )
  )
}

refresh_sesame_resources <- function(enabled) {
  if (!isTRUE(enabled)) {
    messagef("Skipping sesameDataCache() because --skip-sesame-cache was set.")
    return(invisible(NULL))
  }
  tryCatch(
    {
      sesameData::sesameDataCache()
      messagef("sesameData cache refreshed.")
    },
    error = function(e) {
      warning("sesameDataCache() failed: ", conditionMessage(e))
    }
  )
}

list_idat_prefixes <- function(idat_dir) {
  if (!dir.exists(idat_dir)) {
    stop("IDAT directory does not exist: ", idat_dir)
  }

  idat_files <- list.files(
    idat_dir,
    pattern = "[Rr]ed\\.idat$|[Gg]rn\\.idat$",
    full.names = FALSE
  )

  if (length(idat_files) == 0) {
    stop("No IDAT files found in: ", idat_dir)
  }

  prefixes <- unique(sub("_[Rr]ed\\.idat$|_[Gg]rn\\.idat$", "", idat_files))
  prefixes <- sort(prefixes)
  messagef("Detected %d IDAT prefixes in %s", length(prefixes), idat_dir)
  prefixes
}

extract_beta_vector <- function(obj) {
  if (is.numeric(obj)) {
    return(obj)
  }
  if (!is.null(obj$beta)) {
    return(obj$beta)
  }
  if (!is.null(obj$betas)) {
    return(obj$betas)
  }
  NULL
}

process_one_sample <- function(prefix_path, prep_candidates, bpparam, debug_rds = NULL, state = NULL) {
  for (prep in prep_candidates) {
    prep_label <- if (is.na(prep)) "<default>" else prep
    messagef("Trying prep=%s for %s", prep_label, basename(prefix_path))

    sesame_obj <- tryCatch(
      {
        if (is.na(prep)) {
          sesame::openSesame(prefix_path, BPPARAM = bpparam)
        } else {
          sesame::openSesame(prefix_path, prep = prep, BPPARAM = bpparam)
        }
      },
      error = function(e) {
        warning(
          sprintf(
            "openSesame failed for %s with prep=%s: %s",
            prefix_path,
            prep_label,
            conditionMessage(e)
          )
        )
        NULL
      }
    )

    if (is.null(sesame_obj)) {
      next
    }

    if (!is.null(debug_rds) && !isTRUE(state$debug_saved)) {
      saveRDS(sesame_obj, debug_rds)
      state$debug_saved <- TRUE
      messagef("Saved first successful sesame object to %s", debug_rds)
    }

    beta_vec <- extract_beta_vector(sesame_obj)
    if (is.null(beta_vec) || length(beta_vec) == 0) {
      warning(sprintf("No usable beta vector produced for %s with prep=%s.", prefix_path, prep_label))
      next
    }

    attr(beta_vec, "prep_used") <- prep_label
    messagef("Recovered %d beta values for %s", length(beta_vec), basename(prefix_path))
    return(beta_vec)
  }

  NULL
}

build_beta_matrix <- function(prefixes, idat_dir, prep_candidates, debug_rds = NULL) {
  state <- new.env(parent = emptyenv())
  state$debug_saved <- FALSE
  bpparam <- BiocParallel::SerialParam()

  beta_list <- list()
  failed_prefixes <- character(0)

  for (prefix in prefixes) {
    prefix_path <- file.path(idat_dir, prefix)
    beta_vec <- process_one_sample(
      prefix_path = prefix_path,
      prep_candidates = prep_candidates,
      bpparam = bpparam,
      debug_rds = debug_rds,
      state = state
    )
    if (is.null(beta_vec)) {
      failed_prefixes <- c(failed_prefixes, prefix)
      next
    }
    beta_list[[prefix]] <- beta_vec
  }

  if (length(beta_list) == 0) {
    stop(
      paste(
        "No valid beta vectors were produced.",
        "Check sesame configuration, Mammal40 compatibility, and the first saved debug RDS if requested.",
        sep = "\n"
      )
    )
  }

  lengths <- vapply(beta_list, length, integer(1))
  if (length(unique(lengths)) != 1) {
    stop(
      paste0(
        "Samples produced inconsistent beta-vector lengths: ",
        paste(names(lengths), lengths, collapse = "; "),
        ". This usually indicates mixed platforms or corrupted IDATs."
      )
    )
  }

  beta_mat <- do.call(cbind, beta_list)
  colnames(beta_mat) <- names(beta_list)

  list(
    beta_matrix = beta_mat,
    failed_prefixes = failed_prefixes
  )
}

choose_metadata_id_column <- function(meta, requested = "auto") {
  if (!identical(requested, "auto")) {
    if (!requested %in% colnames(meta)) {
      stop("Requested metadata ID column not found: ", requested)
    }
    return(requested)
  }

  for (candidate in c("OriginalSampleName", "OriginalSampleName.1", "sample")) {
    if (candidate %in% colnames(meta)) {
      return(candidate)
    }
  }

  stop("Could not auto-detect a metadata ID column for exploratory partial mapping.")
}

write_partial_mapping <- function(beta_mat, metadata_path, partial_output, metadata_id_column) {
  if (!file.exists(metadata_path)) {
    warning("Metadata file not found for exploratory partial mapping: ", metadata_path)
    return(list(written = FALSE, reason = "Metadata file missing.", column = ""))
  }

  meta <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- choose_metadata_id_column(meta, metadata_id_column)
  n_map <- min(nrow(meta), ncol(beta_mat))
  if (n_map <= 0) {
    warning("Exploratory partial mapping requested, but no rows were available to map.")
    return(list(written = FALSE, reason = "No metadata rows available for mapping.", column = id_col))
  }

  mapped <- beta_mat
  colnames(mapped)[seq_len(n_map)] <- as.character(meta[[id_col]][seq_len(n_map)])
  write.csv(mapped, partial_output, row.names = TRUE, quote = FALSE)
  warning(
    paste0(
      "Exploratory partial mapping wrote ", n_map,
      " columns using metadata order and column '", id_col,
      "'. This file is non-canonical and should not replace the technical-ID matrix."
    )
  )

  list(
    written = TRUE,
    reason = sprintf("Mapped first %d columns using metadata order and column %s.", n_map, id_col),
    column = id_col
  )
}

write_summary <- function(
  summary_output,
  idat_dir,
  metadata_path,
  beta_mat,
  prefixes,
  detected_prefixes,
  failed_prefixes,
  prep_candidates,
  partial_info,
  technical_output,
  partial_output,
  prefix_filter
) {
  summary_df <- data.frame(
    idat_dir = normalizePath(idat_dir, winslash = "/", mustWork = FALSE),
    metadata_path = normalizePath(metadata_path, winslash = "/", mustWork = FALSE),
    n_idat_prefixes_detected = length(detected_prefixes),
    n_idat_prefixes_selected = length(prefixes),
    selected_prefixes = paste(prefixes, collapse = ";"),
    prefix_filter = prefix_filter,
    n_samples_rebuilt = ncol(beta_mat),
    n_failed_prefixes = length(failed_prefixes),
    failed_prefixes = paste(failed_prefixes, collapse = ";"),
    n_features = nrow(beta_mat),
    prep_candidates = paste(ifelse(is.na(prep_candidates), "<default>", prep_candidates), collapse = ","),
    technical_output = normalizePath(technical_output, winslash = "/", mustWork = FALSE),
    partial_mapping_written = isTRUE(partial_info$written),
    partial_mapping_column = partial_info$column,
    partial_mapping_reason = partial_info$reason,
    partial_output = if (isTRUE(partial_info$written)) {
      normalizePath(partial_output, winslash = "/", mustWork = FALSE)
    } else {
      ""
    },
    stringsAsFactors = FALSE
  )
  write.csv(summary_df, summary_output, row.names = FALSE, quote = TRUE)
}

main <- function() {
  defaults <- default_paths()
  opts <- parse_args(commandArgs(trailingOnly = TRUE), defaults)

  dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
  messagef("OMIX007582 rebuild output directory: %s", normalizePath(opts$output_dir, winslash = "/", mustWork = FALSE))

  require_bioc_packages(c("sesame", "sesameData", "BiocParallel"))
  suppressPackageStartupMessages({
    library(sesame)
    library(sesameData)
    library(BiocParallel)
  })

  refresh_sesame_resources(opts$refresh_sesame_cache)
  detected_prefixes <- list_idat_prefixes(opts$idat_dir)
  prefixes <- select_idat_prefixes(
    detected_prefixes = detected_prefixes,
    requested_prefixes = opts$prefixes,
    max_prefixes = opts$max_prefixes
  )
  prefix_filter <- if (length(opts$prefixes) > 0) {
    paste0("prefixes=", paste(opts$prefixes, collapse = ";"))
  } else if (!is.na(opts$max_prefixes)) {
    paste0("max-prefixes=", opts$max_prefixes)
  } else {
    "all"
  }
  messagef(
    "Selected %d of %d detected IDAT prefixes for processing.",
    length(prefixes),
    length(detected_prefixes)
  )
  rebuilt <- build_beta_matrix(
    prefixes = prefixes,
    idat_dir = opts$idat_dir,
    prep_candidates = opts$prep_candidates,
    debug_rds = opts$debug_rds
  )
  beta_mat <- rebuilt$beta_matrix

  write.csv(beta_mat, opts$technical_output, row.names = TRUE, quote = FALSE)
  messagef(
    "Wrote technical-ID beta matrix to %s (%d CpGs x %d samples).",
    opts$technical_output,
    nrow(beta_mat),
    ncol(beta_mat)
  )

  partial_info <- list(
    written = FALSE,
    reason = "Exploratory partial mapping disabled.",
    column = ""
  )
  if (isTRUE(opts$enable_partial_mapping)) {
    partial_info <- write_partial_mapping(
      beta_mat = beta_mat,
      metadata_path = opts$metadata,
      partial_output = opts$partial_output,
      metadata_id_column = opts$metadata_id_column
    )
  }

  write_summary(
    summary_output = opts$summary_output,
    idat_dir = opts$idat_dir,
    metadata_path = opts$metadata,
    beta_mat = beta_mat,
    prefixes = prefixes,
    detected_prefixes = detected_prefixes,
    failed_prefixes = rebuilt$failed_prefixes,
    prep_candidates = opts$prep_candidates,
    partial_info = partial_info,
    technical_output = opts$technical_output,
    partial_output = opts$partial_output,
    prefix_filter = prefix_filter
  )
  messagef("Wrote rebuild summary to %s", opts$summary_output)
}

main()
