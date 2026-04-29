# Core utilities
# ================

#' Log an INFO-level message
#'
#' Respects `CONFIG$info` when `show_info` is `NULL`.
#' @param message Text.
#' @param show_info Whether to print; default from `CONFIG$info`.
#' @return `NULL`, invisibly.
log_info <- function(message, show_info = NULL) {
  if (is.null(show_info)) {
    if (exists("CONFIG") && !is.null(CONFIG$info)) {
      show_info <- CONFIG$info
    } else {
      show_info <- FALSE
    }
  }
  
  if (show_info) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    msg <- sprintf("[%s] [INFO] %s", timestamp, message)
    message(msg)
    flush.console()
  }
}

#' Log a DEBUG-level message
#'
#' Respects `CONFIG$debug` when `show_debug` is `NULL`.
#' @param message Text.
#' @param show_debug Whether to print; default from `CONFIG$debug`.
#' @return `NULL`, invisibly.
log_debug <- function(message, show_debug = NULL) {
  if (is.null(show_debug)) {
    if (exists("CONFIG") && !is.null(CONFIG$debug)) {
      show_debug <- CONFIG$debug
    } else {
      show_debug <- FALSE
    }
  }
  
  if (show_debug) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    msg <- sprintf("[%s] [DEBUG] %s", timestamp, message)
    message(msg)
    flush.console()
  }
}

#' Map sex code to label (pedtools convention).
#' @param x Sex code (1 = male, 2 = female).
#' @return `"M"`, `"F"`, or `"UNK"`.
sex_label <- function(x) {
  if (length(x) == 0 || is.na(x)) return("UNK")
  if (x == 1) return("M")
  if (x == 2) return("F")
  "UNK"
}

#' Test whether object is a non-singleton pedigree.
#' @param obj R object.
#' @return Logical.
is_ped_not_singleton <- function(obj) {
  inherits(obj, "ped") && !inherits(obj, "singleton")
}

# Validation
# ==========

#' Check that a path exists and is readable
#'
#' @param file_path File path.
#' @return `TRUE` invisibly; stops otherwise.
validate_file <- function(file_path) {
  log_info("Validating file")
  log_debug(sprintf("Path: %s", file_path))
  if (is.null(file_path) || file_path == "") {
    stop("No file path provided")
  }
  if (!file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }
  if (!file.access(file_path, mode = 4) == 0) {
    stop(sprintf("File not readable: %s", file_path))
  }
  log_info("File OK")
  return(TRUE)
}

# Pedigree processing
# ====================

#' Reorder individuals in each pedigree by numeric ID order
#'
#' Trims trailing whitespace from individual IDs before sorting (see local helper).
#' @param pedigree_list List of `ped` objects.
#' @param verbose If `TRUE`, emit debug logs.
#' @return List of pedigrees with individuals reordered.
reorder_pedigrees_by_id <- function(pedigree_list, verbose = FALSE) {
  remove_trailing_space <- function(x) {
    sub("\\s$", "", x)
  }

  if (is.null(pedigree_list) || length(pedigree_list) == 0) {
    log_info("Empty pedigree list; skip reorder")
    return(pedigree_list)
  }
  
  log_info("Reordering pedigrees by numeric ID order")
  if (verbose) {
    log_debug(sprintf("Pedigrees to process: %d", length(pedigree_list)))
    log_debug("Using stringr::str_sort(numeric = TRUE)")
  }
  
  original_names <- names(pedigree_list)
  
  result <- lapply(seq_along(pedigree_list), function(i) {
    ped <- pedigree_list[[i]]
    ped_name <- original_names[i] %||% sprintf("Pedigree_%d", i)
    
    if (!is.ped(ped)) {
      if (verbose) {
        log_debug(sprintf("Skip '%s': not a valid pedigree", ped_name))
      }
      return(ped)
    }

    ped$ID <- remove_trailing_space(ped$ID)
    original_ids <- ped$ID
    n_ids <- length(original_ids)
    
    sorted_ids <- stringr::str_sort(original_ids, numeric = TRUE)
    
    order_changed <- !identical(original_ids, sorted_ids)
    
    if (verbose && order_changed) {
        log_debug(sprintf("Pedigree '%s': reordering %d IDs", ped_name, n_ids))
      if (n_ids > 0) {
        n_show <- min(5, n_ids)
        ids_before <- paste(head(original_ids, n_show), collapse = ", ")
        ids_after <- paste(head(sorted_ids, n_show), collapse = ", ")
        if (n_ids > n_show) {
          ids_before <- paste0(ids_before, ", ...")
          ids_after <- paste0(ids_after, ", ...")
        }
        log_debug(sprintf("  Before: [%s]", ids_before))
        log_debug(sprintf("  After: [%s]", ids_after))
      }
    } else if (verbose) {
      log_debug(sprintf("Pedigree '%s': %d IDs already sorted", ped_name, n_ids))
    }
    
    reorderPed(ped, sorted_ids)
  })
  
  names(result) <- original_names
  
  n_reordered <- sum(sapply(seq_along(pedigree_list), function(i) {
    if (!is.ped(pedigree_list[[i]])) return(FALSE)
    original_ids <- remove_trailing_space(pedigree_list[[i]]$ID)
    sorted_ids <- stringr::str_sort(original_ids, numeric = TRUE)
    !identical(original_ids, sorted_ids)
  }))
  
  log_info(sprintf("Reorder done: %d pedigrees", length(result)))
  if (verbose) {
    log_debug(sprintf("Stats: %d reordered, %d unchanged", 
                      n_reordered, length(result) - n_reordered))
  }
  
  return(result)
}

# Load and prepare data
# =====================

#' Read a Familias `.fam` file via `pedFamilias::readFam`.
#' @param file Path to the `.fam` file.
#' @param verbose If `TRUE`, extra logging (not passed to `readFam`).
#' @return List of pedigree objects.
read_famfile <- function(file, verbose = FALSE) {
  log_info("Reading .fam file")
  log_debug(sprintf("Path: %s", file))
  if (!requireNamespace("pedFamilias", quietly = TRUE)) {
    stop("Package 'pedFamilias' is required")
  }
  
  tryCatch({
    result <- pedFamilias::readFam(
      file, 
      useDVI = TRUE, 
      prefixAdded = "EXTRA_", 
      simplify1 = TRUE, 
      deduplicate = TRUE, 
      verbose = FALSE
    )
    log_info(sprintf("Read .fam: %d pedigrees", length(result)))
    return(result)
  }, error = function(e) {
    stop(sprintf("Error reading file %s: %s", file, conditionMessage(e)))
  })
}

#' Extract a single pedigree from a complex POI Component object.
#' @param x Nested structure from `readFam` (paths tried in order).
#' @return A `ped` object, or `NULL`.
extract_ped <- function(x) {
  # Some .fam files are parsed directly as a ped object (not nested).
  if (is_ped_not_singleton(x)) {
    return(x)
  }

  paths <- list(
    c("Reference pedigree", "_comp1"),
    c("Reference pedigree", "_comp2"),
    "_comp1",
    "_comp2"
  )

  for (p in paths) {
    obj <- purrr::pluck(x, !!!p, .default = NULL)

    if (is_ped_not_singleton(obj))
      return(obj)

    if (is.list(obj)) {
      cand <- purrr::keep(obj, is_ped_not_singleton)
      if (length(cand))
        return(cand[[1]])
    }
  }

  # Fallback: walk nested lists and return the first non-singleton ped found.
  if (is.list(x)) {
    queue <- x
    while (length(queue) > 0) {
      obj <- queue[[1]]
      queue <- queue[-1]

      if (is_ped_not_singleton(obj)) {
        return(obj)
      }

      if (is.list(obj)) {
        queue <- c(queue, obj)
      }
    }
  }

  NULL
}

#' Rename `EXTRA_1` to the missing-person label when present
#'
#' `pedFamilias::readFam(..., prefixAdded = "EXTRA_")` may label the missing person as
#' `EXTRA_1`. Before [pedtools::mergePed()], that placeholder is relabelled in two steps:
#' `EXTRA_1` -> `transient` (see `CONFIG$merge_mp_transient_label`) -> `mp_id`, so the merge key
#' is never an `EXTRA_*` label.
#' @param ped A [pedtools::ped()] object.
#' @param mp_id Missing-person label (same as `mergePed(..., by = mp_id)` and UI).
#' @param transient Internal bridge ID; default is `CONFIG$merge_mp_transient_label` or
#'   `"__KINSHIP_MP__"` if `CONFIG` is unavailable.
#' @return The pedigree with `EXTRA_1` (and optionally `transient`) relabelled to `mp_id`, or
#'   `ped` unchanged when neither applies.
rename_extra1 <- function(ped, mp_id = NULL, transient = NULL) {
  if (is.null(mp_id) || !nzchar(as.character(mp_id)[1])) {
    mp_id <- if (exists("CONFIG", inherits = TRUE) &&
                 is.character(CONFIG$mp_id) &&
                 nzchar(CONFIG$mp_id)) {
      CONFIG$mp_id
    } else {
      stop("Missing person ID is required (CONFIG$mp_id not set)", call. = FALSE)
    }
  }
  if (is.null(transient)) {
    transient <- if (exists("CONFIG", inherits = TRUE) &&
                     is.character(CONFIG$merge_mp_transient_label) &&
                     nzchar(CONFIG$merge_mp_transient_label)) {
      CONFIG$merge_mp_transient_label
    } else {
      "__KINSHIP_MP__"
    }
  }
  ped_out <- ped
  if ("EXTRA_1" %in% labels(ped_out)) {
    ped_out <- relabel(
      ped_out,
      new = transient,
      old = "EXTRA_1",
      reorder = FALSE,
      returnLabs = FALSE
    )
  }
  if (transient %in% labels(ped_out)) {
    ped_out <- relabel(
      ped_out,
      new = mp_id,
      old = transient,
      reorder = FALSE,
      returnLabs = FALSE
    )
  }
  ped_out
}

#' Count Mendelian-inconsistent markers per MPI pedigree
#'
#' Uses [pedtools::mendelianCheck()] with `remove = FALSE`. Pedigrees are not modified; the count is
#' the length of failing marker indices (intra-family checks).
#' @param mpi_list Named list of MPI [pedtools::ped()] objects.
#' @return Named integer vector: MPI family id -> number of inconsistent markers (`0` if none).
mpi_mendelian_mismatch_counts <- function(mpi_list) {
  if (is.null(mpi_list) || length(mpi_list) == 0) {
    return(structure(integer(0), names = character(0)))
  }
  n <- length(mpi_list)
  counts <- integer(n)
  nm <- names(mpi_list)
  for (i in seq_len(n)) {
    ped <- mpi_list[[i]]
    err_idx <- mendelianCheck(ped, remove = FALSE, verbose = FALSE)
    counts[i] <- as.integer(length(err_idx))
  }
  if (is.null(nm)) {
    nm <- vapply(seq_len(n), function(i) {
      as.character(famid(mpi_list[[i]]))[1]
    }, character(1))
  }
  names(counts) <- nm
  counts
}

#' Extract locus attributes from MPI pedigrees
#'
#' Uses [pedtools::getLocusAttributes()]. Clears `mutmod` on load; mutation is applied later via
#' [apply_mutation_to_locus_attributes()].
#' @param mpi List of MPI `ped` objects.
#' @param markers Marker names to include, or `NULL` for all.
#' @param verbose If `TRUE`, emit extra debug logs.
#' @return List of locus attribute lists (`mutmod` cleared before hypothesis build).
get_locus_attributes_from_ped <- function(mpi, markers = NULL, verbose = FALSE) {
  if (is.null(mpi) || length(mpi) == 0) {
    stop("MPI pedigree list is empty or NULL")
  }
  
  log_info("Getting locus attributes from pedigrees")
  log_debug(sprintf("Pedigree count: %d", length(mpi)))
  
  tryCatch({
    log_debug("Attributes: alleles, afreq, name, chrom, posMb, mutmod")
    locus_attributes <- getLocusAttributes(
      mpi,
      markers = markers,
      attribs = c("alleles", "afreq", "name", "chrom", "posMb", "mutmod"),
      checkComps = TRUE,
      simplify = TRUE
    )
    
    log_info(sprintf("Locus attributes for %d markers", length(locus_attributes)))
    
    log_debug("Clearing mutmod on load; mutation is attached in apply_mutation_to_locus_attributes")
    locus_attributes <- purrr::map(locus_attributes, ~{
      .x["mutmod"] <- list(NULL)
      .x
    })
    
    log_info(sprintf("Ready: %d markers with attributes", length(locus_attributes)))
    return(locus_attributes)
  }, error = function(e) {
    stop(sprintf("Error getting locus attributes: %s", conditionMessage(e)))
  })
}

# Mutation model and LR helpers
# =============================

#' Apply a mutation model to a pedigree (`pedtools::setMutmod`).
#'
#' @param ped A `ped` object.
#' @param mut_model One of `"none"` (effective rate 0), `"equal"`, or `"stepwise"`.
#' @param mut_rate Per-locus rate for `"equal"` / `"stepwise"` (ignored for `"none"`).
#' @param mut_range Step size parameter for `"stepwise"` (`setMutmod` argument `range`).
#' @param mut_range2 Second rate parameter for `"stepwise"` (`setMutmod` argument `rate2`).
#' @return The pedigree with mutation models attached to each marker.
apply_mutation_model <- function(ped, 
                                 mut_model = "none",
                                 mut_rate = 0.002,
                                 mut_range = NULL,
                                 mut_range2 = NULL) {
  
  if (is.null(mut_model) || mut_model == "NULL" || mut_model == "none") {
    ped <- setMutmod(ped, model = "equal", rate = 0)
  } else if (mut_model == "equal") {
    ped <- setMutmod(ped, model = "equal", rate = mut_rate)
  } else if (mut_model == "extended_stepwise" || mut_model == "stepwise") {
    if (is.null(mut_range) || is.null(mut_range2)) {
      stop("stepwise model requires mut_range and mut_range2")
    }
    ped <- setMutmod(ped, model = "stepwise", 
                     rate = mut_rate, 
                     range = mut_range, 
                     rate2 = mut_range2)
  } else {
    stop(sprintf("Unknown mutation model '%s'. Use 'none', 'equal', or 'stepwise'.", mut_model))
  }
  
  return(ped)
}

#' Attach mutation models to a `locus_attributes` list (for `setMarkers` on the merged pedigree).
#'
#' Uses a singleton `ped`, \code{setMarkers}, \code{apply_mutation_model}, and
#' \code{getLocusAttributes} so merged Ped 1 receives consistent `mutmod` entries.
#' H1 `_comp1` / `_comp2` pedigrees are not rebuilt from this list; they receive
#' \code{setMutmod} via \code{apply_mutation_model} after genotype masking.
#'
#' @param locus_attributes Output of \code{get_locus_attributes_from_ped()} (typically with `mutmod` cleared).
#' @param mut_model,mut_rate,mut_range,mut_range2 Same semantics as \code{apply_mutation_model}.
#' @return List of the same length and marker order, with \code{mutmod} set per locus.
apply_mutation_to_locus_attributes <- function(locus_attributes,
                                               mut_model = "none",
                                               mut_rate = NULL,
                                               mut_range = NULL,
                                               mut_range2 = NULL) {
  if (length(locus_attributes) == 0L) {
    return(locus_attributes)
  }
  mr <- mut_rate %||% 0.002
  tmp <- pedtools::singleton(".__tmpMUT__", sex = 1)
  tmp <- setMarkers(tmp, locusAttributes = locus_attributes, missing = 0, checkCons = TRUE)
  tmp <- apply_mutation_model(
    tmp,
    mut_model = mut_model,
    mut_rate = mr,
    mut_range = mut_range,
    mut_range2 = mut_range2
  )
  nm <- nMarkers(tmp)
  if (nm == 0L) {
    return(locus_attributes)
  }
  mk <- name(tmp, seq_len(nm))
  out <- getLocusAttributes(
    tmp,
    markers = mk,
    attribs = c("alleles", "afreq", "name", "chrom", "posMb", "mutmod"),
    checkComps = TRUE,
    simplify = TRUE
  )
  names_orig <- vapply(locus_attributes, function(x) as.character(x$name)[1], character(1))
  names_new <- vapply(out, function(x) as.character(x$name)[1], character(1))
  out_ordered <- vector("list", length(locus_attributes))
  for (i in seq_along(names_orig)) {
    idx <- match(names_orig[i], names_new)
    if (is.na(idx)) {
      stop(
        sprintf("apply_mutation_to_locus_attributes: marker '%s' missing after setMutmod", names_orig[i]),
        call. = FALSE
      )
    }
    out_ordered[[i]] <- out[[idx]]
  }
  out_ordered
}

#' Extract per-marker LR (column \code{Ped 1:Ped 2}) from \code{forrel::kinshipLR} output.
#'
#' @param lr_obj Object returned by \code{forrel::kinshipLR}.
#' @return Named list with \code{Marker} and \code{LR}, or \code{NULL}.
extract_lr_per_marker <- function(lr_obj) {
  lr_per_marker_obj <- lr_obj[["LRperMarker"]]
  if (is.null(lr_per_marker_obj)) {
    return(NULL)
  }

  marker_names_lr <- rownames(lr_per_marker_obj)
  lr_partial_lr <- NA_real_
  cn <- tryCatch(colnames(lr_per_marker_obj), error = function(e) NULL)
  if (!is.null(cn)) {
    cn_trim <- trimws(cn)
    if ("Ped 1:Ped 2" %in% cn_trim) {
      col_idx <- which(cn_trim == "Ped 1:Ped 2")[1]
      lr_partial_lr <- suppressWarnings(as.numeric(lr_per_marker_obj[, col_idx]))
    } else if (length(cn_trim) >= 2) {
      col_ok <- rep(TRUE, ncol(lr_per_marker_obj))
      for (j in seq_len(ncol(lr_per_marker_obj))) {
        col_j <- suppressWarnings(as.numeric(lr_per_marker_obj[, j]))
        col_ok[j] <- !all((abs(col_j - 1) < 1e-12) | is.infinite(col_j), na.rm = TRUE)
      }
      if (any(col_ok)) {
        col_idx <- which(col_ok)[1]
        lr_partial_lr <- suppressWarnings(as.numeric(lr_per_marker_obj[, col_idx]))
      } else {
        lr_partial_lr <- suppressWarnings(as.numeric(lr_per_marker_obj[, 1L]))
      }
    } else if (length(cn_trim) == 1L && ncol(lr_per_marker_obj) >= 1L) {
      lr_partial_lr <- suppressWarnings(as.numeric(lr_per_marker_obj[, 1L]))
    }
  } else {
    lr_partial_lr <- tryCatch(
      suppressWarnings(as.numeric(lr_per_marker_obj[, 1L])),
      error = function(e) NA_real_
    )
  }

  if (is.null(marker_names_lr) || length(marker_names_lr) == 0) {
    lr_per_marker_df <- as.data.frame(lr_per_marker_obj)
    if ("Marker" %in% names(lr_per_marker_df)) {
      marker_names_lr <- lr_per_marker_df$Marker
    } else {
      marker_names_lr <- seq_len(nrow(lr_per_marker_df))
    }
  }

  marker_names <- as.character(marker_names_lr)
  marker_names <- gsub("\u00A0", " ", marker_names)
  marker_names <- gsub("\u200B", "", marker_names)
  marker_names <- gsub("\u200C", "", marker_names)
  marker_names <- gsub("\uFEFF", "", marker_names)
  marker_names <- trimws(marker_names)

  if (length(lr_partial_lr) != length(marker_names)) {
    lr_partial <- rep(NA_real_, length(marker_names))
  } else {
    lr_partial <- lr_partial_lr
  }

  list(Marker = marker_names, LR = lr_partial)
}

#' Restrict a hypothesis to a marker subset
#'
#' Applies [pedtools::selectMarkers()] to `Ped 1` and to `Ped 2` components (`_comp1`, `_comp2`, or a single `ped`).
#' Preserves `markers_lr_single_source` when present.
#' @param hypothesis List with `Ped 1` and `Ped 2` as built by [build_pedigree_hypotheses()].
#' @param markers Character vector of marker names to keep.
#' @return Updated hypothesis list.
subset_hypothesis_markers <- function(hypothesis, markers) {
  if (length(markers) == 0L) {
    stop("subset_hypothesis_markers: at least one marker is required")
  }
  h <- hypothesis
  h[["Ped 1"]] <- selectMarkers(h[["Ped 1"]], markers = markers)
  p2 <- h[["Ped 2"]]
  if (is.list(p2) && "_comp1" %in% names(p2)) {
    if (!is.null(p2[["_comp1"]]) && is.ped(p2[["_comp1"]])) {
      p2[["_comp1"]] <- selectMarkers(p2[["_comp1"]], markers = markers)
    }
    if (!is.null(p2[["_comp2"]]) && is.ped(p2[["_comp2"]])) {
      p2[["_comp2"]] <- selectMarkers(p2[["_comp2"]], markers = markers)
    }
    h[["Ped 2"]] <- p2
  } else if (is.ped(p2)) {
    h[["Ped 2"]] <- selectMarkers(p2, markers = markers)
  }
  ms <- attr(hypothesis, "markers_lr_single_source")
  if (!is.null(ms)) {
    attr(h, "markers_lr_single_source") <- ms
  }
  h
}

#' Rebuild per-marker partial LRs in full panel order
#'
#' Markers listed in `excluded_markers` are assigned partial LR `1` (neutral factor). Others are taken from
#' `lr_detail` in marker order.
#' @param lr_detail List with `Marker` and `LR` from [extract_lr_per_marker()].
#' @param markers_panel_order Full marker order (e.g. from H0 merged pedigree).
#' @param excluded_markers Markers to force to LR 1.
#' @return List with `Marker` and `LR` vectors, or `NULL` if inputs are invalid.
merge_lr_detail_with_excluded <- function(lr_detail, markers_panel_order, excluded_markers) {
  if (is.null(lr_detail) || length(markers_panel_order) == 0L) {
    return(NULL)
  }
  ex <- unique(as.character(excluded_markers))
  nm <- as.character(lr_detail$Marker)
  val <- as.numeric(lr_detail$LR)
  lr_out <- numeric(length(markers_panel_order))
  for (i in seq_along(markers_panel_order)) {
    m <- as.character(markers_panel_order[i])
    if (m %in% ex) {
      lr_out[i] <- 1
      next
    }
    hit <- which(nm == m)
    lr_out[i] <- if (length(hit) > 0L) val[hit[1L]] else NA_real_
  }
  list(Marker = as.character(markers_panel_order), LR = lr_out)
}

#' Run [forrel::kinshipLR()] on an MPI+POI Component hypothesis
#'
#' Markers in `attr(hypothesis, "markers_lr_single_source")` (typed on one side only) are excluded from
#' the numeric likelihood run and assigned partial LR `1` so the displayed product matches the total LR.
#' @param hypothesis List with `Ped 1` (merged H0) and `Ped 2` (MPI + POI Component components).
#' @return Named list: `lr` (total), `lr_obj` ([forrel::kinshipLR()] output or `NULL`), `lr_detail`
#'   (per-marker `Marker`/`LR`), `lr_per_marker` (numeric vector), `nMarkers_uninformative`,
#'   `nMarkers_compared`.
kinship_lr_mpi_poic_hypothesis <- function(hypothesis) {
  ped1_full <- hypothesis[["Ped 1"]]
  n_m <- tryCatch(nMarkers(ped1_full), error = function(e) 0L)
  markers_panel <- if (n_m > 0L) {
    as.character(name(ped1_full, seq_len(n_m)))
  } else {
    character(0)
  }
  exclude_lr <- attr(hypothesis, "markers_lr_single_source") %||% character(0)
  exclude_lr <- intersect(as.character(exclude_lr), markers_panel)
  keep_lr <- setdiff(markers_panel, exclude_lr)

  if (length(keep_lr) == 0L && length(markers_panel) > 0L) {
    return(list(
      lr = 1,
      lr_obj = NULL,
      lr_detail = list(Marker = markers_panel, LR = rep(1, length(markers_panel))),
      lr_per_marker = rep(1, length(markers_panel)),
      nMarkers_uninformative = length(markers_panel),
      nMarkers_compared = 0L
    ))
  }

  if (length(markers_panel) == 0L) {
    lr_obj <- kinshipLR(hypothesis)
    lr_raw <- lr_obj[["LRtotal"]][["Ped 1:Ped 2"]]
    lr <- suppressWarnings(as.numeric(lr_raw))
    if (length(lr) == 0L) lr <- NA_real_
    else lr <- lr[1L]
    lr_detail <- extract_lr_per_marker(lr_obj)
    lr_per_marker <- if (is.null(lr_detail)) numeric(0) else lr_detail$LR
    nMarkers_uninformative <- sum(round(lr_per_marker, 6) == 1.0, na.rm = TRUE)
    nMarkers_compared <- length(lr_per_marker) - nMarkers_uninformative
    return(list(
      lr = lr, lr_obj = lr_obj, lr_detail = lr_detail,
      lr_per_marker = lr_per_marker,
      nMarkers_uninformative = nMarkers_uninformative,
      nMarkers_compared = nMarkers_compared
    ))
  }

  hip_lr <- if (length(exclude_lr) > 0L) {
    subset_hypothesis_markers(hypothesis, markers = keep_lr)
  } else {
    hypothesis
  }
  lr_obj <- kinshipLR(hip_lr)
  lr_raw <- lr_obj[["LRtotal"]][["Ped 1:Ped 2"]]
  lr <- suppressWarnings(as.numeric(lr_raw))
  if (length(lr) == 0L) lr <- NA_real_
  else lr <- lr[1L]

  lr_detail_core <- extract_lr_per_marker(lr_obj)
  lr_detail <- if (length(exclude_lr) > 0L) {
    merge_lr_detail_with_excluded(lr_detail_core, markers_panel, exclude_lr)
  } else {
    lr_detail_core
  }
  lr_per_marker <- if (is.null(lr_detail)) numeric(0) else lr_detail$LR
  nMarkers_uninformative <- sum(round(lr_per_marker, 6) == 1.0, na.rm = TRUE)
  nMarkers_compared <- length(lr_per_marker) - nMarkers_uninformative

  list(
    lr = lr, lr_obj = lr_obj, lr_detail = lr_detail,
    lr_per_marker = lr_per_marker,
    nMarkers_uninformative = nMarkers_uninformative,
    nMarkers_compared = nMarkers_compared
  )
}

#' Format \code{allele1/allele2} for one individual and marker from \code{getAlleles(H0)}.
#'
#' Columns are \code{<marker>.1} and \code{<marker>.2}; names are matched with whitespace-normalised stems.
#' @param alleles_mat Matrix from \code{getAlleles(ped_H0)}.
#' @param sample_id Row label (individual ID).
#' @param marker_label Marker name as in the per-marker LR table.
#' @return \code{"a/b"} or \code{"-"} if missing or incomplete.
modal_marker_genotype_string <- function(alleles_mat, sample_id, marker_label) {
  if (is.null(alleles_mat) || is.null(colnames(alleles_mat))) return("-")
  rn <- rownames(alleles_mat)
  if (is.null(rn) || !(sample_id %in% rn)) return("-")

  cn <- colnames(alleles_mat)
  mk <- gsub("\u00A0", " ", as.character(marker_label))
  mk <- gsub("\u200B|\u200C|\uFEFF", "", mk)
  mk <- trimws(gsub("\\s+", " ", mk))

  strip_spaces_lower <- function(x) {
    tolower(gsub("[[:space:]]+", "", x))
  }

  fmt <- function(c1, c2) {
    if (!c1 %in% cn || !c2 %in% cn) return("-")
    a1 <- alleles_mat[sample_id, c1, drop = TRUE]
    a2 <- alleles_mat[sample_id, c2, drop = TRUE]
    s1 <- as.character(a1)
    s2 <- as.character(a2)
    miss <- function(a, sc) {
      is.na(a) || is.na(sc) || !nzchar(sc) || sc %in% c("0", "NA")
    }
    if (miss(a1, s1) || miss(a2, s2)) return("-")
    paste0(s1, "/", s2)
  }

  # 1) Literal stem (after normalising label): <marker>.1 / .2
  c1e <- paste0(mk, ".1")
  c2e <- paste0(mk, ".2")
  if (c1e %in% cn && c2e %in% cn) {
    return(fmt(c1e, c2e))
  }

  # 2) Scan *.1 columns: space-stripped stem must match marker stem (e.g. Penta D)
  key <- strip_spaces_lower(mk)
  j1 <- grep("\\.1$", cn, perl = TRUE)
  for (j in j1) {
    stem <- sub("\\.1$", "", cn[j])
    if (strip_spaces_lower(stem) == key) {
      c1 <- cn[j]
      c2 <- sub("\\.1$", ".2", c1, fixed = TRUE)
      return(fmt(c1, c2))
    }
  }

  "-"
}

#' Classify a typed individual for modal column styling (MPI vs POI Component).
#'
#' Uses labels from \code{Ped 2} \code{_comp1} (MPI) and \code{_comp2} (POI Component).
#' IDs present only in POI Component are \code{"poic"}; others default to \code{"mpi"}
#' (including the shared missing-person ID when it appears in both).
#'
#' @param id Individual ID.
#' @param lab_mpi Character vector of labels from \code{Ped 2[["_comp1"]]}.
#' @param lab_poic Character vector of labels from \code{Ped 2[["_comp2"]]}.
#' @return \code{"mpi"} or \code{"poic"}.
classify_modal_individual_side <- function(id, lab_mpi, lab_poic) {
  id <- as.character(id)[1]
  in_m <- id %in% lab_mpi
  in_p <- id %in% lab_poic
  if (in_p && !in_m) return("poic")
  if (in_m && !in_p) return("mpi")
  if (in_m && in_p) return("mpi")
  "mpi"
}

#' Build wide data.frame for per-marker LR modal (same as shown in the modal, including total row).
#'
#' @param key Key \code{"MPI+POIc"} matching \code{mpi_poic_list}.
#' @param mpi_poic_list Hypothesis list from \code{build_pedigree_hypotheses}.
#' @param lr_details_by_key List of stored per-marker LR objects keyed like \code{mpi_poic_list}.
#' @param lr_modal_total Total LR scalar (from results row).
#' @param mut_settings List with \code{mut_model}, \code{mut_rate}, \code{mut_range}, \code{mut_range2}, \code{mp_id}.
#' @param default_mp_id Fallback missing-person ID.
#' @return A \code{data.frame}, or \code{NULL} if inputs are invalid. Attribute
#'   \code{sample_column_side} is a named character vector (\code{mpi} / \code{poic}) for
#'   genotype columns (for styling).
build_lr_modal_detail_df <- function(
    key,
    mpi_poic_list,
    lr_details_by_key,
    lr_modal_total,
    mut_settings,
    default_mp_id = NULL) {

  if (is.null(key) || !nzchar(as.character(key)[1])) return(NULL)
  stored <- lr_details_by_key[[key]]
  if (is.null(stored) || length(stored$Marker) == 0L) return(NULL)

  hypothesis <- mpi_poic_list[[key]]
  if (is.null(hypothesis)) return(NULL)

  ped2 <- hypothesis[["Ped 2"]]
  lab_mpi <- character(0)
  lab_poic <- character(0)
  if (is.list(ped2)) {
    if (!is.null(ped2[["_comp1"]]) && is.ped(ped2[["_comp1"]])) {
      lab_mpi <- labels(ped2[["_comp1"]])
    }
    if (!is.null(ped2[["_comp2"]]) && is.ped(ped2[["_comp2"]])) {
      lab_poic <- labels(ped2[["_comp2"]])
    }
  }

  if (is.null(default_mp_id) || !nzchar(as.character(default_mp_id)[1])) {
    default_mp_id <- if (exists("CONFIG", inherits = TRUE) &&
                         is.character(CONFIG$mp_id) &&
                         nzchar(CONFIG$mp_id)) {
      CONFIG$mp_id
    } else {
      stop("Missing person ID is required (CONFIG$mp_id not set)", call. = FALSE)
    }
  }
  mp_id_val <- mut_settings$mp_id %||% default_mp_id

  ped_h0 <- hypothesis[["Ped 1"]]
  if (is.null(ped_h0)) return(NULL)

  ped_h0 <- apply_mutation_model(
    ped = ped_h0,
    mut_model = mut_settings$mut_model,
    mut_rate = mut_settings$mut_rate,
    mut_range = mut_settings$mut_range,
    mut_range2 = mut_settings$mut_range2
  )

  alleles_h0 <- tryCatch(getAlleles(ped_h0), error = function(e) NULL)

  marker_names <- stored$Marker
  lr_partial <- stored$LR
  lr_total <- suppressWarnings(as.numeric(lr_modal_total[1]))

  typed_ids <- tryCatch(typedMembers(ped_h0), error = function(e) character(0))
  if (length(typed_ids) == 0 && !is.null(mp_id_val) && nzchar(as.character(mp_id_val)[1])) {
    typed_ids <- as.character(mp_id_val)[1]
  }
  sample_cols <- make.unique(as.character(typed_ids), sep = "_")

  detail_df <- data.frame(
    Marker = marker_names,
    LR = lr_partial,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!is.null(alleles_h0) && length(typed_ids) > 0) {
    for (i in seq_along(typed_ids)) {
      sid <- typed_ids[i]
      detail_df[[sample_cols[i]]] <- vapply(
        marker_names,
        function(m) modal_marker_genotype_string(alleles_h0, sid, m),
        character(1)
      )
    }
  }

  sample_side <- NULL
  if (length(sample_cols) > 0L) {
    sample_side <- vapply(
      seq_along(typed_ids),
      function(i) classify_modal_individual_side(typed_ids[i], lab_mpi, lab_poic),
      character(1)
    )
    names(sample_side) <- sample_cols
    attr(detail_df, "sample_column_side") <- sample_side
  }

  genotype_cols <- setdiff(names(detail_df), c("Marker", "LR"))
  if (length(genotype_cols) > 1L) {
    genotype_cols <- sort(genotype_cols, na.last = TRUE)
    detail_df <- detail_df[, c("Marker", "LR", genotype_cols), drop = FALSE]
    if (!is.null(sample_side)) {
      sample_side <- sample_side[genotype_cols]
      attr(detail_df, "sample_column_side") <- sample_side
    }
  }

  total_row <- detail_df[1, , drop = FALSE]
  total_row$Marker <- "Total LR"
  total_row$LR <- lr_total
  for (nm in setdiff(names(detail_df), c("Marker", "LR"))) {
    total_row[[nm]] <- ""
  }
  out <- rbind(detail_df, total_row)
  if (!is.null(sample_side)) {
    attr(out, "sample_column_side") <- sample_side
  }
  out
}

#' Build MPI\u00d7POI Component (POIc) hypothesis list for `kinshipLR` (H0 merged pedigree, H1 two components).
#' @param MPI List of MPI `ped` objects.
#' @param poic List of POI Component `ped` objects.
#' @param locus_attributes Locus list from `get_locus_attributes_from_ped()`.
#' @param mp_id Missing-person label (must match IDs on both sides).
#' @param verbose If TRUE, more `log_debug` output.
#' @param progress Optional Shiny progress object with `$set(value, detail)`.
#' @param mut_model,mut_rate,mut_range,mut_range2 Mutation model applied to `locus_attributes` before building (H0/H1).
#' @return Named list: keys `"MPI_FAM+POIc_FAM"`, values `list(Ped 1 = ..., Ped 2 = list(_comp1, _comp2))`.
build_pedigree_hypotheses <- function(MPI,
                     poic,
                     locus_attributes,
                     mp_id = NULL,
                     verbose = FALSE,
                     progress = NULL,
                     mut_model = "none",
                     mut_rate = NULL,
                     mut_range = NULL,
                     mut_range2 = NULL) {
  if (is.null(mp_id) || !nzchar(as.character(mp_id)[1])) {
    mp_id <- if (exists("CONFIG", inherits = TRUE) &&
                 is.character(CONFIG$mp_id) &&
                 nzchar(CONFIG$mp_id)) {
      CONFIG$mp_id
    } else {
      stop("Missing person ID is required (CONFIG$mp_id not set)", call. = FALSE)
    }
  }
  
  total_comparisons <- length(MPI) * length(poic)
  log_info("Building MPI–POI Component hypotheses")
  log_info(sprintf("Total comparisons: %d", total_comparisons))
  log_debug(sprintf("MPI: %d pedigrees, POI Component: %d pedigrees", length(MPI), length(poic)))
  
  locus_attributes <- tryCatch(
    apply_mutation_to_locus_attributes(
      locus_attributes,
      mut_model = mut_model,
      mut_rate = mut_rate,
      mut_range = mut_range,
      mut_range2 = mut_range2
    ),
    error = function(e) {
      stop(sprintf("Error applying mutation model to locus attributes: %s", conditionMessage(e)), call. = FALSE)
    }
  )
  log_info("Mutation models attached to locus attributes for merged Ped 1")
  
  out <- list()
  merge_failures <- list()
  sex_filtered <- character(0)
  markers <- vapply(locus_attributes, function(x) as.character(x$name)[1], character(1))
  log_debug(sprintf("Markers available: %d", length(markers)))
  
  current_comparison <- 0
  
  for (i in seq_along(MPI)) {
    ped_MPI0 <- MPI[[i]]
    fam1 <- famid(ped_MPI0) %||% paste0("MPI_", i)

    for (j in seq_along(poic)) {
      current_comparison <- current_comparison + 1
      
      if (!is.null(progress) && is.function(progress$set)) {
        tryCatch({
          progress$set(
            value = current_comparison / total_comparisons,
            detail = sprintf("%s + %s (%d/%d)", fam1, 
                            famid(poic[[j]]) %||% paste0("POIc_", j),
                            current_comparison, total_comparisons)
          )
        }, error = function(e) {
          invisible()
        })
      }
      
      ped_poic0 <- poic[[j]]
      fam2 <- famid(ped_poic0) %||% paste0("POIc_", j)

      key <- paste0(fam1, "+", fam2)

      log_debug(sprintf("Comparison %d/%d: %s", current_comparison, total_comparisons, key))

      ped_MPI <- rename_extra1(ped_MPI0, mp_id)
      ped_poic <- rename_extra1(ped_poic0, mp_id)

      ## 1. Missing-person sex compatibility filter
      ## Rule:
      ## - MPI M/F -> only merge with POI Component M/F of the same sex
      ## - MPI UNK   -> merge with both POI Component sexes
      if (mp_id %in% labels(ped_MPI) && mp_id %in% labels(ped_poic)) {
        s_MPI <- as.integer(getSex(ped_MPI, ids = mp_id)[1])
        s_poic <- as.integer(getSex(ped_poic, ids = mp_id)[1])

        mpi_known <- !is.na(s_MPI) && s_MPI %in% c(1L, 2L)
        poic_known <- !is.na(s_poic) && s_poic %in% c(1L, 2L)

        if (mpi_known && (!poic_known || s_MPI != s_poic)) {
          log_debug(sprintf(
            "Skip %s due to MP sex mismatch/unknown (MPI=%s, POIc=%s)",
            key, sex_label(s_MPI), sex_label(s_poic)
          ))
          sex_filtered <- c(sex_filtered, key)
          next
        }
      }

      ## Copies with markers for H1 (_comp1 / _comp2); same labels as mergePed (EXTRA_1 -> transient -> mp_id)
      ped_MPI_h1 <- ped_MPI
      ped_poic_h1 <- ped_poic

      ## 2. Alleles before structural merge (after EXTRA_1 -> transient -> mp_id so rownames match merge)
      alleles_MPI <- if (nMarkers(ped_MPI) > 0) getAlleles(ped_MPI) else NULL
      alleles_poic <- if (nMarkers(ped_poic) > 0) getAlleles(ped_poic) else NULL

      ## 3. Strip markers for structural merge
      if (nMarkers(ped_MPI) > 0)
        ped_MPI <- removeMarkers(ped_MPI, 1:nMarkers(ped_MPI))

      if (nMarkers(ped_poic) > 0)
        ped_poic <- removeMarkers(ped_poic, 1:nMarkers(ped_poic))

      ## 4. Structural merge
      log_debug(sprintf("Structural merge for %s", key))
      ped_merge <- tryCatch(
        mergePed(ped_MPI, ped_poic, by = mp_id),
        error = function(e) {
          msg <- conditionMessage(e)
          log_info(sprintf("mergePed failed for %s: %s", key, msg))
          merge_failures[[key]] <<- msg
          NULL
        }
      )

      if (is.null(ped_merge)) {
        next
      }

      famid(ped_merge) <- fam1

      ## 5. Attach markers and allele frequencies
      log_debug(sprintf("Setting markers for %s", key))
      ped_merge <- setMarkers(
        ped_merge,
        locusAttributes = locus_attributes,
        missing = 0,
        checkCons = TRUE
      )

      ped_merge <- selectMarkers(ped_merge, markers = markers)
      markers <- name(ped_merge, seq_len(nMarkers(ped_merge)))

      ## ---------------------------------
      ## 6. Import merged genotypes
      ## ---------------------------------
      log_debug(sprintf("Importing genotypes for %s", key))
      n_ind <- length(labels(ped_merge))
      n_loc <- length(markers)

      # Column order: marker.1, marker.2, ... (matches ped after selectMarkers)
      col_names <- character(2 * n_loc)
      for (k in seq_along(markers)) {
        col_names[2 * k - 1] <- paste0(markers[k], ".1")
        col_names[2 * k] <- paste0(markers[k], ".2")
      }

      geno_full <- matrix(
        0,
        nrow = n_ind,
        ncol = 2 * n_loc,
        dimnames = list(labels(ped_merge), col_names)
      )

      copy_geno <- function(target, source) {
        common_ids <- intersect(rownames(source), rownames(target))
        if (length(common_ids) == 0) return(target)
        
        common_cols <- intersect(colnames(target), colnames(source))
        if (length(common_cols) > 0) {
          for (id in common_ids) {
            target[id, common_cols] <- source[id, common_cols, drop = FALSE]
          }
        }
        target
      }

      # MPI first, then POI relatives (MP row not overwritten from POI)
      markers_zero_MPI <- character(0)
      if (!is.null(alleles_MPI)) {
        geno_full <- copy_geno(geno_full, alleles_MPI)
        
        mpi_ids_in_geno <- intersect(rownames(alleles_MPI), rownames(geno_full))
        if (length(mpi_ids_in_geno) > 0) {
          for (id in mpi_ids_in_geno) {
            common_cols <- intersect(colnames(alleles_MPI), colnames(geno_full))
            if (length(common_cols) > 0) {
              # Ensure MPI alleles survived the copy
              for (col in common_cols) {
                if (!is.na(alleles_MPI[id, col]) && alleles_MPI[id, col] != 0) {
                  if (is.na(geno_full[id, col]) || geno_full[id, col] != alleles_MPI[id, col]) {
                    if (verbose) {
                      message(sprintf(
                        "[WARN] Allele lost when copying MPI: %s, column %s, original=%s, copied=%s",
                        id, col, alleles_MPI[id, col], geno_full[id, col]
                      ))
                    }
                    geno_full[id, col] <- alleles_MPI[id, col]
                  }
                }
              }
            }
          }
        }
        
        # Markers that are all-zero for MPI rows after the MPI copy
        # (check marker.1 / marker.2 pairs)
        for (k in seq_along(markers)) {
          marker_name <- markers[k]
          col1 <- paste0(marker_name, ".1")
          col2 <- paste0(marker_name, ".2")
          
          if (col1 %in% colnames(geno_full) && col2 %in% colnames(geno_full)) {
            mpi_ids <- intersect(rownames(alleles_MPI), rownames(geno_full))
            if (length(mpi_ids) > 0) {
              # All MPI rows missing for this locus?
              all_zero <- all(
                (geno_full[mpi_ids, col1] == 0 | is.na(geno_full[mpi_ids, col1])) &
                (geno_full[mpi_ids, col2] == 0 | is.na(geno_full[mpi_ids, col2]))
              )
              if (all_zero) {
                markers_zero_MPI <- c(markers_zero_MPI, marker_name)
              }
            }
          }
        }
      }

      # POI Component copy next (does not overwrite the MP row)
      markers_zero_poic <- character(0)
      if (!is.null(alleles_poic)) {
        ids <- setdiff(rownames(alleles_poic), mp_id)
        if (length(ids) > 0) {
          common_cols <- intersect(colnames(geno_full), colnames(alleles_poic))
          if (length(common_cols) > 0) {
            geno_full[ids, common_cols] <- alleles_poic[ids, common_cols, drop = FALSE]
            
            # Markers that are all-zero for POI Component rows after that copy
            for (k in seq_along(markers)) {
              marker_name <- markers[k]
              col1 <- paste0(marker_name, ".1")
              col2 <- paste0(marker_name, ".2")
              
              if (col1 %in% colnames(geno_full) && col2 %in% colnames(geno_full)) {
                # All POI Component rows missing for this locus?
                all_zero <- all(
                  (geno_full[ids, col1] == 0 | is.na(geno_full[ids, col1])) &
                  (geno_full[ids, col2] == 0 | is.na(geno_full[ids, col2]))
                )
                if (all_zero) {
                  markers_zero_poic <- c(markers_zero_poic, marker_name)
                }
              }
            }
          }
        }
      }
      
      # Zero missing only on the subtree that contributes no data for that locus.
      # Previously the whole column was zeroed: if POI Component did not type the locus but MPI did,
      # MPI alleles were lost (modal and tables showed "-" for the whole row).
      mpi_row_ids <- if (!is.null(alleles_MPI)) {
        intersect(rownames(alleles_MPI), rownames(geno_full))
      } else {
        character(0)
      }
      poic_row_ids <- if (!is.null(alleles_poic)) {
        intersect(setdiff(rownames(alleles_poic), mp_id), rownames(geno_full))
      } else {
        character(0)
      }

      only_mpi <- setdiff(markers_zero_MPI, markers_zero_poic)
      only_poic <- setdiff(markers_zero_poic, markers_zero_MPI)
      both_sides <- intersect(markers_zero_MPI, markers_zero_poic)

      zero_geno_rows <- function(rows, col1, col2) {
        if (length(rows) == 0L) return(invisible())
        if (col1 %in% colnames(geno_full) && col2 %in% colnames(geno_full)) {
          geno_full[rows, col1] <- 0
          geno_full[rows, col2] <- 0
        }
      }

      for (marker_name in only_mpi) {
        zero_geno_rows(mpi_row_ids, paste0(marker_name, ".1"), paste0(marker_name, ".2"))
      }
      for (marker_name in only_poic) {
        zero_geno_rows(poic_row_ids, paste0(marker_name, ".1"), paste0(marker_name, ".2"))
      }
      for (marker_name in both_sides) {
        col1 <- paste0(marker_name, ".1")
        col2 <- paste0(marker_name, ".2")
        if (col1 %in% colnames(geno_full) && col2 %in% colnames(geno_full)) {
          geno_full[, col1] <- 0
          geno_full[, col2] <- 0
        }
      }

      # data.frame with columns in the order setAlleles expects (marker.1, marker.2, ...)
      geno_full <- as.data.frame(geno_full, stringsAsFactors = FALSE)
      
      missing_cols <- setdiff(col_names, colnames(geno_full))
      if (length(missing_cols) > 0) {
        for (col in missing_cols) {
          geno_full[[col]] <- 0
        }
      }
      geno_full <- geno_full[, col_names, drop = FALSE]
      
      # Validate alleles: invalid values -> 0 (setAlleles rejects alleles not in locus_attributes)
      for (k in seq_along(markers)) {
        marker_name <- markers[k]
        col1 <- paste0(marker_name, ".1")
        col2 <- paste0(marker_name, ".2")
        
        if (col1 %in% colnames(geno_full) && col2 %in% colnames(geno_full)) {
          marker_idx <- which(vapply(locus_attributes, function(x) as.character(x$name)[1], character(1)) == marker_name)
          if (length(marker_idx) > 0) {
            valid_alleles <- locus_attributes[[marker_idx]]$alleles
            
            validate_allele <- function(allele) {
              if (is.na(allele) || identical(allele, 0) || identical(allele, "0")) {
                return(0)
              }
              
              if (allele %in% valid_alleles) {
                return(allele)
              }
              
              allele_char <- as.character(allele)
              valid_char <- as.character(valid_alleles)
              if (allele_char %in% valid_char) {
                idx <- which(valid_char == allele_char)
                return(valid_alleles[idx[1]])
              }
              
              return(0)
            }
            
            original_col1 <- geno_full[, col1]
            original_col2 <- geno_full[, col2]
            
            geno_full[, col1] <- sapply(geno_full[, col1], validate_allele)
            geno_full[, col2] <- sapply(geno_full[, col2], validate_allele)
            
            # If either allele was invalid, zero both (avoid half-called genotypes)
            for (row_idx in seq_len(nrow(geno_full))) {
              orig_val1 <- original_col1[row_idx]
              orig_val2 <- original_col2[row_idx]
              new_val1 <- geno_full[row_idx, col1]
              new_val2 <- geno_full[row_idx, col2]
              
              allele1_was_invalid <- (!is.na(orig_val1) && orig_val1 != 0 && new_val1 == 0)
              allele2_was_invalid <- (!is.na(orig_val2) && orig_val2 != 0 && new_val2 == 0)
              
              if (allele1_was_invalid || allele2_was_invalid) {
                geno_full[row_idx, col1] <- 0
                geno_full[row_idx, col2] <- 0
              }
            }
          }
        }
      }
      

      log_debug(sprintf("setAlleles for %s", key))
      ped_merge <- setAlleles(ped_merge, alleles = geno_full)

      log_debug(sprintf("Building H0/H1 for %s", key))
      mpi_poic <- list()

      mpi_poic[["Ped 1"]] <- ped_merge

      ped_poic1 <- ped_poic_h1
      ped_poic1 <- relabel(ped_poic1, new = "POI", old = mp_id)

      # H1 masking (same rules as merged geno_full); then setMutmod on each component.
      if (nMarkers(ped_MPI_h1) > 0 && length(c(only_mpi, both_sides)) > 0) {
        alleles_comp1 <- getAlleles(ped_MPI_h1)
        for (marker_name in c(only_mpi, both_sides)) {
          col1 <- paste0(marker_name, ".1")
          col2 <- paste0(marker_name, ".2")
          if (col1 %in% colnames(alleles_comp1) && col2 %in% colnames(alleles_comp1)) {
            alleles_comp1[, col1] <- 0
            alleles_comp1[, col2] <- 0
          }
        }
        ped_MPI_h1 <- setAlleles(ped_MPI_h1, alleles = alleles_comp1)
      }

      if (nMarkers(ped_poic1) > 0 && length(c(only_poic, both_sides)) > 0) {
        alleles_comp2 <- getAlleles(ped_poic1)
        for (marker_name in c(only_poic, both_sides)) {
          col1 <- paste0(marker_name, ".1")
          col2 <- paste0(marker_name, ".2")
          if (col1 %in% colnames(alleles_comp2) && col2 %in% colnames(alleles_comp2)) {
            alleles_comp2[, col1] <- 0
            alleles_comp2[, col2] <- 0
          }
        }
        ped_poic1 <- setAlleles(ped_poic1, alleles = alleles_comp2)
      }

      if (nMarkers(ped_MPI_h1) > 0L) {
        ped_MPI_h1 <- apply_mutation_model(
          ped_MPI_h1,
          mut_model = mut_model,
          mut_rate = mut_rate %||% 0.002,
          mut_range = mut_range,
          mut_range2 = mut_range2
        )
      }
      if (nMarkers(ped_poic1) > 0L) {
        ped_poic1 <- apply_mutation_model(
          ped_poic1,
          mut_model = mut_model,
          mut_rate = mut_rate %||% 0.002,
          mut_range = mut_range,
          mut_range2 = mut_range2
        )
      }

      mpi_poic[["Ped 2"]] <- list(
        "_comp1" = ped_MPI_h1,
        "_comp2" = ped_poic1
      )

      # Markers typed on one family only: kinshipLR factorisation does not yield a clean
      # P_merge / (P_MPI * P_POIc); partial LR can be misleading — exclude from product (factor 1).
      markers_lr_single_source <- unique(c(only_poic, only_mpi))
      attr(mpi_poic, "markers_lr_single_source") <- markers_lr_single_source

      out[[key]] <- mpi_poic
      log_debug(sprintf("Hypothesis %s OK (%d/%d)", 
                       key, current_comparison, total_comparisons))
    }
  }

  if (length(merge_failures) > 0) {
    lines <- paste0("  ", names(merge_failures), ": ", unlist(merge_failures, use.names = FALSE))
    stop(sprintf(
      "mergePed failed for %d of %d MPI–POI pair(s). Check MP ID, sex, and structural compatibility.\n%s",
      length(merge_failures), total_comparisons, paste(lines, collapse = "\n")
    ), call. = FALSE)
  }

  if (length(sex_filtered) > 0) {
    log_info(sprintf(
      "Sex filter skipped %d MPI–POI pair(s) (MPI known sex requires POIc same sex).",
      length(sex_filtered)
    ))
  }

  log_info(sprintf("Hypotheses built: %d", length(out)))
  return(out)
}

# Filtering and summary stats
# ===========================

#' Filter a results table by LR threshold
#'
#' @param results Results `data.frame`.
#' @param threshold Minimum LR (rows with `LR >= threshold` or infinite LR are kept).
#' @return List with `results`, `n_before`, `n_after`, `n_filtered`.
filter_results_by_threshold <- function(results, threshold = 1) {
  log_info("Filtering results by LR threshold")
  log_debug(sprintf("Threshold: LR >= %s", threshold))
  
  if (is.null(results) || nrow(results) == 0) {
    log_info("Nothing to filter")
    return(list(
      results = results,
      n_before = 0,
      n_after = 0,
      n_filtered = 0
    ))
  }
  
  n_before <- nrow(results)
  log_debug(sprintf("Rows before filter: %d", n_before))
  
  if ("LR" %in% names(results)) {
    filtered <- results %>%
      dplyr::filter(LR >= threshold | is.infinite(LR))
  } else {
    log_debug("No LR column; returning all rows")
    filtered <- results
  }
  
  n_after <- nrow(filtered)
  n_filtered <- n_before - n_after
  
  log_info(sprintf("After filter: %d rows", n_after))
  log_debug(sprintf("Rows removed: %d", n_filtered))
  
  return(list(
    results = filtered,
    n_before = n_before,
    n_after = n_after,
    n_filtered = n_filtered
  ))
}

#' Summary statistics for a results table
#'
#' @param result_df Results `data.frame` with at least `LR`, `MPI`, `POIc`, `nMarkers`.
#' @return Named list of scalar summaries (`total_comparisons`, `unique_mpi`, `unique_poic`, `max_lr`,
#'   `min_lr`, `mean_markers`, `median_lr`).
add_summary_stats <- function(result_df) {
  log_info("Summary statistics")
  
  if (is.null(result_df) || nrow(result_df) == 0) {
    log_info("Empty results; zeroed summary")
    return(list(
      total_comparisons = 0,
      unique_mpi = 0,
      unique_poic = 0,
      max_lr = NA,
      min_lr = NA,
      mean_markers = NA,
      median_lr = NA
    ))
  }
  
  stats <- list(
    total_comparisons = nrow(result_df),
    unique_mpi = length(unique(result_df$MPI)),
    unique_poic = length(unique(result_df$POIc)),
    max_lr = max(result_df$LR, na.rm = TRUE),
    min_lr = min(result_df$LR, na.rm = TRUE),
    mean_markers = mean(result_df$nMarkers, na.rm = TRUE),
    median_lr = median(result_df$LR, na.rm = TRUE)
  )
  
  log_info(sprintf("Summary: %d comparisons", stats$total_comparisons))
  log_debug(sprintf("Unique MPI: %d, unique POIc: %d", stats$unique_mpi, stats$unique_poic))
  
  return(stats)
}

# Plotting
# ========

#' Internal neighbor IDs for pedigree plot-branch sizing (parent–child edges only).
#' When at `focal_id`, do not step to children whose co-parent is in `other_spouse_ids`.
#' @keywords internal
ped_plotbranch_neighbors <- function(x, member_int, focal_id, other_spouse_ids) {
  member_id <- x$ID[member_int]
  neighbor_ints <- integer()
  fa <- father(x, member_id)[1]
  mo <- mother(x, member_id)[1]
  if (length(fa) == 1L && !is.na(fa)) neighbor_ints <- c(neighbor_ints, internalID(x, fa))
  if (length(mo) == 1L && !is.na(mo)) neighbor_ints <- c(neighbor_ints, internalID(x, mo))
  kids <- children(x, member_id)
  for (child_id in kids) {
    if (member_id == focal_id) {
      fk <- father(x, child_id)[1]
      mk <- mother(x, child_id)[1]
      co_parent <- NA_character_
      if (length(fk) == 1L && !is.na(fk) && fk == member_id) {
        co_parent <- mother(x, child_id)[1]
      } else if (length(mk) == 1L && !is.na(mk) && mk == member_id) {
        co_parent <- father(x, child_id)[1]
      }
      if (!is.na(co_parent) && co_parent %in% other_spouse_ids)
        next
    }
    neighbor_ints <- c(neighbor_ints, internalID(x, child_id))
  }
  unique.default(neighbor_ints)
}

#' Count individuals reachable from `spouse_id` without crossing into the other spouse
#' line of `focal_id` (used to pick the "main" branch for [pedtools::plot.ped] `spouseOrder`).
#' @keywords internal
ped_branch_size_exclusive <- function(x, focal_id, spouse_id) {
  other_spouse_ids <- setdiff(spouses(x, focal_id, internal = FALSE), spouse_id)
  start_int <- internalID(x, spouse_id)
  n <- pedsize(x)
  visited <- logical(n)
  queue <- start_int
  visited[start_int] <- TRUE
  count <- 0L
  while (length(queue)) {
    visiting_int <- queue[1L]
    queue <- queue[-1L]
    count <- count + 1L
    for (neighbor_int in ped_plotbranch_neighbors(x, visiting_int, focal_id, other_spouse_ids)) {
      if (neighbor_int < 1L || neighbor_int > n) next
      if (!visited[neighbor_int]) {
        visited[neighbor_int] <- TRUE
        queue <- c(queue, neighbor_int)
      }
    }
  }
  count
}

#' Build `spouseOrder` for [pedtools::plot.ped] / [pedtools::plotPedList].
#'
#' For each individual with two or more spouses, orders spouses left-to-right so that
#' the side with the larger exclusive branch (more members when the other marriage
#' line is blocked at the focal person) comes first. Ties break on ID order.
#'
#' @param x A [pedtools::ped()] object.
#' @return `NULL` if no hints are needed, otherwise a list of ID vectors for `spouseOrder`.
spouse_order_main_branch <- function(x) {
  if (!is.ped(x))
    return(NULL)
  order_list <- list()
  for (focal_id in labels(x)) {
    spouse_ids <- spouses(x, focal_id, internal = FALSE)
    if (length(spouse_ids) < 2L)
      next
    if (length(spouse_ids) == 2L) {
      spouse_a <- spouse_ids[1L]
      spouse_b <- spouse_ids[2L]
      branch_size_a <- ped_branch_size_exclusive(x, focal_id, spouse_a)
      branch_size_b <- ped_branch_size_exclusive(x, focal_id, spouse_b)
      if (branch_size_a > branch_size_b ||
          (branch_size_a == branch_size_b && as.character(spouse_a) <= as.character(spouse_b))) {
        order_vec <- c(spouse_a, focal_id, spouse_b)
      } else {
        order_vec <- c(spouse_b, focal_id, spouse_a)
      }
    } else {
      branch_sizes <- vapply(spouse_ids, function(s) ped_branch_size_exclusive(x, focal_id, s), integer(1))
      rank_order <- order(-branch_sizes, spouse_ids)
      spouses_sorted <- spouse_ids[rank_order]
      order_vec <- c(spouses_sorted[1L], focal_id, spouses_sorted[-1L])
    }
    order_list[[length(order_list) + 1L]] <- order_vec
  }
  if (length(order_list) == 0L)
    return(NULL)
  order_list
}

#' Initial kinship2 `hints$order` vector (same idea as [kinship2::autohint] init): within each
#' depth level, subjects get order 1, 2, ... in row order.
#' @keywords internal
kinship2_order_init <- function(x) {
  if (!requireNamespace("kinship2", quietly = TRUE))
    stop("Package 'kinship2' is required for pedigree plot hints", call. = FALSE)
  k2 <- as_kinship2_pedigree(x)
  n <- length(k2$id)
  depth <- kinship2::kindepth(k2, align = TRUE)
  horder <- integer(n)
  for (lev in sort(unique(as.vector(depth)))) {
    who <- which(depth == lev & horder == 0L)
    if (length(who))
      horder[who] <- seq_len(length(who))
  }
  horder
}

#' Internal IDs of `id` plus full siblings (same father and mother in `x`).
#' @keywords internal
full_sibling_internal_ids <- function(x, id) {
  unique.default(c(internalID(x, id), internalID(x, siblings(x, id, half = FALSE))))
}

#' Replicate kinship2 autohint `shift()` for non-twin sibships: move `subject_int` to the
#' leftmost (`go_left = TRUE`) or rightmost (`go_left = FALSE`) position among `sib_ints`.
#' @keywords internal
apply_sib_shift_horizontal <- function(hint, subject_int, sib_ints, go_left) {
  sib_ints <- unique.default(as.integer(sib_ints))
  subject_int <- as.integer(subject_int)
  if (!subject_int %in% sib_ints)
    sib_ints <- c(sib_ints, subject_int)
  if (length(sib_ints) < 2L)
    return(hint)
  if (isTRUE(go_left)) {
    hint[subject_int] <- min(hint[sib_ints]) - 1L
  } else {
    hint[subject_int] <- max(hint[sib_ints]) + 1L
  }
  hint[sib_ints] <- rank(hint[sib_ints])
  hint
}

#' Adjust `hints$order` so the father of `mp_id` is rightmost among his full siblings
#' (paternal aunts/uncles to his left) and the mother of `mp_id` is leftmost among her
#' full siblings (maternal aunts/uncles to her right), matching kinship2 sibling-shift logic.
#' @keywords internal
mp_parent_sibling_order_hints <- function(x, mp_id) {
  n <- pedsize(x)
  if (!mp_id %in% labels(x))
    return(integer(n))
  horder <- kinship2_order_init(x)
  fa <- father(x, mp_id)[1]
  mo <- mother(x, mp_id)[1]
  if (length(fa) == 1L && !is.na(fa)) {
    sfi <- full_sibling_internal_ids(x, fa)
    if (length(sfi) > 1L)
      horder <- apply_sib_shift_horizontal(horder, internalID(x, fa), sfi, go_left = FALSE)
  }
  if (length(mo) == 1L && !is.na(mo)) {
    smi <- full_sibling_internal_ids(x, mo)
    if (length(smi) > 1L)
      horder <- apply_sib_shift_horizontal(horder, internalID(x, mo), smi, go_left = TRUE)
  }
  horder
}

#' Convert [spouse_order_main_branch()] output to a kinship2 `spouse` hint matrix (internal IDs).
#' Same pairing rules as pedtools `.spouseOrder()`.
#' @keywords internal
spouse_order_list_to_spouse_matrix <- function(x, plotorder_list) {
  if (is.null(plotorder_list) || length(plotorder_list) == 0L)
    return(NULL)
  if (!is.list(plotorder_list))
    plotorder_list <- list(plotorder_list)
  all_pairs <- list()
  for (ids in plotorder_list) {
    ids_int <- internalID(x, ids)
    if (anyNA(ids_int))
      next
    if (length(ids_int) < 2L)
      next
    new_pairs <- lapply(seq.int(2L, length(ids_int)), function(i) ids_int[(i - 1L):i])
    all_pairs <- c(all_pairs, new_pairs)
  }
  if (length(all_pairs) == 0L)
    return(NULL)
  valid_pairs <- list()
  for (p in all_pairs) {
    if (anyNA(p) || length(p) < 2L)
      next
    spouse_vec <- spouses(x, p[2L], internal = TRUE)
    is_spouse <- p[1L] %in% spouse_vec
    if (!isTRUE(is_spouse))
      next
    valid_pairs[[length(valid_pairs) + 1L]] <- p
  }
  if (length(valid_pairs) == 0L)
    return(NULL)
  cbind(do.call(rbind, valid_pairs), 0L)
}

#' Build a full `hints` list for [pedtools::plot.ped] (cannot combine `spouseOrder` with `hints`).
#'
#' @param x A [pedtools::ped()].
#' @param spouse_order_list Result of [spouse_order_main_branch()] or `NULL`.
#' @param mp_id Label of the missing person / POI in this plot (parents of this ID get sibling shifts).
#' @param use_spouse_order Whether to encode `spouse_order_list` into `hints$spouse`.
#' @param use_mp_parent_sibling_order Whether to shift father's and mother's full sibships.
#' @return `NULL` if nothing to pass, else `list(order = , spouse = )` for `plot(..., hints = )`.
build_ped_plot_hints <- function(
  x,
  spouse_order_list = NULL,
  mp_id = NULL,
  use_spouse_order = TRUE,
  use_mp_parent_sibling_order = TRUE
) {
  has_spouse_mat <- isTRUE(use_spouse_order) &&
    !is.null(spouse_order_list) && length(spouse_order_list) > 0L
  spouse_mat <- NULL
  if (has_spouse_mat) {
    spouse_mat <- tryCatch(
      spouse_order_list_to_spouse_matrix(x, spouse_order_list),
      error = function(e) {
        log_debug(sprintf("Ignoring spouse hint due to error: %s", conditionMessage(e)))
        NULL
      }
    )
  }
  use_mp <- isTRUE(use_mp_parent_sibling_order) && !is.null(mp_id) && mp_id %in% labels(x)
  if (!has_spouse_mat && !use_mp)
    return(NULL)
  if (use_mp) {
    horder <- mp_parent_sibling_order_hints(x, mp_id)
  } else {
    horder <- seq_along(x$ID)
  }
  list(order = horder, spouse = spouse_mat)
}

#' Plot H0–H1 pedigree layout for MPI–POI Component comparison (missing branch + POI panel).
#' @param hypothesis_elem One element from `build_pedigree_hypotheses()` (`Ped 1` / `Ped 2`).
#' @param labs Individual labels.
#' @param hatched IDs to hatch.
#' @param missing_id Missing-person ID.
#' @param POI_id POI label in H0 plot.
#' @param MP.col Fill colour for missing person.
#' @param POI.col Fill colour for POI.
#' @param titles Plot titles (expressions allowed).
#' @param cex Text scale.
#' @param width Figure width (inches).
#' @param height Figure height (inches).
#' @param mar Per-panel margins only: each of `plot_H0`, `plot_fam`, and `plot_poi` gets this as
#'   `margins =` (pedtools equivalent of `par(mar)` inside each `plot.ped` call). This is **not** passed
#'   to [pedtools::plotPedList()] and does **not** set outer margins for the combined layout—only the
#'   three pedigree panels. Vector of length 4: bottom, left, top, right (same order as [graphics::par()] `mar`).
#' @param fmar `plotPedList` frame margin (0–0.5).
#' @param auto_spouse_order If `TRUE`, encode [spouse_order_main_branch()] into `hints$spouse`
#'   so remarriages plot with the larger branch first.
#' @param auto_mp_parent_sibling_order If `TRUE`, pass `hints$order` so the father of the MP is
#'   rightmost among his full siblings (paternal aunts/uncles to his left) and the mother of the MP
#'   is leftmost among her full siblings (maternal aunts/uncles to her right).
#' @param ... Passed to `plotPedList`.
#' @return Invisible; draws plot.
missing_branch_plot <- function(
  hypothesis_elem,
  labs       = NULL,
  hatched   = NULL,
  missing_id = NULL,
  POI_id     = "POI",
  MP.col     = "#FF9999",
  POI.col    = "lightgreen",
  titles = NULL,
  cex = 1.2,
  width = 8,
  height = 6,
  mar = c(4, 4, 4, 4),
  fmar = 0.010,
  auto_spouse_order = TRUE,
  auto_mp_parent_sibling_order = TRUE,
  ...
) {
  if (is.null(missing_id) || !nzchar(as.character(missing_id)[1])) {
    missing_id <- if (exists("CONFIG", inherits = TRUE) &&
                      is.character(CONFIG$mp_id) &&
                      nzchar(CONFIG$mp_id)) {
      CONFIG$mp_id
    } else {
      stop("Missing person ID is required (CONFIG$mp_id not set)", call. = FALSE)
    }
  }
  missing_id <- as.character(missing_id)[1]

  if (is.null(titles)) {
    titles <- c(
      bquote(H[0] * ": POI component = " * .(missing_id)),
      expression(H[1] * ": POI component unrelated")
    )
  }

  # ------------------------------------------------------------
  # Checks
  # ------------------------------------------------------------
  # Ped 1: merged MPI + POI Component
  ped_fam <- hypothesis_elem[["Ped 1"]]
  
  # Ped 2: _comp1 = original MPI; _comp2 = POI Component relabelled as "POI"
  ped_MPI <- hypothesis_elem[["Ped 2"]][["_comp1"]]
  ped_poi <- hypothesis_elem[["Ped 2"]][["_comp2"]]

  if (!is.ped(ped_fam) || !is.ped(ped_MPI) || !is.ped(ped_poi)) {
    cls <- function(x) paste(class(x), collapse = "/")
    stop(
      sprintf(
        "Invalid hypothesis structure for plotting (Ped 1=%s, Ped 2._comp1=%s, Ped 2._comp2=%s)",
        cls(ped_fam), cls(ped_MPI), cls(ped_poi)
      ),
      call. = FALSE
    )
  }

  if (!missing_id %in% labels(ped_fam))
    stop("Missing person not found in merged pedigree")

  # ------------------------------------------------------------
  # Defaults
  # ------------------------------------------------------------
  if (is.null(labs))
    labs <- setNames(labels(ped_fam), labels(ped_fam))

  if (is.null(hatched))
    hatched <- typedMembers(ped_fam)

  # ------------------------------------------------------------
  # H0: MP == POI
  # ------------------------------------------------------------
  ped_H0 <- relabel(ped_fam, new = POI_id, old = missing_id)

  labs_H0 <- labs[labs != missing_id]
  labs_H0 <- c(labs_H0, setNames(POI_id, POI_id))
  labs_H0 <- labs_H0[labs_H0 %in% labels(ped_H0)]

  hatched_H0 <- setdiff(hatched, missing_id)

  spouse_order_h0 <- if (isTRUE(auto_spouse_order)) spouse_order_main_branch(ped_H0) else NULL
  plot_hints_h0 <- build_ped_plot_hints(
    ped_H0,
    spouse_order_list = spouse_order_h0,
    mp_id = POI_id,
    use_spouse_order = auto_spouse_order,
    use_mp_parent_sibling_order = auto_mp_parent_sibling_order
  )

  # Per-panel `margins` for pedtools::plot.ped (H0 panel only)—not a plotPedList() argument
  plot_H0 <- c(
    list(
      ped_H0,
      labs     = labs_H0,
      hatched  = hatched_H0,
      fill     = setNames(POI.col, POI_id),
      arrows   = FALSE,
      margins  = mar
    ),
    if (!is.null(plot_hints_h0)) list(hints = plot_hints_h0) else NULL
  )

  ped_fam_h1 <- removeIndividuals(
    ped_MPI,
    ids = missing_id,
    remove = "descendants",
    verbose = FALSE
  )

  # Keep addChildren inputs scalar; vectorised parent inputs can return a list instead of a ped.
  mp_father <- father(ped_fam, id = missing_id)[1]
  mp_mother <- mother(ped_fam, id = missing_id)[1]
  mp_sex <- as.integer(getSex(ped_fam, ids = missing_id)[1])

  ped_fam_h1 <- tryCatch(
    addChildren(
      ped_fam_h1,
      father = mp_father,
      mother = mp_mother,
      nch = 1,
      sex = mp_sex,
      ids = missing_id
    ),
    error = function(e) {
      log_debug(sprintf("Fallback to ped_MPI in H1 family panel: %s", conditionMessage(e)))
      ped_MPI
    }
  )
  if (!is.ped(ped_fam_h1)) {
    log_debug("H1 family panel is not a ped after addChildren; fallback to ped_MPI")
    ped_fam_h1 <- ped_MPI
  }

  labs_fam <- labs[labs %in% labels(ped_fam_h1)]
  names(labs_fam) <- names(labs)[labs %in% labels(ped_fam_h1)]

  hatched_fam <- intersect(hatched, labels(ped_fam_h1))
  if (length(hatched_fam) == 0) hatched_fam <- NULL

  spouse_order_fam <- if (isTRUE(auto_spouse_order)) spouse_order_main_branch(ped_fam_h1) else NULL
  plot_hints_fam <- build_ped_plot_hints(
    ped_fam_h1,
    spouse_order_list = spouse_order_fam,
    mp_id = missing_id,
    use_spouse_order = auto_spouse_order,
    use_mp_parent_sibling_order = auto_mp_parent_sibling_order
  )

  # MPI / missing-person branch panel (H1 left)
  plot_fam <- c(
    list(
      ped_fam_h1,
      labs     = labs_fam,
      hatched  = hatched_fam,
      fill     = setNames(MP.col, missing_id),
      arrows   = FALSE,
      margins  = mar
    ),
    if (!is.null(plot_hints_fam)) list(hints = plot_hints_fam) else NULL
  )


  ped_poi_h1 <- ped_poi

  hatched_poi <- intersect(hatched, labels(ped_poi_h1))
  if (length(hatched_poi) == 0) hatched_poi <- NULL

  spouse_order_poi <- if (isTRUE(auto_spouse_order)) spouse_order_main_branch(ped_poi_h1) else NULL
  plot_hints_poi <- build_ped_plot_hints(
    ped_poi_h1,
    spouse_order_list = spouse_order_poi,
    mp_id = POI_id,
    use_spouse_order = auto_spouse_order,
    use_mp_parent_sibling_order = auto_mp_parent_sibling_order
  )

  # POI Component panel (H1 right)
  plot_poi <- c(
    list(
      ped_poi_h1,
      labs     = setNames(labels(ped_poi_h1), labels(ped_poi_h1)),
      hatched  = hatched_poi,
      fill     = setNames(POI.col, POI_id),
      arrows   = FALSE,
      margins  = mar
    ),
    if (!is.null(plot_hints_poi)) list(hints = plot_hints_poi) else NULL
  )

  # Whole-figure plot size (inches); unrelated to `mar` above. plotPedList() resets par(mar) globally.
  old_par <- par(pin = c(width, height))
  on.exit(par(old_par), add = TRUE)

  plotPedList(
    list(plot_H0, plot_fam, plot_poi),
    groups = list(1, 2:3),
    titles = titles,
    cex    = cex,
    fmar   = fmar,
    arrows = FALSE,
    ...
  )
}

# Export
# ======

#' Write results to CSV
#'
#' When `summary_stats` is provided, appends a second table (`Metric`, `Value`) after a blank line.
#' @param result_df Results table.
#' @param file_path Output path.
#' @param summary_stats Optional named list of summary scalars (appended after the main table).
#' @return `TRUE` on success.
export_to_csv <- function(result_df, file_path, summary_stats = NULL) {
  log_info("Exporting to CSV")
  log_debug(sprintf("Path: %s", file_path))
  log_debug(sprintf("Rows to export: %d", nrow(result_df)))
  
  tryCatch({
    if ("LR" %in% names(result_df)) {
      result_df$LR <- signif(as.numeric(result_df$LR), digits = 5)
    }
    write.csv(result_df, file_path, row.names = FALSE, fileEncoding = "UTF-8")
    if (!is.null(summary_stats) && length(summary_stats) > 0L) {
      summary_df <- data.frame(
        Metric = names(summary_stats),
        Value = unlist(summary_stats, use.names = FALSE),
        stringsAsFactors = FALSE
      )
      cat("\n", file = file_path, append = TRUE)
      write.csv(summary_df, file_path, row.names = FALSE, fileEncoding = "UTF-8", append = TRUE)
      log_info("Summary block appended to CSV")
    }
    log_info("CSV written")
    return(TRUE)
  }, error = function(e) {
    stop(sprintf("CSV export error: %s", conditionMessage(e)))
  })
}
