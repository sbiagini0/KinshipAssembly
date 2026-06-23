## ============================================================
## SERVER - SERVER-SIDE LOGIC
## ============================================================

server <- function(input, output, session) {
  
  # Empty results schema (reset / no analysis): keeps DT from showing stale rows
  empty_results_tibble <- function() {
    tibble::tibble(
      MPI = character(0),
      MPI_mismatch = integer(0),
      POIc = character(0),
      POIc_mismatch = integer(0),
      LR = numeric(0),
      POI_sex = character(0),
      nComp = integer(0),
      nMarkers = integer(0)
    )
  }
  
  # Show a success toast after the browser is likely to have finished receiving a download
  # (downloadHandler runs on the server before the file reaches the client; immediate
  # showNotification can appear too early, and CSV had no success message at all.)
  notify_after_download <- function(msg, delay_sec = 0.55) {
    session$onFlushed(function() {
      later::later(function() {
        # Pass `session` explicitly: no reactive domain inside `later()`.
        showNotification(msg, type = "message", duration = 3, session = session)
      }, delay = delay_sec)
    }, once = TRUE)
  }
  
  # Reactive values container used throughout the app
  values <- reactiveValues(
    mpi_data = NULL,
    poic_data = NULL,
    locus_attributes = NULL,
    mpi_poic_list = NULL,  # from build_pedigree_hypotheses()
    results = NULL,        # full LR table (unfiltered); UI filters by input$lr_threshold
    data_loaded = FALSE,
    locus_loaded = FALSE,
    analysis_complete = FALSE,
    analysis_running = FALSE,
    stop_calculation = FALSE,
    lr_modal_key = NULL,
    lr_modal_total = NULL,
    last_lr_modal_mutation = NULL,
    lr_details_by_key = NULL,     # Per MPI+POIc key: Marker + partial LR from last run
    # Canonical "MPI+POIc" key for the detail modal (pedigree + per-marker table); set when a results row is selected
    pedigree_plot_key = NULL,
    mpi_load_error = NULL,
    poic_load_error = NULL,
    mpi_mendelian_mismatch = NULL,  # named int: MPI fam id -> N markers failing mendelianCheck (remove=FALSE)
    poic_mendelian_mismatch = NULL,  # named int: POIc fam id -> N markers failing mendelianCheck (remove=FALSE)
    file_input_revision = 0L      # increment on reset to recreate fileInputs (browser keeps stale names otherwise)
  )
  
  # Initial button state
  shinyjs::disable("run_analysis")
  shinyjs::disable("stop_analysis")
  
  # Download report UI: real downloadButton only when there is a results table.
  # shinyjs::disable does not block Shiny download links (<a>); conditional UI avoids a live href.
  output$download_report_ui <- renderUI({
    has_table <- isTRUE(values$analysis_complete) &&
      !is.null(values$results) &&
      is.data.frame(values$results) &&
      nrow(values$results) > 0L &&
      !isTRUE(values$analysis_running)
    if (has_table) {
      downloadButton(
        "downloadReport",
        label = tagList(icon("file-download"), " Download report"),
        class = "btn-primary btn-action",
        style = "width: 100%;"
      )
    } else {
      tags$button(
        type = "button",
        class = "btn btn-primary btn-action",
        disabled = "disabled",
        style = "width: 100%; opacity: 0.45; cursor: not-allowed;",
        icon("file-download"),
        " Download report"
      )
    }
  })
  
  # ============================================
  # LOCUS ATTRIBUTES FROM MPI PEDIGREES
  # ============================================
  # Locus attributes are automatically derived from MPI pedigrees using
  # getLocusAttributes() when the MPI file is loaded.
  
  # ============================================
  # OBSERVERS AND EVENTS
  # ============================================
  
  # Helper: build joint pedigrees (called from Compute LR, not on file load / config change)
  build_pedigrees <- function() {
    if (is.null(values$mpi_data) || is.null(values$poic_data) || is.null(values$locus_attributes)) {
      return(FALSE)
    }
    
    log_info("Building combined pedigrees")
    
    tryCatch({
      mp_id_val <- if (!is.null(input$mp_id) && input$mp_id != "") {
        input$mp_id
      } else {
        CONFIG$mp_id
      }
      
      # Ensure a mutation model has been selected
      if (is.null(input$mut_model) || input$mut_model == "") {
        showNotification("Select a mutation model to continue.", 
                         type = "error", duration = 5)
        return(FALSE)
      }
      
      mut_model_val <- input$mut_model
      
      # Obtain parameters according to the selected model
      mut_rate_val <- if (mut_model_val %in% c("equal", "stepwise")) {
        if (!is.null(input$mut_rate)) {
          input$mut_rate
        } else {
          CONFIG$mut_rate
        }
      } else {
        NULL
      }
      
      mut_range_val <- if (mut_model_val == "stepwise") {
        if (!is.null(input$mut_range)) {
          input$mut_range
        } else {
          0.5  # default value
        }
      } else {
        NULL
      }
      
      mut_range2_val <- if (mut_model_val == "stepwise") {
        if (!is.null(input$mut_range2)) {
          input$mut_range2
        } else {
          0.1  # default value
        }
      } else {
        NULL
      }
      
      mp_sex_unknown_val <- isTRUE("yes" %in% input$mp_sex_unknown)
      
      comparison_mode_val <- if (!is.null(input$comparison_mode)) {
        input$comparison_mode
      } else {
        "all"
      }
      mpi_for_build <- values$mpi_data
      poic_for_build <- values$poic_data
      comparison_pairs <- NULL
      if (identical(comparison_mode_val, "custom")) {
        if (is.null(input$selected_mpi) || length(input$selected_mpi) == 0L ||
            is.null(input$selected_poic) || length(input$selected_poic) == 0L) {
          showNotification(
            "Select at least one MPI and one POI Component for custom comparisons.",
            type = "error", duration = 5
          )
          return(FALSE)
        }
        mpi_for_build <- subset_named_ped_list(values$mpi_data, input$selected_mpi)
        poic_for_build <- subset_named_ped_list(values$poic_data, input$selected_poic)
        if (length(mpi_for_build) == 0L || length(poic_for_build) == 0L) {
          showNotification(
            "No valid MPI or POI Component selection for custom comparisons.",
            type = "error", duration = 5
          )
          return(FALSE)
        }
        comparison_pairs <- expand.grid(
          mpi = names(mpi_for_build),
          poic = names(poic_for_build),
          stringsAsFactors = FALSE
        )
        log_info(sprintf(
          "Custom comparison mode: building %d MPI x %d POI Component = %d pair(s): %s",
          length(mpi_for_build), length(poic_for_build), nrow(comparison_pairs),
          paste(paste0(comparison_pairs$mpi, "+", comparison_pairs$poic), collapse = ", ")
        ))
      }
      
      withProgress(message = "Building pedigrees...", value = 0, {
        progress_obj <- list(
          set = function(value, detail = "") {
            setProgress(value = value, detail = detail)
          }
        )
        
        # Mutation models are set on locus attributes inside build_pedigree_hypotheses;
        # H1 components receive setMutmod via apply_mutation_model after allele edits.
        setProgress(value = 0.1, detail = "Building pedigrees (10%)...")
        
        values$mpi_poic_list <- tryCatch({
          build_pedigree_hypotheses(
            MPI = mpi_for_build,
            poic = poic_for_build,
            locus_attributes = values$locus_attributes,
            mp_id = mp_id_val,
            verbose = CONFIG$info,
            progress = progress_obj,
            mut_model = mut_model_val,
            mut_rate = mut_rate_val,
            mut_range = mut_range_val,
            mut_range2 = mut_range2_val,
            mp_sex_unknown = mp_sex_unknown_val,
            comparison_pairs = comparison_pairs
          )
        }, error = function(e) {
          stop(sprintf("Error building pedigrees: %s", conditionMessage(e)))
        })
        
        setProgress(value = 1, detail = paste("Done —", length(values$mpi_poic_list), "combinations (100%)"))
      })
      
      values$lr_details_by_key <- NULL
      
      if (is.null(values$mpi_poic_list) || length(values$mpi_poic_list) == 0L) {
        showNotification("No pedigree combinations could be built with the current settings.", 
                         type = "warning", duration = 5)
        return(FALSE)
      }
      
      log_info(sprintf("Pedigree construction completed: %d combinations generated", length(values$mpi_poic_list)))
      TRUE
      
    }, error = function(e) {
      showNotification(paste("Error building pedigrees:", conditionMessage(e)), 
                       type = "error", duration = 5)
      values$mpi_poic_list <- NULL
      FALSE
    })
  }
  
  # Drop cached hypotheses/results when analysis settings change (rebuilt on Compute LR)
  observeEvent(
    list(
      input$mut_model, input$mut_rate, input$mut_range, input$mut_range2,
      input$mp_id, input$mp_sex_unknown, input$comparison_mode,
      input$selected_mpi, input$selected_poic
    ),
    {
      if (!isTRUE(values$data_loaded)) return()
      values$mpi_poic_list <- NULL
      values$lr_details_by_key <- NULL
      values$results <- NULL
      values$analysis_complete <- FALSE
    },
    ignoreInit = TRUE
  )
  
  # Load POI Component file SECOND, with progress reporting and error handling
  # Input id includes file_input_revision so reset can mount fresh controls (clears filename in browser)
  observeEvent(
    {
      rid <- values$file_input_revision
      input[[paste0("file_poic_", rid)]]
    },
    {
      rid <- values$file_input_revision
      poic <- input[[paste0("file_poic_", rid)]]
      req(poic)
      
    tryCatch({
      values$poic_load_error <- NULL
      file_ext <- tools::file_ext(poic$name)
      if (nzchar(file_ext) && tolower(file_ext) == "rds") {
        stop("RDS import is not supported. Please upload a .fam file.")
      }
      
      # Load .fam file with progress bar
        withProgress(message = "Loading POI Component file...", value = 0, {
          # Step 1: validate file
          setProgress(value = 0.1, detail = "Validating file (10%)...")
          validate_file(poic$datapath)
          
          # Step 2: read .fam file
          setProgress(value = 0.3, detail = "Reading file (30%)...")
          log_info("Reading POI Component file")
          poic_raw <- read_famfile(poic$datapath, verbose = CONFIG$info)
          
          # Step 3: extract pedigrees from complex POI Component structure
          setProgress(value = 0.5, detail = "Extracting pedigrees (50%)...")
          poic_extracted <- purrr::map(poic_raw, extract_ped) |> purrr::compact()
          
          # Step 4: assign FAMIDs and trim whitespace in IDs
          setProgress(value = 0.65, detail = "Assigning FAMIDs (65%)...")
          for (i in seq_along(poic_extracted)) {
            famid(poic_extracted[[i]]) <- names(poic_extracted)[i]
          }
          
          # Step 4.5: reorder pedigrees by numeric ID
          setProgress(value = 0.75, detail = "Reordering pedigrees by ID (75%)...")
          log_info("Reordering POI Component pedigrees by ID")
          log_debug(sprintf("Pedigrees to reorder: %d", length(poic_extracted)))
          
          poic_extracted <- reorder_pedigrees_by_id(poic_extracted, verbose = CONFIG$debug)
          
          # Optionally verify reordering by logging the first few IDs
          if (length(poic_extracted) > 0) {
            first_ped <- poic_extracted[[1]]
            if (is.ped(first_ped) && length(first_ped$ID) > 0) {
              first_ids <- head(first_ped$ID, 5)
              log_debug(sprintf("First IDs of first pedigree after reorder: %s", 
                                paste(first_ids, collapse = ", ")))
            }
          }
          
          # Step 5: harmonise markers across pedigrees
          setProgress(value = 0.85, detail = "Harmonising markers (85%)...")
          poic_harmonised <- harmoniseMarkers(poic_extracted)
          
          # Step 7: finish processing
          setProgress(value = 0.95, detail = "Finalising (95%)...")
          
          values$poic_data <- poic_harmonised
          values$poic_mendelian_mismatch <- pedigree_mendelian_mismatch_counts(values$poic_data)
          n_poic_mm <- values$poic_mendelian_mismatch
          log_info(sprintf(
            "Mendelian check (POI Component): %d pedigrees, %d with >=1 inconsistent marker",
            length(n_poic_mm), sum(n_poic_mm > 0L, na.rm = TRUE)
          ))
          for (nm in names(n_poic_mm)) {
            if (n_poic_mm[[nm]] > 0L) {
              log_info(sprintf("  POI Component %s: %d inconsistent marker(s) (intrafamily)", nm, n_poic_mm[[nm]]))
            }
          }
          n_poic <- length(values$poic_data)
          
          # Step 8: done (100%)
          setProgress(value = 1, detail = paste("Done —", n_poic, "pedigrees (100%)"))
        })
        
        showNotification(paste("POI Component file loaded:", n_poic, "pedigrees"), 
                         type = "message", duration = 3)
      
      # Derive locus_attributes from MPI pedigrees using getLocusAttributes
      if (!is.null(values$mpi_data) && length(values$mpi_data) > 0) {
        tryCatch({
          log_info("Getting locus attributes from MPI pedigrees")
          values$locus_attributes <- get_locus_attributes_from_ped(values$mpi_data, markers = NULL, verbose = CONFIG$info)
          values$locus_loaded <- TRUE
        }, error = function(e) {
          stop(sprintf("Error getting locus attributes from pedigrees: %s", conditionMessage(e)))
        })
      }
      
      # When both files are loaded, ready for analysis (pedigrees built on Compute LR)
      if (!is.null(values$mpi_data) && !is.null(values$poic_data) && !is.null(values$locus_attributes)) {
        values$data_loaded <- TRUE
        log_info("Both files loaded; ready for analysis")
      }
      
    }, error = function(e) {
      # Collect full error message
      error_msg <- conditionMessage(e)
      if (length(error_msg) == 0 || nchar(error_msg) == 0) {
        error_msg <- "Unknown error while loading file"
      }
      
      values$poic_load_error <- error_msg
      
      # Notification with full error
      showNotification(
        HTML(paste(
          tags$strong("Error loading POI Component file:"), br(),
          tags$code(error_msg)
        )), 
        type = "error", 
        duration = 15,
        closeButton = TRUE
      )
      
      message("Error loading POI Component file: ", error_msg)
      
      values$poic_data <- NULL
      values$poic_mendelian_mismatch <- NULL
      values$data_loaded <- FALSE
    })
    },
    ignoreInit = TRUE
  )
  
  # Load MPI file FIRST, with progress reporting and error handling
  observeEvent(
    {
      rid <- values$file_input_revision
      input[[paste0("file_mpi_", rid)]]
    },
    {
      rid <- values$file_input_revision
      mpi <- input[[paste0("file_mpi_", rid)]]
      req(mpi)
      
    tryCatch({
      values$mpi_load_error <- NULL
      file_ext <- tools::file_ext(mpi$name)
      if (nzchar(file_ext) && tolower(file_ext) == "rds") {
        stop("RDS import is not supported. Please upload a .fam file.")
      }
      
      # Load .fam file with progress bar
        withProgress(message = "Loading MPI file...", value = 0, {
          # Step 1: validate file (5%)
          setProgress(value = 0.05, detail = "Validating file (5%)...")
          validate_file(mpi$datapath)
          
          # Step 2: read .fam file (20%)
          setProgress(value = 0.2, detail = "Reading .fam (20%)...")
          log_info("Reading MPI file")
          mpi_raw <- read_famfile(mpi$datapath, verbose = CONFIG$info)
          
          # Step 3: filter pedigrees using exclusion pattern
          setProgress(value = 0.3, detail = "Filtering pedigrees (30%)...")
          selected <- !grepl(CONFIG$exclude_patterns, names(mpi_raw))
          mpi_filtered <- mpi_raw[selected]
          
          # Step 4: extract `_comp1` from each pedigree
          setProgress(value = 0.4, detail = "Extracting components (40%)...")
          mpi_extracted <- purrr::map(mpi_filtered, ~.x[['_comp1']])
          mpi_extracted <- purrr::compact(mpi_extracted)
          
          # Step 5: assign FAMIDs and trim whitespace in IDs
          setProgress(value = 0.5, detail = "Assigning FAMIDs (50%)...")
          for (i in seq_along(mpi_extracted)) {
            famid(mpi_extracted[[i]]) <- names(mpi_extracted)[i]
          }
          
          # Step 5.5: reorder pedigrees by numeric ID
          setProgress(value = 0.55, detail = "Reordering pedigrees by ID (55%)...")
          log_info("Reordering MPI pedigrees by ID")
          log_debug(sprintf("Pedigrees to reorder: %d", length(mpi_extracted)))
          
          mpi_extracted <- reorder_pedigrees_by_id(mpi_extracted, verbose = CONFIG$debug)
          
          # Optionally verify reordering by logging the first few IDs
          if (length(mpi_extracted) > 0) {
            first_ped <- mpi_extracted[[1]]
            if (is.ped(first_ped) && length(first_ped$ID) > 0) {
              first_ids <- head(first_ped$ID, 5)
              log_debug(sprintf("First IDs of first pedigree after reorder: %s", 
                                paste(first_ids, collapse = ", ")))
            }
          }
          
          # Step 6: harmonise markers across pedigrees
          setProgress(value = 0.8, detail = "Harmonising markers (80%)...")
          mpi_harmonised <- harmoniseMarkers(mpi_extracted)
          
          # Step 7: finish processing
          setProgress(value = 0.9, detail = "Finalising (90%)...")
          
          values$mpi_data <- mpi_harmonised
          values$mpi_mendelian_mismatch <- pedigree_mendelian_mismatch_counts(values$mpi_data)
          n_mm <- values$mpi_mendelian_mismatch
          log_info(sprintf(
            "Mendelian check (MPI): %d pedigrees, %d with >=1 inconsistent marker",
            length(n_mm), sum(n_mm > 0L, na.rm = TRUE)
          ))
          for (nm in names(n_mm)) {
            if (n_mm[[nm]] > 0L) {
              log_info(sprintf("  MPI %s: %d inconsistent marker(s) (intrafamily)", nm, n_mm[[nm]]))
            }
          }
          n_mpi <- length(values$mpi_data)
          
          # Step 8: done (100%)
          setProgress(value = 1, detail = paste("Done —", n_mpi, "pedigrees (100%)"))
        })
        
        showNotification(paste("MPI file loaded:", n_mpi, "pedigrees"), 
                         type = "message", duration = 3)
      
      # Derive locus_attributes from MPI pedigrees using getLocusAttributes
      if (!is.null(values$mpi_data) && length(values$mpi_data) > 0) {
        tryCatch({
          log_info("Getting locus attributes from MPI pedigrees")
          values$locus_attributes <- get_locus_attributes_from_ped(values$mpi_data, markers = NULL, verbose = CONFIG$info)
          values$locus_loaded <- TRUE
        }, error = function(e) {
          stop(sprintf("Error getting locus attributes from pedigrees: %s", conditionMessage(e)))
        })
      }
      
      # When both files are loaded, ready for analysis (pedigrees built on Compute LR)
      if (!is.null(values$mpi_data) && !is.null(values$poic_data) && !is.null(values$locus_attributes)) {
        values$data_loaded <- TRUE
        log_info("Both files loaded; ready for analysis")
      }
      
    }, error = function(e) {
      # Collect full error message
      error_msg <- conditionMessage(e)
      if (length(error_msg) == 0 || nchar(error_msg) == 0) {
        error_msg <- "Unknown error while loading file"
      }
      
      values$mpi_load_error <- error_msg
      
      # Notification with full error
      showNotification(
        HTML(paste(
          tags$strong("Error loading MPI file:"), br(),
          tags$code(error_msg)
        )), 
        type = "error", 
        duration = 15,
        closeButton = TRUE
      )
      
      message("Error loading MPI file: ", error_msg)
      
      values$mpi_data <- NULL
      values$mpi_mendelian_mismatch <- NULL
      values$data_loaded <- FALSE
    })
    },
    ignoreInit = TRUE
  )
  
  # Observer to request stopping a long-running LR computation
  observeEvent(input$stop_analysis, {
    values$stop_calculation <- TRUE
    values$analysis_running <- FALSE
    shinyjs::disable("stop_analysis")
    showNotification("Calculation stopped.", type = "warning", duration = 3)
  })
  
  # Observe data and state to enable/disable UI buttons
  observe({
    mpi_loaded <- !is.null(values$mpi_data)
    poic_loaded <- !is.null(values$poic_data)
    
    # Validate comparison selection when using custom mode
    comparison_ok <- TRUE
    if (!is.null(input$comparison_mode) && input$comparison_mode == "custom") {
      comparison_ok <- !is.null(input$selected_mpi) && length(input$selected_mpi) > 0 &&
        !is.null(input$selected_poic) && length(input$selected_poic) > 0
    }
    
    mut_ok <- !is.null(input$mut_model) && nzchar(input$mut_model)
    
    if (values$data_loaded && !values$analysis_running && mpi_loaded && poic_loaded && comparison_ok && mut_ok) {
      shinyjs::enable("run_analysis")
    } else {
      shinyjs::disable("run_analysis")
    }
    
    if (values$analysis_running) {
      shinyjs::enable("stop_analysis")
    } else {
      shinyjs::disable("stop_analysis")
    }
  })
  
  # UI for selecting specific MPI × POI Component comparisons when using custom mode
  output$comparison_select_ui <- renderUI({
    req(input$comparison_mode)
    
    if (input$comparison_mode != "custom") {
      return(NULL)
    }
    
    # If we do not have data yet, show a short help text
    if (is.null(values$mpi_data) || is.null(values$poic_data)) {
      return(
        tags$p(class = "help-text",
               "Load MPI and POI Component files to choose specific comparisons.")
      )
    }
    
    mpi_choices <- names(values$mpi_data)
    poic_choices <- names(values$poic_data)
    
    tagList(
      shinyWidgets::pickerInput(
        inputId = "selected_mpi",
        label = tags$span(icon("user"), " Select MPI"),
        choices = mpi_choices,
        multiple = TRUE,
        options = list(
          `actions-box` = TRUE,
          `live-search` = TRUE,
          size = 10
        )
      ),
      shinyWidgets::pickerInput(
        inputId = "selected_poic",
        label = tags$span(icon("users"), " Select POI Component"),
        choices = poic_choices,
        multiple = TRUE,
        options = list(
          `actions-box` = TRUE,
          `live-search` = TRUE,
          size = 10
        )
      ),
      tags$p(class = "help-text",
             "Select at least one MPI and one POI Component; only those pairs are built and compared.")
    )
  })
  
  # Run comparisons and compute LR values
  observeEvent(input$run_analysis, {
    req(values$data_loaded)
    
    log_info("Starting MPI × POI Component comparison analysis")
    
    # Validate that a mutation model has been selected
    if (is.null(input$mut_model) || input$mut_model == "") {
      showNotification("Select a mutation model to continue.", 
                       type = "error", duration = 5)
      return(NULL)
    }
    
    values$stop_calculation <- FALSE
    values$analysis_running <- TRUE
    values$analysis_complete <- FALSE
    values$lr_details_by_key <- NULL
    values$results <- NULL
    
    if (!isTRUE(build_pedigrees())) {
      values$analysis_running <- FALSE
      return(NULL)
    }
    
    # Mutation model parameters
    mut_model_val <- input$mut_model
    if (mut_model_val == "none") {
      mut_rate_val <- 0
      mut_range_val <- NULL
      mut_range2_val <- NULL
    } else if (mut_model_val == "equal") {
      mut_rate_val <- if (!is.null(input$mut_rate)) input$mut_rate else 0.002
      mut_range_val <- NULL
      mut_range2_val <- NULL
    } else if (mut_model_val == "stepwise") {
      mut_rate_val <- if (!is.null(input$mut_rate)) input$mut_rate else 0.002
      mut_range_val <- if (!is.null(input$mut_range)) input$mut_range else 0.1
      mut_range2_val <- if (!is.null(input$mut_range2)) input$mut_range2 else 0.000001
    } else {
      mut_rate_val <- NULL
      mut_range_val <- NULL
      mut_range2_val <- NULL
    }
    
    # mpi_poic_list already contains only the pairs built (all or custom subset)
    mpi_poic_to_compare <- values$mpi_poic_list
    
    if (length(mpi_poic_to_compare) == 0L) {
      showNotification("No comparisons to run.", type = "warning", duration = 5)
      values$analysis_running <- FALSE
      return()
    }
    
    progress <- Progress$new(session, min = 0, max = 1)
    progress$set(message = "Starting analysis...", value = 0)
    
    tryCatch({
      mp_id_val <- if (!is.null(input$mp_id) && input$mp_id != "") {
        input$mp_id
      } else {
        CONFIG$mp_id
      }

        # Snapshot of mutation settings for the per-marker modal
        values$last_lr_modal_mutation <- list(
          mut_model = mut_model_val,
          mut_rate = mut_rate_val,
          mut_range = mut_range_val,
          mut_range2 = mut_range2_val,
          mp_id = mp_id_val
        )
      
      if (length(mpi_poic_to_compare) == 0) {
        progress$close()
        values$analysis_running <- FALSE
        showNotification("No data to compare.", type = "error", duration = 5)
        return()
      }
      
      # Mutation parameters are fixed at pedigree build time (build_pedigree_hypotheses).
      log_info("Computing likelihood ratios")
      progress$set(message = "Computing LR...", value = 0.1, detail = "Step 1/2")
      
      total_comparisons <- length(mpi_poic_to_compare)
      current_comparison <- 0
      
      # Compute LR for each comparison
      results_list <- vector("list", total_comparisons)
      lr_details_by_key <- list()
      lr_compute_errors <- character()
      for (key in names(mpi_poic_to_compare)) {
        if (values$stop_calculation) break
        
        current_comparison <- current_comparison + 1
        tryCatch({
          progress$set(
            value = 0.1 + (current_comparison / total_comparisons) * 0.75,
            detail = sprintf("Computing LR for %s (%d/%d)", key, current_comparison, total_comparisons)
          )
          
          # kinship_lr_mpi_poic_hypothesis() in logic/functions.R
          hypothesis <- mpi_poic_to_compare[[key]]
          kr <- kinship_lr_mpi_poic_hypothesis(hypothesis)
          lr <- kr$lr
          lr_obj <- kr$lr_obj
          lr_detail <- kr$lr_detail
          lr_per_marker <- kr$lr_per_marker
          nMarkers_uninformative <- kr$nMarkers_uninformative
          nMarkers_compared <- kr$nMarkers_compared
          if (!is.null(lr_detail)) {
            lr_details_by_key[[key]] <- lr_detail
          }
          
          # Parse MPI and POIc identifiers from the key
          parts <- strsplit(key, "+", fixed = TRUE)[[1]]
          fam1 <- parts[1]
          fam2 <- paste(parts[-1], collapse = "+")
          
          # Extract metadata from Ped 2 _comp2 (POI Component relabelled to "POI")
          ped2 <- hypothesis[["Ped 2"]]
          poi_sex_val <- NA_integer_
          ncomp_val <- NA_integer_
          if (is.list(ped2) && !is.null(ped2[["_comp2"]]) && is.ped(ped2[["_comp2"]])) {
            comp2 <- ped2[["_comp2"]]  # POI Component relabelled as "POI"
            
            # POI sex: obtained from "POI" in _comp2
            if ("POI" %in% labels(comp2)) {
              poi_sex_val <- getSex(comp2, ids = "POI")
            }
            
            # nComp: number of typed individuals in _comp2 excluding the missing person (POI)
            typed_ids <- typedMembers(comp2)
            typed_ids <- setdiff(typed_ids, "POI")
            ncomp_val <- length(typed_ids)
          }
          
          mpi_mm_map <- values$mpi_mendelian_mismatch
          poic_mm_map <- values$poic_mendelian_mismatch
          fam1c <- as.character(fam1)[1]
          fam2c <- as.character(fam2)[1]
          mpi_mismatch_val <- if (!is.null(mpi_mm_map) && length(mpi_mm_map) &&
              fam1c %in% names(mpi_mm_map)) {
            as.integer(mpi_mm_map[[fam1c]])
          } else {
            NA_integer_
          }
          poic_mismatch_val <- if (!is.null(poic_mm_map) && length(poic_mm_map) &&
              fam2c %in% names(poic_mm_map)) {
            as.integer(poic_mm_map[[fam2c]])
          } else {
            NA_integer_
          }
          
          results_list[[current_comparison]] <- tibble(
            MPI = fam1,
            MPI_mismatch = mpi_mismatch_val,
            POIc = fam2,
            POIc_mismatch = poic_mismatch_val,
            LR = as.numeric(lr),
            POI_sex = sex_label(poi_sex_val),
            nComp = ncomp_val,
            nMarkers = nMarkers_compared
          )

          # Drop large objects from this iteration
          rm(lr_obj, lr, lr_detail, lr_per_marker, nMarkers_uninformative, nMarkers_compared)
          if (current_comparison %% 20 == 0) {
            gc()
          }
          
          # Optionally update progress every 10 comparisons
          if (current_comparison %% 10 == 0 || current_comparison == total_comparisons) {
            # Already updated above via progress$set()
          }
        }, error = function(e) {
          msg <- sprintf("%s: %s", key, conditionMessage(e))
          lr_compute_errors <<- c(lr_compute_errors, msg)
          log_info(paste("LR computation failed:", msg))
        })
      }
      
      if (length(lr_compute_errors) > 0) {
        showNotification(
          paste(
            "LR computation failed for", length(lr_compute_errors),
            "comparison(s). Details were written to the R console."
          ),
          type = "warning",
          duration = 10
        )
        for (m in lr_compute_errors) message("LR computation: ", m)
      }
      
      if (values$stop_calculation) {
        progress$close()
        values$analysis_running <- FALSE
        return()
      }
      
      # Step 3: post-process and filter results
      log_info("Post-processing and filtering results")
      progress$set(message = "Processing results...", value = 0.9, detail = "Step 2/2")
      
      # Combine per-comparison results into a single data frame (omit failed comparisons)
      results_list <- purrr::compact(results_list)
      results <- dplyr::bind_rows(results_list)
      
      if (is.null(results) || nrow(results) == 0) {
        progress$close()
        values$analysis_running <- FALSE
        log_info("No results produced")
        showNotification("No results. Check data and settings.", 
                         type = "warning", duration = 5)
        return()
      }
      
      log_info(sprintf("Results: %d comparisons (stored unfiltered; table uses LR threshold)", nrow(results)))
      
      # Informative snapshot at current threshold (table can be refined without re-running)
      threshold <- ifelse(is.null(input$lr_threshold) || is.na(input$lr_threshold), 
                          1, input$lr_threshold)
      
      log_info("LR threshold snapshot (for notifications)")
      filter_info <- filter_results_by_threshold(results, threshold)
      
      if (filter_info$n_after == 0) {
        log_info("No rows meet the current LR threshold")
        showNotification(paste("No pairs with LR >=", threshold, "at current threshold. Lower the threshold or inspect full results in the console export."), 
                         type = "warning", duration = 6)
      } else if (filter_info$n_filtered > 0) {
        log_info(sprintf("At threshold %s: %d of %d rows listed in table", 
                         threshold, filter_info$n_after, filter_info$n_before))
        showNotification(paste("At LR >=", threshold, ":", filter_info$n_after, "of", filter_info$n_before, "pairs would appear in the table (adjust threshold anytime)."), 
                         type = "message", duration = 4)
      }
      
      values$results <- results
      values$lr_details_by_key <- lr_details_by_key
      values$analysis_complete <- TRUE
      values$analysis_running <- FALSE
      
      log_info("Analysis finished")
      
      progress$close()
      
      showNotification("Analysis complete.", type = "message", duration = 3)
      updateTabsetPanel(session, "main_tabs", selected = "comparisons")
      
    }, error = function(e) {
      progress$close()
      values$analysis_running <- FALSE
      showNotification(
        paste("Analysis error:", conditionMessage(e)),
        type = "error",
        duration = 10
      )
    })
  })
  
  # Hypothesis list entry for the active detail modal (MPI+POIc key from the selected results row only)
  selected_hypothesis_elem <- reactive({
    req(values$mpi_poic_list)
    key <- values$pedigree_plot_key
    req(!is.null(key), nzchar(as.character(key)[1]))
    key <- as.character(key)[1]
    if (key %in% names(values$mpi_poic_list)) {
      return(values$mpi_poic_list[[key]])
    }
    NULL
  })
  
  # ============================================
  # FILE INPUTS (dynamic UI so reset can fully clear uploads)
  # ============================================
  
  output$mpi_file_input <- renderUI({
    rid <- values$file_input_revision
    fileInput(
      paste0("file_mpi_", rid),
      label = NULL,
      accept = ".fam",
      width = "100%",
      buttonLabel = tags$span(icon("folder-open"), " Browse"),
      placeholder = "No file selected"
    )
  })
  
  output$poic_file_input <- renderUI({
    rid <- values$file_input_revision
    fileInput(
      paste0("file_poic_", rid),
      label = NULL,
      accept = ".fam",
      width = "100%",
      buttonLabel = tags$span(icon("folder-open"), " Browse"),
      placeholder = "No file selected"
    )
  })
  
  # ============================================
  # OUTPUTS - FILE INFORMATION
  # ============================================
  
  output$mpi_status <- renderUI({
    rid <- values$file_input_revision
    mpi_up <- input[[paste0("file_mpi_", rid)]]
    err <- values$mpi_load_error
    if (!is.null(err) && nzchar(as.character(err)[1])) {
      return(tagList(
        div(class = "status-badge status-error",
            icon("exclamation-circle"), " Error loading MPI file"),
        div(class = "help-text", style = "color: #721c24; margin-top: 5px;",
            tags$strong("Details:"), br(),
            tags$code(as.character(err)[1], style = "font-size: 11px; word-wrap: break-word;"))
      ))
    }
    if (!is.null(values$mpi_data)) {
      return(div(class = "status-badge status-success",
                 icon("check-circle"), " MPI file loaded"))
    }
    if (!is.null(mpi_up)) {
      return(div(class = "status-badge status-info",
                 icon("info-circle"), " MPI file selected"))
    }
    NULL
  })
  
  output$mpi_info <- renderUI({
    rid <- values$file_input_revision
    mpi_up <- input[[paste0("file_mpi_", rid)]]
    if (!is.null(values$mpi_data)) {
      div(class = "help-text",
          tags$strong("File:"), mpi_up$name %||% "", br(),
          tags$strong("Pedigrees:"), length(values$mpi_data)
      )
    } else if (!is.null(mpi_up)) {
      div(class = "help-text",
          tags$strong("File:"), mpi_up$name
      )
    } else {
      NULL
    }
  })
  
  output$poic_status <- renderUI({
    rid <- values$file_input_revision
    poic_up <- input[[paste0("file_poic_", rid)]]
    err <- values$poic_load_error
    if (!is.null(err) && nzchar(as.character(err)[1])) {
      return(tagList(
        div(class = "status-badge status-error",
            icon("exclamation-circle"), " Error loading POI Component file"),
        div(class = "help-text", style = "color: #721c24; margin-top: 5px;",
            tags$strong("Details:"), br(),
            tags$code(as.character(err)[1], style = "font-size: 11px; word-wrap: break-word;"))
      ))
    }
    if (!is.null(values$poic_data)) {
      return(div(class = "status-badge status-success",
                 icon("check-circle"), " POI Component file loaded"))
    }
    if (!is.null(poic_up)) {
      return(div(class = "status-badge status-info",
                 icon("info-circle"), " POI Component file selected"))
    }
    NULL
  })
  
  output$poic_info <- renderUI({
    rid <- values$file_input_revision
    poic_up <- input[[paste0("file_poic_", rid)]]
    if (!is.null(values$poic_data)) {
      div(class = "help-text",
          tags$strong("File:"), poic_up$name %||% "", br(),
          tags$strong("Pedigrees:"), length(values$poic_data)
      )
    } else if (!is.null(poic_up)) {
      div(class = "help-text",
          tags$strong("File:"), poic_up$name
      )
    } else {
      NULL
    }
  })
  
  # ============================================
  # OUTPUTS - DATA SUMMARY
  # ============================================
  
  output$data_loaded <- reactive({
    values$data_loaded
  })
  outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)
  
  # Status icons for side-panel sections
  output$config_status_icon <- renderUI({
    if (values$data_loaded) {
      tags$span(icon("check-circle", style = "color: #28a745; margin-left: 10px;"), 
                title = "Files loaded — settings enabled")
    } else {
      tags$span(icon("lock", style = "color: #ccc; margin-left: 10px;"), 
                title = "Load files to enable")
    }
  })
  
  output$comparison_status_icon <- renderUI({
    if (values$data_loaded && !is.null(input$mut_model) && input$mut_model != "") {
      tags$span(icon("check-circle", style = "color: #28a745; margin-left: 10px;"), 
                title = "Settings complete — comparison mode enabled")
    } else {
      tags$span(icon("lock", style = "color: #ccc; margin-left: 10px;"), 
                title = "Complete settings to enable")
    }
  })
  
  # Enable/disable side-panel inputs according to current state
  observe({
    if (values$data_loaded) {
      shinyjs::enable("mut_model")
      shinyjs::enable("mut_rate")
      shinyjs::enable("mut_range")
      shinyjs::enable("mut_range2")
      shinyjs::enable("mp_id")
    } else {
      shinyjs::disable("mut_model")
      shinyjs::disable("mut_rate")
      shinyjs::disable("mut_range")
      shinyjs::disable("mut_range2")
      shinyjs::disable("mp_id")
    }
  })
  
  observe({
    if (values$data_loaded && !is.null(input$mut_model) && input$mut_model != "") {
      shinyjs::enable("comparison_mode")
      shinyjs::enable("mp_sex_unknown")
    } else {
      shinyjs::disable("comparison_mode")
      shinyjs::disable("mp_sex_unknown")
    }
  })
  
  # ============================================
  # OUTPUTS - ANALYSIS STATUS
  # ============================================
  
  output$analysis_running <- reactive({
    values$analysis_running
  })
  outputOptions(output, "analysis_running", suspendWhenHidden = FALSE)
  
  output$analysis_complete <- reactive({
    values$analysis_complete
  })
  outputOptions(output, "analysis_complete", suspendWhenHidden = FALSE)
  
  # ============================================
  # OUTPUTS - FILTERED RESULTS
  # ============================================
  
  filtered_results <- reactive({
    if (is.null(values$results) || nrow(values$results) == 0) {
      return(empty_results_tibble())
    }
    
    filtered <- values$results
    
    threshold <- ifelse(is.null(input$lr_threshold) || is.na(input$lr_threshold),
                        1, input$lr_threshold)
    
    if ("LR" %in% names(filtered)) {
      filtered <- filtered %>%
        dplyr::filter(LR >= threshold | is.infinite(LR))
    }
    
    if ("GF" %in% names(filtered)) {
      names(filtered)[names(filtered) == "GF"] <- "MPI"
    }
    
    filtered
  })
  
  output$filtered_results_table <- DT::renderDataTable({
    if (!isTRUE(values$analysis_complete)) {
      fr0 <- empty_results_tibble()
      return(
        DT::datatable(
          fr0,
          options = list(dom = "t", ordering = FALSE),
          rownames = FALSE,
          selection = "none",
          class = "display compact cell-border"
        ) %>%
          DT::formatStyle(columns = colnames(fr0), textAlign = "center", verticalAlign = "middle")
      )
    }
    
    fr <- filtered_results()
    if (nrow(fr) == 0) {
      return(
        DT::datatable(
          fr,
          options = list(
            scrollX = TRUE,
            dom = "lBfrtip",
            buttons = c("copy", "csv", "excel"),
            columnDefs = list(list(className = "dt-center", targets = "_all"))
          ),
          extensions = "Buttons",
          rownames = FALSE,
          filter = "top",
          selection = "single",
          class = "display compact cell-border"
        ) %>%
          DT::formatStyle(columns = colnames(fr), textAlign = "center", verticalAlign = "middle")
      )
    }
    
    dt <- DT::datatable(
      fr,
      options = list(
        pageLength = 25,
        lengthMenu = list(
          c(10, 25, 50, 100, 250, 500, -1),
          c("10", "25", "50", "100", "250", "500", "All")
        ),
        scrollX = TRUE,
        order = list(list(4, "desc")),
        dom = "lBfrtip",
        buttons = c("copy", "csv", "excel"),
        columnDefs = list(list(className = "dt-center", targets = "_all"))
      ),
      extensions = "Buttons",
      rownames = FALSE,
      filter = "top",
      selection = "single",
      class = "display compact cell-border"
    )
    if ("LR" %in% names(fr)) {
      dt <- dt %>%
        DT::formatSignif(columns = "LR", digits = 5)
    }
    dt <- dt %>%
      DT::formatStyle(columns = colnames(fr), textAlign = "center", verticalAlign = "middle")
    
    if ("LR" %in% names(fr)) {
      dt <- dt %>%
        DT::formatStyle(
          columns = colnames(fr),
          valueColumns = "LR",
          backgroundColor = DT::styleInterval(0.999999999, c("#ffebee", "#e8f5e9"))
        )
    }
    
    dt
  })

  # ============================================
  # OUTPUT — PER-MARKER LR TABLE (DETAIL MODAL)
  # Row order matches stored partial LRs; alleles from getAlleles(H0), columns *.1 / *.2
  # ============================================
  output$lr_details_modal_table <- DT::renderDataTable({
    req(values$analysis_complete)
    req(values$lr_modal_key)
    req(values$mpi_poic_list)
    req(values$last_lr_modal_mutation)
    req(values$lr_modal_total)
    
    key <- values$lr_modal_key
    stored <- values$lr_details_by_key[[key]]
    validate(need(
      !is.null(stored) && length(stored$Marker) > 0,
      "No per-marker LRs for this pair. Run «Compute LR» again."
    ))
    
    validate(need(!is.null(values$mpi_poic_list[[key]]), "Hypothesis not found for this MPI × POI Component pair."))
    
    detail_df <- build_lr_modal_detail_df(
      key = key,
      mpi_poic_list = values$mpi_poic_list,
      lr_details_by_key = values$lr_details_by_key,
      lr_modal_total = values$lr_modal_total,
      mut_settings = values$last_lr_modal_mutation,
      default_mp_id = CONFIG$mp_id
    )
    validate(need(!is.null(detail_df), "Invalid hypothesis: missing 'Ped 1' (H0) or could not build table."))
    
    sample_side <- attr(detail_df, "sample_column_side")
    
    # 0-based LR column index for JS (column always exists; fallback = second column)
    lr_col_idx_js <- match("LR", names(detail_df))
    if (is.na(lr_col_idx_js)) lr_col_idx_js <- 2L
    lr_col_idx_js <- as.integer(lr_col_idx_js - 1L)
    
    dt_obj <- DT::datatable(
      detail_df,
      options = list(
        dom = "t",
        paging = FALSE,
        scrollX = FALSE,
        autoWidth = FALSE,
        ordering = FALSE,
        columnDefs = list(list(className = "dt-center", targets = "_all")),
        # Equal column widths + row classes: LR=0 (.lr-partial-zero), partial LR=1 (.lr-partial-one, not Total row).
        # Done in drawCallback (not rowCallback) to avoid invalid JS if the index fails.
        drawCallback = htmlwidgets::JS(paste0(
          "function(settings) {",
          "  var lrColIdx = ", lr_col_idx_js, ";",
          "  var api = new $.fn.dataTable.Api(settings);",
          "  function markPartialLRRowClasses() {",
          "    $(api.table().body()).find('tr').each(function() {",
          "      var $tr = $(this);",
          "      var $tds = $tr.find('td');",
          "      if ($tds.length <= lrColIdx) return;",
          "      var markerTxt = $tds.eq(0).text().replace(/\\s+/g, ' ').trim();",
          "      var isTotalRow = markerTxt === 'Total LR';",
          "      var txt = $tds.eq(lrColIdx).text().replace(/\\s/g, '');",
          "      var v = parseFloat(txt.replace(',', '.'));",
          "      var isZero = !isNaN(v) && v === 0;",
          "      var isOne = !isNaN(v) && Math.abs(v - 1) < 1e-8;",
          "      $tr.toggleClass('lr-partial-zero', isZero);",
          "      $tr.toggleClass('lr-partial-one', !isTotalRow && isOne);",
          "    });",
          "  }",
          "  function runEqual() {",
          "    var n = api.columns().count();",
          "    if (n < 1) return;",
          "    var pct = (100 / n) + '%';",
          "    var css = { width: pct, minWidth: 0, maxWidth: 'none', boxSizing: 'border-box' };",
          "    $(api.table().node()).css({ width: '100%', tableLayout: 'fixed' });",
          "    for (var i = 0; i < n; i++) {",
          "      $(api.column(i).header()).css(css);",
          "      $(api.column(i).nodes()).css(css);",
          "    }",
          "  }",
          "  runEqual();",
          "  markPartialLRRowClasses();",
          "  [0, 50, 150, 300, 500].forEach(function(ms) {",
          "    setTimeout(function() {",
          "      if (window.kinshipEqualizeDetailModalTable) window.kinshipEqualizeDetailModalTable(null);",
          "      else runEqual();",
          "      markPartialLRRowClasses();",
          "    }, ms);",
          "  });",
          "}"
        ))
      ),
      rownames = FALSE
    ) %>%
      DT::formatStyle(
        "Marker",
        fontWeight = "bold",
        textAlign = "center",
        verticalAlign = "middle"
      ) %>%
      DT::formatSignif(columns = "LR", digits = 5) %>%
      DT::formatStyle(
        setdiff(names(detail_df), "LR"),
        textAlign = "center",
        verticalAlign = "middle"
      ) %>%
      DT::formatStyle(
        "LR",
        textAlign = "center",
        verticalAlign = "middle"
      )
    
    # Subtle column backgrounds: MPI (reference family) vs POI Component genotypes
    if (!is.null(sample_side) && length(sample_side) > 0) {
      poic_cols <- names(sample_side)[sample_side == "poic"]
      mpi_cols <- names(sample_side)[sample_side == "mpi"]
      if (length(poic_cols) > 0) {
        dt_obj <- dt_obj %>% DT::formatStyle(
          poic_cols,
          backgroundColor = "#eef6fc",
          textAlign = "center",
          verticalAlign = "middle"
        )
      }
      if (length(mpi_cols) > 0) {
        dt_obj <- dt_obj %>% DT::formatStyle(
          mpi_cols,
          backgroundColor = "#f5f4f0",
          textAlign = "center",
          verticalAlign = "middle"
        )
      }
    }
    
    return(dt_obj)
  })
  
  outputOptions(output, "lr_details_modal_table", suspendWhenHidden = FALSE)
  
  # Pedigree plot inside the results-row modal only (600px height matches the former main-panel plot)
  output$hypothesis_plot_modal <- renderPlot({
    req(selected_hypothesis_elem())
    hyp <- selected_hypothesis_elem()
    validate(need(is.list(hyp), "No valid hypothesis structure available for plotting."))
    
    mp_id_val <- if (!is.null(input$mp_id) && input$mp_id != "") {
      input$mp_id
    } else {
      CONFIG$mp_id
    }
    
    tryCatch({
      missing_branch_plot(hypothesis_elem = hyp,
                        missing_id = mp_id_val)
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error plotting pedigree:", conditionMessage(e)), 
           cex = 1.2, col = "red")
    })
  }, width = reactive({
    w <- session$clientData$output_hypothesis_plot_modal_width
    if (is.null(w) || !is.finite(suppressWarnings(as.numeric(w)))) {
      return(1400L)
    }
    nw <- as.integer(floor(as.numeric(w)))
    if (length(nw) != 1L || nw < 320) 1400L else nw
  }), height = 600)
  
  outputOptions(output, "hypothesis_plot_modal", suspendWhenHidden = FALSE)
  
  # Open detail modal: pedigree (top), downloads, per-marker LR table (bottom); pair is taken only from the row
  observeEvent(input$filtered_results_table_rows_selected, {
    req(input$filtered_results_table_rows_selected)
    req(filtered_results())
    req(values$analysis_complete)
    
    comparison <- filtered_results()[input$filtered_results_table_rows_selected, , drop = FALSE]
    mpi_name <- if ("MPI" %in% names(comparison)) {
      as.character(comparison$MPI[1])
    } else if ("GF" %in% names(comparison)) {
      as.character(comparison$GF[1])
    } else {
      NULL
    }
    
    req(!is.null(comparison$POIc))
    req(!is.null(mpi_name), nzchar(mpi_name))
    
    poic_name <- as.character(comparison$POIc)[1]
    
    modal_key <- paste0(mpi_name, "+", poic_name)
    values$pedigree_plot_key <- modal_key
    values$lr_modal_key <- modal_key
    values$lr_modal_total <- suppressWarnings(as.numeric(comparison$LR))
    
    removeModal()
    showModal(modalDialog(
      title = tags$div(
        tags$div(style = "font-weight: 600;", "Comparison detail"),
        tags$div(tags$small(paste0(mpi_name, " vs ", poic_name)))
      ),
      tagList(
        tags$div(
          class = "kinship-modal-detail-inner",
          tags$h5(style = "margin-top: 0;", icon("project-diagram"), " Pedigree hypothesis plot"),
          tags$div(
            class = "kinship-modal-pedigree-block",
            plotOutput("hypothesis_plot_modal", height = "600px", width = "100%")
          ),
          fluidRow(
            column(
              6,
              downloadButton(
                "download_pedigree_modal",
                label = tags$span(icon("download"), " Download pedigree image"),
                class = "btn-primary",
                style = "width: 100%;"
              )
            ),
            column(
              6,
              downloadButton(
                "download_lr_modal_table",
                label = tags$span(icon("download"), " Download table (CSV)"),
                class = "btn-primary",
                style = "width: 100%;"
              )
            )
          ),
          tags$hr(style = "margin: 16px 0; border-top: 1px solid #e0e0e0;"),
          tags$h5(style = "margin-top: 0;", icon("table"), " Genotypes and likelihood ratios"),
          DT::dataTableOutput("lr_details_modal_table", width = "100%")
        )
      ),
      footer = modalButton("Close"),
      easyClose = TRUE,
      fade = TRUE,
      size = "l",
      # kinship-detail-modal: fixed width in style.css + column refit script in global.R
      class = "modal-fit-content kinship-detail-modal"
    ))
  })
  
  output$download_lr_modal_table <- downloadHandler(
    filename = function() {
      key <- values$lr_modal_key %||% "LR_detail"
      safe <- gsub("[^A-Za-z0-9._+-]", "_", as.character(key)[1])
      paste0("LR_per_marker_", safe, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(values$lr_modal_key)
      req(values$mpi_poic_list)
      req(values$lr_details_by_key)
      req(values$last_lr_modal_mutation)
      df <- build_lr_modal_detail_df(
        key = values$lr_modal_key,
        mpi_poic_list = values$mpi_poic_list,
        lr_details_by_key = values$lr_details_by_key,
        lr_modal_total = values$lr_modal_total,
        mut_settings = values$last_lr_modal_mutation,
        default_mp_id = CONFIG$mp_id
      )
      if (is.null(df)) {
        showNotification(
          "No data to export. Run Compute LR and open the modal from a results row.",
          type = "error"
        )
        stop("No data to export.")
      }
      if ("LR" %in% names(df)) {
        df$LR <- signif(as.numeric(df$LR), digits = 5)
      }
      utils::write.csv(df, file, row.names = FALSE, fileEncoding = "UTF-8")
      notify_after_download("Table downloaded successfully.")
    }
  )
  
  # Download pedigree image from the detail modal (same PNG resolution as before: 1600x1000 @ 150 dpi)
  output$download_pedigree_modal <- downloadHandler(
    filename = function() {
      key <- values$pedigree_plot_key %||% "MPI_POIc"
      safe <- gsub("[^A-Za-z0-9._+-]", "_", as.character(key)[1])
      paste0("pedigree_", safe, "_", format(Sys.Date(), "%Y%m%d"), ".png")
    },
    content = function(file) {
      req(selected_hypothesis_elem())
      
      tryCatch({
        png(file, width = 1600, height = 1000, res = 150)
        
        mp_id_val <- if (!is.null(input$mp_id) && input$mp_id != "") {
          input$mp_id
        } else {
          CONFIG$mp_id
        }
        
        missing_branch_plot(hypothesis_elem = selected_hypothesis_elem(),
                          missing_id = mp_id_val)
        dev.off()
        notify_after_download("Pedigree image downloaded successfully.")
      }, error = function(e) {
        showNotification(paste("Download error:", e$message), 
                         type = "error", duration = 5)
      })
    }
  )
  
  # ============================================
  # DOWNLOAD HANDLERS
  # ============================================
  
  # Download handler for the main CSV report
  output$downloadReport <- downloadHandler(
    filename = function() {
      format(Sys.time(), "KinshipAssembly_report_%Y-%m-%d_%H%M%S.csv")
    },
    content = function(file) {
      req(values$results)
      report_df <- isolate(values$results)
      if (is.null(report_df) || nrow(report_df) == 0L) {
        showNotification("No results to export. Run Compute LR first.", type = "warning")
        stop("No results to export.", call. = FALSE)
      }
      ss <- add_summary_stats(report_df)
      tryCatch({
        export_to_csv(
          result_df = report_df,
          file_path = file,
          summary_stats = ss
        )
        showNotification("Report downloaded.", type = "message")
      }, error = function(e) {
        showNotification(
          paste("Report download error:", conditionMessage(e)),
          type = "error"
        )
      })
    }
  )
  
  # ============================================
  # RESET BUTTON
  # ============================================
  
  # Reset button
  observeEvent(input$reset, {
    # If an analysis is running, request it to stop first
    if (values$analysis_running) {
      values$stop_calculation <- TRUE
      values$analysis_running <- FALSE
    }
    
    # Ask user to confirm full reset
    showModal(modalDialog(
      title = tags$span(icon("exclamation-triangle"), " Confirm reset"),
      tags$div(
        class = "kinship-reset-confirm-marker",
        "This will clear all loaded data and results. Continue?"
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_reset", "Reset", class = "btn-danger")
      ),
      easyClose = FALSE,
      fade = TRUE,
      class = "kinship-reset-confirm-modal"
    ))
  }, ignoreInit = TRUE)
  
  observeEvent(input$confirm_reset, {
    # Close confirm dialog and any open Shiny modal (e.g. comparison detail underneath)
    removeModal()
    removeModal()
    
    if (isTRUE(values$analysis_running)) {
      values$stop_calculation <- TRUE
      values$analysis_running <- FALSE
    }
    
    # Clear modal / LR detail state (was missing and left stale UI)
    values$lr_modal_key <- NULL
    values$lr_modal_total <- NULL
    values$last_lr_modal_mutation <- NULL
    values$lr_details_by_key <- NULL
    values$pedigree_plot_key <- NULL
    
    # Clear data and results
    values$results <- NULL
    values$mpi_poic_list <- NULL
    values$locus_attributes <- NULL
    values$poic_data <- NULL
    values$mpi_data <- NULL
    values$mpi_mendelian_mismatch <- NULL
    values$poic_mendelian_mismatch <- NULL
    values$data_loaded <- FALSE
    values$locus_loaded <- FALSE
    values$analysis_complete <- FALSE
    values$mpi_load_error <- NULL
    values$poic_load_error <- NULL
    values$stop_calculation <- FALSE
    
    # Recreate file inputs so the browser truly clears filenames (resetInput/shinyjs alone is unreliable)
    values$file_input_revision <- values$file_input_revision + 1L
    
    # Reset settings / comparison (not the whole sidebar — avoids fighting dynamic file UI)
    shinyjs::reset("config_section")
    shinyjs::reset("comparison_section")
    
    # Align with CONFIG
    updateSelectInput(session, "mut_model", selected = CONFIG$mut_model)
    updateNumericInput(session, "mut_rate", value = CONFIG$mut_rate)
    updateNumericInput(session, "mut_range", value = 0.1)
    updateNumericInput(session, "mut_range2", value = 0.000001)
    updateTextInput(session, "mp_id", value = CONFIG$mp_id)
    updateCheckboxGroupButtons(session, "mp_sex_unknown", selected = character(0))
    updateNumericInput(session, "lr_threshold", value = 1)
    updateRadioGroupButtons(session, "comparison_mode", selected = "all")
    
    shinyjs::disable("run_analysis")
    shinyjs::disable("stop_analysis")
    
    updateTabsetPanel(session, "main_tabs", selected = "comparisons")
    
    output$hypothesis_plot_modal <- renderPlot(NULL)
    
    # Clear DT row selection after flush (avoids ghost selection on new data)
    session$onFlushed(function() {
      proxy <- DT::dataTableProxy("filtered_results_table")
      DT::selectRows(proxy, NULL)
    }, once = TRUE)
    
    showNotification("Session reset. You can load new files.", type = "message", duration = 3)
  })
}


