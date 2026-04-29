## ============================================================
## UI - USER INTERFACE
## Comparisons: results table only; row click opens modal (pedigree + per-marker LRs + downloads).
## ============================================================

ui <- fluidPage(
  custom_css,
  useShinyjs(),
  
  div(class = "main-header",
      div(class = "main-header-content",
          h1(icon("dna", class = "fa-lg"), " ", UI_CONFIG$app_title),
          h3(UI_CONFIG$app_subtitle)
      )
  ),
  
  fluidRow(
    column(width = 3,
           div(id = "kinship_sidebar", class = "sidebar-panel",
               div(class = "section-title", icon("upload"), " Upload files"),
               
               div(class = "upload-file-block",
                   div(class = "upload-header-row",
                       tags$span(class = "upload-field-label",
                                icon("file"), " MPI/DVI (.fam)"),
                       div(class = "upload-status-inline",
                           uiOutput("mpi_status"))),
                   uiOutput("mpi_file_input")),
               
               uiOutput("mpi_info"),
               
               br(),
               
               div(class = "upload-file-block",
                   div(class = "upload-header-row",
                       tags$span(class = "upload-field-label",
                                 icon("file"), " POI Component (.fam)"),
                       div(class = "upload-status-inline",
                           uiOutput("poic_status"))),
                   uiOutput("poic_file_input")),
               
               uiOutput("poic_info"),
               
               tags$hr(style = "border-top: 1px solid #e0e0e0; margin: 20px 0;"),
               
               div(id = "config_section",
                 div(class = "section-title", icon("cog"), " Settings",
                     uiOutput("config_status_icon")),
                 
                 selectInput("mut_model",
                            label = tags$span(icon("sliders-h"), " Mutation model", tags$span("*", style = "color: red;")),
                            choices = list(
                              "None" = "none",
                              "Equal" = "equal",
                              "Extended stepwise" = "stepwise"
                            ),
                            selected = CONFIG$mut_model,
                            width = "100%"),
                 
                 conditionalPanel(
                   condition = "input.mut_model == 'equal' || input.mut_model == 'stepwise'",
                   numericInput("mut_rate",
                               label = tags$span(icon("percent"), " Mutation rate"),
                               value = 0.002,
                               min = 0,
                               max = 1,
                               step = 0.0001,
                               width = "100%"),
                   tags$p(class = "help-text", 
                          "Default: 0.002")
                 ),
                 
                 conditionalPanel(
                   condition = "input.mut_model == 'stepwise'",
                   numericInput("mut_range",
                               label = tags$span(icon("arrows-alt-h"), " Range"),
                               value = 0.1,
                               min = 0,
                               max = 10,
                               step = 0.1,
                               width = "100%"),
                  tags$p(class = "help-text",
                         "Default: 0.1"),
                   numericInput("mut_range2",
                               label = tags$span(icon("arrows-alt-h"), " Range 2"),
                               value = 0.000001,
                               min = 0,
                               max = 10,
                               step = 0.1,
                              width = "100%"),
                  tags$p(class = "help-text",
                         "Default: 0.000001")
                 ),
                 
                 br(),
                 
                 textInput("mp_id",
                          label = tags$span(icon("user"), " Missing person ID"),
                          value = CONFIG$mp_id,
                          width = "100%"),
                 
                 tags$p(class = "help-text", 
                        "Label used to identify the missing person in pedigrees."),
                 
                 tags$hr(style = "border-top: 1px solid #e0e0e0; margin: 20px 0;")
               ),
               
               div(id = "comparison_section",
                 div(class = "section-title", icon("exchange-alt"), " Comparison mode",
                     uiOutput("comparison_status_icon")),
                 
                 radioGroupButtons(
                   inputId = "comparison_mode",
                  label = NULL,
                  choices = c("All vs all" = "all", 
                              "Choose comparisons" = "custom"),
                   selected = "all",
                   checkIcon = list(
                     yes = icon("check-circle", style = "color: #28a745"),
                     no = icon("circle", style = "color: #ccc")
                   ),
                   justified = TRUE,
                   status = "primary"
                 ),
                 
                 tags$p(class = "help-text", 
                       "Compare every MPI/DVI with every POI Component, or restrict to selected pairs."),
                 
                 uiOutput("comparison_select_ui"),
                 
                 tags$hr(style = "border-top: 1px solid #e0e0e0; margin: 20px 0;")
               ),
               
               numericInput("lr_threshold", 
                          label = tags$span(icon("filter"), " LR threshold"),
                          value = 1,
                          min = 0,
                          step = 0.1,
                          width = "100%"),
               
               tags$p(class = "help-text", 
                      "After «Compute LR», the table shows pairs with LR \u2265 threshold; you can change the threshold without re-running."),
               
               tags$hr(style = "border-top: 1px solid #e0e0e0; margin: 20px 0;"),
               
               div(class = "section-title", icon("play-circle"), " Actions"),
               
               actionButton("run_analysis", 
                           label = tags$span(icon("calculator"), " Compute LR"),
                           class = "btn-success btn-action"),
               
               br(),
               
               actionButton("stop_analysis", 
                           label = tags$span(icon("stop"), " Stop"),
                           class = "btn-danger btn-action"),
               
               br(),
               
               uiOutput("download_report_ui"),
               
               br(),
               
               actionButton("reset", 
                           label = tags$span(icon("redo"), " Reset"),
                           class = "btn-secondary btn-action"),
               
               tags$hr(style = "border-top: 1px solid #e0e0e0; margin: 20px 0;")
           )
    ),
    
    column(width = 9,
           tabsetPanel(
             id = "main_tabs",
             type = "tabs",
             
             tabPanel(
               title = tags$span(icon("chart-line"), " Comparisons"),
               value = "comparisons",
               
               br(),
               
               div(class = "results-only-panel",
                   h4(icon("table"), " Results"),
                   tags$p(class = "help-text",
                         "After ", tags$strong("Compute LR"), ", the table lists each MPI/DVI\u2013POI Component pair whose ",
                          tags$strong("LR"), " is at or above your ", tags$strong("LR threshold"), " (default 1). ",
                          "Change the threshold in the sidebar without running again. ",
                          "Click a row for pedigrees, genotypes, and per-marker LRs."),
                   conditionalPanel(
                     condition = "output.analysis_complete",
                     tagList(
                       tags$p(class = "help-text", style = "font-style: normal; color: #495057; margin-bottom: 6px;",
                              tags$strong("Columns (brief):")),
                       tags$ul(class = "results-column-key",
                              tags$li(tags$strong("MPI/DVI"), " \u2014 Name of the MPI/DVI family."),
                              tags$li(tags$strong("MPI_mismatch"), " \u2014 Mendelian inconsistencies inside the MPI/DVI."),
                               tags$li(tags$strong("POIc"), " \u2014 Name of the POI Component family being compared."),
                               tags$li(tags$strong("LR"), " \u2014 Total likelihood ratio for this pair (support for \u201csame person\u201d vs unrelated)."),
                              tags$li(tags$strong("POI_sex"), " \u2014 Sex recorded for the POI in this comparison (M / F / UNK)."),
                               tags$li(tags$strong("nComp"), " \u2014 Number of typed relatives in the POI Component branch (the POI label excluded)."),
                               tags$li(tags$strong("nMarkers"), " \u2014 Markers that entered the LR after harmonisation.")
                       ),
                       DT::dataTableOutput("filtered_results_table", width = "100%")
                     )
                   ),
                   conditionalPanel(
                     condition = "!output.analysis_complete",
                     div(class = "info-card",
                         h5(icon("info-circle"), " No results yet"),
                       tags$p("Load MPI/DVI and POI Component files from the sidebar, then run the analysis.")
                     )
                   )
               )
             ),
             
             tabPanel(
               title = tags$span(icon("question-circle"), " Help"),
               value = "help",
               
               br(),
                       
               div(class = "info-card",
                   h5(icon("lightbulb"), " In plain language"),
                   tags$p(
                    "This app compares ", tags$strong("MPI/DVI"), " (Missing Person Identification / Disaster Victim Identification) families with ",
                     tags$strong("POI Component"), " (Person Of Interest Component) families using STR markers. ",
                     "For each pair it computes a ", tags$strong("likelihood ratio (LR)"), ": ",
                     "how much the genetic data support the idea that the missing person is the POI ",
                     "versus being unrelated. Higher LR usually means stronger support; your ",
                     tags$strong("LR threshold"), " hides weaker pairs from the table."
                   )
               ),
               
               div(class = "info-card",
                  h5(icon("book"), " Quick start"),
                   tags$ol(
                     tags$li(
                      tags$strong("Upload"), " \u2014 Choose an MPI/DVI ", tags$code(".fam"), " file, then a POI Component ", tags$code(".fam"), ". ",
                      "MPI/DVI first helps the app build the comparison list correctly."
                     ),
                     tags$li(
                       tags$strong("Check settings"), " \u2014 Mutation model and ", tags$strong("Missing person ID"), " must match how your files are built ",
                       "(the ID is the label of the missing individual in the pedigrees)."
                     ),
                     tags$li(
                      tags$strong("Choose comparisons"), " \u2014 ", tags$strong("All vs all"), " runs every MPI/DVI with every POI Component family. ",
                       tags$strong("Choose comparisons"), " lets you pick specific names in the lists (after both files are loaded)."
                     ),
                     tags$li(
                       tags$strong("Run"), " \u2014 Click ", tags$strong("Compute LR"), ". ",
                       "Large jobs can take a while; use ", tags$strong("Stop"), " if needed."
                     ),
                     tags$li(
                       tags$strong("Read the table"), " \u2014 Use the ", tags$strong("LR threshold"), " to show only pairs with LR at or above that value (default ", tags$code("1"), "). ",
                       "You can change the threshold without running again."
                     ),
                     tags$li(
                       tags$strong("Open a row"), " \u2014 Click a row to see pedigrees, genotypes, and per-marker LRs, and to download an image or CSV."
                     ),
                     tags$li(
                       tags$strong("Export"), " \u2014 ", tags$strong("Download report"), " saves a CSV of the rows currently shown, plus a short summary block at the end."
                     ),
                     tags$li(
                       tags$strong("Reset"), " \u2014 ", tags$strong("Reset"), " clears uploads and results for a new session."
                     )
                   )
               ),
               
               div(class = "info-card",
                  h5(icon("exclamation-triangle"), " Please note"),
                  tags$ul(
                    tags$li("Files must be compatible ", tags$code(".fam"), " inputs for this workflow."),
                    tags$li("Alleles that do not match the panel are set to missing so the LR stays consistent."),
                    tags$li("If both alleles at a marker are missing, that marker does not help for that pair."),
                   tags$li("More markers and more MPI/DVI\u00d7POIc pairs mean longer run times."),
                    tags$li("The ", tags$strong("Missing person ID"), " must be spelled exactly as in the pedigrees.")
                  )
              ),
               
               div(class = "info-card",
                  h5(icon("cog"), " Technical notes"),
                   tags$ul(
                     tags$li(
                      tags$strong("Hypotheses:"), " ",
                      "H0: missing person is the POI; H1: unrelated. LR = P(data | H0) / P(data | H1)."
                     ),
                     tags$li(
                      tags$strong("Mutation:"), " ",
                      "Applied on the missing person\u2019s genotype before LR. ",
                      "None \u2192 rate 0; Equal \u2192 your rate; Extended stepwise \u2192 Rate, Range and Range 2."
                     ),
                     tags$li(
                      tags$strong("Plots:"), " ",
                      "H0 shows the merged pedigree; H1 splits MPI/DVI branch and POI Component branch."
                     ),
                     tags$li(
                      tags$strong("Markers:"), " ",
                      "Locus info follows the MPI/DVI file; MPI/DVI and POI Component markers are harmonised; invalid alleles \u2192 ", tags$code("0"), "."
                     ),
                     tags$li(
                       tags$strong("Single-side markers:"), " ",
                       "Markers typed on one family only may be excluded from the product LR but still listed with partial factor 1 where applicable."
                     )
                   )
               ),
               
               div(class = "info-card app-version-banner",
                   tags$p(
                     tags$strong("Version "), UI_CONFIG$version,
                     tags$strong(" \u00b7 Author: "), UI_CONFIG$author
                   )
               )
             )
           )
    )
  )
)
