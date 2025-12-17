forest_tab_ui <- function() {
  tabPanel("Risk by Gene Signature",
           fluidRow(
             column(5,
                    tags$div(class = "main-title", "Risk by Gene Signature"),
                    tags$div(style = "margin-bottom:12px; font-size:1.08em; color:#2c3e50;",
                             "This tab displays risk analysis by predefined gene signature using Cox proportional hazards models. Hazard ratios and confidence intervals are calculated for selected gene signatures and clinical endpoints, using the currently filtered sample set."
                    ),
                    tags$div(
                      style = "margin-bottom:8px;",
                      strong("How to construct a forest plot:"),
                      tags$ul(
                        tags$li("Select a clinical endpoint for survival analysis."),
                        tags$li("Filter samples by subtype as needed."),
                        tags$li("Choose one or more gene signatures to display."),
                        tags$li("Adjust the number of bins for numeric signatures to control grouping granularity."),
                        tags$li("The plot and results table will update automatically for selected samples.")
                      )
                    ),
                    tags$div(class = "section-title", "Risk by Gene Signature (Cox Model)"),
                    selectInput(
                      "forest_endpoint", "Select clinical endpoint",
                      choices = c(
                        "Biological Progression",
                        "Clinical Progression",
                        "Clinical Progression During Follow-up",
                        "BCG Any",
                        "BCG Adequate",
                        "Recurrence"
                      ),
                      selected = "Clinical Progression"
                    ),
                    pickerInput(
                      "forest_filter_subtypes",
                      "Filter samples by subtype",
                      choices = levels(metadata$subtype_5_class),
                      selected = levels(metadata$subtype_5_class),
                      multiple = TRUE,
                      options = list(`actions-box` = TRUE, `live-search` = TRUE)
                    ),
                    pickerInput(
                      "forest_signature_select",
                      "Select signature(s) to show",
                      choices = c("All", forest_signature_names),
                      selected = "All",
                      multiple = TRUE,
                      options = list(`actions-box` = TRUE, `live-search` = TRUE)
                    ),
                    sliderInput(
                      "forest_binsize",
                      "Number of bins for numeric signatures",
                      min = 2, max = 20, value = 10, step = 1
                    )
             ),
             column(7,
                    plotOutput("forest_plot", height = "700px")
             )
           ),
           fluidRow(
             column(12,
                    DT::dataTableOutput("forest_table"),
                    downloadButton("download_forest_table", "Download Table", style = "margin-top:10px;")
             )
           )
  )
}