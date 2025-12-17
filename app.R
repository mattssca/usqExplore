# Load global objects and utilities
source("global.R")

# Load UI modules
source("ui/ui_metadata_tab.R")
source("ui/ui_exploration_tab.R")
source("ui/ui_heatmap_tab.R")
source("ui/ui_survival_tab.R")
source("ui/ui_forest_tab.R")
source("ui/ui_correlation_tab.R")
source("ui/ui_external_tab.R")

# Load server modules  
source("server/server_metadata.R")
source("server/server_exploration.R")
source("server/server_heatmap.R")
source("server/server_survival.R")
source("server/server_forest.R")
source("server/server_correlation.R")

ui <- fluidPage(
  theme = bs_theme(bootswatch = "minty"),
  tags$head(
    tags$style(HTML("
      .annotation-legends-row {
        display: flex;
        flex-wrap: wrap;
        gap: 30px;
        align-items: flex-start;
      }
      .annotation-legend-item {
        min-width: 120px;
        margin-right: 20px;
      }
      .main-title { font-size: 2em; font-weight: bold; margin-bottom: 10px; }
      .plot-container { background: #fff; padding: 20px; border-radius: 8px; box-shadow: 0 2px 8px #ccc; }
      .filter-title { font-size: 1.3em; font-weight: bold; margin-bottom: 10px; }
      .after-plot-text { margin-top: 15px; font-size: 1.1em; }
      .section-title { font-size: 1.2em; font-weight: bold; margin-top: 15px; margin-bottom: 8px; }
      .btn-file, .fileinput-button {
        background-color: #dd6262 !important;
        border-color: #dd6262 !important;
        color: #fff !important;
      }
    "))
  ),
  fluidRow(
    column(
      width = 12,
      tags$img(src = "logo.png", height = "200px", style = "margin-bottom:10px;")
    )
  ),
  tabsetPanel(
    tabPanel("Welcome",
             fluidRow(
               column(12,
                      tags$div(style = "margin-top:20px; font-size:1.15em;",
                               tags$p("Welcome to the UROSCANSEQ Explorer! This app allows you to interactively explore gene expression and clinical metadata from the UROSCANSEQ dataset."),
                               tags$p("Start by filtering samples in the Metadata tab. All downstream analyses will use the selected samples."),
                               tags$ul(
                                 tags$li("Metadata tab: Filter samples and view metadata table."),
                                 tags$li("Exploration Plots: Create interactive plots based on selected samples."),
                                 tags$li("Heatmap: Visualize gene expression heatmaps."),
                                 tags$li("Gene Expression Outcomes: Kaplan-Meier survival analysis."),
                                 tags$li("Risk by Gene Signature: Cox model forest plots for gene signatures.")
                               ),
                               tags$p(
                                 "For more information and source code, visit our ",
                                 tags$a(href = "https://github.com/LundBladderCancerGroup", "GitHub", target = "_blank"),
                                 "."
                               )
                      )
               )
             )
    ),
    metadata_tab_ui(),
    exploration_tab_ui(),
    heatmap_tab_ui(),
    survival_tab_ui(),
    forest_tab_ui(),
    correlation_tab_ui(),
    external_tab_ui()
  )
)

server <- function(input, output, session) {
  # Call metadata server module
  metadata_module <- metadata_server(input, output, session)
  
  # Call other server modules
  exploration_server(input, output, session, metadata_module$data_filtered)
  heatmap_server(input, output, session, metadata_module$data_filtered)
  survival_server(input, output, session, metadata_module$data_filtered)
  forest_server(input, output, session, metadata_module$data_filtered)
  correlation_server(input, output, session, metadata_module$data_filtered)
}

shinyApp(ui, server)