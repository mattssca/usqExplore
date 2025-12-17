heatmap_tab_ui <- function() {
  tabPanel("Heatmap",
           tags$div(class = "main-title", "Heatmap"),
           tags$div(style = "margin-bottom:8px; font-size:1.08em; color:#2c3e50;",
                    "Visualize gene expression heatmaps for selected genes and samples. Add annotation tracks and meta-genes to explore sample-level information and gene signatures in the filtered dataset."
           ),
           tags$div(
             style = "margin-bottom:8px;",
             strong("How to construct a heatmap:"),
             tags$ul(
               tags$li("Upload a custom gene list or select a predefined signature."),
               tags$li("Optionally define meta-genes for annotation tracks."),
               tags$li("Add annotation variables to display sample-level information."),
               tags$li("Click 'Plot Heatmap' to generate the heatmap for selected samples and genes.")
             )
           ),
           tags$div(class = "section-title", "Genes"),
           fileInput("gene_file", "Upload gene list (txt, one gene per line)", accept = ".txt"),
           pickerInput(
             "gene_signature",
             "Select predefined gene list (signature)",
             choices = c("Custom Genes", unique(genes_to_plot_df$signature)),
             multiple = FALSE,
             selected = "Early_CC",
             options = list(`live-search` = TRUE)
           ),
           uiOutput("gene_select_ui"),
           tags$hr(),
           tags$div(class = "section-title", "Define Meta-genes"),
           uiOutput("metagene_ui"),
           actionBttn("add_metagene", "Add Meta-gene", color = "success", style = "simple"),
           tags$hr(),
           tags$div(class = "section-title", "Add Annotation Track"),
           pickerInput(
             "heatmap_ann_vars",
             "Add annotation variable(s) for heatmap (in addition to Subtype 5 Class)",
             choices = setdiff(names(metadata), "subtype_5_class"),
             multiple = TRUE,
             options = list(`live-search` = TRUE)
           ),
           uiOutput("ann_color_pickers"),
           tags$hr(),
           fluidRow(
             column(
               width = 6,
               actionBttn("plot_heatmap", "Plot Heatmap", width = "100%", style = "simple", color = "success"),
               tags$div(style = "height:20px") 
             ),
             column(
               width = 6,
               div(style = "display: flex; align-items: center; height: 100%; justify-content: flex-end;",
                   prettySwitch("sort_heatmap_samples", "Sort samples by Proliferation Score", value = TRUE, status = "info"),
                   prettySwitch("main_heatmap_subtype", "Show 7-class subtype annotation", value = FALSE, status = "info")
               )
             )
           ),
           plotlyOutput("heatmap_plot"),
           tags$div(style = "height:20px"),
           uiOutput("heatmap_legends")
  )
}