metadata_tab_ui <- function() {
  tabPanel("Metadata",
           fluidRow(
             column(5,
                    tags$div(class = "main-title", "Metadata"),
                    # Add predefined sample set selector
                    selectInput("predefined_sample_set", 
                                "Select predefined sample set:", 
                                choices = names(predefined_sample_lists),
                                selected = "All Samples"),
                    tags$div(style = "margin-bottom:15px; font-style:italic; color:#666;",
                             textOutput("predefined_set_info")),
                    
                    tags$div(style = "margin-bottom:8px; font-size:1.15em; color:#2c3e50;",
                             tags$p("You can subset your data in three ways:"),
                             tags$ul(
                               tags$li("Select a predefined sample set from a specific publication above"),
                               tags$li("Apply one or more filters below to interactively select samples based on metadata"),
                               tags$li("Upload a txt file (one sample per line) to subset directly to those samples")
                             ),
                             tags$p("If no filters are added and no file is uploaded, all samples are included.")
                    ),
                    
                    fileInput("sample_id_file", "Upload sample ID list", accept = ".txt"),
                    tags$div(id = "filter_ui_container"),
                    actionBttn("add_filter", "Add Filter", style = "simple", color = "success"),
                    tags$div(style = "height:30px"),
                    tags$div(style = "height:18px"),
                    tags$div(style = "height:18px"),
                    downloadBttn("download_ids", "Download Sample IDs", style = "unite", color = "success")
             ),
             column(7,
                    tags$div(class = "section-title", "Metadata Table"),
                    DT::dataTableOutput("table")
             )
           )
  )
}