survival_tab_ui <- function() {
  tabPanel("Gene Expression Outcomes",
           fluidRow(
             column(5,
                    tags$div(class = "main-title", "Gene Expression Outcomes"),
                    tags$div(style = "margin-bottom:12px; font-size:1.08em; color:#2c3e50;",
                             "This tab allows you to perform survival analysis using gene expression or meta-gene data. Kaplan-Meier curves are generated for selected genes or meta-genes, stratified by expression level and clinical endpoint. All analyses use the currently filtered sample set."
                    ),
                    tags$div(
                      style = "margin-bottom:8px;",
                      strong("How to construct a survival plot:"),
                      tags$ul(
                        tags$li("Select a gene or meta-gene to stratify samples."),
                        tags$li("Choose a clinical endpoint and cutoff method for grouping."),
                        tags$li("Optionally filter samples by subtype."),
                        tags$li("The plot and risk table will update automatically for selected samples.")
                      )
                    ),
                    tags$div(class = "section-title", "Kaplan-Meier Plot"),
                    selectizeInput("km_gene", "Select gene or meta-gene",
                                   choices = expr_gene_names, selected = NULL, options = list(server = TRUE)),
                    selectInput(
                      "km_endpoint", "Select clinical endpoint",
                      choices = c(
                        "Biological Progression",
                        "Clinical Progression",
                        "Clinical Progression During Follow-up",
                        "BCG Any",
                        "BCG Adequate",
                        "Recurrance"
                      ),
                      selected = "Clinical Progression"
                    ),
                    pickerInput(
                      "km_filter_subtypes",
                      "Filter samples by subtype",
                      choices = levels(metadata$subtype_5_class),
                      selected = levels(metadata$subtype_5_class),
                      multiple = TRUE,
                      options = list(`actions-box` = TRUE, `live-search` = TRUE)
                    ),
                    radioButtons("km_cut", "Cutoff method",
                                 choices = c("Median" = "median", "Tertiles" = "tertile", "Quartiles" = "quartile"),
                                 inline = TRUE
                    ),
                    checkboxInput("km_show_confint", "Show confidence intervals", value = TRUE)
             ),
             column(7,
                    plotOutput("km_plot", height = "700px")
             )
           ),
           fluidRow(
             column(12,
                    # The risk table is part of the ggplot output in survminer::ggsurvplot
             )
           )
  )
}