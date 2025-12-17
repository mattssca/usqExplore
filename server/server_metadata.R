metadata_server <- function(input, output, session) {
  
  # Reactive values for filters
  filter_count <- reactiveVal(0)
  filter_ids <- reactiveVal(character(0))
  
  # Predefined sample set info
  output$predefined_set_info <- renderText({
    if (is.null(input$predefined_sample_set)) return("")
    
    set_name <- input$predefined_sample_set
    if (set_name %in% names(predefined_sample_lists)) {
      sample_count <- length(predefined_sample_lists[[set_name]])
      paste("This set contains", sample_count, "samples")
    } else {
      ""
    }
  })
  
  # Add Filter button logic
  observeEvent(input$add_filter, {
    # Store current filter input values before adding new filter
    ids_prev <- isolate(filter_ids())
    filter_values <- lapply(ids_prev, function(fid) {
      list(
        cat = input[[paste0(fid, "_cat")]],
        var = input[[paste0(fid, "_var")]],
        val = input[[paste0(fid, "_val")]]
      )
    })
    names(filter_values) <- ids_prev
    
    count <- filter_count() + 1
    filter_count(count)
    id <- paste0("filter_", count)
    ids <- c(filter_ids(), id)
    filter_ids(ids)
    
    insertUI(
      selector = "#filter_ui_container",
      where = "beforeEnd",
      ui = tags$div(id = id,
                    fluidRow(
                      column(5,
                             pickerInput(paste0(id, "_cat"), "Category", choices = names(var_categories), width = "100%", options = list(`live-search` = TRUE)),
                             uiOutput(paste0(id, "_var_ui")),
                             uiOutput(paste0(id, "_val_ui")),
                             actionBttn(paste0(id, "_remove"), "Remove", style = "simple", color = "danger", size = "md"),
                             tags$div(style = "height:18px")
                      )
                    )
      )
    )
    
    # Restore previous filter selections for all filters after UI is rendered
    later::later(function() {
      for (i in seq_along(ids_prev)) {
        fid <- ids_prev[i]
        vals <- filter_values[[fid]]
        if (!is.null(vals$cat)) updatePickerInput(session, paste0(fid, "_cat"), selected = vals$cat)
        # Wait for var UI to be rendered before updating var and val
        later::later(function() {
          if (!is.null(vals$var)) updatePickerInput(session, paste0(fid, "_var"), selected = vals$var)
          later::later(function() {
            if (!is.null(vals$val)) {
              if (!is.null(vals$var) && is.numeric(metadata[[vals$var]])) {
                updateSliderInput(session, paste0(fid, "_val"), value = vals$val)
              } else {
                updatePickerInput(session, paste0(fid, "_val"), selected = vals$val)
              }
            }
          }, delay = 0.1)
        }, delay = 0.1)
      }
    }, delay = 0.2)
  })
  
  # Variable selection UI for each filter
  observe({
    lapply(filter_ids(), function(id) {
      output[[paste0(id, "_var_ui")]] <- renderUI({
        cat_input <- input[[paste0(id, "_cat")]]
        if (is.null(cat_input)) return(NULL)
        pickerInput(paste0(id, "_var"), "Variable", choices = var_categories[[cat_input]], width = "100%",
                    options = list(`live-search` = TRUE))
      })
    })
  })
  
  # Remove filter observers
  observe({
    lapply(filter_ids(), function(id) {
      observeEvent(input[[paste0(id, "_remove")]], {
        removeUI(selector = paste0("#", id))
        filter_ids(setdiff(filter_ids(), id))
      }, ignoreInit = TRUE)
    })
  })
  
  # Value selection UI for each filter
  observe({
    lapply(filter_ids(), function(id) {
      output[[paste0(id, "_val_ui")]] <- renderUI({
        var_input <- input[[paste0(id, "_var")]]
        if (is.null(var_input)) return(NULL)
        vals <- metadata[[var_input]]
        if (is.numeric(vals)) {
          rng <- range(vals, na.rm = TRUE)
          sliderInput(
            paste0(id, "_val"), "Range",
            min = rng[1], max = rng[2], value = rng, step = diff(rng)/100,
            width = "100%"
          )
        } else {
          pickerInput(paste0(id, "_val"), "Value(s)", choices = unique(vals), multiple = TRUE, width = "100%",
                      options = list(`live-search` = TRUE))
        }
      })
    })
  })
  
  # Main data filtering reactive
  data_filtered <- reactive({
    df <- metadata
    
    # Apply predefined sample set first
    if (!is.null(input$predefined_sample_set) && input$predefined_sample_set != "All Samples") {
      predefined_ids <- predefined_sample_lists[[input$predefined_sample_set]]
      df <- df[df$sample_id %in% predefined_ids, , drop = FALSE]
    }
    
    # If a sample ID file is uploaded, subset to those sample IDs
    if (!is.null(input$sample_id_file)) {
      ids <- tryCatch({
        readLines(input$sample_id_file$datapath)
      }, error = function(e) NULL)
      ids <- trimws(ids)
      ids <- ids[ids != ""]
      df <- df[df$sample_id %in% ids, , drop = FALSE]
    }
    
    # Apply manual filters
    for (id in filter_ids()) {
      var <- input[[paste0(id, "_var")]]
      val <- input[[paste0(id, "_val")]]
      if (!is.null(var) && !is.null(val)) {
        if (is.numeric(metadata[[var]])) {
          df <- df[df[[var]] >= val[1] & df[[var]] <= val[2], , drop = FALSE]
        } else if (length(val) > 0) {
          df <- df[df[[var]] %in% val, , drop = FALSE]
        }
      }
    }
    df
  })
  
  # Render metadata table
  output$table <- DT::renderDataTable({
    data_filtered()
  }, options = list(pageLength = 10))
  
  # Download filtered sample IDs
  output$download_ids <- downloadHandler(
    filename = function() {
      paste0("sample_ids_", Sys.Date(), ".txt")
    },
    content = function(file) {
      writeLines(as.character(data_filtered()$sample_id), file)
    }
  )
  
  # Return the filtered data for other modules
  return(list(data_filtered = data_filtered))
}