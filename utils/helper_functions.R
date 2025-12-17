# Bin Numeric Variables (internal)
int_bin_numeric_variables <- function(this_data = NULL, num_bins = NULL){
  numeric_columns = sapply(this_data, is.numeric)
  bin_column = function(column, num_bins){
    range_min = min(column, na.rm = TRUE)
    range_max = max(column, na.rm = TRUE)
    binned_column = cut(column, 
                        breaks = seq(range_min, range_max, length.out = num_bins + 1), 
                        labels = FALSE, 
                        include.lowest = TRUE)
    return(binned_column)
  }
  this_data[numeric_columns] <- lapply(this_data[numeric_columns], bin_column, num_bins = num_bins)
  return(this_data)
}

# Helper to get current meta-gene input values
get_current_metagene_inputs <- function(mg_list, input) {
  out <- list()
  for (mg_id in names(mg_list)) {
    out[[mg_id]] <- list(
      name = input[[paste0(mg_id, "_name")]],
      genes = input[[paste0(mg_id, "_genes")]],
      filter = input[[paste0(mg_id, "_filter")]]
    )
  }
  out
}
