#' @title Get HM Order.
#'
#' @description Return a sample order based on the late early ratio.
#'
#' @details Modified version of the LundTax2023 heatmap function, this function only returns the 
#' order of samples to be used for downstream plotting functions. The order is dictated by the early 
#' late ratio.
#' 
#' @param these_predictions Required parameter, should be the output from 
#' [LundTax2023Classifier::lundtax_predict_sub()]. Note, this function required the `lundtax_predict_sub` 
#' to be run with `include_data = TRUE` argument.
#' @param norm Boolean parameter. Set to TRUE (default) to normalize the data into Z-scaled values.
#' @param return_this Set to "late_early" to return late/early ratio. Set to "sample_order" to 
#' return the sample order based on late/early ratio (default).
#' 
#' @return Rreturns the sample order.
#' 
#' @import dplyr
#' 
#' @export 
#'
#' @examples
#' #get predictions with data
#' my_predictions = lundtax_predict_sub(this_data = sjodahl_2017,
#'                                      include_data = TRUE, 
#'                                      gene_id = "hgnc_symbol", 
#'                                      impute = TRUE, 
#'                                      adjust = TRUE)
#'                                      
#' #run function                                    
#' get_sample_order(these_predictions = my_predictions)
#'
get_sample_order = function(these_predictions = NULL,
                            return_this = "sample_order"){
  
  #check incoming data and parameter combinations
  this_data = as.matrix(these_predictions$data)
  score_matrix = these_predictions$subtype_scores
  pred_labels5 = these_predictions$predictions_5classes
  pred_labels7 = these_predictions$predictions_7classes
  
  #gene signatures for plotting
  genes_to_plot = list(Early_CC = c(signatures$proliferation[which(signatures$proliferation$signature == "EarlyCellCycle"), 1]),
                       Late_CC = c(signatures$proliferation[which(signatures$proliferation$signature == "LateCellCycle"), 1]),
                       Late_Early = NULL)
  
  #get genes in provided data for downstream filtering steps
  these_genes = row.names(this_data)
  genes_early = intersect(genes_to_plot$Early_CC, these_genes)
  genes_late = intersect(genes_to_plot$Late_CC, these_genes)
  
  #combine genes
  genes_cc = na.omit(c(genes_early, genes_late))
    
  #row split for the heatmap
  row_split = c(rep("Early", length(genes_early)),
                rep("Late", length(genes_late)))
  #row title
  row_title_cc = c("Early Cell Cycle", "Late Cell Cycle")
    
  #late and Early scores
  late_score = apply(this_data[intersect(rownames(this_data),genes_to_plot$Late_CC),], 2, median)
  early_score = apply(this_data[intersect(rownames(this_data),genes_to_plot$Early_CC),], 2, median)
    
  #ratio
  late_early_ratio = late_score - early_score
    
  #add genes to genes_to_plot object
  genes_to_plot$Late_Early = late_early_ratio

  #order samples by late_early cell cycle
  sample_order = order(late_early_ratio)
  
  #end function
  print(sample_order)
  
  if(return_this == "late_early"){
    message("Returning late/early ratio...")
    return(as.data.frame(late_early_ratio))
  }else if(return_this == "sample_order"){
    message("Returning sample order...")
    return(sample_order)
  }
}
