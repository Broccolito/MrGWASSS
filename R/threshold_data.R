#' Threshold the GWAS Summary Statistics
#'
#' This function thresholds the GWAS summary statistics based on the specified p value threshold
#'
#' @param data an data frame object loaded from load_data function
#' @param threshold the p-value threshold specified, usually at 1e-5 or 1e-8
#' @return data_filtered, a data frame of the filtered summary statistics
#' @export
threshold_data = function(data, threshold = 1e-5){
  data_filtered = data %>%
    filter(pvalue <= threshold) %>%
    arrange(pvalue)
  return(data_filtered)
}