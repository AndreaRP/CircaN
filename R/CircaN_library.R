#' circan
#'
#' Nonlinear least squares model for accurate detection of circadian expression patterns.
#' @param data Dataframe containing the expression data. Samples must be in columns and genes in rows.
#' For an example see data(expression_example).
#' @param meta Dataframe containing the metadata for the samples. Must have at least a 'sample' column
#' with the sample name as it appears in the data matrix;
#' a 'time' column with the time point the sample was collected;
#' and an 'ind' column containing information for the individual the sample comes from.
#' For an example see data(metadata_example).
#' @param shiny Is the package running in a shiny app? default to FALSE.
#' @param mode Algorithm to use in the NLS regression. Must be one of 'default' for Gauss-Newton, 'plinear' for the Golub-Pereyra algorithm
#' for partially linear least-squares models and 'port' for the ‘nl2sol’ algorithm from the Port library. Default is default. See nls documentation
#' for extended info.
#' @param init_value Initial value for the period. Default is set to 24.
#' @param max_per Maximum period to regress. Default is set to Inf.
#' @param min_per Minimum period to regress. Default is set to -Inf.
#' @keywords CircaN circadian regression
#' @export
#' @examples
#' # This runs CircaN on the example data with the 'Port' algorithm.
#' circan(data=expression_example, meta=metadata_example, mode="port")
circan <- function(data, meta, shiny=FALSE, mode="default", init_value=24, max_per=Inf, min_per=-Inf){
  print("Fitting sine wave...")
  circan_sine_results <- circan_sine(data=data, s2c=meta, mode=mode, init_value=init_value, min_per = min_per, max_per = max_per)
  print("Fitting triangular wave...")
  circan_triangular_results <- circan_triangular(data=data, s2c=meta, mode=mode, init_value=init_value, min_per = min_per, max_per = max_per)
  print("Fitting peak wave...")
  circan_peak_results <- circan_peak(data=data, s2c=meta, mode=mode, init_value=init_value, min_per = min_per, max_per = max_per)
  print("Fitting linear + sine wave...")
  circan_line_results <- circan_line(data=data, s2c=meta, mode=mode, init_value=init_value, min_per = min_per, max_per = max_per)
  print("Fitting exponential wave...")
  circan_exp_results <- circan_exp(data=data, s2c=meta, mode=mode, init_value=init_value, min_per = min_per, max_per = max_per)
  print("Fitting damped wave wave...")
  circan_damp_results <- circan_damp(data=data, s2c=meta, mode=mode, init_value=init_value, min_per = min_per, max_per = max_per)
  print("Fitting cosine2 wave...")
  circan_cosine2_results <- circan_cosine2(data=data, s2c=meta, mode=mode, init_value=init_value, min_per = min_per, max_per = max_per)
  print("Calculating best fit...")
  results_list <- grep(glob2rx("circan*results"), ls(), value=TRUE)
  # Select best fitting curve
  circan_results <- best_fit_select(results_list=mget(results_list))
  circan_results$combined_pval <- get_fisher(as.numeric(circan_results$p.value.amp), as.numeric(circan_results$p.value.per))
  circan_results$BH_combined <- p.adjust(as.numeric(circan_results$combined_pval), method = "BH")
  return(circan_results)
}
