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


#' full_mode_analysis
#'
#' Runs analysis on CircaN, JTK and MetaCycle and integrates the results.
#' @param data The data matrix with the omics data.
#' @param s2c Data frame with at least the columns sample (which must coincide with the samples in data) and time (recording
#' time corresponding to each sample.
#' @param algorithms character vector indicating which algorithms to test the data with. Must beone, or  a combination
#' of circan, jtk and metacycle.
#' @param min_per Minimum period to search for.
#' @param max_per Maximum period to search for.
#' @param circan_init_value Initial value for the period of CircaN regression.
#' @param circan_mode Algorithm to use in CircaN regression. Must be one of 'default' for Gauss-Newton, 'plinear' for the Golub-Pereyra algorithm
#' for partially linear least-squares models and 'port' for the ‘nl2sol’ algorithm from the Port library.
#' @param mc_cycMethod a character vector(length 1 or 2 or 3) for MetaCycle. User-defined methods for detecting rhythmic signals, must be
#' selected as any one, any two or all three methods(default) from "ARS"(ARSER), "JTK"(JTK_CYCLE) and "LS"(Lomb-Scargle).
#' @keywords circan, jtk, metacycle
#' @export
#' @examples
#' results <- full_mode_analysis(my_data, s2c)
full_mode_analysis <- function(data, s2c, algorithms = c("circan", "jtk", "metacycle"), min_per = 20, max_per = 28
                               , circan_init_value = 24, circan_mode="port", mc_cycMethod = c("LS", "JTK")){
  ################
  # setwd("/data3/arubio/projects/Andrea_CircaN_rebuttal/")
  # algorithms = c("circan", "jtk", "metacycle")
  # data <- read.csv("/data3/arubio/projects/Andrea_CircaN_rebuttal/reports/bio_data/fitzgerald/files/fitzgerald_norm_data.csv", row.names = 1)[1:200,]
  # s2c <- read.csv("/data3/arubio/projects/Andrea_CircaN_rebuttal/reports/bio_data/fitzgerald/files/fitzgerald_meta.csv")
  # min_per = 20
  # max_per = 28
  # circan_init_value = 24
  # circan_mode="port"
  # mc_cycMethod = c("LS", "JTK")
  # circan_files <- paste("/data3/arubio/src/dependencies/CircaN/"
  #                     , list.files(paste("/data3/arubio/src/dependencies/CircaN/", sep ="")), sep="")
  # sapply(circan_files, source)
  library("dplyr")
  ################
  results <- data.frame()
  if(any(grepl("circan", tolower(algorithms), fixed = TRUE))){
    data_aux <- data
    data_aux$features <- rownames(data_aux)
    data_aux <- data_aux %>% dplyr::select("features", everything())
    cat("Running CircaN")
    circan_results <- circan(data = data_aux, meta = s2c, mode = circan_mode
                           , init_value = circan_init_value, min_per = min_per, max_per = max_per)
    rm(data_aux)
    circan_results <- circan_results[order(circan_results[,1]),]
    colnames(circan_results)[-1] <- paste("circan", colnames(circan_results[,-1]), sep="_")
    results <- as.data.frame(cbind.fill(results,circan_results), stringsAsFactors = F)
    rm(circan_results)
  }
  if(any(grepl("jtk", tolower(algorithms), fixed = TRUE))){
    cat("Running JTK")
    jtk_results <- jtk_wrapper(data = data, s2c = s2c, min_per = min_per, max_per = max_per)
    jtk_results <- jtk_results[order(jtk_results[,1]),]
    colnames(jtk_results)[-1] <- paste("jtk", colnames(jtk_results[,-1]), sep="_")
    results <- as.data.frame(cbind.fill(results, jtk_results), stringsAsFactors = F)
    rm(jtk_results)
  }
  if(any(grepl("metacycle", tolower(algorithms), fixed = TRUE))){
    cat("Running MetaCycle")
    library("MetaCycle")
    mc_results <- meta2d_wrapper(data = data, timepoints = s2c$time
                               , min_per = min_per, max_per = max_per, cycMethod = mc_cycMethod)
    colnames(mc_results)[-1] <- paste("mc", colnames(mc_results[,-1]), sep="_")
    results <- as.data.frame(cbind.fill(results, mc_results), stringsAsFactors = F)
    rm(mc_results)
  }
  return(results)
}
