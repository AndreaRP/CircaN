setClass("harmony_gene", slots=list(vec1="numeric", vec2="numeric", r1="numeric", r2="numeric", gs.vec="numeric"))


#############################
setwd("/home/arubio/Dropbox/CircaN/R/")
# list.files()
source("Utils.R")
source("circan_sine.R")
source("circan_peak.R")
source("circan_cosine2.R")
source("circan_damp.R")
source("circan_exp.R")
source("circan_line.R")
source("circan_square.R")
source("circan_triangular.R")
load("../data/expression_example.RData")
load("../data/metadata_example.RData")
expression_example <- expression_example[1:20,]

p <- circan(data=expression_example, meta=metadata_example, mode="port")
#############################

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


#' circan
#'
#' Nonlinear least squares model for accurate detection of circadian expression patterns.
#' @param data Dataframe containing the expression data. Samples must be in columns and genes in rows.
#' For an example see data(expression_example).
#' @param s2c Dataframe containing the metadata for the samples. Must have at least a 'sample' column
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
#' circan(data=expression_example, s2c=metadata_example, mode="port")

circan <- function(data, s2c, shiny=FALSE, mode="default", init_value=24, max_per=Inf, min_per=-Inf){
  s2c$time <- as.numeric(as.character(s2c$time))
  s2c$ind <- as.integer(as.character(s2c$ind))
  s2c$sample <- as.character(s2c$sample)
  # metadata must have sample, ind and time columns. Data is a data frame or matrix
  # with the ids in the first column
  t <- unique(as.numeric(as.character(s2c$time)))
  data0 <- data
  total_genes <- nrow(data)
  rownames(data) <- as.character(data[,1])
  data <- data[,-1]
  data <- as.matrix(data)
  results.cols <- c("feature"
                    ,"estimate.amp","std.error.amp","statistic.amp","p.value.amp"
                    ,"estimate.phase","std.error.phase","statistic.phase","p.value.phase"
                    ,"estimate.per" ,"std.error.per","statistic.per", "p.value.per"
                    , "r")
  results <- setNames(data.frame(matrix(ncol = 14, nrow = 0)), results.cols)
  for(gene in 1:nrow(data)){
    # Standarize data
    data.st <- (data[gene,]-mean(data[gene,]))/sqrt(var(data[gene,]))
    # Group data by timepoint
    df0 <- merge(data.st, s2c, by.x = 0, by.y = "sample")
    df0 <- dplyr::rename(df0, data = x)
    df <- as.data.frame(df0[,c("data", "ind", "time")])
    df$data <- as.numeric(df$data)
    df$time <- as.numeric(df$time)
    gd <- nlme::groupedData(data~time|ind,data=df)
    # plot(gd)
    result = tryCatch({
      # Try fitting all parameters
      nls.model=nls(data~amp*sin(2*pi*time/per+phase)
                    , start=list(amp=1, phase=0, per=init_value)
                    , algorithm = mode
                    , lower = list(per=min_per)
                    , upper = list(per=max_per)
                    , data = gd
                    , control = list(maxiter = 200)
      )
      stats <- as.data.frame(broom::tidy(nls.model))
      rownames(stats) <- stats[,1]
      stats <- stats[,-1]
      vec <- as.numeric(c(t(stats)))
      names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
      # R (goodness of fit)
      r <- cor(gd$data, predict(nls.model))
    }, error = function(e1) {
      convergence_error1 <- grepl("Convergence failure", as.character(e1))
      if(convergence_error1){
        tryCatch({
          # Try with fixed period
          nls.model=nls(  data~amp*sin(2*pi*time/init_value+phase)
                          , start=list(amp=1,phase=0)
                          , algorithm = mode
                          , data=gd
                          , control = list(maxiter = 200
                                           # , warnOnly = TRUE
                          )
          )
          stats <- as.data.frame(broom::tidy(nls.model))
          rownames(stats) <- stats[,1]
          stats <- stats[,-1]
          vec <- as.numeric(c(t(stats)))
          names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
          # As we haven't estimated the period we add the data to vec
          temp <- c(24, "NA", "NA", "NA")
          names(temp) <- c("estimate.period", "std.error.period"
                           , "statistic.period", "p.value.period")
          vec <<- c(vec, temp)
          # R (goodness of fit)
          r <<- cor(gd$data, predict(nls.model))
        }, error = function(e2) {
          convergence_error2 <- grepl("Convergence failure", as.character(e2))
          if(convergence_error2){
            # Can't be adjusted to a circadian curve
            vec <<- rep("NA", times = ncol(results)-1)
            r <<- "NA"
          }
        })
      }
    }, finally = {
      res <- rbind(c(feature=rownames(data)[gene], vec, r = r))
      results[gene,] <- res
    })
    if(shiny) incProgress(amount=1/total_genes)
  }
  # Calc adjusted p.val for period
  results$BH.q.value.per <- p.adjust(as.numeric(results$p.value.per), method="BH")
  results$BH.q.value.amp <- p.adjust(as.numeric(results$p.value.amp), method="BH")
  results$BH.q.value.phase <- p.adjust(as.numeric(results$p.value.phase), method="BH")
  return(results)
}

#' circan_groups
#'
#' Nonlinear least squares model for accurate detection of circadian expression patterns and
#' group comparison.
#' @param data Dataframe containing the expression data of samples from both groups. Samples must be in columns and genes in rows.
#' @param s2c Dataframe containing the metadata for the samples. Must have at least a 'sample' column
#' with the sample name as it appears in the data matrix;
#' a 'time' column with the time point the sample was collected;
#' and an 'ind' column containing information for the individual the sample comes from.
#' @param id_col Column in the expression dataframe containing the feature IDs.
#' @param sample_col Column in the metadata dataframe containing the names of the samples.
#' @param time_col Column in the metadata dataframe containing the time point of each sample.
#' @param group_col Column in the metadata dataframe containing the group to which sample belongs.
#' @param ind_col Column in the metadata dataframe containing the individual from which each sample comes.
#' Please note that the individuals from both groups must be different.
#' @param shiny Is the package running in a shiny app? default to FALSE.
#' @param mode Algorithm to use in the NLS regression. Must be one of 'default' for Gauss-Newton, 'plinear' for the Golub-Pereyra algorithm
#' for partially linear least-squares models and 'port' for the ‘nl2sol’ algorithm from the Port library. Default is default. See nls documentation
#' for extended info.
#' @param init_value Initial value for the period. Default is set to 24.
#' @param max_per Maximum period to regress. Default is set to Inf.
#' @param min_per Minimum period to regress. Default is set to -Inf.
#' @keywords CircaN circadian regression, group comparison.
#' @export


circan_groups <- function(data, s2c, id_col, sample_col, time_col
                          , group_col, ind_col, shiny=FALSE, mode="default"
                          , init_value=24, max_per=Inf, min_per=-Inf){
  
  s2c[,time_col] <- as.numeric(as.character(s2c[,time_col]))
  s2c[,ind_col] <- as.integer(as.character(s2c[,ind_col]))
  s2c[,sample_col]<- as.character(s2c[,sample_col])
  s2c[,group_col] <- as.character(s2c[,group_col])
  # Formatting data
  rown <- data[,id_col]
  data <- data[,-id_col]
  data[,1:ncol(data)] <- sapply(data[,1:ncol(data)], as.numeric)
  rownames(data) <- rown
  
  G1_name <- unique(s2c[,group_col])[1]
  G2_name <- unique(s2c[,group_col])[2]
  data1 <- data[,s2c[which(s2c[,group_col] == G1_name),sample_col]]
  data2 <- data[,s2c[which(s2c[,group_col] == G2_name),sample_col]]
  
  results.cols <- c("estimate.amp","std.error.amp","statistic.amp","p.value.amp"
                    ,"estimate.phase","std.error.phase","statistic.phase","p.value.phase"
                    ,"estimate.per" ,"std.error.per","statistic.per", "p.value.per"
                    , "r")
  group.cols <- c("feature", "estimate.Amp2_Amp1", "std.error.Amp2_Amp1", "statistic.Amp2_Amp1", "p.value.Amp2_Amp1"
                  , "estimate.Ph2_Ph1", "std.error.Ph2_Ph1", "statistic.Ph2_Ph1", "p.value.Ph2_Ph1"
                  , "estimate.Per2_Per1", "std.error.Per2_Per1", "statistic.Per2_Per1", "p.value.Per2_Per1")
  
  G1 <- setNames(data.frame(matrix(ncol = length(results.cols)+1, nrow = 0))
                 , c("feature", paste(results.cols, "1", sep="")))
  G2 <- setNames(data.frame(matrix(ncol = length(results.cols)+1, nrow = 0))
                 , c("feature", paste(results.cols, "2", sep="")))
  comp <- setNames(data.frame(matrix(ncol = length(group.cols), nrow = 0))
                   , group.cols)
  # results <- new("harmony_data", G1=G1, G2=G2, Comparison=comp)
  results <- new("harmony_data", G1=G1, G2=G2, Comparison=comp
                 , run.info=list(G1=G1_name, G2=G2_name, id_col=id_col, sample_col=sample_col
                                 , time_col=time_col, group_col=group_col
                                 , ind_col=ind_col, mode=mode, init_value=init_value
                                 , max_per=max_per, min_per=min_per))
  env=new.env()
  assign("h_gene", new("harmony_gene"), env=env)
  
  for(gene in 1:nrow(data)){
    # Center data
    data.cent.1 <- (data1[gene,]-mean(as.numeric(data1[gene,])))
    data.cent.2 <- (data2[gene,]-mean(as.numeric(data2[gene,])))
    # Merged data
    merged_data <- cbind(data.cent.1, data.cent.2)
    rownames(merged_data) <- "data"
    # Construct df0 (ind must be different between groups)
    df0 <- merge(t(merged_data), s2c, by.x = 0, by.y = sample_col)
    # Change names to match formula
    df0 <- dplyr::rename(df0, group = colnames(s2c)[group_col]
                         , time = colnames(s2c)[time_col]
                         , ind = colnames(s2c)[ind_col])
    df0 <- df0[order(df0$group, df0$time),]
    df0$data <- as.numeric(df0$data)
    df0$time <- as.numeric(df0$time)
    # Group
    gd <- nlme::groupedData(data~time|ind, data=df0)
    # Possible try
    result = tryCatch({
      nls.model=stats::nls(data~
                             as.numeric(group==unique(s2c[,group_col])[1])*amp1*sin(2*pi*time/per1+phase1) +
                             as.numeric(group==unique(s2c[,group_col])[2])*amp2*sin(2*pi*time/per2+phase2)
                           , algorithm = mode
                           , lower = list(per=min_per)
                           , upper = list(per=max_per)
                           , data = gd
                           , start = c(amp1=1, phase1=0, per1=init_value, amp2=1, phase2=0, per2=init_value)
      )
      # Each group' stats separately
      stats <- as.data.frame(broom::tidy(nls.model))
      rownames(stats) <- stats[,1]
      coeff <- stats$term
      stats <- stats[,-1]
      # Populate G1 and G2
      vec <- as.numeric(c(t(stats)))
      names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
      vec1 <- vec[grep("1", names(vec))]
      vec2 <- vec[grep("2", names(vec))]
      # R (goodness of fit)
      r1 <- cor(gd[which(gd$group==unique(s2c[,group_col])[1]),"data"]
                , predict(nls.model)[1:length(which(gd$group==unique(s2c[,group_col])[1]))])
      r2 <- cor(gd[which(gd$group==unique(s2c[,group_col])[2]),"data"]
                , predict(nls.model)[(length(which(gd$group==unique(s2c[,group_col])[1]))+1):length(gd$group)])
      
      # Groups comparison
      K <- rbind(c(-1, 0, 0, 1 , 0, 0)
                 , c( 0, -1, 0, 0, 1, 0)
                 , c( 0, 0, -1, 0, 0, 1)
      )
      rownames(K) <- c("Amp2_Amp1", "Ph2_Ph1", "Per2_Per1")
      colnames(K) <- names(coef(nls.model))
      
      group.stats <- as.data.frame(broom::tidy(summary(multcomp::glht(nls.model, linfct = K))))
      rownames(group.stats) <- group.stats[,1]
      coeff.group <- group.stats$lhs
      group.stats <- group.stats[,-c(1:2)]
      gs.vec <- as.numeric(c(t(group.stats)))
      names(gs.vec) <- c(outer(colnames(group.stats), rownames(group.stats), paste, sep = "."))
      assign("h_gene", new("harmony_gene", vec1=vec1, vec2=vec2, r1=r1, r2=r2, gs.vec=gs.vec), env=env)
    }, error = function(err1) {
      convergence_error1 <- grepl("Convergence failure", as.character(err1))
      if (convergence_error1){
        result2 = tryCatch({ ## Per 1 fixed
          # Try with fixed period1
          nls.model=stats::nls(data~
                                 as.numeric(group==unique(s2c[,group_col])[1])*amp1*sin(2*pi*time/init_value+phase1) +
                                 as.numeric(group==unique(s2c[,group_col])[2])*amp2*sin(2*pi*time/per2+phase2)
                               , algorithm = mode
                               , lower = list(per=min_per)
                               , upper = list(per=max_per)
                               , data = gd
                               , start = c(amp1=1, phase1=0, amp2=1, phase2=0, per2=init_value)
          )
          # statistic values
          stats <- as.data.frame(broom::tidy(nls.model))
          rownames(stats) <- stats[,1]
          stats <- stats[,-1]
          vec <- as.numeric(c(t(stats)))
          names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
          # As we haven't estimated period 1 we add the data to vec1
          vec1 <- vec[grep("1", names(vec))]
          temp <- c(init_value, NA, NA, NA)
          names(temp) <- c("estimate.per1", "std.error.per1"
                           , "statistic.per1", "p.value.per1")
          vec1 <- c(vec1, temp)
          mode(vec1) <- "numeric"
          vec2 <- vec[grep("2", names(vec))]
          # R (goodness of fit)
          r1 <- cor(gd[which(gd$group==unique(s2c[,group_col])[1]),"data"]
                    , predict(nls.model)[1:length(which(gd$group==unique(s2c[,group_col])[1]))])
          r2 <- cor(gd[which(gd$group==unique(s2c[,group_col])[2]),"data"]
                    , predict(nls.model)[(length(which(gd$group==unique(s2c[,group_col])[1]))+1):length(gd$group)])
          
          # Groups comparison
          K <- rbind(c(-1, 0, 1, 0, 0)
                     , c(0, -1, 0, 1, 0)
          )
          rownames(K) <- c("Amp2_Amp1", "Ph2_Ph1")
          colnames(K) <- names(coef(nls.model))
          
          # group.stats <- tidy(summary(glht(nls.model, linfct = K)))
          group.stats <- as.data.frame(broom::tidy(summary(multcomp::glht(nls.model, linfct = K))))
          rownames(group.stats) <- group.stats[,1]
          coeff.group <- group.stats$lhs
          group.stats <- group.stats[,-c(1:2)]
          gs.vec <- as.numeric(c(t(group.stats)))
          names(gs.vec) <- c(outer(colnames(group.stats), rownames(group.stats), paste, sep = "."))
          gs.vec <- c(gs.vec, setNames(as.numeric(c(vec2["estimate.per2"]-vec1["estimate.per1"], NA, NA, NA))
                                       , c("estimate.Per2_Per1", "std.error.Per2_Per1"
                                           , "statistic.Per2_Per1", "p.value.Per2_Per1")))
          assign("h_gene", new("harmony_gene", vec1=vec1, vec2=vec2, r1=r1, r2=r2, gs.vec=gs.vec), env=env)
          
        }, error = function(err2) {
          convergence_error2 <- grepl("Convergence failure", as.character(err2))
          if(convergence_error2){
            result = tryCatch({
              # Fixed per2
              nls.model=stats::nls(data~
                                     as.numeric(group==unique(s2c[,group_col])[1])*amp1*sin(2*pi*time/per1+phase1) +
                                     as.numeric(group==unique(s2c[,group_col])[2])*amp2*sin(2*pi*time/init_value+phase2)
                                   , algorithm = mode
                                   , lower = list(per=min_per)
                                   , upper = list(per=max_per)
                                   , data = gd
                                   , start = c(amp1=1, phase1=0, amp2=1, phase2=0, per1=init_value)
              )
              # statistic values
              stats <- as.data.frame(broom::tidy(nls.model))
              rownames(stats) <- stats[,1]
              stats <- stats[,-1]
              vec <- as.numeric(c(t(stats)))
              names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
              
              vec1 <- vec[grep("1", names(vec))]
              # As we haven't estimated period 2 we add the data to vec2
              vec2 <- vec[grep("2", names(vec))]
              temp <- c(init_value, "NA", "NA", "NA")
              names(temp) <- c("estimate.per2", "std.error.per2"
                               , "statistic.per2", "p.value.per2")
              vec2 <- c(vec2, temp)
              mode(vec2) <- "numeric"
              
              # R (goodness of fit)
              r1 <- cor(gd[which(gd$group==unique(s2c[,group_col])[1]),"data"]
                        , predict(nls.model)[1:length(which(gd$group==unique(s2c[,group_col])[1]))])
              r2 <- cor(gd[which(gd$group==unique(s2c[,group_col])[2]),"data"]
                        , predict(nls.model)[(length(which(gd$group==unique(s2c[,group_col])[1]))+1):length(gd$group)])
              
              # Groups comparison
              K <- rbind(c(-1, 0, 1, 0, 0)
                         , c(0, -1, 0, 1, 0)
              )
              rownames(K) <- c("Amp2_Amp1", "Ph2_Ph1")
              colnames(K) <- names(coef(nls.model))
              
              group.stats <- as.data.frame(broom::tidy(summary(multcomp::glht(nls.model, linfct = K))))
              rownames(group.stats) <- group.stats[,1]
              coeff.group <- group.stats$lhs
              group.stats <- group.stats[,-c(1:2)]
              gs.vec <- as.numeric(c(t(group.stats)))
              names(gs.vec) <- c(outer(colnames(group.stats), rownames(group.stats), paste, sep = "."))
              gs.vec <- c(gs.vec, setNames(as.numeric(c(vec2["estimate.per2"]-vec1["estimate.per1"], NA, NA, NA))
                                           , c("estimate.Per2_Per1", "std.error.Per2_Per1"
                                               , "statistic.Per2_Per1", "p.value.Per2_Per1")))
              assign("h_gene", new("harmony_gene", vec1=vec1, vec2=vec2, r1=r1, r2=r2, gs.vec=gs.vec), env=env)
            }, error = function(err3) {
              convergence_error3 <- grepl("Convergence failure", as.character(err3))
              if(convergence_error3){
                # Fix per 1 and 2
                result = tryCatch({
                  # Fixed per2
                  nls.model=stats::nls(data~
                                         as.numeric(group==unique(s2c[,group_col])[1])*amp1*sin(2*pi*time/init_value+phase1) +
                                         as.numeric(group==unique(s2c[,group_col])[2])*amp2*sin(2*pi*time/init_value+phase2)
                                       , algorithm = mode
                                       , lower = list(per=min_per)
                                       , upper = list(per=max_per)
                                       , data = gd
                                       , start = c(amp1=1, phase1=0, amp2=1, phase2=0)
                  )
                  # statistic values
                  stats <- as.data.frame(broom::tidy(nls.model))
                  rownames(stats) <- stats[,1]
                  stats <- stats[,-1]
                  vec <- as.numeric(c(t(stats)))
                  names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
                  # As we haven't estimated period 1 we add the data to vec2
                  vec1 <- vec[grep("1", names(vec))]
                  temp <- c(init_value, NA, NA, NA)
                  names(temp) <- c("estimate.per1", "std.error.per1"
                                   , "statistic.per1", "p.value.per1")
                  vec1 <- c(vec1, temp)
                  mode(vec1) <- "numeric"
                  # As we haven't estimated period 2 we add the data to vec2
                  vec2 <- vec[grep("2", names(vec))]
                  temp <- c(init_value, NA, NA, NA)
                  names(temp) <- c("estimate.per2", "std.error.per2"
                                   , "statistic.per2", "p.value.per2")
                  vec2 <- c(vec2, temp)
                  mode(vec2) <- "numeric"
                  
                  # R (goodness of fit)
                  r1 <- cor(gd[which(gd$group==unique(s2c[,group_col])[1]),"data"]
                            , predict(nls.model)[1:length(which(gd$group==unique(s2c[,group_col])[1]))])
                  r2 <- cor(gd[which(gd$group==unique(s2c[,group_col])[2]),"data"]
                            , predict(nls.model)[(length(which(gd$group==unique(s2c[,group_col])[1]))+1):length(gd$group)])
                  
                  # Groups comparison
                  K <- rbind(c(-1, 0, 1, 0)
                             , c(0, -1, 0, 1)
                  )
                  rownames(K) <- c("Amp2_Amp1", "Ph2_Ph1")
                  colnames(K) <- names(coef(nls.model))
                  
                  # group.stats <- tidy(summary(glht(nls.model, linfct = K)))
                  group.stats <- as.data.frame(broom::tidy(summary(multcomp::glht(nls.model, linfct = K))))
                  rownames(group.stats) <- group.stats[,1]
                  coeff.group <- group.stats$lhs
                  group.stats <- group.stats[,-c(1:2)]
                  gs.vec <- as.numeric(c(t(group.stats)))
                  names(gs.vec) <- c(outer(colnames(group.stats), rownames(group.stats), paste, sep = "."))
                  gs.vec <- c(gs.vec, setNames(as.numeric(c(vec2["estimate.per2"]-vec1["estimate.per1"], NA, NA, NA))
                                               , c("estimate.Per2_Per1", "std.error.Per2_Per1"
                                                   , "statistic.Per2_Per1", "p.value.Per2_Per1")))
                  assign("h_gene", new("harmony_gene", vec1=vec1, vec2=vec2, r1=r1, r2=r2, gs.vec=gs.vec), env=env)
                }, error = function(err4) {
                  convergence_error4 <- grepl("Convergence failure", as.character(err4))
                  if(convergence_error4){
                    # Can't be adjusted to a circadian curve
                    vec1 = vec2 <- rep(NA, times = length(results.cols))
                    r1 = r2 <- NA
                    gs.vec <- rep(NA, times = length(group.cols)-1)
                    assign("h_gene", new("harmony_gene", vec1=as.numeric(vec1), vec2=as.numeric(vec2)
                                         , r1=as.numeric(r1), r2=as.numeric(r2), gs.vec=as.numeric(gs.vec))
                           , env=env)
                  }# end error 4
                })
              } # end error 3
            })
          }# end error 2
        })
      } # end error 1
    }, finally = {
      temp_gene <- get("h_gene", env=env)
      
      res1 <- rbind(c(feature=rownames(data)[gene], as.numeric(temp_gene@vec1), r1 = as.numeric(temp_gene@r1)))
      results@G1[gene,] <- res1
      res2 <- rbind(c(feature=rownames(data)[gene], as.numeric(temp_gene@vec2), r2 = as.numeric(temp_gene@r2)))
      results@G2[gene,] <- res2
      results@Comparison[gene,] <- as.matrix(rbind(c(feature=rownames(data)[gene], as.numeric(temp_gene@gs.vec))))
    })
    if(shiny) incProgress(amount=1/total_genes)
  }
  
  # G1
  # Calc adjusted p.val for period
  results@G1$BH.q.value.per1 <- p.adjust(as.numeric(results@G1$p.value.per1), method="BH")
  results@G1$BH.q.value.amp1 <- p.adjust(as.numeric(results@G1$p.value.amp1), method="BH")
  results@G1$BH.q.value.phase1 <- p.adjust(as.numeric(results@G1$p.value.phase1), method="BH")
  # G2
  results@G2$BH.q.value.per2 <- p.adjust(as.numeric(results@G2$p.value.per2), method="BH")
  results@G2$BH.q.value.amp2 <- p.adjust(as.numeric(results@G2$p.value.amp2), method="BH")
  results@G2$BH.q.value.phase2 <- p.adjust(as.numeric(results@G2$p.value.phase2), method="BH")
  # Comparison
  results@Comparison$BH.q.value.Per2_Per1 <- p.adjust(as.numeric(results@Comparison$p.value.Per2_Per1), method="BH")
  results@Comparison$BH.q.value.Amp2_Amp1 <- p.adjust(as.numeric(results@Comparison$p.value.Amp2_Amp1), method="BH")
  results@Comparison$BH.q.value.Ph2_Ph1 <- p.adjust(as.numeric(results@Comparison$p.value.Ph2_Ph1), method="BH")
  
  
  # if(shiny) incProgress(total_genes/total_genes)
  return(results)
}