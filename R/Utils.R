'%!in%' <- function(x,y)!('%in%'(x,y))

#' supress_na_warning
#'
#' Supresses warning when BH q value is computed on datasets with NAs
#' @keywords NA, regression
#' @export
#' @examples
#' withCallingHandlers(BH <- p.adjust(as.numeric(results$p.value.per), method = "BH"), warning = supress_na_warning)
supress_na_warning <- function(w) if( any(grepl( "NAs introduced by coercion", w))) invokeRestart( "muffleWarning" )


#' best_fit_select
#'
#' Obtains the curve form from a series of fit results that fit best each feature.
#' @param results_list A list of datasets with different regression results
#' @keywords NA, regression
#' @export
#' @examples
#' results_list <- list.files(path="~/results/", pattern="*.csv")
#' best_fit_results <- do.call("best_fit_select", as.list(results_list))
best_fit_select <- function(results_list){
  args_list <- results_list
  n_features <- NA
  # Checks if the input is file or data
  for (i in 1:length(args_list)){
    a <- as.data.frame(args_list[i])
    colnames(a) <- sapply(strsplit(colnames(a), "_results."), "[[", 2)
    temp_name <- names(args_list[i])
    if(is.na(n_features)){
      n_features <- nrow(a)
    }else{
      if(n_features != nrow(a)){
        stop("The datasets provided have different number of rows")
      }
    }
    assign(temp_name, a, envir = .GlobalEnv)
  }
  # Once everything is loaded, we loop through the features searching for the best AIC
  best_fit_results <- data.frame()
  fit_data_colnames <- c("feature", "estimate.amp", "p.value.amp", "BH.q.value.amp"
                       , "estimate.phase", "p.value.phase", "BH.q.value.phase"
                       , "estimate.per", "p.value.per", "BH.q.value.per"
                       , "AIC", "BIC", "r", "best_curve")
  # global results data
  best_fit_results <- setNames(data.frame(matrix(ncol = length(fit_data_colnames), nrow = 0))
                               , fit_data_colnames
  )
  # print(n_features)
  for(f in 1:n_features){
    minimum_aic <- Inf
    best_adjust <- c(NA)
    feature_name <- NA
    for (i in 1:length(args_list)){
      temp_name <- names(args_list[i])
      temp <- get(temp_name)[f,]
      current_aic <- as.numeric(temp["AIC"])
      feature_name <- temp["feature"]
      if(!is.na(current_aic) & current_aic < minimum_aic){
        algorithm <- paste(unlist(strsplit(temp_name, "_"))[c(1:2)], collapse="_")
        assign("best_adjust", c(temp, "best_curve"=algorithm))
        minimum_aic <- current_aic
      }
      best_adjust <- unlist(best_adjust)
    }
    # Individual gene data
    df <- setNames(data.frame(matrix(ncol = length(fit_data_colnames), nrow = 0))
                   , fit_data_colnames
    )
    # Add the results from the best fit, if found
    if (all(is.na(best_adjust))){
      df <- base::rbind(df, c(feature_name, rep("NA", (length(fit_data_colnames)-1))))
      df <- setNames(df, fit_data_colnames)
      best_fit_results <- base::rbind(best_fit_results, df, stringsAsFactors=F, make.row.names = F)
      best_fit_results <- setNames(best_fit_results, fit_data_colnames)
    }else{
      df <- best_adjust[fit_data_colnames]
      df <- setNames(df, fit_data_colnames)
      df <- t(as.data.frame((df)))
      best_fit_results <- base::rbind(best_fit_results, df, stringsAsFactors=F, make.row.names = F)
    }
  }
  # Remove all datasets from envir
  for (i in 1:length(args_list)){
    temp_name <- names(args_list[i])
    rm(temp_name, envir =.GlobalEnv)
  }
  return(best_fit_results)
}



#' statistics_calc
#'
#' Obtains the number of tp, fp, fn and tn contained in a dataset
#' @param detected_genes the genes detected as positives
#' @param real_positive real circadian genes
#' @param real_negative real non-circadian genes
#' @keywords regression
#' @export
#' @examples
#' results_list <- list.files(path="~/results/", pattern="*.csv")
#' best_fit_results <- do.call("best_fit_select", as.list(results_list))
statistics_calc <- function(detected_genes, real_positive, real_negative){
  # Check how many genes overlap between predicted and real
  tp <-  Reduce(intersect, list(real_positive, detected_genes))
  # fp <- detected_genes[detected_genes %!in% real_positive]
  fp <- base::setdiff(detected_genes, real_positive)
  # fn <- real_positive[real_positive %!in% detected_genes]
  fn <- base::setdiff(real_positive, detected_genes)
  # tn <- real_negative[real_negative %!in% detected_genes]
  tn <- base::setdiff(real_negative, detected_genes)

  df <- data.frame(rbind(tp=length(tp), fp=length(fp),fn=length(fn), tn=length(tn)))
  colnames(df) <- "statistics"
  return(df)
}


#' data_summary
#'
#' Obtains the sd and mean grouping by column. Will return a df with the mean and sd by group of
#' the chosen variable.
#' from: STHDA barplot tutorial
#' @param data dataframe with data to summarize
#' @param varname column name for the variable which we want to summarize
#' @param groupnames grouping variables
#' @keywords sd, mean, summarize
#' @export
#' @examples
#' results_list <- data_summary(ToothGrowth, varname="len", groupnames=c("supp", "dose"))
#' summarized_data <- data_summary(all_data, varname = "meta2d_AMP", groupnames = "amp")
data_summary <- function(data, varname, groupnames){
  library("plyr")
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#' fisher_method
#'
#' Gets combined p.val by Fisher's method
#' @param x numeric vector with the p.values to merge
#' @keywords fisher p.value
#' @export
#' @examples
fisher_method <- function(x){
  if(!anyNA(x)){
    p <- pchisq(-2 * sum(log(x)), df=2*length(x), lower=FALSE)
  }else{
    p <- ifelse(all(is.na(x)), NA, x[!is.na(x)])
  }
  return(p)
}


#' get_fisher
#'
#' Uses fisher_method to calculate the fisher combined p.value of the provided numeric vectors
#' in a pairwise manner.
#' @param ... numeric vectors with the p.values to merge
#' @keywords fisher p.value
#' @export
#' @examples
#' circan_results$combined_pval <- get_fisher(as.numeric(circan_results$p.value.amp), as.numeric(circan_results$p.value.per))
get_fisher <- function(...){
  args_list <- as.data.frame(do.call(cbind, list(...)))
  ps <- c()
  for(i in 1:nrow(args_list)){
    ps <- c(ps, fisher_method(args_list[i,]))
  }
  return(ps)
}


#' meta2d_wrapper
#'
#' Wraps the meta 2d method with JTK + LS and returns a data frame with the results. The period is set between 20 and 28h.
#' @param data The data matrix with the omics data.
#' @param timepoints numeric vector containing the timepoints corresponding to each column in data.
#' @param cycMethod a character vector(length 1 or 2 or 3). User-defined methods for detecting rhythmic signals, must be
#' selected as any one, any two or all three methods(default) from "ARS"(ARSER), "JTK"(JTK_CYCLE) and "LS"(Lomb-Scargle).
#' @keywords metacycle, meta2d
#' @export
#' @examples
#' meta2d_results <- meta2d_wrapper(my_data, s2c$time)
meta2d_wrapper <- function(data, timepoints, min_per=20, max_per=28, cycMethod = c("LS", "JTK")){
  data_name <- deparse(substitute(data))
  write.csv(data, paste("./", data_name, ".csv", sep=""))
  meta_name <- deparse(substitute(meta_file))
  meta2d(infile=paste("./", data_name, ".csv", sep="")
         , outdir = "./"
         , filestyle="csv"
         , timepoints = as.numeric(timepoints)
         , minper = min_per
         , maxper = max_per
         , cycMethod = cycMethod
         , analysisStrategy = "auto"
         , outputFile = TRUE
         , outIntegration = "both"
         , adjustPhase = "predictedPer"
         , combinePvalue = "fisher"
         , weightedPerPha = FALSE
         , ARSmle = "auto"
         , ARSdefaultPer = 24
         , outRawData = FALSE
         , releaseNote = TRUE
         , outSymbol = "")
  # Move results
  file.remove(paste("./LSresult_",data_name,".csv", sep=""))
  meta2d_results <- read.csv(paste("./meta2d_",data_name,".csv", sep=""))
  file.remove(paste("./meta2d_",data_name,".csv", sep=""))
  return(meta2d_results)
}



#' jtk_wrapper
#'
#' Wraps JTK method. The period is set between 20 and 28h.
#' @param data The data matrix with the omics data.
#' @param s2c Data frame with at least the columns sample (which must coincide with the samples in data) and time (recording
#' time corresponding to each sample.
#' @keywords jtk
#' @export
#' @examples
#' jtk_results <- jtk_wrapper(my_data, s2c)
jtk_wrapper <- function(data, s2c, min_per=20, max_per=28){
  annot <- data.frame(Probeset=rownames(data))
  timepoints <- length(unique(s2c$time))
  nrep <- unique(table(s2c$time))
  lag <- unique(diff(as.numeric(unique(s2c$time))))
  data <- data[,s2c$sample]

  jtkdist(timepoints, nrep) # total time points, # replicates per time point

  periods <- ceiling(min_per/lag):round(max_per/lag) # number of time points per cycle. (10/6=; 20/6)
  # cat(timepoints, nrep, lag, periods)
  jtk.init(periods, lag)  # 4 is the number of hours between time points

  res <- apply(data,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results_jtk <- cbind(annot,res,data)
  results_jtk <- results_jtk[order(res$ADJ.P,-res$AMP),]
  return(results_jtk)
}


cbind.fill <- function(...){
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function (x)
    rbind(x, matrix(, n-nrow(x), ncol(x)))))
}
