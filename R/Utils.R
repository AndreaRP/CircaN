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
  fp <- detected_genes[detected_genes %!in% real_positive]
  fn <- real_positive[real_positive %!in% detected_genes]
  tn <- real_negative[real_negative %!in% detected_genes]
  
  df <- data.frame(rbind(tp=length(tp), fp=length(fp),fn=length(fn), tn=length(tn)))
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

# best_fit_select <- function(...){
#   args_list <- list(...)
#   print("#############")
#   print(length(args_list))
#   print(class(args_list))
#   print("#############")
#   n_features <- NA
#   # Checks if the input is file or data
#   for (a in args_list){
#     print("#############")
#     print(names(a))
#     print("#############")
#     if(grepl(".csv", a)){
#       temp_name <- sapply(strsplit(a, "\\."), "[[", 1)
#       # temp <- read.csv(paste("./data/results/", a, sep=""), stringsAsFactors = F)[1:20,]#######
#     }else{
#       temp_name <- a
#       temp <- get(a)
#     }
#     if(is.na(n_features)){
#       n_features <- nrow(temp)
#     }else{
#       if(n_features != nrow(temp)){
#         stop("The datasets provided have different number of rows")
#       }
#     }
#     assign(temp_name, temp, envir = .GlobalEnv)
#   }
#   # Once everything is loaded, we loop through the features searching for the best AIC
#   best_fit_results <- data.frame()
#   fit_data_colnames <- c("feature", "estimate.amp", "p.value.amp", "BH.q.value.amp"
#                          , "estimate.phase", "p.value.phase", "BH.q.value.phase"
#                          , "estimate.per", "p.value.per", "BH.q.value.per"
#                          , "AIC", "BIC", "r", "best_curve")
#   # global results data
#   best_fit_results <- setNames(data.frame(matrix(ncol = length(fit_data_colnames), nrow = 0))
#                       , fit_data_colnames
#   )
#   # print(n_features)
#   for(f in 1:n_features){
#     # print(paste("feature ", f)) ####
#     minimum_aic <- Inf
#     best_adjust <- c(NA)
#     feature_name <- NA
#     # feat_name <- (NA)########
#     for (a in args_list){
#       # print(a) ####
#       temp_name <- sapply(strsplit(a, "\\."), "[[", 1)
#       temp <- get(temp_name)[f,]
#       current_aic <- as.numeric(temp["AIC"])
#       # if(is.na(feat_name) & as.character(temp["feature"])!=feat_name) feat_name <- as.character(temp["feature"])########
#       feature_name <- temp["feature"]
#       if(!is.na(current_aic) & current_aic < minimum_aic){
#         # algorithm <- paste(rev(rev(unlist(strsplit(a, "_")))[-c(1:3)]), collapse="_")
#         algorithm <- paste(unlist(strsplit(a, "_"))[c(1:2)], collapse="_")
#         assign("best_adjust", c(temp, "best_curve"=algorithm))
#         minimum_aic <- current_aic
#       }
#       best_adjust <- unlist(best_adjust)
#     }
#     # Individual gene data
#     df <- setNames(data.frame(matrix(ncol = length(fit_data_colnames), nrow = 0))
#                    , fit_data_colnames
#     )
#     # print(paste("feature ", feature_name))
#     # Add the results from the best fit, if found
#     if (all(is.na(best_adjust))){
#       df <- base::rbind(df, c(feature_name, rep("NA", (length(fit_data_colnames)-1))))
#       df <- setNames(df, fit_data_colnames)
#       best_fit_results <- base::rbind(best_fit_results, df, stringsAsFactors=F, make.row.names = F)
#       best_fit_results <- setNames(best_fit_results, fit_data_colnames)
#       # print(best_fit_results)
#     }else{
#       df <- best_adjust[fit_data_colnames]
#       # print(df)#######
#       df <- setNames(df, fit_data_colnames)
#       df <- t(as.data.frame((df)))
#       best_fit_results <- base::rbind(best_fit_results, df, stringsAsFactors=F, make.row.names = F)
#       # print("best_fit_results Not NA") ####
#       # print(data.frame(best_fit_results), na.print = "NA2") ####
#     }
#   }
#   # Remove all datasets from envir
#   for (a in args_list){
#     temp_name <- sapply(strsplit(a, "\\."), "[[", 1)
#     rm(temp_name, envir =.GlobalEnv)
#   }
#   return(best_fit_results)
# }
# 


