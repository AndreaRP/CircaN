#' circan_damp
#'
#' Fit to a damped sine curve in the frame of a nonlinear least squares model for accurate detection of circadian expression patterns.
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
#' circan_damp(data=expression_example, s2c=metadata_example, mode="port")
circan_damp <- function(data, s2c, shiny = FALSE, mode = "default", init_value = 24,
                        max_per = Inf, min_per = -Inf) {
  #############################
  # data <- subsampled_exprs_nls
  # s2c <- subsampled_s2c
  # mode="port"
  # init_value = 24
  # max_per = Inf
  # min_per = -Inf
  # shiny = FALSE
  #############################

  s2c$time <- as.numeric(as.character(s2c$time))
  s2c$ind <- as.integer(as.character(s2c$ind))
  s2c$sample <- as.character(s2c$sample)
  t <- unique(as.numeric(as.character(s2c$time)))
  data0 <- data
  total_genes <- nrow(data)
  rownames(data) <- as.character(data[, 1])
  data <- data[, -1]
  data <- as.matrix(data)
  results.cols <- c("feature", "estimate.amp", "std.error.amp",
                    "statistic.amp", "p.value.amp", "estimate.phase", "std.error.phase",
                    "statistic.phase", "p.value.phase", "estimate.per", "std.error.per",
                    "statistic.per", "p.value.per", "AIC", "BIC", "r")
  results <- setNames(data.frame(matrix(ncol = 16, nrow = 0)), results.cols)

  for (gene in 1:nrow(data)) {
    #############################
    # gene <- 2
    #############################

    # cat(gene, "\t")
    data.st <- (data[gene, ] - mean(data[gene, ]))/sqrt(var(data[gene,]))
    df0 <- merge(data.st, s2c, by.x = 0, by.y = "sample")
    df0 <- dplyr::rename(df0, data = x)
    df <- as.data.frame(df0[, c("data", "ind", "time")])
    df$data <- as.numeric(df$data)
    df$time <- as.numeric(df$time)
    gd <- nlme::groupedData(data ~ time | ind, data = df)
    # Find likely start for amp
    amp_init <- abs(max(gd$data) - min(gd$data))/2

    result = tryCatch({
      nls.model = nls(data ~ amp * cos(2*pi/per * time - phase*(2*pi/per)) * exp(-damp*time)
                      , start = list(amp = amp_init, phase = 0, per = init_value, damp=0)
                      , algorithm = mode
                      , lower = list(amp=0, phase=0, per = min_per, damp=0)
                      , upper = list(amp=Inf, phase=Inf, per = max_per,damp=Inf)
                      , data = gd
                      , control = list(maxiter = 200)
      )
      stats <- as.data.frame(broom::tidy(nls.model))
      stats <- stats[which(stats$term %in% c("amp", "phase", "per")),]
      rownames(stats) <- stats[, 1]
      stats <- stats[, -1]
      vec <- as.numeric(c(t(stats)))
      names(vec) <- c(outer(colnames(stats), rownames(stats),
                            paste, sep = "."))
      # Un-standarize amplitude estimation
      vec["estimate.amp"] <- (as.numeric(vec["estimate.amp"])*sqrt(var(data[gene,])))
      aic <- AIC(nls.model)
      bic <- BIC(nls.model)
      akaike <- c(aic, bic)
      names(akaike) <- c("AIC", "BIC")
      r <- cor(gd$data, predict(nls.model))
    }, error = function(e1) {
      convergence_error1 <- grepl("Convergence failure|the inverse cannot be computed",
                                  as.character(e1))
      if (convergence_error1) {
        tryCatch({ # Fix period
          nls.model = nls(data ~ amp * cos(2*pi/init_value * time - phase*(2*pi/init_value)) * exp(-damp*time)
                          , start = list(amp = amp_init, phase = 0, damp = 0)
                          , lower = list(amp=0, phase=0, damp=0)
                          , upper = list(amp=Inf, phase=Inf, damp=Inf)
                          , algorithm = mode
                          , data = gd
                          , control = list(maxiter = 200)
          )
          stats <- as.data.frame(broom::tidy(nls.model))
          stats <- stats[which(stats$term %in% c("amp", "phase")),]
          rownames(stats) <- stats[, 1]
          stats <- stats[, -1]
          vec <- as.numeric(c(t(stats)))
          names(vec) <- c(outer(colnames(stats), rownames(stats),
                                paste, sep = "."))
          # Un-standarize amplitude estimation
          vec["estimate.amp"] <- (as.numeric(vec["estimate.amp"])*sqrt(var(data[gene,])))
          temp <- c(init_value, "NA", "NA", "NA")
          names(temp) <- c("estimate.period", "std.error.period",
                           "statistic.period", "p.value.period")
          vec <<- c(vec, temp)

          aic <- AIC(nls.model)
          bic <- BIC(nls.model)
          akaike <<- c(aic, bic)
          names(akaike) <- c("AIC", "BIC")

          r <<- cor(gd$data, predict(nls.model))
        }, error = function(e2) {
          convergence_error2 <- grepl("Convergence failure|the inverse cannot be computed",
                                      as.character(e2))
          if (convergence_error2) {
            vec <<- rep("NA", times = ncol(results) - 4)

            akaike <<- c("NA", "NA")

            r <<- "NA"
          }
        })
      }
    }, finally = {
      res <- rbind(c(feature = rownames(data)[gene], vec, akaike, r = r))
      results[gene, ] <- res
    })
    if (shiny)
      incProgress(amount = 1/total_genes)
  }
  withCallingHandlers(results$BH.q.value.per <- p.adjust(as.numeric(results$p.value.per), method = "BH")
                      , warning = supress_na_warning)
  withCallingHandlers(results$BH.q.value.amp <- p.adjust(as.numeric(results$p.value.amp), method = "BH")
                      , warning = supress_na_warning)
  withCallingHandlers(results$BH.q.value.phase <- p.adjust(as.numeric(results$p.value.phase), method = "BH")
                      , warning = supress_na_warning)
  return(results)
}
