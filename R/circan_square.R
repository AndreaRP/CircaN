arccot <- function(x){
  y = atan(1/x)
  return(y)
}

circan_square <- function(data, s2c, shiny = FALSE, mode = "default", init_value = 24,
                       max_per = Inf, min_per = -Inf) {
  #############################
  # data <- data_file
  # s2c <- meta_file
  # mode="port"
  # init_value = 24
  # max_per = 28
  # min_per = 20
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
    # gene <- 8
    #############################

    # cat(gene, "\t")
    data.st <- (data[gene, ] - mean(data[gene, ]))/sqrt(var(data[gene,]))
    df0 <- merge(data.st, s2c, by.x = 0, by.y = "sample")
    df0 <- dplyr::rename(df0, data = x)
    df <- as.data.frame(df0[, c("data", "ind", "time")])
    df$data <- as.numeric(df$data)
    df$time <- as.numeric(df$time)
    df$ind <- as.numeric(as.character(df$ind))
    df <- df[order(df$time),]
    gd <- nlme::groupedData(data ~ time | ind, data = df)
    # Find likely start for amp
    amp_init <- abs(max(gd$data) - min(gd$data))/2

    result = tryCatch({
      nls.model = nls(data ~ atan(sin(2*pi*time/per + phase))+arccot(sin(2*pi*time/per + phase))
                    , start = list(amp = amp_init, phase = 0, per = init_value)
                    , lower = list(amp = 0, phase = 0, per = min_per)
                    , upper = list(amp = Inf, phase = Inf, per = max_per)
                    , algorithm = mode
                    , data = gd
                    , control = list(maxiter = 200)
                    )
      stats <- as.data.frame(broom::tidy(nls.model))
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
          nls.model = nls(data ~ (atan(sin(2*pi*time/init_value + phase))+arccot(sin(2*pi*time/init_value + phase)))
                          , start = list(amp = amp_init, phase = 0)
                          , lower = list(amp=0, phase=0)
                          , upper = list(amp=Inf, phase=Inf)
                          , algorithm = mode
                          , data = gd
                          , control = list(maxiter = 200))
          stats <- as.data.frame(broom::tidy(nls.model))
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