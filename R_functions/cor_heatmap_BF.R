cor.BF <- function(x, y, iterations = 10000) {
  # Take two vectors as input: Compute the correlation coefficent and the corresponding BF
  # Return a named vector of length two: the Bayes factor and estimated coefficient rho
  library(BayesFactor)
  
  cor <- correlationBF(x, y)
  
  BF <- extractBF(cor)$bf
  coeff <- mean(posterior(cor, iterations = iterations, progress = FALSE)[, "rho"])
  
  out <- c(BF, coeff)
  names(out) <- c("BF", "rho")
  
  return(out)
}

# For example:
# cor.BF(rnorm(50), rnorm(50))

cor_heatmap_BF <- function(X, plot.all = FALSE, labels=NULL, low.col = "lightblue", high.col = "red", title=NULL, subtitle=NULL, digits=2, show.N=TRUE, show.BF=TRUE, legend.position = "right") {
  # X: a data.frame or matrix that can be passed to corr.test()
  #    data.frame is preferred because colnames will be retained for labels
  #  code largely taken from:
  #   http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
  # Load necessary libraries (will throw error if not installed):
  library(ggplot2)
  library(reshape2) # for melt()
  # library(psych) # for corr.test()
  
  if(!is.null(labels)) {
    if(length(labels) != ncol(X)) {
      warning(paste("Number of specified labels does not correspond to number of columns in X!", length(labels), "!=", ncol(X), "\nOriginal colnames will be retained and parameter `labels` will be ignored."))
    } else {
      colnames(X) <- labels
    }
  } 
  
  #Compute correlations:
  # cormat <- corr.test(X, adjust = p.adjust)
  cormat <- list(r = matrix(NA, ncol = ncol(X), nrow = ncol(X)),  # pre-allocate
                 bf = matrix(NA, ncol = ncol(X), nrow = ncol(X)),
                 n = matrix(NA, ncol = ncol(X), nrow = ncol(X))) 
  
  for(i in seq_along(X)) {
    for(j in seq_along(X)) {
      cor.bf <- cor.BF(as.matrix(X)[, i], as.matrix(X)[, j])
      cormat[["r"]][i, j] <- cor.bf["rho"]
      cormat[["bf"]][i, j] <- cor.bf["BF"]
      cormat[["n"]][i, j] <- nrow(na.omit(X[, c(i, j)]))
    }
  }
  
  cormat$r <- round(cormat$r, digits)
  
  # Remove redundant correlations if not otherwise specified:
  if(!plot.all) {
    cormat$r[upper.tri(cormat$r, diag = TRUE)] <- NA
    cormat$bf[upper.tri(cormat$bf, diag = TRUE)] <- NA
    cormat$n[upper.tri(cormat$n, diag = TRUE)] <- NA
  }
  
  colnames(cormat$r) <- colnames(cormat$bf) <- colnames(cormat$n) <- colnames(X)
  row.names(cormat$r) <- row.names(cormat$bf) <- row.names(cormat$n) <- colnames(X)
  
  # Put together data structure:
  cormat.melt <- melt(cormat$r, na.rm = TRUE, value.name = "r")
  cormat.melt <- merge(cormat.melt, melt(cormat$bf, value.name = "bf"))
  cormat.melt <- merge(cormat.melt, melt(cormat$n, value.name = "n"))
  
  if(show.BF) cormat.melt$bf.label <- prettyNum(round(cormat.melt$bf, 2), big.mark = ",", scientific = FALSE)
  
  corplot <- 
    ggplot(cormat.melt, aes(Var1, Var2, fill=r)) + 
    geom_tile(alpha=.8) + 
    geom_text(aes(label=format(r, nsmall=digits))) + 
    theme_minimal() + 
    coord_fixed() + 
    labs(title=title, subtitle=subtitle, x=NULL, y=NULL) + 
    theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust = 1), axis.text.y = element_text(angle = 45,  vjust = 1, hjust = 1), legend.position = legend.position) + 
    scale_fill_gradient2(low = low.col, high = high.col, mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Correlation")
  
  # Conditional additions:
  if(show.N) corplot <- corplot + geom_text(aes(label=paste("N =", n)), size=2, nudge_y = -0.25, color="darkgray") 
  if(show.BF) {
    corplot <- corplot +
      geom_text(aes(label="Bayes factor"), size=2, nudge_y = 0.35, color="darkgray", hjust="center") +
      geom_text(aes(label=bf.label), size=2, nudge_y = 0.25, color="darkgray", hjust="center")
  }
  
  return(corplot)
}