#' Plot Decision Boundaries and Linear Discriminant Scores
#'
#' This function plots the decision boundaries and linear discriminant (LD)
#' scores for a given ULDA model. If it is a binary classification problem, a
#' density plot is created. Otherwise, a scatter plot with decision boundaries
#' is generated.
#'
#' @param x A fitted ULDA model object.
#' @param datX A data frame containing the predictor variables.
#' @param response A factor representing the response variable (training labels)
#'   corresponding to `datX`.
#' @param ... Additional arguments.
#'
#' @return A `ggplot2` plot object, either a density plot or a scatter plot with
#'   decision boundaries.
#'
#' @export
#'
#' @examples
#' fit <- folda(datX = iris[, -5], response = iris[, 5], subsetMethod = "all")
#' plot(fit, iris[, -5], iris[, 5])

plot.ULDA <- function(x, datX, response, ...){

  LD1 <- LD2 <- Y <- density <- NULL

  predictInternal <- function(x, LDscores){ # predict class based on LDscores
    x$groupMeans <- x$groupMeans[, seq_len(ncol(LDscores)), drop = FALSE]
    LDscores <- as.matrix(LDscores)
    loglikelihood <- LDscores %*% t(x$groupMeans) + matrix(log(x$prior) - 0.5 * rowSums(x$groupMeans^2), nrow(LDscores), length(x$prior), byrow = TRUE)
    likelihood <- exp(loglikelihood - apply(loglikelihood, 1, max))
    posterior <- likelihood / apply(likelihood, 1, sum)
    return(rownames(x$groupMeans)[max.col(-posterior %*% t(x$misClassCost), ties.method = "first")])
  }

  if(missing(datX) || missing(response)) stop("Please input the training X and Y for plotting")
  response <- droplevels(as.factor(response))
  if(!identical(x$misClassCost, checkPriorAndMisClassCost(NULL, NULL, response)$misClassCost))
    warning("With customized misClassCost, plot may not reflect real decision boundaries")
  LDscores <- getLDscores(modelLDA = x, data = datX, nScores = min(2, ncol(x$scaling)))
  datPlot <- cbind.data.frame(response = response, LDscores)

  if(dim(x$scaling)[2] == 1){ # Only one LD is available, draw the histogram
    estimatedPrior <- table(response, dnn = NULL) / length(response)
    datPlot <- do.call(rbind, lapply(seq_along(estimatedPrior), function(i) cbind(with(stats::density(datPlot$LD1[response == names(estimatedPrior)[i]]), data.frame(LD1 = x, density = y * estimatedPrior[i])), response = names(estimatedPrior)[i])))
    positionX <- as.vector(-diff(log(x$prior) - 0.5 * rowSums(x$groupMeans^2)) / diff(x$groupMeans))
    positionY <- max(datPlot$density) - 0.1
    p <- ggplot2::ggplot(data = datPlot)+
      ggplot2::geom_ribbon(ggplot2::aes(x = LD1, ymin = 0, ymax = density, fill = response), alpha = 0.8)+
      ggplot2::geom_vline(xintercept = positionX, color = grDevices::hcl.colors(3)[2])+
      ggplot2::geom_label(x = positionX, y = positionY, label = "Decision\nBoundary")+
      ggplot2::scale_fill_manual(values = grDevices::hcl.colors(3)[c(1,3)])+
      ggplot2::labs(title = "Density plot of the first linear discriminant score")+
      ggplot2::theme_bw()
  } else{
    LDranges <- apply(LDscores, 2, range)
    fakeLDs <- expand.grid(list(LD1 = seq(from = LDranges[1,1], to = LDranges[2,1], length.out = 100),
                                LD2 = seq(from = LDranges[1,2], to = LDranges[2,2], length.out = 100)))
    fakeLDs$Y <- factor(predictInternal(x = x, LDscores = fakeLDs), levels = levels(response))
    p <- ggplot2::ggplot(data = datPlot)+
      ggplot2::geom_raster(data = fakeLDs, ggplot2::aes(x = LD1, y = LD2, fill = Y), alpha = 0.4, show.legend = FALSE)+
      ggplot2::geom_point(ggplot2::aes(x = LD1, y = LD2, color = response), alpha = 0.5)+
      ggplot2::geom_contour(data = fakeLDs, ggplot2::aes(x = LD1, y = LD2, z = as.integer(Y)), color = "white", linetype = "dotted")+
      ggplot2::scale_color_manual(values = grDevices::hcl.colors(nlevels(response)))+
      ggplot2::scale_fill_manual(values = grDevices::hcl.colors(nlevels(response)))+
      ggplot2::labs(title = "Scatter plot of the first two linear discriminant scores")+
      ggplot2::theme_bw()
  }
  return(p)
}


#' @export
print.ULDA <- function(x, ...) {
  cat("\nOverall Pillai's trace: ", format(x$statPillai, digits = 4), "\n", sep = "")
  cat("Associated p-value: ", format(x$pValue, digits = 4), "\n\n", sep = "")

  cat("Prediction Results on Training Data:\n")
  accuracy <- sum(diag(x$confusionMatrix)) / sum(x$confusionMatrix)
  cat("Refitting Accuracy: ", format(accuracy, digits = 4), "\n", sep = "")
  cat("Gini Index: ", format(x$predGini, digits = 4), "\n", sep = "")
  cat("\nConfusion Matrix:\n")
  print(x$confusionMatrix)

  cat("\nGroup means of LD scores:\n")
  print(x$groupMeans)

  if (!is.null(x$forwardInfo)) {
    cat("\nForward Selection Results:\n")
    print(x$forwardInfo)
  }

  invisible(x)
}


