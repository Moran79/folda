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

  LD1 <- LD2 <- Y <- xText <- yPred <- density <- NULL

  predictInternal <- function(x, LDscores){ # predict class based on LDscores
    x$groupMeans <- x$groupMeans[, seq_len(ncol(LDscores)), drop = FALSE]
    LDscores <- as.matrix(LDscores)
    loglikelihood <- LDscores %*% t(x$groupMeans) + matrix(log(x$prior) - 0.5 * rowSums(x$groupMeans^2), nrow(LDscores), length(x$prior), byrow = TRUE)
    likelihood <- exp(loglikelihood - apply(loglikelihood, 1, max))
    posterior <- likelihood / apply(likelihood, 1, sum)
    return(rownames(x$groupMeans)[max.col(-posterior %*% t(x$misClassCost), ties.method = "first")])
  }

  getBoundary1D <- function(prior, groupMeans, xRange, tol = 1e-10){
    levelName <- names(prior); prior <- as.vector(prior)
    groupMeans <- as.vector(groupMeans)
    targetNum <- log(prior) - 0.5 * groupMeans^2
    cutoffCandidates <- c(); cutOffFake <- c()

    for(i in seq_len(length(prior) - 1)){
      for(j in seq(i+1, length(prior), by = 1)){
        cutoffCandidates <- c(cutoffCandidates, -diff(targetNum[c(i,j)]) / diff(groupMeans[c(i,j)]))
      }
    }
    cutoffCandidates <- sort(stats::na.omit(cutoffCandidates))
    cutoffCandidates <- cutoffCandidates[cutoffCandidates >= xRange[1] & cutoffCandidates <= xRange[2]]
    yPred <- c(); leftClass <- rightClass <- levelName[which.max(prior)]

    for(cutoff in cutoffCandidates){
      leftClass <- levelName[which.max((cutoff - tol) * groupMeans + targetNum)]
      rightClass <- levelName[which.max((cutoff + tol) * groupMeans + targetNum)]
      if(leftClass != rightClass){
        yPred <- c(yPred, leftClass)
      }else cutOffFake <- c(cutOffFake, cutoff)
    }
    yPred <- c(yPred, rightClass)
    cutoffCandidates <- setdiff(cutoffCandidates, cutOffFake)
    xText <- c(xRange[1], cutoffCandidates, xRange[2])
    xText <- diff(xText) / 2 + xText[seq_len(length(xText) - 1)] # get mid point
    return(list(positionX = cutoffCandidates,
                textInfo = data.frame(xText = xText, yPred = yPred)))
  }

  if(missing(datX) || missing(response)) stop("Please input the training X and Y for plotting")
  response <- droplevels(as.factor(response))
  if(!identical(x$misClassCost, checkPriorAndMisClassCost(NULL, NULL, response)$misClassCost))
    warning("With customized misClassCost, plot may not reflect real decision boundaries")
  LDscores <- round(getLDscores(modelULDA = x, data = datX, nScores = min(2, ncol(x$scaling))), 10) # increased stability
  datPlot <- cbind.data.frame(response = response, LDscores)
  colorManual <- grDevices::hcl.colors(nlevels(response)); names(colorManual) <- levels(response)

  if(dim(x$scaling)[2] == 1){ # Only one LD is available, draw the histogram
    estimatedPrior <- table(response, dnn = NULL) / length(response)
    datPlot <- rbind(datPlot, datPlot) # In case some groups only have one observation
    datPlot <- do.call(rbind, lapply(seq_along(estimatedPrior), function(i) cbind(with(stats::density(datPlot$LD1[datPlot$response == names(estimatedPrior)[i]]), data.frame(LD1 = x, density = y * estimatedPrior[i])), response = names(estimatedPrior)[i])))
    datBounds <- getBoundary1D(prior = x$prior, groupMeans = x$groupMeans, xRange = range(LDscores))
    datBounds$textInfo$density <- max(datPlot$density) # write the y position for text
    p <- ggplot2::ggplot(data = datPlot)+
      ggplot2::geom_ribbon(ggplot2::aes(x = LD1, ymin = 0, ymax = density, fill = response), alpha = 0.8)+
      ggplot2::geom_vline(xintercept = datBounds$positionX, color = "black", linetype = "dotted")+
      ggplot2::geom_text(data = datBounds$textInfo, ggplot2::aes(x = xText, y = density, label = yPred, color = yPred), show.legend = FALSE)+
      ggplot2::scale_fill_manual(values = colorManual)+
      ggplot2::scale_color_manual(values = colorManual)+
      ggplot2::labs(title = "Density plot of the first linear discriminant score", subtitle = "Dotted line is the decision boundary")+
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
      ggplot2::scale_color_manual(values = colorManual)+
      ggplot2::scale_fill_manual(values = colorManual)+
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


