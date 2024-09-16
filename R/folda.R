#' Forward Uncorrelated Linear Discriminant Analysis
#'
#' This function fits a ULDA (Uncorrelated Linear Discriminant Analysis) model
#' to the provided data, with an option for forward selection of variables based
#' on Pillai's trace or Wilks' Lambda. It can also handle missing values,
#' perform downsampling, and compute the linear discriminant scores and group
#' means for classification. The function returns a fitted ULDA model object.
#'
#' @param datX A data frame of predictor variables.
#' @param response A factor representing the response variable with multiple
#'   classes.
#' @param subsetMethod A character string specifying the method for variable
#'   selection. Options are `"forward"` for forward selection or `"all"` for
#'   using all variables. Default is `"forward"`.
#' @param testStat A character string specifying the test statistic to use for
#'   forward selection. Options are `"Pillai"` or `"Wilks"`. Default is
#'   `"Pillai"`.
#' @param correction A logical value indicating whether to apply a multiple
#'   comparison correction during forward selection. Default is `TRUE`.
#' @param alpha A numeric value between 0 and 1 specifying the significance
#'   level for the test statistic during forward selection. Default is 0.1.
#' @param prior A numeric vector representing the prior probabilities for each
#'   class in the response variable. If `NULL`, the observed class frequencies
#'   are used as the prior. Default is `NULL`.
#' @param misClassCost A square matrix \eqn{C}, where each element \eqn{C_{ij}}
#'   represents the cost of classifying an observation into class \eqn{i} given
#'   that it truly belongs to class \eqn{j}. If `NULL`, a default matrix with
#'   equal misclassification costs for all class pairs is used. Default is
#'   `NULL`.
#' @param missingMethod A character vector of length 2 specifying how to handle
#'   missing values for numerical and categorical variables, respectively.
#'   Default is `c("medianFlag", "newLevel")`.
#' @param downSampling A logical value indicating whether to perform
#'   downsampling to balance the class distribution in the training data or to
#'   improve computational efficiency. Default is `FALSE`. Note that if
#'   downsampling is applied and the `prior` is `NULL`, the class prior will be
#'   calculated based on the downsampled data. To retain the original prior,
#'   please specify it explicitly using the `prior` parameter.
#' @param kSample An integer specifying the maximum number of samples to take
#'   from each class during downsampling. If `NULL`, the number of samples is
#'   limited to the size of the smallest class. Default is `NULL`.
#'
#' @return A list of class `ULDA` containing the following components:
#'   \item{scaling}{The matrix of scaling coefficients for the linear
#'   discriminants.} \item{groupMeans}{The group means of the linear
#'   discriminant scores.} \item{prior}{The prior probabilities for each class.}
#'   \item{misClassCost}{The misclassification cost matrix.}
#'   \item{misReference}{A reference for handling missing values.}
#'   \item{terms}{The terms used in the model formula.} \item{xlevels}{The
#'   levels of the factors used in the model.} \item{varIdx}{The indices of the
#'   selected variables.} \item{varSD}{The standard deviations of the selected
#'   variables.} \item{varCenter}{The means of the selected variables.}
#'   \item{statPillai}{The Pillai's trace statistic.} \item{pValue}{The p-value
#'   associated with Pillai's trace.} \item{predGini}{The Gini index of the
#'   predictions on the training data.} \item{confusionMatrix}{The confusion
#'   matrix for the training data predictions.} \item{forwardInfo}{Information
#'   about the forward selection process, if applicable.} \item{stopInfo}{A
#'   message indicating why forward selection stopped, if applicable.}
#'
#' @export
#'
#' @references Howland, P., Jeon, M., & Park, H. (2003). \emph{Structure
#'   preserving dimension reduction for clustered text data based on the
#'   generalized singular value decomposition}. SIAM Journal on Matrix Analysis
#'   and Applications
#'
#'   Wang, S. (2024). A New Forward Discriminant Analysis Framework Based On
#'   Pillai's Trace and ULDA. \emph{arXiv preprint arXiv:2409.03136}. Available
#'   at \url{https://arxiv.org/abs/2409.03136}.
#'
#' @examples
#' # Fit the ULDA model
#' fit <- folda(datX = iris[, -5], response = iris[, 5], subsetMethod = "all")
#'
#' # Fit the ULDA model with forward selection
#' fit <- folda(datX = iris[, -5], response = iris[, 5], subsetMethod = "forward")
folda <- function(datX,
                  response,
                  subsetMethod = c("forward", "all"),
                  testStat = c("Pillai", "Wilks"),
                  correction = TRUE,
                  alpha = 0.1,
                  prior = NULL,
                  misClassCost = NULL,
                  missingMethod = c("medianFlag", "newLevel"),
                  downSampling = FALSE,
                  kSample = NULL){

  # Pre-processing: Arguments ----------------------------------

  datX <- data.frame(datX) # change to data.frame, remove the potential tibble attribute
  for(i in seq_along(datX)){ # remove ordered factors
    if(inherits(datX[,i], c("ordered"))) class(datX[,i]) <- "factor"
  }
  response <- droplevels(factor(response, ordered = FALSE)) # remove ordered factors
  subsetMethod <- match.arg(subsetMethod, c("forward", "all"))

  # Pre-processing: Data Cleaning -----------------------------------------------

  idxTrain <- getDownSampleInd(response = response,
                               downSampling = downSampling,
                               kSample = kSample)
  response <- droplevels(response[idxTrain])
  priorAndMisClassCost <- checkPriorAndMisClassCost(prior = prior, misClassCost = misClassCost, response = response)

  imputedSummary <- missingFix(data = datX[idxTrain, , drop = FALSE], missingMethod = missingMethod)
  datX <- imputedSummary$data # this step also removes some constant columns
  if(any(dim(datX) == 0)) stop("No available observations or features, which maybe due to preprocessing steps.")

  # Pre-processing: Data Transformation -----------------------------------------------

  modelFrame <- stats::model.frame(formula = ~.-1, datX, na.action = "na.fail")
  Terms <- stats::terms(modelFrame)
  m <- scale(stats::model.matrix(Terms, modelFrame)) # constant cols would be changed to NaN in this step
  currentVarList <- seq_len(ncol(m))

  # Forward Selection -----------------------------------------------

  if(subsetMethod == "forward"){
    forwardRes <- forwardSel(m = m,
                             response = response,
                             testStat = testStat,
                             alpha = alpha,
                             correction = correction)

    # When no variable is selected, use the full model
    if(length(forwardRes$currentVarList) != 0){
      # modify the design matrix to make it more compact
      selectedVarRawIdx <- unique(sort(attributes(m)$assign[forwardRes$currentVarList]))
      modelFrame <- stats::model.frame(formula = ~.-1, datX[, selectedVarRawIdx, drop = FALSE], na.action = "na.fail")
      Terms <- stats::terms(modelFrame)
      m <- scale(stats::model.matrix(Terms, modelFrame))
      currentVarList <- which(colnames(m) %in% forwardRes$forwardInfo$var)
    }
  }

  varSD <- attr(m,"scaled:scale")[currentVarList]
  varCenter <- attr(m,"scaled:center")[currentVarList]
  m <- m[, currentVarList, drop = FALSE]

  # ULDA -----------------------------------------------

  # Step 1: SVD on the combined matrix H
  groupMeans <- tapply(c(m), list(rep(response, dim(m)[2]), col(m)), function(x) mean(x, na.rm = TRUE))
  Hb <- sqrt(tabulate(response)) * groupMeans
  Hw <- m - groupMeans[response, , drop = FALSE]
  if(diff(dim(m)) < 0){ # More rows than columns
    qrRes <- qrEigen(Hw)
    fitSVD <- saferSVD(rbind(Hb, qrRes$R))
  }else fitSVD <- saferSVD(rbind(Hb, Hw))

  # Step 2: SVD on the P matrix
  N <- nrow(m); J <- nlevels(response)
  rankT <- sum(fitSVD$d >= max(dim(fitSVD$u), dim(fitSVD$v)) * .Machine$double.eps * fitSVD$d[1])
  fitSVDp <- saferSVD(fitSVD$u[seq_len(J), seq_len(rankT), drop = FALSE], uFlag = FALSE)
  rankAll <- min(J - 1, sum(fitSVDp$d >= max(J, rankT) * .Machine$double.eps * fitSVDp$d[1]))

  # Step 3: Transform Sw into identity matrix
  unitSD <- diag(sqrt((N - J) / abs(1 - fitSVDp$d^2 + 1e-5)), nrow = rankAll) # Scale to unit var
  scalingFinal <- (fitSVD$v[, seq_len(rankT), drop = FALSE] %*% diag(1 / fitSVD$d[seq_len(rankT)], nrow = rankT) %*% fitSVDp$v[, seq_len(rankAll), drop = FALSE]) %*% unitSD
  rownames(scalingFinal) <- colnames(m)
  groupMeans <- groupMeans %*% scalingFinal
  rownames(groupMeans) <- levels(response)
  colnames(groupMeans) <- colnames(scalingFinal) <- paste("LD", seq_len(ncol(groupMeans)), sep = "")

  # Summary and outputs -----------------------------------------------

  statPillai <- sum(fitSVDp$d[seq_len(rankAll)]^2)
  p <- rankT; s <- rankAll; numF <- N - J - p + s; denF <- abs(p - J + 1) + s
  pValue <- ifelse(numF > 0, stats::pbeta(1 - statPillai / s, shape1 = numF * s / 2, shape2 = denF * s / 2), 0)

  res <- list(scaling = scalingFinal, groupMeans = groupMeans, prior = priorAndMisClassCost$prior,
              misClassCost = priorAndMisClassCost$misClassCost, misReference = imputedSummary$ref,
              terms = Terms, xlevels = stats::.getXlevels(Terms, modelFrame), varIdx = currentVarList,
              varSD = varSD, varCenter = varCenter, statPillai = statPillai, pValue = pValue)

  if(subsetMethod == "forward"){
    res$forwardInfo = forwardRes$forwardInfo
    res$stopInfo <- forwardRes$stopInfo
  }

  class(res) <- "ULDA"
  pred <- factor(stats::predict(res, datX), levels = levels(response))
  res$predGini <- 1 - sum(unname(table(pred) / dim(datX)[1])^2)
  res$confusionMatrix <- table(Predicted = pred, Actual = response)
  return(res)
}


#' Predict Method for ULDA Model
#'
#' This function predicts the class labels or class probabilities for new data
#' using a fitted ULDA model. The prediction can return either the most likely
#' class (`"response"`) or the posterior probabilities for each class
#' (`"prob"`).
#'
#' @param object A fitted `ULDA` model object.
#' @param newdata A data frame containing the new predictor variables for which
#'   predictions are to be made.
#' @param type A character string specifying the type of prediction to return.
#'   `"response"` returns the predicted class labels, while `"prob"` returns the
#'   posterior probabilities for each class. Default is `"response"`.
#' @param ... Additional arguments.
#'
#' @return If `type = "response"`, the function returns a vector of predicted
#'   class labels. If `type = "prob"`, it returns a matrix of posterior
#'   probabilities, where each row corresponds to a sample and each column to a
#'   class.
#'
#' @export
#'
#' @examples
#' fit <- folda(datX = iris[, -5], response = iris[, 5], subsetMethod = "all")
#'
#' # Predict class labels
#' predictions <- predict(fit, iris, type = "response")
#'
#' # Predict class probabilities
#' prob_predictions <- predict(fit, iris, type = "prob")
predict.ULDA <- function(object, newdata, type = c("response", "prob"), ...){
  type <- match.arg(type, c("response", "prob"))
  LDscores <- getLDscores(modelULDA = object, data = newdata)
  loglikelihood <- LDscores %*% t(object$groupMeans) + matrix(log(object$prior) - 0.5 * rowSums(object$groupMeans^2), nrow(LDscores), length(object$prior), byrow = TRUE)
  likelihood <- exp(loglikelihood - apply(loglikelihood, 1, max))
  posterior <- likelihood / apply(likelihood, 1, sum)
  if(type == "prob") return(posterior)
  return(rownames(object$groupMeans)[max.col(-posterior %*% t(object$misClassCost), ties.method = "first")])
}


