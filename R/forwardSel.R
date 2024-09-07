#' Forward Selection via Multivariate Test Statistics
#'
#' This function performs forward selection on a dataset based on multivariate
#' test statistics (`Pillai` or `Wilks`). It iteratively adds variables that
#' contribute the most to the test statistic until no significant variables are
#' found or a stopping criterion is met.
#'
#' @param m A numeric matrix containing the predictor variables. Rows represent
#'   observations and columns represent variables.
#' @param response A factor representing the response variable with multiple
#'   levels (groups).
#' @param testStat A character string specifying the test statistic to use. Can
#'   be `"Pillai"` (default) or `"Wilks"`.
#' @param alpha A numeric value between 0 and 1 specifying the significance
#'   level for the test. Default is 0.1.
#' @param correction A logical value indicating whether to apply a multiple
#'   comparison correction. Default is `TRUE`.
#'
#' @return A list with three components: \item{currentVarList}{A vector of
#'   selected variable indices based on the forward selection process.}
#'   \item{forwardInfo}{A data frame containing detailed information about the
#'   forward selection process, including the selected variables, test
#'   statistics, and thresholds.} \item{stopInfo}{A character string describing
#'   why the selection process stopped.}
#'
#' @references Wang, S. (2024). A New Forward Discriminant Analysis Framework
#' Based On Pillai's Trace and ULDA. \emph{arXiv preprint arXiv:2409.03136}.
#' Available at \url{https://arxiv.org/abs/2409.03136}.
forwardSel <- function(m,
                       response,
                       testStat = "Pillai",
                       alpha = 0.1,
                       correction = TRUE){

  testStat <- match.arg(testStat, c("Pillai", "Wilks"))
  if (alpha < 0 || alpha > 1) stop("Parameter 'alpha' must be between 0 and 1")
  if (!is.logical(correction) || length(correction) != 1) stop("correction must be either TRUE or FALSE")

  # Initialization -----------------------------------------------

  N = nrow(m); J = nlevels(response); p = 0
  currentCandidates <- seq_len(ncol(m)); currentVarList = c()
  statSaved <- numeric(ncol(m) + 1) + ifelse(testStat == "Wilks", 1, 0) # (i+1)-th place saves i-th result
  statDiff <- numeric(ncol(m)); stopInfo <- "Normal"

  groupMeans <- tapply(c(m), list(rep(response, dim(m)[2]), col(m)), function(x) mean(x, na.rm = TRUE))
  mW <- m - groupMeans[response, , drop = FALSE]
  Sw <- St <- matrix(NA, nrow = ncol(m), ncol = ncol(m))
  diag(Sw) <- apply(mW^2, 2, sum) / N; diag(St) <- apply(m^2, 2, sum) / N

  forwardInfo <- data.frame(var = character(ncol(m)),
                            statOverall = NA,
                            statDiff = NA)

  if(N <= J) return(list(currentVarList = currentVarList, forwardInfo = forwardInfo, stopInfo = "No enough observations"))

  # Threshold for Wilks' Lambda
  maxVar <- min(ncol(m), N - J); maxVarSeq <- seq_len(maxVar) - 1
  if(correction){
    correctionFactor <- ncol(m) - maxVarSeq
  }else correctionFactor <- 1
  threshold <- (N - J - maxVarSeq) / (N - J - maxVarSeq + stats::qf(1 - alpha / correctionFactor, df1 = J - 1, N - J - maxVarSeq) * (J - 1))

  # Main Loop -----------------------------------------------

  while(length(currentCandidates) != 0 && p < maxVar){
    nCandidates <- length(currentCandidates)
    p = p + 1

    if(testStat == "Pillai"){ # Update threshold for Pillai
      Jnow <- J - statSaved[p]
      if(Jnow <= 1){
        stopInfo <- "Perfect separation"
        break
      }
      correctionFactor <- ifelse(correction, 1 / nCandidates, 1)
      threshold[p] <- stats::qbeta((1 - alpha)^correctionFactor, shape1 = (Jnow - 1) / 2, shape2 = (N - Jnow) / 2)
    }

    bestVarInfo <- getBestVar(currentVar = currentVarList,
                               newVar = currentCandidates,
                               Sw = Sw,
                               St = St,
                               testStat = testStat)
    bestVar <- bestVarInfo$varIdx

    if(bestVarInfo$stopflag){ # St = 0
      stopInfo <- "Perfect linear dependency"
      break
    }

    if(testStat == "Pillai"){
      statDiff[p] <- bestVarInfo$stat - statSaved[p]
    }else statDiff[p] <- bestVarInfo$stat / statSaved[p]
    stopFlag <- ifelse(testStat == "Pillai",
                       statDiff[p] < threshold[p],
                       statDiff[p] > threshold[p])

    if(stopFlag){
      stopInfo <- "No significant variables"
      break
    }

    # Save the info
    currentVarList <- c(currentVarList, bestVar)
    currentCandidates <- setdiff(currentCandidates, c(bestVar, bestVarInfo$collinearVar))
    forwardInfo$var[p] <- colnames(m)[bestVar]
    forwardInfo$statDiff[p] <- statDiff[p]
    forwardInfo$statOverall[p] <- statSaved[p+1] <- bestVarInfo$stat

    if(testStat == "Wilks" && bestVarInfo$stat == 0){ # Lambda = 0
      stopInfo <- "Wilks' Lambda = 0"
      break
    }

    # Update the Sw and St on the new added column
    Sw[currentCandidates, bestVar] <- Sw[bestVar, currentCandidates] <- as.vector(crossprod(mW[, currentCandidates, drop = FALSE], mW[, bestVar, drop = FALSE])) / N
    St[currentCandidates, bestVar] <- St[bestVar, currentCandidates] <- as.vector(crossprod(m[, currentCandidates, drop = FALSE], m[, bestVar, drop = FALSE])) / N
  }

  forwardInfo <- cbind.data.frame(forwardInfo[seq_along(currentVarList), , drop = FALSE],
                                  threshold = threshold[seq_along(currentVarList)])
  return(list(currentVarList = currentVarList, forwardInfo = forwardInfo, stopInfo = stopInfo))
}


#' Select Best Variable at Current Step Based on Multivariate Test Statistics
#'
#' This function selects the best variable based on the specified multivariate
#' test statistic (`Pillai` or `Wilks`). It evaluates the statistic for each
#' candidate variable in `newVar` when combined with `currentVar`, and returns
#' the index and test statistic of the best variable. It also identifies
#' collinear variables.
#'
#' @param currentVar A numeric vector indicating the indices of currently
#'   selected variables.
#' @param newVar A numeric vector indicating the indices of candidate variables
#'   to be tested.
#' @param Sw A matrix representing the within-class scatter matrix.
#' @param St A matrix representing the total scatter matrix.
#' @param testStat A character string specifying the test statistic to use. Can
#'   be either `"Pillai"` or `"Wilks"`. Default is `"Pillai"`.
#'
#' @return A list containing: \item{stopflag}{A logical value indicating whether
#'   the best variable is collinear (i.e., should the selection stop).}
#'   \item{varIdx}{The index of the selected variable from `newVar` based on the
#'   test statistic.} \item{stat}{The value of the test statistic for the
#'   selected variable.} \item{collinearVar}{A vector of indices from `newVar`
#'   representing collinear variables.}
getBestVar <- function(currentVar, newVar, Sw, St, testStat = "Pillai"){
  if(testStat == "Pillai"){
    statAll <- sapply(newVar, function(i) getPillai(Sw[c(i,currentVar),c(i,currentVar), drop = FALSE], St[c(i,currentVar),c(i,currentVar), drop = FALSE]))
    currentVarIdx <- which.max(statAll)
    collinearVar <- which(statAll == -1)
  }else if(testStat == "Wilks"){
    statAll <- sapply(newVar, function(i) getWilks(Sw[c(i,currentVar),c(i,currentVar), drop = FALSE], St[c(i,currentVar),c(i,currentVar), drop = FALSE]))
    currentVarIdx <- which.min(statAll)
    collinearVar <- which(statAll == 2)
  }

  return(list(stopflag = (currentVarIdx %in% collinearVar),
              varIdx = newVar[currentVarIdx],
              stat = statAll[currentVarIdx],
              collinearVar = newVar[collinearVar]))
}
