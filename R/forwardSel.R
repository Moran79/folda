#' Title
#'
#' @param datX
#' @param response
#' @param method
#' @param firstOnly
#' @param alpha
#' @param correction
#' @param fixThreshold
#'
#' @return
#' @export
#'
#' @examples
forwardSel <- function(datX, response, method = "Pillai", firstOnly = FALSE, alpha = 0.05, correction = TRUE, fixThreshold = FALSE){
  #> Given a matrix, output the important vars
  #> method = Wilks / Pillai

  # Pre-processing
  response <- droplevels(as.factor(response)) # some levels are branched out

  # constant factors
  levelObs <- sapply(datX, nlevels)
  if(any(levelObs == 1)) datX[,levelObs == 1] <- 0

  if(is.data.frame(datX)){
    modelFrame <- model.frame(formula = ~.-1, datX, na.action = "na.fail")
    Terms <- terms(modelFrame)
    m <- scale(model.matrix(Terms, modelFrame)) # constant cols would be changed to NaN in this step
  }else if (is.matrix(datX)){
    m <- scale(datX); colnames(m) <- paste0("X", seq_len(ncol(m)))
  }
  idxOriginal <- currentCandidates <- as.vector(which(apply(m, 2, function(x) !any(is.nan(x))))) # remove constant columns and intercept

  # Program starts
  m <- m[, currentCandidates, drop = FALSE] # all columns should be useful
  groupMeans <- tapply(c(m), list(rep(response, dim(m)[2]), col(m)), function(x) mean(x, na.rm = TRUE))
  mW <- m - groupMeans[response, , drop = FALSE]

  # Initialize
  n = nrow(m); g = nlevels(response); p = 0; currentVarList = c()
  previousPillai <- previousDiff <- previousDiffDiff <- numeric(ncol(m)+1) + ifelse(method == "Wilks", 1, 0);
  previousDiff[1] <- Inf; diffChecker <- 0
  kRes <- 1; currentCandidates <- seq_len(ncol(m))
  Sw <- St <- matrix(NA, nrow = ncol(m), ncol = ncol(m))
  diag(Sw) <- apply(mW^2,2,sum) / n; diag(St) <- apply(m^2,2,sum) / n
  stopFlag <- 0;pillaiThreshold <- numeric(ncol(m))

  if(method == "Wilks"){
    pillaiThreshold <- numeric(ncol(m))
    maxVar <- seq_len(min(ncol(m), n-g-1))
    if(correction){
      correctionPower <- ncol(m) + 1 - maxVar
    }else correctionPower <- 1

    if(fixThreshold){
      pillaiThreshold <- (n-g-maxVar) / (n-g-maxVar + 4 * (g-1)) # 4 as threshold
    }else pillaiThreshold[maxVar] <- (n-g-maxVar) / (n-g-maxVar + qf(1 - alpha / correctionPower, df1 = g - 1, n-g-maxVar) * (g-1))
  }

  stepInfo <- data.frame(var = character(2*ncol(m)),
                         pillaiToEnter = 0,
                         partialWilks = 0,
                         threshold = pillaiThreshold,
                         pillai = 0)

  #> If n <= g, which means there are too few observations,
  #> we output all columns, and leave that problem to outside function
  if(anyNA(pillaiThreshold)){
    currentVarList <- currentCandidates; currentCandidates <- c()
    stepInfo$var[seq_along(currentVarList)] <- colnames(m)[currentVarList]
    stopFlag <- 4
  }

  # Stepwise selection starts!
  while(length(currentCandidates) != 0){
    nCandidates <- length(currentCandidates)
    p = p + 1
    if(method == "Pillai"){
      gNow <- g - previousPillai[p]
      correctionPower <- ifelse(correction, 1 / nCandidates, 1)
      pillaiThreshold[p] <- qbeta((1 - alpha)^correctionPower, shape1 = (gNow - 1) / 2, shape2 = (n - gNow) / 2)
    }

    if(firstOnly & p > 1){ # Only first variable needed to check type I error
      stopFlag <- 6
      break
    }

    # if(length(currentVarList) >= 5){ ### TESTING ON IRIS ONLY ###
    #   stopFlag <- 7
    #   break
    # }

    selectVarInfo <- selectVar(currentVar = currentVarList,
                               newVar = currentCandidates,
                               Sw = Sw,
                               St = St,
                               method = method)
    bestVar <- selectVarInfo$varIdx
    if(selectVarInfo$stopflag){ # If St = 0, stop. [Might never happens, since there are other variables to choose]
      stopFlag <- 1
      break
    }

    # get the difference in Pillai's trace
    if(method == "Pillai"){
      previousDiff[p+1] <- selectVarInfo$statistics - previousPillai[p]
    }else{
      previousDiff[p+1] <- selectVarInfo$statistics / previousPillai[p]
      stepInfo$partialWilks[p] <- (n-g-p+1) / (g-1) * (1-previousDiff[p+1]) / previousDiff[p+1]
    }

    # Check the stopping rule
    stopFlag <- ifelse(method == "Pillai",
                       previousDiff[p+1] < pillaiThreshold[p],
                       previousDiff[p+1] > pillaiThreshold[p])

    if(is.na(stopFlag)){ # Lambda = 0
      stopFlag <- 5
      break
    }

    if(stopFlag){ # If no significant variable selected, stop
      stopFlag <- 2
      break
    }

    # if(diffChecker == 10){ # converge
    #   stopFlag <- 3
    #   break
    # }

    # Add the variable into the model
    previousPillai[p+1] <- selectVarInfo$statistics
    currentVarList <- c(currentVarList, bestVar)
    currentCandidates <- setdiff(currentCandidates, bestVar)
    stepInfo$var[kRes] <- colnames(m)[bestVar]
    stepInfo$pillaiToEnter[kRes] <- previousDiff[p+1]
    stepInfo$pillai[kRes] <- previousPillai[p+1]
    kRes <- kRes + 1

    # Update the Sw and St on the new added column
    Sw[currentCandidates, bestVar] <- Sw[bestVar, currentCandidates] <- as.vector(t(mW[, currentCandidates, drop = FALSE]) %*% mW[,bestVar, drop = FALSE]) / n
    St[currentCandidates, bestVar] <- St[bestVar, currentCandidates] <- as.vector(t(m[, currentCandidates, drop = FALSE]) %*% m[,bestVar, drop = FALSE]) / n

    # Update Pillai Threshold
    # getBetaThresholdSim <- function(m, response, alpha = 0.05, nSim = 100){
    #   fakeCols <- matrix(rnorm(nSim * nrow(m)), ncol = 100)
    #   apply(fakeCols, 2, function(x) 1)
    # }
  }

  # Remove the empty rows in the stepInfo if stepLDA does not select all variables
  stepInfo$threshold = pillaiThreshold # Update the threshold used
  stepInfo <- stepInfo[seq_along(currentVarList),]

  return(list(currentVarList = idxOriginal[currentVarList], stepInfo = stepInfo, stopFlag = stopFlag))
}
