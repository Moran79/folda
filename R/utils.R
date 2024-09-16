#' Impute Missing Values and Add Missing Flags to a Data Frame
#'
#' This function imputes missing values in a data frame based on specified
#' methods for numerical and categorical variables. Additionally, it can add
#' flag columns to indicate missing values. For numerical variables, missing
#' values can be imputed using the mean or median. For categorical variables,
#' missing values can be imputed using the mode or a new level. This function
#' also removes constant columns (all NAs or all observed but the same value).
#'
#' @param data A data frame containing the data to be processed. Missing values
#'   (`NA`) will be imputed based on the methods provided in `missingMethod`.
#' @param missingMethod A character vector of length 2 specifying the methods
#'   for imputing missing values. The first element specifies the method for
#'   numerical variables (`"mean"`, `"median"`, `"meanFlag"`, or
#'   `"medianFlag"`), and the second element specifies the method for
#'   categorical variables (`"mode"`, `"modeFlag"`, or `"newLevel"`). If
#'   `"Flag"` is included, a flag column will be added for the corresponding
#'   variable type.
#'
#' @return A list with two elements: \item{data}{The original data frame with
#'   missing values imputed, and flag columns added if applicable.} \item{ref}{A
#'   reference row containing the imputed values and flag levels, which can be
#'   used for future predictions or reference.}
#'
#' @export
#'
#' @examples
#' dat <- data.frame(
#'   X1 = rep(NA, 5),
#'   X2 = factor(rep(NA, 5), levels = LETTERS[1:3]),
#'   X3 = 1:5,
#'   X4 = LETTERS[1:5],
#'   X5 = c(NA, 2, 3, 10, NA),
#'   X6 = factor(c("A", NA, NA, "B", "B"), levels = LETTERS[1:3])
#' )
#' missingFix(dat)
missingFix <- function(data, missingMethod = c("medianFlag", "newLevel")){

  # Preprocessing: Argument #
  if (length(missingMethod) != 2) stop("missingMethod must be a vector of length 2")
  numMethod <- match.arg(missingMethod[1], c("mean", "median", "meanFlag", "medianFlag"))
  catMethod <- match.arg(missingMethod[2], c("mode", "modeFlag", "newLevel"))

  # remove constant columns
  data <- data[, nonConstInd(data = data), drop = FALSE]
  if(any(dim(data) == 0)) return(list(data = data, ref = data)) # all columns are constant

  # Add Missing Columns #
  numOrNot <- getNumFlag(data); NAorNot <- sapply(data, anyNA)

  createFlagColumns <- function(addFlag, targetBool){
    if(addFlag & sum(targetBool) > 0){
      NAcols <- which(targetBool & NAorNot)
      if(length(NAcols) > 0){
        flagCols <- is.na(data[, NAcols, drop = FALSE]) + 0
        colnames(flagCols) <- paste(colnames(flagCols), "FLAG", sep = "_")
        data <- cbind(data, flagCols)
      }
    }
    return(data)
  }
  data <- createFlagColumns(addFlag = grepl("Flag", numMethod), targetBool = numOrNot)
  data <- createFlagColumns(addFlag = grepl("Flag", catMethod), targetBool = !numOrNot)

  # Imputation Starts #
  dataNRef <- rbind(data, NA) # add ref to the last row of data, initialize using NAs
  dataNRef[nrow(dataNRef), seq_len(ncol(dataNRef) - length(numOrNot)) + length(numOrNot)] <- 1 # flags

  for(i in seq_along(numOrNot)){
    #> This for loop is inevitable since we have to add factor levels to ref,
    #> and add the new level to the original data
    if(numOrNot[i]){ # numerical variables
      targetValue <- do.call(ifelse(grepl("mean", numMethod), "mean", "median"),
                             list(dataNRef[,i], na.rm = TRUE))
    } else{ # categorical vars

      if(is.factor(dataNRef[,i])){
        currentLevels <- levels(droplevels(dataNRef[,i])) # drop extra levels
      }else{
        currentLevels <- sort(unique(stats::na.omit(dataNRef[,i])))
      }

      if(grepl("newLevel", catMethod) & NAorNot[i]){
        currentLevels <- c(currentLevels, "new0_0Level")
        targetValue <- "new0_0Level"
      } else targetValue <- getMode(dataNRef[,i])

      dataNRef[, i] <- factor(dataNRef[,i], levels = currentLevels)
    }
    dataNRef[is.na(dataNRef[,i]), i] <- targetValue
  }

  # remove constant columns after imputation
  dataNRef <- dataNRef[, nonConstInd(data = dataNRef), drop = FALSE]
  missingRef <- dataNRef[nrow(dataNRef),,drop = FALSE]
  rownames(missingRef) <- NULL

  return(list(data = dataNRef[-nrow(dataNRef),,drop = FALSE],
              ref = missingRef))
}


#' Align Data with a Missing Reference
#'
#' This function aligns a given dataset (`data`) with a reference dataset
#' (`missingReference`). It ensures that the structure, column names, and factor
#' levels in `data` match the structure of `missingReference`. If necessary,
#' missing columns are initialized with `NA`, and factor levels are adjusted to
#' match the reference. Additionally, it handles the imputation of missing
#' values based on the reference and manages flag variables for categorical or
#' numerical columns.
#'
#' @param data A data frame to be aligned and adjusted according to the
#'   `missingReference`.
#' @param missingReference A reference data frame that provides the structure
#'   (column names, factor levels, and missing value reference) for aligning
#'   `data`.
#'
#' @return A data frame where the structure, column names, and factor levels of
#'   `data` are aligned with `missingReference`. Missing values in `data` are
#'   imputed based on the first row of the `missingReference`, and flag
#'   variables are updated accordingly.
#'
#' @export
#'
#' @examples
#' data <- data.frame(
#'   X1_FLAG = c(0, 0, 0),
#'   X1 = factor(c(NA, "C", "B"), levels = LETTERS[2:3]),
#'   X2_FLAG = c(NA, 0, 1),
#'   X2 = c(2, NA, 3)
#' )
#'
#' missingReference <- data.frame(
#'   X1_FLAG = 1,
#'   X1 = factor("A", levels = LETTERS[1:2]),
#'   X2 = 1,
#'   X2_FLAG = 1
#' )
#'
#' getDataInShape(data, missingReference)
getDataInShape <- function(data, missingReference){
  cname <- colnames(missingReference)
  refIdx <- match(cname, colnames(data))
  if(anyNA(refIdx)){ # initialize with NA if missing columns
    data[, cname[which(is.na(refIdx))]] <- NA
    refIdx <- match(cname, colnames(data))
  }

  for(currentIdx in seq_len(ncol(missingReference))){ # Main program starts
    #> The tricky part is the iterator is based on the missingReference, NOT the data.
    dataIdx <- refIdx[currentIdx]; currentVarName <- cname[currentIdx]

    ### New-level Fix for Categorical Variable ###
    numOrNot <- getNumFlag(missingReference[, currentIdx])
    if(!numOrNot & !identical(levels(data[, dataIdx]), levels(missingReference[, currentIdx])))
      data[, dataIdx] <- factor(data[, dataIdx], levels = levels(missingReference[, currentIdx]))

    missingOrNot <- is.na(data[, dataIdx])
    if(!any(missingOrNot)) next

    ## Scenario 1: It has a related flag variable in the data ##
    #> Only modify those flags where the original variable is missing
    #> Keep other parts still, since there could already be imputed values
    #> in the original variable that have been taken care of
    currentFlagIdx <- which(cname == paste(currentVarName,"FLAG",sep = "_"))
    if(length(currentFlagIdx) == 1) data[which(missingOrNot), refIdx[currentFlagIdx]] <- 1

    ## Scenario 2: It is a flag and it has an original variable in (or not in) the data ##
    #> Only impute those NAs in the flag, but keep the values that are already in the flag
    if(grepl("_FLAG$", currentVarName)){
      orginalVarIdx <- which(cname == sub("_FLAG$", "", currentVarName))
      if(length(orginalVarIdx) == 1){
        data[which(missingOrNot), dataIdx] <- is.na(data[which(missingOrNot), refIdx[orginalVarIdx]]) + 0
      } else data[which(missingOrNot), dataIdx] <- 1 # The original data is NOT found
      next
    }

    ### For numerical & categorical variables ###
    data[which(missingOrNot), dataIdx] <- missingReference[1, currentIdx]
  }

  return(data[, refIdx, drop = FALSE])
}


#' Identify Non-Constant Columns in a Data Frame
#'
#' @noRd
#'
#' @param data A data frame in which columns will be checked for constant
#'   values. Columns can be of any type (numeric, integer, logical, or factor).
#' @param tol A numeric tolerance value (default is `1e-8`) that applies to
#'   numerical columns.
#' @param na.rm A logical value. If `FALSE` (default), missing values are
#'   retained during the check.
#'
#' @return An integer vector containing the indices of the non-constant columns
#'   in the data frame. If all columns are constant, an empty vector is
#'   returned.
nonConstInd <- function(data, tol = 1e-8, na.rm = FALSE){

  if (!is.data.frame(data)) stop("data must be a data.frame")

  nonConstIndHelper <- function(x, tol){
    if(getNumFlag(x)) x <- round(x, digits = -log(tol, base = 10))
    if(na.rm) x <- stats::na.omit(x)
    return(length(unique(x)) > 1)
  }

  idxNotConst <- which(sapply(data, function(x) nonConstIndHelper(x = x, tol = tol)))

  return(idxNotConst)
}


#' Helper Function to Generate Training Set Indices Through Downsampling
#'
#' This function selects the indices for the training set based on the class
#' vector `response`. It allows for optional downsampling to balance the class
#' distribution by limiting the number of samples per class.
#'
#' @noRd
#'
#' @param response A factor vector representing the class labels.
#' @param downSampling A logical value indicating whether downsampling should be
#'   applied. If `TRUE`, downsampling is performed to limit the number of
#'   samples per class based on `kSample`. Note that this may not result in
#'   equal class frequencies, as `kSample` defines an upper limit for each
#'   class, not a lower limit.
#' @param kSample An integer specifying the maximum number of samples to be
#'   selected per class. If `NULL`, the number of samples is limited to the size
#'   of the smallest class.
#'
#' @return An integer vector of indices representing the selected samples from
#'   the original `response` vector.
getDownSampleInd <- function(response,
                             downSampling = FALSE,
                             kSample = NULL){
  if (!is.logical(downSampling) || length(downSampling) != 1) stop("downSampling must be either TRUE or FALSE")

  idxSelected <- which(!is.na(response)); response <- droplevels(response[idxSelected])

  if(downSampling){
    idxForEachClass <- split(seq_along(response), response)
    resTable <- sapply(idxForEachClass, length)
    if(is.null(kSample)) kSample <- min(resTable)
    nForEachClass <- pmin(resTable, kSample)
    finalIdx <- do.call(c, lapply(seq_along(resTable), function(i) idxForEachClass[[i]][sample(resTable[i], nForEachClass[i])]))
    idxSelected <- idxSelected[finalIdx]
  }

  return(idxSelected)
}


#' Identify Numeric, Integer, or Logical Columns in a Data Frame
#'
#' This function checks whether the columns in a data frame (or a vector) are of
#' type numeric, integer, or logical. It can return a logical vector indicating
#' whether each column matches these types, or, if `index = TRUE`, it returns
#' the indices of the matching columns.
#'
#' @noRd
#'
#' @param data A data frame or a vector. The function will check the data types
#'   of the columns (if `data` is a data frame) or the type of the vector.
#' @param index A logical value. If `FALSE` (default), the function returns a
#'   logical vector indicating which columns are numeric, integer, or logical.
#'   If `TRUE`, it returns the indices of these columns.
#' @return If `index = FALSE` (default), the function returns a logical vector
#'   with one element for each column (or the vector itself), where `TRUE`
#'   indicates that the column is of type numeric, integer, or logical, and
#'   `FALSE` indicates it is not. If `index = TRUE`, the function returns an
#'   integer vector containing the indices of the columns that are numeric,
#'   integer, or logical.
getNumFlag <- function(data, index = FALSE){
  #> This function should be replaced (or deleted) with caution,
  #> since the function in sapply framework output wrong classes.
  if(is.null(dim(data))) return(inherits(data, c("numeric", "integer", "logical")))

  if (!is.data.frame(data)) stop("data must be a data.frame")

  numOrNot <- sapply(data, function(col) inherits(col, c("numeric", "integer", "logical")))

  if(index){numOrNot <- which(numOrNot)}

  return(numOrNot)
}


#' Calculate the Mode of a Factor Variable with Optional Priors
#'
#' @noRd
#'
#' @param v A factor or vector that can be coerced into a factor. The mode will
#'   be calculated from the levels of this factor.
#' @param prior A numeric vector of prior weights for each level of the factor.
#'
#' @return The mode of the factor `v` as a character string. If all values are
#'   `NA`, the function returns `NA`.
getMode <- function(v, prior){
  #> NA will be ignored
  v <- as.factor(v)
  if(missing(prior)){
    prior = rep(1, nlevels(v)) # equal prior
  }else{
    if (is.null(names(prior))){
      stopifnot(length(prior) == nlevels(v))
      names(prior) <- levels(v)
    } else prior <- prior[match(levels(v), names(prior))]
  }

  summary_table <- table(v) * prior
  if(length(summary_table) == 0 | sum(summary_table) == 0) return(NA) # all NA
  return(names(which.max(summary_table)))
}


#' Check and Normalize Prior Probabilities and Misclassification Costs
#'
#' This function verifies and normalizes the provided prior probabilities and
#' misclassification cost matrix for a given response variable. It ensures that
#' the lengths of the prior and the dimensions of the misclassification cost
#' matrix match the number of levels in the response variable. If `prior` or
#' `misClassCost` are not provided, default values are used: the prior is set to
#' the observed frequencies of the response, and the misclassification cost
#' matrix is set to 1 for all misclassifications and 0 for correct
#' classifications.
#'
#' @noRd
#'
#' @param prior A numeric vector representing the prior probabilities for each
#'   class in the response variable. If `NULL`, the observed frequencies of the
#'   response are used as the default prior.
#' @param misClassCost A square matrix representing the misclassification costs
#'   for each pair of classes in the response variable. If `NULL`, a default
#'   misclassification matrix is created where all misclassifications have a
#'   cost of 1 and correct classifications have a cost of 0.
#' @param response A factor representing the response variable with multiple
#'   classes.
#'
#' @return A list containing: \item{prior}{A normalized vector of prior
#'   probabilities for each class.} \item{misClassCost}{A square matrix
#'   representing the misclassification costs, with rows and columns labeled by
#'   the levels of the response variable.}
checkPriorAndMisClassCost <- function(prior, misClassCost, response){

  matchWrapper <- function(nameObj, nameTarget){
    if (length(nameObj) != length(nameTarget))
      stop("Length mismatch: expected ", length(nameTarget), " names but got ", length(nameObj))
    targetIndex <- match(nameTarget, nameObj)
    if (anyNA(targetIndex))
      stop("Names in '", paste(nameTarget[is.na(targetIndex)], collapse = ", "), "' do not match with response levels.")
    return(targetIndex)
  }

  freqObs <- table(response, dnn = NULL) / length(response) # Default: Estimated Prior

  if (is.null(prior)) { # prior fix
    prior <- as.vector(freqObs)
  } else {
    if (!is.numeric(prior) || any(prior < 0)) stop("prior must be non-negative numbers")
    if (length(prior) != nlevels(response))
      stop("The length of 'prior' must be equal to the number of response levels (", nlevels(response), ").")
    if (!is.null(names(prior))) prior <- prior[matchWrapper(names(prior), levels(response))]
    prior <- as.vector(prior) / sum(prior)
  }

  if (is.null(misClassCost)){ # fix misClassCost fix
    misClassCost <- matrix(1, nrow = nlevels(response), ncol = nlevels(response)) - diag(1, nrow = nlevels(response))
  } else{
    if (!is.matrix(misClassCost) || !getNumFlag(as.vector(misClassCost)))
      stop("'misClassCost' must be a numeric matrix.")

    if(length(unique(c(dim(misClassCost), nlevels(response)))) != 1)
      stop("misclassification costs matrix has wrong dimension")

    if(!all.equal(colnames(misClassCost), rownames(misClassCost)))
      stop("misClassCost: colnames should be the same as rownames")

    if (!is.null(colnames(misClassCost))){
      misClassCost <- misClassCost[matchWrapper(colnames(misClassCost), levels(response)),
                                   matchWrapper(colnames(misClassCost), levels(response))]
    }
  }

  names(prior) <- colnames(misClassCost) <- rownames(misClassCost) <- levels(response)
  return(list(prior = prior, misClassCost = misClassCost))
}


#' Compute Linear Discriminant Scores
#'
#' @noRd
#'
#' @param modelLDA A fitted LDA model object containing the scaling matrix and
#'   the reference structure for missing data.
#' @param data A data frame containing the predictor variables for which to
#'   compute the linear discriminant scores.
#' @param nScores An integer specifying the number of discriminant scores to
#'   compute. If `-1` (default), all scores are computed.
#' @return A matrix of linear discriminant scores, where rows correspond to
#'   observations and columns correspond to the computed discriminant scores.
#'   If `nScores > 0`, only the specified number of scores is returned; otherwise,
#'   all scores are computed and returned.
getLDscores <- function(modelLDA, data, nScores = -1){
  data <- getDataInShape(data = data, missingReference = modelLDA$misReference)
  modelX <- getDesignMatrix(modelLDA = modelLDA, data = data, scale = TRUE)
  if(nScores > 0){
    nScores <- min(nScores, ncol(modelLDA$scaling))
    modelLDA$scaling <- modelLDA$scaling[, seq_len(nScores), drop = FALSE]
  }
  LDscores <- modelX %*% modelLDA$scaling
  return(LDscores)
}


#' Generate the Design Matrix for LDA Model
#'
#' @noRd
#'
#' @param modelLDA A fitted LDA model object containing the terms, variable
#'   indices, variable centers, and scaling factors.
#' @param data A data frame containing the predictor variables that are used to
#'   create the design matrix.
#' @param scale A logical value indicating whether to scale the design matrix
#'   based on the mean and standard deviation of the variables (default is
#'   `FALSE`).
#' @return A design matrix where each row corresponds to an observation and
#'   each column to a predictor variable. If `scale = TRUE`, the variables are
#'   centered and scaled based on the means and standard deviations provided in
#'   the LDA model object.
getDesignMatrix <- function(modelLDA, data, scale = FALSE){
  Terms <- stats::delete.response(modelLDA$terms)
  modelX <- stats::model.matrix(Terms, data = data, xlev = modelLDA$xlevels)
  if(scale){ # Reserved for internal usage from getLDscores
    modelX <- sweep(modelX[, modelLDA$varIdx, drop = FALSE], 2, modelLDA$varCenter, "-")
    modelX <- sweep(modelX, 2, modelLDA$varSD, "/")
  }
  return(modelX)
}


saferSVD <- function(x, uFlag = TRUE) {
  # For users that does not have the latest RcppEigen module

  fitSVD <- tryCatch({
    if (uFlag) {
      svdEigen(x)
    } else svdEigen(x, uFlag = FALSE)
  }, error = function(e) {NULL})

  if (!is.null(fitSVD) && !anyNA(fitSVD$v)) return(fitSVD)

  # if svdEigen is NULL / NA in the fitSVD -> fallback to base::svd
  fitSVD <- tryCatch({
    if (uFlag) {
      svd(x)
    } else svd(x, nu = 0L)
  }, error = function(e) {stop("Both svd modules failed, check your input data.")})

  return(fitSVD)
}





