#' Extracting specific columns from a data\cr
#'
#' This function extracts specific columns from \code{X} based on \code{varind}.
#' It provides an efficient procedure for conducting ABC evaluation,
#' especially when working with high-dimensional data.
#'
#' @details Please be aware that this function automatically renames column names
#'  into a designated format (e.g., X.1, X.2 for main effects, and X.1X.2 for
#'  interaction effect, etc), regardless of the original column names in \code{X}.
#'
#'  Under no heredity condition, this function can be applied in the context of
#'  interaction only linear regression models. See Example section for details.
#'
#' @param X Input data. An optional data frame, or numeric matrix of dimension
#' \code{n} by \code{nmain.p}. Note that the two-way interaction effects should not
#' be included in \code{X} because this function automatically generates the
#' corresponding two-way interaction effects if needed.
#' @param varind A numeric vector of class \code{c()} that specifies the indices
#' of variables to be extracted from \code{X}. Duplicated values are not allowed.
#' See Example section for details.
#' @param interaction.ind A two-column numeric matrix containing all possible
#' two-way interaction effects. It must be generated outside of this function using
#' \code{t(utils::combn())}. See Example section for details.
#'
#' @return A numeric matrix is returned.
#' @export
#'
#' @seealso \code{\link{ABC}}, \code{\link{initial}}.
#'
#' @examples # Extract main effect X1 and X2 from X1,...X4
#' set.seed(0)
#' X1 <- matrix(rnorm(20), ncol = 4)
#' y1 <- X1[, 2] + rnorm(5)
#' interaction.ind <- t(combn(4,2))
#'
#' @examples # Extract main effect X1 and interaction effect X1X2 from X1,..X4
#' Extract(X1, varind = c(1,5), interaction.ind)
#'
#' @examples # Extract interaction effect X1X2 from X1,...X4
#' Extract(X1, varind = 5, interaction.ind)
#'
#' @examples # Extract using duplicated values in varind.
#' try(Extract(X1, varind = c(1,1), interaction.ind)) # this will not run

Extract <- function(X, varind, interaction.ind = NULL){
  if (is.null(interaction.ind)) stop("Interaction.ind is missing.
                                       Use t(utils::combn()) to generate interaction matrix.")
  if (as.logical(any(duplicated(varind[which(varind!=0)])))){
    stop("There cannot be duplicated values in varind.")
  }
  ncoln <- length(varind)
  nrown <- nrow(X)
  nmain.p <- ncol(X)
  mainind <- varind[which(varind%in%1:nmain.p)]
  mainvars.matrix <- X[, mainind]
  interind <- varind[which(varind>nmain.p)]

  a <- interaction.ind[interind-nmain.p,]
  if (length(a) >1){
    intervars.matrix <- matrix(0, nrow = nrown,ncol = mydim(a)[1])
    for (i in 1:mydim(a)[1]) {
      if (mydim(a)[1]==1) intervars.matrix[,i] <- X[, a[1]]*X[,a[2]]
      else {
        intervars.matrix[,i] <- X[, a[i,1]]*X[,a[i,2]]
      }
    }
    colnames(intervars.matrix) <- paste0("X.", interind)
    data_extract <- cbind(mainvars.matrix, intervars.matrix)
  }
  else{
    data_extract <- mainvars.matrix
  }
  if (length(mainind) == 1){
    data_extract <- as.matrix(data_extract)
    colnames(data_extract)[which(colnames(data_extract)=="mainvars.matrix" )] <- paste0("X.", mainind)
  }
  if (length(mainind) == 1 && length(interind) ==0){
    colnames(data_extract) <- paste0("X.", mainind)
  }
  if (length(mainind) == 0){
    data_extract <- data_extract
  }
  if (length(mainind)>1){
    colnames(data_extract)[1:dim(mainvars.matrix)[2]] <- paste0("X.", mainind)
  }
  return(data_extract)
}

mychoose <- function(k1){
  if (k1 == 1){
    return(1)
  }else{
    return(choose(k1,2))
  }
}

mydim <- function(x){
  if (is.matrix(x)) dim(x)
  else return(1)
}

mysample <- function(x, size, replace = F, prob = NULL){
  if (length(x) == 1) return(x)
  if (length(x) > 1) return(sample(x, size, replace, prob))
}

sort_zeros <- function(vec){
  non_zeros <- vec[vec != 0]
  zeros <- vec[vec == 0]
  if (length(zeros) > 0){
    return(c(sort(non_zeros), zeros))
  }else{
    return(c(sort(non_zeros)))
  }
}

unique_rows <- function(matrix) {
  unique_matrix <- matrix[!duplicated(matrix[,ncol(matrix)]), ]
  return(unique_matrix)
}
