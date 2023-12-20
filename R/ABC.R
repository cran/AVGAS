#' Evaluating ABC for each fitted model\cr
#'
#' This function evaluates ABC score for fitted model, one model at a time. For a model I,
#' the ABC is defined as
#' \deqn{ABC(I)=\sum\limits_{i=1}^n\bigg(Y_i-\hat{Y}_i^{I}\bigg)^2+2r_I\sigma^2+\lambda\sigma^2C_I.}
#' When comparing ABC of fitted models to the same data set, the smaller
#' the ABC, the better fit.
#'
#' @details
#' \itemize{
#' \item For inputs \code{pi1}, \code{pi2}, and \code{pi3}, the number needs to
#' satisfy the condition: \eqn{\pi_1+\pi_2+\pi_3=1-\pi_0} where \eqn{\pi_0}
#' is a numeric value between 0 and 1, the smaller the better.
#' \item For input \code{lambda}, the number needs to satisfy the condition:
#' \eqn{\lambda\geq 5.1/log(2)}.
#' }
#'
#' @param X Input data. An optional data frame, or numeric matrix of dimension
#' \code{n} by \code{nmain.p}. Note that the two-way interaction effects should not
#' be included in \code{X} because this function automatically generates the
#' corresponding two-way interaction effects if needed.
#' @param y Response variable. A \code{n}-dimensional vector, where \code{n} is the number
#' of observations in \code{X}.
#' @param heredity Whether to enforce Strong, Weak, or No heredity. Default is "Strong".
#' @param nmain.p A numeric value that represents the total number of main effects
#' in \code{X}.
#' @param sigma The standard deviation of the noise term. In practice, sigma is usually
#' unknown. In such case, this function automatically estimate sigma using root mean
#' square error (RMSE). Default is NULL. Otherwise, users need to enter a numeric value.
#' @param extract A either "Yes" or "No" logical vector that represents whether or not
#' to extract specific columns from \code{X}. Default is "No".
#' @param varind Only used when \code{extract = "Yes"}. A numeric vector of class
#' \code{c()} that specifies the indices of variables to be extracted from \code{X}.
#' If \code{varind} contains indices of two-way interaction effects, then this function
#' automatically generates corresponding two-way interaction effects from \code{X}.
#' @param interaction.ind Only used when \code{extract = "Yes"}. A two-column numeric
#' matrix containing all possible two-way interaction effects. It must be generated
#' outside of this function using \code{t(utils::combn())} or \code{indchunked()}. See Example section for
#' details.
#' @param pi1 A numeric value between 0 and 1, defined by users. Default is 0.32.
#' For guidance on selecting an appropriate value, please refer to the Details section.
#' @param pi2 A numeric value between 0 and 1, defined by users. Default is 0.32.
#' For guidance on selecting an appropriate value, please refer to the Details section.
#' @param pi3 A numeric value between 0 and 1, defined by users. Default is 0.32.
#' For guidance on selecting an appropriate value, please refer to the Details section.
#' @param lambda A numeric value defined by users. Default is 10.
#' For guidance on selecting an appropriate value, please refer to the Details section.
#'
#' @return A numeric value is returned. It represents the ABC score of the fitted model.
#'
#' @export
#'
#' @seealso \code{\link{Extract}}, \code{\link{initial}}.
#' @references
#' Ye, C. and Yang, Y., 2019. \emph{High-dimensional adaptive minimax sparse estimation with interactions.}
#'
#' @examples # sigma is unknown
#' set.seed(0)
#' nmain.p <- 4
#' interaction.ind <- t(combn(4,2))
#' X <- matrix(rnorm(50*4,1,0.1), 50, 4)
#' epl <- rnorm(50,0,0.01)
#' y<- 1+X[,1]+X[,2]+X[,3]+X[,4]+epl
#' ABC(X, y, nmain.p = 4, interaction.ind = interaction.ind)
#' ABC(X, y, nmain.p = 4, extract = "Yes",
#'     varind = c(1,2,5), interaction.ind = interaction.ind)
#' #'
#' @examples # users want to enter a suggested value for sigma
# ABC(X, y, nmain.p = 4, sigma = 0.01)
# ABC(X, y, nmain.p = 4, sigma = 0.01, extract = "Yes",
#     varind = c(1,2,5), interaction.ind = interaction.ind)
#'
#' @examples # model with only one predictor
#' try(ABC(X, y, nmain.p = 4, extract = "Yes",
#'   varind = 1, interaction.ind = interaction.ind)) # warning message


ABC <- function(X, y, heredity = "Strong", nmain.p, sigma = NULL,
                extract = "No", varind = NULL, interaction.ind = NULL,
                pi1 = 0.32, pi2 = 0.32, pi3 = 0.32, lambda = 10){
  colnames(X) <- make.names(rep("","X",ncol(X)+1),unique=TRUE)[-1]
  n <- dim(X)[1]

  if (is.null(sigma)){
    full <- Extract(X, varind = c(1:(dim(X)[2]+dim(interaction.ind)[1])), interaction.ind)
    sigma <- estimateSigma(full, y)$sigmahat
  }else{
    sigma <- sigma
  }

  if (extract == "Yes"){
    if (is.null(varind)) stop("You must specify the variables to be extracted")
    if (is.null(interaction.ind)) stop("Interaction.ind is missing.
                                       Use t(utils::combn()) to generate interaction matrix.")
    data_extract <- Extract(X, varind, interaction.ind)
    data <- data_extract
    r.I <- Matrix::rankMatrix(data)[1]
    allpar <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(data)))
  }else{
    r.I <- Matrix::rankMatrix(X)[1]
    allpar <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(X)))
    data <- X
  }
  k2 <- sum(allpar > nmain.p)
  k1 <- length(allpar) - k2
  if (r.I < ncol(data)){
    yhat <- pracma::orth(data)%*%solve(crossprod(pracma::orth(data)),
                                       t(pracma::orth(data))%*%y, tol = 1e-50)
  }else{
    yhat <- (data)%*%solve(crossprod(data), t(data)%*%y, tol = 1e-50)
  }
  SSE <- sum((y - yhat)^2)

  if (heredity == "Strong"){
    C.I.strong <- -log(pi1)+log(min(nmain.p,n))+log(min(mychoose(k1),n))
    + log(choose(nmain.p,k1))+ log(choose(choose(k1,2),k2))
    ABC <- SSE+2*r.I*sigma^2+lambda*sigma^2*C.I.strong
    if(k1 == 1) warning("This model contains only one predictor")
  }
  else if (heredity == "Weak"){
    K <- k1*nmain.p-choose(k1,2)-k1
    C.I.weak <- -log(pi2)+log(min(nmain.p,n))+log(min(K,n))
    + log(choose(nmain.p,k1))+ log(choose(K,k2))
    ABC <- SSE+2*r.I*sigma^2+lambda*sigma^2*C.I.weak
    if(k1 == 1) warning("This model contains only one predictor")
  }
  else if (heredity == "No"){
    C.I.no <- -log(pi3)+log(min(nmain.p,n))+log(min(choose(nmain.p,2),n))
    + log(choose(nmain.p,k1))+ log(choose(choose(nmain.p,2),k2))
    ABC <- SSE+2*r.I*sigma^2+lambda*sigma^2*C.I.no
    if(k1 == 1) warning("This model contains only one predictor")
  }
  return(ABC)
}

