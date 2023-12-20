#' Suggesting values for \code{r2}
#'
#' This function suggests the values for \code{r2}.
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
#' @param r1 A numeric value indicating the maximum number of main effects.
#' @param r2 A numeric value indicating the maximum number of interaction effects.
#' @param interaction.ind A two-column numeric matrix containing all possible
#' two-way interaction effects. It must be generated outside of this function
#' using \code{t(utils::combn())}. See Example section for details.
#' @param pi1 A numeric value between 0 and 1, defined by users. Default is 0.32.
#' For guidance on selecting an appropriate value, please refer to \code{\link{ABC}}.
#' @param pi2 A numeric value between 0 and 1, defined by users. Default is 0.32.
#' For guidance on selecting an appropriate value, please refer to \code{\link{ABC}}.
#' @param pi3 A numeric value between 0 and 1, defined by users. Default is 0.32.
#' For guidance on selecting an appropriate value, please refer to \code{\link{ABC}}.
#' @param lambda A numeric value defined by users. Default is 10.
#' For guidance on selecting an appropriate value, please refer to \code{\link{ABC}}.
#' @param q A numeric value indicating the number of models in each generation (e.g.,
#' the population size). Default is 40.
#'
#' @return A \code{list} of output. The components are:
#' \item{InterRank}{Rank of all candidate interaction effects. A two-column numeric
#' matrix. The first column contains indices of ranked two-way interaction effects, and the
#' second column contains its corresponding ABC score.}
#' \item{mainind.sel}{Selected main effects.  A \code{r1}-dimensional vector.}
#' \item{mainpool}{Ranked main effects in \code{X}.}
#' \item{plot}{Plot of potential interaction effects and their corresponding ABC scores.}
#' @export
#'
#' @seealso \code{\link{initial}}.
#'
#' @examples # under Strong heredity
# set.seed(0)
# nmain.p <- 4
# interaction.ind <- t(combn(4,2))
# X <- matrix(rnorm(50*4,1,0.1), 50, 4)
# epl <- rnorm(50,0,0.01)
# y<- 1+X[,1]+X[,2]+X[,1]*X[,2]+epl
# d1 <- detect(X, y, nmain.p = 4, r1 = 3, r2 = 3,
#     interaction.ind = interaction.ind, q = 5)
#'
#' @examples # under No heredity
#' set.seed(0)
#' nmain.p <- 4
#' interaction.ind <- t(combn(4,2))
#' X <- matrix(rnorm(50*4,1,0.1), 50, 4)
#' epl <- rnorm(50,0,0.01)
#' y<- 1+X[,1]+X[,2]+X[,1]*X[,2]+epl
#' d2 <- detect(X, y, heredity = "No", nmain.p = 4, r1 = 3, r2 = 3,
#'     interaction.ind = interaction.ind, q = 5)
#'

detect <- function (X, y, heredity = "Strong",
                    nmain.p, sigma = NULL,  r1, r2,
                    interaction.ind = NULL,
                    pi1 = 0.32, pi2 = 0.32, pi3 = 0.32,
                    lambda = 10, q = 40){

  bbb <- int(X, y, heredity = heredity,
             nmain.p = nmain.p, sigma = sigma,  r1 = r1, r2 = r2,
             interaction.ind = interaction.ind,
             pi1 = pi1, pi2 = pi2, pi3 = pi3,
             lambda = lambda, q = q)

  interpool <- bbb$InterRank
  ccc <- as.data.frame(interpool)
  inter <- ccc[,1]
  scores <- ccc[,2]
  if (dim(interpool)[1] <= 50){
    gp <- ggplot2::ggplot(ccc,
                          ggplot2::aes(x = stats::reorder(as.character(inter),
                                                   +as.numeric(scores)), y = as.numeric(scores))) +
      geom_point() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
  }else{
    gp <- ggplot2::ggplot(as.data.frame(interpool)[1:50,],
                          ggplot2::aes(x = stats::reorder(as.character(inter),
                                                   +as.numeric(scores)), y = as.numeric(scores))) +
      geom_point() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
  }
  return(plot = gp)
}

int <- function(X, y, heredity = "Strong",
                 nmain.p, sigma = NULL,  r1, r2,
                 interaction.ind = NULL,
                 pi1 = 0.32, pi2 = 0.32, pi3 = 0.32,
                 lambda = 10, q = 40){
  if (is.null(interaction.ind)) stop("Interaction.ind is missing. Use t(utils::combn()) to generate interaction matrix.")
  colnames(X) <- make.names(rep("","X",ncol(X)+1),unique=TRUE)[-1]
  max_model_size <- r1 + r2

  DCINFO <- DCSIS(X,y,nsis=(dim(X)[1])/log(dim(X)[1]))
  mainind <- as.numeric(gsub(".*?([0-9]+).*", "\\1",  colnames(X)[order(DCINFO$rankedallVar)]))
  Shattempind <- mainind[1:r1]

  # no heredity pool
  if (heredity == "No"){
    interpooltemp <- interaction.ind
  }

  # weak heredity pool
  if (heredity =="Weak"){
    for (i in 1:q) {
      df <- rbind(interaction.ind[interaction.ind[,1] %in% Shattempind,][order(stats::na.exclude(match(interaction.ind[,1], Shattempind))),],
                  interaction.ind[interaction.ind[,2] %in% Shattempind,][order(stats::na.exclude(match(interaction.ind[,2], Shattempind))),])
      interpooltemp <- df[!duplicated(df),]
    }
  }

  # strong heredity pool
  if (heredity =="Strong"){
    interpooltemp <- t(utils::combn(sort(Shattempind),2))
  }

  intercandidates.ind <- match(do.call(paste, as.data.frame(interpooltemp)), do.call(paste, as.data.frame(interaction.ind)))+nmain.p

  interscoreind <- list()
  for (i in 1:length(intercandidates.ind)) {
    interscoreind[[i]] <- ABC(X, y, heredity = heredity, nmain.p = nmain.p, sigma = sigma,
                              extract = "Yes", varind = c(Shattempind,intercandidates.ind[i]),
                              interaction.ind = interaction.ind,
                              pi1 = pi1, pi2 = pi2, pi3= pi3, lambda = lambda)
  }
  interscoreind <- interscoreind
  MA <- as.matrix(cbind(inter = intercandidates.ind, scores = interscoreind))
  if (dim(MA)[1] == 1){
    MB <- MA
  }else{
    MB <- MA[order(unlist(MA[,2]), na.last = TRUE),]
  }
  interpool <- MB
  return(list(InterRank = interpool,
              mainind.sel = Shattempind,
              mainpool = mainind))
}

DCSIS <- function(X,Y,nsis=(dim(X)[1])/log(dim(X)[1])){
  if (dim(X)[1]!=length(Y)) {
    stop("X and Y should have same number of rows")
  }
  if (missing(X)|missing(Y)) {
    stop("The data is missing")
  }
  if (TRUE%in%(is.na(X)|is.na(Y)|is.na(nsis))) {
    stop("The input vector or matrix cannot have NA")
  }
  n=dim(X)[1]
  p=dim(X)[2]
  B=matrix(1,n,1)
  C=matrix(1,1,p)
  sxy1=matrix(0,n,p)
  sxy2=matrix(0,n,p)
  sxy3=matrix(0,n,1)
  sxx1=matrix(0,n,p)
  syy1=matrix(0,n,1)
  for (i in 1:n){
    XX1=abs(X-B%*%X[i,])
    YY1=sqrt(apply((Y-B%*%Y[i])^2,1,sum))
    sxy1[i,]=apply(XX1*(YY1%*%C),2,mean)
    sxy2[i,]=apply(XX1,2,mean)
    sxy3[i,]=mean(YY1)
    XX2=XX1^2
    sxx1[i,]=apply(XX2,2,mean)
    YY2=YY1^2
    syy1[i,]=mean(YY2)
  }
  SXY1=apply(sxy1,2,mean)
  SXY2=apply(sxy2,2,mean)*apply(sxy3,2,mean)
  SXY3=apply(sxy2*(sxy3%*%C),2,mean)
  SXX1=apply(sxx1,2,mean)
  SXX2=apply(sxy2,2,mean)^2
  SXX3=apply(sxy2^2,2,mean)
  SYY1=apply(syy1,2,mean)
  SYY2=apply(sxy3,2,mean)^2
  SYY3=apply(sxy3^2,2,mean)
  dcovXY=sqrt(SXY1+SXY2-2*SXY3)
  dvarXX=sqrt(SXX1+SXX2-2*SXX3)
  dvarYY=sqrt(SYY1+SYY2-2*SYY3)
  dcorrXY=dcovXY/sqrt(dvarXX*dvarYY)
  A=order(dcorrXY,decreasing=TRUE)
  return (list(rankedallVar = A,
               rankednsisVar = A[1:min(length(order(dcorrXY,decreasing=TRUE)),nsis)],
               scoreallVar = dcorrXY)
  )
}

indchunked <- function(n, chunk_size) {
  chunk <- matrix(0, chunk_size, 2)
  idx <- 1
  chunk_idx <- 1
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      chunk[chunk_idx,] <- c(i, j)
      chunk_idx <- chunk_idx + 1
      if(chunk_idx > chunk_size) {
        chunk <- chunk
        chunk_idx <- 1
      }
      idx <- idx + 1
    }
  }
  if(chunk_idx > 1) {
    chunk[1:(chunk_idx-1),]
    chunk <- chunk
  }
  return(chunk)
}

