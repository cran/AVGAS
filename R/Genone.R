#' Gathering useful information for first generation
#'
#' This function automatically ranks all candidate interaction effects under
#' Strong, Weak, or No heredity condition, compare and obtain first generation
#' candidate models. The selected models will be re-ordered so that main effects
#' come first, followed by interaction effects. Only two-way interaction effects
#' will be considered.
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
#' @param r1 A numeric value indicating the maximum number of main effects.
#' @param r2 A numeric value indicating the maximum number of interaction effects.
#' @param sigma The standard deviation of the noise term. In practice, sigma is usually
#' unknown. In such case, this function automatically estimate sigma using root mean
#' square error (RMSE). Default is NULL. Otherwise, users need to enter a numeric value.
#' @param interaction.ind A two-column numeric matrix containing all possible
#' two-way interaction effects. It must be generated outside of this function
#' using \code{t(utils::combn())}. See Example section for details.
#' @param lambda A numeric value defined by users. Default is 10.
#' For guidance on selecting an appropriate value, please refer to the Details section.
#' @param q A numeric value indicating the number of models in each generation (e.g.,
#' the population size). Default is 40.
#' @param allout Whether to print all outputs from this function. A "Yes" or "No"
#' logical vector. Default is "No". See Value section for details.
#' @param interonly Whether or not to consider fitted models with only two-way
#' interaction effects. A “Yes" or "No" logical vector. Default is "No".
#' @param pi1 A numeric value between 0 and 1, defined by users. Default is 0.32.
#' For guidance on selecting an appropriate value, please refer to \code{\link{ABC}}.
#' @param pi2 A numeric value between 0 and 1, defined by users. Default is 0.32.
#' For guidance on selecting an appropriate value, please refer to \code{\link{ABC}}.
#' @param pi3 A numeric value between 0 and 1, defined by users. Default is 0.32.
#' For guidance on selecting an appropriate value, please refer to \code{\link{ABC}}.
#' @param aprob A numeric value between 0 and 1, defined by users.
#' The addition probability during mutation. Default is 0.9.
#' @param dprob A numeric value between 0 and 1, defined by users.
#' The deletion probability during mutation. Default is 0.9.
#' @param aprobm A numeric value between 0 and 1, defined by users.
#' The main effect addition probability during addition. Default is 0.1.
#' @param aprobi A numeric value between 0 and 1, defined by users.
#' The interaction effect addition probability during addition. Default is 0.9.
#' @param dprobm A numeric value between 0 and 1, defined by users.
#' The main effect deletion probability during deletion. Default is 0.9.
#' @param dprobi A numeric value between 0 and 1, defined by users.
#' The interaction effect deletion probability during deletion. Default is 0.1.
#'
#' @return A list of output. The components are:
#' \item{newparents}{ New parents models used for t+1-th generation. A numeric matrix
#' of dimension \code{q} by \code{r1+r2} where each row represents a fitted model.
#' Duplicated models are allowed.}
#' \item{parents_models}{ A numeric matrix containing all fitted models from
#' \code{\link{initial}}, \code{\link{cross}}, and \code{\link{mut}} where each
#' row corresponding to a fitted model and each column representing the predictor
#' index in that model. Duplicated models are allowed.}
#' \item{parents_models_cleaned}{ A numeric matrix containing fitted models from
#' \code{\link{initial}}, \code{\link{cross}}, and \code{\link{mut}} with ABC scores.
#' Each row corresponding to a fitted model; the first 1 to \code{r1 + r2} columns
#' representing the predictor indices in that model, and the last column is a numeric value
#' representing the ABC score of that fitted model. Duplicated models are not allowed.}
#' \item{InterRank}{ Rank of all candidate interaction effects. A two-column numeric
#' matrix. The first column contains indices of ranked two-way interaction effects, and the
#' second column contains its corresponding ABC score.}
#'
#' @export
#' @seealso \code{\link{initial}}, \code{\link{cross}}, \code{\link{mut}}, \code{\link{ABC}}, and \code{\link{Extract}}.
#'
#' @examples # allout = "No"
#' set.seed(0)
#' nmain.p <- 4
#' interaction.ind <- t(combn(4,2))
#' X <- matrix(rnorm(50*4,1,0.1), 50, 4)
#' epl <- rnorm(50,0,0.01)
#' y <- 1+X[,1]+X[,2]+X[,1]*X[,2]+epl
#' g1 <- Genone(X, y, nmain.p = 4, r1= 3, r2=3,
#'     interaction.ind = interaction.ind, q = 5)
#'
#' @examples # allout = "Yes"
#' g2 <- Genone(X, y, nmain.p = 4, r1= 3, r2=3,
#'     interaction.ind = interaction.ind, q = 5, allout = "Yes")
Genone <- function(X, y, heredity = "Strong", nmain.p, r1, r2,
                   sigma = NULL, interaction.ind = NULL,
                   lambda = 10, q = 40, allout = "No",
                   interonly = "No",  pi1 = 0.32, pi2 = 0.32, pi3 = 0.32,
                   aprob = 0.9, dprob = 0.9, aprobm = 0.1, aprobi=0.9, dprobm = 0.9, dprobi = 0.1){
  if (is.null(interaction.ind)) stop("Interaction.ind is missing.
                                       Use t(utils::combn()) to generate interaction matrix.")

  initial_parents <- initial(X, y, heredity = heredity,
                             nmain.p = nmain.p, sigma = sigma, r1 = r1, r2 = r2,
                             interaction.ind = interaction.ind,
                             pi1 = pi1, pi2 = pi2, pi3= pi3, lambda = lambda, q = q)
  InterRank <- initial_parents$InterRank
  MatrixA <- initial_parents$initialize
  MatrixB <- cross(initial_parents, heredity = heredity,
                   nmain.p = nmain.p, r1 = r1, r2 = r2, interaction.ind = interaction.ind)
  MatrixC <-  mut(initial_parents, heredity = heredity, nmain.p = nmain.p, r1 = r1, r2 = r2,
                  interaction.ind = interaction.ind,
                  interonly = interonly, aprob = aprob, dprob = dprob,
                  aprobm = aprobm, aprobi=aprobi, dprobm = dprobm, dprobi = dprobi)
  UnionABC <- rbind(MatrixA, MatrixB, MatrixC)

  a <- t(apply(UnionABC, 1, function(x) {x <- sort(x, decreasing = TRUE);x}))
  b <- t(apply(a, 1, function(x) {x[x != 0] <- sort(x[x != 0]);x}))
  UnionABC_cleaned <- dplyr::distinct(as.data.frame(b))
  ABCscore.1 <- list()
  for (i in 1:dim(UnionABC_cleaned)[1]) {
    ABCscore.1[[i]] <- ABC(X, y, heredity = heredity, nmain.p = nmain.p, sigma = sigma,
                           extract = "Yes", varind = c(as.numeric(UnionABC_cleaned[i,][which(!UnionABC_cleaned[i,]==0)])),
                           interaction.ind = interaction.ind,
                           pi1 = pi1, pi2 = pi2, pi3= pi3, lambda = lambda)
  }
  ABCscore.1 <- as.matrix(ABCscore.1)

  if (dim(UnionABC_cleaned)[1] < q) {
    MatrixD <- stats::na.omit(UnionABC_cleaned[order(as.numeric(ABCscore.1))[1:q],])
  }else{
    MatrixD <- UnionABC_cleaned[order(as.numeric(ABCscore.1))[1:q],]
  }
  MatrixD <- as.matrix(MatrixD)
  rownames(MatrixD) <- colnames(MatrixD) <- NULL
  temp <- as.matrix(cbind(UnionABC_cleaned, m.scores = unlist(ABCscore.1)))
  modelw.score <- temp[order(unlist(temp[,((r1+r2)+1)])),]
  colnames(modelw.score) <- NULL
  if (allout == "Yes"){
    return(list(newparents = MatrixD, # Generation 1 parents matrix used for next generation
                parents_models = UnionABC, # Gen 1 initialize, crossover, mutation not cleaned
                parents_models_cleaned = modelw.score, # Gen 1 initialize, crossover, mutation cleaned
                InterRank = InterRank
    ))
  }
  else{
    return(list(newparents = MatrixD, # Generation 1 parents matrix used for next generation
                InterRank = InterRank
    ))
  }
}
