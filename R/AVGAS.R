#' A Variable selection using Genetic AlgorithmS
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
#' @param r1 A numeric value indicating the maximum number of main effects. This number
#' can be different from the \code{r1} defined in \code{\link{detect}}.
#' @param r2 A numeric value indicating the maximum number of interaction effects. This number
#' can be different from the \code{r1} defined in \code{\link{detect}}.
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
#' @param take Only used when \code{allout = "No"}. Number of top candidate models
#' to display. Default is 3.
#'
#' @return A list of output. The components are:
#' \item{final_model}{The final selected model.}
#' \item{cleaned_candidate_model}{All candidate models where each row corresponding
#' to a fitted model; the first 1 to \code{r1 + r2} columns representing the predictor
#' indices in that model, and the last column is a numeric value representing the
#' ABC score of that fitted model. Duplicated models are not allowed.}
#' \item{InterRank}{Rank of all candidate interaction effects. A two-column numeric
#' matrix. The first column contains indices of ranked two-way interaction effects, and the
#' second column contains its corresponding ABC score.}
#' @export
#' @seealso \code{\link{initial}}, \code{\link{cross}}, \code{\link{mut}}, \code{\link{ABC}}, \code{\link{Genone}}, and \code{\link{Extract}}.
#' @importFrom utils combn
#' @importFrom selectiveInference estimateSigma
#' @importFrom Matrix rankMatrix
#' @importFrom pracma orth
#' @importFrom stats rnorm
#' @importFrom stats reorder
#' @importFrom VariableScreening screenIID
#' @importFrom stats na.exclude
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 labs
#' @importFrom dplyr distinct
#' @importFrom stats na.omit
#'
#' @examples # allout = "No"
#'
# set.seed(0)
# nmain.p <- 4
# interaction.ind <- t(combn(4,2))
# X <- matrix(rnorm(50*4,1,0.1), 50, 4)
# epl <- rnorm(50,0,0.01)
# y <- 1+X[,1]+X[,2]+X[,1]*X[,2]+epl
#
# a1 <- AVGAS(X, y, nmain.p=4, r1=3, r2=3,
#     interaction.ind = interaction.ind, q=5)
#'
#' @examples # allout = "Yes"
# a2 <- AVGAS(X, y, nmain.p=4, r1=3, r2=3,
#     interaction.ind = interaction.ind, q=5, allout = "Yes")
#'
AVGAS <- function(X, y, heredity = "Strong", nmain.p, r1, r2, sigma = NULL,
                  interaction.ind = NULL, lambda = 10, q = 40, allout = "No",
                  interonly = "No",  pi1 = 0.32, pi2 = 0.32, pi3 = 0.32,
                  aprob = 0.9, dprob = 0.9, aprobm = 0.1, aprobi=0.9,
                  dprobm = 0.9, dprobi = 0.1, take = 3){

  first <- Genone(X, y, heredity = heredity, nmain.p = nmain.p, r1 = r1, r2 = r2, sigma = sigma,
                  interaction.ind = interaction.ind, lambda = lambda, q = q, allout = "Yes",
                  interonly = interonly, pi1 = pi1, pi2 = pi2, pi3 = pi3, aprob = aprob, dprob = dprob,
                  aprobm=aprobm, aprobi = aprobi, dprobm = dprobm, dprobi = dprobi)
  parents <- first
  parents$initialize <- first$newparents
  oldABCscore <- as.numeric(first$parents_models_cleaned[1,(r1+r2)+1])
  InterRank = first$InterRank
  count <- 0
  repeat{
    count <- count + 1
    NewMatrixB <- cross(parents, heredity = heredity,
                        nmain.p = nmain.p, r1 = r1, r2 = r2, interaction.ind = interaction.ind)
    NewmatrixC <- mut(parents, heredity = heredity, nmain.p = nmain.p, r1 = r1, r2 = r2,
                      interaction.ind = interaction.ind,
                      interonly = interonly, aprob = aprob, dprob = dprob,
                      aprobm = aprobm, aprobi=aprobi, dprobm = dprobm, dprobi = dprobi)
    MatrixE <- rbind(NewMatrixB, NewmatrixC)

    if (interonly == "Yes"){
      MatrixEF <- MatrixE
      for (i in 1:dim(MatrixEF)[1]) {
        if (length(which(MatrixEF[i,]%in% 1:nmain.p))>0){
          tempind <- which(MatrixEF[i,]%in% 1:nmain.p)
          MatrixEF[i,tempind] <- 0
        }
      }
      MatrixE <- MatrixEF
    }else{
      MatrixE <- MatrixE
    }
    NewABCscore <- list()
    for (i in 1:dim(MatrixE)[1]){
      NewABCscore[[i]] <-  ABC(X, y, heredity = heredity, nmain.p = nmain.p, sigma = sigma,
                               extract = "Yes", varind = c(as.numeric(MatrixE[i,][which(!MatrixE[i,]==0)])),
                               interaction.ind = interaction.ind,
                               pi1 = pi1, pi2 = pi2, pi3= pi3, lambda = lambda)
    }
    temp2 <- as.matrix(cbind(MatrixE, New.scores = unlist(NewABCscore)))
    temp3 <- as.matrix(rbind(first$parents_models_cleaned, temp2))
    new.modelw.score <- temp3[order(unlist(temp3[,((r1+r2)+1)])),]

    new.modelw.score <- cbind(t(apply(t(apply(new.modelw.score[,-((r1+r2)+1)], 1, function(x){x <- sort(x, decreasing = TRUE);x})), 1, function(x) {x[x != 0] <- sort(x[x != 0]);x})),new.modelw.score[,(r1+r2)+1])

    new.modelw.score <- unique_rows(new.modelw.score)
    if (nrow(new.modelw.score) < q){
      MatrixD <- new.modelw.score[,1:(r1+r2)]
      MatrixD <- rbind(MatrixD, matrix(0, nrow = q - nrow(MatrixD), ncol = (r1+r2)))
    }else{
      MatrixD <- new.modelw.score[1:q, 1:(r1+r2)]
    }
    NewABCscore <- as.numeric(new.modelw.score[1,(r1+r2)+1])
    if (round(NewABCscore,10) >= round(oldABCscore,10)){
      break
    }else{
      oldABCscore <- NewABCscore
      parents <- first
      parents$initialize <- MatrixD
    }
  }
  selected_model <- as.numeric(MatrixD[1,])
  model.match <- predictor_match(new.modelw.score, r1 = r1, r2 = r2 , nmain.p = nmain.p, interaction.ind = interaction.ind)
  model.match.matrix <- matrix(0, nrow = length(model.match), ncol = (r1+r2))
  for (b in 1: nrow(model.match.matrix)) {
    model.match.matrix[b,1:length(model.match[[b]])] <- model.match[[b]]
  }
  model.match.matrix <- as.matrix(cbind(model.match.matrix, ABCscore = round(unlist(new.modelw.score[,(r1+r2)+1]), 4)))
  if (allout == "Yes"){
    return(list(
      final_model = model.match[1],
      cleaned_candidate_model = model.match.matrix,
      InterRank = InterRank
    ))
  }else{
    return(list(
      final_model = model.match[1],
      cleaned_candidate_model = model.match.matrix[1:take,]
    ))
  }
}

