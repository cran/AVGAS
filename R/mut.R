#' Performing mutation\cr
#'
#' This function performs mutation which only stores all fitted models without
#' making any comparison. The selected indices in each fitted model will be
#' automatically re-ordered so that main effects comes first, followed by
#' two-way interaction effects, and zero reservation spaces.
#'
#' @param parents  A numeric matrix of dimension \code{q} by \code{r1+r2},
#' obtained from \code{initial} or previous generation where each row corresponding
#' a fitted model and each column representing the predictor index in the fitted model.
#' @param heredity Whether to enforce Strong, Weak, or No heredity. Default is "Strong".
#' @param nmain.p A numeric value that represents the total number of main effects
#' in \code{X}.
#' @param r1 A numeric value indicating the maximum number of main effects.
#' @param r2 A numeric value indicating the maximum number of interaction effects.
#' @param interaction.ind A two-column numeric matrix containing all possible
#' two-way interaction effects. It must be generated outside of this function
#' using \code{t(utils::combn())}. See Example section for details.
#' @param interonly Whether or not to consider fitted models with only two-way
#' interaction effects. A “Yes" or "No" logical vector. Default is "No".
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
#' @return  A numeric matrix \code{single.child.mutated} is returned. Each row representing
#' a fitted model, and each column corresponding to the predictor index in the fitted model.
#' Duplicated models are allowed.
#' @export
#' @seealso \code{\link{initial}}.
#'
#' @examples # Under Strong heredity, interonly = "No"
#' set.seed(0)
#' nmain.p <- 4
#' interaction.ind <- t(combn(4,2))
#' X <- matrix(rnorm(50*4,1,0.1), 50, 4)
#' epl <- rnorm(50,0,0.01)
#' y <- 1+X[,1]+X[,2]+X[,1]*X[,2]+epl
#' p1 <- initial(X, y, nmain.p = 4, r1 = 3, r2 = 3,
#'     interaction.ind = interaction.ind, q = 5)
#' m1 <- mut(p1, nmain.p = 4, r1 = 3, r2 = 3,
#'     interaction.ind =interaction.ind)
#' @examples # Under Strong heredity, interonly = "Yes"
#' m2 <- mut(p1, heredity = "No", nmain.p = 4, r1 = 3, r2 = 3,
#'     interaction.ind =interaction.ind, interonly = "Yes")

mut <- function(parents, heredity = "Strong", nmain.p,
                r1, r2, interaction.ind = NULL, interonly = "No",
                aprob = 0.9, dprob = 0.9, aprobm = 0.1, aprobi=0.9, dprobm = 0.9, dprobi = 0.1){
  if (is.null(interaction.ind)) stop("Interaction.ind is missing.
                                       Use t(utils::combn()) to generate interaction matrix.")

  single.child.mutated <- parents$initialize

  for (i in 1: nrow(single.child.mutated)) {
    if (length(single.child.mutated[i,][which(single.child.mutated[i,] == 0)])>0){
      addition <- stats::rbinom(1, 1, prob = aprob)
      if (as.logical(addition)){
        additionindpool <- setdiff(union(as.numeric(parents$mainpool),as.numeric(parents$InterRank[,1]) ),
                                   single.child.mutated[i,][which(!(single.child.mutated[i,]) == 0)])
        additionindpool.main <- additionindpool[additionindpool%in% 1:nmain.p]
        additionindpool.inter <- additionindpool[additionindpool>nmain.p]

        aamain <- length(single.child.mutated[i,][which(single.child.mutated[i,]%in%1:nmain.p)])
        aainter <- length(single.child.mutated[i,][which(single.child.mutated[i,]>nmain.p)])

        if (!length(additionindpool.inter)==0){
          additionind <- mysample(stats::na.omit(c(additionindpool.main[1], additionindpool.inter[1])), 1 , prob=c(aprobm,aprobi))
          if (additionind<=nmain.p & !aamain==0 & aamain<r1){
            additionind <- mysample(additionindpool.main[1],1)
            single.child.mutated[i,max(which(!single.child.mutated[i,] == 0))+1] <- additionind
          }else{
            single.child.mutated[i,] <- single.child.mutated[i,]
          }
          if (additionind>nmain.p & !aainter==0 & aainter<r2){
            if (heredity == "Strong" | heredity == "Weak"){
              check <- Heredity(x = c(single.child.mutated[i,][which(!single.child.mutated[i,]==0)],additionind),
                                nmain.p = nmain.p, interaction.ind = interaction.ind, heredity = heredity)
              if (check == TRUE){
                single.child.mutated[i,max(which(!single.child.mutated[i,] == 0))+1] <- additionind
              }
            }
          }else{
            single.child.mutated[i,] <- single.child.mutated[i,]
          }
        }
      }
      else{
        single.child.mutated[i,] <- single.child.mutated[i,]
      }
    }

    deletion <- stats::rbinom(1, 1, prob = dprob)
    if (as.logical(deletion)){
      if (sum(!single.child.mutated[i,]==0)>1){
        bbb <- as.numeric(single.child.mutated[i,][which(single.child.mutated[i,]%in% 1:nmain.p)])
        ccc <- as.numeric(single.child.mutated[i,][which(single.child.mutated[i,]>nmain.p)])

        dmain <- stats::rbinom(1, 1, prob = dprobm)
        if (as.logical(dmain) & !length(bbb)==0){
          sample_index <- as.numeric(mysample(bbb,1))
          deletionind <- sample_index
        }
        dinter <- stats::rbinom(1, 1, prob = dprobi)
        if (dmain == FALSE & as.logical(dinter) & !length(dinter)==0){
          sample_index <- as.numeric(mysample(ccc,1))
          deletionind <- sample_index
        }

        if (dmain == FALSE & dinter == FALSE){
          single.child.mutated[i,] <- single.child.mutated[i,]
        }else{
          if (heredity == "No"){
            if (interonly == "Yes"){
              single.child.mutated[i,] <- replace(single.child.mutated[i,], which(single.child.mutated[i,] < nmain.p+1), 0)
            }else{
              single.child.mutated[i,] <- replace(single.child.mutated[i,], which(single.child.mutated[i,] == deletionind), 0)
            }
          }

          if (heredity == "Strong"){
            if (deletionind %in% 1:nmain.p){
              mutate.inter <- single.child.mutated[i,][single.child.mutated[i,] >nmain.p]
              for (j in 1:length(mutate.inter)) {
                if (any(interaction.ind[mutate.inter[j]-nmain.p,] %in% deletionind)){
                  single.child.mutated[i,] <- replace(single.child.mutated[i,], which(single.child.mutated[i,] == mutate.inter[j]), 0)
                }
              }
            }
            single.child.mutated[i,] <- replace(single.child.mutated[i,], which(single.child.mutated[i,] == deletionind), 0)
          }

          if (heredity =="Weak"){
            if (deletionind %in% 1:nmain.p){
              mutate.inter <- single.child.mutated[i,][single.child.mutated[i,] > nmain.p]
              for (j in 1:length(mutate.inter)) {
                if (!any(interaction.ind[mutate.inter[j]-nmain.p,]%in% setdiff(single.child.mutated[i,],deletionind))){
                  single.child.mutated[i,] <- replace(single.child.mutated[i,], which(single.child.mutated[i,] == mutate.inter[j]), 0)
                }
              }
            }
            single.child.mutated[i,] <- replace(single.child.mutated[i,], which(single.child.mutated[i,] == deletionind), 0)
          }
        }
      }
    }
    else{
      single.child.mutated[i,] <- single.child.mutated[i,]
    }
  }
  single.child.mutated <- as.matrix(single.child.mutated)
  for (i in 1:nrow(single.child.mutated)) {
    single.child.mutated[i,] <- sort_zeros(single.child.mutated[i,])
  }
  return(single.child.mutated)
}
