#' Performing crossover\cr
#'
#' This function performs crossover which only stores all fitted models without
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
#'
#' @return A numeric matrix \code{single.child.bit} is returned. Each row representing
#' a fitted model, and each column corresponding to the predictor index in the fitted model.
#' Duplicated models are allowed.
#' @export
#' @seealso \code{\link{initial}}.
#'
#' @examples # Under Strong heredity
#' set.seed(0)
#' nmain.p <- 4
#' interaction.ind <- t(combn(4,2))
#' X <- matrix(rnorm(50*4,1,0.1), 50, 4)
#' epl <- rnorm(50,0,0.01)
#' y<- 1+X[,1]+X[,2]+X[,1]*X[,2]+epl
#' p1 <- initial(X, y, nmain.p = 4, r1 = 3, r2 = 3,
#'     interaction.ind = interaction.ind, q = 5)
#' c1 <- cross(p1, nmain.p=4, r1 = 3, r2 = 3,
#'     interaction.ind = interaction.ind)

cross <- function(parents, heredity = "Strong", nmain.p, r1, r2, interaction.ind = NULL){
  if (is.null(interaction.ind)) stop("Interaction.ind is missing.
                                       Use t(utils::combn()) to generate interaction matrix.")
  max_model_size  <- length(parents$initialize[1,])
  parentsMB <- parents$InterRank[,1]

  single.child.bit <- matrix(0,nrow=choose(dim(parents$initialize)[1],2),ncol=max_model_size)
  tempcount <- 0
  for (i in 1:(dim(parents$initialize)[1]-1)) {
    for (j in ((i+1):(dim(parents$initialize)[1]))) {
      tempcount <- tempcount + 1
      crossind <- union(parents$initialize[i,][which(!parents$initialize[i,]==0)],
                        parents$initialize[j,][which(!parents$initialize[j,]==0)])
      crossind <- as.numeric(unlist(crossind))

      crossindmain <- as.numeric(unique(crossind[which(crossind%in% 1:nmain.p)]))
      crossindinter <- as.numeric(unique(crossind[which(crossind>nmain.p)]))

      if (length(crossindmain)<=r1 & length(crossindinter)<=r2){
        if (length(crossindmain)>0){
          single.child.bit[tempcount, c(1:length(crossindmain))] <- crossindmain
        }else{
          single.child.bit[tempcount, c(1:r1)] <- 0
        }
        if (length(crossindinter)>0){
          single.child.bit[tempcount,
                           c((max(which(!single.child.bit[tempcount,] == 0))+1):((max(which(!single.child.bit[tempcount,] == 0))+1)+length(crossindinter)-1))] <- crossindinter
        }else{
          single.child.bit[tempcount, c((length(crossindmain)+1):max_model_size)] <- 0
        }
      }else{
        single.child.bit[tempcount, c(1:min(r1, length(crossindmain)))] <- sort(mysample(crossindmain, min(r1,length(crossindmain))))
        single.child.bit[tempcount, (max(which(!single.child.bit[tempcount,] == 0))+1)] <- as.numeric(unlist(parentsMB))[1]
      }
    }
  }
  single.child.bit <- as.matrix(single.child.bit)
  return(single.child.bit)
}


