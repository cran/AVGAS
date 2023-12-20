predictor_match <- function(candidate.model, r1,r2, nmain.p, interaction.ind){
  model.match <- list()
  if (mydim(candidate.model)[1]==1){
    ee <- candidate.model[1:(r1+r2)]
    ee1 <- ee[which(ee%in% 1:nmain.p)]
    ee2 <- ee[which(ee > nmain.p)]
    if (length(ee2) == 0){
      model.match[[1]] <-c(paste0("X.", ee1))
    }
    if (length(ee1) == 0){
      model.match[[1]] <- c( paste0("X.", interaction.ind[ee2-nmain.p,1], "X.", interaction.ind[ee2-nmain.p,2]))
    }
    if (!length(ee2) ==0 & !length(ee1) ==0){
      model.match[[1]] <-c(paste0("X.", ee1), paste0("X.", interaction.ind[ee2-nmain.p,1], "X.", interaction.ind[ee2-nmain.p,2]))
    }
  }else{
    for (i in 1:nrow(candidate.model)) {
      ee <- candidate.model[i, 1:(r1+r2)]
      ee1 <- ee[which(ee%in% 1:nmain.p)]
      ee2 <- ee[which(ee > nmain.p)]
      if (length(ee2) == 0){
        model.match[[i]] <-c(paste0("X.", ee1))
      }
      if (length(ee1) == 0){
        model.match[[i]] <- c( paste0("X.", interaction.ind[ee2-nmain.p,1], "X.", interaction.ind[ee2-nmain.p,2]))
      }
      if (!length(ee2) ==0 & !length(ee1) ==0){
        model.match[[i]] <-c(paste0("X.", ee1), paste0("X.", interaction.ind[ee2-nmain.p,1], "X.", interaction.ind[ee2-nmain.p,2]))
      }
    }
  }
  return(model.match = model.match)
}
