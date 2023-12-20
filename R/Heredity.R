Heredity <- function(x, nmain.p, interaction.ind, heredity = "Strong"){
  if (is.null(heredity)) stop("You must define a heredity condition")
  intereffect <- x[which(x>nmain.p)]
  if (length(intereffect)== 0) stop("This model contains no interaction effect")
  maineffect <- x[which(x%in% 1:nmain.p)]
  check <- c()
  if (heredity == "Strong"){
    for (i in 1:length(intereffect)) {
      check <- c(check, interaction.ind[as.numeric(intereffect[i]-nmain.p),])
    }
    result <- all(check%in% maineffect)
  }
  else if (heredity == "Weak"){
    for (i in 1:length(intereffect)) {
      ee <- interaction.ind[as.numeric(intereffect[i]-nmain.p),]
      check[i] <- any(ee%in% maineffect)
    }
    result <- all(check)
  }
  return(result)
}
