PopeProjection <- function(fileAbundancia, fileCaptura, scenario, a, b, k, Linf, vectorM, freq, sp, Ts){
  abundancia <- readAtLength(fileAbundancia, sp = sp)

  if(ncol(abundancia) == 1){
    abundancia <- matrix(rep(abundancia, each = 2), ncol = 2, byrow = TRUE)
  }

  captura <- readAtLength(fileCaptura, sp = sp)

  species[sp, "a"] <- a
  species[sp, "b"] <- b
  LHT$anchoveta$G[scenario] <- c(k, Linf, -0.210)
  LHT$anchoveta$M[scenario] <- vectorM

  Nsim <- array(dim = c(ncol(captura)+1, dim(abundancia)))
  Nsim[1, , ] <- abundancia
  Bsim <- array(dim=c(ncol(captura), ncol(abundancia)))

  for(i in 1:ncol(captura)){
    sim <- projectPOPE(N = Nsim[i, ,], captura[, i], scenario = scenario, freq = freq, sp = sp, Ts = Ts)
    Nsim[i+1, , ] <- sim$raw$N[Ts+1, , ]
    Bsim[i, ]     <- sim$raw$B[Ts+1, ]
  }
  raw <- list(Abundancia = Nsim, Biomasa = Bsim)

  return(list(Abundancia = apply(Nsim, c(1,2), median), Biomasa = apply(Bsim, 1, median), raw = raw))
}

an <- as.numeric
ac <- as.character
anc <- function(x, ...) as.numeric(as.character(x, ...))
