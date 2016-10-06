getSpeciesInfo <- function(sp, data=NULL){
  otros <- NULL
  out <- if(sp %in% rownames(species)) as.list(species[sp, ]) else otros
  return(out)
}


.getScenarios <- function(scenario, n){

  scenarios <- c("neutro", "desfavorable", "favorable")

  if(is.character(scenario)){
    if(length(scenario)!=1) stop("You can only provide one scenario")
    return(rep(scenario, n))
  }

  if(length(scenario)!=length(scenarios))
    stop("You must indicate one probability per scenario.")
  set.seed(880820)
  scenario <- sample(x=scenarios, size=n, prob = scenario, replace=TRUE)

  return(scenario)
}


.createMarks <- function(specie, phi=FALSE){
  marks <- seq(from=specie$Lmin, to=specie$Lmax, by=specie$bin)
  if(isTRUE(phi)){
    marks_inf <- marks - 0.5*specie$bin
    marks_sup <- marks + 0.5*specie$bin
    marks <- sort(unique(c(marks_inf, marks_sup)))
  }
  return(marks)
}


brody <- function(l, Linf, k){
  out <- Linf - (Linf-l)*exp(-k)
  return(out)
}


.lengthProj <- function(l, marcas){

  x <- sort(c(l, marcas))
  dif <- diff(x)
  pos <- findInterval(l, marcas) + 1
  pos <- seq(from=pos[1], to=pos[2])
  props <- dif[pos]
  props <- props/sum(props)
  out <- numeric(length(marcas)-1)
  out[pos-1] <- props
  return(out)
}


lengthProjMatrix <- function(sp, scenario, freq){

  if(!(scenario %in% c("neutro", "favorable", "desfavorable")))
    stop("Wrong scenario.")

  dt <- 1/freq
  species <- getSpeciesInfo(sp)
  bin <- species$bin

  marcas <- .createMarks(species)
  marcas_phi <- .createMarks(species, phi=TRUE)
  marcas_inf <- marcas - 0.5*bin
  marcas_sup <- marcas + 0.5*bin

  k    <- LHT[[sp]]$G["k", scenario]*dt
  Linf <- LHT[[sp]]$G["Linf", scenario]

  l_inf <- brody(marcas_inf, Linf=Linf, k=k)
  l_sup <- brody(marcas_sup, Linf=Linf, k=k)

  newMarcas <- cbind(l_inf, l_sup)
  A <- t(apply(newMarcas, 1, .lengthProj, marcas=marcas_phi))
  return(A)
}


naturalMortality <- function(sp, scenario, freq){

  if(!(scenario %in% c("neutro", "favorable", "desfavorable")))
    stop("Wrong scenario.")
  dt <- 1/freq
  species <- getSpeciesInfo(sp)
  bin <- species$bin

  marcas <- .createMarks(species)

  M_table <- LHT[[sp]]$M

  mPos <- findInterval(marcas, M_table$size)

  M <- M_table[mPos, scenario]*dt
  names(M) <- marcas

  return(M)
}


.projectPOPE <- function (N0, catch, scenario, freq, sp, Ts){

  A <- lengthProjMatrix(sp=sp, scenario=scenario, freq=freq)
  M <- naturalMortality(sp=sp, scenario=scenario, freq=freq)

  N <- matrix(ncol=length(M), nrow=Ts+1)

  N[1, ] <- N0

  for(t in seq_len(Ts)){
    nNew <- (N[t,]*exp(-(M/2)) - catch)*exp(-(M/2))
    N[t+1, ] <- as.numeric(nNew %*% A)
  }

  species <- getSpeciesInfo(sp)
  marcas <- .createMarks(species)

  weights <- species$a*marcas^species$b

  maturity <- .maturity.ojive(sp)

  B   <- N %*% weights
  BD  <- N %*% (maturity*weights)
  BDR <- tail(BD, 1)

  return(list(N=N, C=C, B=B, BD=BD, BDR=BDR))
}

readAtLength = function(file, sp = "anchoveta"){
  if(is.null(file)) return(NULL)

  base =  read.csv(file, stringsAsFactors = FALSE)
  colnames(base) = tolower(colnames(base))
  specie = getSpeciesInfo(sp)
  marcas = .createMarks(specie)
  newBase = expand.grid(x = marcas)
  base = merge(base, newBase, all = T)
  base = as.matrix(base[,-1])
  base[is.na(base)] = 0

  return(base)
}


.maturity.ojive = function(sp) {

  species = getSpeciesInfo(sp)
  marcas = .createMarks(species)

  #Use "mat1" and "mat2" values of the current survey.
  #Change the values in the species.csv (auxiliar folder,
  #then run the script species.R and build & reload the package.)
  out = 1/(1+exp(species$mat1+species$mat2*marcas))
  return(out)

}

getDates <- function(x){
  timeVector <- sapply(sub("^([[:alpha:]]*).*", "\\1", x), grep, x = month.abb_spanish)
  timeVector <- paste(15, timeVector, sep = "-",
                      an(gsub(pattern = "[^\\d]+", replacement = "", x = x, perl = TRUE)))
  timeVector <- as.Date(timeVector, format = "%d-%m-%y")

  return(timeVector)
}
