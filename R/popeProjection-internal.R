getSpeciesInfo <- function(sp, data=NULL){
  otros <- NULL
  out <- if(sp %in% rownames(species)) as.list(species[sp, ]) else otros
  return(out)
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


lengthProjMatrix <- function(sp, k, Linf, freq){

  dt      <- 1/freq
  species <- getSpeciesInfo(sp)
  bin     <- species$bin

  marcas     <- .createMarks(species)
  marcas_phi <- .createMarks(species, phi=TRUE)
  marcas_inf <- marcas - 0.5*bin
  marcas_sup <- marcas + 0.5*bin

  k    <- k*dt
  Linf <- Linf

  l_inf <- brody(marcas_inf, Linf=Linf, k=k)
  l_sup <- brody(marcas_sup, Linf=Linf, k=k)

  newMarcas <- cbind(l_inf, l_sup)
  A <- t(apply(newMarcas, 1, .lengthProj, marcas=marcas_phi))
  return(A)
}


naturalMortality <- function(sp, vectorM, sizeM, freq){

  dt <- 1/freq
  species <- getSpeciesInfo(sp)
  bin <- species$bin

  marcas <- .createMarks(species)

  M_table <- data.frame(size = sizeM, vectorM = vectorM)

  mPos <- findInterval(marcas, M_table$size)

  M <- M_table[mPos, "vectorM"]*dt
  names(M) <- marcas

  return(M)
}


.projectPOPE <- function (N0, catch, a, b, k, Linf, sizeM, vectorM, freq, sp, Ts){

  A <- lengthProjMatrix(sp=sp, k=k, Linf=Linf, freq=freq)
  M <- naturalMortality(sp=sp, sizeM=sizeM, vectorM=vectorM, freq=freq)

  N <- matrix(ncol=length(M), nrow=Ts+1)

  N[1, ] <- N0

  for(t in seq_len(Ts)){
    nNew <- (N[t,]*exp(-(M/2)) - catch)*exp(-(M/2))
    N[t+1, ] <- as.numeric(nNew %*% A)
  }

  species <- getSpeciesInfo(sp)
  marcas <- .createMarks(species)

  weights <- a*marcas^b

  maturity <- .maturity.ojive(sp)

  B   <- N %*% weights
  BD  <- N %*% (maturity*weights)
  BDR <- tail(BD, 1)

  return(list(N=N, C=C, B=B, BD=BD, BDR=BDR))
}


.maturity.ojive = function(sp) {
  species = getSpeciesInfo(sp)
  marcas = .createMarks(species)

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


makeLengthStairs <- function(projectionData, surveyData, catchData, absolute = FALSE,
                             xlim = c(2, 20, 0.5), ylimProj = NULL, ylimCatch = c(0, 50, 10),
                             ylimBiomass = c(0, 10e6, 5e6), ylimBiomassFactor = 1e6,
                             cols = c("green4", "red", "blue3"), ltys = c("solid", "solid", "dotted"),
                             lwds = 2, main = NA, ...){

  require(ruisu)

  thousands <- TRUE

  # Read data base
  if(!any(is.element(c("data.frame", "matrix"), class(projectionData)))){
    projectionData <- read.csv(file = projectionData, stringsAsFactors = FALSE, check.names = FALSE,
                               row.names = 1)
  }

  if(!any(is.element(c("data.frame", "matrix"), class(surveyData)))){
    surveyData <- read.csv(file = surveyData, stringsAsFactors = FALSE, check.names = FALSE,
                           row.names = 1)
  }

  if(!any(is.element(c("data.frame", "matrix"), class(catchData)))){
    catchData <- read.csv(file = catchData, stringsAsFactors = FALSE, check.names = FALSE,
                          row.names = 1)
  }

  projectionData[projectionData < 0.0001] <- NA
  surveyData[surveyData < 0.0001] <- NA
  catchData[catchData < 0.0001] <- NA

  if(isTRUE(absolute)){

    if(is.null(ylimProj)){
      newYlim <- pretty(c(0, apply(projectionData, 2, max, na.rm = TRUE)))
      ylimProj <- c(range(newYlim), unique(diff(newYlim)))
    }

    ylimCatch <- ylimProj
  }else{
    ylimProj <- c(0, 20, 5)
  }

  # Check if ncol of data bases are the same
  if(ncol(projectionData)*ncol(surveyData)*ncol(catchData) != ncol(projectionData)^3){
    stop("All data must have the same number of columns.")
  }

  # Get length marks
  allMarks <- an(rownames(projectionData))

  # In surveyData, keep colnames only for columns where colSums > 0
  index <- colSums(surveyData, na.rm = TRUE) > 0
  colnames(surveyData)[!index] <- ""

  # Vector of three data
  allBases <- c("projectionData", "surveyData", "catchData")

  # Standardize cols, ltys and lwds
  cols <- rep(cols, length.out = 3)
  ltys <- rep(ltys, length.out = 3)
  lwds <- rep(lwds, length.out = 3)

  # Build plot layout
  groupFactor <- ncol(projectionData)
  layoutMatrix <- matrix(rep(x = seq(groupFactor), each = 4), nrow = groupFactor, byrow = TRUE)
  layoutMatrix <- cbind(layoutMatrix, seq(groupFactor + 1, length.out = groupFactor))

  # Set graphic parameters
  x11()

  layout(layoutMatrix)
  par(mar = rep(0, 4), oma = c(5, 5, 4, 5), xaxs = "i", yaxs = "i")

  # Loop for each column (steps of stairs)
  for(j in seq(ncol(projectionData))){

    # Make an empty canvas
    plot(1, 1, pch = NA, axes = FALSE, xlab = NA, ylab = NA, xlim = xlim[1:2], ylim = ylimProj[1:2],
         main = main)

    # Loop for each data base
    for(i in seq_along(allBases)){

      # Temporal data
      tempData <- get(allBases[i])

      # Convert to percentage
      lengthVector <- tempData[,j]

      if(!isTRUE(absolute)){
        lengthVector <- lengthVector/sum(lengthVector, na.rm = TRUE)*100
      }

      # With catchData, it'll use ylimCatch limits
      if(allBases[i] == "catchData"){


        if(isTRUE(absolute)){
          lengthVector <- lengthVector
        }else{
          lengthVector <- lengthVector/ylimCatch[2]*ylimProj[2]
        }
      }

      # Make lines
      lines(allMarks, lengthVector, col = cols[i], lty = ltys[i], lwd = lwds[i])
    }

    # Juvenile line
    abline(v = 12, lty = "dotted", col = "red")

    # Make step labels
    mtext(text = colnames(projectionData)[j], side = 3, line = -1.5, adj = 0.99)
    mtext(text = colnames(surveyData)[j], side = 3, line = -1.5, adj = 0.01)

    # Make X axis
    if(j == ncol(projectionData)){
      axis(side = 1, at = seq(xlim[1], xlim[2], xlim[3]), labels = NA)
      axis(side = 1, at = seq(xlim[1], xlim[2]), tick = FALSE)
    }

    # Make Y axis
    ylimFactor <- ifelse(isTRUE(thousands), 1e3, 1)

    if(j %% 2 > 0){
      axis(side = 2, at = seq(ylimProj[1], ylimProj[2], ylimProj[3]), las = 2,
           labels = seq(ylimProj[1], ylimProj[2], ylimProj[3])/ylimFactor)
    }else{
      axis(side = 4, at = seq(ylimProj[1], ylimProj[2],
                              length.out = (ylimCatch[2] - ylimCatch[1])/ylimCatch[3] + 1),
           las = 2, labels = seq(ylimCatch[1], ylimCatch[2], ylimCatch[3])/ylimFactor)
    }

    box()
  }

  # Vector of three data
  allBases <- c("surveyData", "catchData", "projectionData")

  par(mar = c(0, 4, 0, 0), xaxs = "r")

  for(j in seq(ncol(projectionData))){

    biomassBases <- NULL
    for(i in seq_along(allBases)){
      # Temporal data
      tempData <- get(allBases[i])

      # Calculate biomass
      lengthVector <- (growthParams[j, "a"]*an(rownames(tempData))^growthParams[j, "b"])*tempData[,j]
      biomassBases <- c(biomassBases, sum(lengthVector, na.rm = TRUE))
    }

    if(j == ncol(projectionData)){
      names.arg <- c("Crucero", "Capturas", "Projección")
    }else{
      names.arg <- rep(NA, 3)
    }

    barplot(biomassBases, names.arg = names.arg, col = cols, ylim = ylimBiomass[1:2], axes = FALSE, las = 3)

    xAxisLab <- seq(ylimBiomass[1], ylimBiomass[2], ylimBiomass[3])/ylimBiomassFactor
    xAxisLab[1] <- NA
    axis(side = 4, at = seq(ylimBiomass[1], ylimBiomass[2], ylimBiomass[3]), labels = xAxisLab, las = 2)
    box()

  }

  # X and Y axis text
  mtext(text = paste0(ifelse(isTRUE(absolute), "", "Frecuencia ("),
                      ifelse(isTRUE(absolute), "Millones individuos", "%"),
                      ifelse(isTRUE(absolute), "", ")")),
        side = 2, outer = TRUE, line = 3)
  mtext(text = "Longitud total (cm)", side = 1, outer = TRUE, line = 3)
  mtext(text = "Biomasa (millones de t)", side = 4, outer = TRUE, line = 3)

  return(invisible())
}


# Inverse functions -------------------------------------------------------

brodyInverse <- function(l, Linf, k){
  out <- (l - Linf*(1-exp(-k)))/exp(-k)
  out[out <= 2] <- 2 ####
  return(out)
}


.lengthProjInverse <- function(l, marcas){

  x     <- sort(c(l, marcas))
  dif   <- diff(x)
  pos   <- findInterval(l, marcas) + 1
  pos   <- seq(from=pos[1], to=pos[2])
  props <- dif[pos]
  props <- props/sum(props)
  out   <- numeric(length(marcas)-1)
  out[pos-1] <- props
  return(out)
}


lengthProjMatrixInverse <- function(sp, k, Linf, freq){

  dt      <- 1/freq
  species <- getSpeciesInfo(sp)
  bin     <- species$bin

  marcas     <- .createMarks(species)
  marcas_phi <- .createMarks(species, phi=TRUE)
  marcas_inf <- marcas - 0.5*bin
  marcas_sup <- marcas + 0.5*bin

  k    <- k*dt
  Linf <- Linf

  l_inf <- brodyInverse(marcas_inf, Linf=Linf, k=k)
  l_sup <- brodyInverse(marcas_sup, Linf=Linf, k=k)
  l_sup[length(l_sup)] <- rev(tail(marcas))[1] ###

  newMarcas <- cbind(l_inf, l_sup)
  A <- t(apply(newMarcas, 1, .lengthProjInverse, marcas=marcas_phi))
  A[is.na(A)] <- 0
  return(A)
}


naturalMortalityInverse <- function(sp, vectorM, sizeM, freq){

  dt <- 1/freq
  species <- getSpeciesInfo(sp)
  bin <- species$bin

  marcas  <- .createMarks(species)

  M_table <- data.frame(size = sizeM, vectorM = vectorM)

  mPos    <- findInterval(marcas, M_table$size)

  M <- M_table[mPos, "vectorM"]*dt
  names(M) <- marcas

  return(M)
}


.projectPOPEInverse <- function(N0, catch, a, b, k, Linf, sizeM, vectorM, freq, sp, Ts){

  A <- lengthProjMatrixInverse(sp = sp, k = k, Linf = Linf, freq = freq)
  M <- naturalMortalityInverse(sp = sp, sizeM = sizeM, vectorM = vectorM, freq = freq)

  N <- matrix(ncol = length(M), nrow = Ts + 1)

  N[1, ] <- N0

  for(t in seq_len(Ts)){
    nNew <- (N[t,] + catch*exp(-M/2))/exp(-M)
    N[t+1, ] <- as.numeric(nNew %*% A)
  }

  species <- getSpeciesInfo(sp)
  marcas  <- .createMarks(species)

  weights <- a*marcas^b

  maturity <- .maturity.ojive(sp)

  B   <- N %*% weights
  BD  <- N %*% (maturity*weights)
  BDR <- tail(BD, 1)

  return(list(N = N, C = C, B = B, BD = BD, BDR = BDR))
}
