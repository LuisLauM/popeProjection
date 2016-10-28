#' Title
#'
#' @param allCatches
#' @param allSurveys
#' @param surveyNames
#' @param allProjections
#' @param catchDataFactor
#'
#' @return
#' @export
#'
#' @examples
autoProjector <- function(allCatches, allSurveys, surveyNames, growthParams, allProjections = NULL,
                          catchDataFactor = 1e-3, sp = "anchoveta",
                          mortalityInfo = list(sizeM = c(0, 8, 12), vectorM = c(1.29, 0.92, 0.83)),
                          freq = 12, Ts = 1){

  myFolder <- dirname(allCatches)

  allCatches <- readAtLength(file = allCatches, sp = sp, check.names = FALSE)
  allSurveys <- read.csv(file = allSurveys, stringsAsFactors = FALSE,
                         check.names = FALSE, row.names = 1, na.strings = 0)

  allCatches <- allCatches*catchDataFactor

  allCatches[is.na(allCatches)] <- 0
  allSurveys[is.na(allSurveys)] <- 0

  if(length(surveyNames) != ncol(allSurveys)){
    stop("Lengths of surveyNames and ncol(allSurveys) must be the same.")
  }

  spInfo <- getSpeciesInfo(sp = sp)
  allMarks <- seq(spInfo$Lmin, spInfo$Lmax, spInfo$bin)

  catchDates <- getDates(colnames(allCatches))
  surveyDates <- getDates(colnames(allSurveys))

  if(is.null(allProjections)){
    cat("\nCreating projection table...\n")

    allProjections <- matrix(data = NA, nrow = length(allMarks), ncol = ncol(allCatches),
                             dimnames = list(allMarks,
                                             paste0(substr(as.yearmon(catchDates), 1, 3), "-",
                                                    substr(as.yearmon(catchDates), 8, 9))))

    write.csv(x = allProjections, file = file.path(myFolder, "allProjections.csv"), na = "")

    allProjections <- file.path(myFolder, "allProjections.csv")

    cat(paste("\n'allProjections' file created in", allProjections, "\n"))
  }

  allProjections <- read.csv(file = allProjections, stringsAsFactors = FALSE,
                             check.names = FALSE, row.names = 1)

  projDates <- getDates(colnames(allProjections))

  growthParams <- read.csv(file = growthParams, stringsAsFactors = FALSE, row.names = 1)

  index <- as.numeric(na.omit(match(c(1, colnames(allProjections), 2), rownames(growthParams))))
  growthParams <- growthParams[index,]

  simulation_biomass <- NULL
  for(i in seq(ncol(allProjections) - 1)){

    if(is.element(projDates[i], surveyDates)){
      seedLengthVector <- allSurveys[,match(projDates[i], surveyDates)]

      if(i == 1){
        allProjections[,i] <- seedLengthVector
      }
    }else{
      seedLengthVector <- allProjections[,i]
    }

    tempLengths <- cbind(seedLengthVector, seedLengthVector)
    tempCatches <- allCatches[,i]

    simulation <- projectPOPE(N = tempLengths, catch = tempCatches,
                              a = growthParams$a[i], b = growthParams$b[i],
                              k = growthParams$k[i], Linf = growthParams$Linf[i],
                              sizeM = mortalityInfo$sizeM, vectorM = mortalityInfo$vectorM,
                              freq = freq, sp = sp, Ts = Ts)
    simulation_abundance <- simulation$raw$N[Ts + 1,,][,1]
    simulation_biomass <- c(simulation_biomass, simulation$raw$B[Ts + 1, 1])

    simulation_abundance[simulation_abundance <= 0] <- 0
    allProjections[,i + 1] <- simulation_abundance

  }

  return(list(capturas = allCatches,
              cruceros = allSurveys,
              proyecciones = round(allProjections, 0),
              crecimiento = growthParams,
              nombres_cruceros = surveyNames,
              biomasas = simulation_biomass))
}



#' Title
#'
#' @param N
#' @param catch
#' @param a
#' @param b
#' @param k
#' @param Linf
#' @param sizeM
#' @param vectorM
#' @param freq
#' @param sp
#' @param Ts
#'
#' @return
#' @export
#'
#' @examples
projectPOPE <- function(N, catch, a, b, k, Linf, sizeM, vectorM, freq, sp, Ts){

  matrixN    <- array(dim=c(Ts+1, dim(N)))
  matrixB    <- matrix(nrow=Ts+1, ncol=ncol(N))
  matrixBD   <- matrix(nrow=Ts+1, ncol=ncol(N))
  matrixBDR  <- numeric(ncol(N))

  for(i in seq_len(ncol(N))){

    N0 <- N[, i]
    sim <- .projectPOPE(N=N0, catch=catch, a=a, b=b, k=k, Linf=Linf, sizeM=sizeM, vectorM=vectorM, freq=freq, sp=sp, Ts=Ts)

    matrixN[, ,i] <- sim$N
    matrixB[, i]  <- sim$B
    matrixBD[, i] <- sim$BD
    matrixBDR[i]  <- sim$BDR
  }

  N    <- apply(matrixN, c(1,2), median)
  B    <- apply(matrixB, 1, median)
  BD   <- apply(matrixBD, 1, median)
  BDR  <- median(matrixBDR)

  rawData <- list(N=matrixN, B=matrixB, BD=matrixBD, BDR=matrixBDR)

  output <- list(N=N, B=B, BD=BD, BDR=BDR,
                 raw=rawData)

  attr(output, which="sp") <-  sp
  attr(output, which="freq") <-  freq
  attr(output, which="Ts") <-  Ts

  class(output) <- "surveyProj"

  return(output)
}


#' Title
#'
#' @param allData
#' @param rangeDate
#' @param nStairs
#' @param outputDir
#' @param outputPattern
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
autoProjectorPlot <- function(allData, rangeDate, nStairs, addBiomassBar = TRUE, absolute = FALSE,
                              outputDir = "./", outputPattern = "Fig_",
                              width = 1200, height = nStairs*120 + ifelse(isTRUE(addBiomassBar), 60, 0), res = 180,
                              ...){

  allProjections <- allData$proyecciones
  allSurveys <- allData$cruceros
  allCatches <- allData$capturas

  index <- min(unique(an(gsub(x = rangeDate, pattern = "[[:alpha:] \\& [:punct:]]", replacement = ""))))
  allDates <- expand.grid(month.abb_spanish, seq(index, index + 50))
  allDates <- apply(allDates, 1, paste, collapse = "-")

  if(diff(match(rangeDate, allDates)) < nStairs){
    index <- seq(match(rangeDate[1], allDates), match(rangeDate[2], allDates))
  }else{
    index <- seq(match(rangeDate[1], allDates),
                 length.out = ceiling(diff(match(rangeDate, allDates))/nStairs)*nStairs)
  }

  allDates <- allDates[index]

  for(i in seq(length(allDates)/nStairs)){

    tempDates <- allDates[seq((nStairs - 1)*(i - 1) + 1, length.out = nStairs)]

    projectionData <- surveyData <- catchData <- matrix(data = NA, nrow = nrow(allData$proyecciones), ncol = length(tempDates),
                                                        dimnames = list(dimnames(allData$proyecciones)[[1]], tempDates))

    index0 <- match(rangeDate, colnames(allData$proyecciones))
    index0 <- colnames(allData$proyecciones)[seq(index0[1], index0[2])]

    index1 <- an(na.omit(match(index0, tempDates)))
    index2 <- an(na.omit(match(tempDates, index0)))
    projectionData[,index1] <- as.matrix(allData$proyecciones[,index2])

    index1 <- an(na.omit(match(colnames(allData$capturas), tempDates)))
    index2 <- an(na.omit(match(tempDates, colnames(allData$capturas))))
    catchData[,index1] <- as.matrix(allData$capturas[,index2])

    index1 <- an(na.omit(match(colnames(allData$cruceros), tempDates)))
    index2 <- an(na.omit(match(tempDates, colnames(allData$cruceros))))
    surveyData[,index1] <- as.matrix(allData$cruceros[,index2])

    index <- an(na.omit(match(rownames(allData$crecimiento), tempDates)))
    growthParams <- allData$crecimiento[index,]

    colnames(surveyData)[index1] <- allData$nombres_cruceros[index2]

    makeLengthStairs(projectionData = projectionData, surveyData = surveyData, catchData = catchData,
                     growthParams = growthParams, addBiomassBar = addBiomassBar, absolute = absolute, ...)

    filename <- paste0(outputDir, outputPattern, tempDates[1], "-", tempDates[length(tempDates)], ".png")

    dev.copy(png, filename = filename, width = width, height = height, res = res)
    dev.off()
  }

  return(invisible())
}


#' Title
#'
#' @param N
#' @param catch
#' @param a
#' @param b
#' @param k
#' @param Linf
#' @param sizeM
#' @param vectorM
#' @param freq
#' @param sp
#' @param Ts
#'
#' @return
#' @export
#'
#' @examples
projectPOPEInverse <- function(N, catch, a, b, k, Linf, sizeM, vectorM, freq, sp, Ts){

  matrixN   <- array(dim = c(Ts + 1, dim(N)))
  matrixB   <- matrix(nrow = Ts + 1, ncol = ncol(N))
  matrixBD  <- matrix(nrow = Ts + 1, ncol = ncol(N))
  matrixBDR <- numeric(ncol(N))

  for(i in seq_len(ncol(N))){

    N0  <- N[, i]
    sim <- .projectPOPEInverse(N = N0, catch = catch, a = a, b = b, k = k, Linf = Linf, sizeM = sizeM,
                               vectorM = vectorM, freq = freq, sp = sp, Ts = Ts)

    matrixN[, ,i] <- sim$N
    matrixB[, i]  <- sim$B
    matrixBD[, i] <- sim$BD
    matrixBDR[i]  <- sim$BDR
  }

  N    <- apply(matrixN, c(1,2), median)
  B    <- apply(matrixB, 1, median)
  BD   <- apply(matrixBD, 1, median)
  BDR  <- median(matrixBDR)

  rawData <- list(N = matrixN, B = matrixB, BD = matrixBD, BDR = matrixBDR)

  output <- list(N = N, B = B, BD = BD, BDR = BDR, raw = rawData)

  attr(output, which = "sp")   <- sp
  attr(output, which = "freq") <- freq
  attr(output, which = "Ts")   <- Ts

  class(output) <- "surveyProj"

  return(output)
}


#' Title
#'
#' @param file
#' @param sp
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
readAtLength = function(file, sp = "anchoveta", ...){
  if(is.null(file)) return(NULL)

  base <- read.csv(file, stringsAsFactors = FALSE, ...)
  colnames(base)[1] <- "x"
  specie = getSpeciesInfo(sp)
  marcas = .createMarks(specie)
  newBase = expand.grid(x = marcas)
  base = merge(base, newBase, all = T)
  base = as.matrix(base[,-1])
  base[is.na(base)] = 0

  return(base)
}
