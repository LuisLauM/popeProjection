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
autoProjector <- function(allCatches, allSurveys, surveyNames, allProjections = NULL,
                          catchDataFactor = 1e-3){

  myFolder <- dirname(allCatches)

  allCatches <- read.csv(file = allCatches, stringsAsFactors = FALSE,
                         check.names = FALSE, row.names = 1, na.strings = 0)
  allSurveys <- read.csv(file = allSurveys, stringsAsFactors = FALSE,
                         check.names = FALSE, row.names = 1, na.strings = 0)

  allCatches <- allCatches*catchDataFactor

  allCatches[is.na(allCatches)] <- 0
  allSurveys[is.na(allSurveys)] <- 0

  if(length(surveyNames) != ncol(allSurveys)){
    stop("Lengths of surveyNames and ncol(allSurveys) must be the same.")
  }

  growthParams <- allSurveys[1:4,]
  allSurveys <-  allSurveys[5:nrow(allSurveys),]

  allMarks <- seq(2, 20, 0.5)

  catchDates <- getDates(colnames(allCatches))
  surveyDates <- getDates(colnames(allSurveys))

  if(is.null(allProjections)){
    cat("\nCreating projection table...\n")

    allProjections <- matrix(data = NA, nrow = length(allMarks) + 1, ncol = ncol(allCatches),
                             dimnames = list(c("fixed", allMarks),
                                             paste0(substr(as.yearmon(catchDates), 1, 3), "-",
                                                    substr(as.yearmon(catchDates), 8, 9))))
    allProjections[1,] <- 0

    write.csv(x = allProjections, file = file.path(myFolder, "allProjections.csv"), na = "")

    allProjections <- file.path(myFolder, "allProjections.csv")

    cat(paste("\n'allProjections' file created in", allProjections, "\n"))
  }

  allProjections <- read.csv(file = allProjections, stringsAsFactors = FALSE,
                             check.names = FALSE, row.names = 1)

  projDates <- getDates(colnames(allProjections))

  projCondition <- an(allProjections[1,])
  allProjections <- allProjections[-1,]

  diffIndex <- an(round(diff(surveyDates)/30, 0))
  monthDiffs <- c(0, cumsum(diffIndex))
  allGrowth <- matrix(data = NA, nrow = 4, ncol = ncol(allProjections),
                      dimnames = list(rownames(growthParams), colnames(allProjections)))
  for(i in seq(ncol(allSurveys) - 1)){
    index <- seq(monthDiffs[i] + 1, monthDiffs[i + 1])

    allGrowth[,index] <- matrix(rep(growthParams[,i], times = diffIndex[i]), nrow = 4)
  }

  for(i in seq(ncol(allProjections) - 1)){

    if(is.element(projDates[i], surveyDates)){
      seedLengthVector <- allSurveys[,match(projDates[i], surveyDates)]
    }else{
      seedLengthVector <- allProjections[,i]
    }

    species["anchoveta", "a"] <- allGrowth["a", i]
    species["anchoveta", "b"] <- allGrowth["b", i]
    LHT$anchoveta$G["neutro"] <- c(allGrowth["k", i], allGrowth["Linf", i], -0.210)
    LHT$anchoveta$M["neutro"] <- rep(0.8, 3)

    tempLengths <- cbind(seedLengthVector, seedLengthVector)

    simulation <- projectPOPE(N = tempLengths, catch = allCatches[,i], scenario = "neutro",
                                    freq = 12, sp = "anchoveta", Ts = 1)

    if(!is.na(projCondition[i + 1]) && projCondition[i + 1] == 0){
      allProjections[,i + 1] <- simulation$raw$N[2,,1]
    }

  }

  return(list(capturas = allCatches,
              cruceros = allSurveys,
              proyecciones = allProjections,
              nombres_cruceros = surveyNames))
}

#' Title
#'
#' @param N
#' @param catch
#' @param scenario
#' @param freq
#' @param sp
#' @param Ts
#'
#' @return
#' @export
#'
#' @examples
projectPOPE <- function(N, catch, scenario, freq, sp, Ts){

  matrixN    <- array(dim=c(Ts+1, dim(N)))
  matrixB    <- matrix(nrow=Ts+1, ncol=ncol(N))
  matrixBD   <- matrix(nrow=Ts+1, ncol=ncol(N))
  matrixBDR  <- numeric(ncol(N))

  scenario <- .getScenarios(scenario, ncol(N))

  for(i in seq_len(ncol(N))){

    N0 <- N[, i]
    sim <- .projectPOPE(N=N0, catch=catch, scenario=scenario[i], freq=freq, sp=sp, Ts=Ts)

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
