allCatches <- "../auxiliar/catches_may.csv"
allCatches <- "../auxiliar/allCatches.csv"
allCatches <- "../auxiliar/capturas1.1.csv"

allSurveys <- "../auxiliar/abundance_may.csv"
allSurveys <- "../auxiliar/allSurveys.csv"
allSurveys <- "../auxiliar/abundancia160304.csv"

growthParams <- "../auxiliar/growthParams.csv"

allProjections <- NULL

catchDataFactor <- 1e-3

surveyNames <- paste("Cr.", c("160304"))
# surveyNames <- paste("Cr.", c("110204", "110809", "120204", "120809",
#                               "130204", "130809", "140204", "140810-Aj",
#                               "150204", "150810", "160304", "160506"))
# surveyNames <- paste("Cr.", c("160304", "160506"))



# Hacer proyecciones ------------------------------------------------------

allData <- autoProjector(allCatches = allCatches, allSurveys = allSurveys, surveyNames = surveyNames,
                         growthParams = growthParams, catchDataFactor = 1)
write.csv(x = allData$proyecciones, file = "../auxiliar/allProjections.csv", na = "")


# Hacer gráficos ----------------------------------------------------------

rangeDate <- c("May-16", "Ene-17")

nStairs <- 4

allData <- allData

ylimProj <- c(0, 150000, 50000) # Left
ylimCatch <- c(0, 150000, 50000) # Right
cols <- c("green4", "red", "blue")
absolute <- TRUE

autoProjectorPlot(allData = allData, rangeDate = rangeDate, nStairs = nStairs, outputDir = "../auxiliar/",
                  ylimProj = ylimProj, ylimCatch = ylimCatch, cols = cols, absolute = absolute)
