
# Proyectar hacia adelante ------------------------------------------------

sp <- "anchoveta"

fileAbundancia <- "../auxiliar/abundancia160304.csv"
fileCaptura    <- "../auxiliar/capturas1.1.csv"

marcas <- seq(2, 20, 0.5)

N1 <- readAtLength(file = fileAbundancia, sp = sp)
N1 <- matrix(rep(N1, each = 2), ncol = 2, byrow = TRUE)

captura <- readAtLength(fileCaptura, sp = sp)
Ts <- 1

sim1 <- projectPOPE(N = N1, captura[,4], a = 0.003925, b = 3.2178, k = 0.83, Linf = 19.21, sizeM = c(0, 8, 12),
                    vectorM = c(0.8, 0.8, 0.8), freq = 12, sp = "anchoveta", Ts = Ts)

abundancia1 <- cbind(sim1$N[1, ], sim1$N[2, ])

par(mfrow = c(1,2))
matplot(marcas, abundancia1, type = "l", axes = FALSE, ylim = c(0, 70000))
axis(1, seq(2,20, 0.5), las = 2)
axis(2)
box()



# Proyectar hacia atrás ---------------------------------------------------

sim1.1 <- projectPOPEInverse(N = sim1$raw$N[Ts+1,,], captura[, 4], a = 0.003925, b = 3.2178, k = 0.83, Linf = 19.21,
                             sizeM = c(0, 8, 12), vectorM = c(0.8, 0.8, 0.8), freq = 12, sp = "anchoveta", Ts = Ts)

abundancia2 <- cbind(sim1.1$N[1, ], sim1.1$N[2, ])

matplot(marcas, abundancia2, type = "l", axes = FALSE, ylim = c(0, 70000))
axis(1, seq(2,20,0.5), las = 2)
axis(2)
box()



# Hacer proyecciones ------------------------------------------------------
allCatches <- "../auxiliar/capturas1.1.csv"

allSurveys <- "../auxiliar/abundancia160304.csv"

growthParams <- "../auxiliar/growthParams.csv"

allProjections <- NULL

catchDataFactor <- 1e-3

surveyNames <- paste("Cr.", c("160304"))

allData <- autoProjector(allCatches = allCatches, allSurveys = allSurveys, surveyNames = surveyNames,
                         growthParams = growthParams, catchDataFactor = 1)
write.csv(x = allData$proyecciones, file = "../auxiliar/allProjections.csv", na = "")


# Hacer gráficos ----------------------------------------------------------

rangeDate <- c("May-16", "Ene-17")

nStairs <- 13

allData <- allData

ylimProj <- c(0, 150000, 50000) # Left
ylimCatch <- c(0, 150000, 50000) # Right
cols <- c("green4", "red", "blue")
absolute <- TRUE

autoProjectorPlot(allData = allData, rangeDate = rangeDate, nStairs = nStairs, outputDir = "../auxiliar/",
                  ylimProj = ylimProj, ylimCatch = ylimCatch, cols = cols, absolute = absolute)
