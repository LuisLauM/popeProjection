exampleProjection <- PopeProjection(fileAbundancia, fileCaptura,
                                    scenario = "neutro", a = 0.003925, b = 3.2178,
                                    k = 0.83, Linf = 19.21, vectorM = c(0.8, 0.8, 0.8),
                                    freq = 12, sp = "anchoveta", Ts = 1)

exampleProjection$Biomasa
exampleProjection$Abundancia
matplot(marcas, t(exampleProjection$Abundancia), type = "l")
