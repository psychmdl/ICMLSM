library(psych)
library(retimes)
library(gamlss.dist)
source('ICM_Simulations_for_Dissertation.R')


ICM_Fit <- function(
  Observed_Data1 = c(0, -38, -54, 94, 83, 39),
  Observed_Data2 = c(0, -39, -34, 75, 82, 5),
  Num_Sims = 10000
){
  Simulated_Data <- list()
  FitData <- list()
  for(i in 1:Num_Sims){
    competition = runif(1, 0, .75)
    L1Strength = runif(1, 1, 10)
    L1Record <- vector()
    h = runif(1, 1, 10) #Language Activation parameter on switch trials = h*3
    hrecord <- vector()
    cor1 <- vector()
    cor2 <- vector()
    x = ICMSimulation(Simulations = 1, Participants = 2, Num_SubBlocks = 8, Trial_Comparisons = c(1,2,3,4,5,6), Comp = competition, h1 = h, h2 = h, L1_Strength = L1Strength, L2_Strength = L1Strength)
    Simulated_Data1 <- x$Simulation.Results[seq(2, 12, 2),]$Mean.RT
    Simulated_Data2 <- x$Simulation.Results[seq(1, 11, 2),]$Mean.RT
    Simulated_Data[i] <- data.frame(Simulated_Data1, Simulated_Data2)

    cor1[i] <- cor(Simulated_Data1[i], Observed_Data1)
    cor2[i] <- cor(Simulated_Data2[i], Observed_Data2)
    hrecord[i] = h
    L1Record[i] = L1Strength
  }
  df1 <- data.frame(cor1, cor2, hrecord, L1Record)

  df2 <- do.call("rbind", Simulated_Data)

  resultsList <- list(fits = df1, SimData = df2)

}
