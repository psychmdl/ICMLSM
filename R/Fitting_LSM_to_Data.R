library(psych)
library(retimes)
library(gamlss.dist)
source('LSM_Simulations_for_Dissertation.R')


LSM_Fit <- function(
  Observed_Data1 = c(0, -38, -54, 94, 83, 39),
  Observed_Data2 = c(0, -39, -34, 75, 82, 5),
  Num_Sims = 10000
){
  Simulated_Data <- list()
  FitData <- list()
  hrecord <- vector()
  cor1 <- vector()
  cor2 <- vector()
  chi1 <- vector()
  pval1 <- vector()
  chi2 <- vector()
  pval2 <- vector()
  L1Record <- vector()
  comprecord <- vector()
  RMS1 <- vector()
  RMS2 <- vector()
  for(i in 1:Num_Sims){
    competition = runif(1, 0, .54)
    L1Strength = runif(1, 1, 10)
    h = runif(1, .01, 10) #Language Activation parameter on switch trials = h*3*Inh2
    x = LSMSimulation(Simulations = 1, Participants = 1, Num_SubBlocks = 8, Trial_Comparisons = c(1,2,3,4,5,6), Comp = competition, h1 = h, h2 = h, L1_Strength = L1Strength, L2_Strength = L1Strength, NoiseMu = .0001, NoiseTau = .0001, NoiseSigma = .001)
    Simulated_Data1 <- x$Simulation.Results[seq(1, 11, 2),]$Mean.RT
    Simulated_Data2 <- x$Simulation.Results[seq(2, 12, 2),]$Mean.RT
    Simulated_Data3 <- data.frame(Simulated_Data1, Simulated_Data2, Simulation = i)
    Simulated_Data[[i]] <- Simulated_Data3

    chidata1 <- data.frame(Simulated_Data1 +500 , Observed_Data1 +500) #Make Sure Chi Square isn't negative
    chidata2 <- data.frame(Simulated_Data2 +500 , Observed_Data2 + 500)

    cor1[i] <- cor(Simulated_Data1, Observed_Data1)
    cor2[i] <- cor(Simulated_Data2, Observed_Data2)
    hrecord[i] = h
    L1Record[i] = L1Strength
    comprecord[i] = competition
    chi1[i] <- suppressWarnings(as.numeric(chisq.test(chidata1)[1]))
    pval1[i] <- suppressWarnings(as.numeric(chisq.test(chidata1)[3]))

    chi2[i] <- suppressWarnings(as.numeric(chisq.test(chidata2)[1]))
    pval2[i] <- suppressWarnings(as.numeric(chisq.test(chidata2)[3]))

    RMS1[i] <- sqrt(mean((Simulated_Data1 - Observed_Data1)^2))
    RMS2[i] <- sqrt(mean((Simulated_Data2 - Observed_Data2)^2))

  }
  df1 <- data.frame(cor1, cor2, hrecord, L1Record, comprecord, chi1, pval1, chi2, pval2, RMS1, RMS2)
  df2 <- do.call("rbind", Simulated_Data)

  resultsList <- list(fits = df1, SimData = df2)

}
