source("LSM_Single_Trial.R")
library(psych)
library(retimes)
library(gamlss.dist)

LSMSimulation <- function(
  Participants = 40,
  Num_SubBlocks = 128,
  Simulations = 1,
  Length_Sub_Block = 6,
  Trial_Comparisons = c(5, 6),
  LangOrder = c("L1", "L1", "L1", "L2", "L1", "L1"),
  SwitchOrder = c("Stay", "Stay", "Stay", "Switch", "Switch", "Stay"),
  SemantOrder1 = c("False", "True", "True", "True", "True", "True"),
  SemantOrder2 = c("False", rep("True", 3), "False", "True"),
  Comp =.53,
  L1_Activation = 3, #rnorm(1, 4, .5),
  L2_Activation = 1.5, #rnorm(1, 2, .5)
  L1_Distractors = 3, #Most Active distractor
  L2_Distractors = 1.5,
  L1_Initial = L1_Activation,
  L1_Initial_D = L1_Distractors,
  L2_Initial = L2_Activation,
  L2_Initial_D = L2_Distractors,
  L1_Others = 3, #other distractors that get spreading activation as well
  L2_Others = 1.5,
  L1_Initial_Others = L1_Others,
  L2_Initial_Others = L2_Others,
  L1_Strength = 1.5,
  L2_Strength = 5,
  t = 0.0,
  u = .01,
  p1 = .75/2,
  p2 = ((1-p1)/5)/2,
  p3 = ((1-p1)*(4/5))/2,
  ISI = 1000/20,
  c = 0.01,
  h1 = 1.5,
  h2 = 1.5,
  Inh1 = 3.0,
  Inh2 = 3.0,
  NoiseMu = rnorm(1, 720, 50),
  NoiseTau = rnorm(1, 400, 50),
  NoiseSigma = rnorm(1, 100, 5)
){
  listofdfs <- list()
  listofdfs2 <- list()
  listofdfs3 <- list()

  pb0 = txtProgressBar(min = 0, max = Simulations, style = 3)
  for(h in 1:Simulations){
    #print(h)
    pb = txtProgressBar(min = 0, max = Participants, style = 3)

    for(i in 1:Participants){
      Noise1 = Noise1
      Noise2 = Noise2
      Noise3 = Noise3
      x = rep(LangOrder, Num_SubBlocks)
      y = rep(SwitchOrder, Num_SubBlocks)
      z = c(rep(SemantOrder1, Num_SubBlocks/2), rep(SemantOrder2, Num_SubBlocks/2))
      Condition = c(rep("SemantOrder1", length(SemantOrder1)*Num_SubBlocks/2), rep("SemantOrder2", length(SemantOrder1)*Num_SubBlocks/2))
      Nx = length(x)
      Cond = vector("numeric", length = Nx)
      Simulated_RT = vector("numeric", length = Nx)
      L1Act = vector("numeric", length = Nx)
      L2Act = vector("numeric", length = Nx)
      L1_Dist = vector("numeric", length = Nx)
      L2_Dist = vector("numeric", length = Nx)
      L1_Max = vector("numeric", length = Nx)
      L2_Max = vector("numeric", length = Nx)
      L1_Max_Dist = vector("numeric", length = Nx)
      L2_Max_Dist = vector("numeric", length = Nx)
      L1_Init = vector("numeric", length = Nx)
      L2_Init = vector("numeric", length = Nx)
      L1_Init_D = vector("numeric", length = Nx)
      L2_Init_D = vector("numeric", length = Nx)
      L1_O = vector("numeric", length = Nx)
      L1_Init_O = vector("numeric", length = Nx)
      L2_O = vector("numeric", length = Nx)
      L2_Init_O = vector("numeric", length = Nx)
      Trial = rep(1:Length_Sub_Block, Num_SubBlocks)

      for(k in 1:Nx){
        #Noise1 = rexGAUS(1,NoiseMu1,NoiseTau1, NoiseSigma1)
        if(k == 1){
          resp1 = sim_single_trial(Language = x[k], Switching = y[k], Semantic = z[k], Comp = Comp, L1_Activation = L1_Activation, L2_Activation = L2_Activation, L1_Distractors = L1_Distractors, L2_Distractors = L2_Distractors, L1_Initial = L1_Initial, L1_Initial_D = L1_Initial_D, L2_Initial = L2_Initial, L2_Initial_D = L2_Initial_D, L1_Others = L1_Others, L2_Others = L2_Others, L1_Initial_Others = L1_Initial_Others, L2_Initial_Others = L2_Initial_Others, L1_Strength = L1_Strength, L2_Strength = L2_Strength, t = t, u = u, p1 = p1, p2 = p2, p3 = p3, ISI = ISI,c = c, h1 = h1, h2 = h2, Inh1 = Inh1, Inh2 = Inh2, NoiseMu = Noise1, NoiseTau = Noise2, NoiseSigma = Noise3 )
        } else{
          resp1 = sim_single_trial(Language = x[k], Switching = y[k], Semantic = z[k], L1_Activation = resp1$L1_Activation, L2_Activation = resp1$L2_Activation, L1_Distractors = resp1$L1_Distractors, L2_Distractors = resp1$L2_Distractors, L2_Others = resp1$L2_Others, L1_Others = resp1$L1_Others, t = t, u = u, p1 = p1, p2 = p2, p3 = p3, ISI = ISI,c = c, h1 = h1, h2 = h2, Inh1 = Inh1, Inh2 = Inh2, Comp = Comp, NoiseMu = Noise1, NoiseTau = Noise2, NoiseSigma = Noise3)
          }

        Simulated_RT[k] = resp1$RT
        Cond[k] = Condition[k]
        L1Act[k] = resp1$L1_Activation
        L2Act[k] = resp1$L2_Activation
        L1_Dist[k] = resp1$L1_Distractors
        L2_Dist[k] = resp1$L2_Distractors
        L1_Max[k] = resp1$L1_Max_Activation
        L2_Max[k] = resp1$L2_Max_Activation
        L1_Max_Dist[k] = resp1$L1_Dist_Max
        L2_Max_Dist[k] = resp1$L2_Dist_Max
        L1_Init[k] = resp1$L1_Initial
        L2_Init[k] = resp1$L2_Initial
        L1_Init_D[k] = resp1$L1_Initial_D
        L2_Init_D[k] = resp1$L2_Initial_D
        L1_O[k] = resp1$L1_Others
        L2_O[k] = resp1$L2_Others
        L1_Init_O[k] = resp1$L1_Initial_Others
        L2_Init_O[k] = resp1$L2_Initial_Others


      }
      df = data.frame(Trial, Simulated_RT, x, y, z, L1Act, L2Act, L1_Dist, L2_Dist, L1_Max, L2_Max, L1_Max_Dist, L2_Max_Dist, L1_Init, L2_Init, L1_Init_D, L2_Init_D, L1_O, L2_O, L1_Init_O, L2_Init_O, Condition2 = Cond)
      df$Participant <- i
      listofdfs[[i]]<- df #Each Participant gets a df
      setTxtProgressBar(pb,i)
    }
    close(pb)
    All_Blocks2 = do.call("rbind", listofdfs) #Bind all the participants together into a big df

    comparison = vector()
    Sem = vector()
    TrialNo = vector()
    condlist = c("SemantOrder1", "SemantOrder2")
    o = 0

    for(m in 1:length(Trial_Comparisons)){
      for(n in 1:length(condlist)){
        o = o + 1
        dfcond <- subset(All_Blocks2, All_Blocks2$Trial == Trial_Comparisons[m])
        dfcond <- subset(dfcond, dfcond$Condition2 == condlist[n])
        exgausmean <- timefit(dfcond$Simulated_RT)
        comparison[o] <- as.numeric(exgausmean@par[1] + exgausmean@par[2])
        Sem[o] <- paste0(condlist[n])
        TrialNo[o] <- paste0("Trial ", Trial_Comparisons[m])
        }
      }


    df2<- data.frame(Mean.RT = comparison, Condition = Sem, TrialNo)
    df2$Simulation = h

    listofdfs2[[h]] <- df2
    Sys.sleep(.1)
    setTxtProgressBar(pb0,h)
  }
  close(pb)
  close(pb0)
  df3 <- do.call("rbind", listofdfs2)
  resultList <- list(Last_Simulation = All_Blocks2, Simulation.Results = df3)

  }
