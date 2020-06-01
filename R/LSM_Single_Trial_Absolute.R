library(gamlss.dist)
library(retimes)

ifelse(isNamespaceLoaded("gamlss.dist") == FALSE,
  paste0("Load Package gamlss.dist in order to run simulation"),


sim_single_trial_Abs <- function(
  Language = c("L1","L2"),
  Switching = c("Switch","Stay"),
  Semantic = c("True", "False"),
  NoiseMu = 100,
  NoiseTau = 400,
  NoiseSigma = 100,
  Noise = rexGAUS(1,NoiseMu,NoiseTau, NoiseSigma),
  Comp =.53,
  L1_Activation = 3, #rnorm(1, 4, .5),
  L2_Activation = 3.0, #rnorm(1, 2, .5)
  L1_Distractors = 3, #Most Active distractor
  L2_Distractors = 3.0,
  L1_Initial = L1_Activation,
  L1_Initial_D = L1_Distractors,
  L2_Initial = L2_Activation,
  L2_Initial_D = L2_Distractors,
  L1_Others = 3, #other distractors that get spreading activation as well
  L2_Others = 3.0,
  L1_Initial_Others = L1_Others,
  L2_Initial_Others = L2_Others,
  L1_Strength = 1.5,
  L2_Strength = 1.5,
  RT_L1 = Noise,
  RT_L2 = Noise,
  t = 0.0,
  u = 1.0, #changed from .1 to increase effect
  p1 = .75/2,
  p2 = ((1-p1)/5)/2,
  p3 = ((1-p1)*(4/5))/2,
  ISI = 1000/20,
  c = 0.01,
  h1 = 1.5,
  h2 = 1.5,
  Inh1 = 3.0,
  Inh2 = 3.0
){
  if((Semantic == "True") && (Switching == "Stay")){
    L1_Initial = L1_Activation
    L2_Initial = L2_Activation
    L1_Initial_D = L1_Distractors
    L2_Initial_D = L2_Distractors
    L1_Initial_Others = L1_Others
    L2_Initial_Others = L2_Others
    if(Language == "L1"){
      safe = 0
      repeat{
        if(( L1_Activation >= L1_Activation && L1_Activation >= Comp)  || safe > 1000){
          break
        }
        SA2 = 1/(1+L2_Strength*exp(1)^(-t)) #semantic activation
        DK2 = (SA2 - (L2_Strength*exp(1)^(-t))/(1+L2_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L2_Distractors = L2_Initial_D + p3*SA2 #+ p3*DK
        L2_Others = L2_Initial_Others + p2*SA2 #+ p2*DK
        L2_Activation = L2_Initial + p1*SA2 #+ p1*DK

        SA1 = 1/(1+ L1_Strength*exp(1)^(-t)) #semantic activation
        DK1 = (SA1 - (L1_Strength*exp(1)^(-t))/(1+L1_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L1_Distractors = L1_Initial_D + p3*SA1 #+ p3*DK1
        L1_Others = L1_Initial_Others + p2*SA1 #+ p2*DK1
        L1_Activation = L1_Initial + p1*SA1 #+ p1*DK1

        t = t + u
        RT_L1 = RT_L1 + 20*u*t
        safe = safe + 1
      }

      RT = RT_L1

      L1_Max_Activation = L1_Activation
      L1_Dist_Max = L1_Distractors
      L1_Distractors = (L1_Activation - L1_Initial)*exp(1)^(-c*ISI) + L1_Initial
      L1_Activation = (L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others
      L1_Others =(L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others

      L2_Max_Activation = L2_Activation
      L2_Dist_Max = L2_Distractors
      L2_Distractors = (L2_Activation - L2_Initial)*exp(1)^(-c*ISI) + L2_Initial
      L2_Activation = (L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
      L2_Others =(L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others

    }
    if(Language == "L2"){
      safe = 0
      repeat{
        if(( L2_Activation >= L2_Activation && L2_Activation >= Comp)  || safe > 1000){
          break
        }
        SA2 = 1/(1+L2_Strength*exp(1)^(-t)) #semantic activation
        DK2 = (SA2 - (L2_Strength*exp(1)^(-t))/(1+L2_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L2_Distractors = L2_Initial_D + p3*SA2 #+ p3*DK
        L2_Others = L2_Initial_Others + p2*SA2 #+ p2*DK
        L2_Activation = L2_Initial + p1*SA2 #+ p1*DK

        SA1 = 1/(1+ L1_Strength*exp(1)^(-t)) #semantic activation
        DK1 = (SA1 - (L1_Strength*exp(1)^(-t))/(1+L1_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L1_Distractors = L1_Initial_D + p3*SA1 #+ p3*DK1
        L1_Others = L1_Initial_Others + p2*SA1 #+ p2*DK1
        L1_Activation = L1_Initial + p1*SA1 #+ p1*DK1
        t = t + u
        RT_L2 = RT_L2 + 20*u*t
        safe = safe +1
      }

      RT = RT_L2


      L1_Max_Activation = L1_Activation
      L1_Dist_Max = L1_Distractors
      L1_Distractors = (L1_Activation - L1_Initial)*exp(1)^(-c*ISI) + L1_Initial
      L1_Activation = (L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others
      L1_Others =(L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others

      L2_Max_Activation = L2_Activation
      L2_Dist_Max = L2_Distractors
      L2_Distractors = (L2_Activation - L2_Initial)*exp(1)^(-c*ISI) + L2_Initial
      L2_Activation = (L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
      L2_Others =(L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
    }
  }
  if(Semantic == "False" && Switching == "Stay"){
    L1_Initial = 3
    L1_Initial_D = 3
    L1_Activation = 3
    L1_Distractors = 3
    L1_Initial_Others = 3
    L1_Others = 3

    L2_Initial = 3.0
    L2_Initial_D = 3.0
    L2_Initial_Others = 3.0
    if(Language == "L1"){
      safe = 0
      repeat{
        if(( L1_Activation >= L1_Activation && L1_Activation >= Comp) || safe >1000){
          break
        }
        SA2 = 1/(1+L2_Strength*exp(1)^(-t)) #semantic activation
        DK2 = (SA2 - (L2_Strength*exp(1)^(-t))/(1+L2_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L2_Distractors = L2_Initial_D + p3*SA2 #+ p3*DK
        L2_Others = L2_Initial_Others + p2*SA2 #+ p2*DK
        L2_Activation = L2_Initial + p1*SA2 #+ p1*DK

        SA1 = 1/(1+ L1_Strength*exp(1)^(-t)) #semantic activation
        DK1 = (SA1 - (L1_Strength*exp(1)^(-t))/(1+L1_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L1_Distractors = L1_Initial_D + p3*SA1 #+ p3*DK1
        L1_Others = L1_Initial_Others + p2*SA1 #+ p2*DK1
        L1_Activation = L1_Initial + p1*SA1 #+ p1*DK1

        t = t + u
        RT_L1 = RT_L1 + 20*u*t
        safe = safe + 1
      }

      RT = RT_L1


      L1_Max_Activation = L1_Activation
      L1_Dist_Max = L1_Distractors
      L1_Distractors = (L1_Activation - L1_Initial)*exp(1)^(-c*ISI) + L1_Initial
      L1_Activation = (L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others
      L1_Others =(L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others

      L2_Max_Activation = L2_Activation
      L2_Dist_Max = L2_Distractors
      L2_Distractors = (L2_Activation - L2_Initial)*exp(1)^(-c*ISI) + L2_Initial
      L2_Activation = (L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
      L2_Others =(L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
    }
    if(Language == "L2"){
      safe = 0
      repeat{
        if(( L2_Activation >= L2_Activation && L2_Activation) >= Comp  || safe > 1000){
          break
        }
        SA2 = 1/(1+L2_Strength*exp(1)^(-t)) #semantic activation
        DK2 = (SA2 - (L2_Strength*exp(1)^(-t))/(1+L2_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L2_Distractors = L2_Initial_D + p3*SA2 #+ p3*DK
        L2_Others = L2_Initial_Others + p2*SA2 #+ p2*DK
        L2_Activation = L2_Initial + p1*SA2 #+ p1*DK

        SA1 = 1/(1+ L1_Strength*exp(1)^(-t)) #semantic activation
        DK1 = (SA1 - (L1_Strength*exp(1)^(-t))/(1+L1_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L1_Distractors = L1_Initial_D + p3*SA1 #+ p3*DK1
        L1_Others = L1_Initial_Others + p2*SA1 #+ p2*DK1
        L1_Activation = L1_Initial + p1*SA1 #+ p1*DK1
        t = t + u
        RT_L2 = RT_L2 + 20*u*t
        safe = safe + 1
      }

      RT = RT_L2


      L1_Max_Activation = L1_Activation
      L1_Dist_Max = L1_Distractors
      L1_Distractors = (L1_Activation - L1_Initial)*exp(1)^(-c*ISI) + L1_Initial
      L1_Activation = (L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others
      L1_Others =(L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others

      L2_Max_Activation = L2_Activation
      L2_Dist_Max = L2_Distractors
      L2_Distractors = (L2_Activation - L2_Initial)*exp(1)^(-c*ISI) + L2_Initial
      L2_Activation = (L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
      L2_Others =(L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
    }
  }
  if(Semantic == "True" && Switching == "Switch"){
    L1_Initial = L1_Activation
    L2_Initial = L2_Activation
    L1_Initial_Others = L1_Others
    L2_Initial_Others = L2_Others
    L1_Initial_D = L1_Distractors
    L2_Initial_D = L2_Distractors
    if(Language == "L1"){
      safe = 0
      repeat{

        if(( L1_Activation >= L1_Activation && L1_Activation >= Comp)  || safe > 1000){
          break
        }
        SA2 = 1/(1+(Inh2*h2)*L2_Strength*exp(1)^(-t)) #semantic activation
        DK2 = (SA2 - ((Inh2*h2)*L2_Strength*exp(1)^(-t))/(1+L2_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L2_Distractors = L2_Initial_D + p3*SA2 #+ p3*DK
        L2_Others = L2_Initial_Others + p2*SA2 #+ p2*DK
        L2_Activation = L2_Initial + p1*SA2 #+ p1*DK

        SA1 = 1/(1+(Inh1*h1)*L1_Strength*exp(1)^(-t)) #semantic activation
        DK1 = (SA1 - ((Inh1*h1)*L1_Strength*exp(1)^(-t))/(1+L1_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L1_Distractors = L1_Initial_D + p3*SA1 #+ p3*DK1
        L1_Others = L1_Initial_Others + p2*SA1 #+ p2*DK1
        L1_Activation = L1_Initial + p1*SA1 #+ p1*DK1
        t = t + u
        RT_L1 = RT_L1 + 20*u*t
        safe = safe + 1
      }

      RT = RT_L1


      L1_Max_Activation = L1_Activation
      L1_Dist_Max = L1_Distractors
      L1_Distractors = (L1_Activation - L1_Initial)*exp(1)^(-c*ISI) + L1_Initial
      L1_Activation = (L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others
      L1_Others =(L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others

      L2_Max_Activation = L2_Activation
      L2_Dist_Max = L2_Distractors
      L2_Distractors = (L2_Activation - L2_Initial)*exp(1)^(-c*ISI) + L2_Initial
      L2_Activation = (L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
      L2_Others =(L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
    }
    if(Language == "L2"){
    safe = 0
      repeat{
        if(( L2_Activation >= L2_Activation && L2_Activation >= Comp)  || safe > 1000){
          break
        }
        SA2 = 1/(1+(Inh2*h2)*L2_Strength*exp(1)^(-t)) #semantic activation
        DK2 = (SA2 - ((Inh2*h2)*L2_Strength*exp(1)^(-t))/(1+L2_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L2_Distractors = L2_Initial_D + p3*SA2 #+ p3*DK
        L2_Others = L2_Initial_Others + p2*SA2 #+ p2*DK
        L2_Activation = L2_Initial + p1*SA2 #+ p1*DK

        SA1 = 1/(1+(Inh1*h1)*L1_Strength*exp(1)^(-t)) #semantic activation
        DK1 = (SA1 - ((Inh1*h1)*L1_Strength*exp(1)^(-t))/(1+L1_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L1_Distractors = L1_Initial_D + p3*SA1 #+ p3*DK1
        L1_Others = L1_Initial_Others + p2*SA1 #+ p2*DK1
        L1_Activation = L1_Initial + p1*SA1 #+ p1*DK1
        t = t + u
        RT_L2 = RT_L2 + 20*u*t
        safe = safe + 1
      }

      RT = RT_L2


      L1_Max_Activation = L1_Activation
      L1_Dist_Max = L1_Distractors
      L1_Distractors = (L1_Activation - L1_Initial)*exp(1)^(-c*ISI) + L1_Initial
      L1_Activation = (L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others
      L1_Others =(L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others

      L2_Max_Activation = L2_Activation
      L2_Dist_Max = L2_Distractors
      L2_Distractors = (L2_Activation - L2_Initial)*exp(1)^(-c*ISI) + L2_Initial
      L2_Activation = (L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
      L2_Others =(L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others

    }
  }
  if(Semantic == "False" && Switching == "Switch"){
    L1_Initial = 3
    L2_Initial = 3.0
    L1_Initial_D = 3
    L2_Initial_D = 3.0
    L1_Initial_Others = 3
    L2_Initial_Others = 3.0
    if(Language == "L1"){
      safe = 0

      repeat{
        if(( L1_Activation >= L1_Activation && L1_Activation >= Comp)  || safe > 1000 ){
          break
        }
        SA2 = 1/(1+(Inh2*h2)*L2_Strength*exp(1)^(-t)) #semantic activation
        DK2 = (SA2 - ((Inh2*h2)*L2_Strength*exp(1)^(-t))/(1+L2_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L2_Distractors = L2_Initial_D + p3*SA2 #+ p3*DK
        L2_Others = L2_Initial_Others + p2*SA2 #+ p2*DK
        L2_Activation = L2_Initial + p1*SA2 #+ p1*DK

        SA1 = 1/(1+(Inh1*h1)*L1_Strength*exp(1)^(-t)) #semantic activation
        DK1 = (SA1 - ((Inh1*h1)*L1_Strength*exp(1)^(-t))/(1+L1_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L1_Distractors = L1_Initial_D + p3*SA1 #+ p3*DK1
        L1_Others = L1_Initial_Others + p2*SA1 #+ p2*DK1
        L1_Activation = L1_Initial + p1*SA1 #+ p1*DK1

        t = t + u
        RT_L1 = RT_L1 + 20*u*t
        safe = safe +1
      }

      RT = RT_L1


      L1_Max_Activation = L1_Activation
      L1_Dist_Max = L1_Distractors
      L1_Distractors = (L1_Activation - L1_Initial)*exp(1)^(-c*ISI) + L1_Initial
      L1_Activation = (L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others
      L1_Others =(L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others

      L2_Max_Activation = L2_Activation
      L2_Dist_Max = L2_Distractors
      L2_Distractors = (L2_Activation - L2_Initial)*exp(1)^(-c*ISI) + L2_Initial
      L2_Activation = (L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
      L2_Others =(L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
    }
    if(Language == "L2"){
      safe = 0
      repeat{
        if(( L2_Activation >= L2_Activation && L2_Activation >= Comp)  || safe > 1000){
          break
        }
        SA2 = 1/(1+(Inh2*h2)*L2_Strength*exp(1)^(-t)) #semantic activation
        DK2 = (SA2 - ((Inh2*h2)*L2_Strength*exp(1)^(-t))/(1+L2_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L2_Distractors = L2_Initial_D + p3*SA2 #+ p3*DK
        L2_Others = L2_Initial_Others + p2*SA2 #+ p2*DK
        L2_Activation = L2_Initial + p1*SA2 #+ p1*DK

        SA1 = 1/(1+(Inh1*h1)*L1_Strength*exp(1)^(-t)) #semantic activation
        DK1 = (SA1 - ((Inh1*h1)*L1_Strength*exp(1)^(-t))/(1+L1_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)
        L1_Distractors = L1_Initial_D + p3*SA1 #+ p3*DK1
        L1_Others = L1_Initial_Others + p2*SA1 #+ p2*DK1
        L1_Activation = L1_Initial + p1*SA1 #+ p1*DK1
        t = t + u
        RT_L2 = RT_L2 + 20*u*t
        safe = safe + 1
      }
      RT = RT_L2


      L1_Max_Activation = L1_Activation
      L1_Dist_Max = L1_Distractors
      L1_Distractors = (L1_Activation - L1_Initial)*exp(1)^(-c*ISI) + L1_Initial
      L1_Activation = (L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others
      L1_Others =(L1_Others - L1_Initial_Others)*exp(1)^(-c*ISI) + L1_Initial_Others

      L2_Max_Activation = L2_Activation
      L2_Dist_Max = L2_Distractors
      L2_Distractors = (L2_Activation - L2_Initial)*exp(1)^(-c*ISI) + L2_Initial
      L2_Activation = (L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
      L2_Others =(L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
    }
  }
  response = data.frame(RT, L1_Activation, L2_Activation, L1_Distractors, L2_Distractors, L1_Max_Activation, L2_Max_Activation, L1_Dist_Max, L2_Dist_Max, L1_Initial, L2_Initial, L1_Initial_D, L2_Initial_D, L1_Others, L1_Initial_Others, L2_Others, L2_Initial_Others, Comp)
}
)
