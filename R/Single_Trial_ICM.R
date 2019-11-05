library(gamlss.dist)
library(retimes)

sim_single_trial_ICM <- function(
  Language = c("L1","L2"),
  Switching = c("Switch","Stay"),
  Semantic = c("True", "False"),
  NoiseMu = 100,
  NoiseTau = 400,
  NoiseSigma = 100,
  Noise = rexGAUS(1,NoiseMu,NoiseTau, NoiseSigma),
  Comp = 0.60,
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
  RT_L1 = Noise,
  RT_L2 = Noise,
  t = 0.0,
  u = 1, #changed from .01 to increase effect
  p1 = .75,
  p2 = ((1-p1)/5),
  p3 = ((1-p1)*(4/5)),
  ISI = 1000/20,
  c = .01,
  h1 = 1, #inhibition parameters
  h2 = 1,
  Y1 = 5*(h1 + L1_Strength), #Reactivation Parameters
  Y2 = 5*(h2 + L2_Strength)
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
        if((Comp <= L1_Activation/sum(L1_Others + L1_Distractors+L2_Activation+L2_Others+L2_Distractors)) | safe > 1000){
          break
        }
        SA = 1/(1+L1_Strength*exp(1)^(-t)) #semantic activation
        DK = (SA - (L1_Strength*exp(1)^(-t))/(1+L1_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)

        L1_Distractors = L1_Initial_D + p3*SA + p3*DK
        L1_Others = L1_Initial_Others + p2*SA + p2*DK
        L1_Activation = L1_Initial + p1*SA + p1*DK

        L2_Distractors = L2_Initial_D*exp(1)^(-h2*t)
        L2_Others =  L2_Initial_Others*exp(1)^(-h2*t)
        L2_Activation =  L2_Initial*exp(1)^(-h2*t)

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
      L2_Distractors = L2_Max_Activation*exp(1)^(-c*ISI)
      L2_Dist_Max = L2_Distractors
      L2_Others = L2_Others*exp(1)^(-c*ISI)
      L2_Activation = L2_Others*exp(1)^(-c*ISI)

    }
    if(Language == "L2"){
      safe = 0
      repeat{
        if((Comp <= L2_Activation/sum(L2_Distractors + L2_Others+L1_Activation+L1_Distractors+L1_Others)) | safe > 1000){
          break
        }
        SA = 1/(1+L2_Strength*exp(1)^(-t)) #semantic activation
        DK = (SA - (L2_Strength*exp(1)^(-t))/(1+L2_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)

        L2_Distractors = L2_Initial_D + p3*SA + p3*DK
        L2_Others = L2_Initial_Others + p2*SA + p2*DK
        L2_Activation = L2_Initial + p1*SA + p1*DK

        L1_Distractors = L1_Initial_D*exp(1)^(-h1*t)
        L1_Others =  L1_Initial_Others*exp(1)^(-h1*t)
        L1_Activation =  L1_Initial*exp(1)^(-h1*t)
        t = t + u
        RT_L2 = RT_L2 + 20*u*t
        safe = safe + 1
      }

      RT = RT_L2


      L2_Max_Activation = L2_Activation
      L2_Dist_Max = L2_Distractors
      L2_Distractors = (L2_Activation - L2_Initial)*exp(1)^(-c*ISI) + L2_Initial
      L2_Activation = (L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
      L2_Others =(L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others

      L1_Max_Activation = L1_Activation
      L1_Distractors = L1_Max_Activation*exp(1)^(-c*ISI)
      L1_Dist_Max = L1_Distractors
      L1_Others = L1_Others*exp(1)^(-c*ISI)
      L1_Activation = L1_Others*exp(1)^(-c*ISI)
    }
  }
  if(Semantic == "False" && Switching == "Stay"){
    if(Language == "L1"){
      L1_Initial = 3
      L1_Initial_D = 3
      L1_Activation = 3
      L1_Distractors = 3
      L1_Initial_Others = 3
      L1_Others = 3

      L2_Initial = L2_Initial
      L2_Initial_D = L2_Initial_D
      L2_Initial_Others = L2_Others
      safe = 0
      repeat{
        if((Comp <= L1_Activation/sum(L1_Others + L1_Distractors+L2_Activation+L2_Others+L2_Distractors)) | safe > 1000){
          break
        }
        SA = 1/(1+L1_Strength*exp(1)^(-t)) #semantic activation
        DK = (SA - (L1_Strength*exp(1)^(-t))/(1+L1_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)

        L1_Distractors = L1_Initial_D + p3*SA + p3*DK
        L1_Others = L1_Initial_Others + p2*SA + p2*DK
        L1_Activation = L1_Initial + p1*SA + p1*DK

        L2_Distractors = L2_Initial_D*exp(1)^(-h2*t)
        L2_Others =  L2_Initial_Others*exp(1)^(-h2*t)
        L2_Activation =  L2_Initial*exp(1)^(-h2*t)

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
      L2_Distractors = L2_Max_Activation*exp(1)^(-c*ISI)
      L2_Dist_Max = L2_Distractors
      L2_Others = L2_Others*exp(1)^(-c*ISI)
      L2_Activation = L2_Others*exp(1)^(-c*ISI)
    }
    if(Language == "L2"){
      safe = 0
      L2_Activation = 1.5
      L2_Distractors = 1.5
      L2_Initial = 1.5
      L2_Initial_D = 1.5
      L2_Others = 1.5
      L2_Initial_Others = 1.5

      L1_Initial = L1_Initial
      L1_Initial_D = L1_Initial_D
      L1_Initial_Others = L1_Initial_Others
      repeat{
        if((Comp <= L2_Activation/sum(L2_Distractors + L2_Others+L1_Activation+L1_Distractors+L1_Others)) | safe > 1000){
          break
        }
        SA = 1/(1+L2_Strength*exp(1)^(-t)) #semantic activation
        DK = (SA - (L2_Strength*exp(1)^(-t))/(1+L2_Strength*exp(1)^(-t))^2)*exp(1)^(c*u)

        L2_Distractors = L2_Initial_D + p3*SA + p3*DK
        L2_Others = L2_Initial_Others + p2*SA + p2*DK
        L2_Activation = L2_Initial + p1*SA + p1*DK

        L1_Distractors = L1_Initial_D*exp(1)^(-h1*t)
        L1_Others =  L1_Initial_Others*exp(1)^(-h1*t)
        L1_Activation =  L1_Initial*exp(1)^(-h1*t)
        t = t + u
        RT_L2 = RT_L2 + 20*u*t
        safe = safe + 1
      }

      RT = RT_L2


      L2_Max_Activation = L2_Activation
      L2_Dist_Max = L2_Distractors
      L2_Distractors = (L2_Activation - L2_Initial)*exp(1)^(-c*ISI) + L2_Initial
      L2_Activation = (L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
      L2_Others =(L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others

      L1_Max_Activation = L1_Activation
      L1_Distractors = L1_Max_Activation*exp(1)^(-c*ISI)
      L1_Dist_Max = L1_Distractors
      L1_Others = L1_Others*exp(1)^(-c*ISI)
      L1_Activation = L1_Others*exp(1)^(-c*ISI)
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
        if(L1_Initial >=3 && L1_Initial_D >= 3 && L1_Initial_Others >= 3 | safe > 1000){
          break
        }
        RA = 1/(1+(Y1)*exp(1)^(-t)) #RA = Reactivation
        L1_Initial = L1_Initial + RA
        L1_Initial_D  = L1_Initial_D + RA
        L1_Initial_Others = L1_Initial_Others + RA

        L2_Distractors = L2_Initial_D*exp(1)^(-h2*t)
        L2_Others =  L2_Initial_Others*exp(1)^(-h2*t)
        L2_Activation =  L2_Initial*exp(1)^(-h2*t)
        t = t + u
        safe = safe + 1
      }
      t1= t-t
      safe = 0

      repeat{
        if((Comp <= L1_Activation/sum(L1_Others + L1_Distractors+L2_Activation+L2_Others+L2_Distractors)) | safe > 1000){
          break
        }
        SA = 1/(1+(Y1)*exp(1)^(-t1)) #semantic activation
        DK = (SA - ((Y1)*exp(1)^(-t1))/(1+(Y1)*exp(1)^(-t1))^2)*exp(1)^(c*u)
        L1_Distractors = L1_Initial_D + p3*SA + p3*DK
        L1_Others = L1_Initial_Others + p2*SA + p2*DK
        L1_Activation = L1_Initial + p1*SA + p1*DK

        L2_Distractors = L2_Initial_D*exp(1)^(-h2*t)
        L2_Others =  L2_Initial_Others*exp(1)^(-h2*t)
        L2_Activation =  L2_Initial*exp(1)^(-h2*t)
        t = t + u
        t1 = t1+u
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
      L2_Distractors = L2_Max_Activation*exp(1)^(-c*ISI)
      L2_Dist_Max = L2_Distractors
      L2_Others = L2_Others*exp(1)^(-c*ISI)
      L2_Activation = L2_Others*exp(1)^(-c*ISI)
    }
    if(Language == "L2"){
      safe = 0
      repeat{
        if(L2_Initial >=1.5 && L2_Initial_D >= 1.5 && L2_Initial_Others >= 1.5 | safe > 1000){
          break
        }
        RA = 1/(1+(Y2)*exp(1)^(-t)) #RA = Reactivation
        L2_Initial = L2_Initial + RA
        L2_Initial_D  = L2_Initial_D + RA
        L2_Initial_Others = L2_Initial_Others + RA

        L1_Distractors = L1_Initial_D*exp(1)^(-h1*t)
        L1_Others =  L1_Initial_Others*exp(1)^(-h1*t)
        L1_Activation =  L1_Initial*exp(1)^(-h1*t)
        t = t + u
        safe = safe + 1
      }
      t1= t-t
      safe = 0
      repeat{
        if((Comp <= L2_Activation/sum(L2_Distractors + L2_Others+L1_Activation+L1_Distractors+L1_Others)) | safe > 1000){
          break
        }
        SA = 1/(1+(Y2)*exp(1)^(-t1)) #semantic activation
        DK = (SA - ((Y2)*exp(1)^(-t1))/(1+Y2*exp(1)^(-t1))^2)*exp(1)^(c*u)

        L2_Distractors = L2_Initial_D + p3*SA + p3*DK
        L2_Others = L2_Initial_Others + p2*SA + p2*DK
        L2_Activation = L2_Initial + p1*SA + p1*DK

        L1_Distractors = L1_Initial_D*exp(1)^(-h1*t)
        L1_Others =  L1_Initial_Others*exp(1)^(-h1*t)
        L1_Activation =  L1_Initial*exp(1)^(-h1*t)
        t = t + u
        t1 = t1+u
        RT_L2 = RT_L2 + 20*u*t
        safe = safe + 1
      }

      RT = RT_L2


      L2_Max_Activation = L2_Activation
      L2_Dist_Max = L2_Distractors
      L2_Distractors = (L2_Activation - L2_Initial)*exp(1)^(-c*ISI) + L2_Initial
      L2_Activation = (L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
      L2_Others =(L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others

      L1_Max_Activation = L1_Activation
      L1_Distractors = L1_Max_Activation*exp(1)^(-c*ISI)
      L1_Dist_Max = L1_Distractors
      L1_Others = L1_Others*exp(1)^(-c*ISI)
      L1_Activation = L1_Others*exp(1)^(-c*ISI)

    }
  }
  if(Semantic == "False" && Switching == "Switch"){

    if(Language == "L1"){
      safe = 0
      L2_Activation = 1.5
      L2_Distractors = 1.5
      L2_Initial = 1.5
      L2_Initial_D = 1.5
      L2_Others = 1.5
      L2_Initial_Others = 1.5

      L1_Initial = L1_Initial
      L1_Initial_D = L1_Initial_D
      L1_Initial_Others = L1_Initial_Others

      repeat{
        if(L1_Initial >=3 && L1_Initial_D >= 3 && L1_Initial_Others >= 3 | safe > 1000){
          break
        }
        RA = 1/(1+(Y1)*exp(1)^(-t)) #RA = Reactivation
        L1_Initial = L1_Initial + RA
        L1_Initial_D  = L1_Initial_D + RA
        L1_Initial_Others = L1_Initial_Others + RA

        L2_Distractors = L2_Initial_D*exp(1)^(-h2*t)
        L2_Others =  L2_Initial_Others*exp(1)^(-h2*t)
        L2_Activation =  L2_Initial*exp(1)^(-h2*t)
        t = t + u
        safe = safe + 1
      }
      t1 = t - t
      safe = 0
      repeat{
        if((Comp <= L1_Activation/sum(L1_Others + L1_Distractors+L2_Activation+L2_Others+L2_Distractors)) | safe > 1000){
          break
        }
        SA = 1/(1+(Y1)*exp(1)^(-t1)) #semantic activation
        DK = (SA - ((Y1)*exp(1)^(-t1))/(1+(Y1)*exp(1)^(-t1))^2)*exp(1)^(c*u)
        L1_Distractors = L1_Initial_D + p3*SA + p3*DK
        L1_Others = L1_Initial_Others + p2*SA + p2*DK
        L1_Activation = L1_Initial + p1*SA + p1*DK

        L2_Distractors = L2_Initial_D*exp(1)^(-h2*t)
        L2_Others =  L2_Initial_Others*exp(1)^(-h2*t)
        L2_Activation =  L2_Initial*exp(1)^(-h2*t)
        t1 = t1 +u
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
      L2_Distractors = L2_Max_Activation*exp(1)^(-c*ISI)
      L2_Dist_Max = L2_Distractors
      L2_Others = L2_Others*exp(1)^(-c*ISI)
      L2_Activation = L2_Others*exp(1)^(-c*ISI)
    }
    if(Language == "L2"){
      safe = 0
      L1_Initial = 3
      L1_Initial_D = 3
      L1_Activation = 3
      L1_Distractors = 3
      L1_Initial_Others = 3
      L1_Others = 3

      L2_Initial = L2_Initial
      L2_Initial_D = L2_Initial_D
      L2_Initial_Others = L2_Others

      repeat{
        if(L2_Initial >=1.5 && L2_Initial_D >= 1.5 && L2_Initial_Others >= 1.5 | safe > 1000){
          break
        }
        RA = 1/(1+(Y2)*exp(1)^(-t)) #RA = Reactivation
        L2_Initial = L2_Initial + RA
        L2_Initial_D  = L2_Initial_D + RA
        L2_Initial_Others = L2_Initial_Others + RA

        L1_Distractors = L1_Initial_D*exp(1)^(-h1*t)
        L1_Others =  L1_Initial_Others*exp(1)^(-h1*t)
        L1_Activation =  L1_Initial*exp(1)^(-h1*t)
        t = t + u
        safe = safe + 1
      }
      t1= t-t
      safe = 0

      repeat{
        if((Comp <= L2_Activation/sum(L2_Distractors + L2_Others+L1_Activation+L1_Distractors+L1_Others)) | safe > 1000){
          break
        }
        SA = 1/(1+(Y2)*exp(1)^(-t1)) #semantic activation
        DK = (SA - ((Y2)*exp(1)^(-t1))/(1+Y2*exp(1)^(-t1))^2)*exp(1)^(c*u)

        L2_Distractors = L2_Initial_D + p3*SA + p3*DK
        L2_Others = L2_Initial_Others + p2*SA + p2*DK
        L2_Activation = L2_Initial + p1*SA + p1*DK

        L1_Distractors = L1_Initial_D*exp(1)^(-h1*t)
        L1_Others =  L1_Initial_Others*exp(1)^(-h1*t)
        L1_Activation =  L1_Initial*exp(1)^(-h1*t)
        t = t + u
        t1 = t1+u
        RT_L2 = RT_L2 + 20*u*t
        safe = safe + 1
      }
      RT = RT_L2


      L2_Max_Activation = L2_Activation
      L2_Dist_Max = L2_Distractors
      L2_Distractors = (L2_Activation - L2_Initial)*exp(1)^(-c*ISI) + L2_Initial
      L2_Activation = (L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others
      L2_Others =(L2_Others - L2_Initial_Others)*exp(1)^(-c*ISI) + L2_Initial_Others

      L1_Max_Activation = L1_Activation
      L1_Distractors = L1_Max_Activation*exp(1)^(-c*ISI)
      L1_Dist_Max = L1_Distractors
      L1_Others = L1_Others*exp(1)^(-c*ISI)
      L1_Activation = L1_Others*exp(1)^(-c*ISI)
    }
  }
  response = data.frame(RT, L1_Activation, L2_Activation, L1_Distractors, L2_Distractors, L1_Max_Activation, L2_Max_Activation, L1_Dist_Max, L2_Dist_Max, L1_Initial, L2_Initial, L1_Initial_D, L2_Initial_D, L1_Others, L1_Initial_Others, L2_Others, L2_Initial_Others, Comp)
}
