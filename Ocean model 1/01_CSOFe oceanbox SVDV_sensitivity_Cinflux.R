################################################################################
# 3 box model of ocean iron and sulphur chemistry
# Sebastiaan van de Velde & Filip Meysman
# sensitivity sulfate influx
################################################################################

# Loading required packages

require(AquaEnv)
require(rootSolve)
source(file="CSOFe oceanbox SVDV v05_model only.R")

##################################################################################
# O2 20% PAL
##################################################################################

#==============================================================================
# varying carbon flux
# O2 20% PAL
#==============================================================================

load("00_CSOFe_runsteady.Rdata")

#------------------
# run tests
#------------------

PL        <- sens00$PL
PP_ref    <- 4*1E+18
PL$k_PyFe <- 0
PL$O2_eq  <- sens00$sens[2401]

PP.frac <- c(seq(1,0.95,-0.01),seq(0.94999,0.75,-0.00001),seq(0.74,0.01,-0.01))
#PP.frac <- c(seq(1,0.1,-0.05))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens00$output$C[2401,]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                   <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance           <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                          "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                          "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01a_baseline_20percO2",type="png")

#------------------
# save output
#------------------

sens01a <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01a,file="01a_baseline_20percO2.Rdata")

#==============================================================================
# varying carbon flux
# O2 20% PAL
#==============================================================================

load("01a_baseline_20percO2.Rdata")

#------------------
# run tests
#------------------

PL <- sens01a$PL
PP_ref   <- 4*1E+18

PP.frac <- c(seq(0.01,0.74,0.01),seq(0.75,0.94999,0.00001),seq(0.95,1,0.01))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens01a$output$C[nrow(sens01a$output$C),]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                  <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance          <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                         "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                         "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01b_baseline_20percO2",type="png")

#------------------
# save output
#------------------

sens01b <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01b,file="01b_baseline_20percO2.Rdata")

##################################################################################
# O2 15% PAL
##################################################################################

#==============================================================================
# varying carbon flux
# O2 15% PAL
#==============================================================================

load("00_CSOFe_runsteady.Rdata")

#------------------
# run tests
#------------------

PL        <- sens00$PL
PP_ref    <- 4*1E+18
PL$k_PyFe <- 0
PL$O2_eq  <- sens00$sens[2551]

PP.frac <- c(seq(1,0.8,-0.1),0.75,seq(0.74999,0.60,-0.00001),seq(0.59,0.01,-0.01))
#PP.frac <- c(seq(1,0.1,-0.05))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens00$output$C[256,]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                   <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance           <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                          "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                          "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01a_baseline_15percO2",type="png")

#------------------
# save output
#------------------

sens01a <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01a,file="01a_baseline_15percO2.Rdata")

#==============================================================================
# varying carbon flux
# O2 15% PAL
#==============================================================================

load("01a_baseline_15percO2.Rdata")

#------------------
# run tests
#------------------

PL <- sens01a$PL
PP_ref   <- 4*1E+18

PP.frac <- c(seq(0.01,0.59,0.01),seq(0.60,0.74999,0.00001),0.75,seq(0.8,1,0.1))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens01a$output$C[nrow(sens01a$output$C),]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                  <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance          <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                         "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                         "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01b_baseline_15percO2",type="png")

#------------------
# save output
#------------------

sens01b <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01b,file="01b_baseline_15percO2.Rdata")

##################################################################################
# O2 10% PAL
##################################################################################

#==============================================================================
# varying carbon flux
# O2 10% PAL
#==============================================================================

load("00_CSOFe_runsteady.Rdata")

#------------------
# run tests
#------------------

PL        <- sens00$PL
PP_ref    <- 4*1E+18
PL$k_PyFe <- 0
PL$O2_eq  <- sens00$sens[2701]

PP.frac <- c(seq(1,0.5,-0.1),seq(0.49999,0.4,-0.00001),seq(0.39,0.01,-0.01))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens00$output$C[2701,]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                   <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance           <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                          "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                          "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01a_baseline_10percO2",type="png")

#------------------
# save output
#------------------

sens01a <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01a,file="01a_baseline_10percO2.Rdata")

#==============================================================================
# varying carbon flux
# O2 10% PAL
#==============================================================================

load("01a_baseline_10percO2.Rdata")

#------------------
# run tests
#------------------

PL <- sens01a$PL
PP_ref   <- 4*1E+18

PP.frac <- c(seq(0.01,0.39,0.01),seq(0.4,0.49999,0.00001),seq(0.5,1,0.1))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens01a$output$C[nrow(sens01a$output$C),]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                  <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance          <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                          "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                          "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01b_baseline_10percO2",type="png")

#------------------
# save output
#------------------

sens01b <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01b,file="01b_baseline_10percO2.Rdata")

##################################################################################
# O2 9% PAL
##################################################################################

#==============================================================================
# varying carbon flux
# O2 9% PAL
#==============================================================================

load("00_CSOFe_runsteady.Rdata")

#------------------
# run tests
#------------------

PL        <- sens00$PL
PP_ref    <- 4*1E+18
PL$k_PyFe <- 0
PL$O2_eq  <- sens00$sens[2731]

PP.frac <- c(seq(1,0.5,-0.1),0.45,seq(0.44999,0.35,-0.00001),seq(0.34,0.01,-0.01))
#PP.frac <- c(seq(1,0.1,-0.05))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens00$output$C[2731,]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                   <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance           <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                          "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                          "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01a_baseline_9percO2",type="png")

#------------------
# save output
#------------------

sens01a <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01a,file="01a_baseline_9percO2.Rdata")

#==============================================================================
# varying carbon flux
# O2 9% PAL
#==============================================================================

load("01a_baseline_9percO2.Rdata")

#------------------
# run tests
#------------------

PL <- sens01a$PL
PP_ref   <- 4*1E+18

PP.frac <- c(seq(0.01,0.34,0.01),seq(0.35,0.449999,0.00001),0.45,seq(0.5,1,0.1))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens01a$output$C[nrow(sens01a$output$C),]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                  <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance          <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                         "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                         "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01b_baseline_9percO2",type="png")

#------------------
# save output
#------------------

sens01b <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01b,file="01b_baseline_9percO2.Rdata")

##################################################################################
# O2 8% PAL
##################################################################################

#==============================================================================
# varying carbon flux
# O2 8% PAL
#==============================================================================

load("00_CSOFe_runsteady.Rdata")

#------------------
# run tests
#------------------

PL        <- sens00$PL
PP_ref    <- 4*1E+18
PL$k_PyFe <- 0
PL$O2_eq  <- sens00$sens[2761]

PP.frac <- c(seq(1,0.4,-0.1),seq(0.39999,0.3,-0.00001),seq(0.29,0.01,-0.01))
#PP.frac <- c(seq(1,0.1,-0.05))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens00$output$C[2761,]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                   <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance           <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                          "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                          "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01a_baseline_8percO2",type="png")

#------------------
# save output
#------------------

sens01a <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01a,file="01a_baseline_8percO2.Rdata")

#==============================================================================
# varying carbon flux
# O2 8% PAL
#==============================================================================

load("01a_baseline_8percO2.Rdata")

#------------------
# run tests
#------------------

PL <- sens01a$PL
PP_ref   <- 4*1E+18

PP.frac <- c(seq(0.01,0.29,0.01),seq(0.3,0.39999,0.00001),seq(0.4,1,0.1))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens01a$output$C[nrow(sens01a$output$C),]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                  <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance          <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                         "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                         "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01b_baseline_8percO2",type="png")

#------------------
# save output
#------------------

sens01b <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01b,file="01b_baseline_8percO2.Rdata")

##################################################################################
# O2 7% PAL
##################################################################################

#==============================================================================
# varying carbon flux
# O2 7% PAL
#==============================================================================

load("00_CSOFe_runsteady.Rdata")

#------------------
# run tests
#------------------

PL        <- sens00$PL
PP_ref    <- 4*1E+18
PL$k_PyFe <- 0
PL$O2_eq  <- sens00$sens[2791]

PP.frac <- c(seq(1,0.4,-0.1),0.35,seq(0.34999,0.25,-0.00001),seq(0.24,0.01,-0.01))
#PP.frac <- c(seq(1,0.1,-0.05))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens00$output$C[2791,]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                   <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance           <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                          "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                          "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01a_baseline_7percO2",type="png")

#------------------
# save output
#------------------

sens01a <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01a,file="01a_baseline_7percO2.Rdata")

#==============================================================================
# varying carbon flux
# O2 7% PAL
#==============================================================================

load("01a_baseline_7percO2.Rdata")

#------------------
# run tests
#------------------

PL <- sens01a$PL
PP_ref   <- 4*1E+18

PP.frac <- c(seq(0.01,0.24,0.01),seq(0.25,0.34999,0.00001),0.35,seq(0.4,1,0.1))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens01a$output$C[nrow(sens01a$output$C),]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                  <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance          <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                         "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                         "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01b_baseline_7percO2",type="png")

#------------------
# save output
#------------------

sens01b <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01b,file="01b_baseline_7percO2.Rdata")

##################################################################################
# O2 6% PAL
##################################################################################

#==============================================================================
# varying carbon flux
# O2 6% PAL
#==============================================================================

load("00_CSOFe_runsteady.Rdata")

#------------------
# run tests
#------------------

PL        <- sens00$PL
PP_ref    <- 4*1E+18
PL$k_PyFe <- 0
PL$O2_eq  <- sens00$sens[2821]

PP.frac <- c(seq(1,0.3,-0.1),seq(0.29999,0.2,-0.00001),seq(0.19,0.01,-0.01))
#PP.frac <- c(seq(1,0.1,-0.05))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens00$output$C[2821,]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                   <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance           <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                          "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                          "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01a_baseline_6percO2",type="png")

#------------------
# save output
#------------------

sens01a <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01a,file="01a_baseline_6percO2.Rdata")

#==============================================================================
# varying carbon flux
# O2 6% PAL
#==============================================================================

load("01a_baseline_6percO2.Rdata")

#------------------
# run tests
#------------------

PL <- sens01a$PL
PP_ref   <- 4*1E+18

PP.frac <- c(seq(0.01,0.19,0.01),seq(0.2,0.29999,0.00001),seq(0.3,1,0.1))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens01a$output$C[nrow(sens01a$output$C),]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                  <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance          <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                         "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                         "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01b_baseline_6percO2",type="png")

#------------------
# save output
#------------------

sens01b <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01b,file="01b_baseline_6percO2.Rdata")

##################################################################################
# O2 5% PAL
##################################################################################

#==============================================================================
# varying carbon flux
# O2 5% PAL
#==============================================================================

load("00_CSOFe_runsteady.Rdata")

#------------------
# run tests
#------------------

PL        <- sens00$PL
PP_ref    <- 4*1E+18
PL$k_PyFe <- 0
PL$O2_eq  <- sens00$sens[2851]

PP.frac <- c(seq(1,0.3,-0.1),0.25,seq(0.24999,0.15,-0.00001),seq(0.14,0.01,-0.01))
#PP.frac <- c(seq(1,0.1,-0.05))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens00$output$C[2851,]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                   <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance           <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                          "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                          "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01a_baseline_5percO2",type="png")

#------------------
# save output
#------------------

sens01a <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01a,file="01a_baseline_5percO2.Rdata")

#==============================================================================
# varying carbon flux
# O2 5% PAL
#==============================================================================

load("01a_baseline_5percO2.Rdata")

#------------------
# run tests
#------------------

PL <- sens01a$PL
PP_ref   <- 4*1E+18

PP.frac <- c(seq(0.01,0.14,0.01),seq(0.15,0.24999,0.00001),0.25,seq(0.3,1,0.1))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens01a$output$C[nrow(sens01a$output$C),]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                  <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance          <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                         "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                         "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")


savePlot(filename="01b_baseline_5percO2",type="png")

#------------------
# save output
#------------------

sens01b <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01b,file="01b_baseline_5percO2.Rdata")

##################################################################################
# O2 4% PAL
##################################################################################

#==============================================================================
# varying carbon flux
# O2 4% PAL
#==============================================================================

load("00_CSOFe_runsteady.Rdata")

#------------------
# run tests
#------------------

PL        <- sens00$PL
PP_ref    <- 4*1E+18
PL$k_PyFe <- 0
PL$O2_eq  <- sens00$sens[2881]

PP.frac <- c(seq(1,0.2,-0.1),seq(0.19999,0.05,-0.00001),seq(0.04,0.01,-0.01))
#PP.frac <- c(seq(1,0.1,-0.05))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens00$output$C[2881,]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                   <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance           <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                          "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                          "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01a_baseline_4percO2",type="png")

#------------------
# save output
#------------------

sens01a <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01a,file="01a_baseline_4percO2.Rdata")

#==============================================================================
# varying carbon flux
# O2 4% PAL
#==============================================================================

load("01a_baseline_4percO2.Rdata")

#------------------
# run tests
#------------------

PL <- sens01a$PL
PP_ref   <- 4*1E+18

PP.frac <- c(seq(0.01,0.04,0.01),seq(0.05,0.19999,0.00001),seq(0.2,1,0.1))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens01a$output$C[nrow(sens01a$output$C),]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                  <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance          <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                         "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                         "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01b_baseline_4percO2",type="png")

#------------------
# save output
#------------------

sens01b <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01b,file="01b_baseline_4percO2.Rdata")

##################################################################################
# O2 3% PAL
##################################################################################

#==============================================================================
# varying carbon flux
# O2 3% PAL
#==============================================================================

load("00_CSOFe_runsteady.Rdata")

#------------------
# run tests
#------------------

PL        <- sens00$PL
PP_ref    <- 4*1E+18
PL$k_PyFe <- 0
PL$O2_eq  <- sens00$sens[2911]

PP.frac <- c(seq(1,0.2,-0.1),seq(0.19999,0.05,-0.00001),seq(0.04,0.01,-0.01))
#PP.frac <- c(seq(1,0.1,-0.1),seq(0.1,0.01,-0.01))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens00$output$C[2911,]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                   <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance           <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                          "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                          "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01a_baseline_3percO2",type="png")

#------------------
# save output
#------------------

sens01a <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01a,file="01a_baseline_3percO2.Rdata")

#==============================================================================
# varying carbon flux
# O2 3% PAL
#==============================================================================

load("01a_baseline_3percO2.Rdata")

#------------------
# run tests
#------------------

PL <- sens01a$PL
PP_ref   <- 4*1E+18

PP.frac <- c(seq(0.01,0.04,0.01),seq(0.05,0.19999,0.00001),seq(0.2,1,0.1))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens01a$output$C[nrow(sens01a$output$C),]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                  <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance          <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                         "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                         "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01b_baseline_3percO2",type="png")

#------------------
# save output
#------------------

sens01b <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01b,file="01b_baseline_3percO2.Rdata")

##################################################################################
# O2 2% PAL
##################################################################################

#==============================================================================
# varying carbon flux
# O2 2% PAL
#==============================================================================

load("00_CSOFe_runsteady.Rdata")

#------------------
# run tests
#------------------

PL        <- sens00$PL
PP_ref    <- 4*1E+18
PL$k_PyFe <- 0
PL$O2_eq  <- sens00$sens[2941]

PP.frac <- c(seq(1,0.09,-0.001),seq(0.0899999,0.04,-0.00001),seq(0.03,0.01,-0.01))
#PP.frac <- c(seq(1,0.1,-0.1),seq(0.1,0.01,-0.01))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens00$output$C[2941,]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                   <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance           <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                          "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                          "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste("timestep=",i,i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01a_baseline_2percO2",type="png")

#------------------
# save output
#------------------

sens01a <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01a,file="01a_baseline_2percO2.Rdata")

#==============================================================================
# varying carbon flux
# O2 2% PAL
#==============================================================================

load("01a_baseline_2percO2.Rdata")

#------------------
# run tests
#------------------

PL <- sens01a$PL
PP_ref   <- 4*1E+18

PP.frac <- c(seq(0.01,0.03,0.01),seq(0.04,0.0899999,0.00001),seq(0.09,1,0.001))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens01a$output$C[nrow(sens01a$output$C),]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                  <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance          <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                         "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                         "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01b_baseline_2percO2",type="png")

#------------------
# save output
#------------------

sens01b <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01b,file="01b_baseline_1percO2.Rdata")

##################################################################################
# O2 1% PAL
##################################################################################

load("00_CSOFe_runsteady.Rdata")

#------------------
# run tests
#------------------

PL        <- sens00$PL
PP_ref    <- 4*1E+18
PL$k_PyFe <- 0
PL$O2_eq  <- sens00$sens[2971]

PP.frac <- c(seq(1,0.07,-0.001),seq(0.0699999,0.03,-0.00001),seq(0.02,0.01,-0.01))
#PP.frac <- c(seq(1,0.1,-0.1),seq(0.1,0.01,-0.01))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens00$output$C[2971,]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                   <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance           <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                          "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                          "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01a_baseline_1percO2",type="png")

#------------------
# save output
#------------------

sens01a <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01a,file="01a_baseline_1percO2.Rdata")

#==============================================================================
# varying carbon flux
# O2 1% PAL
#==============================================================================

load("01a_baseline_1percO2.Rdata")

#------------------
# run tests
#------------------

PL <- sens01a$PL
PP_ref   <- 4*1E+18

PP.frac <- c(seq(0.01,0.02,0.01),seq(0.03,0.0699999,0.00001),seq(0.07,1,0.001))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens01a$output$C[nrow(sens01a$output$C),]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                  <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance          <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                         "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                         "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01b_baseline_1percO2",type="png")

#------------------
# save output
#------------------

sens01b <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01b,file="01b_baseline_1percO2.Rdata")

##################################################################################
# O2 0.1% PAL
##################################################################################

load("00_CSOFe_runsteady.Rdata")

#------------------
# run tests
#------------------

PL        <- sens00$PL
PP_ref    <- 4*1E+18
PL$k_PyFe <- 0
PL$O2_eq  <- sens00$sens[2999]

PP.frac <- c(seq(1,0.01,-0.01),seq(0.009999,0.000099,-0.00001))
#PP.frac <- c(seq(1,0.1,-0.1),seq(0.1,0.001,-0.001),0.0001)

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens00$output$C[2999,]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                   <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance           <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                          "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                          "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")


savePlot(filename="01a_baseline_01percO2",type="png")

#------------------
# save output
#------------------

sens01a <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01a,file="01a_baseline_01percO2.Rdata")

#==============================================================================
# varying carbon flux
# O2 0.1% PAL
#==============================================================================

load("01a_baseline_01percO2.Rdata")

#------------------
# run tests
#------------------

PL <- sens01a$PL
PP_ref   <- 4*1E+18

PP.frac <- c(seq(0.000099,0.009999,0.00001),seq(0.01,1,0.01))

PP.sens.output       <- NULL
PP.sens.output.rates <- NULL
dy_y                  <- NULL
mass.balance          <- NULL

SV_new <- sens01a$output$C[nrow(sens01a$output$C),]
for (i in 1:length(PP.frac)){
  
  PL$PP_ref <- PP.frac[i]*PP_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,maxsteps = 1e+9,stol = 1e-10)
  output <- steady.1D(y=output$y,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  
  PP.sens.output        <- rbind(PP.sens.output,output$y)
  PP.sens.output.rates  <- rbind(PP.sens.output.rates,output$other$rate)
  dy_y                  <- rbind(dy_y,(output$other$dy/output$y))
  mass.balance          <- rbind(mass.balance,data.frame("S" =(mass.balance.function(species="S",result=output,PL=PL)/PL$F_SO4_ref),
                                                         "C" =(mass.balance.function(species="C",result=output,PL=PL)/PL$PP_ref),
                                                         "Fe"=(mass.balance.function(species="Fe",result=output,PL=PL)/PL$F_FeOOH_ref)))
  SV_new <- output$y
  
  print(paste(i/length(PP.frac)*100,"% complete"))
  #if(i == length(PP.frac)){beep()}
}

# x11()
# plot(x=mass.balance[,"S"],ylim=c(-0.005,0.005),pch=16)
# points(x=mass.balance[,"Fe"],col="red",pch=16)
# points(x=mass.balance[,"C"],col="dodgerblue",pch=16)

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

xname <- "fraction of PP"

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"O2_D"],H=PP.sens.output[,"O2_H"],L=PP.sens.output[,"O2_L"]),
                 xname=xname,yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"HS_D"],H=PP.sens.output[,"HS_H"],L=PP.sens.output[,"HS_L"]),
                 xname=xname,yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"Fe_D"],H=PP.sens.output[,"Fe_H"],L=PP.sens.output[,"Fe_L"]),
                 xname=xname,yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"SO4_D"],H=PP.sens.output[,"SO4_H"],L=PP.sens.output[,"SO4_L"]),
                 xname=xname,yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output[,"FeOOH_D"],H=PP.sens.output[,"FeOOH_H"],L=PP.sens.output[,"FeOOH_L"]),
                 xname=xname,yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"AR_D"],H=PP.sens.output.rates[,"AR_H"],L=PP.sens.output.rates[,"AR_L"]),
                 xname=xname,yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"SR_D"],H=PP.sens.output.rates[,"SR_H"],L=PP.sens.output.rates[,"SR_L"]),
                 xname=xname,yname="SR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"DIR_D"],H=PP.sens.output.rates[,"DIR_H"],L=PP.sens.output.rates[,"DIR_L"]),
                 xname=xname,yname="DIR  [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"MG_D"],H=PP.sens.output.rates[,"MG_H"],L=PP.sens.output.rates[,"MG_L"]),
                 xname=xname,yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_CH2O"],H=PP.sens.output.rates[,"RR_CH2O_H"],L=PP.sens.output.rates[,"RR_CH2O_L"]),
                 xname=xname,yname="Carbon Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeOOH"],H=PP.sens.output.rates[,"RR_FeOOH_H"],L=PP.sens.output.rates[,"RR_FeOOH_L"]),
                 xname=xname,yname="FeOOH Rain Rate/ Burial in D [mmol yr-1]")

sensitivity.plot(sensitivity=PP.frac,output=data.frame(D=PP.sens.output.rates[,"BR_FeS2"],H=PP.sens.output.rates[,"RR_FeS2_H"],L=PP.sens.output.rates[,"RR_FeS2_L"]),
                 xname=xname,yname="FeS2 Rain Rate/ Burial in D [mmol yr-1]")

savePlot(filename="01b_baseline_01percO2",type="png")

#------------------
# save output
#------------------

sens01b <- list(output=list("C"=PP.sens.output,"R"=PP.sens.output.rates),sens=PP.frac,PL=PL)
save(sens01b,file="01b_baseline_01percO2.Rdata")