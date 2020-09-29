################################################################################
# 3 box model of ocean iron and sulphur chemistry
# Sebastiaan van de Velde & Filip Meysman
# run steady
################################################################################

# Loading required packages

require(AquaEnv)
require(rootSolve)
source(file="CSOFe oceanbox SVDV v05_model only.R")

#==============================================================================
# Run to steady proterozoic (anoxic ocean -> O2_atm = 1% PAL)
#==============================================================================

output$y <- steady.output[1,]

#------------------
# run tests
#------------------

O2.eq.sens <- c(seq(0.3,0.003,-0.0001),seq(0.003,0.0003,-0.0001))
steady.output       <- NULL
steady.output.rates <- NULL
PL$v_sink_mineral <- 63.79

SV_new <- output$y
for (i in 1:length(O2.eq.sens)){
  
  PL$O2_eq <- O2.eq.sens[i]
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,stol = 1e-9)#,nspec=N.vars)
  SV_new <- output$y
  output <- steady.1D(y=SV_new,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  steady.output       <- rbind(steady.output,output$y)
  steady.output.rates <- rbind(steady.output.rates,output$other$rate)
  SV_new <- output$y
  
  print(paste(i/length(O2.eq.sens)*100,"% finished"))
  #print(SV_new)
}

#------------------
# plot output
#------------------

length <- nrow(steady.output)

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

sensitivity.plot(sensitivity=O2.eq.sens[1:length],output=data.frame(D=steady.output[,"O2_D"][1:length],H=steady.output[,"O2_H"][1:length],L=steady.output[,"O2_L"][1:length]),
                 xname="[O2]eq mmol kg-1",yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens[1:length],output=data.frame(D=steady.output[,"HS_D"][1:length],H=steady.output[,"HS_H"][1:length],L=steady.output[,"HS_L"][1:length]),
                 xname="[O2]eq mmol kg-1",yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens[1:length],output=data.frame(D=steady.output[,"Fe_D"][1:length],H=steady.output[,"Fe_H"][1:length],L=steady.output[,"Fe_L"][1:length]),
                 xname="[O2]eq mmol kg-1",yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens[1:length],output=data.frame(D=steady.output[,"SO4_D"][1:length],H=steady.output[,"SO4_H"][1:length],L=steady.output[,"SO4_L"][1:length]),
                 xname="[O2]eq mmol kg-1",yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens[1:length],output=data.frame(D=steady.output[,"FeOOH_D"][1:length],H=steady.output[,"FeOOH_H"][1:length],L=steady.output[,"FeOOH_L"][1:length]),
                 xname="[O2]eq mmol kg-1",yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens[1:length],output=data.frame(D=steady.output[,"FeS2_D"][1:length],H=steady.output[,"FeS2_H"][1:length],L=steady.output[,"FeS2_L"][1:length]),
                 xname="[O2]eq mmol kg-1",yname="FeS2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens[1:length],output=data.frame(D=steady.output.rates[,"AR_D"][1:length],H=steady.output.rates[,"AR_H"][1:length],L=steady.output.rates[,"AR_L"][1:length]),
                 xname="[O2]eq mmol kg-1",yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens[1:length],output=data.frame(D=steady.output.rates[,"SR_D"][1:length],H=steady.output.rates[,"SR_H"][1:length],L=steady.output.rates[,"SR_L"][1:length]),
                 xname="[O2]eq mmol kg-1",yname="SR  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens[1:length],output=data.frame(D=steady.output.rates[,"DIR_D"][1:length],H=steady.output.rates[,"DIR_H"][1:length],L=steady.output.rates[,"DIR_L"][1:length]),
                 xname="[O2]eq mmol kg-1",yname="DIR  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens[1:length],output=data.frame(D=steady.output.rates[,"BR_CH2O"][1:length],H=steady.output.rates[,"RR_CH2O_H"][1:length],L=steady.output.rates[,"RR_CH2O_L"][1:length]),
                 xname="[O2]eq mmol kg-1",yname="Carbon Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens[1:length],output=data.frame(D=steady.output.rates[,"BR_FeOOH"][1:length],H=steady.output.rates[,"RR_FeOOH_H"][1:length],L=steady.output.rates[,"RR_FeOOH_L"][1:length]),
                 xname="[O2]eq mmol kg-1",yname="FeOOH Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens[1:length],output=data.frame(D=steady.output.rates[,"BR_FeS2"][1:length],H=steady.output.rates[,"RR_FeS2_H"][1:length],L=steady.output.rates[,"RR_FeS2_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeS2 Rain Rate/ Burial in D [mol yr-1]")

savePlot(filename="figures/00_CSOFe_runsteady",type="png")

#------------------
# save output
#------------------

sens00 <- list(output=list("C"=steady.output,"R"=steady.output.rates),sens=O2.eq.sens,PL=PL)
save(sens00,file="00_CSOFe_runsteady.Rdata")

#==============================================================================
# Run to steady at O2_atm = 10 % - lower SO4 flux to 100 uM steady-state conc
#==============================================================================

load("00_CSOFe_runsteady.Rdata")

#------------------
# run tests
#------------------

PL        <- sens00$PL
PP_ref    <- 4*1E+18
PL$k_PyFe <- 0
PL$O2_eq  <- sens00$sens[2701]
PL$v_sink_mineral <- 63.79

F_SO4_ref <- 3.3e+15
F.SO4.sens <- seq(1.0,0.35,-0.001)
steady.output       <- NULL
steady.output.rates <- NULL

SV_new <- sens00$output$C[2701,]
for (i in 1:length(F.SO4.sens)){
  
  PL$F_SO4_ref <- F.SO4.sens[i]*F_SO4_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,stol = 1e-9)#,nspec=N.vars)
  SV_new <- output$y
  output <- steady.1D(y=SV_new,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  steady.output       <- rbind(steady.output,output$y)
  steady.output.rates <- rbind(steady.output.rates,output$other$rate)
  SV_new <- output$y
  
  print(paste(i/length(F.SO4.sens)*100,"% finished"))
  #print(SV_new)
}

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output[,"O2_D"],H=steady.output[,"O2_H"],L=steady.output[,"O2_L"]),
                 xname="fraction of FSO4",yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output[,"HS_D"],H=steady.output[,"HS_H"],L=steady.output[,"HS_L"]),
                 xname="fraction of FSO4",yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output[,"Fe_D"],H=steady.output[,"Fe_H"],L=steady.output[,"Fe_L"]),
                 xname="fraction of FSO4",yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output[,"SO4_D"],H=steady.output[,"SO4_H"],L=steady.output[,"SO4_L"]),
                 xname="fraction of FSO4",yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output[,"FeOOH_D"],H=steady.output[,"FeOOH_H"],L=steady.output[,"FeOOH_L"]),
                 xname="fraction of FSO4",yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output.rates[,"AR_D"],H=steady.output.rates[,"AR_H"],L=steady.output.rates[,"AR_L"]),
                 xname="fraction of FSO4",yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output.rates[,"SR_D"],H=steady.output.rates[,"SR_H"],L=steady.output.rates[,"SR_L"]),
                 xname="fraction of FSO4",yname="SR  [mol yr-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output.rates[,"DIR_D"],H=steady.output.rates[,"DIR_H"],L=steady.output.rates[,"DIR_L"]),
                 xname="fraction of FSO4",yname="DIR  [mol yr-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output.rates[,"MG_D"],H=steady.output.rates[,"MG_H"],L=steady.output.rates[,"MG_L"]),
                 xname="fraction of FSO4",yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output.rates[,"BR_CH2O"],H=steady.output.rates[,"RR_CH2O_H"],L=steady.output.rates[,"RR_CH2O_L"]),
                 xname="fraction of FSO4",yname="Carbon Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output.rates[,"BR_FeOOH"],H=steady.output.rates[,"RR_FeOOH_H"],L=steady.output.rates[,"RR_FeOOH_L"]),
                 xname="fraction of FSO4",yname="FeOOH Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output.rates[,"BR_FeS2"],H=steady.output.rates[,"RR_FeS2_H"],L=steady.output.rates[,"RR_FeS2_L"]),
                 xname="fraction of FSO4",yname="FeS2 Rain Rate/ Burial in D [mol yr-1]")

savePlot(filename="figures/00b_CSOFe_runsteady",type="png")

#------------------
# save output
#------------------

sens00 <- list(output=list("C"=steady.output,"R"=steady.output.rates),sens=F.SO4.sens,PL=PL)
save(sens00,file="00b_CSOFe_runsteady.Rdata")

#==============================================================================
# Run to steady at O2_atm = 1 % - lower SO4 flux to 100 uM steady-state conc
#==============================================================================

load("00_CSOFe_runsteady.Rdata")

#------------------
# run tests
#------------------

PL        <- sens00$PL
PL$k_PyFe <- 0
PL$O2_eq  <- sens00$sens[2971]
PL$v_sink_mineral <- 63.79

F_SO4_ref <- 3.3e+15
F.SO4.sens <- seq(1.0,0.35,-0.0001)
steady.output       <- NULL
steady.output.rates <- NULL

SV_new <- sens00$output$C[2971,]
for (i in 1:length(F.SO4.sens)){
  
  PL$F_SO4_ref <- F.SO4.sens[i]*F_SO4_ref
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,stol = 1e-9)#,nspec=N.vars)
  SV_new <- output$y
  output <- steady.1D(y=SV_new,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  steady.output       <- rbind(steady.output,output$y)
  steady.output.rates <- rbind(steady.output.rates,output$other$rate)
  SV_new <- output$y
  
  print(paste(i/length(F.SO4.sens)*100,"% finished"))
  #print(SV_new)
}

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output[,"O2_D"],H=steady.output[,"O2_H"],L=steady.output[,"O2_L"]),
                 xname="fraction of FSO4",yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output[,"HS_D"],H=steady.output[,"HS_H"],L=steady.output[,"HS_L"]),
                 xname="fraction of FSO4",yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output[,"Fe_D"],H=steady.output[,"Fe_H"],L=steady.output[,"Fe_L"]),
                 xname="fraction of FSO4",yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output[,"SO4_D"],H=steady.output[,"SO4_H"],L=steady.output[,"SO4_L"]),
                 xname="fraction of FSO4",yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output[,"FeOOH_D"],H=steady.output[,"FeOOH_H"],L=steady.output[,"FeOOH_L"]),
                 xname="fraction of FSO4",yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output.rates[,"AR_D"],H=steady.output.rates[,"AR_H"],L=steady.output.rates[,"AR_L"]),
                 xname="fraction of FSO4",yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output.rates[,"SR_D"],H=steady.output.rates[,"SR_H"],L=steady.output.rates[,"SR_L"]),
                 xname="fraction of FSO4",yname="SR  [mol yr-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output.rates[,"DIR_D"],H=steady.output.rates[,"DIR_H"],L=steady.output.rates[,"DIR_L"]),
                 xname="fraction of FSO4",yname="DIR  [mol yr-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output.rates[,"MG_D"],H=steady.output.rates[,"MG_H"],L=steady.output.rates[,"MG_L"]),
                 xname="fraction of FSO4",yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output.rates[,"BR_CH2O"],H=steady.output.rates[,"RR_CH2O_H"],L=steady.output.rates[,"RR_CH2O_L"]),
                 xname="fraction of FSO4",yname="Carbon Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output.rates[,"BR_FeOOH"],H=steady.output.rates[,"RR_FeOOH_H"],L=steady.output.rates[,"RR_FeOOH_L"]),
                 xname="fraction of FSO4",yname="FeOOH Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=F.SO4.sens,output=data.frame(D=steady.output.rates[,"BR_FeS2"],H=steady.output.rates[,"RR_FeS2_H"],L=steady.output.rates[,"RR_FeS2_L"]),
                 xname="fraction of FSO41",yname="FeS2 Rain Rate/ Burial in D [mol yr-1]")

savePlot(filename="figures/00cc_CSOFe_runsteady",type="png")

#------------------
# save output
#------------------

sens00 <- list(output=list("C"=steady.output,"R"=steady.output.rates),sens=F.SO4.sens,PL=PL)
save(sens00,file="00cc_CSOFe_runsteady.Rdata")

#==============================================================================
# Run to steady proterozoic (anoxic ocean -> O2_atm = 1% PAL) at Circ = 50perc
#==============================================================================

#------------------
# run tests
#------------------

O2.eq.sens <- c(seq(0.3,0.003,-0.001),seq(0.003,0.0003,-0.0001))
steady.output       <- NULL
steady.output.rates <- NULL
PL$v_sink_mineral <- 63.79
PL$k_PyFe <- 0
PL$UT <- 0.50*25*3.1536e13
PL$UM <- 0.50*30*3.1536e13
PL$UT_L <- PL$UT / PL$V_L  
PL$UT_H <- PL$UT / PL$V_H  
PL$UT_D <- PL$UT / PL$V_D  
PL$UM_H <- PL$UM / PL$V_H  
PL$UM_D <- PL$UM / PL$V_D

SV_new <- output$y
for (i in 1:length(O2.eq.sens)){
  
  PL$O2_eq <- O2.eq.sens[i]
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,stol = 1e-9)#,nspec=N.vars)
  SV_new <- output$y
  output <- steady.1D(y=SV_new,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  steady.output       <- rbind(steady.output,output$y)
  steady.output.rates <- rbind(steady.output.rates,output$other$rate)
  SV_new <- output$y
  
  print(paste(i/length(O2.eq.sens)*100,"% finished"))
  #print(SV_new)
}

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"O2_D"],H=steady.output[,"O2_H"],L=steady.output[,"O2_L"]),
                 xname="[O2]eq mmol kg-1",yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"HS_D"],H=steady.output[,"HS_H"],L=steady.output[,"HS_L"]),
                 xname="[O2]eq mmol kg-1",yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"Fe_D"],H=steady.output[,"Fe_H"],L=steady.output[,"Fe_L"]),
                 xname="[O2]eq mmol kg-1",yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"SO4_D"],H=steady.output[,"SO4_H"],L=steady.output[,"SO4_L"]),
                 xname="[O2]eq mmol kg-1",yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"FeOOH_D"],H=steady.output[,"FeOOH_H"],L=steady.output[,"FeOOH_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"AR_D"],H=steady.output.rates[,"AR_H"],L=steady.output.rates[,"AR_L"]),
                 xname="[O2]eq mmol kg-1",yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"SR_D"],H=steady.output.rates[,"SR_H"],L=steady.output.rates[,"SR_L"]),
                 xname="[O2]eq mmol kg-1",yname="SR  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"DIR_D"],H=steady.output.rates[,"DIR_H"],L=steady.output.rates[,"DIR_L"]),
                 xname="[O2]eq mmol kg-1",yname="DIR  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"MG_D"],H=steady.output.rates[,"MG_H"],L=steady.output.rates[,"MG_L"]),
                 xname="[O2]eq mmol kg-1",yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_CH2O"],H=steady.output.rates[,"RR_CH2O_H"],L=steady.output.rates[,"RR_CH2O_L"]),
                 xname="[O2]eq mmol kg-1",yname="Carbon Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_FeOOH"],H=steady.output.rates[,"RR_FeOOH_H"],L=steady.output.rates[,"RR_FeOOH_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeOOH Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_FeS2"],H=steady.output.rates[,"RR_FeS2_H"],L=steady.output.rates[,"RR_FeS2_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeS2 Rain Rate/ Burial in D [mol yr-1]")

savePlot(filename="00d_CSOFe_runsteady",type="png")

#------------------
# save output
#------------------

sens00 <- list(output=list("C"=steady.output,"R"=steady.output.rates),sens=O2.eq.sens,PL=PL)
save(sens00,file="00d_CSOFe_runsteady.Rdata")

#==============================================================================
# Run to steady proterozoic (anoxic ocean -> O2_atm = 1% PAL) at Circ = 150perc
#==============================================================================

#------------------
# run tests
#------------------

O2.eq.sens <- c(seq(0.3,0.003,-0.001),seq(0.003,0.0003,-0.0001))
steady.output       <- NULL
steady.output.rates <- NULL
PL$v_sink_mineral <- 63.79
PL$k_PyFe <- 0
PL$UT <- 1.5*25*3.1536e13
PL$UM <- 1.5*30*3.1536e13
PL$UT_L <- PL$UT / PL$V_L  
PL$UT_H <- PL$UT / PL$V_H  
PL$UT_D <- PL$UT / PL$V_D  
PL$UM_H <- PL$UM / PL$V_H  
PL$UM_D <- PL$UM / PL$V_D

SV_new <- output$y
for (i in 1:length(O2.eq.sens)){
  
  PL$O2_eq <- O2.eq.sens[i]
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,stol = 1e-9)#,nspec=N.vars)
  SV_new <- output$y
  output <- steady.1D(y=SV_new,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  steady.output       <- rbind(steady.output,output$y)
  steady.output.rates <- rbind(steady.output.rates,output$other$rate)
  SV_new <- output$y
  
  print(paste(i/length(O2.eq.sens)*100,"% finished"))
  #print(SV_new)
}

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"O2_D"],H=steady.output[,"O2_H"],L=steady.output[,"O2_L"]),
                 xname="[O2]eq mmol kg-1",yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"HS_D"],H=steady.output[,"HS_H"],L=steady.output[,"HS_L"]),
                 xname="[O2]eq mmol kg-1",yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"Fe_D"],H=steady.output[,"Fe_H"],L=steady.output[,"Fe_L"]),
                 xname="[O2]eq mmol kg-1",yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"SO4_D"],H=steady.output[,"SO4_H"],L=steady.output[,"SO4_L"]),
                 xname="[O2]eq mmol kg-1",yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"FeOOH_D"],H=steady.output[,"FeOOH_H"],L=steady.output[,"FeOOH_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"AR_D"],H=steady.output.rates[,"AR_H"],L=steady.output.rates[,"AR_L"]),
                 xname="[O2]eq mmol kg-1",yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"SR_D"],H=steady.output.rates[,"SR_H"],L=steady.output.rates[,"SR_L"]),
                 xname="[O2]eq mmol kg-1",yname="SR  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"DIR_D"],H=steady.output.rates[,"DIR_H"],L=steady.output.rates[,"DIR_L"]),
                 xname="[O2]eq mmol kg-1",yname="DIR  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"MG_D"],H=steady.output.rates[,"MG_H"],L=steady.output.rates[,"MG_L"]),
                 xname="[O2]eq mmol kg-1",yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_CH2O"],H=steady.output.rates[,"RR_CH2O_H"],L=steady.output.rates[,"RR_CH2O_L"]),
                 xname="[O2]eq mmol kg-1",yname="Carbon Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_FeOOH"],H=steady.output.rates[,"RR_FeOOH_H"],L=steady.output.rates[,"RR_FeOOH_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeOOH Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_FeS2"],H=steady.output.rates[,"RR_FeS2_H"],L=steady.output.rates[,"RR_FeS2_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeS2 Rain Rate/ Burial in D [mol yr-1]")

savePlot(filename="00e_CSOFe_runsteady",type="png")

#------------------
# save output
#------------------

sens00 <- list(output=list("C"=steady.output,"R"=steady.output.rates),sens=O2.eq.sens,PL=PL)
save(sens00,file="00e_CSOFe_runsteady.Rdata")

#==============================================================================
# Run to steady proterozoic SLOW SETTLING (anoxic ocean -> O2_atm = 1% PAL)
#==============================================================================

#------------------
# run tests
#------------------

O2.eq.sens <- c(seq(0.3,0.003,-0.0001),seq(0.003,0.0003,-0.0001))
steady.output       <- NULL
steady.output.rates <- NULL
PL$k_PyFe <- 0
PL$v_sink_mineral <- 182.625/10

SV_new <- output$y
for (i in 1:length(O2.eq.sens)){
  
  PL$O2_eq <- O2.eq.sens[i]
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,stol = 1e-9)#,nspec=N.vars)
  SV_new <- output$y
  output <- steady.1D(y=SV_new,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  steady.output       <- rbind(steady.output,output$y)
  steady.output.rates <- rbind(steady.output.rates,output$other$rate)
  SV_new <- output$y
  
  print(paste(i/length(O2.eq.sens)*100,"% finished"))
  #print(SV_new)
}

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"O2_D"],H=steady.output[,"O2_H"],L=steady.output[,"O2_L"]),
                 xname="[O2]eq mmol kg-1",yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"HS_D"],H=steady.output[,"HS_H"],L=steady.output[,"HS_L"]),
                 xname="[O2]eq mmol kg-1",yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"Fe_D"],H=steady.output[,"Fe_H"],L=steady.output[,"Fe_L"]),
                 xname="[O2]eq mmol kg-1",yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"SO4_D"],H=steady.output[,"SO4_H"],L=steady.output[,"SO4_L"]),
                 xname="[O2]eq mmol kg-1",yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"FeOOH_D"],H=steady.output[,"FeOOH_H"],L=steady.output[,"FeOOH_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"AR_D"],H=steady.output.rates[,"AR_H"],L=steady.output.rates[,"AR_L"]),
                 xname="[O2]eq mmol kg-1",yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"SR_D"],H=steady.output.rates[,"SR_H"],L=steady.output.rates[,"SR_L"]),
                 xname="[O2]eq mmol kg-1",yname="SR  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"DIR_D"],H=steady.output.rates[,"DIR_H"],L=steady.output.rates[,"DIR_L"]),
                 xname="[O2]eq mmol kg-1",yname="DIR  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"MG_D"],H=steady.output.rates[,"MG_H"],L=steady.output.rates[,"MG_L"]),
                 xname="[O2]eq mmol kg-1",yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_CH2O"],H=steady.output.rates[,"RR_CH2O_H"],L=steady.output.rates[,"RR_CH2O_L"]),
                 xname="[O2]eq mmol kg-1",yname="Carbon Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_FeOOH"],H=steady.output.rates[,"RR_FeOOH_H"],L=steady.output.rates[,"RR_FeOOH_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeOOH Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_FeS2"],H=steady.output.rates[,"RR_FeS2_H"],L=steady.output.rates[,"RR_FeS2_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeS2 Rain Rate/ Burial in D [mol yr-1]")

savePlot(filename="figures/00_CSOFe_runsteady_slowsettling",type="png")

#------------------
# save output
#------------------

sens00 <- list(output=list("C"=steady.output,"R"=steady.output.rates),sens=O2.eq.sens,PL=PL)
save(sens00,file="00_CSOFe_runsteady_slowsettling.Rdata")

#==============================================================================
# Run to steady proterozoic SLOW SETTLING (anoxic ocean -> O2_atm = 1% PAL)
#==============================================================================

#------------------
# run tests
#------------------

O2.eq.sens <- c(seq(0.3,0.003,-0.0001),seq(0.003,0.0003,-0.0001))
steady.output       <- NULL
steady.output.rates <- NULL
PL$k_PyFe <- 0
PL$v_sink_mineral <- 182.625

SV_new <- output$y
for (i in 1:length(O2.eq.sens)){
  
  PL$O2_eq <- O2.eq.sens[i]
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,stol = 1e-9)#,nspec=N.vars)
  SV_new <- output$y
  output <- steady.1D(y=SV_new,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  steady.output       <- rbind(steady.output,output$y)
  steady.output.rates <- rbind(steady.output.rates,output$other$rate)
  SV_new <- output$y
  
  print(paste(i/length(O2.eq.sens)*100,"% finished"))
  #print(SV_new)
}

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"O2_D"],H=steady.output[,"O2_H"],L=steady.output[,"O2_L"]),
                 xname="[O2]eq mmol kg-1",yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"HS_D"],H=steady.output[,"HS_H"],L=steady.output[,"HS_L"]),
                 xname="[O2]eq mmol kg-1",yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"Fe_D"],H=steady.output[,"Fe_H"],L=steady.output[,"Fe_L"]),
                 xname="[O2]eq mmol kg-1",yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"SO4_D"],H=steady.output[,"SO4_H"],L=steady.output[,"SO4_L"]),
                 xname="[O2]eq mmol kg-1",yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"FeOOH_D"],H=steady.output[,"FeOOH_H"],L=steady.output[,"FeOOH_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"AR_D"],H=steady.output.rates[,"AR_H"],L=steady.output.rates[,"AR_L"]),
                 xname="[O2]eq mmol kg-1",yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"SR_D"],H=steady.output.rates[,"SR_H"],L=steady.output.rates[,"SR_L"]),
                 xname="[O2]eq mmol kg-1",yname="SR  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"DIR_D"],H=steady.output.rates[,"DIR_H"],L=steady.output.rates[,"DIR_L"]),
                 xname="[O2]eq mmol kg-1",yname="DIR  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"MG_D"],H=steady.output.rates[,"MG_H"],L=steady.output.rates[,"MG_L"]),
                 xname="[O2]eq mmol kg-1",yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_CH2O"],H=steady.output.rates[,"RR_CH2O_H"],L=steady.output.rates[,"RR_CH2O_L"]),
                 xname="[O2]eq mmol kg-1",yname="Carbon Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_FeOOH"],H=steady.output.rates[,"RR_FeOOH_H"],L=steady.output.rates[,"RR_FeOOH_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeOOH Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_FeS2"],H=steady.output.rates[,"RR_FeS2_H"],L=steady.output.rates[,"RR_FeS2_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeS2 Rain Rate/ Burial in D [mol yr-1]")

savePlot(filename="figures/00_CSOFe_runsteady_fastsettling",type="png")

#------------------
# save output
#------------------

sens00 <- list(output=list("C"=steady.output,"R"=steady.output.rates),sens=O2.eq.sens,PL=PL)
save(sens00,file="00_CSOFe_runsteady_fastsettling.Rdata")

#==============================================================================
# Run to steady proterozoic F_FeOOH/3 (anoxic ocean -> O2_atm = 1% PAL)
#==============================================================================

#------------------
# run tests
#------------------

O2.eq.sens <- c(seq(0.3,0.003,-0.0001),seq(0.003,0.0003,-0.0001))
steady.output       <- NULL
steady.output.rates <- NULL
PL$k_PyFe <- 0
PL$v_sink_mineral <- 63.79

PL$F_FeOOH_ref <- 0.3*0.6*1E15

SV_new <- output$y
for (i in 1:length(O2.eq.sens)){
  
  PL$O2_eq <- O2.eq.sens[i]
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,stol = 1e-9)#,nspec=N.vars)
  SV_new <- output$y
  output <- steady.1D(y=SV_new,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  steady.output       <- rbind(steady.output,output$y)
  steady.output.rates <- rbind(steady.output.rates,output$other$rate)
  SV_new <- output$y
  
  print(paste(i/length(O2.eq.sens)*100,"% finished"))
  #print(SV_new)
}

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"O2_D"],H=steady.output[,"O2_H"],L=steady.output[,"O2_L"]),
                 xname="[O2]eq mmol kg-1",yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"HS_D"],H=steady.output[,"HS_H"],L=steady.output[,"HS_L"]),
                 xname="[O2]eq mmol kg-1",yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"Fe_D"],H=steady.output[,"Fe_H"],L=steady.output[,"Fe_L"]),
                 xname="[O2]eq mmol kg-1",yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"SO4_D"],H=steady.output[,"SO4_H"],L=steady.output[,"SO4_L"]),
                 xname="[O2]eq mmol kg-1",yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"FeOOH_D"],H=steady.output[,"FeOOH_H"],L=steady.output[,"FeOOH_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"AR_D"],H=steady.output.rates[,"AR_H"],L=steady.output.rates[,"AR_L"]),
                 xname="[O2]eq mmol kg-1",yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"SR_D"],H=steady.output.rates[,"SR_H"],L=steady.output.rates[,"SR_L"]),
                 xname="[O2]eq mmol kg-1",yname="SR  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"DIR_D"],H=steady.output.rates[,"DIR_H"],L=steady.output.rates[,"DIR_L"]),
                 xname="[O2]eq mmol kg-1",yname="DIR  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"MG_D"],H=steady.output.rates[,"MG_H"],L=steady.output.rates[,"MG_L"]),
                 xname="[O2]eq mmol kg-1",yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_CH2O"],H=steady.output.rates[,"RR_CH2O_H"],L=steady.output.rates[,"RR_CH2O_L"]),
                 xname="[O2]eq mmol kg-1",yname="Carbon Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_FeOOH"],H=steady.output.rates[,"RR_FeOOH_H"],L=steady.output.rates[,"RR_FeOOH_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeOOH Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_FeS2"],H=steady.output.rates[,"RR_FeS2_H"],L=steady.output.rates[,"RR_FeS2_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeS2 Rain Rate/ Burial in D [mol yr-1]")

savePlot(filename="figures/00_CSOFe_runsteady_FFeOOHdec",type="png")

#------------------
# save output
#------------------

sens00 <- list(output=list("C"=steady.output,"R"=steady.output.rates),sens=O2.eq.sens,PL=PL)
save(sens00,file="00_CSOFe_runsteady_FFeOOHdec.Rdata")

#==============================================================================
# Run to steady proterozoic F_FeOOHx2 (anoxic ocean -> O2_atm = 1% PAL)
#==============================================================================

#------------------
# run tests
#------------------

O2.eq.sens <- c(seq(0.3,0.003,-0.0001),seq(0.003,0.0003,-0.0001))
steady.output       <- NULL
steady.output.rates <- NULL
PL$k_PyFe <- 0
PL$v_sink_mineral <- 63.79

PL$F_FeOOH_ref <- 2.0*0.6*1E15

SV_new <- output$y
for (i in 1:length(O2.eq.sens)){
  
  PL$O2_eq <- O2.eq.sens[i]
  output <- runsteady(y=SV_new,func=boxmodel.CSOFe,parms=PL,stol = 1e-9)#,nspec=N.vars)
  SV_new <- output$y
  output <- steady.1D(y=SV_new,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
  steady.output       <- rbind(steady.output,output$y)
  steady.output.rates <- rbind(steady.output.rates,output$other$rate)
  SV_new <- output$y
  
  print(paste(i/length(O2.eq.sens)*100,"% finished"))
  #print(SV_new)
}

#------------------
# plot output
#------------------

win.graph(width=8/3*4,height=8,pointsize=10)

par(mfrow=c(3,4),oma=c(2,1,0,1))

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"O2_D"],H=steady.output[,"O2_H"],L=steady.output[,"O2_L"]),
                 xname="[O2]eq mmol kg-1",yname="O2 conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"HS_D"],H=steady.output[,"HS_H"],L=steady.output[,"HS_L"]),
                 xname="[O2]eq mmol kg-1",yname="HS conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"Fe_D"],H=steady.output[,"Fe_H"],L=steady.output[,"Fe_L"]),
                 xname="[O2]eq mmol kg-1",yname="Fe conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"SO4_D"],H=steady.output[,"SO4_H"],L=steady.output[,"SO4_L"]),
                 xname="[O2]eq mmol kg-1",yname="SO4 conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output[,"FeOOH_D"],H=steady.output[,"FeOOH_H"],L=steady.output[,"FeOOH_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeOOH conc [mmol kg-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"AR_D"],H=steady.output.rates[,"AR_H"],L=steady.output.rates[,"AR_L"]),
                 xname="[O2]eq mmol kg-1",yname="AR [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"SR_D"],H=steady.output.rates[,"SR_H"],L=steady.output.rates[,"SR_L"]),
                 xname="[O2]eq mmol kg-1",yname="SR  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"DIR_D"],H=steady.output.rates[,"DIR_H"],L=steady.output.rates[,"DIR_L"]),
                 xname="[O2]eq mmol kg-1",yname="DIR  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"MG_D"],H=steady.output.rates[,"MG_H"],L=steady.output.rates[,"MG_L"]),
                 xname="[O2]eq mmol kg-1",yname="MG  [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_CH2O"],H=steady.output.rates[,"RR_CH2O_H"],L=steady.output.rates[,"RR_CH2O_L"]),
                 xname="[O2]eq mmol kg-1",yname="Carbon Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_FeOOH"],H=steady.output.rates[,"RR_FeOOH_H"],L=steady.output.rates[,"RR_FeOOH_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeOOH Rain Rate/ Burial in D [mol yr-1]")

sensitivity.plot(sensitivity=O2.eq.sens,output=data.frame(D=steady.output.rates[,"BR_FeS2"],H=steady.output.rates[,"RR_FeS2_H"],L=steady.output.rates[,"RR_FeS2_L"]),
                 xname="[O2]eq mmol kg-1",yname="FeS2 Rain Rate/ Burial in D [mol yr-1]")

savePlot(filename="figures/00_CSOFe_runsteady_FFeOOHinc",type="png")

#------------------
# save output
#------------------

sens00 <- list(output=list("C"=steady.output,"R"=steady.output.rates),sens=O2.eq.sens,PL=PL)
save(sens00,file="00_CSOFe_runsteady_FFeOOHinc.Rdata")
