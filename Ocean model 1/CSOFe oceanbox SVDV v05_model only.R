################################################################################
# 3 box model of ocean iron and sulphur chemistry
# Sebastiaan van de Velde & Filip Meysman
################################################################################

# Loading required packages

require(AquaEnv)
require(rootSolve)

N.vars <- 24

#===============================================================================
# Units
#===============================================================================

# mmol
# m
# yr
# kg

#===============================================================================
# Parameters
#===============================================================================

# L = low latitude surface ocean
# H = high latitude surface ocean
# D = deep ocean

PL <- list()

# -------------------------------------------
# Geometry of the ocean box model
# -------------------------------------------

PL$Aoc <- 3.35*1E14   # [m^2] surface area of the oceans (excluding continental shelves)  ETOPO1
PL$Voc <- 1.3361*1E18 # [m^3] volume of the oceans (excluding continental shelves)  ETOPO1
PL$z_b <- 10898       # [m] "bottom horizon": 10898 m in ETOPO1                                                                   
PL$a_L <- 0.85        # [-] fraction of total ocean surface in L box Toggweiler2003 (Boudreau2009, Sarmiento1984, many old refs: 0.85)
PL$h_L <- 100         # [m] depth of box L Sarmiento1984
PL$h_H <- 350         # [m] depth of box H Sarmiento1984

PL$V_L <- PL$h_L*PL$a_L*PL$Aoc      # [m^3] volume of box L
PL$V_H <- PL$h_H*(1-PL$a_L)*PL$Aoc  # [m^3] volume of box H
PL$V_D <- PL$Voc - PL$V_L - PL$V_H  # [m^3] volume of box D                                
PL$A_D <- PL$Aoc                    # [m^2] surface area of box D
PL$A_L <- PL$a_L*PL$Aoc             # [m^2] surface area of box L
PL$A_H <- (1-PL$a_L)*PL$Aoc         # [m^2] surface area of box H
PL$h_D <- PL$V_D/PL$A_D             # [m] mean depth of box D

# -------------------------------------------
# Environmental parameters in each ocean box 
# -------------------------------------------

# Temperature and salinity

PL$t_L <- 21.5 # [deg C]   temperature of box L
PL$t_H <- 2    # [deg C]   temperature of box H
PL$t_D <- 2    # [deg C]   temperature of box D
PL$S_L <- 35   #  -      salinity of box L
PL$S_H <- 35   #  -      salinity of box H
PL$S_D <- 35   #  -      salinity of box D

# Pressure [bar] (via function gauge_p from AquaEnv package) 
# Note: gauge pressure = hydrostatic pressure - atmospheric pressure

PL$p_L <- gauge_p(PL$h_L/2)  # [bar] mean gauge pressure of box L
PL$p_H <- gauge_p(PL$h_H/2)  # [bar] mean gauge pressure of box H                           
PL$p_D <- gauge_p(PL$h_D/2 + max(PL$h_H,PL$h_L))  # [bar] mean gauge pressure of box D                           

# Density [kg m-3] (via function aquaenv from AquaEnv package) 

PL$dens_L <- aquaenv(S=PL$S_L, t=PL$t_L, p=PL$p_L)$density[[1]]                           
PL$dens_H <- aquaenv(S=PL$S_H, t=PL$t_H, p=PL$p_H)$density[[1]]  
PL$dens_D <- aquaenv(S=PL$S_D, t=PL$t_D, p=PL$p_D)$density[[1]]  

# Total mass in each box [kg]

PL$m_L <- PL$dens_L * PL$V_L
PL$m_H <- PL$dens_H * PL$V_H
PL$m_D <- PL$dens_D * PL$V_D

# -------------------------------------------
# Biogeochemical rate parameters 
# -------------------------------------------

# Circulation parameters 

PL$UT     <- 25 *3.1536e13 # [Sv]  flow rate of thermohaline circulation     Boudreau2009: extracted from code; paper: 24.1
PL$UM     <- 30 *3.1536e13 # [Sv]  mixing coefficient (H <-> D)              Boudreau2009

PL$v_sink <- 182.625       # [kg m-2 yr-1] (from hypoxia problem of modelling of the environment)
PL$v_sink_mineral <- 182.625  

PL$kl      <- 4.8 *365      # [m yr-1] piston velocity CO2 air-sea exchange (Sarmiento&Gruber 2006)
PL$O2_eq   <- 0.3           # [mmol kg-1] O2 equilibration concentration
PL$CH4_eq  <- 0.0           # [mmol kg-1] O2 equilibration concentration

# Unit conversions 

PL$UT_L <- PL$UT / PL$V_L  # [1/yr]  (transport) proportionality factor for thermohaline circulation, box L
PL$UT_H <- PL$UT / PL$V_H  # [1/yr]  (transport) proportionality factor for thermohaline circulation, box H
PL$UT_D <- PL$UT / PL$V_D  # [1/yr]  (transport) proportionality factor for thermohaline circulation, box D

PL$UM_H <- PL$UM / PL$V_H  # [1/yr]  (transport) proportionality factor for higher latitude mixing, box H
PL$UM_D <- PL$UM / PL$V_D  # [1/yr]  (transport) proportionality factor for higher latitude mixing, box D

# Supply of iron, sulphate, organic matter and oxygen to box L

PL$F_FeOOH_ref <- 0.6*1E15 # [mmol yr-1] (Canfield and Farquhar, 2009): based on the 1.2*1E12 mol yr-1 modern pyrite burial rate 
# matches 1.3*1E12 mol yr-1 input in ocean, assuming only 50% is reactive and rest is refractory
# [mmol yr-1] (Poulton and Raiswell, 2002): total 6.8*1e12 -> 1.3*1E12 mol yr-1 in ocean without inne shore areas  
PL$F_SO4_ref   <- 3.3*1E15 # [mmol yr-1] (Poulton and Canfield, 2011)
PL$PP_ref      <- 4*1E+18  # [mmol yr-1] (Longhurst et al., 1995) & (Antoine et al., 1996)
PL$f_photo     <- 1        # - 1 mole C yields 1 mole oxygen (other estimates are 1.4 mol O2 per 1 mol C)
PL$fPP_L       <- 0.85     # - fraction of primary production that occurs in box L
PL$fPP_H       <- 0.15     # - fraction of primary production that occurs in box H
PL$SO4_D_ref   <- 26       # [mmol kg-1] present day sulphate concentration 
PL$gypsum      <- 2.1*1E15 # [mmol yr-1] present day evaporite formation
PL$k_gypsum    <- PL$F_SO4_ref/PL$SO4_D_ref # [mmol yr-1] calibrated as F_gypsum=k_gypsum*[SO4]

# reaction rate parameters -> need to be calibrated, too high values lead to numerical instability

PL$k_min   <- 0.05        # [yr-1] (from previous sensitivity test, to get 200*1e15 mmol yr-1 OM export from box L to box D = 5 % of total PP)
PL$K_O2    <- 0.001       # [mmol kg-1] (van de Velde and Meysman, 2016)
PL$K_FeOOH <- 10.4        # [mmol kg-1] (van de Velde and Meysman, 2016)
PL$K_SO4   <- 0.9         # [mmol kg-1] (van de Velde and Meysman, 2016)

PL$k_FIO   <- 1e2#1e7   # [mmol-1 m3 yr-1] (Millero et al., 1987)
PL$k_CSO   <- 1e2#1e7   # [mmol-1 kg yr-1] (Millero et al., 1987)
PL$k_SIR   <- 1e2#494   # [mmol-1 m3 yr-1] (Poulton et al., 2004)
PL$k_PyO   <- 1e2#1e7   # [mmol-1 m3 yr-1] (Meysman et al., 2015)

PL$k_PyFe  <- 0#1e2   # [mmol-1 m3 yr-1] zelf uitgevonden
PL$k_PyS   <- 1e2#1e4   # [mmol-1 m3 yr-1] zelf uitgevonden

PL$k_AMO   <- 1e2 # model stability?
PL$k_AnMO  <- 1e2 # model stability?
  
  
#===============================================================================
# Auxilliary functions
#===============================================================================

FSAT <- function(C,K=1e-6,n) (C/K)^n/((C/K)^n + 1)

# -------------------------------------------
# plot function
# -------------------------------------------

sensitivity.plot <- function(sensitivity,output,xname,xlim=NULL,yname,ylim=NULL,logy=FALSE){
  
  if (is.null(ylim)){
    ylim <- range(output)
  }
  
  if (is.null(xlim)){
    xlim <- range(sensitivity)
  }
  
  if (logy==F){
  par(mar=c(2,3,2,2))
  
  plot(x=sensitivity,y=output[,"D"],
       axes=FALSE, frame.plot=FALSE, 
       ylim=ylim, ylab="",cex.lab=1.2,
       xlim=xlim, xlab="",type="p",
       col= "black", lwd=2, cex=1,pch=16)
  axis(1, cex.axis=1);  mtext(side=1, text=xname, line=3,cex=1)
  axis(2, cex.axis=1);  mtext(side=2, text=yname, line=2,cex=1)
  
  points(x=sensitivity,y=output[,"L"],col="dodgerblue",pch=15)
  points(x=sensitivity,y=output[,"H"],col="red",pch=17)
  
  legend("right",bty="n",legend=c("box D","box L","box H"),pch=c(16,15,17),col=c("black","dodgerblue","red"))
  }
  
  if (logy==T){
    par(mar=c(2,3,2,2))
    
    plot(x=sensitivity,y=output[,"D"],
         axes=FALSE, frame.plot=FALSE, 
         ylim=ylim, ylab="",cex.lab=1.2,
         xlim=xlim, xlab="",type="p",
         col= "black", lwd=2, cex=1,pch=16,log="y")
    axis(1, cex.axis=1);  mtext(side=1, text=xname, line=3,cex=1)
    axis(2, cex.axis=1);  mtext(side=2, text=yname, line=2,cex=1)
    
    points(x=sensitivity,y=output[,"L"],col="dodgerblue",pch=15,log="y")
    points(x=sensitivity,y=output[,"H"],col="red",pch=17,log="y")
    
    legend("right",bty="n",legend=c("box D","box L","box H"),pch=c(16,15,17),col=c("black","dodgerblue","red"))
  }
}

# -------------------------------------------
# mass balance check
# -------------------------------------------

mass.balance.function <- function(species,result,PL){
  if(species=="S"){
    
    input  <- PL$F_SO4_ref
    output <- result$other$rate["BR_SO4"][[1]] + 2*result$other$rate["BR_FeS2"][[1]]
  
    deficit <- input-output
    return(deficit)
    }
  
  if(species=="Fe"){
    
    input  <- PL$F_FeOOH_ref
    output <- result$other$rate["BR_FeOOH"][[1]] + result$other$rate["BR_FeS2"][[1]]
    
    deficit <- input-output
    return(deficit)
  }
  
  if(species=="C"){
    
    input  <- PL$PP_ref
    output <- result$other$rate["BR_CH2O"][[1]] + 
      (result$other$rate["AR_L"][[1]]  + result$other$rate["AR_H"][[1]]  + result$other$rate["AR_D"][[1]]  + 
       result$other$rate["DIR_L"][[1]] + result$other$rate["DIR_H"][[1]] + result$other$rate["DIR_D"][[1]] +
       result$other$rate["SR_L"][[1]]  + result$other$rate["SR_H"][[1]]  + result$other$rate["SR_D"][[1]]  )
    
    deficit <- input-output
    return(deficit)
  }
  
  }

#===============================================================================
# Model formulation OCEANBOX
#===============================================================================

boxmodel.CSOFe <- function(time, initialstate, parameters)
{
 with(as.list(c(initialstate, parameters)),{       

#===============================================================================
# Transport equations
#===============================================================================

# O2 air-sea exchange
   
 E_L <- kl*A_L*dens_L*(O2_eq - O2_L) # [mmol yr-1]
 E_H <- kl*A_H*dens_H*(O2_eq - O2_H) # [mmol yr-1]
 
 CH4_E_L <- kl*A_L*dens_L*(CH4_eq - CH4_L) # [mmol yr-1]
 CH4_E_H <- kl*A_H*dens_H*(CH4_eq - CH4_H) # [mmol yr-1]
 
# thermohaline circulation
 
 R_UTO2_L <- UT_L * (O2_D*(dens_D/dens_L) - O2_L)
 R_UTO2_H <- UT_H * (O2_L*(dens_L/dens_H) - O2_H)
 R_UTO2_D <- UT_D * (O2_H*(dens_H/dens_D) - O2_D) 

 R_UTSO4_L <- UT_L * (SO4_D*(dens_D/dens_L) - SO4_L)
 R_UTSO4_H <- UT_H * (SO4_L*(dens_L/dens_H) - SO4_H)
 R_UTSO4_D <- UT_D * (SO4_H*(dens_H/dens_D) - SO4_D)
 
 R_UTFe_L <- UT_L * (Fe_D*(dens_D/dens_L) - Fe_L)
 R_UTFe_H <- UT_H * (Fe_L*(dens_L/dens_H) - Fe_H)
 R_UTFe_D <- UT_D * (Fe_H*(dens_H/dens_D) - Fe_D)

 R_UTHS_L <- UT_L * (HS_D*(dens_D/dens_L) - HS_L)
 R_UTHS_H <- UT_H * (HS_L*(dens_L/dens_H) - HS_H)
 R_UTHS_D <- UT_D * (HS_H*(dens_H/dens_D) - HS_D)

 R_UTCH4_L <- UT_L * (CH4_D*(dens_D/dens_L) - CH4_L)
 R_UTCH4_H <- UT_H * (CH4_L*(dens_L/dens_H) - CH4_H)
 R_UTCH4_D <- UT_D * (CH4_H*(dens_H/dens_D) - CH4_D)
 
# high latitude mixing
 
 R_UMO2_H <- UM_H * (O2_D*(dens_D/dens_H) - O2_H)
 R_UMO2_D <- UM_D * (O2_H*(dens_H/dens_D) - O2_D)
 
 R_UMSO4_H <- UM_H * (SO4_D*(dens_D/dens_H) - SO4_H)
 R_UMSO4_D <- UM_D * (SO4_H*(dens_H/dens_D) - SO4_D)
 
 R_UMFe_H <- UM_H * (Fe_D*(dens_D/dens_H) - Fe_H)
 R_UMFe_D <- UM_D * (Fe_H*(dens_H/dens_D) - Fe_D)
 
 R_UMHS_H <- UM_H * (HS_D*(dens_D/dens_H) - HS_H)
 R_UMHS_D <- UM_D * (HS_H*(dens_H/dens_D) - HS_D)

 R_UMCH4_H <- UM_H * (CH4_D*(dens_D/dens_H) - CH4_H)
 R_UMCH4_D <- UM_D * (CH4_H*(dens_H/dens_D) - CH4_D)

 # Riverine input of FeOOH and SO into L box
 
 F_FeOOH_L <- F_FeOOH_ref          
 F_SO4_L   <- F_SO4_ref          
 
 # Primary production simulated as CH2O and O2 flux into L box
 
 F_CH2O_L <- fPP_L*PP_ref          
 F_O2_L   <- f_photo*fPP_L*PP_ref        
 
 F_CH2O_H <- fPP_H*PP_ref          
 F_O2_H   <- f_photo*fPP_H*PP_ref 
 
 # rain of solid particles
 
 F_CH2O_LD  <- v_sink           *CH2O_L *A_L*(CH2O_L>0)
 F_FeOOH_LD <- v_sink_mineral   *FeOOH_L*A_L*(FeOOH_L>0)
 F_FeS2_LD  <- v_sink_mineral*FeS2_L *A_L*(FeS2_L>0)

 F_CH2O_HD  <- v_sink           *CH2O_H *A_H*(CH2O_H>0)
 F_FeOOH_HD <- v_sink_mineral*FeOOH_H*A_H*(FeOOH_H>0)
 F_FeS2_HD  <- v_sink_mineral*FeS2_H *A_H*(FeS2_H>0)
 
 # sediment burial of solid particles
 
 F_CH2O_burial  <- v_sink           *CH2O_D *A_D*(CH2O_D>0)
 F_FeOOH_burial <- v_sink_mineral*FeOOH_D*A_D*(FeOOH_D>0)
 F_FeS2_burial  <- v_sink_mineral*FeS2_D *A_D*(FeS2_D>0)
 
 # gypsum formation
 
 F_SO4_gyp <- k_gypsum*SO4_D*(SO4_D>0) # lumps gypsum formation and sedimentary pyrite formation to get SO4.conc to present day values 
 
#===============================================================================
# Reaction equations
#===============================================================================

# Mineralisation reactions 

 R.min_L <- k_min*CH2O_L
 R.min_D <- k_min*CH2O_D
 R.min_H <- k_min*CH2O_H
 
 O2_lim_L <- O2_L/(O2_L+K_O2)#*(O2_L>0)
 O2_inh_L <- K_O2/(O2_L+K_O2)
 O2_lim_D <- O2_D/(O2_D+K_O2)#*(O2_D>0)
 O2_inh_D <- K_O2/(O2_D+K_O2)
 O2_lim_H <- O2_H/(O2_H+K_O2)#*(O2_H>0)
 O2_inh_H <- K_O2/(O2_H+K_O2)
 
 FeOOH_lim_L <- FeOOH_L/(K_FeOOH+FeOOH_L)#*FSAT(C=FeOOH_L,n=5)
 FeOOH_inh_L <- K_FeOOH/(K_FeOOH+FeOOH_L)#*FSAT(C=FeOOH_L,n=5)
 FeOOH_lim_D <- FeOOH_D/(K_FeOOH+FeOOH_D)#*FSAT(C=FeOOH_D,n=5)
 FeOOH_inh_D <- K_FeOOH/(K_FeOOH+FeOOH_D)#*FSAT(C=FeOOH_D,n=5)
 FeOOH_lim_H <- FeOOH_H/(K_FeOOH+FeOOH_H)#*FSAT(C=FeOOH_H,n=5)
 FeOOH_inh_H <- K_FeOOH/(K_FeOOH+FeOOH_H)#*FSAT(C=FeOOH_H,n=5)
 
 SO4_lim_L <- SO4_L/(K_SO4+SO4_L)#*FSAT(C=SO4_L,n=5)
 SO4_inh_L <- K_SO4/(K_SO4+SO4_L)#*FSAT(C=SO4_L,n=5)
 SO4_lim_D <- SO4_D/(K_SO4+SO4_D)#*FSAT(C=SO4_D,n=5)
 SO4_inh_D <- K_SO4/(K_SO4+SO4_D)#*FSAT(C=SO4_D,n=5)
 SO4_lim_H <- SO4_H/(K_SO4+SO4_H)#*FSAT(C=SO4_H,n=5)
 SO4_inh_H <- K_SO4/(K_SO4+SO4_H)#*FSAT(C=SO4_H,n=5)
 
 f_AR_L  <- O2_lim_L/(O2_lim_L + FeOOH_lim_L*O2_inh_L + SO4_lim_L*FeOOH_inh_L*O2_inh_L + SO4_inh_L*FeOOH_inh_L*O2_inh_L)
 f_AR_D  <- O2_lim_D/(O2_lim_D + FeOOH_lim_D*O2_inh_D + SO4_lim_D*FeOOH_inh_D*O2_inh_D + SO4_inh_D*FeOOH_inh_D*O2_inh_D)
 f_AR_H  <- O2_lim_H/(O2_lim_H + FeOOH_lim_H*O2_inh_H + SO4_lim_H*FeOOH_inh_H*O2_inh_H + SO4_inh_H*FeOOH_inh_H*O2_inh_H)
 f_DIR_L <- (FeOOH_lim_L*O2_inh_L)/(O2_lim_L + FeOOH_lim_L*O2_inh_L + SO4_lim_L*FeOOH_inh_L*O2_inh_L + SO4_inh_L*FeOOH_inh_L*O2_inh_L)
 f_DIR_D <- (FeOOH_lim_D*O2_inh_D)/(O2_lim_D + FeOOH_lim_D*O2_inh_D + SO4_lim_D*FeOOH_inh_D*O2_inh_D + SO4_inh_D*FeOOH_inh_D*O2_inh_D)
 f_DIR_H <- (FeOOH_lim_H*O2_inh_H)/(O2_lim_H + FeOOH_lim_H*O2_inh_H + SO4_lim_H*FeOOH_inh_H*O2_inh_H + SO4_inh_H*FeOOH_inh_H*O2_inh_H)
 f_SR_L  <- (SO4_lim_L*FeOOH_inh_L*O2_inh_L)/(O2_lim_L + FeOOH_lim_L*O2_inh_L + SO4_lim_L*FeOOH_inh_L*O2_inh_L + SO4_inh_L*FeOOH_inh_L*O2_inh_L)
 f_SR_D  <- (SO4_lim_D*FeOOH_inh_D*O2_inh_D)/(O2_lim_D + FeOOH_lim_D*O2_inh_D + SO4_lim_D*FeOOH_inh_D*O2_inh_D + SO4_inh_D*FeOOH_inh_D*O2_inh_D)
 f_SR_H  <- (SO4_lim_H*FeOOH_inh_H*O2_inh_H)/(O2_lim_H + FeOOH_lim_H*O2_inh_H + SO4_lim_H*FeOOH_inh_H*O2_inh_H + SO4_inh_H*FeOOH_inh_H*O2_inh_H)
 f_MG_L  <- (SO4_inh_L*FeOOH_inh_L*O2_inh_L)/(O2_lim_L + FeOOH_lim_L*O2_inh_L + SO4_lim_L*FeOOH_inh_L*O2_inh_L + SO4_inh_L*FeOOH_inh_L*O2_inh_L)
 f_MG_D  <- (SO4_inh_D*FeOOH_inh_D*O2_inh_D)/(O2_lim_D + FeOOH_lim_D*O2_inh_D + SO4_lim_D*FeOOH_inh_D*O2_inh_D + SO4_inh_D*FeOOH_inh_D*O2_inh_D)
 f_MG_H  <- (SO4_inh_H*FeOOH_inh_H*O2_inh_H)/(O2_lim_H + FeOOH_lim_H*O2_inh_H + SO4_lim_H*FeOOH_inh_H*O2_inh_H + SO4_inh_H*FeOOH_inh_H*O2_inh_H)
 
 AR_L <- f_AR_L*R.min_L*(O2_L>0)
 AR_D <- f_AR_D*R.min_D*(O2_D>0)
 AR_H <- f_AR_H*R.min_H*(O2_H>0)
 
 DIR_L <-  f_DIR_L*R.min_L*(FeOOH_L>0)
 DIR_D <-  f_DIR_D*R.min_D*(FeOOH_D>0)
 DIR_H <-  f_DIR_H*R.min_H*(FeOOH_H>0)
 
 SR_L <-  f_SR_L*R.min_L*(SO4_L>0)
 SR_D <-  f_SR_D*R.min_D*(SO4_D>0)
 SR_H <-  f_SR_H*R.min_H*(SO4_H>0)
 
 MG_L <-  f_MG_L*R.min_L*(CH2O_L>0)
 MG_D <-  f_MG_D*R.min_D*(CH2O_D>0)
 MG_H <-  f_MG_H*R.min_H*(CH2O_H>0)
 
 # Reoxidation reactions 
 
 FIO_L <- k_FIO*O2_L*Fe_L*(Fe_L>0)*(O2_L>0)
 FIO_D <- k_FIO*O2_D*Fe_D*(Fe_D>0)*(O2_D>0)
 FIO_H <- k_FIO*O2_H*Fe_H*(Fe_H>0)*(O2_H>0)
 
 CSO_L <- k_CSO*O2_L*HS_L*(HS_L>0)*(O2_L>0)
 CSO_D <- k_CSO*O2_D*HS_D*(HS_D>0)*(O2_D>0)
 CSO_H <- k_CSO*O2_H*HS_H*(HS_H>0)*(O2_H>0)
 
 AMO_L <- k_AMO*O2_L*CH4_L*(CH4_L>0)*(O2_L>0)
 AMO_D <- k_AMO*O2_D*CH4_D*(CH4_D>0)*(O2_D>0)
 AMO_H <- k_AMO*O2_H*CH4_H*(CH4_H>0)*(O2_H>0)
 
 AnMO_L <- k_AnMO*SO4_L*CH4_L*(CH4_L>0)*(O2_L>0)
 AnMO_D <- k_AnMO*SO4_D*CH4_D*(CH4_D>0)*(O2_D>0)
 AnMO_H <- k_AnMO*SO4_H*CH4_H*(CH4_H>0)*(O2_H>0)

 # Fe-S interaction reactions 
 
 SIR_L <- k_SIR*FeOOH_L*HS_L*(HS_L>0)*(FeOOH_L>0)
 SIR_D <- k_SIR*FeOOH_D*HS_D*(HS_D>0)*(FeOOH_D>0)
 SIR_H <- k_SIR*FeOOH_H*HS_H*(HS_H>0)*(FeOOH_H>0)
 
 PyO_L <- k_PyO*FeS2_L*O2_L*(O2_L>0)*(FeS2_L>0)
 PyO_D <- k_PyO*FeS2_D*O2_D*(O2_D>0)*(FeS2_D>0)
 PyO_H <- k_PyO*FeS2_H*O2_H*(O2_H>0)*(FeS2_H>0)
 
 PyFe_L <- k_PyFe*FeOOH_L*HS_L*(HS_L>0)*(FeOOH_L>0)
 PyFe_D <- k_PyFe*FeOOH_D*HS_D*(HS_D>0)*(FeOOH_D>0)
 PyFe_H <- k_PyFe*FeOOH_H*HS_H*(HS_H>0)*(FeOOH_H>0)
 
 PyS_L <- k_PyS*Fe_L*HS_L*(HS_L>0)*(Fe_L>0)*(SO4_L>0) #*SO4_L
 PyS_D <- k_PyS*Fe_D*HS_D*(HS_D>0)*(Fe_D>0)*(SO4_D>0) #*SO4_D
 PyS_H <- k_PyS*Fe_H*HS_H*(HS_H>0)*(Fe_H>0)*(SO4_H>0) #*SO4_H
 
 
#===============================================================================
# Mass balances 
#===============================================================================
 
dCH2O_L <- F_CH2O_L/m_L - F_CH2O_LD/m_L                      - AR_L - SR_L - DIR_L - MG_L
dCH2O_H <- F_CH2O_H/m_H - F_CH2O_HD/m_H                      - AR_H - SR_H - DIR_H - MG_H
dCH2O_D <- F_CH2O_LD/m_D*dens_L/dens_D + F_CH2O_HD/m_D*dens_H/dens_D - F_CH2O_burial/m_D - AR_D - SR_D - DIR_D - MG_D      

dFeOOH_L <- F_FeOOH_L/m_L - F_FeOOH_LD/m_L                       - 4*DIR_L + FIO_L - 8*SIR_L - 2*PyFe_L + PyO_L
dFeOOH_H <- - F_FeOOH_HD/m_H                                     - 4*DIR_H + FIO_H - 8*SIR_H - 2*PyFe_H + PyO_H
dFeOOH_D <- F_FeOOH_LD/m_D*dens_L/dens_D + F_FeOOH_HD/m_D*dens_H/dens_D - F_FeOOH_burial/m_D - 4*DIR_D + FIO_D - 8*SIR_D - 2*PyFe_D + PyO_D

dFeS2_L <- - F_FeS2_LD/m_L                                   + PyFe_L + 4*PyS_L - PyO_L
dFeS2_H <- - F_FeS2_HD/m_H                                   + PyFe_H + 4*PyS_H - PyO_H
dFeS2_D <- F_FeS2_LD/m_D*dens_L/dens_D + F_FeS2_HD/m_D*dens_H/dens_D - F_FeS2_burial/m_D + PyFe_D + 4*PyS_D - PyO_D

dO2_L <- f_photo*F_CH2O_L/m_L + E_L/m_L + R_UTO2_L              - AR_L - 2*CSO_L - 1/4*FIO_L - 15/4*PyO_L - 2*AMO_L
dO2_H <- f_photo*F_CH2O_H/m_H + E_H/m_H + R_UTO2_H + R_UMO2_H   - AR_H - 2*CSO_H - 1/4*FIO_H - 15/4*PyO_H - 2*AMO_H
dO2_D <- R_UTO2_D + R_UMO2_D                                    - AR_D - 2*CSO_D - 1/4*FIO_D - 15/4*PyO_D - 2*AMO_D

dSO4_L <- F_SO4_L/m_L + R_UTSO4_L                  - 1/2*SR_L + CSO_L + SIR_L - PyS_L + 2*PyO_L - AnMO_L
dSO4_H <- R_UTSO4_H + R_UMSO4_H                    - 1/2*SR_H + CSO_H + SIR_H - PyS_H + 2*PyO_H - AnMO_H
dSO4_D <- R_UTSO4_D + R_UMSO4_D - F_SO4_gyp/m_D    - 1/2*SR_D + CSO_D + SIR_D - PyS_D + 2*PyO_D - AnMO_D

dFe_L <- R_UTFe_L            + 4*DIR_L - FIO_L + 8*SIR_L - 4*PyS_L + PyFe_L
dFe_H <- R_UTFe_H + R_UMFe_H + 4*DIR_H - FIO_H + 8*SIR_H - 4*PyS_H + PyFe_H
dFe_D <- R_UTFe_D + R_UMFe_D + 4*DIR_D - FIO_D + 8*SIR_D - 4*PyS_D + PyFe_D

dHS_L <- R_UTHS_L            + 1/2*SR_L - CSO_L - SIR_L - 7*PyS_L - 2*PyFe_L + AnMO_L
dHS_H <- R_UTHS_H + R_UMHS_H + 1/2*SR_H - CSO_H - SIR_H - 7*PyS_H - 2*PyFe_H + AnMO_H
dHS_D <- R_UTHS_D + R_UMHS_D + 1/2*SR_D - CSO_D - SIR_D - 7*PyS_D - 2*PyFe_D + AnMO_D

dCH4_L <- CH4_E_L/m_L + R_UTCH4_L              + 1./2.*MG_L - AMO_L - AnMO_L
dCH4_H <- CH4_E_H/m_H + R_UTCH4_H + R_UMCH4_H  + 1./2.*MG_H - AMO_H - AnMO_H
dCH4_D <-               R_UTCH4_D + R_UMCH4_D  + 1./2.*MG_D - AMO_D - AnMO_D                                  


#===============================================================================
# Return statement
#===============================================================================

ratesofchanges <- c(dCH2O_L  = dCH2O_L , dCH2O_H  = dCH2O_H , dCH2O_D  = dCH2O_D ,
		                dFeOOH_L = dFeOOH_L, dFeOOH_H = dFeOOH_H, dFeOOH_D = dFeOOH_D,
		                dFeS2_L  = dFeS2_L , dFeS2_H  = dFeS2_H , dFeS2_D  = dFeS2_D ,
		                dO2_L    = dO2_L   , dO2_H    = dO2_H   , dO2_D    = dO2_D   ,
		                dSO4_L   = dSO4_L  , dSO4_H   = dSO4_H  , dSO4_D   = dSO4_D  ,
		                dFe_L    = dFe_L   , dFe_H    = dFe_H   , dFe_D    = dFe_D   ,
		                dHS_L    = dHS_L   , dHS_H    = dHS_H   , dHS_D    = dHS_D   ,
		                dCH4_L   = dCH4_L  , dCH4_H   = dCH4_H  , dCH4_D   = dCH4_D   )
		  
concentrations <- c(CH2O_L  = CH2O_L , CH2O_H  = CH2O_H , CH2O_D  = CH2O_D ,
                    FeOOH_L = FeOOH_L, FeOOH_H = FeOOH_H, FeOOH_D = FeOOH_D,
                    FeS2_L  = FeS2_L , FeS2_H  = FeS2_H , FeS2_D  = FeS2_D ,
                    O2_L    = O2_L   , O2_H    = O2_H   , O2_D    = O2_D   ,
                    SO4_L   = SO4_L  , SO4_H   = SO4_H  , SO4_D   = SO4_D  ,
                    Fe_L    = Fe_L   , Fe_H    = Fe_H   , Fe_D    = Fe_D   ,
                    HS_L    = HS_L   , HS_H    = HS_H   , HS_D    = HS_D   ,
                    CH4_L   = CH4_L  , CH4_H   = CH4_H  , CH4_D   = CH4_D   )

rates  <- c(RR_CH2O_L  = F_CH2O_LD,      # [mol yr-1]  organic matter rain
            RR_CH2O_H  = F_CH2O_HD,      # [mol yr-1]  organic matter rain
            RR_FeOOH_L = F_FeOOH_LD,     # [mol yr-1]  iron oxide rain
		        RR_FeOOH_H = F_FeOOH_HD,     # [mol yr-1]  iron oxide rain
		        RR_FeS2_L  = F_FeS2_LD,       # [mol yr-1]  iron sulphide rain
		        RR_FeS2_H  = F_FeS2_HD,       # [mol yr-1]  iron sulphide rain
		        BR_CH2O    = F_CH2O_burial,  # [mol yr-1]  organic matter burial
		        BR_FeOOH   = F_FeOOH_burial, # [mol yr-1]  iron oxide burial
		        BR_FeS2    = F_FeS2_burial,   # [mol yr-1]  iron sulphide burial
		        BR_SO4     = F_SO4_gyp,   # [mol yr-1]  iron sulphide burial
		        
		        AR_L  = AR_L*m_L,  # [mol yr-1]
		        AR_H  = AR_H*m_H,
		        AR_D  = AR_D*m_D,
		        DIR_L = DIR_L*m_L,
		        DIR_H = DIR_H*m_H,
		        DIR_D = DIR_D*m_D,
		        SR_L  = SR_L*m_L,
		        SR_H  = SR_H*m_H,
		        SR_D  = SR_D*m_D,
		        MG_L  = MG_L*m_L,
		        MG_H  = MG_H*m_H,
		        MG_D  = MG_D*m_D,
		    
		        FIO_L  = FIO_L*m_L,
		        FIO_H  = FIO_H*m_H,
		        FIO_D  = FIO_D*m_D,
		        CSO_L  = CSO_L*m_L,
		        CSO_H  = CSO_H*m_H,
		        CSO_D  = CSO_D*m_D, 
		        AMO_L  = AMO_L*m_L,
		        AMO_H  = AMO_H*m_H,
		        AMO_D  = AMO_D*m_D,
		        AnMO_L = AnMO_L*m_L,
		        AnMO_H = AnMO_H*m_H,
		        AnMO_D = AnMO_D*m_D,
		        SIR_L  = SIR_L*m_L,
		        SIR_H  = SIR_H*m_H,
		        SIR_D  = SIR_D*m_D,
		        PyFe_L = PyFe_L*m_L,
		        PyFe_H = PyFe_H*m_H,
		        PyFe_D = PyFe_D*m_D,
		        PyS_L  = PyS_L*m_L,
		        PyS_H  = PyS_H*m_H,
		        PyS_D  = PyS_D*m_D,
		        PyO_L  = PyO_L*m_L,
		        PyO_H  = PyO_H*m_H,
		        PyO_D  = PyO_D*m_D
		        )
		  
		return(list(ratesofchanges=ratesofchanges, other=list(conc=concentrations,rate=rates)))
 }
 )
}

#==============================================================================
# Initial conditions 
#==============================================================================

SV <- list(
  
  CH2Oi_L  = 2.251251e+00, # mmol/kg-soln    initial CH2O value box L 
  CH2Oi_H  = 6.571390e-01, # mmol/kg-soln    initial CH2O value box H 
  CH2Oi_D  = 1.855115e-03, # mmol/kg-soln    initial CH2O value box D 
  
  FeOOHi_L = 2.499875e-02, # mmol/kg-soln    initial FeOOH value box L 
  FeOOHi_H = 1.404909e-07, # mmol/kg-soln    initial FeOOH value box H 
  FeOOHi_D = 1.681526e-02, # mmol/kg-soln    initial FeOOH value box D 
  
  FeS2i_L  = 3.188341e-16, # mmol/kg-soln    initial total FeS concentration 
  FeS2i_D  = 1.130908e-21, # mmol/kg-soln    initial total FeS concentration 
  FeS2i_H  = 1.561664e-21, # mmol/kg-soln    initial total FeS concentration 
  
  O2i_L    = 3.001201e-01, # mmol/kg-soln    initial total O2 concentration (fully oxygenated)
  O2i_D    = 2.303450e-01, # mmol/kg-soln    initial total O2 concentration (fully oxygenated)
  O2i_H    = 2.993231e-01, # mmol/kg-soln    initial total O2 concentration (fully oxygenated)
  
  SO4i_L   = 2.609580e+01, # mmol/kg-soln    initial total SO4 concentration (present day)
  SO4i_D   = 2.600000e+01, # mmol/kg-soln    initial total SO4 concentration (present day)
  SO4i_H   = 2.600185e+01, # mmol/kg-soln    initial total SO4 concentration (present day)
  
  Fei_L    = 1.594552e-09, # mmol/kg-soln    initial total Fe concentration (fully oxygenated)
  Fei_D    = 1.614070e-12, # mmol/kg-soln    initial total Fe concentration (fully oxygenated)
  Fei_H    = 2.646287e-14, # mmol/kg-soln    initial total Fe concentration (fully oxygenated)
  
  HSi_L    = 6.000607e-08, # mmol/kg-soln    initial total HS concentration (fully oxygenated)
  HSi_D    = 8.393194e-11, # mmol/kg-soln    initial total HS concentration (fully oxygenated)
  HSi_H    = 1.766437e-08, # mmol/kg-soln    initial total HS concentration (fully oxygenated)
  
  CH4i_L   = 2.342227e-09, # mmol/kg-soln    initial total HS concentration (fully oxygenated)
  CH4i_D   = 6.944137e-10, # mmol/kg-soln    initial total HS concentration (fully oxygenated)
  CH4i_H   = 2.543184e-12  # mmol/kg-soln    initial total HS concentration (fully oxygenated)
  
)

with (c(PL,SV),{
			
   		SV <<- c(CH2O_L  = CH2Oi_L , CH2O_H  = CH2Oi_H , CH2O_D  = CH2Oi_D ,
					     FeOOH_L = FeOOHi_L, FeOOH_H = FeOOHi_H, FeOOH_D = FeOOHi_D,
					     FeS2_L  = FeS2i_L , FeS2_H  = FeS2i_H , FeS2_D  = FeS2i_D ,
					     O2_L    = O2i_L   , O2_H    = O2i_H   , O2_D    = O2i_D   ,
					     SO4_L   = SO4i_L  , SO4_H   = SO4i_H  , SO4_D   = SO4i_D  ,
					     Fe_L    = Fei_L   , Fe_H    = Fei_H   , Fe_D    = Fei_D   ,
					     HS_L    = HSi_L   , HS_H    = HSi_H   , HS_D    = HSi_D   ,
					     CH4_L   = CH4i_L  , CH4_H   = CH4i_H  , CH4_D   = CH4i_D   )    
		})

#==============================================================================
# Test trial of the boxmodel to see if it works 
#==============================================================================

PL$v_sink_mineral <- 63.79

#output <- runsteady(y=SV,func=boxmodel.CSOFe,parms=PL)
output <- steady.1D(y=SV,func=boxmodel.CSOFe,parms=PL,nspec=N.vars)
