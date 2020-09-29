###############################################################################
# Model of early diagenesis
# Type: Carbon,sulfur and iron cycling  
# Author: Filip Meysman (filip.meysman@nioz.nl)
#         Sebastiaan van de Velde (sevdevel@vub.ac.be) 
# Affiliation: NIOZ, Korringaweg 7, 4401 NT, Yerseke 
#              VUB, Pleinlaan 2, 1050 BE, Brussel
###############################################################################

# public packages
require(ReacTran)
require(marelac)
#require(AquaEnv)

# non-public packages
require(diagenesis)

#=============================================================================
# Definition of specific functions
#=============================================================================

#-----------------------------------------------------------------------------
# Saturation function
#-----------------------------------------------------------------------------

FSAT <- function(C,K,n) (C/K)^n/((C/K)^n + 1)

#-----------------------------------------------------------------------------
# Function: tortuosity 
#-----------------------------------------------------------------------------

tortuosity <- function(por) 1-2*log(por)

#-----------------------------------------------------------------------------
# Function: init.diffusion.coefficient 
#-----------------------------------------------------------------------------

init.diffusion.coefficient <- function (species,S,TC,P,conv.fac,grid,tort.grid)
{
  Dmol <- conv.fac*diffcoeff(S=S,t=TC,P=P,species=species)[[species]]
  D.grid <- setup.prop.1D(value=Dmol,grid=grid)
  D.grid$mid <- D.grid$mid/tort.grid$mid
  D.grid$int <- D.grid$int/tort.grid$int
  return(D.grid)
}

#=============================================================================
# Definition of objects: species, reactions, elements, simulations
#=============================================================================

# SL = list of species incorporated in the model 
SL <- list()

# chemical species incorporated in the model 

SL <- initialise.list.species(SL,"1D",name="OC",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="FeOOH",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="FeS",phase="solid",type="chemical")

SL <- initialise.list.species(SL,"1D",name="O2",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="SO4",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="HCO3",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="H2O",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="H",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="HS",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Fe",phase="solute",type="chemical")


# Composite species incorporated in the model 

#SL <- initialise.list.species(SL,"1D",name="OC",phase="solid",type="composite")

# Artificial species incorporated in the model 

#SL <- initialise.list.species(SL,"1D",name="Omega.FeS",phase="solute",type="artificial")

# Initialise species matrix 

spec.mat <- initialise.species.matrix(SL)
#spec.mat["OC",c("OC.f","OC.s")] <- c(1,1)

spec.table <- as.data.frame(spec.mat)
remove(spec.mat)

# List of reactions included 

RL <- list()

RL <- initialise.list.reaction(RL,"1D",name="Cmin",type="kinetic")

RL <- initialise.list.reaction(RL,"1D",name="AR",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="DIR",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SR",type="kinetic")

RL <- initialise.list.reaction(RL,"1D",name="CSO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="FIO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SIR",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="ISP",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="ISO",type="kinetic")


# List of elements included 

EL <- list() # create empty element object list

EL <- initialise.list.element(EL,"1D",name="C")
EL <- initialise.list.element(EL,"1D",name="O")
EL <- initialise.list.element(EL,"1D",name="H")
EL <- initialise.list.element(EL,"1D",name="S")
EL <- initialise.list.element(EL,"1D",name="Fe")
EL <- initialise.list.element(EL,"1D",name="e")

# Initialise element matrix 

elt.mat <- initialise.element.matrix(SL, EL)

elt.mat["OC",c("C","H","O")] <- c(1,2,1)

elt.mat["O2",c("O")] <- c(2)
elt.mat["HCO3",c("H","C","O","e")] <- c(1,1,3,1)
elt.mat["H",c("H","e")] <- c(1,-1)
elt.mat["H2O",c("H","O")] <- c(2,1)

elt.mat["FeS",c("S","Fe")] <- c(1,1)
elt.mat["FeOOH",c("Fe","O","H")] <- c(1,2,1)
elt.mat["Fe",c("Fe","e")] <- c(1,-2)

elt.mat["HS",c("H","S","e")] <- c(1,1,1)
elt.mat["SO4",c("S","O","e")] <- c(1,4,2)


elt.table <- as.data.frame(elt.mat)
remove(elt.mat)

################################################################################
# initialise.parameters function
# This function initialises a default set of parameters (baseline simulation),
# but this function is also used when updating the parameter set.
################################################################################

initialise.parameters <- function(PL=NULL)
  
  # Model parameters are intialised on the first call (PL = NULL), and updated 
  # when a parmeter list is supplied (PL not NULL)
  
  # begin initialise.parameters
{
  
  #=============================================================================
  # Physical units 
  #=============================================================================
  
  if (is.null(PL$units$M)) PL$units$M <- "umol" # Mass unit
  if (is.null(PL$units$T)) PL$units$T <- "yr"   # Time unit
  if (is.null(PL$units$L)) PL$units$L <- "cm"   # Length unit
  
  # Limiting concentration: numerical limit on consumption in kinetic reactions
  PL$C.lim <- 1E-3
  
  #=============================================================================
  # Model domain and grid definition
  #=============================================================================
  
  if (is.null (PL$L)) PL$L <- 30   # depth of sediment domain [cm]
  if (is.null (PL$N)) PL$N <- 200  # number of grid layers
  
  # Call "setup.grid.1D" from ReacTran
  if (is.null (PL$grid.type)) PL$grid.type <- "uniform"  # number of grid layers
  if (is.null (PL$dx.1)) PL$dx.1 <- PL$L/(10*PL$N)
  #if (PL$grid.type == "uniform")
  #{
    # even grid - all depth layers have same thickness 
    #PL$grid <- setup.grid.1D(x.up=0,L=PL$L,N=PL$N)
  #} else {
    # uneven grid - higher resolution near the sediment water interface
    PL$grid <- setup.grid.1D(x.up=0,L=PL$L,N=PL$N,dx.1=PL$dx.1)
  #}
  
  #=============================================================================
  # Sediment parameters: 
  #=============================================================================
  
  # Environmental parameters
  
  if (is.null(PL$S))  PL$S  <- 35    # salinity
  if (is.null(PL$TC)) PL$TC <- 10    # temperature [deg C]
  if (is.null(PL$P))  PL$P  <- 1.013 # pressure [bar]
  
  # Porosity profile 
  
  if (is.null(PL$por.0)) PL$por.0 <- 0.8     # porosity at the sediment-water interface
  if (is.null(PL$por.inf)) PL$por.inf  <- 0.8   # asymptotic porosity at depth
  if (is.null(PL$por.x.att)) PL$por.x.att <- 4    # attenuation depth [cm]

  
  PL$por.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$por.0,y.inf=PL$por.inf,x.L=0,x.att=PL$por.x.att)
  PL$svf.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=(1-PL$por.0),y.inf=(1-PL$por.inf),x.L=0,x.att=PL$por.x.att)
  
  # Initialisation tortuosity 

  PL$tort.grid <- setup.prop.1D(value=0,grid=PL$grid)
  PL$tort.grid$mid <- tortuosity(PL$por.grid$mid)
  PL$tort.grid$int <- tortuosity(PL$por.grid$int)

  # Fixed profiles: pH
  
  PL$pH      <- 7.5
  PL$pH.grid <-  setup.prop.1D(value=PL$pH,grid=PL$grid)
  
  #=============================================================================
  # Diffusion coefficients 
  # Uses routine 'diffcoeff' from 'marelac' package
  #=============================================================================
  
  c.fac <- 10000*(3600*24*365.25) # conversion from [m2 s-1] to [cm2 yr-1]
  
  # Diffusion coefficients from marelac package
  PL$D.O2.grid   <- init.diffusion.coefficient(species="O2",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.SO4.grid  <- init.diffusion.coefficient(species="SO4",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.HCO3.grid <- init.diffusion.coefficient(species="HCO3",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.H2O.grid  <- init.diffusion.coefficient(species="H2O",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.H.grid    <- init.diffusion.coefficient(species="H",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.HS.grid   <- init.diffusion.coefficient(species="HS",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.Fe.grid   <- init.diffusion.coefficient(species="Fe",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  
  #=============================================================================
  # Transport parameters
  #=============================================================================
  
  # Diffusive boundary layer
  
  if (is.null (PL$x.DBL)) PL$x.DBL <- 0 # thickness of DBL [cm]
  
  # Advective velocities 
  
  if (is.null(PL$rho.sed)) PL$rho.sed <- 2.6  # density solid sediment [g cm-3]
  if (is.null(PL$u.0)) PL$u.0 <- 0.2  # sedimentation velocity pore water [cm yr-1]
  if (is.null(PL$v.0)) PL$v.0 <- 0.2  # sedimentation velocity solids [cm yr-1]
  
  PL$SedFlux <- PL$rho.sed*(1-PL$por.0)*PL$v.0  # sedimentation flux [g cm-2 yr-1]
  
  PL$v.grid <- setup.prop.1D(value=PL$v.0,grid=PL$grid) # advection velocity solid
  PL$u.grid <- setup.prop.1D(value=PL$u.0,grid=PL$grid) # advection velocity pore water
  
  # Bioturbation profile
  
  if (is.null(PL$x.L))  PL$x.L    <- 5
  if (is.null(PL$x.att))PL$x.att  <- 5
  if (is.null(PL$y.inf))PL$y.inf  <- 0
  if (is.null(PL$Db.0)) PL$Db.0   <- 10     # biodiffusion coefficient Db at SWI [cm2 yr-1]
  PL$Db.grid <- setup.prop.1D(func=p.sig,y.0=PL$Db.0,y.inf=PL$y.inf,x.att=PL$x.att,x.L=PL$x.L, grid = PL$grid)
  
  # Irrigation profile
  
  if (is.null(PL$irr.0)) PL$irr.0 <- 0     # irrigation rate at SWI [yr-1]
  PL$irr.grid <- setup.prop.1D(func=p.exp,y.0=PL$irr.0,y.inf=0,x.L=0,x.att=3,grid=PL$grid)
  
  #=============================================================================
  # Reaction parameters
  #=============================================================================
  
  # Decay constants organic matter
  
  if (is.null(PL$k.min)) PL$k.min    <- 0.1  # decay constant of organic matter [yr-1]
  
  if (is.null(PL$K_O2))   PL$K_O2    <- 0.001      # Monod constant O2 consumption [umol cm-3 or mM] (Meysman et al., 2003)
  if (is.null(PL$K_FeOOH))PL$K_FeOOH <- 4*PL$rho.sed  # Monod constant FeIII reduction [umol cm-3 or mM] (Meysman et al., 2003)
  if (is.null(PL$K_SO4))  PL$K_SO4   <- 0.9     # Monod constant SO4 reduction [umol cm-3 or mM] (Meysman et al., 2003)
  
  # Secondary reactions  
  
  if (is.null(PL$k.CSO))  PL$k.CSO  <- 1E+07   # kinetic constant sulfide oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996; Meysman et al., 2003)
  if (is.null(PL$k.FIO))  PL$k.FIO  <- 1E+07  # kinetic constant Ferrous iron oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996; Meysman et al., 2003)
  if (is.null(PL$k.CSFO)) PL$k.SIR  <- 494     # kinetic constant sulfide oxidation by FeIII [umol-1 cm3 yr-1] (Berg et al., 2003)
  if (is.null(PL$k.ISP))  PL$k.ISP  <- 1E+04   # kinetic constant FeS precipitation [umol-1 cm3 yr-1] (Meysman et al., 2015)
  if (is.null(PL$k.ISO))  PL$k.ISO  <- 1E+07    # kinetic constant iron sulfide oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996)
     
  
  #=============================================================================
  # Boundary conditions
  #=============================================================================
  
  # Upper boundary fluxes

  F.OC <- 15 # mmol m-2 d-1
  if (is.null(PL$F.OC)) PL$F.OC <- F.OC*365.25/10 # conversion to umol cm-2 yr-1 

  F.FeS <- 0 # mmol m-2 d-1
  if (is.null(PL$F.FeS)) PL$F.FeS <- F.FeS*365.25/10 # conversion to umol cm-2 yr-1 
  
  # 1% Fe content of mineral deposition 
  F.FeOOH <- 1E06*0.01*PL$SedFlux/(55.847) # mmol m-2 d-1
  if (is.null(PL$F.FeOOH)) PL$F.FeOOH <- F.FeOOH*365.25/10 # conversion to umol cm-2 yr-1 

  # Upper boundary concentrations

  if (is.null(PL$O2.ow))   PL$O2.ow   <- 0.28 # O2 concentration bottom water [umol cm-3 or mM]
  if (is.null(PL$HCO3.ow)) PL$HCO3.ow <- 2.2 # SumCO2 concentration bottom water [umol cm-3 or mM]
  if (is.null(PL$HS.ow))   PL$HS.ow   <- 0   # SumH2S concentration bottom water [umol cm-3 or mM]
  if (is.null(PL$H2O.ow))  PL$H2O.ow  <- 0
  if (is.null(PL$SO4.ow))  PL$SO4.ow  <- 140.0*PL$S/((32.065+4*15.999)*1.80655)
  if (is.null(PL$Fe.ow))   PL$Fe.ow   <- 0
  if (is.null(PL$H.ow))    PL$H.ow    <- 0 # Proton concentration bottom water, without equilibrium with H2O
    
  # Boundary condition lower boundary 
  # All species: no gradient
  
  if (is.null(PL$O2.ds))   PL$O2.ds   <- NA
  if (is.null(PL$HCO3.ds)) PL$HCO3.ds <- NA
  if (is.null(PL$HS.ds))   PL$HS.ds   <- NA
  if (is.null(PL$H2O.ds))  PL$H2O.ds  <- NA
  if (is.null(PL$H.ds))    PL$H.ds    <- NA
  if (is.null(PL$SO4.ds))  PL$SO4.ds  <- NA
  if (is.null(PL$Fe.ds))   PL$Fe.ds   <- NA
  
   
  #=============================================================================
  # Flags 
  #=============================================================================
  
  if (is.null(PL$simulation.type)) PL$simulation.type <- "direct.steady.state"
  #simulation.type <- "dynamic.steady.state"
  #simulation.type <- "time.dependent"
  
  ##############################################################################
  return(PL)
  # end initialise.parameters
}
##############################################################################

#=============================================================================
# Model formulation
#=============================================================================

CSFe.model <- function (t,state,parameters,summary.call=FALSE) 
{
  
  with(as.list(c(parameters)),{
    
    #---------------------------------------------------------------------------
    # Initialisation of 
    # SV : matrix of state variables
    # OL : Object List
    #---------------------------------------------------------------------------
    
    SV <- matrix(nrow=N,ncol=N.var,data=state)
    for (i in 1:N.var) OL[[var.names[i]]]$C <- SV[,i]
    
    #---------------------------------------------------------------------------
    OL <- within(OL,{
    #---------------------------------------------------------------------------
      
    #=========================================================================
    # Initialisation of concentration depth profiles
    #=========================================================================
    #-------------------------------------------------------------------------
    # Initialisation of species with fixed concentration depth profiles 
    #-------------------------------------------------------------------------
    
    if (!("OC" %in% var.names)) OC$C <- OC.fix
    
    #-------------------------------------------------------------------------
    # Upper and lower boundary conditions
    #-------------------------------------------------------------------------

    O2$C.up   <- O2.ow
    HCO3$C.up <- HCO3.ow
    HS$C.up   <- HS.ow
    SO4$C.up  <- SO4.ow
    H2O$C.up  <- H2O.ow
    H$C.up    <- H.ow
    Fe$C.up   <- Fe.ow
    
    if (is.na(O2.ds))   O2$C.down   <- O2$C[N]   else O2$C.down   <- O2.ds 
    if (is.na(HCO3.ds)) HCO3$C.down <- HCO3$C[N] else HCO3$C.down <- HCO3.ds 
    if (is.na(HS.ds))   HS$C.down   <- HS$C[N]   else HS$C.down   <- HS.ds 
    if (is.na(SO4.ds))  SO4$C.down  <- SO4$C[N]  else SO4$C.down  <- SO4.ds 
    if (is.na(H2O.ds))  H2O$C.down  <- H2O$C[N]  else H2O$C.down  <- H2O.ds 
    if (is.na(H.ds))    H$C.down    <- H$C[N]    else H$C.down    <- H.ds 
    if (is.na(Fe.ds))   Fe$C.down   <- Fe$C[N]   else Fe$C.down   <- Fe.ds 
    
    
    #=========================================================================
    # TRANSPORT terms using tran.1D from ReacTran
    #=========================================================================
    
    if (summary.call) {ind <- c("tran","C.up","C.down","dif.flux","adv.flux","flux","flux.up","flux.down")} else {ind <- c("tran","flux.up","flux.down")}
    
    # Transport terms: Solids [umol cm-3 solid yr-1]
    
    OC[ind] <- tran.1D(C=OC$C,flux.up=F.OC,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    FeOOH[ind]<- tran.1D(C=FeOOH$C,flux.up=F.FeOOH,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    FeS[ind]  <- tran.1D(C=FeS$C,flux.up=F.FeS,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    
    # Transport terms: Solutes [umol cm-3 pore water yr-1]
    
    O2[ind]   <- tran.1D(C=O2$C,C.up=O2$C.up,C.down=O2$C.down,v=u.grid,D=D.O2.grid,AFDW=0.5,VF=por.grid,dx=grid,full.output=summary.call)
    HCO3[ind] <- tran.1D(C=HCO3$C,C.up=HCO3$C.up,C.down=HCO3$C.down,v=u.grid,D=D.HCO3.grid,AFDW=0.5,VF=por.grid,dx=grid,full.output=summary.call)
    H2O[ind]  <- tran.1D(C=H2O$C,C.up=H2O$C.up,C.down=H2O$C.down,v=u.grid,D=D.H2O.grid,AFDW=0.5,VF=por.grid,dx=grid,full.output=summary.call)
    H[ind]    <- tran.1D(C=H$C,C.up=H$C.up,C.down=H$C.down,v=u.grid,D=D.H.grid,AFDW=0.5,VF=por.grid,dx=grid,full.output=summary.call)
    HS[ind]   <- tran.1D(C=HS$C,C.up=HS$C.up,C.down=HS$C.down,v=u.grid,D=D.HS.grid,AFDW=0.5,VF=por.grid,dx=grid,full.output=summary.call)
    SO4[ind]  <- tran.1D(C=SO4$C,C.up=SO4$C.up,C.down=SO4$C.down,v=u.grid,D=D.SO4.grid,AFDW=0.5,VF=por.grid,dx=grid,full.output=summary.call)
    Fe[ind]   <- tran.1D(C=Fe$C,C.up=Fe$C.up,C.down=Fe$C.down,v=u.grid,D=D.Fe.grid,AFDW=0.5,VF=por.grid,dx=grid,full.output=summary.call)
        
    #-----------------------------------------------------------------------------
    # Saturation state calculation
    #-----------------------------------------------------------------------------
    
    # Saturation state FeS 
    
    #Omega.FeS$C <- Fe$C*HS$C/(K_IS)
    #Omega.FeS$C.up <- Fe.ow*HS$C.up/K_IS
    #Omega.FeS$C.down <- Fe$C.down*HS$C.down/K_IS
    
    #=========================================================================
    # Rates of Kinetic reactions 
    # Expressed per volume of bulk sediment [umol cm-3 yr-1]
    #=========================================================================
    
    # Mineralisation reactions
    
    Cmin$R <- svf.grid$mid*k.min*OC$C*FSAT(OC$C,C.lim,5)
    
    FeOOH.lim <- FeOOH$C/(FeOOH$C+K_FeOOH)*FSAT(FeOOH$C,C.lim,5)
    FeOOH.inh <- K_FeOOH/(FeOOH$C+K_FeOOH)
    SO4.lim   <- SO4$C/(SO4$C+K_SO4)*FSAT(SO4$C,C.lim,5)
    
    f.f <- FeOOH.lim/(FeOOH.lim + SO4.lim*FeOOH.inh)
    s.f <- (SO4.lim*FeOOH.inh)/(FeOOH.lim + SO4.lim*FeOOH.inh)
    
#     AR.f$R <- O2.lim*Cmin.f$R  # Aerobic respiration
#     FR.f$R <- FeOOH.lim*O2.inh*Cmin.f$R # FeIII reduction
#     SR.f$R <- SO4.lim*O2.inh*Cmin.f$R # Sulfate reduction
#     
#     AR.s$R <- O2.lim*Cmin.s$R  # Aerobic respiration
#     FR.s$R <- FeOOH.lim*O2.inh*Cmin.s$R # FeIII reduction
#     SR.s$R <- SO4.lim*O2.inh*Cmin.s$R # Sulfate reduction
    
    DIR$R <- f.f*Cmin$R # FeIII reduction
    SR$R  <- s.f*Cmin$R # Sulfate reduction
    
    Cmin$R <- DIR$R + SR$R
    
    # Re-oxidation reactions
    
    CSO$R  <- por.grid$mid*k.CSO*O2$C*FSAT(O2$C,C.lim/100,5)*HS$C*FSAT(HS$C,C.lim/100,5)        # Sulfide re-oxidation
    FIO$R  <- por.grid$mid*k.FIO*O2$C*FSAT(O2$C,C.lim,5)*Fe$C*FSAT(Fe$C,C.lim/100,5)            # iron re-oxidation
    ISP$R  <- por.grid$mid*k.ISP*Fe$C*FSAT(Fe$C,C.lim,5)*HS$C*FSAT(HS$C,C.lim,5)                # iron sulfide precipitation
    
    #=========================================================================
    # Consumption terms: kinetic reactions
    # Solutes [umol cm-3 pore water yr-1] Solids [umol cm-3 solid yr-1] 
    #=========================================================================
    
    OC$kin.reac    <- - (DIR$R + SR$R)/svf.grid$mid
    FeOOH$kin.reac <- (- 4*DIR$R + 4*FIO$R)/svf.grid$mid
    FeS$kin.reac   <- (ISP$R)/svf.grid$mid
        
    O2$kin.reac  <- (- FIO$R - 2*CSO$R)/por.grid$mid
    Fe$kin.reac  <- (4*DIR$R - 4*FIO$R  - ISP$R)/por.grid$mid
    SO4$kin.reac <- (- 0.5*SR$R + CSO$R)/por.grid$mid
    HS$kin.reac  <- (0.5*SR$R - CSO$R - ISP$R)/por.grid$mid
    HCO3$kin.reac<- (DIR$R + SR$R)/por.grid$mid
    H2O$kin.reac <- (6*DIR$R - 6*FIO$R)/por.grid$mid
    H$kin.reac   <- (7*DIR$R + 1/2*SR$R + 8*FIO$R + CSO$R + ISP$R)/por.grid$mid
    
    #=========================================================================
    # Equilibrium reactions
    # Solutes [umol cm-3 pore water yr-1] Solids [umol cm-3 solid yr-1] 
    #=========================================================================

    #=========================================================================
    # Irrigation 
    # Solutes [umol cm-3 pore water yr-1] Solids [umol cm-3 solid yr-1] 
    #=========================================================================

    Fe$irr   <- irr.grid$mid*(Fe.ow - Fe$C)
    O2$irr   <- irr.grid$mid*(O2.ow - O2$C)
    HS$irr   <- irr.grid$mid*(HS.ow - HS$C)
    HCO3$irr <- irr.grid$mid*(HCO3.ow - HCO3$C)
    SO4$irr  <- irr.grid$mid*(SO4.ow - SO4$C)
    
    #=========================================================================
    # Rate of changes composite species 
    #=========================================================================
    
    #=========================================================================
    # Rate of changes of all species
    #=========================================================================
    
  # Fe$irr <- alfa*(PL$Fe.ow - Fe$C)

    OC$ddt    <- OC$tran + OC$kin.reac + OC$eq.reac + OC$irr
    FeOOH$ddt <- FeOOH$tran + FeOOH$kin.reac + FeOOH$eq.reac + FeOOH$irr
    FeS$ddt   <- FeS$tran + FeS$kin.reac + FeS$eq.reac + FeS$irr
  
    O2$ddt    <- O2$tran + O2$kin.reac + O2$eq.reac + O2$irr
    Fe$ddt    <- Fe$tran + Fe$kin.reac + Fe$eq.reac + Fe$irr
    HCO3$ddt  <- HCO3$tran + HCO3$kin.reac + HCO3$eq.reac + HCO3$irr
    H2O$ddt   <- H2O$tran + H2O$kin.reac + H2O$eq.reac + H2O$irr
    H$ddt     <- H$tran + H$kin.reac + H$eq.reac + H$irr
    HS$ddt    <- HS$tran + HS$kin.reac + HS$eq.reac + HS$irr
    SO4$ddt   <- SO4$tran + SO4$kin.reac + SO4$eq.reac + SO4$irr
    
    #=========================================================================
    }) # end within call
  #=========================================================================
  
  #-----------------------------------------------------------------------------
  # Assemble matrix with rate of changes (RC) 
  #-----------------------------------------------------------------------------
  
  RC <- 0*SV
  for (i in 1:N.var) RC[,i] <- OL[[var.names[i]]]$ddt
  
  #-----------------------------------------------------------------------------
  # Monitoring of simulation progession during dynamic simulation 
  #-----------------------------------------------------------------------------
  
  if ((simulation.type %in% c("dynamic.steady.state","time.dependent")) && (t > sim.info$tc)){
    
    # sim.info$fc = keeps track of number of function calls
    # sim.info$tc = keeps track of incremental simulation time
    
    sim.info$fc <<- sim.info$fc + 1  
    sim.info$tc <<- t*1
    
    # Calculate for each state variable 
    #  ACC = absolute concentration change 
    #  RCC = relative concentration change 
    
    SV.ref <- rep(1E-03,N.var)
    ACC.t <- ACC(RC)
    RCC.t <- RCC(RC,SV,SV.ref)
    
    # Output of simulation progress to the screen 
    screenprint.simulation.progress(sim.info$fc,t,SV,ACC.t,RCC.t,var.names)
    
    # Store simulation progress results 
    sim.info$progress <<- rbind(sim.info$progress,c(t,ACC.t,RCC.t))
    points(log10(t),log10(sqrt(sum(RCC.t^2)/length(RCC.t))),pch=16,col="red")
    
  }
  
  #===========================================================================
  # return statement 
  #===========================================================================
  
  if (!summary.call) {return(list(RC))}
  else {
    
    # Initialisation of lists
    
    SL <- OL[names(SL)]
    RL <- OL[names(RL)]
    EL <- EL
    
    # Updating species information
    
    SL <- updating.chemical.species(SL,PL) 
    SL <- updating.composite.species(spec.table,SL) 
    
    # Updating reaction information
    
    for (i in 1:length(RL)) RL[[i]]$R.int <- sum(RL[[i]]$R*grid$dx) 
    
    # Updating element information
    
    EL <- updating.elements(elt.table,EL,SL) 
    
    return(list(SL=SL,RL=RL,EL=EL,PL=parameters))
    
   }
   
   #-----------------------------------------------------------------------------
 }) # end with call
}
#-----------------------------------------------------------------------------

################################################################################
# Initialisation of state variable matrix
################################################################################

initialise.state.variables <- function(PL)
{
  
  PL$var.names <- c("OC","FeOOH","FeS","O2","Fe","SO4","HCO3","HS","H2O","H")
  PL$non.var.names <- NULL
  
  PL$N <- PL$grid$N
  
  PL$N.var <- length(PL$var.names)
  PL$N.non.var <- length(PL$non.var.names)
  
  return(PL)
}

initialise.species.arrays <- function(PL,SL)
{
  
  for (i in 1:length(SL))
  {
    SL[[i]]$C        <- rep(0,PL$N)
    SL[[i]]$tran     <- rep(0,PL$N)
    SL[[i]]$kin.reac <- rep(0,PL$N)
    SL[[i]]$eq.reac  <- rep(0,PL$N)
    SL[[i]]$irr      <- rep(0,PL$N)
  }
  return(SL)
}

# Function that creates a rough guess for the initial state

crude.initial.state <- function(PL)
{
  with(PL,{
    
    # Initialisation State Variable matrix (SV) 
    
    SV.guess <- matrix(nrow=N,ncol=N.var,dimnames=list(1:N,var.names))
    
    if ("OC" %in% PL$var.names) SV.guess[,"OC"]        <- rep(0,length.out=N)
    
    if ("FeOOH" %in% PL$var.names)  SV.guess[,"FeOOH"] <- rep(0,length.out=N)
    if ("FeS" %in% PL$var.names)  SV.guess[,"FeS"]     <- rep(0,length.out=N)
    if ("Fe" %in% PL$var.names) SV.guess[,"Fe"]        <- rep(Fe.ow,length.out=N)
    
    if ("O2" %in% PL$var.names) SV.guess[,"O2"]        <- rep(O2.ow,length.out=N)
    if ("HCO3" %in% PL$var.names) SV.guess[,"HCO3"]    <- rep(HCO3.ow,length.out=N)
    if ("HS" %in% PL$var.names) SV.guess[,"HS"]        <- rep(HS.ow,length.out=N)
    if ("SO4" %in% PL$var.names) SV.guess[,"SO4"]      <- rep(SO4.ow,length.out=N)
    if ("H2O" %in% PL$var.names) SV.guess[,"H2O"]      <- rep(H2O.ow,length.out=N)
    if ("H" %in% PL$var.names) SV.guess[,"H"]          <- rep(H.ow,length.out=N)
    
    return(SV.guess)
    
  })
}


################################################################################
# Test trial model runs
################################################################################

model <- CSFe.model
sim.info <- list(tc= 0,fc=0)

# Initialisation of object list (OL)  

PL <- initialise.parameters(PL=list())
PL$simulation.type <- "direct.steady.state"
SL <- initialise.species.arrays(PL,SL) # default action
OL <- c(SL,RL) # default action

# Test trial 

PL <- initialise.state.variables(PL)
SV <- crude.initial.state(PL)

system.time({
  result <- model(t=0,state=SV,parameters=PL,summary.call=FALSE)
  result <- model(t=0,state=SV,parameters=PL,summary.call=TRUE)
})

remove(result)
remove(SV)

