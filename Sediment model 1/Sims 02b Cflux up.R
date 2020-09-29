# CSFe - Fe isotope model 
# Author: Sebastiaan van de Velde
###############################################################################
is.local <- TRUE
FOLDER   <- "Sensitivityrun"
BASENAME <- "02b 05 Cflux"
#=============================================================================
# Compile packages for current node
if (!is.local){
  
  dir.create(BASENAME)
  
  install.packages("../packages/gsw_1.0-5.tar.gz", repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/oce_1.2-0.tar.gz", repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/rootSolve_1.8.2.tar.gz",repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/deSolve_1.27.1.tar.gz",repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/shape_1.4.4.tar.gz",repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/ReacTran_1.4.3.1.tar.gz",repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/AquaEnv_1.0-4.tar.gz",repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/seacarb_3.2.12.tar.gz",repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/marelac_2.1.10.tar.gz",repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  
  # public packages
  require(AquaEnv,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(rootSolve,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(deSolve,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(shape,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(ReacTran,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(gsw,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(oce,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(seacarb,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(marelac,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  
  # non-public packages
  require(diagenesis,lib.loc=getwd())
}
#=============================================================================
# Source file containing model function
if (!is.local){source("../Model_CSFe_Feiso_HPCversion_v03.R")}
if (is.local){source("Model_CSFe_Feiso_localversion_v04.R")}

#=============================================================================
# parameters to be tested
F.POC.sens   <- 1000.
# Number of simulations
model <- CSFe.model
#=============================================================================
# Dynamic simulation 1 (start from previous initial conditions)
# add FeOOH
#=============================================================================
load("02b 04 Cflux SteadySta 1.Rdata")
# Initialisation simulation type 
sim.info$index <- 1
sim.info$code <- paste(BASENAME,"Transient")
sim.info$name <- paste(sim.info$code,sim.info$index)
# Initialisation parameter list
sim.start <- simulation.list[[length(simulation.list)]]
PL <- sim.start$PL
# Adapt parameter list
PL$simulation.type <- "time.dependent"
#PL$simulation.type <- "direct.steady.state"
PL$F.POC    <- F.POC.sens
PL$F.POC.1  <- PL$F.POC*0.021746
PL$F.POC.2  <- PL$F.POC*0.00725275
PL$F.POC.3  <- PL$F.POC*0.0096717
PL$F.POC.4  <- PL$F.POC*0.0128974
PL$F.POC.5  <- PL$F.POC*0.017199
PL$F.POC.6  <- PL$F.POC*0.0229352
PL$F.POC.7  <- PL$F.POC*0.0305846
PL$F.POC.8  <- PL$F.POC*0.0407852
PL$F.POC.9  <- PL$F.POC*0.0543879
PL$F.POC.10 <- PL$F.POC*0.0725265
PL$F.POC.11 <- PL$F.POC*0.0967046
PL$F.POC.12 <- PL$F.POC*0.12881
PL$F.POC.13 <- PL$F.POC*0.169822
PL$F.POC.14 <- PL$F.POC*0.314677
PL <- initialise.parameters(PL)
# Initialise set of state variables
PL <- initialise.state.variables(PL)
# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)
# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)
# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)
# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)
#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------
sim.info$time.seq <- c(0,1,6,12,24,1*24*365.25,1000*24*365.25,10000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","1 yr","1000 yr","10000 yr")
sim.info$N.out <- length(sim.info$time.seq)
#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------
# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix
# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV.init, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)
# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)
#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------
# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation
save(sim.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Dynamic simulation 1 (start from previous initial conditions)
# add FeOOH
#=============================================================================
load(paste(BASENAME,"Transient","1.Rdata"))
# Initialisation simulation type 
sim.info$index <- 1
sim.info$code <- paste(BASENAME,"Transient")
sim.info$name <- paste(sim.info$code,sim.info$index)
# Initialisation parameter list
sim.start <- simulation.list[[length(simulation.list)]]
PL <- sim.start$PL
# Adapt parameter list
PL$simulation.type <- "time.dependent"
#PL$simulation.type <- "direct.steady.state"
PL <- initialise.parameters(PL)
# Initialise set of state variables
PL <- initialise.state.variables(PL)
# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)
# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)
# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)
# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)
#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------
sim.info$time.seq <- c(0,1,6,12,24,1*24*365.25,1000*24*365.25,10000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","1 yr","1000 yr","10000 yr")
sim.info$N.out <- length(sim.info$time.seq)
#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------
# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix
# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV.init, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)
# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)
#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------
# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation
save(sim.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Dynamic simulation 3 (start from previous initial conditions)
# check steady state
#=============================================================================
load(paste(BASENAME,"Transient","1.Rdata"))
# Initialisation simulation type 
sim.info$index <- 1
sim.info$code <- paste(BASENAME,"SteadySta")
sim.info$name <- paste(sim.info$code,sim.info$index)
# Initialisation parameter list
sim.start <- simulation.list[[length(simulation.list)]]
PL <- sim.start$PL
# Adapt parameter list
#PL$simulation.type <- "time.dependent"
PL$simulation.type <- "direct.steady.state"
PL <- initialise.parameters(PL)
# Initialise set of state variables
PL <- initialise.state.variables(PL)
# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)
# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)
# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)
# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)
#-------------------------------------------------------------------------------
# Steady simulation
#-------------------------------------------------------------------------------
#Actual steady state simulation
sim.info$run.time <- system.time({
  out <- steady.1D(y=sim.info$SV, func=model, parms=PL, names=PL$var.names, nspec=PL$N.var, positive=TRUE) #
  steady.state.reached <- attributes(out)$steady
  if (steady.state.reached) {SV <- out$y; print("steady state reached")} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)
#prepare matrix out
out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2
#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------
# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation
save(sim.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 
#unlink(paste(BASENAME,"Transient","1.Rdata"))
#=============================================================================
# Send email to say you are finished (and clean up after yourself)
#=============================================================================

if (!is.local){
  unlink(BASENAME,recursive=T)
}

if (is.local){
  
  require(mailR)
  
  send.mail(from = "rsebsebr@gmail.com",
            to = c("sebastiv@ucr.edu"),
            subject = paste("I think ", BASENAME, "is done!"),
            body = "Hello, I think the code you were running has finished. Best come check what mistakes you made ...",
            smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = "rsebsebr@gmail.com", passwd = "Rmail888", ssl = TRUE),
            authenticate = TRUE,  send = TRUE)
}
