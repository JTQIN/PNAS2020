###############################################################################
# Simulation file 
# CS model 
# Author: Filip Meysman
###############################################################################

# Source file containing model function
source("simple 1D model_v01.R")

# Source file containing plotting info function
source("plotting info.r")

# Number of simulations

model <- CSFe.model

#=============================================================================
# Z = 1 m - 1
#=============================================================================

load("00 SS-3 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 01 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 1
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 1 m - 2
#=============================================================================

load("01 01 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 01 Cflux sens-2"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 10 m - 1
#=============================================================================

load("01 01 Cflux sens-2 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 02 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 10
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 10 m - 2
#=============================================================================

load("01 02 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 02 Cflux sens-2"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 100 m - 1
#=============================================================================

load("01 02 Cflux sens-2 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 03 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 100
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 100 m - 2
#=============================================================================

load("01 03 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 03 Cflux sens-2"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 200 m - 1
#=============================================================================

load("01 03 Cflux sens-2 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 04 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 200
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 200 m - 2
#=============================================================================

load("01 04 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 04 Cflux sens-2"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 300 m - 1
#=============================================================================

load("01 04 Cflux sens-2 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 05 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 300
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 300 m - 2
#=============================================================================

load("01 05 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 05 Cflux sens-2"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 400 m - 1
#=============================================================================

load("01 05 Cflux sens-2 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 06 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 400
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 400 m - 2
#=============================================================================

load("01 06 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 06 Cflux sens-2"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 500 m - 1
#=============================================================================

load("01 06 Cflux sens-2 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 07 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 500
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 500 m - 2
#=============================================================================

load("01 07 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 07 Cflux sens-2"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 600 m - 1
#=============================================================================

load("01 07 Cflux sens-2 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 08 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 600
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 600 m - 2
#=============================================================================

load("01 08 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 08 Cflux sens-2"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 600 m - 1
#=============================================================================

load("01 08 Cflux sens-2 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 E3 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 650
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 650 m - 2
#=============================================================================

load("01 E3 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 E3 Cflux sens-2"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Z = 750 m - 1
#=============================================================================

load("01 08 Cflux sens-2 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 09 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 750
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 750 m - 2
#=============================================================================

load("01 09 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 09 Cflux sens-2"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 1000 m - 1
#=============================================================================

load("01 09 Cflux sens-2 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 10 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 1000
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 1000 m - 2
#=============================================================================

load("01 10 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 10 Cflux sens-2"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 1250 m - 1
#=============================================================================

load("01 10 Cflux sens-2 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 11 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 1250
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 1250 m - 2
#=============================================================================

load("01 11 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 11 Cflux sens-2"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Z = 1250 m - 3
#=============================================================================

load("01 11 Cflux sens-2 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 11 Cflux sens-3"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 1500 m - 1
#=============================================================================

load("01 11 Cflux sens-3 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 12 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 1500
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 1500 m - 2
#=============================================================================

load("01 12 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 12 Cflux sens-2"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Z = 1500 m - 3
#=============================================================================

load("01 12 Cflux sens-2 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 12 Cflux sens-3"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 



#=============================================================================
# Z = 1500 m - 5
#=============================================================================

load("01 12 Cflux sens-3 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 12 Cflux sens-4"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 1500 m - 3 - doesn't work
#=============================================================================

load("01 12 Cflux sens-4 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 12 Cflux sens-5"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 1550 m - 1
#=============================================================================

load("01 12 Cflux sens-2 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 E5 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 1550
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,10000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","10000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Z = 1550 m - 4
#=============================================================================

load("01 E5 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 E5 Cflux sens-2"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Z = 1600 m - 1
#=============================================================================

load("01 12 Cflux sens-2 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 E4 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 1600
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Z = 1600 m - 2
#=============================================================================

load("01 E4 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 E4 Cflux sens-2"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","3000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 1600 m - 3
#=============================================================================

load("01 E4 Cflux sens-2 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 E4 Cflux sens-3"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,10000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","10000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Z = 1600 m - 4
#=============================================================================

load("01 E4 Cflux sens-3 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 E4 Cflux sens-4"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 2000 m - 1
#=============================================================================

load("01 12 Cflux sens-2 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 13 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 2000
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Z = 2000 m - 2
#=============================================================================

load("01 13 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 13 Cflux sens-2"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","3000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 2000 m - 3
#=============================================================================

load("01 13 Cflux sens-2 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 13 Cflux sens-3"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","3000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Z = 2000 m - 4
#=============================================================================

load("01 13 Cflux sens-3 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 13 Cflux sens-4"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Z = 3000 m - 1
#=============================================================================

load("01 13 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 14 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 3000
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 



#=============================================================================
# Z = 3000 m - 2
#=============================================================================

load("01 14 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 14 Cflux sens-2"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 




#=============================================================================
# Z = 3000 m - 3
#=============================================================================

load("01 14 Cflux sens-2 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 14 Cflux sens-3"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 



#=============================================================================
# Z = 4000 m - 1
#=============================================================================

load("01 14 Cflux sens-3 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 15 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 4000
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 4000 m - 2
#=============================================================================

load("01 15 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 15 Cflux sens-2"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 



#=============================================================================
# Z = 5000 m - 1
#=============================================================================

load("01 15 Cflux sens-2 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 16 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 5000
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 5000 m - 2
#=============================================================================

load("01 16 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 16 Cflux sens-2"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 



#=============================================================================
# Z = 1300 m - 1
#=============================================================================

load("01 11 Cflux sens-3 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 17 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 1300
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,2000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Z = 1300 m - 2
#=============================================================================

load("01 17 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 17 Cflux sens-2"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,2000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Z = 1300 m - 3
#=============================================================================

load("01 17 Cflux sens-2 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 17 Cflux sens-3"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Z = 1325 m - 1
#=============================================================================

load("01 17 Cflux sens-3 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 18 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 1325
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,5000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","5000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 



#=============================================================================
# Z = 1325 m - 2
#=============================================================================

load("01 18 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 18 Cflux sens-2"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,5000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","2000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 1325 m - 3
#=============================================================================

load("01 18 Cflux sens-2 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 18 Cflux sens-3"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 



#=============================================================================
# Z = 1330 m - 1
#=============================================================================

load("01 18 Cflux sens-3 Steady.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 19 Cflux sens-1"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

Z <- 1330
PL$F.OC <- 10^(-0.50860503 - 0.00038900*Z)*1.8*1000

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,5000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","5000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 1330 m - 2
#=============================================================================

load("01 19 Cflux sens-1 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 19 Cflux sens-2"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,5000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","5000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Z = 1330 m - 3
#=============================================================================

load("01 19 Cflux sens-2 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 19 Cflux sens-3"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,5000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","5000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Z = 1330 m - 4
#=============================================================================

load("01 19 Cflux sens-3 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 19 Cflux sens-4"
sim.info$code <- "Transient"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"

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

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,15000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","15000 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Z = 1330 m - 5 
#=============================================================================

load("01 19 Cflux sens-4 Transient.Rdata")

# Initialisation simulation type 

sim.info$index <- "01 19 Cflux sens-5"
sim.info$code <- "Steady"
sim.info$name <- paste(sim.info$index,sim.info$code)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

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
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[2],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[2], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 



