initialise.plotting.info <- function(plot.info,PL)
{
  #=============================================================================
  # Preparation of the way the simulation output is plotted
  #=============================================================================
  
  plot.info$spec.mfrow <- c(4,4)
  plot.info$spec.mar <- c(2,4,5,2)+0.1
  
  plot.info$reac.mfrow <- c(4,4)
  plot.info$reac.mar <- c(2,4,5,2)+0.1
  
  # conversion from mol C cm-3 of solids to percentage of solids [%]
  
  MW.C <- 12
  MW.CaCO3 <- 100.09
  C.fac <- 100/(1E06*PL$rho.sed) # from [umol cm-3 of solids] to percentage [g C / g total sed %]
  FeS.fac <- 1/(PL$rho.sed) # from [umol cm-3 of solids] to [umol g-1 dry wt]
  FeOOH.fac <- 1/(PL$rho.sed) # from [umol cm-3 of solids] to [umol g-1 dry wt]
  X.Fe.fac <- 1/(PL$rho.sed) # from [umol cm-3 of solids] to [umol g-1 dry wt]
   
  plot.info$SL[["OC"]]$value.fac   <- C.fac*MW.C 
  plot.info$SL[["OC"]]$value.label <- "OC [%]"
  
  plot.info$SL[["FeOOH"]]$value.fac   <- FeOOH.fac 
  plot.info$SL[["FeOOH"]]$value.label <- "FeOOH [µmol g-1]"
  
  plot.info$SL[["FeS"]]$value.fac   <- FeS.fac 
  plot.info$SL[["FeS"]]$value.label <- "FeS [µmol g-1]"
  
  plot.info$SL[["O2"]]$depth.lim <- c(0,1)
  
  plot.info$SL[["OC"]]$visible    <- TRUE
  plot.info$SL[["FeOOH"]]$visible <- TRUE
  plot.info$SL[["FeS"]]$visible   <- TRUE
  
  plot.info$SL[["H2O"]]$visible  <- TRUE
  plot.info$SL[["H"]]$visible    <- TRUE
  plot.info$SL[["SO4"]]$visible  <- TRUE
  plot.info$SL[["HCO3"]]$visible <- TRUE
  plot.info$SL[["HS"]]$visible   <- TRUE
  plot.info$SL[["Fe"]]$visible   <- TRUE
  
  plot.info$RL[["Cmin"]]$visible <- TRUE
  plot.info$RL[["DIR"]]$visible  <- TRUE
  plot.info$RL[["SR"]]$visible   <- TRUE
  plot.info$RL[["CSO"]]$visible  <- TRUE
  plot.info$RL[["FIO"]]$visible  <- TRUE
  plot.info$RL[["ISP"]]$visible  <- TRUE
  
  
  return(plot.info)
}  
