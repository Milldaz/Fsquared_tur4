# model_fsquared.R - runs steps to estimate ICES ref pts using Fsquared
# Fsquared/model_fsquared.R

#===============================================================================
# Fsquared - estimation of F reference points for the ICES advice rule
# Authors: participants of ICES MGWG workshop on the validation of new tool for refpts
#	  Iago Mosqueira (WMR) <iago.mosqueira@wur.nl>
#	  Ernesto Jardim (IPMA) <ernesto.jardim@ipma.pt>
#	  John Trochta (IMR) <john.tyler.trochta@hi.no>
#	  Arni Magnusson (SPC) <arnim@spc.int>
#	  Max Cardinale (SLU) <massimiliano.cardinale@slu.se>
#	  Dorleta Garcia (ICES) <dorleta.garcia@ices.dk>
#	  Colin Millar (ICES) <colin.millar@ices.dk>
#
# Distributed under the terms of the EUPL 1.2
#===============================================================================

# Clean start
rm(list=ls())

# Tested using R version 4.5.2
# load libraries
library(TAF)
library(mse)
library(msemodules)
library(FLasher)
library(FLSRTMB)  # to estimate SR parameters and add estimation uncertainty

# Check correct versions are installed:
check.software(full=TRUE)

# Load performance stats function
source("utilities_fsquared.R")

# LOAD stock assessment results, 'run' is output FLStock
# if you have your own assessment load it, for simplicity rename the output
# to "run", but can use other name, just make sure all calls to "run" 
# are updated
load("model/model.rda")
# rename
run <- stk_fit; rm(stk_fit)

#DM Check FLstock object
# Check for NAs in ns and wts etc.
catch(run)
catch.n(run)
catch.wt(run)
discards(run)
discards.n(run)
discards.wt(run)
landings(run)
landings.n(run)
landings.wt(run)
stock(run)
stock.n(run)
stock.wt(run)
m(run)
mat(run)
harvest(run)
# trim back to 2024 (no 2025 catch info), add stock name
#run <- trim(run, year=1978:2024)
name(run) <- "whg.27.47d"

#===============================================================================
# SETUP
#===============================================================================

# Name of the stock
stkname <- name(run)
# TODO: Recruitment models to be used in the OM conditioning
srmodels <- c("segreg") # segreg, bevholt, ricker
# Initial year of projections
iy <- dims(run)$maxyear
# Years to be used to compute SPR0 for stock-recruitment model, last 5
spryrs <- seq(iy - 5, iy)
# Data lag 
dl <- 2
# Management lag
ml <- 1 
# Assessment frequency
af <- 1
# Data year
dy <- iy - dl
# Final year
fy <- iy + 30 #DM changed after check on stable biomass below
# Years to compute probability metrics
pys <- seq(fy - 9, fy) #DM changed to last ten years
# How many years from the past to condition the future
conditioning_ny <- 5
# CV for SSB to add uncertainty in the shortcut estimator
bcv_sa <- 0.5
# CV for F to add uncertainty in the shortcut estimator
fcv_sa <- 0.5
# Years for geometric mean in short term forecast
recyrs_mp <- -40
# TODO: Blim and Btrigger
Blim <- 150477
Btrigger <- 209100
refpts <- FLPar(c(Blim = Blim, Btrigger = Btrigger, Fmsy = NA))  #DM added Fmsy
# TODO: no. of cores to use in parallel, defauls to 2/3 of those in machine
cores <- round(availableCores() * 0.6)
# TODO: F search grid
#DM Start with rough grid, then fine tune
#fg_mp <- seq(0.25, 0.70, length=cores)
#fg_mp <- seq(0.1, 0.5, by=0.01)
fg_mp <- seq(0.23, 0.70, length=6*cores)
# Number of iterations (minimum of 50 for testing, 500 for final)
#DM it <- cores
it <- 50
# it <- max(25, cores * 25)
# Random seed
set.seed(987)

# PARALLEL setup via doFuture
if(os.linux()) {
  plan(multicore, workers=cores)
} else {
  plan(multisession, workers=cores)
}

options(doFuture.rng.onMisuse="ignore")

#===============================================================================
# OM conditioning
#===============================================================================

# Stock-recruitment relationship(s)
# The file utilities_fsquared.R has code examples to condition the OM
# using stock recruitment parameters estimated by the stock assessment
# model, like SS3, SAM and a4a.

# BOOTSTRAP and SELECT model by largest logLik
srpars <- bootstrapSR(run, iters=it, spr0=mean(spr0y(run)[, ac(spryrs)]),
  models=srmodels, verbose=FALSE)

# GENERATE future deviances: lognormal autocorrelated
srdevs <- rlnormar1(sdlog=srpars$sigmaR, rho=srpars$rho, years=seq(dy, fy),
  bias.correct=FALSE)

# BUILD FLom, OM FLR object
om <- FLom(stock=propagate(run, it), refpts=refpts, model="segreg",
  params=srpars, deviances=srdevs, name=stkname)

# TODO: SETUP om future: average of most recent years set by conditioning_ny
om <- fwdWindow(om, end=fy, nsq=conditioning_ny)

#===============================================================================
# diagnostic(s)
#===============================================================================

# TEST Blim with F=0 projection
f0 <- fwd(om, control=fwdControl(quant='fbar', value=0, year=seq(iy + 1, fy)))

# COMPUTE P(SB<Blim) in time
# when F=0 this probability should be 0, otherwise is a sign BLim is not 
# coherent with the stock-recruitment model and parameters 
performance(f0, statistics=icestats["PBlim"])[year %in% seq(iy, fy, by=5),
  .(PBlim=mean(data)), by=year]

# COMPUTE Inter-annual variability of biomass
# to check when biomass stabilizes and set number of years for projections
# leave about 10 years of stable biomass to compute metrics
# review MSE setup section to avoid projecting for longer than necessary
# rerun setup and objects if you changed the final year of projections
dt0 <- performance(f0, statistics=icestats["IACB"])[year %in% iy:fy,
  .(IACB=mean(data)), by=year]

plot(dt0, type="l", main="Inter-annual changes in biomass")

#===============================================================================
# MP
#===============================================================================

# SET intermediate year + start of runs, lags and frequency
mseargs <- list(iy=iy, fy=fy, data_lag=dl, management_lag=ml, frq=af)

# SET shortcut estimator uncertainty: F and SSB deviances and auto-correlation
# Note your SSB deviances and auto-correlation have very little impact on the P(SB<Blim)
sdevs <- shortcut_devs(om, SSBcv=bcv_sa, Fcv=fcv_sa, Fphi=0)

# SETUP constant F rule
frule <- mpCtrl(
  
  # (est)imation method: shortcut.sa + SSB deviances
  est = mseCtrl(method=shortcut.sa,
                args=list(SSBdevs=sdevs$SSB)),
  
  # hcr: constant F (Btrigger=0)
  hcr = mseCtrl(method=hockeystick.hcr,
                args=list(lim=0, trigger=0, target=refpts(om)$Fmsy, min=0,   #DM should be with refpts(om)$Fmsy (not defined yet)
                          metric="ssb", output="fbar")),
  
  # (i)mplementation (sys)tem: tac.is (C ~ F)
  isys = mseCtrl(method=tac.is, args=list(recyrs=recyrs_mp, Fdevs=sdevs$F))
)

# SETUP standard ICES advice rule
arule <- mpCtrl(

  # (est)imation method: shortcut.sa + SSB deviances
  est = mseCtrl(method=shortcut.sa,
    args=list(SSBdevs=sdevs$SSB)),

  # hcr: hockeystick (fbar ~ ssb | lim, trigger, target, min)
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0, trigger=refpts(om)$Btrigger, target=refpts(om)$Fmsy, min=0,   #DM refpts(om)$Fmsy (not defined yet)
    metric="ssb", output="fbar")),

  # (i)mplementation (sys)tem: tac.is (C ~ F)
  isys = mseCtrl(method=tac.is, args=list(recyrs=recyrs_mp, Fdevs=sdevs$F))
)

#===============================================================================
# Run simulations
#===============================================================================

# RUN constant F rule over Ftarget grid
fgrid <- mps(om, ctrl=frule, args=mseargs, hcr=list(target=fg_mp),
             names=paste0("F", fg_mp))
# PLOT
plot(om, fgrid)
# COMPUTE average performance over pys
performance(fgrid) <- performance(fgrid, statistics=icestats["PBlim"], year=pys,
                                  type="frule")

# OR ... RUN over Ftarget grid and return only performance stats
fgridstat <- mps(om, ctrl=frule, args=mseargs, hcr=list(target=fg_mp),
             names=paste0("F", fg_mp), statistics=icestats, type="frule")

# RUN ICES advice rule over Ftarget grid
#hcrfgrid <- mps(om, ctrl=arule, args=mseargs, hcr=list(target=fg_mp),
#  names=paste0("F", fg_mp))
# PLOT
#plot(om, hcrfgrid)
# COMPUTE average performance over pys
#performance(hcrfgrid) <- performance(hcrfgrid, statistics=icestats["PBlim"], year=pys,
#  type="arule")

# FIND Fpa: Ftarget that gives mean P(B < Blim) = 5%
tune <- tunebisect(om, control=arule, args=mseargs,
  tune=list(target=0.3 * c(0.5, 1.5)),
  statistic=icestats["PBlim"], prob=0.05, tol=0.005, years=pys)

# PLOT
plot(om, tune)

# COMPUTE performance
performance(tune) <- performance(tune, statistics=icestats, type="arule", run="tune")

# CHECK Fpa value
args(control(tune, "hcr"))$target

# SAVE
save("fgrid", "fgridstat","om", "tune", file="model/fsquared.rda")

