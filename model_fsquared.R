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
# catch(run)
# catch.n(run)
# catch.wt(run)
# discards(run)
# discards.n(run)
# discards.wt(run)
# landings(run)
# landings.n(run)
# landings.wt(run)
# stock(run)
# stock.n(run)
# stock.wt(run)
# m(run)
# mat(run)
# harvest(run)
# trim back to 2024 (no 2025 catch info)?, add stock name
#run <- trim(run, year=1978:2024)
name(run) <- "whg.27.47d"

#===============================================================================
# SETUP  #DM when all changes implement, reorganise so all settings that users need to specify are at the start, with no later changes needed
#===============================================================================

# Name of the stock  #DM Should also have a run name to allow users to test different settings and save results
stkname <- name(run)
# TODO: Recruitment models to be used in the OM conditioning
srmodels <- c("segreg") # segreg, bevholt, ricker
# Initial year of projections
iy <- dims(run)$maxyear
# Years to be used to compute SPR0 for stock-recruitment model, last 5  #DM Need some explanation here
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
fy <- iy + 30 #DM useful check on stable biomass below, but maybe need some guidance on minimum period (e.g. 1-2 x lifespan + summary period length)
# Years to compute probability metrics
pys <- seq(fy - 9, fy) #DM changed to last ten years (was 6). Need ACOM decision on best practice
# How many years from the past to condition the future
conditioning_ny <- 5 #DM In EQSIm have bio and sel separate 
#DM Currently only average used as a constant for whole forecast period. Add resample option?

# CV for SSB to add uncertainty in the shortcut estimator
bcv_sa <- 0.5 #DM Need informed defaults and ways to estimate this for specific stocks. How does this interact with management/data lag?
# CV for F to add uncertainty in the shortcut estimator
fcv_sa <- 0.5 #DM EQSIm is 0.243 default. Need informed defaults and ways to estimate this for specific stocks.  How does this interact with management/data lag?
Fphi_sa <- 0 #DM added parameter here instead of hardwiring 0 later in the code. Shouldn't be zero.
# Years for geometric mean in short term forecast
recyrs_mp <- -40 #DM need alternatives (e.g. whole time series)

# TODO: Blim and Btrigger
#DM add a section for Blim and Bpa estimation
#DM Add spasmodic test in the code 
#DM Add Bemp estimation
Blim <- 150477
Btrigger <- 209100 #DM Need to be able to compare Bpa with 5th percentile of SSB at Fmsy
refpts <- FLPar(c(Blim = Blim, Btrigger = Btrigger, Fmsy = NA))  #DM added Fmsy


# TODO: no. of cores to use in parallel, defauls to 2/3 of those in machine
cores <- round(availableCores() * 0.6)
# TODO: F search grid  #DM this will need to be done in steps, starting broad then finetuning around Fmsy (steps of 0.01)
#DM Start with rough grid, then fine tune
#fg_mp <- seq(0.25, 0.70, length=cores)
#fg_mp <- seq(0.1, 0.5, by=0.01)
fg_mp <- seq(0.23, 0.78, length=7*cores)
# Number of iterations (minimum of 50 for testing, 500 for final)
#DM it <- cores
it <- 104 #Simualtion time is long, so use a few for initial work, but run final using a lot more (need ACOM guidance)
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
# The file utilities_fsquared.R has code examples to condition the OM #DM it does not
# using stock recruitment parameters estimated by the stock assessment
# model, like SS3, SAM and a4a.

# BOOTSTRAP and SELECT model by largest logLik
srpars <- bootstrapSR(run, iters=it, spr0=mean(spr0y(run)[, ac(spryrs)]),
  models=srmodels, verbose=FALSE)

# GENERATE future deviances: lognormal autocorrelated
srdevs <- rlnormar1(sdlog=srpars$sigmaR, rho=srpars$rho, years=seq(dy, fy),
  bias.correct=FALSE)

#DM need to be able to implement segreg with a specified breakpoint (e.g. Blim or Bloss)
#DM Need plots of SR (and spread of deviances)
#DM what was wrong with eqsr, which had both of the above?

# BUILD FLom, OM FLR object
om <- FLom(stock=propagate(run, it), refpts=refpts, model=srmodels,  #DM should be srmodels not "segreg"
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
# Note your SSB deviances and auto-correlation have very little impact on the P(SB<Blim) #DM not universally true
sdevs <- shortcut_devs(om, SSBcv=bcv_sa, Fcv=fcv_sa, Fphi=Fphi_sa)  #DM added Fphi_sa here

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

#DM Need to extract Fmsy, Flow and Fupper (i.e. 95% of max yield) from fgrid 

medC_conF <- performance(fgrid, statistics=icestats["medC"], year=pys, type="frule")
aggregate(data ~ run, data = medC_conF, FUN = median, na.rm = TRUE)

meanC_conF <- performance(fgrid, statistics=icestats["meanC"], year=pys, type="frule")
aggregate(data ~ run, data = meanC_conF, FUN = mean, na.rm = TRUE)

# OR ... RUN over Ftarget grid and return only performance stats
#fgridstat <- mps(om, ctrl=frule, args=mseargs, hcr=list(target=fg_mp),
#             names=paste0("F", fg_mp), statistics=icestats, type="frule")

# RUN ICES advice rule over Ftarget grid #DM not sure this is needed for all of fg_mp (probably just Fmsy)
hcrfgrid <- mps(om, ctrl=arule, args=mseargs, hcr=list(target=fg_mp),
  names=paste0("F", fg_mp))
# PLOT
plot(om, hcrfgrid)
# COMPUTE average performance over pys
performance(hcrfgrid) <- performance(hcrfgrid, statistics=icestats["PBlim"], year=pys,
  type="arule")

#DM Need to extract 5th percentile of from hcrfgrid (technically shoul dcome form a run with assessment error, check with ACOM)

medC_AR <- performance(hcrfgrid, statistics=icestats["medC"], year=pys, type="arule")
aggregate(data ~ run, data = medC_AR, FUN = median, na.rm = TRUE)

meanC_AR <- performance(hcrfgrid, statistics=icestats["meanC"], year=pys, type="arule")
aggregate(data ~ run, data = meanC_AR, FUN = mean, na.rm = TRUE)


# FIND Fpa: Ftarget that gives mean P(B < Blim) = 5%
tune <- tunebisect(om, control=arule, args=mseargs,
  tune=list(target=0.3 * c(0.5, 1.5)),  #DM change to range around Fmsy estimated above?
  statistic=icestats["PBlim"], prob=0.05, tol=0.005, years=pys)
#DM - this is the quickest way to get Fpa, but need to think about default settings

# PLOT
plot(om, tune)

# COMPUTE performance
performance(tune) <- performance(tune, statistics=icestats, type="arule", run="tune")

# CHECK Fpa value
args(control(tune, "hcr"))$target

# SAVE
#save("fgrid", "hcrfgrid","om", "tune", file="model/fsquared.rda")
save("fgrid", "fgridstat","hcrfgrid","om", "tune", file="model/fsquared.rda")



