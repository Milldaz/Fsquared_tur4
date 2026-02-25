# report.R - DESC
# fsquared/report.R

#===============================================================================
# Fsquared - estimation of F reference points for the ICES advice rule
# Authors: participants of workshop on the validation of new tool for refpts
#   Iago Mosqueira (WMR) <iago.mosqueira@wur.nl>
#	  Ernesto Jardim (IPMA) <ernesto.jardim@ipma.pt>
#	  John Trochta (IMR) <john.tyler.trochta@hi.no>
#	  Arni Magnusson (SPC) <arnim@spc.int>
#	  Max Cardinale (SLU) <massimiliano.cardinale@slu.se>
#	  Dorleta Garcia (ICES) <dorleta.garcia@ices.dk>
#	  Colin Millar (ICES) <colin.millar@ices.dk>
#
# Distributed under the terms of the EUPL 1.2
#===============================================================================


library(mse)
library(msemodules)
library(mseviz)

library(FLCore)
library(ggplotFL)

# model_fsquaredf.R {{{ ----

load("model/fsquared.rda")

#DM need:
# Table of the decisions (setting values chosen) made for the runs to provide a table to include in benchmarks
# Need yield curve plots

# SR plot
# FIT models
fits <- srrTMB(as.FLSRs(run, models=c("segreg")), 
               spr0=mean(spr0y(run)))
# PLOT
taf.png(file="data_srrfits.png")
plotsrs(fits)
dev.off()

taf.png(file="data_srbootstrap.png", width=1400)
plot_bootstrapSR(fits, srpars)
dev.off()


# 1. Extract the FLSR (Stock-Recruitment) object
stock_recruit <- sr(om)
plot(stock_recruit)

ssb <- seq(0, 5e5, length = 1e3)
srplot <- data.frame(ssb = ssb, rec = NA)
srdata <- data.frame(ssb = fits_all$default$segreg@ssb[drop=T], 
                     rec = fits_all$default$segreg@rec[drop=T], scenario = 'all')

  a <- median(srpars['a',drop=T])
  b <- median(srpars['b',drop=T])
  srplot[, 'rec'] <- ifelse(ssb <= b, a * ssb, a * b)

#taf.png(file=paste0('report/', stock, "_srr",".png"))
ggplot(srplot, aes(ssb, rec)) + geom_line(linewidth = 2) + 
  geom_point(aes(ssb, rec), srdata, col = 'black', size = 2) + 
  ggtitle(paste0("Stock recruitment relationships: ", stock))
#dev.off()



ggplot(aes(ssb, rec), data = model.frame(FLQuants(stock_recruit, "ssb", "rec"))) +
  geom_point() +
  geom_smooth(method = "loess") + 
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
  labs(x = "Spawning Stock Biomass (SSB)", y = "Recruitment", title = "Stock-Recruit Relationship")

##########

# plot HCR
plot_hockeystick.hcr(arule$hcr, labels=c(lim="Blim", trigger="Btrigger",
                                         min="", target="Ftarget")) +
  xlab(expression(hat(SSB))) + ylab(expression(bar(F)))


####################
# PLOT PBlim over F values
dat <- performance(fgrid)[, .(data=mean(data)), by=.(mp)]
dat <- performance(hcrfgrid)[, .(data=mean(data)), by=.(mp)]

ggplot(dat, aes(x=mp, y=data)) +
  geom_point() +
  scale_x_discrete(name="F", labels=fg_mp,
    breaks=function(x) x[c(TRUE, FALSE)]) +
    ggtitle("Average probability of SSB falling below Blim over selected years by F level") +
   geom_hline(yintercept = 0.05, linewidth = 0.5, linetype = "dashed", color="gray50")

# PLOT TOs for tune: C ~ IACC, PBlim, riosk2, Pbtrigger

plotTOs(performance(tune), x="C", y=c("IACC", "PBlim", "risk2", "PBtrigger")) +
  theme(legend.position="none")

plotTOs(performance(tune), x="F", y=c("IACC", "PBlim", "risk2", "PBtrigger")) +
  theme(legend.position="none")

# }}}
