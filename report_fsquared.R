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

# model_fsquaredf.R {{{ ----

load("model/fsquared.rda")

# PLOT PBlim over F values
dat <- performance(fgrid)[, .(data=mean(data)), by=.(mp)]

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
