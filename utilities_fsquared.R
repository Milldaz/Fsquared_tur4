# utilities.R - Extra functions
# /utilities.R

# Distributed under the terms of the EUPL-1.2

# ICES performance statistics {{{ ----

icestats <- list(

  # C
  C=list(~C, name="Catch (t)",
    desc="Catch inn tonnes"),

  # L
  L=list(~L, name="Landings (t)",
    desc="Landings in tonnes"),

  # D
  D=list(~L, name="Discards (t)",
    desc="Discards in tonnes"),

  # SB
  SB=list(~SB, name="SB (t)",
    desc="Spawning stock biomass"),

  # R
  R=list(~R, name="Recruits (1e3)",
    desc="Recruitment in numbers"),
  
  # F
  F=list(~F, name="bar(F)",
    desc="Fishing mortality"),

  # cv(C)
  cvC=list(~sqrt(iterVars(C)) / iterMeans(C), name="cv(C)",
    desc="CV of catch per year"),

  # AVVC
  AAVC=list(~abs(C[, -1] - C[, -dim(C)[2]]) / C[, -dim(C)[2]],
    name="AAV(C)", desc="Average annual variability in catch"),

  # IACC
  IACC=list(~100 * (C[, -1] - C[, -dim(C)[2]]) / C[, -dim(C)[2]],
    name="IAC(C)",
    desc="Percentage inter-annual change in catch"),

  # IACB
  IACB=list(~100 * (SB[, -1] - SB[, -dim(SB)[2]]) / SB[, -dim(SB)[2]], name="IAC(B)", desc="Percentage inter-annual change in biomass"),

  # P(SB<SBlim)
  PBlim=list(~iterMeans((SB/Blim) < 1)  , name="P(SB<SB[lim])",
    desc="Probability that spawner biomass is below Blim"),

  # P(SB>SBtrigger)
  PBtrigger=list(~iterMeans((SB/Btrigger) > 1), name="P(SB>B[trigger])",
    desc="Probability that spawner biomass is above Btrigger"),

  # P(SB < SBlim) at least once
  risk2 = list(~iterMeans((SB / Blim) < 1) > 0,
    name="once(P(SB<B[limit]))",
    desc="ICES Risk 2, probability that spawner biomass is above Blim once")
)

# }}}

### Define new stock recruit relationships
## determine segreg model with Blim breakpoint and (roughly) geomean rec above this
segregBlim  <- function(ab, ssb) log(ifelse(ssb >= Blim, ab$a * Blim, ab$a * ssb))

