# toyepidemic
An R package for simulating within-host and between-host dynamics during epidemics

This package provides methods for simulating the evolution of pathogen traits. The simulation method assumes a classical additive model for a trait in which the trait-value is formed as a sum of single-locus and epistatic effects from the pathogen genotype, the host immune system type and a white noise term representing the combination of unknown non-heritable factors such as unknown pathogen and host immune system loci, host environment and measurement error. It is assumed that the genotype of the pathogen can mutate during an infection as a result of neutral genetic drift or within-host selection for an optimal trait-value. At the between-host level various events, such as transmission to susceptible hosts, diagnosis and immunity are simulated according to a Susceptible-Infected-Recovered (SIR) epidemiological model. The rates at which these events occur as well as contact rates between infected and susceptible hosts can be specified as constants or as functions of the pathogen trait. The resulting epidemic is represented as a data.table and a rooted phylogenetic tree (a phylo object). The data contains the history of within- and between- host events. The package provides plotting functions to visualize the epidemic over time.

# How to use it?
================

Here is a quick example on how to use the package:

``` r
library(toyepidemic)

numsAllelesAtSites <- c(4, 3, 2)
# genotype encodings with their allele contents
genotypes <- generateGenotypes(numsAllelesAtSites)

# number of host-types
n <- 6
pe <- runif(n)
pe <- pe/sum(pe)

# General host-type x strain effects (expected phenotypes for genotype by environment combinations)
GEVs <- matrix(NA, nrow=n, ncol=nrow(genotypes))

for(g in 1:nrow(genotypes)) {
  GEVs[, g] <- rnorm(n=n, mean=2+2/nrow(genotypes)*(g), sd=0.4)
}

matplot(GEVs, type='l', lty=1, col=1:6)

sigmae <- .6

timeStep <- 0.05

pg.init <- rep(0, nrow(genotypes))
pg.init[1] <- 1

sde <- rep(sigmae, n)

# death-rate as a function of viral load and natural death rate mu
rateDie <- function(z, mu) {
  V <- 10^z
  Dmin <- 2
  Dmax <- 25*12
  D50 <- 10^3
  Dk <- 1.4
  (V^Dk+D50^Dk)/(Dmin*(V^Dk+D50^Dk)+((Dmax-Dmin)*D50^Dk)) + mu
}

meanRateDie <- function(z, mu) {
  rep(0.01, length(z))
}

# infectionrate as a function of viral load and rate of risky contacts
rateInfect <- function(z, rateContact) {
  V <- 10^z
  Emin <- .3
  Emax <- .6
  E50 <- 10^3
  Ek <- 1.4
  E <- Emin+(Emax-Emin)*V^Ek/(V^Ek+E50^Ek)
  E*rateContact
}

meanRateInfect <- function(z, rateContact) {
  rep(0.45*rateContact, length(z))
}

# per locus mutation rates
rateMutate <- function(GEValues, es, envs, genes) {
  z <- es+GEValues[cbind(envs, genes)]
  V <- 10^z
  Mmin <- 0.00
  Mmax <- 0.2
  M50 <- 10^3
  Mk <- 1.4
  Mmin+(Mmax-Mmin)*V^Mk/(V^Mk+M50^Mk)
}

meanRateMutate <- function(GEValues, es, envs, genes) {
  rep(0.01, length(genes))
}

rateSampleNeutral <- 1/(4*12)

# maximum time before starting graceful fadeout of the epidemic
maxTime <- 1000

# continue the epidemic outbreak until reaching ...
maxNTips <- 2000

# is the special environmental effect unique for each pathogen genoetype in an individual
eUniqForEachG <- TRUE

# are only beneficial (i.e. increasing the trait-value) mutations allowed
selectWithinHost <- TRUE

# time to continue the simulation of transmission after reaching maxNTips (this was introduced
# in order to study post-outbreak dynamics, i.e. epidemic waves after exhaustion of the
# susceptible pool)
expandTimeAfterMaxNTips <- 4

#initial population size (equilibrium in the absence of disease)
N <- 1e5

# setting the between-host and within-host dynamics of the simulation:
# natural death rate
mu <- 1/850

# constant birth rate that maintains this equilibrium
nu <- ifelse(is.finite(N), mu*N, 0)

rateContact <- 1

rateSample <- rateSampleNeutral

selectWithinHost <- TRUE

epidemic <- simulateEpidemic(
  Ninit=N, nu=nu, mu=mu, pe=pe, sde=sde, pg.init=pg.init, GEValues=GEVs,
  rateContact=1/6, rateInfect=rateInfect, rateDie=rateDie,
  rateSample=rateSampleNeutral,
  rateMutate=meanRateMutate,
  numsAllelesAtSites=numsAllelesAtSites, eUniqForEachG=eUniqForEachG,
  selectWithinHost=TRUE,
  timeStep=timeStep, maxTime=maxTime, maxNTips=maxNTips,
  expandTimeAfterMaxNTips=expandTimeAfterMaxNTips,
  process="select/select")

tr1 <- extractTree(epidemic, tips=1:1000)

plot(ape::ladderize(tr1), show.tip.label = FALSE)

# extract the population of all sampled individuals
pop1 <- extractPop(epidemic)
pop1[, z:=calcValue(env, gene, e, GEValues = GEVs)]

drc1 <- extractDRCouples(epidemic)

```

# License

The toyepidemic R-package is licensed under version 3 of the [GNU General Public License](http://www.gnu.org/licenses/.
