# toyepidemic
An R package for simulating within-host and between-host dynamics during epidemics

This package provides methods for simulating the evolution of pathogen traits. The simulation method assumes a classical additive model for a trait in which the trait-value is formed as a sum of single-locus and epistatic effects from the pathogen genotype, the host immune system type and a white noise term representing the combination of unknown non-heritable factors such as unknown pathogen and host immune system loci, host environment and measurement error. It is assumed that the genotype of the pathogen can mutate during an infection as a result of neutral genetic drift or within-host selection for an optimal trait-value. At the between-host level various events, such as transmission to susceptible hosts, diagnosis and immunity are simulated according to a Susceptible-Infected-Recovered (SIR) epidemiological model. The rates at which these events occur as well as contact rates between infected and susceptible hosts can be specified as constants or as functions of the pathogen trait. The resulting epidemic is represented as a data.table and a rooted phylogenetic tree (a phylo object). The data contains the history of within- and between- host events. The package provides plotting functions to visualize the epidemic over time.

# License

The toyepidemic R-package is licensed under version 3 of the [GNU General Public License](http://www.gnu.org/licenses/.
