#' ---
#' title: "R Stan examples: Treatment for Migraines"
#' author: 
#' - name: Georges Monette
#'   affiliation: York University
#' date: "`r format(Sys.time(), '%B %d, %Y at %H:%M')`"
#' output:
#'    html_document:
#'      toc: true
#'      toc_depth: 6
#'      toc_float: true
#' ---
#' 
#+ knitr_setup, include=FALSE
knitr::opts_knit$set()
knitr::opts_chunk$set(cache = TRUE, eval = FALSE)
#'
#' 
#' ## References
#'  
#' - [Monte Carlo for ODE model](https://www.mathworks.com/matlabcentral/fileexchange/35725-monte-carlo-markov-chain-for-inferring-parameters-for-an-ordinary-differential-equation-model?focused=5226125&tab=example)
#' 
#' 
#' The model:
#' 
#' \[\begin{aligned}
#' \frac{{d}}{{dt}} T(t)=  & \lambda  - d_T T(t) - \beta T(t)V(t) \\ 
#' \frac{{d}}{{dt}} I(t)=  & \beta T(t)V(t) -d_I I(t) \\ 
#' \frac{{d}}{{dt}} V(t)=  & pI(t) - cV(t) \\ 
#' \end{aligned} \]
#' 
#' where
#' 
#' - $T$ is expected uninfected cell count
#' - $I$ is expected infected cell count
#' - $V$ is the expected virus count
#' - $\lambda$ is the autonomous rate of cell generation
#' - $d_T$: death rate for uninfected cells
#' - $d_I$: death rate for infected cells
#' - $p$ production rate of viruses by infected cells
#' - $c$ clearance rate for viruses
#' 
#' In the simple model used here, we assume that the process operates on
#' a very large number of cells and the expected values follows the ODE
#' exactly. We incorporate measurement error due to measurements being based
#' on a small sample of blood. If the size of the sample is consistent and if
#' the measured quantities are cell counts in the sample, then the observed
#' value at time $t$, $T_t^o$ has a Poisson distribution with mean $T_t$ and,
#' thus, standard deviation $\sqrt{T_t}$. 
#' 
#' If the observed value is a concentration per volume $C$ measured from
#' a count of a sample of volume $s$, where $T_t$ is the expected 
#' concentration per volume $C$, then we can use
#' $$T+t^o | T_t \sim Normal(T_t, \sqrt{\frac{C}{s}T_t})$$
#' where $\sqrt{\frac{C}{s}T_t})$ is the standard deviation of the normal
#' distribution.  Note that if $s$ varies from observation to observation, it
#' can be included in the data.
#' 
#' The model below can be adapted to include a Gaussian process to add noise to
#' the ODE. Perhaps more interestingly, it would be equally easy to include a
#' thick-tailed Cauchy process that would reduce the influence of occasional
#' 'wild' data.  How well this might work in practice can be explored with simulations.
#' 
#'
#' ## Load packages 
#+ 
library(spida2) # devtools::install_github('gmonette/spida2')

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
if(.Platform$OS.type == 'windows') windowsFonts(Arial=windowsFont("TT Arial")) 

library(lattice)
library(latticeExtra)


#'
#' ## Simulation function
#'
#' The following function simulates `J`` mice that start with no infected cells
#' and an injection of virus varying from `Vmax/J`` to `Vmax`.

sim <- function(){
  Nt <- 2000    # mean initial
  Vmax <- 200;  # maximum initial viral load
  J <- 8       # number of 'mice'
  lambda <- 10   # autonomous new cell generation
  dT <- .02      # death rate of uninfected cells
  dI <- .04      # death rate of infected cells
  Tmean <- lambda /dT # initial cell count at equilibrium
  beta <- .0002  # infection rate   
  p <- .0005    # production rate  
  c <- .003     # viral clearance rate  
  Imat <- Vmat <- Tmat <- matrix(NA,Nt,J)
  Tmat[1,] <- Tmean
  Imat[1,] <- 0  # start uninfected
  Vmat[1,] <- Vmax *(1:J)/J
  for(t in 1:(Nt-1)){
    Tmat[t+1,] <- pmax(Tmat[t,] + lambda - dT * Tmat[t,] - beta * Tmat[t,]*Vmat[t,],0)
    Imat[t+1,] <- pmax(Imat[t,] + beta * Tmat[t,] * Vmat[t,] - dI * Imat[t,],0)
    Vmat[t+1,] <- pmax(Vmat[t,] + p * Imat[t,] - c * Vmat[t,],0)
  }
  ret <- data.frame(CellUI = c(Tmat), CellInf = c(Imat), virus = c(Vmat),
                    time = rep(1:Nt,J), id = rep(1:J, each = Nt),
                    virus0 = rep(Vmax*(1:J)/J, each = Nt ))
  ret
}
gd(lwd=2)
xyplot(CellUI + CellInf + virus ~ time| as.factor(virus0),
                    sim(), type = 'l', 
       auto.key = list(space = 'top',columns=3,
                       text = c('uninfected','infected','virus'),
                       lines=T,points=F),
       ylab = 'cell counts',
       scales = list(x = list(rot=45,alternating = F)),
       sub = 'infection dynamics with varying initial virus counts')
dd <- sim()
head(dd)
tail(dd)
library(latticeExtra)



#+
knitr::knit_exit()
#+
print(readLines('Cell.stan'),quote =F)
system.time(cell_dso <- stan_model('Cell.stan'))
#'
#' ## Prepare data list
#'

dim(dd)
# Suppose counts obtained every 10 time units
# dobs <- dd[seq(1,nrow(dd),by = 10),]
dobs <- dd
head(dobs)
dobs <- within(dobs, {
  Tc <- CellUI + sqrt(CellUI)*rnorm(CellUI)
  Ic <- CellInf + sqrt(CellInf)*rnorm(CellInf)
  Vc <- virus + sqrt(virus)*rnorm(virus)
})
dat <- with(dobs, 
            list(N = length(time),
                 id = nid <- as.numeric(as.factor(id)),
                 J = max(nid),
                 Nt = max(time),
                 time = time,
                 Tc = Tc,
                 Ic = Ic,
                 Vc = Vc,
                 s = 1,
                 C = 1))
                 

# row_vector<lower=0>[J] T1; // uninfected cells 
# row_vector<lower=0>[J] I1; // infected cells
# row_vector<lower=0>[J] V1; // viruses
# real <lower=0> lambda; // autonomous cell generation rate 
# real <lower=0> dT;     // death rate of uninfected cells
# real <lower=0> dI;     // death rate of infected cells
# real <lower=0> beta;   // infection rate 
# real <lower=0> p;      // production rate of viruses per infected cell
# real <lower=0> c;      // 
  
dat$J
inits <- list(T1 = rep(2000,dat$J), 
              I1 = rep(1, dat$J),
              V1 = rep(1, dat$J),
              lambda = 10,
              dT=.0001,
              dI=.0001,
              beta = .0000001,
              p = .0001,
              c = .0001)
initslist <- list(inits,inits,inits,inits)

system.time(
  mod <- sampling(cell_dso, dat, chains=4,
                  init = initslist)
) # 23871 sec = 6.63 hours = overnight

pars <- grepv('^T|^I|^V',names(mod), invert = T)
traceplot(mod, pars = pars)
pairs(mod, pars = pars)
save(mod, file = 'mod.rda')
