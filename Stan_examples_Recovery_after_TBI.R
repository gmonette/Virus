#' ---
#' title: "R Stan examples: IQ recovery after TBI"
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
knitr::opts_chunk$set(comment = '')
#'
#' ## Stan for Asymptotic Recovery Curves
#' 
#' We can use non-linear mixed effects models for normally
#' distributed responses to estimate asymptotic recovery
#' curves for VIQ and PIQ separately. It's feasible
#' to extend these models to handle VIQ and PIQ as a 
#' multivariate response, but doing so is not straighforward.
#' 
#' Here, we use Stan to build a model for a multivariate
#' response which will allow us to  compare the  parameters
#' for the two models in a way that exploits the covariance
#' between the responses.
#' 
#' We start with univariate models which we then combine
#' into a multivariate model.
#' 
#+ include=FALSE
# note that include=FALSE will evaluate but not show code
interactive <- FALSE  # do not run interactive code
knitr::opts_chunk$set(comment='')
if(.Platform$OS.type == 'windows') windowsFonts(Arial=windowsFont("TT Arial")) 
#+ include=FALSE, eval=FALSE
# we don't want to evaluate this when running the code for Markdown
interactive <- TRUE  # run this when running code by hand
#'
#' Loading packages:
#' 
#+ packages
library(spida2)
library(magrittr, pos = 1000) 
library(lattice)
library(latticeExtra)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#'
#' ## Explore data
#'

data(iq)
?iq
head(iq)
names(iq) <- tolower(names(iq))

dd <-iq

(p <- xyplot(piq ~ dayspc | sex, dd, groups = id, type = 'b'))
update(p, xlim = c(0,4000))
if(interactive) {
  rows <- numeric(0)
  # can repeat:
  trellis.focus()
  rows <- c(rows, panel.identify(labels=dd$id))
  # end
  trellis.unfocus()
  rows
  save(rows, file = 'rows.rda')
}

load('rows.rda', verbose = T)
dd[rows,] %>% sortdf(~dcoma+dayspc)
# id = 2600 retested 4 days apart
dd[dd$id %in% dd[rows,'id'],]  %>%  sortdf(~dcoma+dayspc)

if(interactive){
  library(p3d)
  Init3d()
  dd$dcoma.cat <- cut(dd$dcoma, c(-1,2,5,10,20,50,Inf))
  Plot3d( viq ~ piq + log(dayspc) |
            dcoma.cat, dd, groups = id,
          col = heat.colors(6))
  Plot3d( viq ~ log(dcoma+2) + log(dayspc) |
            dcoma.cat, dd, groups = id,
          col = heat.colors(6))
  Plot3d( piq ~ log(dcoma+2) + log(dayspc) |
            dcoma.cat, dd, groups = id,
          col = heat.colors(12)[1:6])
  fg()
  Id3d()
}
#'
#' ## Stan program
#'
stanfile <- 'asymp_model_1.stan' 
pr <- function(x) print(readLines(x), quote = FALSE) # write a function for repetitive tasks
pr(stanfile)
  system.time(
  asymp_model_dso <- stan_model(stanfile)
)
#'
#' ### Prepare data list
#'
names(dd)
pred <- expand.grid(time=seq(0, 3*365, 30), 
                    coma = seq(0,30,5))
dat <- list(
  N = nrow(dd),
  id = nid <- as.numeric(as.factor(dd$id)),
  J = max(nid),
  y = dd$piq,
  time = dd$dayspc,
  coma = sqrt(dd$dcoma),
  Np = nrow(pred),
  time_p = pred$time,
  coma_p = sqrt(pred$coma)
)
set.seed(789723)
mod <- sampling(asymp_model_dso, dat, chains = 4, iter = 2000)
# Note: chains and iter are defaults
traceplot(mod)
names(mod)
pars <- grepv('^u|y_p', names(mod), invert = T)
pars
traceplot(mod, pars = pars)
pairs(mod, pars = pars)
if(interactive) {
  library(shinystan)
  mod_sso <- launch_shinystan(mod)
}
zd <- as.data.frame(mod, pars = pars)
head(zd)
xyplot(hrt~asymp, zd)
print(mod, pars = pars)
get_posterior_mean(mod, pars = pars)[,5]
#'
#' Note that you can do extensive plotting of 
#' the posterior sample with 'shinystan' including
#' some not highly satisfying 3d plots.
#'
#' ## Random initial deficit
#'
#' It wasn't feasible to fit a model that included
#' random initial deficits using 'lme' so this
#' will be an interesting experiment.
#'  
#' Note that the data will be identical.
#'
#'
stanfile <- 'asymp_model_2.stan' 
pr(stanfile)
system.time(
  asymp_model_2_dso <- stan_model(stanfile)
)
#'
#' Fitting the model. The data is the same.
#' 
mod2 <- sampling(asymp_model_2_dso, dat, chains = 4, iter = 2000)

pars <- grepv('^u|y_p', names(mod2), invert = T)
pars
traceplot(mod2, pars = pars)
pairs(mod2, pars = pars)
#'
#'
#' ### Poor convergence: possible remedies
#' 
#' The variability in the sd of initial deficit is
#' poorly estimated.
#' 
#' We can try to improve this in a least two ways:
#' 
#' 1.  Reparametrize the model: Taking the initial point
#'     at time 0 puts it outside the range of the data.
#'     Using a later time as time 0 could improve the 
#'     estimate.  This can be done just by changing the
#'     data and subtracting a constant from time. We then 
#'     need to remember the changed interpretation of 
#'     parameters although most parameters are not
#'     affected, e.g. 'asymp' and 'hrt'
#' 2.  Use an informative prior for 'sigma_idef'. So far
#'     all models have used flat priors for hyper-parameters.
#'     Gelman would not hesitate to use much more informative
#'     priors.
#'     
#' We'll try step 1 because it's easy. We'll subtract 180 
#' from time so that time 0 will be relocated to approximately
#' 6 months.
#' 
#' #### Reparametrize
#' 
dat_6 <- list(
  N = nrow(dd),
  id = nid <- as.numeric(as.factor(dd$id)),
  J = max(nid),
  y = dd$piq,
  time = dd$dayspc - 180,
  coma = sqrt(dd$dcoma),
  Np = nrow(pred),
  time_p = pred$time - 180,
  coma_p = sqrt(pred$coma)
)

mod2_6 <- sampling(asymp_model_2_dso, dat_6, chains = 4, iter = 2000)

pars <- grepv('^u|y_p', names(mod2_6), invert = T)
pars
traceplot(mod2_6, pars = pars)
pairs(mod2_6, pars = pars)
#'
#' #### Informative prior
#'
#' We will try a Gamma prior on 'sigma_idef'
#'
stanfile <- 'asymp_model_3.stan' 
pr(stanfile)
system.time(
  asymp_model_3_dso <- stan_model(stanfile)
)
#'
#' Fitting the model. The data is the same plus the
#' two parameters for the gamma prior.
#' 
mod3_6 <- sampling(asymp_model_3_dso, 
                 c(dat_6, gamma_df = 4, gamma_scale = .5), 
                 chains = 4, iter = 2000,
                 control = list(adapt_delta = 0.85))

pars <- grepv('^u|y_p', names(mod3_6), invert = T)
pars
traceplot(mod3_6, pars = pars)
pairs(mod3_6, pars = pars)
print(mod3_6, pars = pars)
get_posterior_mean(mod3_6, pars = pars)[,5] %>% cbind
get_posterior_mean(mod2_6, pars = pars)[,5] %>% cbind
get_posterior_mean(mod2, pars = pars)[,5] %>% cbind
#'
#' #### Baysian uniform prior 
#'
#' We'll just use our intuition about IQs. The standard deviation
#' couldn't 'possibly' be less than 3 and couldn't 'possibly'
#' be greater than 30!
#' 
stanfile <- 'asymp_model_3b.stan' 
pr(stanfile)
system.time(
  asymp_model_3b_dso <- stan_model(stanfile)
)
#'
#' Fitting the model. The data is the same plus the
#' two parameters for the gamma prior.
#' 
mod3b_6 <- sampling(asymp_model_3b_dso, 
                   dat_6,  
                   chains = 4, iter = 2000,
                   control = list(adapt_delta = 0.85))

pars <- grepv('^u|y_fit', names(mod3b_6), invert = T)
pars
traceplot(mod3b_6, pars = pars)
pairs(mod3b_6, pars = pars)
print(mod3b_6, pars = pars)
get_posterior_mean(mod3b_6, pars = pars)[,5] %>% cbind
get_posterior_mean(mod2_6, pars = pars)[,5] %>% cbind
get_posterior_mean(mod2, pars = pars)[,5] %>% cbind
#'
#' Alas, this didn't work either.
#' 
#' ### Sensitivity analysis
#' 
#' Since it looks like we can't estimate 'sigma_idef' from the data with
#' our model but we are concerned that it might have an uncertain
#' impact on our other results, we can resort to a sensitivity analysis.
#' Try a variety of plausible values and see how it affects other results.
#'

stanfile <- 'asymp_model_3c.stan' 
pr(stanfile)
system.time(
  asymp_model_3c_dso <- stan_model(stanfile)
)
#'
#' Fitting the model. The data is the same plus the
#' two parameters for the gamma prior.
#' 
mod3c_6 <- sampling(asymp_model_3c_dso, 
                    c(dat_6, sigma_idef = .1),  
                    chains = 4, iter = 2000)

pars <- grepv('^u|y_fit', names(mod3c_6), invert = T)
pars
traceplot(mod3c_6, pars = pars)
pairs(mod3c_6, pars = pars)
print(mod3c_6, pars = pars)
get_posterior_mean(mod3c_6, pars = pars)[,5] %>% cbind
get_posterior_mean(mod3c_6, pars = pars)[,5] %>% cbind

#'
#' ## Multivariate model for VIQ and PIQ
#'
#' Instead of simple standard deviations for each VIQ and PIQ, we need
#' covariance matrices to describe the within- and the between-subject
#' covariance between VIQ and PIQ.
#' 
#' The first model uses a uniform prior on covariance matrices using
#' constraints generated by Stan on covariance matrices.  
#' 
#' It is more efficient to use a combination of an LKJ prior on correlation
#' together with priors of variances. The second model below use this
#' form of parametrization.
#' 

head(dd)
nrow(dd)

tab(c(Tab(dd, ~id)))

stanfile <- 'asymp_model_4.stan' 
pr(stanfile)
system.time(
  asymp_model_4_dso <- stan_model(stanfile)
)
#'
#' ### Data for multivariate model
#'
#'

dat <- list(
  N = nrow(dd),
  id = nid <- as.numeric(as.factor(dd$id)),
  J = max(nid),
  iq = cbind(dd$viq,dd$piq),
  time = dd$dayspc,
  coma = sqrt(dd$dcoma)
)

mult_mod <- sampling( asymp_model_4_dso, dat)

pars <- grepv('^u' ,names(mult_mod), invert = T)
traceplot(mult_mod, pars = pars)
pairs(mult_mod, pars = pars)
print(mult_mod, pars = pars)

#'
#' ## Using a Cholesky parametrization
#'

