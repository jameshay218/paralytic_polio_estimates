library(epidemia)
library(lubridate)
library(rstanarm)
library(dplyr)
options(mc.cores = parallel::detectCores())

paralysis <- c(rep(0, 48),1, rep(0,21))
times <- as.Date(seq_along(paralysis),origin="2022-06-01")

data <- data.frame(location="NY",paralysis=paralysis,date=times,pop=1e6)

si <- extraDistr::ddgamma(0:25,5,1)
si <- si/sum(si)
plot(si)
incu_period <- extraDistr::ddgamma(0:25,5,1)
incu_period <- incu_period/sum(incu_period)
rt <- epirt(formula = R(location,date) ~ 1 + rw(time = date, prior_scale=0.05))

inf <- epiinf(gen = si, seed_days=20L, pop_adjust = TRUE, pops = pop,latent = TRUE,
              prior_susc = normal(0.5, 0.1), prior_seeds = normal(5, 1))
length(times)
obs <-  epiobs(formula = paralysis ~ 1, link = "identity",
               prior_intercept = normal((1/2000),0.00025),i2o = incu_period)

fm <- epim(rt = rt, inf = inf, obs = obs, data = data,
           iter = 4e3, control = list(max_treedepth = 12), seed=12345)
fm$stan_summary %>% tail
spaghetti_rt(fm)
spaghetti_infections(fm)
spaghetti_obs(fm,type = "paralysis")
rstan::traceplot(fm$stanfit,pars="paralysis|(Intercept)")
rstan::stan_dens(fm$stanfit,pars="paralysis|(Intercept)")
