#######################################################################
# load packages
# 
library(fitdistrplus)
library(abind)
library(parallel)
library(extraDistr)

# set working directory to source file location
# 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# preliminaries

load("NY census 2019")  

counties <- colnames(age)

load("relative infectivity by age")

vax_profile <- read.csv("age vaccination profiles.csv")

time_dist <- function(duration, weights, n = 2e5){ # length(weights) = number of stages
  
  k<-length(weights)
  times <- NULL
  
  r <- matrix( rexp(k*n, rate=k/duration), nrow=n)
  t <- runif(n, min=0, max=rowSums(r))
  x <- runif(n)
  
  for(i in 1:n){
    
    breaks <- rep(NA,k-1); for(j in 1:(k-1)){breaks[j] <- sum(r[i,1:j])}
    
    interval <- if(t[i]>max(breaks)){k}else{min(which(t[i]<breaks))}
    
    if(x[i]<=weights[interval]){times<-c(times,t[i])}else{}
  }
  
  pars <- fitdist(times,"gamma", method="mle")$estimate
  
  return(pars)
}

beta_pars <- function(m,c){
  a <- (c/m)/((1/m)-1)
  b<-a*(1/m - 1)
  return(c(a,b))
}

immune_stage_dist<-function(t,d,n=5){ # t = time since vaccination; d = duration of each immune stage; n = number of stages
  dist<-rep(0,n)
  s<-min(1+floor(t/d),n)
  if(s>=1){dist[s]<-1}
  return(dist)
}

wane <- function(effect, residual, nsteps=5, shape=5){ 
  # takes vaccine effect in first immune stage, residual effect in last immune stage, number of steps and shape parameter and gives value for all steps, including first and last
  x_initial <- 1-effect
  x_final <- 1-residual*effect
  x <- x_final - (x_final-x_initial)*((nsteps-(1:nsteps))/(nsteps-1))^shape #### Changed from x <- x_final - (x_final-x_initial)*((nsteps-(1:nsteps)+1)/nsteps)^shape 
  return(1-x)
}

#######################################################################

generate_priors <- function(specify_counties, # c("Rockland", "Orange", "Sullivan", "Bronx", "Brooklyn", "Manhattan", "Queens", "Staten"
                            coverage_1, # age 2 (IPV)
                            coverage_2, # age 5 (IPV), pre-mandate
                            coverage_3, # age 5 (IPV), post-mandate
                            coverage_4, # age 2 (OPV)
                            N = 5e3, # size drawn from priors on Kalkowska et al parameters
                            M = 1e3 # # size drawn to estimate generation interval distribution parameters
                            ){
set.seed(1)

n.cores<-detectCores()

# AGE DISTRIBUTION

age_totals <- rowSums(matrix(age[,specify_counties],ncol=length(specify_counties))) # total in each age group

age_dist <- age_totals/sum(age_totals) # age distribution for the counties specified

# PARAMETER VALUES FROM KALKOWSKA ET AL.

immune_stage_duration <- 1
fecal_shed_duration <- 27.8
oral_shed_duration <- 13.4
latent_period <- 3
prop_oral <- 0.83

shedding_weights <- c(12,40,12,4)/68

# FIT DISTRIBUTIONS FOR TIME BETWEEN ONSET OF INFECTIOUSNESS AND TRANSMISSION (FOR SUSCEPTIBLE)

pars_fecal <- time_dist(fecal_shed_duration, shedding_weights) 

pars_oral <- time_dist(oral_shed_duration, shedding_weights) 

pars_latent <- c(shape = 2, rate = 2/3) # sum of two Exp(2/3) RVs

# initial (first immune stage) reductions in susceptibility, infectiousness, and duration of shedding relative to unvaccinated)

vax_effect <- matrix(c(0.28,0.78, 
                       0.41,0.77,
                       0.88,0.79,
                       0.36,0.68,
                       0.54,0.72),
                     ncol=2,byrow=T,dimnames=list(c("susceptibility","infectiousness_fecal","infectiousness_oral","shed_duration_fecal","shed_duration_oral"),c("IPV","OPV"))) 

# fraction of initial vaccine effect that remains in the last immune stage))

vax_residual_effect <- matrix(c(1e-10,0.38, #### Updated to 1e-10 from 0.02
                                0.37,0.91,
                                0.99,0.89,
                                0.42,0.87,
                                0.93,0.97),ncol=2,byrow=T,dimnames=dimnames(vax_effect))

# ----- drawing from priors on input parameters (vaccination coverage, shedding duration, waning rate, etc.)
 
cover_pars <- sapply(list(coverage_1,coverage_2,coverage_3, coverage_4),FUN=beta_pars, c=4)
cover_1 <- rbeta(N, cover_pars[1,1], cover_pars[2,1])
cover_2 <- mapply(rbeta(N, cover_pars[1,2], cover_pars[2,2]),FUN=max,cover_1) #### max(random proportion for IPV at 5, random proportion for IPV at 2)
cover_3 <- mapply(rbeta(N, cover_pars[1,3], cover_pars[2,3]),FUN=max,cover_2)
cover_4 <- rbeta(N, cover_pars[1,4], cover_pars[2,4])
 
immune_duration <- rgamma(N, 1.5*immune_stage_duration, 1.5)
 
duration_fecal <- rtnorm(N, mean=fecal_shed_duration, sd=fecal_shed_duration/6, a=0, b=Inf)
 
duration_oral <- rtnorm(N, mean=oral_shed_duration, sd=oral_shed_duration/6, a=0, b=Inf)
 
latent <- rgamma(N, 1.5*latent_period, 1.5)
 
prop_oral_pars <- beta_pars(m=prop_oral, c=2)
p_oral <- rbeta(N, prop_oral_pars[1], prop_oral_pars[2])
 
vax_effect_pars <- sapply(vax_effect, FUN=beta_pars, c=2)
vax_effect <- array(rbeta(10*N, vax_effect_pars[1,], vax_effect_pars[2,]), dim=c(5,2,N))
 
vax_resid_effect_pars <- sapply(vax_residual_effect, FUN=beta_pars, c=2)
vax_resid_effect <- array(rbeta(10*N, vax_resid_effect_pars[1,], vax_resid_effect_pars[2,]), dim=c(5,2,N))
 
# --- sorting population by age, vaccine received, and waning stage
 
ages <- 0:99

immune_dist_birth <- array(NA, dim=c(N,5,length(ages)))
immune_dist_2 <- array(NA, dim=c(N,5,length(ages)))
immune_dist_5 <- array(NA, dim=c(N,5,length(ages)))
immune_dist_2019 <- array(NA, dim=c(N,5,length(ages)))

for(i in 1:N){
  immune_dist_birth[i,,] <- sapply(ages, immune_stage_dist, d=immune_duration[i])
  immune_dist_2[i,,] <- sapply(ages-2, immune_stage_dist, d=immune_duration[i])
  immune_dist_5[i,,] <- sapply(ages-5, immune_stage_dist, d=immune_duration[i])
  immune_dist_2019[i,,] <- matrix(rep(immune_stage_dist(3, d=immune_duration[i]), length(ages)), nrow=5)
}

LPV_birth <- sweep(immune_dist_birth, MARGIN=3, STATS=vax_profile$LPV.at.birth,FUN="*")
OPV_cov_2 <- sweep(cover_4*immune_dist_2,MARGIN=3,STATS=vax_profile$OPV.at.age.2,FUN="*")
OPV_cov_5 <- sweep(cover_4*immune_dist_2,MARGIN=3,STATS=vax_profile$OPV.at.age.5,FUN="*") #### Added
IPV_cov_2 <- sweep(cover_1*immune_dist_2,MARGIN=3,STATS=vax_profile$IPV.at.age.2,FUN="*")
IPV_cov_5_pre <- sweep(cover_2*immune_dist_5,MARGIN=3,STATS=vax_profile$IPV.at.age.5.premandate,FUN="*") #### Changed from sweep((cover_2-cover_1)*immune_dist_5,MARGIN=3,STATS=vax_profile$IPV.at.age.5.premandate,FUN="*")      #### Pre-mandate IPV coverage at 5 for those with no IPV at 2
IPV_cov_5_post <- sweep(cover_3*immune_dist_5,MARGIN=3,STATS=vax_profile$IPV.at.age.5.postmandate,FUN="*") #### Changed from sweep((cover_3-cover_1)*immune_dist_5,MARGIN=3,STATS=vax_profile$IPV.at.age.5.postmandate,FUN="*")
IPV_cov_2019 <- sweep(cover_3*immune_dist_2019,MARGIN=3,STATS=vax_profile$IPV.in.2019,FUN="*") #### Changed from sweep((cover_3-cover_2)*immune_dist_2019,MARGIN=3,STATS=vax_profile$IPV.in.2019,FUN="*")

OPV <- LPV_birth + OPV_cov_2 + OPV_cov_5 

IPV <- IPV_cov_2 + IPV_cov_5_pre + IPV_cov_5_post + IPV_cov_2019    

IPV[1,,]

colSums(OPV[1,,]+IPV[1,,])

stage_age_vax <- abind(IPV, OPV,along=4)

stage_age_vax[1,,,1]+stage_age_vax[1,,,2]

stage_age_vax <- sweep(stage_age_vax,MARGIN=3,STATS=age_dist,FUN="*")

dimnames(stage_age_vax) <- list(NULL, NULL, NULL, c("IPV", "OPV"))

# --- stage-specific vaccine effectiveness

vax_effect_by.stage <- array(mapply(vax_effect,FUN=wane, residual=vax_resid_effect),dim=c(5,5,2,N), dimnames = list(NULL, c("susceptibility","infectiousness_fecal","infectiousness_oral","shed_duration_fecal","shed_duration_oral"), c("IPV","OPV"), NULL))

vax_effect_by.stage <- aperm(vax_effect_by.stage, perm=c(4,1,2,3))

# --- susceptible, partially immune, and refractory proportions of population

susceptible <- 1 - apply(stage_age_vax, MARGIN = 1, FUN = sum)

stage_vax <- apply(stage_age_vax, MARGIN = c(1,2,4), FUN = sum)

ve_susc_stage <- vax_effect_by.stage[,,"susceptibility",]

refractory <- apply(stage_vax*ve_susc_stage, MARGIN = 1, FUN = sum)

partial <- 1 - susceptible - refractory

# --- susceptible generation interval

# shedding duration specified by Kalkowska et al: fecal_shed_duration, oral_shed_duration
# gamma dist. parameters for time-to-transmission (from onset of infectiousness) fit to the above durations: pars_fecal, pars_oral
# shedding durations drawn from priors: duration_fecal, duration_oral

ratio_duration_fecal <- duration_fecal/fecal_shed_duration
# value drawn from prior divided by original (point) value
ratio_duration_oral <- duration_oral/oral_shed_duration # similar
ratio_latent <- latent/latent_period

pars.fecal <- cbind(rep(pars_fecal[1],N), pars_fecal[2]/ratio_duration_fecal) #### rate par. equals to (pars_fecal[2]/duration_fecal)*fecal_shed_duration
pars.oral <- cbind(rep(pars_oral[1], N), pars_oral[2]/ratio_duration_oral)
pars.latent <- cbind(rep(pars_latent[1],N), pars_latent[2]/ratio_latent)

n_oral <- rbinom(N, size = M, prob=p_oral)
n_fecal <- M - n_oral

n_fecal_oral <- rbind(n_fecal, n_oral)
shape_fecal_oral <- rbind(pars.fecal[,1], pars.oral[,1])
rate_fecal_oral <- rbind(pars.fecal[,2], pars.oral[,2])

s_gen_time <- array(NA, dim=c(M,N))
for(i in 1:N){
  s_gen_time[,i] <- rgamma(M, shape=rep(shape_fecal_oral[,i],n_fecal_oral[,i]), rate=rep(rate_fecal_oral[,i],n_fecal_oral[,i])) + rgamma(M, shape=pars.latent[i,1], rate=pars.latent[i,2])
}

list_s_gen_time <- lapply(asplit(s_gen_time, MARGIN=2), FUN=c)

clust <- makeCluster(n.cores)

s_gen_pars <- parLapply(clust, X=list_s_gen_time, fun=fitdist, distr="gamma", method="mle")

s_gen_pars_matrix <- matrix(sapply(s_gen_pars,FUN=function(x){unname(x$estimate)}), nrow=N, ncol=2, byrow=T)

# --- relative infectiousness - susceptible

age_summed <- apply(stage_age_vax, MARGIN = c(1,3), FUN = sum) 

susc_age <- sweep(-age_summed, MARGIN=2, STATS=age_dist, FUN ="+") # subtracting vaccinated from the age-stratified population 
#### to avoid negative values in `susc_age`, suggest `susc_age[which(susc_age<0)]=0`
susc_infectivity <- rowSums(sweep(susc_age, MARGIN=2, STATS=infectivity_age, FUN="*")) 

s_infectiousness <- susc_infectivity/susceptible

# --- relative infectiousness - partially immune

p_stage_age_vax <- sweep(stage_age_vax, MARGIN=c(1,2,4), STATS= (1 - vax_effect_by.stage[,,"susceptibility",]), FUN="*")/partial # sort partially immune population by age, vaccine received, and waning stage

relative_duration <- 1 - vax_effect_by.stage[,,c("shed_duration_fecal", "shed_duration_oral"),]
# relative duration of shedding for each waning stage, vaccine, and route of transmission (measured against the duration of shedding by the same route in susceptible individuals)

relative_infectiousness <- 1 - vax_effect_by.stage[,,c("infectiousness_fecal", "infectiousness_oral"),]
# relative infectiousness for each waning stage, vaccine, and route of transmission (measured against susceptible)

relative_total_shedding <- relative_duration*relative_infectiousness
dimnames(relative_total_shedding) <- list(NULL, NULL, c("fecal", "oral"), c("IPV", "OPV"))
# relative total shedding for each waning stage, vaccine and route of transmission (measured against transmission via the same route by susceptible individuals)

relative_transmission_prob_by.route <- array(NA, dim = c(N,5,2,2), dimnames=list(NULL, NULL, c("IPV", "OPV"), c("fecal", "oral")))

relative_transmission_prob_by.route[,,,1] <- sweep(relative_total_shedding[,,"fecal",],MARGIN=1, STATS=(1-p_oral),FUN="*")

relative_transmission_prob_by.route[,,,2] <- sweep(relative_total_shedding[,,"oral",],MARGIN=1, STATS=p_oral,FUN="*")

relative_transmission_prob <- apply(relative_transmission_prob_by.route, MARGIN=1:3, FUN=sum)

relative_transmission_prob_age <- relative_transmission_prob%o%infectivity_age
relative_transmission_prob_age <- aperm(relative_transmission_prob_age, perm=c(1,2,4,3))

p_infectiousness <- apply(relative_transmission_prob_age*p_stage_age_vax, MARGIN=1, FUN=sum)/s_infectiousness

# --- partially immune generation interval

shape_params <- array(NA, dim = c(N,5,2,2), dimnames = list(NULL, NULL, c("IPV", "OPV"), c("fecal", "oral")))
rate_params <- array(NA, dim = c(N,5,2,2), dimnames = list(NULL, NULL, c("IPV", "OPV"), c("fecal", "oral")))

shape_params[,,,1] <- pars_fecal[1]
shape_params[,,,2] <- pars_oral[1]

rate_params[,,,1] <- pars_fecal[2]/relative_duration[,,"shed_duration_fecal",]
rate_params[,,,2] <- pars_oral[2]/relative_duration[,,"shed_duration_oral",]

shape_params_rearr <- aperm(shape_params, perm=c(2:4,1))
rate_params_rearr <- aperm(rate_params, perm=c(2:4,1))
# dims = stage, IPV/OPV, fecal/oral, 1:N

p_stage_vax <- apply(p_stage_age_vax, MARGIN=c(1,2,4), FUN=sum)

relative_infectiousness <- sweep(relative_transmission_prob_by.route, MARGIN=1:3, STATS=p_stage_vax, FUN="*")
# dims = 1:N, stage, IPV/OPV, fecal/oral

n_transmit <- mapply(rep(1, N), FUN=rmultinom, size = M, prob = asplit(relative_infectiousness, MARGIN = 1))
# dims = 1:20, 1:N
 
p_gen_time <- array(NA,dim=c(M,N))
  
for(i in 1:N){
  p_gen_time[,i] <- rgamma(M, shape=rep(shape_params_rearr[,,,i], n_transmit[,i]), rate=rep(rate_params_rearr[,,,i], n_transmit[,i])) + rgamma(M, shape=pars.latent[i,1], rate=pars.latent[i,2])
} 

list_p_gen_time <- lapply(asplit(p_gen_time, MARGIN=2), FUN=c)

clust <- makeCluster(n.cores)

p_gen_pars <- parLapply(clust, X=list_p_gen_time, fun=fitdist, distr="gamma", method="mle")

p_gen_pars_matrix <- matrix(sapply(p_gen_pars,FUN=function(x){unname(x$estimate)}), nrow=N, ncol=2, byrow=T)

out <- cbind(susceptible, refractory, p_infectiousness, s_gen_pars_matrix, p_gen_pars_matrix)

return(out)

}


