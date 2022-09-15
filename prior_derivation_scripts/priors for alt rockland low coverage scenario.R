setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("Generate priors.R")

############################

rockland_low <- generate_priors(specify_counties = "Rockland",
                                    coverage_1 = 0.40, 
                                    coverage_2 = 0.9, 
                                    coverage_3 = 0.94,
                                    coverage_4 = 0.9)

pars_rockland_low <- data.frame(matrix(nrow=7,ncol=4,dimnames=list(NULL,c("model parameter","distribution","par1","par2"))))

pars_rockland_low[,1] <- c("prop_susceptible","prop_refractory","relative_infectiousness","susceptible_generation_interval_shape","susceptible_generation_interval_rate","partial_generation_interval_shape","partial_generation_interval_rate")

x<-rockland_low[,6] 
hist(x,breaks=100,freq=F)
fit <- fitdist(x,"norm","mge",gof="CvM")$estimate
curve(dnorm(x, fit[1],fit[2]),from=min(x),to=max(x),add=T, lwd=2)

pars_rockland_low[,2] <- c("beta", #1
                           "beta", #2
                           "beta", #3
                           "norm", #4
                           "beta", #5
                           "norm", #6
                           "beta" #7
)

for(i in 1:7){pars_rockland_low[i,3:4] <- unname(fitdist(rockland_low[,i], distr=pars_rockland_low[i,2], method="mge", gof="CvM")$estimate)}

pars_rockland_low[,3]/(pars_rockland_low[,3]+pars_rockland_low[,4])

scenario <- rep("rockland_low_coverage", 7)

priors <- cbind(scenario,pars_rockland_low)

priors

write.csv(priors, file = "priors - alt low coverage scenario.csv",row.names=F)
