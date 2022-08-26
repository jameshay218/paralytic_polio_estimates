library(Hmisc)
library(ggplot2)
library(patchwork)
library(dplyr)
library(paletteer)
setwd("~/Documents/GitHub/paralytic_polio_estimates/")
source("simulation_functions_twoimmune.R")

setwd("~/Documents/GitHub/paralytic_polio_estimates/sims/")
all_res <-  NULL
for(i in 1:1000){
    if(file.exists(paste0("simulation_",i,".RData"))){
        load(paste0("simulation_",i,".RData"))
        all_res[[i]] <- pars
    }
}
res <- as_tibble(do.call("bind_rows",all_res))
res <- res %>% select(-c(Rt, inc_s,inc_ps,para,para_s,para_ps))

setwd("~/Documents/GitHub/paralytic_polio_estimates/sims_traj/")
all_traj <-  NULL
for(i in 1:1000){
    if(file.exists(paste0("traj_",i,".RData"))){
        load(paste0("traj_",i,".RData"))
        all_traj[[i]] <- res2
    }
}
traj <- as_tibble(do.call("bind_rows",all_traj))

setwd("~/Documents/GitHub/paralytic_polio_estimates/sims_nyc/")
all_nyc <-  NULL
for(i in 1:1000){
    if(file.exists(paste0("nyc_",i,".RData"))){
        load(paste0("nyc_",i,".RData"))
        all_nyc[[i]] <- nyc
    }
}
traj_nyc <- as_tibble(do.call("bind_rows",all_nyc))


traj <- traj %>% left_join(res %>% select(sim, date_start)) %>%
    mutate(t = as.Date(date_start + t)) %>%
    drop_na() 

traj_nyc <- traj_nyc %>% left_join(res %>% select(sim, date_start)) %>%
    mutate(t = as.Date(date_start + t)) %>%
    drop_na() 

traj$vacc_prop <- as.factor(traj$vacc_prop)
trajectories <- traj %>% filter(vacc_prop==1) %>% dplyr::select(-c(vacc_prop,rep))
save(trajectories,file="../outputs/fitted_trajectories.RData")
save(traj_nyc,file="../outputs/nyc_trajectories.RData")
save(res,file="../outputs/filtered_parameters.RData")

#res <- res %>% filter(date_start > "2022-05-01")

samps <- sample(unique(res$sim),1000,replace=TRUE)
samps1 <- sample(unique(res$sim),100,replace=TRUE)
ggplot(traj %>% filter(vacc_prop == 1) %>%
           filter(t <= "2022-08-01") %>%
           filter(sim %in% samps1) %>%
           group_by(sim) %>%
           filter(t != min(t)) %>%
           left_join(res %>%  filter(sim %in% samps1))) + 
    geom_line(aes(x=t,y=inc,group=sim,col=Re>1),size=0.25,alpha=0.5) +
    geom_vline(xintercept=as.Date("2022-07-18"))

ggplot(traj_nyc %>% 
           filter(t <= "2023-08-01") %>%
           filter(sim %in% samps1)) + 
    geom_line(aes(x=t,y=inc,group=sim),size=0.25,alpha=0.5) +
    geom_vline(xintercept=as.Date("2022-07-18"))

ggplot(traj_nyc %>% 
           filter(t <= "2023-08-01") %>%
           filter(sim %in% samps1)) + 
    geom_line(aes(x=t,y=para,group=sim),size=0.25,alpha=0.5) +
    geom_vline(xintercept=as.Date("2022-07-18"))

ggplot(traj %>% filter(vacc_prop == 1) %>%
           filter(sim %in% samps)) + geom_line(aes(x=t,y=para,group=interaction(sim,rep)))

traj_prop <- traj %>% mutate(y=inc>0) %>% group_by(t,vacc_prop) %>% dplyr::summarize(y=sum(y),N=n()) %>% 
    group_by(t,vacc_prop) %>%
    mutate(prop=y/N,lower=binconf(y,N,0.05)[2], upper=binconf(y,N,0.05)[3]) 

traj_para_prop <- traj %>% mutate(y=para>0) %>% group_by(t,vacc_prop) %>% dplyr::summarize(y=sum(y),N=n()) %>% 
    group_by(t,vacc_prop) %>%
    mutate(prop=y/N,lower=binconf(y,N,0.05)[2], upper=binconf(y,N,0.05)[3])

p1 <- traj_prop %>% ggplot() + 
    geom_ribbon(aes(x=t,ymin=lower,ymax=upper,fill=vacc_prop),alpha=0.25) +
    geom_line(aes(x=t,y=prop,col=vacc_prop),size=1) +
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1)) +
    scale_x_date(limits=as.Date(c("2022-08-22","2023-07-01")),
                 breaks="1 month") +
    theme_classic() +
    ylab("Proportion of simulations with\n ongoing infections") +
    theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          panel.grid.major=element_line(size=0.25,color="grey70"),
          legend.position="top") +
    labs(tag="A",fill="Prop vaccinated",color="Prop vaccinated") +
    scale_color_viridis_d() + scale_fill_viridis_d()
p2 <- traj_para_prop %>% 
    ggplot() + 
    geom_ribbon(aes(x=t,ymin=lower,ymax=upper,fill=vacc_prop),alpha=0.25) +
    geom_line(aes(x=t,y=prop,col=vacc_prop),size=1) +
    scale_y_continuous(limits=c(0,0.1),breaks=seq(0,1,by=0.02))+
    scale_x_date(limits=as.Date(c("2022-08-22","2023-07-01")),
                 breaks="1 month")+
    theme_classic()+
    ylab("Proportion of simulations with\n ongoing cases of paralysis") +
    xlab("Date") + 
    theme(panel.grid.major=element_line(size=0.25,color="grey70"),
          axis.text.x=element_text(angle=45,hjust=1),legend.position="none")+
    labs(tag="B",fill="Prop vaccinated",color="Prop vaccinated") +
    scale_color_viridis_d() + scale_fill_viridis_d()

fig1 <- p1/p2

## Generate prior draws for density plots

priors <- simulate_priors(n=100000,
                        incu_mean_prior_mean=16,
                        incu_mean_prior_var=5,
                        incu_var_prior_mean=10,
                        incu_var_prior_var=3,
                        infect_mean_prior_mean=1.28/0.19,
                        infect_mean_prior_var=1,
                        infect_var_prior_mean=1.28/(0.19*0.19),
                        infect_var_prior_var=3,
                        latent_period=3,
                        prob_paralysis_mean=0.0005,
                        prob_paralysis_ps_mean = 0.01*0.0005,
                        prob_paralysis_var = 1e-8,
                        prob_paralysis_ps_var = 1e-11,
                        R0_par1=4.9,R0_par2=2,R0_dist="truncnorm",
                        prop_immune_pars = c(23.4,59.3,17.3),
                        rel_R0_mean = 0.18,rel_R0_var=0.001)

priors <- priors %>% mutate(Re = R0 * prop_immune_groups.3 + R0*rel_R0s*prop_immune_groups.2)
res <- res %>% mutate(Re = R0 * prop_immune_groups.3 + R0*rel_R0*prop_immune_groups.2)
p1 <- ggplot(res) +
    geom_density(data=priors,aes(x=R0,fill="Prior"),alpha=0.1) +
    geom_density(aes(x=R0,fill="Post-filter"),alpha=0.25) +
    theme_classic() +
    geom_vline(xintercept=1,col="red",linetype="dashed") +
    scale_fill_manual(name="",
                       values=c("Prior"="black","Post-filter"="blue")) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits=c(0,10),breaks=seq(0,10,by=1)) +
    xlab("Basic reproductive number in fully susceptible R0") +
    ylab("Density") +
    labs(tag="A") 

p2 <- ggplot(res) +
    geom_density(data=priors,aes(x=R0*rel_R0s,fill="Prior"),alpha=0.1) +
    geom_density(aes(x=R0*rel_R0,fill="Post-filter"),alpha=0.25) +
    theme_classic() +
    geom_vline(xintercept=1,col="red",linetype="dashed") +
    scale_fill_manual(name="",
                      values=c("Prior"="black","Post-filter"="blue")) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits=c(0,5),breaks=seq(0,10,by=1)) +
    xlab("Basic reproductive number in partially immune") +
    ylab("Density") +
    labs(tag="B") +
    theme(legend.position="none")
    
p3 <- ggplot(res) + 
    geom_density(data=priors,aes(x=Re,fill="Prior"),alpha=0.1) +
    geom_density(aes(x=Re,fill="Post-filter"),alpha=0.25) +
    theme_classic() +
    scale_fill_manual(name="",values=c("Prior"="black","Post-filter"="blue")) +
    geom_vline(xintercept=1,col="red",linetype="dashed") +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits=c(0,5),breaks=seq(0,5,by=1))+
    xlab("Effective reproductive number weighted by\n proportion partially immune")+
    ylab("Density")+
    labs(tag="C")+
    theme(legend.position="none")


p4 <- ggplot(res) + 
    geom_density(data=priors,aes(x=prop_immune_groups.1,fill="Prior"),alpha=0.1) +
    geom_density(aes(x=prop_immune_groups.1,fill="Post-filter"),alpha=0.25) +
    theme_classic() +
    scale_fill_manual(name="",values=c("Prior"="black","Post-filter"="blue")) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits=c(0,1)) +
    xlab("Proportion fully immune")+
    ylab("Density")+
    labs(tag="D")+
    theme(legend.position="none")

p5 <- ggplot(res) + 
    geom_density(data=priors,aes(x=prop_immune_groups.2,fill="Prior"),alpha=0.1) +
    geom_density(aes(x=prop_immune_groups.2,fill="Post-filter"),alpha=0.25) +
    theme_classic() +
    scale_fill_manual(name="",values=c("Prior"="black","Post-filter"="blue")) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits=c(0,1)) +
    xlab("Proportion partially immune")+
    ylab("Density")+
    labs(tag="E")+
    theme(legend.position="none")

p6 <- ggplot(res) + 
    geom_density(data=priors,aes(x=prop_immune_groups.3,fill="Prior"),alpha=0.1) +
    geom_density(aes(x=prop_immune_groups.3,fill="Post-filter"),alpha=0.25) +
    theme_classic() +
    scale_fill_manual(name="",values=c("Prior"="black","Post-filter"="blue")) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits=c(0,1)) +
    xlab("Proportion fully susceptible")+
    ylab("Density")+
    labs(tag="F")+
    theme(legend.position="none")

p7 <- ggplot(res) + 
    geom_density(data=priors,aes(x=prob_paralysis_s,fill="Prior"),alpha=0.1) +
    geom_density(aes(x=prob_paralysis_s,fill="Post-filter"),alpha=0.25) +
    theme_classic() +
    scale_fill_manual(name="",values=c("Prior"="black","Post-filter"="blue")) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits=c(0,0.001))+
    xlab("Infection-paralysis-ratio in susceptibles")+
    ylab("Density")+
    labs(tag="G")+
    theme(legend.position="none")

p8 <- ggplot(res) + 
    geom_density(data=priors,aes(x=prob_paralysis_ps,fill="Prior"),alpha=0.1) +
    geom_density(aes(x=prob_paralysis_ps,fill="Post-filter"),alpha=0.25) +
    theme_classic() +
    scale_fill_manual(name="",values=c("Prior"="black","Post-filter"="blue")) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits=c(0,0.00005))+
    xlab("Infection-paralysis-ratio in partially immune")+
    ylab("Density")+
    labs(tag="H")+
    theme(legend.position="none")


fig2 <- (p1+theme(legend.position="none")) + p2 + p3 + p4 + p5 + p6 + p7 + p8 + ggpubr::get_legend(p1) + plot_layout(ncol=3)

## Proportion with Rt>1 by 2022-08-20
y <- traj %>% 
    ungroup() %>%
    filter(t == "2022-08-20", vacc_prop==1) %>% 
    mutate(y=Rt>1) %>%
    pull(y)
sum(y)/nrow(res)

## Proportion with Incidence > 0 by 2022-08-20
y <- traj %>% 
    ungroup() %>%
    filter(t == "2022-08-20", vacc_prop==1) %>% 
    mutate(y=inc>0) %>%
    pull(y)
sum(y)/nrow(res)

## Some estimates
## Paralysis cases by end of outbreak
table1 <- traj %>% 
    group_by(sim, vacc_prop) %>% 
    dplyr::summarize(y=sum(para)) %>%
    group_by(vacc_prop) %>%
    dplyr::summarize(mean_para=mean(y),
              median_para=median(y),
              lower95=quantile(y,0.025), lower80=quantile(y,0.1),
              upper80=quantile(y,0.9),upper95=quantile(y,0.975))

## Proportion of outbreaks with more than the initial case of paralysis
table2 <- traj %>% 
    group_by(sim, vacc_prop) %>% 
    dplyr::summarize(y=sum(para)) %>%
    ungroup() %>%
    mutate(x=y>1) %>%
    group_by(vacc_prop) %>%
    dplyr::summarize(y=sum(x)/n())


## Number of infections by end of simulation
table3 <- traj %>% 
    group_by(sim, vacc_prop) %>% 
    dplyr::summarize(y=sum(inc)) %>%
    group_by(vacc_prop) %>%
    dplyr::summarize(mean_inf=mean(y),
                     median_inf=median(y),
                     lower95=quantile(y,0.025), lower80=quantile(y,0.1),
                     upper80=quantile(y,0.9),upper95=quantile(y,0.975))

## Proportion of outbreaks with more than 100 further infections
table4 <- traj %>% 
    group_by(sim, vacc_prop) %>% 
    dplyr::summarize(y=sum(inc)) %>%
    ungroup() %>%
    mutate(x=y>100) %>%
    group_by(vacc_prop) %>%
    dplyr::summarize(y=sum(x)/n())


setwd("~/Documents/GitHub/paralytic_polio_estimates/")

write_csv(table1,"figures/table1.csv")
write_csv(table2,"figures/table2.csv")
write_csv(table3,"figures/table3.csv")
write_csv(table4,"figures/table4.csv")
ggplot2::ggsave("figures/fig1.pdf",fig1,width=7,height=8)
ggplot2::ggsave("figures/fig2.pdf",fig2,width=7,height=6)


## Generation interval plots
incubation_period <- function(t, incu_scale, incu_shape, max_incu_period=50){extraDistr::ddgamma(t,scale=incu_scale,shape=incu_shape)/extraDistr::pdgamma(max_incu_period,scale=incu_scale,shape=incu_shape)}
infectiousness <- function(t,latent_period=3,infect_scale,infect_shape,max_infectious_period=50){extraDistr::ddgamma(t-latent_period, scale=infect_scale, shape=infect_shape)/extraDistr::pdgamma(max_infectious_period, shape=infect_shape,scale=infect_scale)}
#infectiousness_ps <- function(t){c(rep(0, latent_period), extraDistr::ddgamma(t, scale=infect_scale_ps, shape=infect_shape_ps))}
infectiousness_ps <- function(t, latent_period=3,infect_ps_par,max_infectious_period=50){dexp(t-latent_period, infect_ps_par)/pexp(max_infectious_period, infect_ps_par)}

ts <- 0:25
res1 <- res %>% filter(sim %in% samps)
incu_periods_dat <- matrix(0,nrow=nrow(res1),ncol=length(ts))
infectiousness_dat <- matrix(0,nrow=nrow(res1),ncol=length(ts))
infectiousness_ps_dat <- matrix(0,nrow=nrow(res1),ncol=length(ts))
for(x in 1:nrow(res1)){
    incu_periods_dat[x,] <- incubation_period(ts, res1$incu_scale[x],res1$incu_shape[x])
    infectiousness_dat[x,] <- infectiousness(ts, res1$latent_period[x], infect_scale=res1$infect_scale[x],infect_shape=res1$infect_shape[x])
    infectiousness_ps_dat[x,] <- incubation_period(ts, res1$latent_period[x], res1$infect_ps_par[x])
}
incu_periods_dat <- reshape2::melt(incu_periods_dat)
colnames(incu_periods_dat) <- c("sim","time_since_infection","probability")
incu_periods_dat$model <- "Paralysis incubation period"
infectiousness_dat <- reshape2::melt(infectiousness_dat)
colnames(infectiousness_dat) <- c("sim","time_since_infection","probability")
infectiousness_dat$model <- "Susceptible generation interval"

infectiousness_ps_dat <- reshape2::melt(infectiousness_ps_dat)
colnames(infectiousness_ps_dat) <- c("sim","time_since_infection","probability")
infectiousness_ps_dat$model <- "Partially immune generation interval"
dists <- bind_rows(incu_periods_dat,infectiousness_dat,infectiousness_ps_dat)
dists_summary <- dists %>% group_by(model, time_since_infection)
ggplot(dists) + geom_line(aes(x=time_since_infection,y=probability,group=sim),alpha=0.1) +
    facet_wrap(~model,scales="free_y",ncol=1)
