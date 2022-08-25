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

setwd("~/Documents/GitHub/paralytic_polio_estimates/sims_traj/")
all_traj <-  NULL
for(i in 1:1000){
    if(file.exists(paste0("traj_",i,".RData"))){
        load(paste0("traj_",i,".RData"))
        all_traj[[i]] <- res2
    }
}
traj <- as_tibble(do.call("bind_rows",all_traj))

traj <- traj %>% left_join(res %>% select(sim, date_start)) %>%
    mutate(t = as.Date(date_start + t)) %>%
    drop_na() 

traj$vacc_prop <- as.factor(traj$vacc_prop)

#res <- res %>% filter(date_start > "2022-05-01")

samps <- sample(unique(res$sim),1000,replace=TRUE)
samps1 <- sample(unique(res$sim),100,replace=TRUE)
ggplot(traj %>% filter(vacc_prop == 1) %>%
           filter(t <= "2022-08-01") %>%
           filter(sim %in% samps1) %>%
           group_by(sim) %>%
           filter(t != min(t)) %>%
           left_join(res)) + 
    geom_line(aes(x=t,y=inc,group=sim,col=Re>1),size=0.25,alpha=0.5) +
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
          legend.position=c(0.5,0.8),legend.direction = "horizontal") +
    labs(tag="A",fill="Prop vaccinated",color="Prop vaccinated") +
    scale_color_viridis_d() + scale_fill_viridis_d()
p2 <- traj_para_prop %>% 
    ggplot() + 
    geom_ribbon(aes(x=t,ymin=lower,ymax=upper,fill=vacc_prop),alpha=0.25) +
    geom_line(aes(x=t,y=prop,col=vacc_prop),size=1) +
    scale_y_continuous(limits=c(0,0.75),breaks=seq(0,1,by=0.1))+
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
                        R0_par1=0,R0_par2=10,
                        prop_immune_pars = c(23.4,59.3,17.3),
                        rel_R0_mean = 0.18,rel_R0_var=0.001)

p1 <- ggplot(res) +
    geom_density(data=priors,aes(x=R0,fill="Prior"),alpha=0.1) +
    geom_density(aes(x=R0,fill="Post-filter"),alpha=0.25) +
    theme_classic() +
    geom_vline(xintercept=1,col="red",linetype="dashed") +
    scale_fill_manual(name="",
                       values=c("Prior"="black","Post-filter"="blue")) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits=c(0,8),breaks=seq(0,8,by=1)) +
    xlab("Basic reproductive number R0") +
    ylab("Density") +
    labs(tag="A") +
    theme(legend.position=c(0.8,0.8))
    
p2 <- ggplot(res) + 
    geom_density(data=priors,aes(x=Re,fill="Prior"),alpha=0.1) +
    geom_density(aes(x=Re,fill="Post-filter"),alpha=0.25) +
    theme_classic() +
    scale_fill_manual(name="",values=c("Prior"="black","Post-filter"="blue")) +
    geom_vline(xintercept=1,col="red",linetype="dashed") +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits=c(0,8),breaks=seq(0,8,by=1))+
    xlab("Effective reproductive number Re")+
    ylab("Density")+
    labs(tag="B")+
    theme(legend.position="none")

p3 <- ggplot(res) + 
    geom_density(data=priors,aes(x=prob_paralysis,fill="Prior"),alpha=0.1) +
    geom_density(aes(x=prob_paralysis_s,fill="Post-filter"),alpha=0.25) +
    theme_classic() +
    scale_fill_manual(name="",values=c("Prior"="black","Post-filter"="blue")) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits=c(0,0.001))+
    xlab("Infection-paralysis-ratio")+
    ylab("Density")+
    labs(tag="C")+
    theme(legend.position="none")

p4 <- ggplot(res) + 
    geom_density(data=priors,aes(x=prop_immune,fill="Prior"),alpha=0.1) +
    geom_density(aes(x=prop_immune_groups.1,fill="Post-filter"),alpha=0.25) +
    theme_classic() +
    scale_fill_manual(name="",values=c("Prior"="black","Post-filter"="blue")) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits=c(0,1)) +
    xlab("Proportion initially immune")+
    ylab("Density")+
    labs(tag="D")+
    theme(legend.position="none")

fig2 <- p1 + p2 + p3 + p4 + plot_layout(ncol=2)

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


## Proportion of trajectories that have ended by now
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
