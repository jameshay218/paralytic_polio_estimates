library(Hmisc)
library(ggplot2)
library(patchwork)
library(dplyr)
library(paletteer)
setwd("~/Documents/GitHub/paralytic_polio_estimates/")
source("simulation_functions.R")

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
        all_traj[[i]] <- traj_future
    }
}
traj <- as_tibble(do.call("bind_rows",all_traj))

setwd("~/Documents/GitHub/paralytic_polio_estimates/sims_para/")
all_para <-  NULL
for(i in 1:1000){
    if(file.exists(paste0("para_",i,".RData"))){
        load(paste0("para_",i,".RData"))
        all_para[[i]] <- traj_future_par
    }
}
para <- as_tibble(do.call("bind_rows",all_para))

setwd("~/Documents/GitHub/paralytic_polio_estimates/sims_traj_vacc//")
all_traj_vacc <-  NULL
for(i in 1:1000){
    if(file.exists(paste0("traj_vacc_",i,".RData"))){
        load(paste0("traj_vacc_",i,".RData"))
        all_traj_vacc[[i]] <- traj_future_vacc
    }
}
traj_vacc <- as_tibble(do.call("bind_rows",all_traj_vacc))

setwd("~/Documents/GitHub/paralytic_polio_estimates/sims_para_vacc/")
all_para_vacc <-  NULL
for(i in 1:1000){
    if(file.exists(paste0("para_vacc_",i,".RData"))){
        load(paste0("para_vacc_",i,".RData"))
        all_para_vacc[[i]] <- traj_future_par_vacc
    }
}
para_vacc <- as_tibble(do.call("bind_rows",all_para_vacc))

traj <- traj %>% left_join(res %>% select(sim, date_start)) %>%
    mutate(t = as.Date(date_start + t)) %>%
    drop_na() %>%
    mutate(vacc_prop = 0)
traj_vacc <- traj_vacc %>% left_join(res %>% select(sim, date_start)) %>%
    mutate(t = as.Date(date_start + t)) %>%
    drop_na()
para <- para %>% left_join(res %>% select(sim, date_start))%>%
    mutate(t = as.Date(date_start + t))%>%
    drop_na()%>%
    mutate(vacc_prop = 0)
para_vacc <- para_vacc %>% left_join(res %>% select(sim, date_start))%>%
    mutate(t = as.Date(date_start + t))%>%
    drop_na()

traj <- bind_rows(traj, traj_vacc)
para <- bind_rows(para, para_vacc)

traj$vacc_prop <- as.factor(traj$vacc_prop)
para$vacc_prop <- as.factor(para$vacc_prop)

res <- res %>% filter(date_start > "2022-05-01")

samps <- sample(unique(res$sim),1000,replace=TRUE)
samps1 <- sample(unique(res$sim),100,replace=TRUE)
ggplot(traj %>% filter(vacc_prop == 0) %>%
           filter(t <= "2022-08-01") %>%
           filter(sim %in% samps1) %>%
           group_by(sim) %>%
           filter(t != min(t)) %>%
           left_join(res)) + 
    geom_line(aes(x=t,y=inc,group=sim,col=Re>1),size=0.25,alpha=0.5) +
    geom_vline(xintercept=as.Date("2022-07-18")) +
    scale_y_log10()


ggplot(para %>% filter(vacc_prop == 0) %>%
           filter(sim %in% samps)) + geom_line(aes(x=t,y=para,group=interaction(sim,rep)))

traj_prop <- traj %>% mutate(y=inc>0) %>% group_by(t,vacc_prop) %>% dplyr::summarize(y=sum(y),N=n()) %>% 
    group_by(t,vacc_prop) %>%
    mutate(prop=y/N,lower=binconf(y,N,0.05)[2], upper=binconf(y,N,0.05)[3]) 

traj_para_prop <- para %>% mutate(y=para>0) %>% group_by(t,vacc_prop) %>% dplyr::summarize(y=sum(y),N=n()) %>% 
    group_by(t,vacc_prop) %>%
    mutate(prop=y/N,lower=binconf(y,N,0.05)[2], upper=binconf(y,N,0.05)[3])

p1 <- traj_prop %>% ggplot() + 
    geom_ribbon(aes(x=t,ymin=lower,ymax=upper,fill=vacc_prop),alpha=0.25) +
    geom_line(aes(x=t,y=prop,col=vacc_prop),size=1) +
    scale_y_continuous(limits=c(0,0.75),breaks=seq(0,1,by=0.1)) +
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
n_prior <- 100000
immune_prop_pars <- get_beta_pars(0.8, 0.05)
prop_immune <- rbeta(n_prior, immune_prop_pars[1],immune_prop_pars[2])

prob_paralysis_mean <- 1/2000
prob_paralysis_var <- 1/2000/10
prob_para <- rnorm(n_prior, prob_paralysis_mean, prob_paralysis_var)
prob_para[prob_para < 0] <- 0
prob_para[prob_para >1] <- 1

priors <- data.frame(R0 = rlnorm(n_prior,log(1),1),
                     prop_immune=prop_immune,
                     prob_paralysis=prob_para) %>%
    mutate(Re = R0 * (1-prop_immune))

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
    geom_density(aes(x=prob_paralysis,fill="Post-filter"),alpha=0.25) +
    theme_classic() +
    scale_fill_manual(name="",values=c("Prior"="black","Post-filter"="blue")) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits=c(0,0.03))+
    xlab("Infection-paralysis-ratio")+
    ylab("Density")+
    labs(tag="C")+
    theme(legend.position="none")

p4 <- ggplot(res) + 
    geom_density(data=priors,aes(x=prop_immune,fill="Prior"),alpha=0.1) +
    geom_density(aes(x=prop_immune,fill="Post-filter"),alpha=0.25) +
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
table1 <- para %>% 
    group_by(sim, vacc_prop) %>% 
    dplyr::summarize(y=sum(para)) %>%
    group_by(vacc_prop) %>%
    dplyr::summarize(mean_para=mean(y),
              median_para=median(y),
              lower95=quantile(y,0.025), lower80=quantile(y,0.1),
              upper80=quantile(y,0.9),upper95=quantile(y,0.975))
## Proportion of outbreaks with more than the initial case of paralysis
table2 <- para %>% 
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

## Proportion of outbreaks with more than 100 further cases
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
