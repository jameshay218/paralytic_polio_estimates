library(Hmisc)
library(ggplot2)
library(patchwork)
library(dplyr)
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

traj <- traj %>% left_join(res %>% select(sim, date_start)) %>%
    mutate(t = as.Date(date_start + t)) %>%
    drop_na()
para <- para %>% left_join(res %>% select(sim, date_start))%>%
    mutate(t = as.Date(date_start + t))%>%
    drop_na()

res <- res %>% filter(date_start > "2022-05-01")

samps <- sample(unique(res$sim),1000,replace=TRUE)
ggplot(traj %>% filter(sim %in% samps)) + geom_line(aes(x=t,y=inc,group=sim))
ggplot(para %>% filter(sim %in% samps)) + geom_line(aes(x=t,y=para,group=interaction(sim,rep)))

traj_prop <- traj %>% mutate(y=inc>0) %>% group_by(t) %>% dplyr::summarize(y=sum(y),N=n()) %>% 
    group_by(t) %>%
    mutate(prop=y/N,lower=binconf(y,N,0.05)[2], upper=binconf(y,N,0.05)[3]) 

traj_para_prop <- para %>% mutate(y=para>0) %>% group_by(t) %>% dplyr::summarize(y=sum(y),N=n()) %>% 
    group_by(t) %>%
    mutate(prop=y/N,lower=binconf(y,N,0.05)[2], upper=binconf(y,N,0.05)[3])

p1 <- traj_prop %>% ggplot() + 
    geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25,fill="purple") +
    geom_line(aes(x=t,y=prop),color="purple",size=1) +
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25)) +
    scale_x_date(limits=as.Date(c("2022-08-22","2023-07-01")),
                 breaks="1 month") +
    theme_classic() +
    ylab("Proportion of simulations with\n ongoing infections") +
    theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          panel.grid.major=element_line(size=0.25,color="grey70")) +
    labs(tag="A")
p2 <- traj_para_prop %>% 
    ggplot() + 
    geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25,fill="orange") +
    geom_line(aes(x=t,y=prop),color="orange",size=1)+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+
    scale_x_date(limits=as.Date(c("2022-08-22","2023-07-01")),
                 breaks="1 month")+
    theme_classic()+
    ylab("Proportion of simulations with\n ongoing cases of paralysis") +
    xlab("Date") + 
    theme(panel.grid.major=element_line(size=0.25,color="grey70"),
          axis.text.x=element_text(angle=45,hjust=1))+
    labs(tag="B")
fig1 <- p1/p2

## Generate prior draws for density plots
n_prior <- 100000
immune_prop_pars <- get_beta_pars(0.8, 0.05)
prop_immune <- rbeta(n_prior, immune_prop_pars[1],immune_prop_pars[2])

prob_paralysis_mean <- 1/2000
prob_paralysis_var <- 1/2000/10
prob_para <- rnorm(n_prior, prob_paralysis_mean, sqrt(prob_paralysis_var))
prob_para[prob_para < 0] <- 0
prob_para[prob_para >1] <- 1

priors <- data.frame(R0 = rlnorm(n_prior,log(1),1),
                     prop_immune=prop_immune,
                     prob_paralysis=prob_para) %>%
    mutate(Re = R0 * prop_immune)

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
para %>% 
    group_by(sim) %>% 
    dplyr::summarize(y=sum(para)) %>%
    ungroup() %>% pull(y) -> x
mean(x)
quantile(x,c(0.025,0.1,0.5,0.9,0.975))

## Proportion of outbreaks with more than the initial case of paralysis
length(x[x > 1])/length(x)

## Proportion of trajectories that have ended by now
traj %>% filter(sim %in% samps) %>% 
    group_by(sim) %>% 
    filter(t >= date_start) %>%
    dplyr::summarize(y=sum(inc)) %>%
    ungroup() %>% pull(y) -> y
mean(y)
quantile(y,c(0.025,0.1,0.5,0.9,0.975))

## Proportion of outbreaks with more than 100 further cases
length(y[y > 100])/length(y)
setwd("~/Documents/GitHub/paralytic_polio_estimates/")
ggplot2::ggsave("figures/fig1.pdf",fig1,width=7,height=6)
ggplot2::ggsave("figures/fig2.pdf",fig2,width=7,height=6)
