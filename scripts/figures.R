library(Hmisc)
library(ggplot2)
library(patchwork)
library(dplyr)

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
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.2)) +
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
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.2))+
    scale_x_date(limits=as.Date(c("2022-08-22","2023-07-01")),
                 breaks="1 month")+
    theme_classic()+
    ylab("Proportion of simulations with\n ongoing cases of paralysis") +
    xlab("Date") + 
    theme(panel.grid.major=element_line(size=0.25,color="grey70"),
          axis.text.x=element_text(angle=45,hjust=1))+
    labs(tag="B")
p1/p2


ggplot(res) + 
    geom_density(aes(x=Re),fill="blue",alpha=0.25) +
    theme_classic() +
    geom_vline(xintercept=1) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(breaks=seq(0,8,by=1))

para %>% filter(sim %in% samps) %>% group_by(sim) %>% dplyr::summarize(y=sum(para)) %>%
    ungroup() %>% pull(y) -> x
mean(x)
quantile(x,c(0.025,0.1,0.5,0.9,0.975))

length(x[x > 1])/length(x)
