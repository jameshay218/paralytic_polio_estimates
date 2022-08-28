library(ggsci)
library(Hmisc)
library(ggplot2)
library(patchwork)
library(dplyr)
library(paletteer)
setwd("~/Documents/GitHub/paralytic_polio_estimates/")
source("simulation_functions_twoimmune.R")
max_date <- "2023-04-01"
summarize <- dplyr::summarize
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

## Get cumulative incidenc
traj <- traj %>% group_by(sim) %>% mutate(cumu_inc = cumsum(inc),
                                          cumu_para = cumsum(para))
traj <- traj %>% ungroup() %>% mutate(month = round_date(t, "month"))

## Rockland county trajectories only valid if they have some cases in May, June, July and August
traj %>% filter(month %in% as.Date(c("2022-05-01","2022-06-01","2022-07-01","2022-08-01"))) %>% group_by(sim, month) %>% dplyr::summarize(monthly_inc=sum(inc)) %>% mutate(keep=monthly_inc > 0) %>% select(sim, keep, month) %>% pivot_wider(id_cols = sim,values_from="keep",names_from="month") %>%
    filter(`2022-05-01`==TRUE,`2022-06-01`==TRUE,`2022-07-01`==TRUE,`2022-08-01`==TRUE) %>% pull(sim) -> wastewater_sims
traj <- traj %>% filter(sim %in% wastewater_sims)
res <- res %>% filter(sim %in% wastewater_sims)
## Assume NYC is seeded 4 weeks later
traj_nyc <- traj_nyc %>% 
    left_join(res %>% select(sim, date_start)) %>%
    mutate(t = as.Date(date_start + t) + 28) %>%
    drop_na() 

traj_nyc <- traj_nyc %>% group_by(sim) %>% mutate(cumu_inc = cumsum(inc),
                                          cumu_para = cumsum(para))
traj_nyc <- traj_nyc %>% ungroup() %>% mutate(month = round_date(t, "month"))
## NYC trajectories only acceptable if they have cases in June and July
traj_nyc %>% filter(month %in% as.Date(c("2022-06-01","2022-07-01"))) %>% group_by(sim, month) %>% dplyr::summarize(monthly_inc=sum(inc)) %>% mutate(keep=monthly_inc > 0) %>% select(sim, keep, month) %>% pivot_wider(id_cols = sim,values_from="keep",names_from="month") %>%
    filter(`2022-06-01`==TRUE,`2022-07-01`==TRUE) %>% pull(sim) -> wastewater_sims_nyc
traj_nyc <- traj_nyc %>% filter(sim %in% wastewater_sims_nyc)


traj$vacc_prop <- as.factor(traj$vacc_prop)
trajectories <- traj %>% filter(vacc_prop==1) %>% dplyr::select(-c(vacc_prop,rep))
save(trajectories,file="../outputs/fitted_trajectories.RData")
save(traj_nyc,file="../outputs/nyc_trajectories.RData")
save(res,file="../outputs/filtered_parameters.RData")

## Pad prematurely ended trajectories with zeros
N_tot <- length(unique(trajectories$sim))
## Pad up to this date
end_date <- max(trajectories$t)
## Create dummy rows from this date onward, unique to each sim
max_t1 <- trajectories %>% group_by(sim) %>% filter(t==max(t)) %>% select(sim, t, date_start) %>% mutate(t = t + 1) %>% rename(max_t=t) %>% mutate(max_date=end_date)

## Create the dummy rows
extra_rows <- max_t1 %>%
    filter(max_t < max_date) %>%
    group_by(sim) %>%
    distinct() %>%
    do(data.frame(sim=.$sim, t=seq(.$max_t,.$max_date,by="1 day")))
extra_rows <- extra_rows %>% mutate(inc=0,inc_s=0,inc_ps=0,para=0,para_s=0,para_ps=0, Rt=0)

## Add the dummy rows in
extra_rows <- extra_rows %>% left_join(trajectories %>% select(sim, date_start) %>% distinct())
trajectories <- bind_rows(trajectories,extra_rows) %>% arrange(sim, t)
trajectories <- trajectories %>% mutate(Rt = ifelse(is.na(Rt),0, Rt))
trajectories <- trajectories %>% group_by(sim) %>% mutate(cumu_para=cumsum(para)) %>% ungroup()
trajectories <- trajectories %>% group_by(sim) %>% 
    ## If cases aren't what they were 7 days ago, then the epidemic is still ongoing
    mutate(cumu_inc=cumsum(inc)) %>% 
    mutate(ongoing_7= cumu_inc != lag(cumu_inc,7,nafill(0))) %>%
    ungroup()

## No further cases observed
tmp_flags <- trajectories %>% filter(t <= "2022-10-01") %>% group_by(sim) %>% dplyr::summarise(y=sum(para)) %>% mutate(no_further_cases=y==1)

traj_summary1 <- trajectories %>% filter(t <= max_date) %>% 
    #mutate(y=Rt<1) %>% 
    mutate(y=ongoing_7 == TRUE) %>%
    group_by(t) %>%
    dplyr::summarize(
        median_para=median(cumu_para),lower1=quantile(cumu_para,0.1),upper1=quantile(cumu_para,0.9),
        y=sum(y),N=n(),prop=sum(y))%>% group_by(t) %>% 
    mutate(
        prop=y/N,lower=binconf(y,N)[2],upper=binconf(y,N)[3]) %>%
    mutate(model="Current data")
traj_summary2 <- trajectories %>% left_join(tmp_flags) %>%
    filter(no_further_cases==TRUE) %>%
    filter(t <= max_date) %>% 
    #mutate(y=Rt<1) %>% 
    mutate(y=ongoing_7 == TRUE) %>%
    group_by(t) %>% 
    dplyr::summarize(
        median_para=median(cumu_para),lower1=quantile(cumu_para,0.1),upper1=quantile(cumu_para,0.9),
        y=sum(y),N=n())%>% group_by(t) %>% 
    mutate(prop=y/N,lower=binconf(y,N)[2],upper=binconf(y,N)[3])%>%
    mutate(model="No cases by October 1st")
traj_summary3 <- trajectories %>% left_join(tmp_flags) %>%
    filter(no_further_cases==FALSE) %>%
    filter(t <= max_date) %>% 
    #mutate(y=Rt<1) %>% 
    mutate(y=ongoing_7 == TRUE) %>%
    group_by(t) %>% 
    dplyr::summarize(
        median_para=median(cumu_para),lower1=quantile(cumu_para,0.1),upper1=quantile(cumu_para,0.9),
                     y=sum(y),N=n())%>% group_by(t) %>% 
    mutate(prop=y/N,lower=binconf(y,N)[2],upper=binconf(y,N)[3])%>%
    mutate(model="Further cases reported\n by October 1st")


traj_summary <- bind_rows(traj_summary1,traj_summary2,traj_summary3)
p1 <- ggplot(traj_summary) + 
    geom_ribbon(aes(x=t,ymin=lower,ymax=upper,fill=model),alpha=0.25) +
    geom_line(aes(x=t,y=prop,col=model)) +
    geom_text(data=data.frame(x=as.Date(c("2022-06-22","2022-08-20")),
                              y=1.12,label=c("First case of\n paralysis","Last observation")),aes(x=x,y=y,label=label),size=2.5) +
    scale_y_continuous(limits=c(0,1.15),breaks=seq(0,1,by=0.2)) +
    
    geom_segment(data=data.frame(x=as.Date(c("2022-06-22","2022-08-20")),
                                 xend=as.Date(c("2022-06-22","2022-08-20")),
                                 y=0,yend=1),
                 aes(x=x,xend=xend,y=y,yend=yend),
               linetype="dashed") +    ylab("Probability epidemic is ongoing\n (cases in previous 7 days)")+
    theme_classic() +
    scale_fill_nejm() +
    scale_color_nejm() +
    
    scale_x_date(limits=as.Date(c("2022-06-01",max_date)),
                 breaks="1 month") +    
    theme(legend.position="none",axis.text.x=element_blank(),
          axis.line.x=element_blank(),axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),axis.text.y=element_text(size=6),
          axis.title=element_text(size=8),
          legend.text=element_text(size=6),
          legend.title=element_text(size=6),
          panel.grid.major = element_line(size=0.1,color="grey70")) +
    labs(tag="A")

p2 <- ggplot(traj_summary) + 
    geom_ribbon(aes(x=t,ymin=lower1,ymax=upper1,fill=model),alpha=0.25) +
    geom_line(aes(x=t,y=median_para,col=model)) +
    geom_vline(xintercept=as.Date(c("2022-06-22","2022-08-20")),
linetype="dashed") +
    ylab("Cumulative paralysis cases") +
    scale_x_date(limits=as.Date(c("2022-06-01",max_date)),
                 breaks="1 month") +
    theme_classic() +
    xlab("Date") +
    scale_fill_nejm(name="Data assumption")+
    scale_color_nejm(name="Data assumption")+
    theme(legend.position="bottom",
          panel.grid.major = element_line(size=0.1,color="grey70"),
          axis.text.x=element_text(angle=45,hjust=1,size=6),
          axis.text.y=element_text(size=6),
          axis.title=element_text(size=8),
          legend.text=element_text(size=6),
          legend.title=element_text(size=6))+
    labs(tag="C")




p_final_size <- trajectories %>% group_by(sim) %>%
    mutate(y=cumsum(inc)) %>% 
    filter(t == "2022-08-20") %>% 
    mutate(model = "Current data") %>%
    ggplot() + geom_histogram(aes(x=y,fill=model),binwidth=10000,
                              color="black") +
    scale_x_continuous(breaks=seq(0,250000,by=100000))+
    scale_y_continuous(expand=c(0,0)) +
    theme_classic() +
    xlab("Estimated infections by 2022-08-20") + ylab("Simulation count") +
    theme(axis.text=element_text(size=6),
          axis.title=element_text(size=8),
          legend.position="none",
          legend.text=element_text(size=6),
          legend.title=element_text(size=6)) +
    scale_fill_nejm()+
    labs(tag="D")

p_re <- trajectories %>% 
    group_by(sim) %>%
    filter(t == min(t)) %>%
    mutate(model = "Current data") %>%
    ggplot() + geom_density(aes(x=Rt,fill=model),color="black",alpha=0.5) +
    geom_vline(xintercept=1,linetype="dashed") +
    scale_y_continuous(expand=c(0,0)) +
    theme_classic() +
    xlab("Initial effective reproductive number") + ylab("Density") +
    theme(axis.text=element_text(size=6),
          axis.title=element_text(size=8),
          legend.position="none",
          legend.text=element_text(size=6),
          legend.title=element_text(size=6)) +
    scale_fill_nejm() +
    labs(tag="B")



fig1 <- p1 + p2 + p_re +p_final_size +plot_layout(ncol=2,widths=c(2.5,1),byrow=FALSE)
fig1

## For every day after 2022-08-20, if we've seen another 0, 1, 2 or 3 
## cases of paralysis, flag

flag_traj <- trajectories %>% group_by(sim) %>% 
    filter(t > "2022-08-20", t <= max_date) %>%
    group_by(sim) %>% 
    mutate(cumu_para=cumsum(para)) %>%
    ungroup() %>%
    ## Flag when we've seen 1, 2 or 3 total cases after 20th August
    mutate(
            no_new_para = cumu_para == 0, ## Have still no cases
            one_new_para=cumu_para == 1, ## I have only seen one case
           two_new_para=cumu_para == 2, ## I have only seen two cases
           three_new_para=cumu_para == 3 ## I have only seen three cases
           ) %>%
    select(sim, t, no_new_para, one_new_para,two_new_para, three_new_para)


## Find total paralysis cases by end date
traj_cumu_para <- trajectories %>% 
    group_by(sim) %>% mutate(para=cumsum(para)) %>% 
    filter(t == max_date) %>%
    select(sim, para) %>%
    rename(total_para=para)
    
tmp <- flag_traj %>% left_join(traj_cumu_para)


omg1 <- tmp %>% group_by(sim) %>% filter(no_new_para == TRUE) %>%
    filter(t == max(t)) %>%
    mutate(last_date_no_cases = t) %>%
    select(sim, last_date_no_cases,total_para)

omg2 <- tmp %>% group_by(sim) %>% filter(one_new_para == TRUE) %>%
    filter(t == max(t)) %>%
    mutate(last_date_one_case = t) %>%
    select(sim, last_date_one_case,total_para)

omg3 <- tmp %>% group_by(sim) %>% filter(two_new_para == TRUE) %>%
    filter(t == max(t)) %>%
    mutate(last_date_two_cases = t) %>%
    select(sim, last_date_two_cases,total_para)

omg4 <- tmp %>% group_by(sim) %>% filter(three_new_para == TRUE) %>%
    filter(t == max(t)) %>%
    mutate(last_date_three_cases = t) %>%
    select(sim, last_date_three_cases,total_para)

tmp_all <- omg1 %>% left_join(omg2) %>% left_join(omg3) %>% left_join(omg4) %>% ungroup()

dates <- seq(as.Date("2022-08-20"),as.Date(max_date),by="1 day")

ns <- numeric(length(dates))
subset_flags1 <- tmp_all
subset_flags2 <- tmp_all
subset_flags3 <- tmp_all
subset_flags4 <- tmp_all

all_no_para <- NULL
all_one_para <- NULL
all_two_para <- NULL
all_three_para <- NULL
quantile_lower <- 0.1
quantile_upper <- 0.9
for(i in seq_along(dates)){
    print(i)
    subset_flags1 <- subset_flags1 %>% filter(
        !is.na(last_date_no_cases) &
        last_date_no_cases >= dates[i])
    subset_flags2 <- subset_flags2 %>% filter(
        !is.na(last_date_one_case) &last_date_one_case >= dates[i])
    subset_flags3 <- subset_flags3 %>% filter(
        !is.na(last_date_two_cases) &last_date_two_cases >= dates[i])
    subset_flags4 <- subset_flags4 %>% filter(
        !is.na(last_date_three_cases) &last_date_three_cases >= dates[i])
    ns[i] <- nrow(subset_flags1)
    all_no_para[[i]] <- subset_flags1 %>% 
        dplyr::summarize(median_para=median(total_para),
                         mean_para=mean(total_para),
                         lower=quantile(total_para,quantile_lower),
                  upper=quantile(total_para,quantile_upper)) %>%
        mutate(model="0", t =dates[i])
    
    all_one_para[[i]] <- subset_flags2 %>% 
        dplyr::summarize(median_para=median(total_para),
                         mean_para=mean(total_para),
                         lower=quantile(total_para,quantile_lower),
                  upper=quantile(total_para,quantile_upper)) %>%
        mutate(model="1", t =dates[i])
    
    all_two_para[[i]] <- subset_flags3 %>% 
        dplyr::summarize(median_para=median(total_para),
                         mean_para=mean(total_para),
                         lower=quantile(total_para,quantile_lower),
                  upper=quantile(total_para,quantile_upper)) %>%
        mutate(model="2", t =dates[i])
    
    all_three_para[[i]] <- subset_flags4 %>% 
        dplyr::summarize(median_para=median(total_para),
                         mean_para=mean(total_para),
                         lower=quantile(total_para,quantile_lower),
                  upper=quantile(total_para,quantile_upper)) %>%
        mutate(model="3", t =dates[i])
    
}
all_no_para <- do.call("bind_rows",all_no_para)
all_one_para <- do.call("bind_rows",all_one_para)
all_two_para <- do.call("bind_rows",all_two_para)
all_three_para <- do.call("bind_rows",all_three_para)

tmp_comb <- bind_rows(all_no_para,all_one_para, all_two_para, all_three_para)
fig2A <- tmp_comb %>%
    mutate(model = as.factor(model)) %>%
    mutate(model=paste0("Further paralysis cases observed: ", model)) %>%
    ggplot() + 
    geom_ribbon(aes(x=t,ymin=lower,ymax=upper,fill=model),alpha=0.05) +
    geom_line(aes(x=t,y=upper,col=model,linetype="90% quantiles"),size=0.25) +
    geom_line(aes(x=t,y=lower,col=model,linetype="90% quantiles"),size=0.25) +
    #geom_line(aes(x=t,y=upper,col=model),size=0.25) +
    #geom_line(aes(x=t,y=lower,col=model),size=0.25) +
    #geom_line(aes(x=t,y=lower,col=model),linetype="dotted") +
    #geom_line(aes(x=t,y=mean_para,col=model,linetype="Mean"),size=1) +
    
    geom_line(aes(x=t,y=median_para,col=model,linetype="Median"),size=0.8) +
    scale_y_continuous(expand=c(0,0),
                       limits=c(0,50)) +
    scale_fill_nejm(name="Additional paralysis\ncases observed") +
    scale_color_nejm(name="Additional paralysis\ncases observed") +
    scale_linetype_manual(name="",values=c("90% quantiles"="dashed",
                                    "Median"="solid")) +
    ylab(paste0("Projected incidence of paralysis by ",max_date)) +
    #scale_x_date(breaks="1 month") +
    #scale_x_date(breaks="1 month") +
    xlab("Date") +
    theme_classic() +
    theme(panel.grid.major = element_line(size=0.25,color="grey70"),
          panel.grid.minor = element_line(size=0.25,color="grey70"),
          legend.position="none",
          strip.background = element_blank(),
          axis.text=element_text(size=6),
          axis.title=element_text(size=8),
          strip.text=element_text(size=6)) +
    facet_wrap(~model,nrow=1)
fig2A

date_key <- c("2022-09-01"="Sep","2022-10-01"="Oct","2022-11-01"="Nov","2022-12-01"="Dec")
tmp_comb$date1 <- date_key[as.character(tmp_comb$t)]

tmp_comb1 <- tmp_comb %>% filter(t %in% as.Date(c("2022-09-01","2022-10-01","2022-11-01","2022-12-01"))) %>%
    mutate(date2=paste0("Cases observed by ", date1))
tmp_comb1$date2 <- factor(tmp_comb1$date1,levels=c("Cases observed by Sep","Cases observed by Oct","Cases observed by Nov","Cases observed by Dec"))
tmp_comb1$date1 <- factor(tmp_comb1$date1, levels=c("Sep","Oct","Nov","Dec"))
fig2B <- tmp_comb1 %>%   
    mutate(model=paste0("k=",model)) %>%
    ggplot() +
    geom_bar(aes(x=model,y=median_para,fill=model),stat="identity",col="black") +
    geom_errorbar(aes(x=model,ymin=lower,ymax=upper),width=0.5) +
    facet_wrap(~date1,nrow=1,strip.position="bottom") +
    xlab("Month by which k additional cases observed") +
    ylab(paste0("Projected paralysis cases\n by ", max_date)) +
    theme_classic() +
    scale_fill_nejm() +
    theme(legend.position="none",strip.background = element_blank(),
          panel.grid.major=element_line(size=0.1,color="grey70"),
          strip.placement = "outside"
         ) +
    scale_y_continuous(expand=c(0,0),limits=c(0,45))

ggsave(filename="~/Documents/GitHub/paralytic_polio_estimates/figures/fig1.pdf", fig1, height=6,width=8)
ggsave(filename="~/Documents/GitHub/paralytic_polio_estimates/figures/fig1.png", fig1, height=6,width=8,units='in',dpi=300)
ggsave(filename="~/Documents/GitHub/paralytic_polio_estimates/figures/fig2.pdf", fig2B, height=3,width=7)

fig2_alt <- (fig2B+labs(tag="A"))/(fig2B_nyc + labs(tag="B"))

ggsave(filename="~/Documents/GitHub/paralytic_polio_estimates/figures/fig2_alt.pdf", fig2_alt, height=5,width=7)
ggsave(filename="~/Documents/GitHub/paralytic_polio_estimates/figures/fig2_alt.png", fig2_alt, height=5,width=7,units="in",dpi=300)


## Numbers for paper
## % showing continued spread by August 20th
trajectories %>% filter(t == "2022-08-20") %>%
    group_by(ongoing_7) %>% tally() %>%
    pivot_wider(values_from=n,names_from=ongoing_7) %>%
    mutate(prop=`TRUE`/(`TRUE` + `FALSE`))

## Effective reproduction number at seeding
trajectories %>% group_by(sim) %>%
    filter(t == min(t)) %>%
    ungroup() %>%
    dplyr::summarize(mean_re=mean(Rt),median_re=median(Rt),
              lower90=quantile(Rt,0.1),upper90=quantile(Rt,0.9))

## Generation interval estimate
res$infectious_period_mean %>% mean
quantile(res$infectious_period_mean,c(0.1,0.9))

## Doubling time

## 22nd June number of infections
trajectories %>% filter(t == as.Date("2022-06-22")) %>% 
    ungroup() %>% summarize(mean_inf = mean(cumu_inc),
                            median_inf = median(cumu_inc),
                            lower90=quantile(cumu_inc,0.1),
                            upper90=quantile(cumu_inc,0.9))

## 18nd July number of infections
trajectories %>% filter(t == as.Date("2022-07-18")) %>% 
    ungroup() %>% summarize(mean_inf = mean(cumu_inc),
                            median_inf = median(cumu_inc),
                            lower90=quantile(cumu_inc,0.1),
                            upper90=quantile(cumu_inc,0.9))

## 20th August number of infections
trajectories %>% filter(t == as.Date("2022-08-20")) %>% 
    ungroup() %>% summarize(mean_inf = mean(cumu_inc),
                            median_inf = median(cumu_inc),
                            lower90=quantile(cumu_inc,0.1),
                            upper90=quantile(cumu_inc,0.9))

## 1st April 2023 number of infections
## 22nd June number of infections
trajectories %>% filter(t == as.Date(max_date)) %>% 
    ungroup() %>% summarize(mean_inf = mean(cumu_inc),
                            median_inf = median(cumu_inc),
                            lower90=quantile(cumu_inc,0.1),
                            upper90=quantile(cumu_inc,0.9))

## Number of paralytic polio by 1st April
trajectories %>% filter(t == as.Date(max_date)) %>% 
    ungroup() %>% summarize(mean_para = mean(cumu_para),
                            median_para = median(cumu_para),
                            lower90=quantile(cumu_para,0.1),
                            upper90=quantile(cumu_para,0.9))

## Number of paralytic polio by 20th August if no further cases by 1st Oct
traj_summary %>% filter(t == as.Date("2022-08-20")) %>%
    mutate(prop=1-prop,lower=1-lower,upper=1-upper)


## Number of paralytic polio by 1st April if 1 further case by 1st Oct
traj_summary %>% filter(t == max_date)


## Probability of not seeing any cases in NYC at all

## Number of paralysis in NYC if we see 1 case on or after 1st October by 1st April 2023
