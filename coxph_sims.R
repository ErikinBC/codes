# This script simulates the relative efficiency of the Cox-PH model
rm(list=ls())
ll.cran <- c('magrittr','glmnet','survival','flexsurv')
for (k in ll.cran) { library(k,character.only = T,warn.conflicts = F,quietly = T)}

# Directory
dir.base <- "C:/Users/erikinwest/Documents/Research/BL/surv_ph"
setwd(dir.base)

# Call in scripts
source('coxph_funs.R')

##########################################################
## ------ (1) HOW DOES COX DO ON ITS OWN TERMS? ------- ##

# Only Exponential, Weibull, and Gumpertz share the PH assumption (lognormal or loglogistic do not)

# Set up parameters
nsim <- 250
n=100
p=1
beta <- 0.5
dist.list <- c('exponential','weibull','gompertz')
lam=1; alpha=0.5
param.list <- list(exponential=list(lam=lam,alpha=NULL),
                   weibull=list(lam=lam,alpha=alpha),gompertz=list(lam=lam,alpha=alpha))
cens.list <- seq(0.0,0.8,0.2)

sims.store <- list()
i=0
for (didx in 1:length(dist.list)) {
  dist <- dist.list[didx]
  lam <- param.list[[didx]]$lam
  alpha <- param.list[[didx]]$alpha
  for (cens in cens.list) {
    i=i+1
    # Run simulation
    temp.sim <- sim.coxph.para(nsim=nsim,n=n,p=p,beta=beta,dist=dist,lam=lam,alpha=alpha,cens=cens)
    # Add on the information
    temp.sim <- data.frame(temp.sim,dist=dist,cens0=cens,lam,
                           alpha=ifelse(is.null(alpha),'NULL',alpha),n,beta0=beta) 
    # Store
    sims.store[[i]] <- temp.sim
    print(i)
  }
}

# Combine
sims.tbl <- do.call('rbind',lapply(sims.store,function(df) 
  mutate(df,alpha=ifelse(alpha=='NULL',NA,alpha)))) %>% tbl_df
# Add on the coverage
# with(sims.tbl,(beta0 > lb) & (beta0 < ub)) %>% mean

# Summary stats
sims.sum <- sims.tbl %>% 
  group_by(mdl,dist,cens0) %>% 
  mutate(bias=beta-beta0,coverage=((beta0 > lb) & (beta0 < ub)) ) %>% 
  summarise(coverage=mean(coverage),bmu=mean(beta),bvar=var(beta)) %>% tbl_df 


 # Get long format
sims.long <-
sims.sum %>% filter(as.character(mdl)==as.character(dist) | mdl=='coxph') %>% 
  mutate(mdl=ifelse(mdl=='coxph','coxph','parametric'),bias=abs(bmu-beta),
         mse=bias^2+bvar) %>% dplyr::select(-c(bmu,bvar)) %>% 
  gather(measure,value,-(mdl:cens0))

ggplot(sims.long,aes(x=cens0,y=value,color=mdl)) + 
  geom_point(size=4) + 
  facet_grid(measure~dist,scales='free_y') + 
  background_grid(major='xy',minor='none')

sims.dens <- sims.tbl %>% filter(as.character(mdl)==as.character(dist) | mdl=='coxph') %>%
  mutate(mdl=ifelse(mdl=='coxph','coxph','parametric'))

ggplot(sims.dens,aes(x=beta,y=..density..,fill=mdl)) + 
  geom_density(alpha=0.5,color='black') + 
  facet_grid(cens0 ~ dist,scales='free')


###########################################
## ------ (2) NON-PH SIMULATIONS ------- ##


###############################################
## ------ (3) CURE MODEL SIMULATIONS ------- ##




