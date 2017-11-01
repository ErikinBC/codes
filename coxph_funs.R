# This script has the survival dgp properties

# # Notes on censoring, assuming C_i ~ Exp(gam) then P(T_i > C_i) = censoring %:
# lam=1;gam=2;alpha=0.5
# # (i) Exponential
# # gam/(lam+gam) if T_i ~ Exp(lam)
# integrand <- function(x,lam,gam,alpha) { exp(-(lam*x + gam*x)) }
# gam * integrate(integrand,lower=0,upper=Inf,lam=lam,gam=gam,alpha=alpha)$value
# mean( rexp(n=10000,lam) > rexp(n=10000,gam) )
# # (ii) Weibull
# # gam * int[0,infty] exp(-(lam*c_i^alpha + gam*ci)) dci if T_i ~ Weibull(lam,alpha)
# # Note R's scale param is equiv to 1/lam
# integrand <- function(x,lam,gam,alpha) { exp(-(lam*x^alpha + gam*x)) }
# gam * integrate(integrand,lower=0,upper=Inf,lam=lam,gam=gam,alpha=alpha)$value
# mean( rweibull(n=10000,scale=1/lam,shape=alpha) > rexp(n=10000,gam) )
# # (iii) Gompertz
# # gam * int[0,infty] exp( (lam/alpha)*(exp(alpha*c) -1) - gam*c ) dc
# integrand <- function(x,lam,gam,alpha) { exp( -(lam/alpha)*(exp(alpha*x)-1) -gam*x) }
# gam * integrate(integrand,lower=0,upper=Inf,lam=lam,gam=gam,alpha=alpha)$value
# mean( rgompertz(n=10000,rate=lam,shape=alpha) > rexp(n=10000,gam) )

# -------- DATA GENERATING PROCESS FOR PH MODEL WITH COVARIATES ---------- #

# ss=1;n=100;p=1;beta=0.5;dist='weibull';lam=1;alpha=2;cens=0.2
ph.dgp <- function(ss,n,p,beta,dist,lam,alpha=NULL,cens) {
  if ( !( dist %in% c('exponential','weibull','gompertz') ) ) {
    print('Please pick a valid distribution!'); break
  }
  # begin dgp
  set.seed(ss)
  # Design matrix
  X <- matrix(rnorm(n * p),ncol=p)
  eta <- as.vector( X %*%  beta)
  exp.eta <- exp(eta)
  # Sample one point
  exp.eta1 <- median(exp.eta) #sample(exp.eta,1)
  # Survival times
  logu <- log(runif(n))
  if (dist == 'exponential') {
    TT <- -logu/(lam * exp.eta )
  } else if (dist == 'weibull' ) {
    TT <- ( -logu/(lam * exp.eta ) )^(1/alpha)
  } else if (dist == 'gompertz') {
    TT <- (1/alpha)*log(1 - (alpha*logu)/(lam*exp.eta ) )
  } else { print('huh?!')}
  # Use optim() to find the value of the rate parameter of the exponential distribution
  #   for random censoring
  if (dist == 'exponential') {
    integrand <- function(x,lam,gam,alpha,eeta) { exp(-(eeta*lam*x + gam*x)) }
  } else if (dist == 'weibull' ) {
    integrand <- function(x,lam,gam,alpha,eeta) { exp(-(eeta*lam*x^alpha + gam*x)) }
  } else if (dist == 'gompertz') {
    integrand <- function(x,lam,gam,alpha,eeta) { exp( -eeta*(lam/alpha)*(exp(alpha*x)-1) -gam*x) }
  } else { print('huh?!')}
  # Function to minimize
  fr <- function(gam,lam,alpha,eeta,ig,cens) {
    ((gam * integrate(ig,lower=0,upper=10,lam=lam,gam=gam,alpha=alpha,eeta=eeta)$value) - cens)^2
  }
  # Define the line interval to search over
  uint <- c(0,10)
  # Get the gamma for each
  gam.star <- sapply(exp.eta,function(ee) 
    optimize(f=fr,interval=uint,lam=lam,alpha=alpha,eeta=ee,ig=integrand,cens=cens)$minimum )
  # Check
  if ( any(abs(gam.star - max(uint))<1e-5) ) { print('Increase the search interval'); break }
  # Now generate the random censoring times
  CC <- rexp(n,rate=gam.star)
  # Now calculate the observed time point
  is.cens <- (CC<TT)
  Tobs <- TT
  Tobs[is.cens] <- CC[is.cens]
  # Return the data frame
  dat <- data.frame(time=Tobs,event=ifelse(is.cens,0,1),X)
  return(dat)
}


# ------ COMPARE EFFICIENCY OF COXPH MODEL TO PARAMETRIC SURVIVAL MODEL -------- #

# nsim=10;n=100;p=1;beta=0.5;dist='weibull';lam=1;alpha=2;cens=0.2
sim.coxph.para <- function(nsim,n,p,beta,dist,lam,alpha=NULL,cens) {
  # Storae
  results.store <- list()
  for (k in 1:nsim) {
    # Generate the data
    dat <- ph.dgp(ss=k,n=n,p=p,beta=beta,dist=dist,lam=lam,alpha=alpha,cens=cens)
    # Get survival object
    So <- with(dat,Surv(time=time,event=event))
    # Get design matrix
    X <- as.matrix(dat[,-(1:2)])
    # Fit each model and get the coefficient and standard errors
    # - (i) Cox PH - #
    mdl.cox <- coxph(So ~ X)
    beta.cox <- summary(mdl.cox)$coef[1]
    se.cox <- summary(mdl.cox)$coef[3]
    # - (ii) Exponential - #
    mdl.exp <- flexsurvreg(So ~ X,dist='exponential')
    beta.exp <- coef(mdl.exp)['X']
    se.exp <- sqrt(diag(mdl.exp$cov)['X'])
    # - (iii) Weibull - #
    if (class(lapply(1,function(ll) try(flexsurvreg(So ~ X,dist='weibullPH')))[[1]])=='try-error') {
      beta.weibull <- NA
      se.weibull <- NA
    } else {
      mdl.weibull <- flexsurvreg(So ~ X,dist='weibullPH')
      beta.weibull <- coef(mdl.weibull)['X']
      se.weibull <- sqrt(diag(mdl.weibull$cov)['X'])
    }
    # - (iv) Gompertz - #
    if (class(lapply(1,function(ll) try(flexsurvreg(So ~ X,dist='gompertz')))[[1]])=='try-error') {
      beta.gompertz <- NA
      se.gomperts <- NA
    } else {
      mdl.gompertz <- flexsurvreg(So ~ X,dist='gompertz')
      beta.gompertz <- mdl.gompertz$coef['X']
      se.gomperts <- sqrt(diag(mdl.gompertz$cov)['X'])
    }
    # Store the slices and calculate CI
    temp.results <- data.frame(mdl=c('coxph','exponential','weibull','gompertz'),
                               rbind(c(beta.cox,se.cox),c(beta.exp,se.exp),
                                     c(beta.weibull,se.weibull),c(beta.gompertz,se.gomperts)))
    colnames(temp.results)[2:3] <- c('beta','se')
    temp.results <- mutate(temp.results,lb=beta-se*qnorm(0.975),ub=beta+se*qnorm(0.975))
    # Add on censoring information
    temp.results <- cbind(temp.results,cens=1-mean(dat$event))
    # Store and update
    results.store[[k]] <- temp.results
    if (mod(k,25)==0) { print(k) }
  }
  # Combine and save
  results.tbl <- do.call('rbind',results.store)
  # Return
  return(results.tbl)
}

