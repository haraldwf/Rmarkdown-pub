# R functions (http://www.r-project.org/) for estimating DCIS progression rates
#   Harald Weedon-Fekjær <harald.weedon-fekjar@medisin.uio.no>, 2019

scrdata <- data.frame(t      = c(0.67,2,4,6.5),            
                      est    = c(0.54, 0.66, 0.74, 1.10),  
                      est.se = c(1.23, 1.15, 1.18, 1.15))
  # t      = “Time since screening”
  # est    = “Estimate relative risk compared to initial screening”
  # est.se = ”Estimated standard error of est (se previous variable)”

  
  
# -----------------------------------------------------------------
# ----- Estimate mean sojourn time and screening sensitivity: -----
# -----------------------------------------------------------------
estM <- function(scrdata) {
  # Estimate DCIS mean sojourn time and screening sensitivity from screening data
  # 
  # Define log likelihood:
  find.logLik <- function(M,S) {
    exp.valesI <- rep(NA,dim(scrdata)[1]) 
    for (i in 1:length(exp.valesI)) {
      exp.valesI[i] <- ((1/S)*(1-S))*(1-pexp(scrdata$t[i],rate=1/M))*S + 
                       pexp(scrdata$t[i],rate=1/M)
    }
    loglik   <- 0
    for (i in 1:length(exp.valesI)) {
      loglik <- loglik + log(x=dnorm(log(scrdata$est[i]),
                             mean=log(exp.valesI[i]),sd=log(scrdata$est.se[i])))
    }
    return(loglik)
  } 
  # Function for use in optimization of likelihood:
  tmp <- function(par) {
    ret <- - find.logLik(M=par[1],S=par[2])
    return(ret)
  }
  # Optimization of likelihood:
  return(optim(c(2,.9),fn=tmp, lower=c(0,0), upper=c(Inf,1),method="L-BFGS-B")$par)
}  



# ---------------------------------
# ----- Supportive function: -----
# ---------------------------------
ExpScrByTime <- function(pint,M,S,timeLscr=NA,noscr=1) {
  # Calculates expected number at screening
  # pint	= Insidence of pre-clinical screening detectable DCIS at initial screening
  # M        = Mean sojourn time
  # S        = Screening test sensitivity
  # timeLscr = Time since last screening
  # mpscr    = No screened
  if (is.na(timeLscr)) {
    scr.rate <- pint*S     # (No earlier screening)
  } else {
    scr.rate <- (pint*(1-S))*(1-pexp(timeLscr,rate=1/M))*S + 
                 pint*pexp(timeLscr,rate=1/M)*S
  }
  return(scr.rate*noscr)
}

