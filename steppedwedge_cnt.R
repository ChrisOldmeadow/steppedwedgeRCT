# Name: sim_power_steppedwedge_cnt.R
# Author: chris Oldmeadow
# Date: 7/1/2021
# Summary:  simulate power for a stepped wedge trial with a continuous outcome



library(plyr) # data functions
library(lme4) # simulating the data and fitting the mixed models
library(pbapply)  # progress bar for simulations

# using a mixed effects, quick method



designSWT <-function(n.site ,
                  n.id ,
                  min.base, # basleine time for first site
                  n.period,
                  step) {    #time between steps, eg every 2 months 
   
    expdat <- expand.grid( id = seq(n.id),time = seq(n.period),site =seq(n.site)) 
    sites <- data.frame(site = seq(n.site), step = seq(from = (min.base+1), 
                                                       to = (n.period-(min.base+1)) ,
                                                       by = step))
    df<-merge(expdat,sites,by="site")
    df$period <- ifelse(df$time <  df$step,0,
                        ifelse(df$time < (df$step + step), NA, 1 ))  # period is missing during the intervention period and 1 thereafter

    df <- na.exclude(df)
    return(df)
}


simSWRCT <- function(df,
                     nreps ,
                     m1 , # mean of pre intervetion
                     m2 , # mean of post intervetion
                     s, # standard deviation of the outcome
                     rho) { # intra-class correlation
    beta <- c(m1, (m2-m1))
    names(beta)<-c("(Intercept)","period")
    sigma2 <- s*s 
    theta <-sqrt((rho*sigma2)/(1-rho))/s
    names(theta)<-c("site.(Intercept)")
    ss <- simulate(~ period +  (1 | site), nsim = nreps, 
                   family = gaussian,
                   newdata = df, 
                   newparams = list(theta = theta,sigma = s,   beta = beta))

    return(ss)
}


# fits a linear mixed model to the simulated data

fitsim <- function(i,dat = simdat,out = ss) {
    # dat is the data structure
    # out is the list of simulated outcomes

    return <- tryCatch({
        dat$y <- out[[i]]
        res<-coef(summary(lmerTest::lmer(y ~ period + time + (1|site), data = dat)))["period", ]
        names(res)<-c("est","se","df","t","p")
        return(res)
    },
    warning =function(e) {
        #message(e)  # print error message
        return(c(est=NA, se=NA, df = NA, t = NA, p=NA))
    },

    error=function(e) {
        #message(e)  # print error message
        return(c(est=NA, se=NA, df = NA, t=NA, p=NA))
    })

    return(return)
}


# now the power calculations
simdat <- designSWT(n.site = 5,
              n.id = 10,
              min.base = 3,
              n.period = 16,
              step = 2)

table(simdat$site,simdat$time)
table(simdat$site,simdat$period)
 


ss <- simSWRCT(simdat, m1 = 1 ,m2 = 0.55 ,s = 1, rho = .2,nreps = 5000)
fitAll <- t(pbsapply(seq(5000), function(i) fitsim(i)))

pow <- mean( fitAll[,"p"] < 0.05, na.rm = TRUE)
pow    

save.image(file=here('power.RData'))

# load(file = here('power.RData'))



 
