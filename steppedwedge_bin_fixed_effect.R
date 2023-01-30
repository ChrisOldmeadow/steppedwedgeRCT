# Name: sim_power_steppedwedge_bin.R
# Author: chris Oldmeadow
# Date: 18/1/2021
# Summary:  simulate power for a stepped wedge trial with a dichotomous outcome



library(plyr) # data functions
library(lme4) # simulating the data and fitting the mixed models
library(pbapply)  # progress bar for simulations


# using a mixed effects, quick method
designSWT <-function(n.site ,
                     n.seq,  #The number of sequeces, a sequence is a collection of sites that step together
                     n.id ,  # The number of subjects in a site
                     n.period, # The number of data collection periods
                     t.step) {     #time between steps, eg every 2 months

    site.per.seq <- n.site/n.seq
    expdat <- expand.grid( id = seq(n.id),time = seq(n.period),site =seq(n.site) ) 
    sites <- data.frame(site = rep(seq(n.site),eaech = site.per.seq), step =
                        rep(seq(from = (t.step+1), 
                                to = (n.period-1) ,
                                by = t.step), each = site.per.seq))

    df<-merge(expdat,sites,by="site")
    df$period <- ifelse(df$time <  df$step, 0,
                        ifelse(df$time == df$step, NA, 1 ))  # period is missing during the intervention period and 1 thereafter

    df <- na.exclude(df)
    return(df)
}

simSWRCT <- function(df,
                     nreps ,
                     p1 , # prev of outcome at pre intervetion
                     p2 , # prev of outcome at post intervetion
                     rho) { # intra-class correlation
    o2<-p2/(1-p2)
    o1<- p1/(1-p1)

    beta <- c(log(o2), log(o1/o2))
    names(beta)<-c("(Intercept)","period")



    sigma2 <-(pi ^ 2) / 3
    theta <-sqrt((rho*sigma2)/(1-rho))
    names(theta)<-c("site.(Intercept)")

    ss <- simulate(~ period  + (1 | site), nsim = nreps, family = binomial, 
                   newdata = df, newparams = list(theta = theta,   beta = beta))


    return(ss)
}

# fits a logistic mixed model to the simulated data
 
fitsim <- function(i,dat = simdat,out = ss) {
    # dat is the data structure
    # out is the list of simulated outcomes

    return <- tryCatch({
        dat$resp <- out[[i]]
        res<-coef(summary(glm(resp ~ period + factor(site),
                                        family = binomial, 
                                        data=dat,
                                        verbose = FALSE)))["period", ]
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


simdat <- designSWT(n.site = 10,
                    n.seq = 5,
                    n.id = 40,
                    n.period = 36,
                    t.step = 5)

table(simdat$site,simdat$time)
table(simdat$site,simdat$period)
 


ss <- simSWRCT(simdat, p1 = .6 ,p2 = .7 , rho = .2,nreps = 1000)
fitAll <- t(pbsapply(seq(1000), function(i) fitsim(i)))

pow <- mean( fitAll[,"p"] < 0.008, na.rm = TRUE)
pow    




