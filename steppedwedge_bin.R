
library(plyr)
library(lme4)


# using a mixed effects, quick method



sim_designSWT <-function(n.site ,
                  n.id ,
                  n.period,
                  step)     #time between steps, eg every 2 months
{
   
  
  
 
expdat <- expand.grid( id = seq(n.id),time = seq(n.period),site =seq(n.site)) 


sites <- data.frame(site = seq(n.site), step = seq(from = (step), 
                                                   to = (n.period-1) ,
                                                   by = step))


df<-merge(expdat,sites,by="site")

df$period <- ifelse(df$time <=  df$step,0,
                    ifelse(df$time <= (df$step + step), NA, 1 ))  # period is missing during the intervention period and 1 thereafter



df <- na.exclude(df)

  
  
  return(df)
 
  
}

simdat <- sim_designSWT(n.site = 3,
              n.id = 12,
              n.period = 5,
              step = 1)

  table(simdat$site,simdat$time)
table(simdat$site,simdat$period)
 


pow_simSWRCT <- function(df,
                           nreps ,
                  p1 ,
                  p2 ,
                  rho)
  {

  
  o2<-p2/(1-p2)
  o1<- p1/(1-p1)
  
  beta <- c(log(o2), log(o1/o2))
  names(beta)<-c("(Intercept)","period")
  
  
  
  sigma2 <-(pi ^ 2) / 3
  theta <-sqrt((rho*sigma2)/(1-rho))
  names(theta)<-c("site.(Intercept)")
  
  ss <- simulate(~ period  + (1 | site), nsim = nreps, family = binomial, 
                 newdata = df, newparams = list(theta = theta,   beta = beta))
  
  
  
  df$resp <- ss[, 1]
  fit1 <- MASS::glmmPQL(resp ~ period + time , random = ~ 1 | site, family = binomial, data=df, verbose = FALSE)
  
  fitsim <- function(i) {
    df$resp <- ss[, i]
    return <- tryCatch({
      res<-coef(summary(MASS::glmmPQL(resp ~ period + time , random = ~ 1 | site, family = binomial, data=df)))["period", ]
      names(res)<-c("est","se","df","t","p")
      return(res)
    },
    warning =function(e) {
      #message(e)  # print error message
      return(c(est=NA, se=NA, df=NA , t=NA, p=NA))
    },
    
    error=function(e) {
      #message(e)  # print error message
      return(c(est=NA, se=NA, df=NA , t=NA, p=NA))
    })
    
    return(return)
  }
  
  fitAll <- ldply(seq(nreps), function(i) fitsim(i))
  
  pow<-with(fitAll, mean( p < 0.05, na.rm = TRUE))
  
  pow

}


pow_simSWRCT(simdat,
                  nreps = 1200,
                  p1 = 0.15,
                  p2 = 0.35,
                  rho= 0.05)

