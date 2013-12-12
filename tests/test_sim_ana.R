library(spatsurv)

set.seed(10)

dat <- simsurv(mcmc.control=mcmcpars(nits=100,burn=10,thin=11),dist="exp")

X <- as.data.frame(dat$X) # covariates

survtimes <- dat$survtimes
n <- length(survtimes)
censtimes <- runif(n,min(survtimes),max(survtimes))                                    
survdat <- gensens(survtimes,censtimes)                                    

plot(survfit(survdat~1))

su <- survexpon(formula=survdat~age+sex+cancer,
                data=X,
                mcmc.control=mcmcpars(nits=1000,burn=100,thin=9),
                betapriorsd=10,
                omegapriorsd=10)
                
dat1 <- simsurv(mcmc.control=mcmcpars(nits=100,burn=10,thin=11),dist="weibull",theta=c(1,0.5))                