print.scrfit <-
function(x, burn=NULL, ...) {
mcmc<- x$mcmchist
p2k<- x$parms2report
if(is.null(burn))  burn<- max( round(nrow(mcmc)*.05,0), 500)
mcmc<- mcmc[(burn +1):nrow(mcmc),]
smy<-function(v){
     c(mean=mean(v),SD=sqrt(var(v)), quantile(v,c(0.025,0.50,0.975)))
 }
mcmc<-mcmc[,p2k]
tab1<-apply(mcmc,2,smy)
cat("MCMC iterations (total): ", nrow(mcmc),fill=TRUE)
cat(" ",fill=TRUE)
cat(" burn-in: ", burn, fill=TRUE)
cat(" total posterior samples: ", nrow(mcmc), fill=TRUE)
cat("-------------------------------------------------------",fill=TRUE)
cat("Posterior summaries of model parameters:","\n")
print(t(tab1),...)


}
