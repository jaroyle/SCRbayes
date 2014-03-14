print.scrfit <-
function(x, ...) {
mcmc<- x$mcmchist
p2k<- x$parms2report
smy<-function(v){
     c(mean=mean(v),SD=sqrt(var(v)), quantile(v,c(0.025,0.50,0.975)))
 }
mcmc<-mcmc[,p2k]
tab1<-apply(mcmc,2,smy)
cat("Posterior summaries of model parameters:","\n")
print(t(tab1),...)


}
