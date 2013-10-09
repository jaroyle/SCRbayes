###################################################
# SCR.bf calculates bayes factors for a model set based upon 
# Gelfand and Dey (1994) approximation using the posterior samples.
# Inputs: requires the output from a set of SCRx.fn models passed 
# passed to the function in a named list

SCR.bf = function(modelsin, refmodel=NULL){
	require(mvtnorm)
	nmodels = length(modelsin)
	model.lik <- maxL <-avgL<- medL<- rep(NA, nmodels)
	for (i in 1:nmodels){
		model = modelsin[[i]]
		model.post = model$out[,model$parms2report]
		model.post = model.post[,colnames(model.post)!="sigma"&colnames(model.post)!="sigma2" &colnames(model.post)!="D"]
		
		post.avg = apply(model.post,2,mean)
		post.var = var(model.post)
		emp.dens = apply(model.post,2,density)
		emp.dens = lapply(emp.dens,function(x){
			prob<-x$y/sum(x$y)
			return(cbind(x=x$x,prob))})
		chain.dens = matrix(NA, nrow=nrow(model.post), ncol=ncol(model.post))
		for (c in 1:ncol(model.post)){
			chain.dens[,c]=emp.dens[[c]][,2][findInterval(model.post[,c],emp.dens[[c]][,1])]
		}
		f.samples = dmvt(model.post,delta=post.avg, sigma=cov2cor(post.var), df=nrow(model.post),log=F)
		model.lik[i] = (1/nrow(model.post)*sum(f.samples/(model$Lout*apply(chain.dens,1,prod))))^(-1)
		maxL[i] = max(model$Lout)
		avgL[i] = mean(model$Lout) 
		medL[i] = median(model$Lout)
	}
	post.prob = model.lik/sum(model.lik)
	names(model.lik) <- names(modelsin)
	if(!is.null(refmodel)){
	bf = model.lik/model.lik[refmodel]
	Jeff.breaks = c(1,3,10,30,100)
	Jeff.evid = c("Supports Null", "Barely worth mentioning", "Substantial","Strong","Very Strong","Decisive")
	KR.breaks = c(1,3,20,150)
	KR.evid = c("Supports Null","Barely worth mentioning", "Positive","Strong","Very Strong")
	out.table = data.frame("Model" = names(modelsin), "Post.Prob" = post.prob, "BayesFactor"=bf, "log10BF" = log10(bf), "Two*lnBF" = 2*log(bf), "JeffreysEvidence" = Jeff.evid[findInterval(bf, Jeff.breaks)+1],"KRevidence" = KR.evid[findInterval(bf, KR.breaks)+1])
	out.table$JeffreysEvidence = factor(out.table$JeffreysEvidence,levels = c(levels(out.table$JeffreysEvidence),"Null Model"))	
	out.table$KRevidence = factor(out.table$KRevidence,levels = c(levels(out.table$KRevidence),"Null Model"))	
	out.table[refmodel,c("JeffreysEvidence","KRevidence")]="Null Model"
	return(out.table)
	} else {
		out.table = data.frame("Model" = names(modelsin), "Post.Prob" = post.prob)
		return (out.table)
	}
}