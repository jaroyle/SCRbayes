###################################################
# SCR.bic calculates the bayesian information criterion (BIC) for a model set based upon 
# results reported with deltaBIC, measures of evidence against a null/reference model, and BIC weights
# Inputs: requires the output from a set of SCRx.fn models passed 
# to the function in a named list.
# I assume that posterior probability density of each parameter 
# is independent... this may require some adjustment, but I have not found 
# an easy way of computing multivariate empirical PDFs. 

SCR.bic = function(modelsin,ncaps, sort=TRUE){
	
	nmodels = length(modelsin)
	model.loglik <- model.npts <- nparam <- rep(NA, nmodels)
	n.nolik <- list()
	for (i in 1:nmodels){
		model = modelsin[[i]]
		model.loglik[i] = max(model$likelihood[,1])
		model.npts[i] = model$likelihood[which(model$likelihood[,1]==max(model$likelihood[,1])),2]
#		hist(model$likelihood, xlab = "Log Likelihood", main = names(modelsin)[i])
		model.post = model$mcmchist[,model$parms2report]
		model.post = model.post[,colnames(model.post)!="sigma"&colnames(model.post)!="sigma2" &colnames(model.post)!="D"&colnames(model.post)!="lam0"]
		nparam[i] = ncol(model.post)
		n.nolik[[i]] <- c("sum"=sum(model$likelihood[,2]), "mean"=mean(model$likelihood[,2]), "median"=median(model$likelihood[,2]), min = min(model$likelihood[,2]), "max"=max(model$likelihood[,2]))
		}
	BIC = -2*model.loglik + nparam*log(ncaps)
	names(BIC) <- names(modelsin)
	deltaBIC = BIC-min(BIC)
	wBIC = exp(-0.5*deltaBIC)/sum(exp(-0.5*deltaBIC))
	
	n.nolik.df = do.call("rbind", n.nolik)
	

	
	out.table = data.frame("Model" = names(modelsin), "MaxlnL" = model.loglik, "BIC"=BIC, "deltaBIC" = deltaBIC, "wBIC" = wBIC,"No.Obs.Unexplained" = model.npts)
	rownames(out.table)=NULL
	if(sort){
		out.table = out.table[order(deltaBIC),]
	}
	out = list("ModelSelection"=out.table, "LowLikelihoodObservations"=n.nolik.df)
	return(out)
	
}