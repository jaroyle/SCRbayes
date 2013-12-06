###################################################
# SCR.bic calculates the bayesian information criterion (BIC) for a model set based upon 
# results reported with deltaBIC, measures of evidence against a null/reference model, and BIC weights
# Inputs: requires the output from a set of SCRx.fn models passed 
# to the function in a named list.
# I assume that posterior probability density of each parameter 
# is independent... this may require some adjustment, but I have not found 
# an easy way of computing multivariate empirical PDFs. 

SCR.bic = function(modelsin,ncaps, refmodel=NULL){
	
	nmodels = length(modelsin)
	model.loglik <- avg.loglik<- min.loglik <- nparam <- model.npts <- rep(NA, nmodels)
	for (i in 1:nmodels){
		model = modelsin[[i]]
		model.loglik[i] = max(model$likelihood[,1])
		model.npts[i] = model$likelihood[which(model$likelihood[,1]==max(model$likelihood[,1])),2]
		avg.loglik[i] = mean(model$likelihood[,1])
		min.loglik[i] = min(model$likelihood[,1])
		hist(model$likelihood, xlab = "Log Likelihood", main = names(modelsin)[i])
		model.post = model$mcmchist[,model$parms2report]
		model.post = model.post[,colnames(model.post)!="sigma"&colnames(model.post)!="sigma2" &colnames(model.post)!="D"]
		nparam[i] = ncol(model.post)
		}
	BIC = -2*model.loglik + nparam*log(ncaps)
	deltaBIC = ifelse(rep(is.null(refmodel),nmodels), BIC-min(BIC), BIC-BIC[refmodel])
	wBIC = exp(-0.5*deltaBIC)/exp(sum(-0.5*deltaBIC))
	
	KR.breaks = c(1,3,20,150)
	KR.evid = c("Supports Null","Barely worth mentioning", "Positive","Strong","Very Strong")
	
	out.table = data.frame("Model" = names(modelsin), "MaxlnL" = model.loglik, "BIC"=BIC, "deltaBIC" = deltaBIC, "wBIC" = wBIC,"No.Obs.Unexplained" = model.npts,"KRevidence" = KR.evid[findInterval(-deltaBIC, KR.breaks)+1])
	out.table$KRevidence = factor(out.table$KRevidence,levels = c(levels(out.table$KRevidence),"Reference Model"))	
	out.table[ifelse(is.null(refmodel), which(BIC==min(BIC)),refmodel),c("KRevidence")]="Reference Model"
	return(out.table)
	
}