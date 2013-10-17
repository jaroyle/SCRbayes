SCR.area = function(obj, SO, Mb=0, Mbvalue=NULL, Xeff=NULL, Xd=NULL, Xsex = NULL, iter=NULL, scalein=1000,scaleout=100000, useSnowfall=FALSE, nprocs = 1, con.type="SOCK"){
	# Some basic error checking on behavioral covariates
	# The expectations on this input likely need some further thought. 
	# I allow either a global or trap-specific value	
	if (Mb==1){ 
		if(!is.numeric(Mbvalue)) {
			cat("Error: Please supply a numeric value for behavioral capture effect", fill=T)
			return()}
		if (!(length(Mbvalue) %in% c(1,nrow(obj$G)))){
			cat("Error: Behavioral effect should take on a single value or a trap-specific value", fill=T)
			return()
		}
	}
	# Compute activity centers for each iteration of the chain
	# This takes some time and it may be better to incorporate this 
	# into the if statements below for computing param.values to improve calculation time
	centers= (obj$Sout*obj$zout)
	gridchain = matrix(nrow=nrow(centers),ncol=nrow(obj$statespace))
	
	for (k in 1:nrow(centers)){
		for (j in 1:nrow(obj$statespace)){
			gridchain[k,j] = sum(centers[k,]==j)
		}
	}
	#n.animals=apply(gridchain,2,sum)
	
	# if iter is a character defining a functioning passed to apply...
	if (is.character(iter)&iter!="all"){
		param.values = apply(obj$mcmchist,2,iter)
		param.values = array(param.values,c(ifelse(is.matrix(param.values),nrow(param.values),1),ncol(obj$mcmchist)))
		colnames(param.values) = colnames(obj$mcmchist)
		# These lines select the location of the activity centers from the iteration
		# that is closest to the desired total number of animals
		# In case of a tie, the selection is made randomly from the equivalent candidates
		Nsuper.iter = c(apply(matrix(rowSums(gridchain),ncol=1),2,iter))
		n.animals.iter = matrix(NA, nrow=nrow(param.values), ncol = ncol(gridchain))
		for (N in 1:length(Nsuper.iter)){
			n.animals.iter[N,] = gridchain[which(abs(rowSums(gridchain)-Nsuper.iter[N])==min(abs(rowSums(gridchain)-Nsuper.iter[N])))[sample(sum(abs(rowSums(gridchain)-Nsuper.iter[N])==min(abs(rowSums(gridchain)-Nsuper.iter[N]))), 1)],]
		}	
		n.animals.iter = array(n.animals.iter,c(ifelse(is.matrix(n.animals.iter), nrow(n.animals.iter),1),ncol(gridchain)))
	}
	
	# If iter is a numeric giving which iterations of the chain to work with...
	if (is.numeric(iter)){
		param.values = obj$mcmchist[iter,]
		n.animals.iter = gridchain[,iter]
	}
	# Finally, iter can specify all values from the MCMC chain
	# I use this as the null/default option
	if (iter =="all"|is.null(iter)){
		param.values - obj$mcmchist
		n.animals.iter = gridchain
	}
	# String together necessary info to be passed to internal function (below)
	inputs = list("obj"=obj,"param.values"=param.values, "n.animals.iter"=n.animals.iter,"SO"=SO, "Xeff" = Xeff, "Xsex"=Xsex,"Mb"=Mb, "Mbvalue"=Mbvalue, "Xd"=Xd)
	
	# postProb calculates the effective areas, population sizes and densities 
	# for given parameter values, locations of activity centers, and covariates
	postProb = function(i,inputs){
		obj = inputs$obj
		param.values = inputs$param.values
		n.animals.iter = inputs$n.animals.iter
		SO = inputs$SO
		Xeff = inputs$Xeff
		Xsex = inputs$Xsex
		Mb = inputs$Mb
		Mbvalue = inputs$Mbvalue
		Xd =inputs$Xd
	
	# Calculate sigma and baseline detection probability
	# Use weighted mean when no/NULL Xsex specified  
	if (is.null(Xsex)){
	sigma <- weighted.mean(c(param.values[i,1],param.values[i,3]),c(1-param.values[i,"psi.sex"],param.values[i,"psi.sex"]))
	# Calculate baseline detection probability
	# Again, I use a weighted avg for sex-based differences
	lam0 = exp(log(param.values[i,"lam0"]) + param.values[i,"psi.sex"]*param.values[i,"beta.sex"])
	} else { # When Xsex is specified calculate exposure probability for given sex
	sigma <- c(param.values[i,1],param.values[i,3])[Xsex+1]
	# Calculate baseline detection probability
	# Again, I use a weighted avg for sex-based differences
	lam0 = exp(log(param.values[i,"lam0"]) + Xsex*param.values[i,"beta.sex"])	
	}
	# Reconstruct the trap locations on the original scale
	ss.x <- unique(obj$statespace$X_coord)
	ss.y <- rev(unique(obj$statespace$Y_coord))
	coord.scale <- ((obj$Gunscaled[,1]-min(obj$Gunscaled[,1]))/(obj$G[,1]-min(obj$G[,1])))[nrow(obj$Gunscaled)]
	traplocs.utm = data.frame(t((t(as.matrix(obj$traplocs))-apply(obj$G,2,min))*coord.scale + apply(obj$Gunscaled,2,min)))
	
	# Create matrix to store the probability of capture in any trap given the 
	# location of activity center
	prob.anycap <- matrix(NA, nrow=length(ss.y), ncol = length(ss.x))
		
		for (x in ss.x){ # For each x coord of statespace...
		  for (y in ss.y){ # For each y coord of statespace...
			# Distance to all traps		
		    temp.dist <- sqrt((traplocs.utm[,1]-x)^2+(traplocs.utm[,2]-y)^2) # Find distance from central node
		    # Compute the linear predictor components
		    lp.eff = ifelse(is.null(Xeff),0, param.values[i,"beta1(effort)"]*Xeff)
		    lp.behav = ifelse(Mb==0,0, Mb*Mbvalue*param.values[i,"beta.behave"])
		    lp.dens = ifelse(is.null(Xd),0,(Xd[obj$statespace$X_coord==x&obj$statespace$Y_coord==y]*param.values[i,"beta.density"]/(log(sum(exp(Xd*param.values[i,"beta.density"]))))))
		    # Add all the lps, convert to probability of capture in a trap
		    prob.cap <-1 - exp(-lam0*exp(-sigma*(temp.dist/coord.scale)^(2*param.values[i,"theta"])+lp.eff+lp.behav+lp.dens ))
		    # Store probability of capture in ANY trap 
		    prob.anycap[which(ss.y==y), which(ss.x==x)] <- 1- prod(1-prob.cap)^SO
		  }
		}
		# Find Resolution of statespace
		ss.res = c(as.numeric(names(table(diff(ss.x)))),as.numeric(names(table(diff(rev(ss.y))))))
		# Calculate density scaling factor
		# Effective area
		area.iter <- sum(prob.anycap)*(prod(ss.res/scalein))
		popcontrib.iter <- matrix(n.animals.iter[i,],nrow = length(ss.y),ncol=length(ss.x),byrow=F)*apply(prob.anycap,2,rev)
		# Effective number of exposed individuals
		ngrid.iter <- sum(popcontrib.iter)
		# Compute density, scaling to /100 km^2
		density.iter <- ngrid.iter/area.iter *(scaleout/scalein)
		out1 = list("chain"=c(param.values[i,], "area"=area.iter,"Ngrid"=ngrid.iter,"density"=density.iter),"prob.anycap"=prob.anycap)
		return(out1)
	}

# If using snowfall for multicore computations, useful for many interations
	if (useSnowfall){
		require(snowfall)
		sfInit(parallel=TRUE, cpus=nprocs,type=con.type)
		sfExport("postProb")
		
		out1 = sfLapply(1:nrow(param.values),postProb,inputs)
		sfStop()
	} else { # Otherwise series computations using lapply
		out1 = lapply(1:nrow(param.values),postProb,inputs)
	}

# Re-arrange output
	params = t(sapply(out1,"[[",1) )
	probs = sapply(out1,"[[",2)
# Return the updated mcmchist and 
# the corresponding exposure probabilities (can be mapped)	
	return(list("mcmchist" = params, "prob.cap"=probs))
	
}