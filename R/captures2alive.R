captures2alive <- function(scrCaptures, status, n.occasion){
	alive = matrix(1,nrow=max(scrCaptures$individual),ncol=n.occasion)
	dead.col = scrCaptures$occasion[status=="Dead"]+1
	dead.row = scrCaptures$individual[(status=="Dead")]
	rep.deadr = rep(dead.row[dead.col<(n.occasion+1)],times = ((n.occasion+1)-dead.col[dead.col<(n.occasion+1)]))
	rep.deadc = unlist(lapply(dead.col[dead.col<(n.occasion+1)],function(x){seq(x,n.occasion,1)}))
	alive[cbind(rep.deadr,rep.deadc)]=0
	return(alive)
}