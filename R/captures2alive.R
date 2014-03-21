captures2alive <- function(scrCaptures, status){
	alive = matrix(1,nrow=max(scrCaptures$individual),ncol=max(scrCaptures$occasion))
	dead.col = scrCaptures$occasion[status=="Dead"]+1
	dead.row = scrCaptures$individual[(status=="Dead")]
	rep.deadr = rep(dead.row[dead.col<(max(scrCaptures$occasion)+1)],times = ((max(scrCaptures$occasion)+1)-dead.col[dead.col<(max(scrCaptures$occasion)+1)]))
	rep.deadc = unlist(lapply(dead.col[dead.col<(max(scrCaptures$occasion)+1)],function(x){seq(x,max(scrCaptures$occasion),1)}))
	alive[cbind(rep.deadr,rep.deadc)]=0
	return(alive)
}