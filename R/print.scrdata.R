print.scrdata <-
function(x, ...) {
cat("SCR Data Summary ----------------",fill=TRUE)
ntraps<-nrow(x$traps)
nx<-table(x$captures[,"individual"])
cat("\nTotal encounters of each individual:","\n")
print(nx,...)
cat("\nSummary encounter frequencies:","\n")
print(table(nx),...)
a<-table(x$captures[,"individual"],x$captures[,"trapid"])
a<-table(apply(a>0,1,sum))
cat("\nNum. unique traps of capture:","\n")
print(a,...)

if(!is.null(scrobj$Ytel)){
cat("Telemetry Data Summary  ----------------",fill=TRUE)

    cat("\nTelemetry data on ", length(unique(scrobj$Ytel[,3])), " individuals.","\n")




}




}
