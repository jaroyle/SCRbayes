make.ss <-
function(traps, buffer, grid.density)
{y<-seq(min(traps$x)-buffer, max(traps$x)+buffer, by=grid.density)
x<-seq(min(traps$y)-buffer, max(traps$y)+buffer, by=grid.density)
ss.y<-rep(y, length(x))
ss.x<-rep(x[1], length(y))

for (i in 2:length(x))
{
ss.x<-c(ss.x,rep(x[i],length(y)))
}

ss.animal<-data.frame(X=ss.x, Y_coord=ss.y)
list(ss.animal=ss.animal)
}
