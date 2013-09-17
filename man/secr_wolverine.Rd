\name{secr_wolverine}
\alias{secr_wolverine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
analysis of the wolverine data using secr
}
\description{
analysis using secr %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
secr_wolverine()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##
## The function executes a set of commands to conduct a
## likelihood analysis of the wolverine camera trapping data
## the function returns the results for the 2 x 2 habitat mask
##

## The function is currently defined by the following set of commands:


library("secr")
library("scrbook")
data("wolverine")
traps<-as.matrix(wolverine$wtraps)   #[,1:3]
dimnames(traps)<-list(NULL,c("trapID","x","y",paste("day",1:165,sep="")))

traps1<-as.data.frame(traps[,1:3])
traps1$x<-as.numeric(as.character(traps1$x))
traps1$y<-as.numeric(as.character(traps1$y))

# This seems to ignore the trap operation information
trapfile1<-read.traps(data=traps1,detector="proximity")

trapfile2<-scr2secr(scrtraps=traps,type="proximity")

wolv.dat<-wolverine$wcaps
dimnames(wolv.dat)<-list(NULL,c("Session","ID","Occasion","trapID"))
wolv.dat<-as.data.frame(wolv.dat)
wolvcapt1<-make.capthist(wolv.dat,trapfile1,fmt="trapID",noccasions=165)
wolvcapt2<-make.capthist(wolv.dat,trapfile2,fmt="trapID",noccasions=165)

gr<-(as.matrix(wolverine$grid2))
dimnames(gr)<-list(NULL,c("x","y"))
gr2<-read.mask(data=gr)

gr<-(as.matrix(wolverine$grid4))
dimnames(gr)<-list(NULL,c("x","y"))
gr4<-read.mask(data=gr)

gr<-(as.matrix(wolverine$grid8))
dimnames(gr)<-list(NULL,c("x","y"))
gr8<-read.mask(data=gr)

# run model without trap operation information
wolv.secr<-secr.fit(wolvcapt1,model=list(D~1, g0~1, sigma~1), buffer=20000)

# now use wolvcapt2 which has trap operation information
wolv.secr0<-secr.fit(wolvcapt2,model=list(D~1, g0~1, sigma~1), buffer=20000)
wolv.secr2<-secr.fit(wolvcapt2,model=list(D~1, g0~1, sigma~1), buffer=20000,mask=gr2)
wolv.secr4<-secr.fit(wolvcapt2,model=list(D~1, g0~1, sigma~1), buffer=20000,mask=gr4)
wolv.secr8<-secr.fit(wolvcapt2,model=list(D~1, g0~1, sigma~1), buffer=20000,mask=gr8)
 

# reported in the book chapter:
wolv.secr2

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
