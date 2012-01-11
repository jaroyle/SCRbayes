\name{array3d2SCR.fn}
\alias{array3d2SCR.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
utility function to convert 3-d encounter history array to SCR format 
}
\description{
%%  ~~ array3d2SCR.fn is used to transform a 3 dimensional array with dimensions equal to 
	individuals x sampling occasion x traps into a 2 dimensional matrix where column 1=traps, column 2=individuals
	and column 3=sample occasion which can passed into the function SCRf.fn~~
}
\usage{
array3d2SCR.fn(y = object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~object = three dimensional array with dimensions (i,k,j) where 
	i=number of individuals, k=number of sampling occasion, and j=number of traps ~~
}
}
\details{
A typical object is constructed of an array containing zero and 1 variables
indicating if individual i was observed on sampling occasion k in trap j
}
\value{
matrix of nrows and 3 columns is returned where columns contain trap number, individual, and sampling occasions
 of observed animals. Object can be used as a input for the function SCRf.fn
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Andy Royle aroyle@usgs.gov
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~See also SCRf.fn, e2dist, make.ss,expit, image.scale,SCR.gof , ~~~
}
\examples{
##---- test.data<-array(dim=c(200,3,50)) 
##-- for(i in 1:200)
	{
	for(j in 1:3)
	{
	for(k in 1:50)
	{
	test.data[i,j,k]<-rbinom(1,1,0.2)
	}
	}
	}
SCR.data<-array3d2SCR.fn(test.data)

##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y = tiger.data3d) 
{
    dd <- dim(y)
    nind <- dd[1]
    nrep <- dd[2]
    ntrap <- dd[3]
    indid <- repid <- trapid <- NULL
    for (i in 1:nind) {
        yi <- y[i, , ]
        indid <- c(indid, rep(i, sum(yi > 0, na.rm = TRUE)))
        repid <- c(repid, row(yi)[!is.na(yi) & yi > 0])
        trapid <- c(trapid, col(yi)[!is.na(yi) & yi > 0])
    }
    o <- cbind(trapid, indid, repid)
    dimnames(o) <- list(NULL, c("trapid", "individual", "sample"))
    o <- o[order(o[, 1]), ]
    o
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% 
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

