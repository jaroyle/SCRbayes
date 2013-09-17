\name{SCRf.fn}
\alias{SCRf.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
conducts spatial capture recapture analysis 
}
\description{
Estimates density of animals from mark-recapture surveys where individual animals were unique identfied
 and spatial locations of captured animals were recorded.Estimation is conducted in a Bayesian framework.  
}
\usage{
SCRf.fn(scrobj,ni = 1100, burn = 100, skip = 2, nz = 200,theta=2, Msigma = 1, Mb = 0, Msex = 0, Msexsigma = 0, Xd = NULL, Xeff = NULL, Xsex = NULL, coord.scale = 5000, thinstatespace = 1, maxNN = 20, dumprate = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{scrobj}{

  a list of type "scrobj" which has elements "traps",
  "captures" and "statespace"

traps:
traps is a matrix containing a seq of numbers (1,2,3....j) where j= the
number of traps in column 1.  Traps can be cameras, hair snares, or
centerpoints of grid locations.  In columns two and three of the trap
matrix are the X and Y or longitude and latitude of the trap locations,
and the survey dates in columns 4 to column 3+k where k=the number of
survey dates. Each row of the matrix will contain either a 1 if the trap
was was open on the survey data or a 0 if the trap was closed.
Alternately, if survey effort was recorded, each row of the matrix will
contain an index of survey effort on each date associated with each trap
j during the specified survey k 

captures: a matrix
where the number of rows=the number of captures.  Each row contains 3
columns labeled LOC_ID, ANIMAL_ID, and SO; where LOC_ID =the number of
the trap where the animal was captured, ANIMAL_ID = the animal's unique
number (must be a sequence 1,2,,,,number of animals), and SO=the
sampling occasion on which the capture took place 

statespace: contains the coordinate of potential centers of activity
including those from outside the trapping area, the statespace is a
matrix where the X or longitudinal coordinate is in the first column and
the Y or latitudinal coordinate is in the second column, and HABITAT= to
1 or 0 is the third column.  The third column indicates whether the
point defined by the X,Y locations is considered habitat (normally)
habitat=1; (see function make.ss) 

}

\item{ni}{ ni=number of iterations }

\item{burn}{ burn=the number of iterations that are considered "burn-in"
  and will be discarded from the final estimation }

\item{skip}{ skip=the
number of iterations in a series that will be retained.  Skip=3 will
retain every third iteration of the MCMC chain.  }

\item{nz}{ nz=the
number of all zero encounter histories; the number of animals to augment
the data with; experimentation with this number is recommended- the
distribution of the parameter estimate of animal density should not be
skewed or "truncated" by the number of augmented animals }

\item{theta}{ the exponent parameter of the distance function. theta = 1
  is exponential; theta=2 is Gaussian; theta = NA then theta is
  estimated.
  }

\item{Msigma}{
Msigma is an indicator variable that =1 if you are estimating the
spatial model, normally Msigma=1
}
  \item{Mb}{
Mb is an indicator variable that=1 if behavioral effects on capture histories are expected.  For example trap-happiness or trap shyness.  Setting Mb=1 is analagous to running
model Mb described by Otis et al. 1978
}
  \item{Msex}{
Msex is an indicator variable that=1 if detection probability is to be estimated seperately by sex, if Msex=1, you must provide a vector Xsex (see Xsex for more details)
}
  \item{Msexsigma}{
Msex sigma is an indicatore variable that = 1 if sigma is to be estimated seperately by sex, if Msexsigma =1 you must provide a vector Xsex
}
  \item{Xd}{
this option is not yet implemented
}
  \item{Xeff}{
Xeff=j by k matrix where k=the number of surveys and j=the number of traps.  For each row indicate the amount of survey effort conducted in that trap location during survey k,
if effort was equal across all surveyed areas Xeff=NULL; 
If a trap x occasion was not sampled then any value can be used for
effort because those elements of Xeff will not be extracted for the
likelihood calculations. (the file "traps" contains information about
which traps x occasions were sampled, see above). 
}
  \item{Xsex}{
Xsex= a vector of length n where n= the number of animals;Xsex[n] is either a zero or 1 indicating the sex of animal [n] in the capture file,
must be provided if Msex or Msexsigma=1
}
  \item{coord.scale}{
coord.scale scales the coordinate system by the input value.  The default is 5000 and will scale all X and Y locations by 5000
}
  \item{thinstatespace}{
thinstate allows the user to retain every yth point in the state space, fewer points can add with computational time, but two few point will bias density estimates 
}
  \item{maxNN}{
used to determine the "neighborhood" for each grid cell; getNN() function in SCRf.fn; information is used in updating the activity centers in the MCMC algorithm
}
  \item{dumprate}{
number of iterations to run before writing a file of results
}
}
\details{
User is advised to read referenced literature, and experimentation with multiple values to determine the appropriate nz, and the size and density of the statespace.
SCRf.fn will provide a summary at the beginning telling the user how many animals were captured and how many were recaptured.  User's should be aware of limitations of 
their data and determine for themselves whether SCR methods are appropriate for their study. 
}
\value{
returns an object with several features.  the function names can be used to identify the features associated with the returned object.
features that can be accessed using $ include, 
out- a matrix containing the values of the parameter estimates for all retained iterations.

G-a matrix containing the scaled coordinates of the statespace, traplocs-a matrix containing the scaled coordinates of the traps, 


Sout-an i by n matrix containing the number corresponding to the row in the statespace that contains the coordinates of the estimated activity center for animal n in iteration i

zout and i by n matrix indicating where the animal was estimated to be a member of the population zout=1 or not zout=0 for each i iteration of the model; zout combined with Sout can be
used to determine the locations of the animals that were included in the population.see functions spatial.plot and image.scale

}
\references{
Otis, D. L., K. P. Burnham, G. C. White, and D. R. Anderson. 1978. Statistical inference from capture
data on closed animal populations. Wildlife Monographs 62.

Royle, J. A., A. J. Magoun, B. Gardner, P. Valkenburg, and R. E. Lowell. 2011. 
Density estimation in a wolverine population using spatial capture-recapture models.  
Journal of Wildlife Management 75:604-611.

Gardner, B., J. A. Royle, and M. T. Wegan. 2009. Hierarchical models for estimating density from DNA mark–recapture studies.
 Ecology 90:1106-1115.


Thompson, C. M., J. A. Royle, J. D. Gardner. in press. A framework for inference about carnivore density from unstructured spatial sampling of scat using detector dogs.
  Journal of Wildlife Management.
}
\author{
Andy Royle
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{



data("lions")
scrobj<-list(traps=traps.lions,captures=captures.lions,statespace=statespace.lions)

test<-

SCRf.fn(scrobj,ni=220, burn=20,
 skip=2,nz=200,Msigma=1, Mb=0, Msex=0, Msexsigma=0, Xd=NULL,
 Xeff=NULL, Xsex=NULL, thinstatespace=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line




