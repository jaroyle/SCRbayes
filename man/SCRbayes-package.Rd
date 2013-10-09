\name{SCRbayes-package}
\alias{SCRbayes-package}
\alias{SCRbayes}
\docType{package}
\title{
  Bayesian analysis of spatial capture-recapture models.
}
\description{
  Bayesian analysis of spatial capture-recapture models with a discrete
  state-space. The aim of this package is to provide more efficient
  analysis than WinBUGS/JAGS and the R package "secr",
  at least for a limited class of models. 
}
\details{
\tabular{ll}{
Package: \tab SCRbayes\cr
Type: \tab Package\cr
Version: \tab xxxx\cr
Date: \tab  xxxxx\cr
License: \tab What license is it under?\cr
}

The main fitting function is SCRh.fn -- see below for an example and
look at the help file after reading this document. 

QUICK START GUIDE. 

This is the helpfile ?SCRbayes which is "how to get started". If you
have some SCR data laying around and you NEED AN ANSWER FAST(!!!!), then read
this document. We're here to help. 

Motivation: SCRbayes exists to do Bayesian analysis of a certain class
of SCR models (to be described shortly) when the state-space is defined
by a DISCRETE set of points. In secr terminology the discrete
state-space grid is called a "habitat mask".

Why should you use SCRbayes over WinBUGS or JAGS or the R package secr
(which does likelihood analysis)?  (1) The main reason is that SCRbayes
is relatively more robust (in an operational sense) to the "dimension"
of the problem. If you have a lot of traps, sample occassions, or a
large state-space grid, SCRbayes is less likely to run out of memory or
take forever to produce estimates.  (2) SCRbayes has a limited number of
features that secr doesn't have. It handles unknown sex information in
models with sex-specific parameters by putting a prior distribution on
the individual covariate "sex".  (3) In my view, SCRbayes packing and
handling of the data is easier and more intuitive than secr, although
secr handles much more general models.  (4) SCRbayes produces the
posterior distribution of the underlying point pattern of individual
activity centers. This is enormously useful for producing density maps,
or estimates of population size in any arbitrary area. Note: summary
functions for doing all of these things are in development so you have
to do these things by yourself right now.
Note: BUGS/JAGS and secr all have a host of advantages over SCRbayes but
this package isn't about them. 

Discrete state-space: SCRbayes works strictly with a DISCRETE
state-space grid. IF you can define your state-space by a rectangle or a
constrained rectangle then JAGS or WinBUGS can be extremely fast working
with a continuous state-space.  Note: secr can be very efficient too if
you do not define the state-space grid ("habitat mask" in secr
terminology) or if you define it to be relatively coarse. Its only when
the state-space grid becomes enormous that secr starts to break down
[note: this is subject to change we suppose].

Future developments:
(1) Multiple core capabilities will be developed. 
(2) Use the "dumprate" argument (not yet implemented) to get
intermediate output for really huge problems.
(3) Functions to assess convergence of the MCMC algorithm and how long
the MCMC algorithm should be run will be implemented soon. 
(4) Multi-session models: Right now SCRbayes does not handle multi-year,
multi-site or multi-session models in any generality. However, certain
specific cases can be analyzed by "stacking" the data sets. Details to
follow. 


Question: I NEED TO DO AN ANALYSIS OF SOME SCR DATA NOW GODDAMMIT, CAN YOU JUST
CUT THE CRAP AND TELL ME WHAT TO DO?
Answer: Don't worry, we get this question all the time, we're here to help you.

1. You need some trap locations, places where you collected data. These
should be in some Cartesian coordinate system. Typically UTMs.  They
should exist as a MATRIX in R having dimensions "ntraps x 2" where the
columns are (x,y), (easting, northing), (UTME, UTMN) or whatever you
call them. These should be real numbers so you can plot() them.

From the trap coordinate information you
need to produce the "trap deployment file", abbreviated "TDF". 
For more information see ?traps or ?TDF (not written yet)
%XXXXX need as.TDF
%XXXXX class "TDF"
If your study was operated for K sampling occasions then you need to
create a matrix that is "ntraps x K" defining which traps were
operational in which periods. This matrix should have entries 1 if the
trap was operational and 0 if not. 
Using these two pieces of information you can create the TDF easily by
doing this:

newtraps<- cbind(traplocs, trapop)
%my.TDF<- as.TDF(newtraps)
%This 2nd line just adds an attribute declaring the data to be a "TDF"
%data structure. IF you don't have names for your traps it also makes up
%some fake names so there is an extra leading column on the result. 


2. To use SCRbayes you need to define the state-space as a grid of coordinates. What, you
ask, is the "state-space"?
The state-space defines the possible locations of individual activity
centers. It is, in effect, a prior distribution for
"where a guy might live" in the vicinity of your study area. 
The state-space grid should be a grid of regularly-spaced points, it is
a matrix of real numbers having dimensions "nG x 2" where the columns
are (x,y), (easting, northing), etc... 

You need to add a 3rd column which defines whether each point is
available habitat (1) or not (0). We recommend defining everything to
be habitat even if you think it isn't.  The reason is this: Even if a
point corresponds to a parking lot somewhere in an urban area, if the
individual is using the space around the parking lot, then its center of
activity might well be on the parking lot even if it is subjectively
believed to be "non-habitat" by a biologist. 

You should convert your statespace grid to an object of class
"statespace" using the as.statespace() function.
XXXXX NOT YET IMPLEMENTED XXXXXXXX

If you don't have a state-space grid you can easily make one using the
make.statespace function. See ?statespace or ?make.statespace
XXXXX NOT FLESHED-OUT YET XXXXXXXX 

3. 
You need some encounter data from an actual study that you did which is
packaged into a data file called the encounter data files (EDF). 
This has exactly 4 columns, in any order, labeled: session, individual,
occasion and trapid.
These should all be integer data types where session is a "primary
period" (season or year), individual is the integer identity of
individuals in the sample (1 to n), occasion is the integer identify of
sample period and trapid is the integer trap identity 1:ntraps. Trapid
must correspond to ROWS of the trap deployment file (TDF). 


4. If you have these things you can do an analysis easily. See
below. See also the help file ?wolverine

}
\author{
Andy Royle  and Robin Russell  and Chris Sutherland and OTHERS to be identified

Maintainer: Andy Royle, aroyle@usgs.gov

SIGN UP FOR THIS GOOGLE EMAIL LIST: 
}
\references{
The code originates from Royle et al. (2009) which was developed a
little bit and given a GUI and called SPACECAP (Gopalaswamy et
al. 2011???). In parrallel some additional things were done to the code
and this produced the primiative SCRbayes package from Russell et
al. (2012) from which the current SCRbayes package arose and is under
active development. 

You should see the SCR book Royle et al. (2013) for general concepts and
  methods. The R package "scrbook" is not directly related to the
  SCRbayes package although some functions and data sets are shared
  between the two.  

The R package secr (Efford 2011) fits a lot of different SCR models by
likelihood methods. 

}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{

library("SCRbayes")
data(lions)

### This makes the encounter data file "EDF" for the mountain lion data
### from Russell et al. (2012)
### It is best to use exactly the column names given although the
### function should work without those names. 

dimnames(captures.lions)<-
list(1:nrow(captures.lions),c("trapid","individual","occasion"))

### Now we package up the trap information, the EDF and the state-space
### information into a single data structure of class "scrdata"

lion.scrdata <- scrData(traps.lions,captures.lions,statespace.lions)

### Now we do a very short MCMC run to see how things look

test <- SCRh.fn(lion.scrdata,ni=220, burn=20,
 skip=2,nz=200,Msigma=1, Mb=0, Msex=0, Msexsigma=0, 
 Xeff=NULL, Xsex=NULL, thinstatespace=2)

### summarize the output

print(test, digits = 2)


}
