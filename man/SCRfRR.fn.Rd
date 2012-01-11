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
SCRf.fn(ni = 1100, burn = 100, skip = 2, nz = 200, traps = traps1, captures = captures1, statespace = grid900, Msigma = 1, Mb = 0, Msex = 0, Msexsigma = 0, Xd = NULL, Xeff = NULL, Xsex = NULL, coord.scale = 5000, thinstatespace = 1, maxNN = 20, dumprate = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ni}{
ni=number of iterations
}
  \item{burn}{
burn=the number of iterations that are considered "burn-in" and will be discarded from the final estimation
}
  \item{skip}{
skip=the number of iterations in a series that will be retained.  Skip=3 will retain every third iteration of the MCMC chain.
}
  \item{nz}{
nz=the number of all zero encounter histories; the number of animals to augment the data with; experimentation with this number is recommended-
the distribution of the parameter estimate of animal density should not be skewed or "truncated" by the number of augmented animals 
}
  \item{traps}{
traps is a matrix containing a seq of numbers (1,2,3....j) where j= the number of traps in column 1.  Traps can be cameras, hair snares, or centerpoints of grid locations.
In columns two and three of the trap matrix are the X and Y or longitude and latitude of the trap locations, and the survey dates in columns 4 to column 3+k where
 k=the number of survey dates. Each row of the matrix will contain either a 1 if the trap was was open on the survey data or a 0 if the trap was closed.  Alternately, if survey 
effort was recorded, each row of the matrix will contain an index of survey effort on each date associated with each trap j during the specified survey k  
}
  \item{captures}{
captures is a matrix where the number of rows=the number of captures.  Each row contains 3 columns labeled LOC_ID, ANIMAL_ID, and SO; where LOC_ID =the number of the trap
where the animal was captured, ANIMAL_ID = the animal's unique number (must be a sequence 1,2,,,,number of animals), and SO=the sampling occasion on which the capture took
place
}
  \item{statespace}{
statespace= contains the coordinate of potential centers of activity including those from outside the trapping area,
the statespace is a matrix where the X or longitudinal coordinate is in the first column and the Y or latitudinal coordinate is in the second column, and HABITAT= to 1 or 0 is the third column.
The third column indicates whether the point defined by the X,Y locations is considered habitat (normally) habitat=1; (see function make.ss)
}
  \item{Msigma}{
Msigma is an indicator variable that =1 if you are estimating the spatial model, normally Msigma=1
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
J. Andrew Royle
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ni = 1100, burn = 100, skip = 2, nz = 200, traps = traps1, 
    captures = captures1, statespace = grid900, Msigma = 1, Mb = 0, 
    Msex = 0, Msexsigma = 0, Xd = NULL, Xeff = NULL, Xsex = NULL, 
    coord.scale = 5000, thinstatespace = 1, maxNN = 20, dumprate = 1000) 
{
    Y <- captures
    traplocs <- traps[, 2:3]
    MASK <- as.matrix(traps[, 4:ncol(traps)])
    nind <- max(Y[, 2])
    T <- dim(MASK)[2]
    M <- nind + nz
    ntraps <- nrow(traplocs)
    checkdata.fn <- function(Y) {
        encfreqs <- table(table(Y[, 2]))
        totrecaps <- sum(encfreqs * (as.numeric(names(encfreqs)) - 
            1))
        trapfreqs <- table(apply(table(Y[, 2], Y[, 1]) > 0, 1, 
            sum))
        multitraps <- sum(trapfreqs * (as.numeric(names(trapfreqs)) - 
            1))
        c(totrecaps, multitraps)
    }
    smydata <- checkdata.fn(Y)
    cat("You had ", smydata[1], " recaptures in your study", 
        fill = TRUE)
    cat("You had ", smydata[2], " multi-trap recaptures in your study", 
        fill = TRUE)
    if ((smydata[1] < 10) | (smydata[2] < 5)) {
        cat("You should go into the field and obtain more data.", 
            fill = TRUE)
        cat("Try again later.", fill = TRUE)
        return(NULL)
    }
    statespace <- statespace[seq(1, nrow(statespace), thinstatespace), 
        ]
    goodbad <- statespace[, 3]
    G <- statespace[, 1:2]
    G <- G[goodbad == 1, ]
    Gunscaled <- G
    nG <- nrow(G)
    if (is.null(Xd)) {
        Xd <- rep(1, nG)
        beta.den <- 0
    }
    mgx <- min(traplocs[, 1])
    mgy <- min(traplocs[, 2])
    traplocs[, 1] <- (traplocs[, 1] - min(traplocs[, 1]))/coord.scale
    traplocs[, 2] <- (traplocs[, 2] - min(traplocs[, 2]))/coord.scale
    G[, 1] <- (G[, 1] - mgx)/coord.scale
    G[, 2] <- (G[, 2] - mgy)/coord.scale
    msk2 <- array(NA, c(nind + nz, T, ntraps))
    for (i in 1:(nind + nz)) {
        msk2[i, 1:T, 1:ntraps] <- t(MASK[1:ntraps, 1:T])
    }
    msk2 <- as.vector(msk2)
    Ynew <- array(0, dim = c(nind, T, ntraps))
    Ynew[cbind(Y[, 2], Y[, 3], Y[, 1])] <- 1
    Y <- Ynew
    Yaug <- array(0, dim = c(nind + nz, T, ntraps))
    for (j in 1:nind) {
        Yaug[j, 1:T, 1:ntraps] <- Y[j, 1:T, 1:ntraps]
    }
    if (!is.null(Xeff)) {
        Xeffnew <- array(0, dim = c(nind + nz, T, ntraps))
        for (j in 1:M) {
            Xeffnew[j, 1:T, 1:ntraps] <- t(Xeff)
        }
        Xeff.tf <- TRUE
    }
    if (is.null(Xeff)) {
        Xeffnew <- array(0, dim = c(nind + nz, T, ntraps))
        Xeff.tf <- FALSE
    }
    Xeff <- Xeffnew
    if (!is.null(Xsex)) {
        Xsexnew <- c(Xsex, rep(NA, nz))
    }
    if (is.null(Xsex)) {
        Xsexnew <- rep(0, nind + nz)
    }
    Xsex <- Xsexnew
    sex.naflag <- is.na(Xsex)
    prevcap <- array(0, c(nind + nz, T, ntraps))
    for (i in 1:(nind)) {
        for (j in 1:ntraps) {
            tmp <- Yaug[i, 1:T, j]
            if (any(tmp == 1)) {
                fst <- min((1:T)[tmp == 1])
                if (fst < T) 
                  prevcap[i, (fst + 1):T, j] <- 1
            }
        }
    }
    prevcap <- as.vector(prevcap)
    arr.trues <- array(TRUE, c(nind + nz, T, ntraps))
    idx <- which(arr.trues, arr.ind = TRUE)
    y <- as.vector(Yaug)
    y <- y[msk2 == 1]
    Xeff <- as.vector(Xeff)
    Xeff <- Xeff[msk2 == 1]
    prevcap <- prevcap[msk2 == 1]
    indid <- idx[msk2 == 1, 1]
    repid <- idx[msk2 == 1, 2]
    trapid <- idx[msk2 == 1, 3]
    getNN <- function(maxNN, G) {
        nG <- nrow(G)
        NN <- matrix(NA, nrow = nG, ncol = maxNN)
        for (i in 1:nG) {
            od <- sqrt((G[i, 1] - G[, 1])^2 + (G[i, 2] - G[, 
                2])^2)
            NN[i, 1:maxNN] <- order(od)[1:maxNN]
        }
        numnn <- rep(0, nrow(NN))
        for (i in 1:nG) {
            for (j in 1:ncol(NN)) {
                if (any(NN[NN[i, j], ] == i, na.rm = TRUE)) 
                  next
                else NN[i, j] <- NA
            }
            numnn[i] <- sum(!is.na(NN[i, ]))
            NN[i, 1:numnn[i]] <- NN[i, !is.na(NN[i, ])]
        }
        if (min(numnn) < 3) 
            cat("State-space grid has isolated or nearly-isolated cells increase maxNN or modify state-space", 
                fill = TRUE)
        out <- list(NN = NN, numnn = numnn)
        return(out)
    }
    hld <- getNN(maxNN, G)
    NN <- hld$NN
    numnn <- hld$numnn
    centers1 <- rep(NA, nind)
    for (i in 1:nind) {
        tt <- t(as.matrix(Yaug[i, , ]))
        tt <- col(tt)[tt == 1]
        xxx <- traplocs[tt, ]
        av.coord <- colSums(xxx)/nrow(xxx)
        dvec <- as.vector(e2dist(matrix(av.coord, ncol = 2), 
            G))
        centers1[i] <- (1:length(dvec))[dvec == min(dvec)][1]
    }
    centers2 <- sample(1:nG, M - nind, replace = TRUE)
    centers <- c(centers1, centers2)
    S <- G[centers, ]
    if (Msexsigma == 0) 
        bsigma <- 1
    if (Msexsigma == 1) 
        bsigma <- c(3, 3)
    loglam0 <- log(0.018)
    beta <- 0
    beta.sex <- 0
    beta1 <- 0
    lam0 <- exp(loglam0)
    psi <- 0.5
    psi.sex <- mean(Xsex, na.rm = TRUE)
    z <- c(rep(1, nind), rbinom(nz, 1, psi))
    if (sum(sex.naflag) > 0) 
        Xsex[sex.naflag] <- rbinom(sum(sex.naflag), 1, 0.65)
    lik.fn <- function(lpc, y1) {
        llvector.new <- -1 * exp(lpc)
        part2 <- exp(exp(lpc[y1])) - 1
        part2[part2 == 0] <- min(part2[part2 != 0])
        llvector.new[y1] <- llvector.new[y1] + log(part2)
        llvector.new
    }
    trapgridbig <- traplocs[trapid, ]
    y1 <- y == 1
    c1 <- (S[indid, 1] - trapgridbig[, 1])^2
    c2 <- (S[indid, 2] - trapgridbig[, 2])^2
    gof.new <- gof.data <- rep(NA, (ni - burn)/skip)
    out <- matrix(NA, nrow = (ni - burn)/skip, ncol = 12)
    dimnames(out) <- list(NULL, c("bsigma", "sigma", "bsigma2", 
        "sigma2", "lam0", "beta", "beta1(effort)", "beta.sex", 
        "psi", "psi.sex", "Nsuper", "Nin"))
    zout <- matrix(NA, nrow = (ni - burn)/skip, ncol = M)
    Sout <- matrix(NA, nrow = (ni - burn)/skip, ncol = M)
    m <- 1
    LM1 <- LM2 <- matrix(NA, nrow = M, ncol = length(indid)/M)
    ones <- rep(1, ncol(LM1))
    if (Msexsigma == 0) 
        lp.sigma <- Msigma * bsigma * (c1 + c2)
    if (Msexsigma == 1) 
        lp.sigma <- bsigma[Xsex[indid] + 1] * (c1 + c2)
    acc.count <- 0
    delta <- 0.05
    for (i in 1:ni) {
        cat("iter: ", i, fill = TRUE)
        lp <- loglam0 + Mb * beta * prevcap - lp.sigma + beta1 * 
            Xeff + Msex * beta.sex * Xsex[indid]
        loglam0c <- rnorm(1, loglam0, 0.1)
        lpc <- loglam0c + Mb * beta * prevcap - lp.sigma + beta1 * 
            Xeff + Msex * beta.sex * Xsex[indid]
        llvector <- lik.fn(lp, y1)
        llvector.new <- lik.fn(lpc, y1)
        LM1[1:length(LM1)] <- llvector.new
        LM2[1:length(LM2)] <- llvector
        if (runif(1) < exp(sum(((LM1[z == 1, ] - LM2[z == 1, 
            ]) \%*\% ones)))) {
            loglam0 <- loglam0c
            lam0 <- exp(loglam0)
            llvector <- llvector.new
            lp <- lpc
            LM2 <- LM1
        }
        if (Msexsigma == 0) {
            bsigmac <- exp(rnorm(1, log(bsigma), delta))
            lp.sigmac <- Msigma * bsigmac * (c1 + c2)
            lpc <- loglam0 + Mb * beta * prevcap - lp.sigmac + 
                beta1 * Xeff + Msex * beta.sex * Xsex[indid]
            llvector.new <- lik.fn(lpc, y1)
            LM1[1:length(LM1)] <- llvector.new
            if (runif(1) < exp(sum(((LM1[z == 1, ] - LM2[z == 
                1, ]) \%*\% ones)))) {
                lp.sigma <- lp.sigmac
                bsigma <- bsigmac
                llvector <- llvector.new
                lp <- lpc
                LM2 <- LM1
                acc.count <- acc.count - 1
            }
            else {
                acc.count <- acc.count + 1
            }
        }
        if (Msexsigma == 1) {
            bsigmac <- c(exp(rnorm(1, log(bsigma[1]), 2 * delta)), 
                bsigma[2])
            lp.sigmac <- bsigmac[Xsex[indid] + 1] * (c1 + c2)
            lpc <- loglam0 + Mb * beta * prevcap - lp.sigmac + 
                beta1 * Xeff + Msex * beta.sex * Xsex[indid]
            llvector.new <- lik.fn(lpc, y1)
            LM1[1:length(LM1)] <- llvector.new
            if (runif(1) < exp(sum(((LM1[z == 1, ] - LM2[z == 
                1, ]) \%*\% ones)))) {
                lp.sigma <- lp.sigmac
                bsigma <- bsigmac
                llvector <- llvector.new
                lp <- lpc
                LM2 <- LM1
            }
            else {
            }
            bsigmac <- c(bsigma[1], exp(rnorm(1, log(bsigma[2]), 
                2 * delta)))
            lp.sigmac <- bsigmac[Xsex[indid] + 1] * (c1 + c2)
            lpc <- loglam0 + Mb * beta * prevcap - lp.sigmac + 
                beta1 * Xeff + Msex * beta.sex * Xsex[indid]
            llvector.new <- lik.fn(lpc, y1)
            LM1[1:length(LM1)] <- llvector.new
            if (runif(1) < exp(sum(((LM1[z == 1, ] - LM2[z == 
                1, ]) \%*\% ones)))) {
                lp.sigma <- lp.sigmac
                bsigma[2] <- bsigmac[2]
                llvector <- llvector.new
                lp <- lpc
                LM2 <- LM1
            }
        }
        cat("DELTA: ", delta, fill = TRUE)
        if (any(bsigma < 0)) {
            cat("negative bsigma....", fill = TRUE)
            return(0)
        }
        if (Msex == 1) {
            beta.sexc <- rnorm(1, beta.sex, 0.1)
            lpc <- loglam0 + Mb * beta * prevcap - lp.sigma + 
                beta1 * Xeff + Msex * beta.sexc * Xsex[indid]
            llvector.new <- lik.fn(lpc, y1)
            LM1[1:length(LM1)] <- llvector.new
            if (runif(1) < exp(sum(((LM1[z == 1, ] - LM2[z == 
                1, ]) \%*\% ones)))) {
                beta.sex <- beta.sexc
                llvector <- llvector.new
                lp <- lpc
                LM2 <- LM1
            }
        }
        if (Mb == 1) {
            betac <- rnorm(1, beta, 0.1)
            lpc <- loglam0 + Mb * betac * prevcap - lp.sigma + 
                beta1 * Xeff + Msex * beta.sex * Xsex[indid]
            llvector.new <- lik.fn(lpc, y1)
            LM1[1:length(LM1)] <- llvector.new
            if (runif(1) < exp(sum(((LM1[z == 1, ] - LM2[z == 
                1, ]) \%*\% ones)))) {
                beta <- betac
                llvector <- llvector.new
                lp <- lpc
                LM2 <- LM1
            }
        }
        if (Xeff.tf) {
            beta1c <- rnorm(1, beta1, 0.1)
            lpc <- loglam0 + Mb * beta * prevcap - lp.sigma + 
                beta1c * Xeff + Msex * beta.sex * Xsex[indid]
            llvector.new <- lik.fn(lpc, y1)
            LM1[1:length(LM1)] <- llvector.new
            if (runif(1) < exp(sum(((LM1[z == 1, ] - LM2[z == 
                1, ]) \%*\% ones)))) {
                beta1 <- beta1c
                llvector <- llvector.new
                lp <- lpc
                LM2 <- LM1
            }
        }
        probz <- exp(rowsum(llvector[indid > nind], indid[indid > 
            nind]))
        probz <- (probz * psi)/(probz * psi + (1 - psi))
        z[(nind + 1):M] <- rbinom(M - nind, 1, probz)
        psi <- rbeta(1, 1 + sum(z), 1 + M - sum(z))
        if (!is.null(Xsex)) {
            tmp.sex <- Xsex
            tmp.sex[sex.naflag] <- 1 - Xsex[sex.naflag]
            if (Msexsigma == 0) 
                lp.sigmac <- Msigma * bsigma * (c1 + c2)
            if (Msexsigma == 1) 
                lp.sigmac <- bsigma[tmp.sex[indid] + 1] * (c1 + 
                  c2)
            lpc <- loglam0 + Mb * beta * prevcap - lp.sigmac + 
                beta1 * Xeff + Msex * beta.sex * tmp.sex[indid]
            llvector.new <- lik.fn(lpc, y1)
            lik.othersex <- exp(rowsum(llvector.new, indid))
            lik.sex <- exp(rowsum(llvector, indid))
            prior.curr <- (psi.sex^Xsex) * ((1 - psi.sex)^(1 - 
                Xsex))
            prior.cand <- (psi.sex^tmp.sex) * ((1 - psi.sex)^(1 - 
                tmp.sex))
            swtch <- sex.naflag & (runif(M, 0, 1) < ((lik.othersex * 
                prior.cand)/(lik.sex * prior.curr)))
            Xsex[swtch] <- 1 - Xsex[swtch]
            psi.sex <- rbeta(1, 0.1 + sum(Xsex[z == 1]), 0.1 + 
                sum(z) - sum(Xsex[z == 1]))
            if (Msexsigma == 0) 
                lp.sigma <- Msigma * bsigma * (c1 + c2)
            if (Msexsigma == 1) 
                lp.sigma <- bsigma[Xsex[indid] + 1] * (c1 + c2)
            lp <- loglam0 + Mb * beta * prevcap - lp.sigma + 
                beta1 * Xeff + Msex * beta.sex * Xsex[indid]
            llvector.new <- lik.fn(lp, y1)
            LM1[1:length(LM1)] <- llvector.new
            llvector <- llvector.new
            LM2 <- LM1
        }
        newcenters <- trunc(runif(M, 0, numnn[centers])) + 1
        newcenters <- NN[cbind(centers, newcenters)]
        qnew <- 1/numnn[centers]
        qold <- 1/numnn[newcenters]
        Sc <- G[newcenters, ]
        c1c <- (Sc[indid, 1] - trapgridbig[, 1])^2
        c2c <- (Sc[indid, 2] - trapgridbig[, 2])^2
        if (Msexsigma == 0) 
            lp.sigmac <- Msigma * bsigma * (c1c + c2c)
        if (Msexsigma == 1) 
            lp.sigmac <- bsigma[Xsex[indid] + 1] * (c1c + c2c)
        lpc <- loglam0 + Mb * beta * prevcap - lp.sigmac + beta1 * 
            Xeff + Msex * beta.sex * Xsex[indid]
        llvector.new <- lik.fn(lpc, y1)
        LM1[1:length(LM1)] <- llvector.new
        likdiff <- (LM1 - LM2) \%*\% ones
        likdiff[z == 0] <- 0
        likdiff <- likdiff + log(qold/qnew)
        accept <- runif(M) < exp(likdiff)
        cat("accept rate: ", mean(accept), fill = TRUE)
        S[accept, ] <- Sc[accept, ]
        centers[accept] <- newcenters[accept]
        c1 <- (S[indid, 1] - trapgridbig[, 1])^2
        c2 <- (S[indid, 2] - trapgridbig[, 2])^2
        LM2[accept, ] <- LM1[accept, ]
        if (Msexsigma == 0) 
            lp.sigma <- Msigma * bsigma * (c1 + c2)
        if (Msexsigma == 1) 
            lp.sigma <- bsigma[Xsex[indid] + 1] * (c1 + c2)
        lp <- loglam0 + Mb * beta * prevcap - lp.sigma + beta1 * 
            Xeff + Msex * beta.sex * Xsex[indid]
        llvector <- lik.fn(lp, y1)
        LM2 <- LM1
        if ((i > burn) & (i\%\%skip == 0)) {
            sigma <- sqrt(1/(2 * bsigma))
            if (Msexsigma == 0) {
                sigmatmp <- c(sigma, sigma)
                bsigmatmp <- c(bsigma, bsigma)
            }
            else {
                sigmatmp <- sigma
                bsigmatmp <- bsigma
            }
            logmu <- loglam0 + Mb * beta * prevcap - lp.sigma + 
                beta1 * Xeff + Msex * beta.sex * Xsex[indid]
            mu <- (1 - exp(-exp(logmu))) * z[indid]
            newy <- rbinom(length(mu), 1, mu)
            gof.stats <- cbind(y, newy, mu)
            gof.stats <- aggregate(gof.stats, list(indid), sum)
            gof.data[m] <- sum((sqrt(gof.stats[, 2]) - sqrt(gof.stats[, 
                4]))[z == 1]^2)
            gof.new[m] <- sum((sqrt(gof.stats[, 3]) - sqrt(gof.stats[, 
                4]))[z == 1]^2)
            zout[m, ] <- z
            Sout[m, ] <- centers
            out[m, ] <- c(bsigmatmp[1], sigmatmp[1], bsigmatmp[2], 
                sigmatmp[2], lam0, beta, beta1, beta.sex, psi, 
                psi.sex, sum(z), NA)
            print(out[m, ])
            if (m\%\%dumprate == 0) {
            }
            m <- m + 1
        }
    }
    list(out = out, G = G, Gunscaled = Gunscaled, traplocs = traplocs, 
        Sout = Sout, zout = zout, statespace = statespace, gof.data = gof.data, 
        gof.new = gof.new)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line




