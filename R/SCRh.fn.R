SCRh.fn <-
function (scrobj, ni = 1100, burn = 100, skip = 2, nz = 200,
Msigma = 1,
Mb = 0,
Msex = 0,
Msexsigma = 0,
Mss=0,
Meff = 0,
Mtel = 0,
theta = 1,
coord.scale = 1000, area.per.pixel = 1,
    thinstatespace = 1, maxNN = 20, dumprate = 1000)
{
    call <- match.call()
    traps <- scrobj$traps
    captures <- scrobj$captures
    statespace <- scrobj$statespace
    alive <- scrobj$alive
# Mss = fit "Xd" covariate
    Xd <- scrobj$Xd
# Meff = fit "Xeff" covariate
    Xeff <- scrobj$Xeff
    Xsex <- scrobj$Xsex
    Xtel <- scrobj$Xtel
    if(is.null(scrobj$Ytel)){
        loglik.tel<-function(bsigma,Ytel,Msexsigma){
            return(0)
            }
        }
    if(!is.null(scrobj$Ytel)){
       Ytel<- scrobj$Ytel
    if(is.null(Xtel)){
        loglik.tel<-function(bsigma,Ytel,Msexsigma){
            sbar<-attributes(Ytel)$sbar
            sex.tel<-attributes(Ytel)$sex

            sigma <- sqrt(1/(2*bsigma) )
            if(length(bsigma)==1){
return(            sum(    dnorm(Ytel[,1],sbar[,1][Ytel[,3]],sigma,log=TRUE) +
           dnorm(Ytel[,2],sbar[,2][Ytel[,3]],sigma,log=TRUE)            ))
            }
            if(length(bsigma==2)){
              tmp<- sex.tel[Ytel[,3]]
             return(sum(dnorm(Ytel[,1],sbar[,1][Ytel[,3]],sigma[tmp+1],log=TRUE) +
           dnorm(Ytel[,2],sbar[,2][Ytel[,3]],sigma[tmp+1],log=TRUE)            ))
          }
        }
        ##tel.gr<-    make.statespace(ll = statespace, buffer = 0.01, nx = 40)



    }

if(!is.null(Xtel)){
# extract coordinates and covariate values here

}


    }


    captures <- cbind(captures[, 4], captures[, 2], captures[,
        3])
    if (length(unique(captures[, 2])) != length(min(captures[,
        2]):max(captures[, 2]))) {
        cat("Error: individuals not numbered sequentially, renumbering them now",
            fill = TRUE)
        captures[, 2] <- as.numeric(factor(captures[, 2]))
    }
    Y <- captures
    traplocs <- traps[, 2:3]

    xrange <- max(traplocs[,1])-min(traplocs[,1])
    yrange <- max(traplocs[,2])-min(traplocs[,2])
if( (xrange + yrange) < coord.scale)
    {
cat("coord.scale appears to be too large, check this out, maybe set to 1",fill=TRUE)
         return(0)
    }
    MASK <- as.matrix(traps[, 4:ncol(traps)])
    nind <- max(Y[, 2])
    T <- dim(MASK)[2]
    M <- nind + nz
    ntraps <- nrow(traplocs)
    totalarea <- nrow(statespace) * area.per.pixel
    thinned <- seq(1, nrow(statespace), thinstatespace)
    statespace <- statespace[thinned, ]
    Xd <- Xd[thinned]
    goodbad <- statespace[, 3]
    G <- statespace[, 1:2]
    G <- G[goodbad == 1, ]
    Gunscaled <- G
    Xd <- Xd[goodbad == 1]
    nG <- nrow(G)
    new.area.per.pixel <- totalarea/nG
    mgx <- min(traplocs[, 1])
    mgy <- min(traplocs[, 2])
    traplocs[, 1] <- (traplocs[, 1] - mgx)/coord.scale
    traplocs[, 2] <- (traplocs[, 2] - mgy)/coord.scale
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
    if (Meff ==1 ) {
        #####################!is.null(Xeff)) {
        Xeffnew <- array(0, dim = c(nind + nz, T, ntraps))
        for (j in 1:M) {
            Xeffnew[j, 1:T, 1:ntraps] <- t(Xeff)
        }
       #####################  Xeff.tf <- TRUE
    }
    if (Meff ==0){
        Xeffnew <- array(0, dim = c(nind + nz, T, ntraps))
       ####################### Xeff.tf <- FALSE
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
    alive.trues <- array(1, c(nind + nz, T, ntraps))
    for (i in 1:nind) {
        for (t in 1:T) {
            alive.trues[i, T, 1:ntraps] <- alive[i, t]
        }
    }
    alive.trues <- as.vector(alive.trues)
    aliveid <- alive.trues[msk2 == 1]
    arr.trues <- array(TRUE, c(nind + nz, T, ntraps))
    idx <- which(arr.trues, arr.ind = TRUE)
    y <- as.vector(Yaug)
    y <- y[msk2 == 1 & alive.trues == 1]
    Xeff <- as.vector(Xeff)
    Xeff <- Xeff[msk2 == 1 & alive.trues == 1]
    prevcap <- prevcap[msk2 == 1 & alive.trues == 1]
    indid.LM <- idx[msk2 == 1, 1]
    indid <- idx[msk2 == 1 & alive.trues == 1, 1]
    repid <- idx[msk2 == 1 & alive.trues == 1, 2]
    trapid <- idx[msk2 == 1 & alive.trues == 1, 3]
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
        tt <- row(tt)[tt == 1]
        xxx <- as.matrix(traplocs[tt, ], ncol = 2, byrow = FALSE)
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
    update.theta <- FALSE
    if (is.na(theta)) {
        theta <- 0.75
        update.theta <- TRUE
    }
    loglam0 <- log(0.018)
    beta.behave <- 0
    beta.sex <- 0
    beta1 <- 0
    lam0 <- exp(loglam0)
    psi <- 0.5
    psi.sex <- mean(Xsex, na.rm = TRUE)
    z <- c(rep(1, nind), rbinom(nz, 1, psi))
    if (sum(sex.naflag) > 0)
        Xsex[sex.naflag] <- rbinom(sum(sex.naflag), 1, 0.65)
 # obsolete
 #   if (is.null(Xd)) {
 #       Xd <- rep(1, nG)
 #   }
    beta.den <- 0


    lik.fn <- function(lpc, y1) {
        llvector.new <- -1 * exp(lpc)
        part2 <- exp(exp(lpc[y1])) - 1
        part2[part2 == 0] <- .Machine$double.eps
        llvector.new[y1] <- llvector.new[y1] + log(part2)
        llvector.new
    }
    trapgridbig <- traplocs[trapid, ]
    y1 <- y == 1
    c1 <- (S[indid, 1] - trapgridbig[, 1])^2
    c2 <- (S[indid, 2] - trapgridbig[, 2])^2
    gof.new <- gof.data <- rep(NA, (ni - burn)/skip)
    out <- matrix(NA, nrow = (ni - burn)/skip, ncol = 16)
    dimnames(out) <- list(NULL, c("bsigma", "sigma", "bsigma2",
        "sigma2", "loglam0","lam0", "beta.behave", "beta1(effort)", "beta.sex",
        "psi", "psi.sex", "Nsuper", "theta", "beta.density",
        "D","beta(habitat)"))
    zout <- matrix(NA, nrow = (ni - burn)/skip, ncol = M)
    Sout <- matrix(NA, nrow = (ni - burn)/skip, ncol = M)
    m <- 1
    LM1 <- LM2 <- matrix(0, nrow = M, ncol = length(indid.LM)/M)
    ones <- rep(1, ncol(LM1))
    if (Msexsigma == 0)
        lp.sigma <- Msigma * bsigma * (c1 + c2)^theta
    if (Msexsigma == 1)
        lp.sigma <- bsigma[Xsex[indid] + 1] * (c1 + c2)^theta
    acc.count <- 0
    delta <- 0.05


    for (i in 1:ni) {
        cat("iter: ", i, fill = TRUE)
        lp <- loglam0 + Mb * beta.behave * prevcap - lp.sigma +
            Meff* beta1 * Xeff + Msex * beta.sex * Xsex[indid]
        loglam0c <- rnorm(1, loglam0, 0.1)
        lpc <- loglam0c + Mb * beta.behave * prevcap - lp.sigma +
            Meff*beta1 * Xeff + Msex * beta.sex * Xsex[indid]
        llvector <- lik.fn(lp, y1)
        llvector.new <- lik.fn(lpc, y1)
        LM1[aliveid == 1] <- llvector.new
        LM2[aliveid == 1] <- llvector

        if (runif(1) < exp(sum(((LM1[z == 1, ] - LM2[z == 1,
            ]) %*% ones)))) {
            loglam0 <- loglam0c
            lam0 <- exp(loglam0)
            llvector <- llvector.new
            lp <- lpc
            LM2 <- LM1
        }

        if (update.theta) {
            lp <- loglam0 + Mb * beta.behave * prevcap - lp.sigma +
                Meff* beta1 * Xeff + Msex * beta.sex * Xsex[indid]
            thetac <- rnorm(1, theta, 0.02)
            if (thetac >= 0.5 & thetac <= 1) {
                if (Msexsigma == 0)
                  lp.sigmac <- Msigma * bsigma * (c1 + c2)^thetac
                if (Msexsigma == 1)
                  lp.sigmac <- bsigma[Xsex[indid] + 1] * (c1 +
                    c2)^thetac
                lpc <- loglam0 + Mb * beta.behave * prevcap -
                  lp.sigmac + Meff*beta1 * Xeff + Msex * beta.sex *
                  Xsex[indid]
                llvector <- lik.fn(lp, y1)
                llvector.new <- lik.fn(lpc, y1)
                LM1[aliveid == 1] <- llvector.new
                LM2[aliveid == 1] <- llvector
                if (runif(1) < exp(sum(((LM1[z == 1, ] - LM2[z ==
                  1, ]) %*% ones)))) {
                  theta <- thetac
                  lp.sigma <- lp.sigmac
                  llvector <- llvector.new
                  lp <- lpc
                  LM2 <- LM1
                }
            }
        }
        if (Msexsigma == 0) {
            bsigmac <- exp(rnorm(1, log(bsigma), delta))
            lp.sigmac <- Msigma * bsigmac * (c1 + c2)^theta
            lpc <- loglam0 + Mb * beta.behave * prevcap - lp.sigmac +
                Meff*beta1 * Xeff + Msex * beta.sex * Xsex[indid]
            llvector.new <- lik.fn(lpc, y1)
            LM1[aliveid == 1] <- llvector.new

            mh.ratio<-
            exp(sum(((LM1[z == 1, ] - LM2[z ==1, ]) %*% ones)) +
                Mtel*(loglik.tel(bsigmac,Ytel,Msexsigma) - loglik.tel(bsigma,Ytel,Msexsigma))  )
            if (runif(1) < mh.ratio){
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
            lp.sigmac <- bsigmac[Xsex[indid] + 1] * (c1 + c2)^theta
            lpc <- loglam0 + Mb * beta.behave * prevcap - lp.sigmac +
                Meff*beta1 * Xeff + Msex * beta.sex * Xsex[indid]
            llvector.new <- lik.fn(lpc, y1)
            LM1[aliveid == 1] <- llvector.new
            mh.ratio<-             exp(sum(((LM1[z == 1, ] - LM2[z ==
                                             1, ]) %*% ones)) +
               Mtel*(loglik.tel(bsigmac,Ytel,Msexsigma) - loglik.tel(bsigma,Ytel,Msexsigma))  )

            if (runif(1) < mh.ratio){
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
            lp.sigmac <- bsigmac[Xsex[indid] + 1] * (c1 + c2)^theta
            lpc <- loglam0 + Mb * beta.behave * prevcap - lp.sigmac +
                Meff*beta1 * Xeff + Msex * beta.sex * Xsex[indid]
            llvector.new <- lik.fn(lpc, y1)
            LM1[aliveid == 1] <- llvector.new
 mh.ratio<-  exp(sum(((LM1[z == 1, ] - LM2[z ==1, ]) %*% ones)) +
             Mtel*(loglik.tel(bsigmac,Ytel,Msexsigma) - loglik.tel(bsigma,Ytel,Msexsigma))  )
            if (runif(1) < mh.ratio) {
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
            lpc <- loglam0 + Mb * beta.behave * prevcap - lp.sigma +
                Meff*beta1 * Xeff + Msex * beta.sexc * Xsex[indid]
            llvector.new <- lik.fn(lpc, y1)
            LM1[aliveid == 1] <- llvector.new
            if (runif(1) < exp(sum(((LM1[z == 1, ] - LM2[z ==
                1, ]) %*% ones)))) {
                beta.sex <- beta.sexc
                llvector <- llvector.new
                lp <- lpc
                LM2 <- LM1
            }
        }
        if (Mb == 1) {
            beta.behave.c <- rnorm(1, beta.behave, 0.1)
            lpc <- loglam0 + Mb * beta.behave.c * prevcap - lp.sigma +
                Meff*beta1 * Xeff + Msex * beta.sex * Xsex[indid]
            llvector.new <- lik.fn(lpc, y1)
            LM1[aliveid == 1] <- llvector.new
            if (runif(1) < exp(sum(((LM1[z == 1, ] - LM2[z ==
                1, ]) %*% ones)))) {
                beta.behave <- beta.behave.c
                llvector <- llvector.new
                lp <- lpc
                LM2 <- LM1
            }
        }
        if (Meff == 1){
            beta1c <- rnorm(1, beta1, 0.1)
            lpc <- loglam0 + Mb * beta.behave * prevcap - lp.sigma +
                Meff*beta1c * Xeff + Msex * beta.sex * Xsex[indid]
            llvector.new <- lik.fn(lpc, y1)
            LM1[aliveid == 1] <- llvector.new
            if (runif(1) < exp(sum(((LM1[z == 1, ] - LM2[z ==
                1, ]) %*% ones)))) {
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
                lp.sigmac <- Msigma * bsigma * (c1 + c2)^theta
            if (Msexsigma == 1)
                lp.sigmac <- bsigma[tmp.sex[indid] + 1] * (c1 +
                  c2)^theta
            lpc <- loglam0 + Mb * beta.behave * prevcap - lp.sigmac +
                Meff*beta1 * Xeff + Msex * beta.sex * tmp.sex[indid]
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
                lp.sigma <- Msigma * bsigma * (c1 + c2)^theta
            if (Msexsigma == 1)
                lp.sigma <- bsigma[Xsex[indid] + 1] * (c1 + c2)^theta
            lp <- loglam0 + Mb * beta.behave * prevcap - lp.sigma +
                Meff*beta1 * Xeff + Msex * beta.sex * Xsex[indid]
            llvector.new <- lik.fn(lp, y1)
            LM1[aliveid == 1] <- llvector.new
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
            lp.sigmac <- Msigma * bsigma * (c1c + c2c)^theta
        if (Msexsigma == 1)
            lp.sigmac <- bsigma[Xsex[indid] + 1] * (c1c + c2c)^theta
        lpc <- loglam0 + Mb * beta.behave * prevcap - lp.sigmac +
            Meff*beta1 * Xeff + Msex * beta.sex * Xsex[indid]
        llvector.new <- lik.fn(lpc, y1)
        LM1[aliveid == 1] <- llvector.new
        likdiff <- (LM1 - LM2) %*% ones
        likdiff[z == 0] <- 0
        logprior.new <- Xd[newcenters] * beta.den
        logprior.old <- Xd[centers] * beta.den
        likdiff <- likdiff + log(qold/qnew) + (logprior.new -
            logprior.old)
        accept <- runif(M) < exp(likdiff)
        cat("accept rate: ", mean(accept), fill = TRUE)
        S[accept, ] <- Sc[accept, ]
        centers[accept] <- newcenters[accept]
        c1 <- (S[indid, 1] - trapgridbig[, 1])^2
        c2 <- (S[indid, 2] - trapgridbig[, 2])^2
        LM2[accept, ] <- LM1[accept, ]
   if(Mss==1){
        beta.den.c <- rnorm(1, beta.den, 0.25)
        numerator.new <- exp(Xd * beta.den.c)
        loglik.new <- Xd[centers] * beta.den.c - log(sum(numerator.new))
        numerator.old <- exp(Xd * beta.den)
        loglik.old <- Xd[centers] * beta.den - log(sum(numerator.old))
        if (runif(1) < exp(sum(loglik.new) - sum(loglik.old))) {
            beta.den <- beta.den.c
        }
        }
        if (Msexsigma == 0)
            lp.sigma <- Msigma * bsigma * (c1 + c2)^theta
        if (Msexsigma == 1)
            lp.sigma <- bsigma[Xsex[indid] + 1] * (c1 + c2)^theta
        lp <- loglam0 + Mb * beta.behave * prevcap - lp.sigma +
            Meff*beta1 * Xeff + Msex * beta.sex * Xsex[indid]
        llvector <- lik.fn(lp, y1)
        LM2 <- LM1


        ###
        ###
        ###
        ###
        if ((i > burn) & (i%%skip == 0)) {
            sigma <- sqrt(1/(2 * bsigma))
            if (Msexsigma == 0) {
                sigmatmp <- c(sigma, sigma)
                bsigmatmp <- c(bsigma, bsigma)
            }
            else {
                sigmatmp <- sigma
                bsigmatmp <- bsigma
            }
            logmu <- loglam0 + Mb * beta.behave * prevcap - lp.sigma +
                Meff*beta1 * Xeff + Msex * beta.sex * Xsex[indid]
            mu <- (1 - exp(-exp(logmu))) * z[indid]
            newy <- rbinom(length(mu), 1, mu)
            gof.stats <- cbind(y, newy, mu)
            gof.stats <- aggregate(gof.stats, list(indid), sum)
            gof.data[m] <- sum((sqrt(gof.stats[, 2]) - sqrt(gof.stats[,
                4]))[z == 1]^2)
            gof.new[m] <- sum((sqrt(gof.stats[, 3]) - sqrt(gof.stats[,
                4]))[z == 1]^2)
            density <- sum(z)/totalarea
            zout[m, ] <- z
            Sout[m, ] <- centers
            out[m, ] <- c(bsigmatmp[1], sigmatmp[1], bsigmatmp[2],
                sigmatmp[2], log(lam0), lam0, beta.behave, beta1, beta.sex,
                psi, psi.sex, sum(z), theta, beta.den, density,NA)
            print(out[m, ])
            if (m%%dumprate == 0) {
            }
            m <- m + 1
        }
    }
    parms.2.report <- c(TRUE, TRUE,
                        Msexsigma == 1,
                        Msexsigma == 1, TRUE, TRUE, Mb == 1,
                        Meff ==1, Msex == 1, TRUE, Msex == 1,
        TRUE, update.theta, Mss == 1, TRUE,FALSE)
    out <- list(mcmchist = out, G = G, Gunscaled = Gunscaled,
        traplocs = traplocs, Sout = Sout, zout = zout, statespace = statespace,
        gof.data = gof.data, gof.new = gof.new, call = call,
        parms2report = parms.2.report)
    class(out) <- c("scrfit", "list")
    return(out)
}
