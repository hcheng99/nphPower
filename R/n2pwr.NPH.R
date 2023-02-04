

###input parameters
n2pwr.NPH<- function(method = "MaxLR"
                     ,entry   = 1
                     ,fup      = 1
                     ,maxfup = entry+fup
                     ,CtrlHaz
                     ,hazR
                     ,transP1
                     ,transP0
                     ,Wlist
                     ,entry_pdf0=function(x){(1/entry)*(x>=0&x<=entry)}
                     ,entry_pdf1=entry_pdf0
                     ,eventN = NULL
                     ,totalN = NULL
                     ,ratio    = 1
                     ,alpha    = 0.05
                     ,alternative=c("two.sided")
                     ,k        = 100
                     ,nreps = 10
){
  if (is.null(eventN)) {
    mce <- 1
    old.event <- NA
  }
  else {
    mce <- 0
    old.event <- eventN
  }
  if (is.null(totalN)) {
    mct <- 1
    old.total <- NA
  }else {
    mct <- 0
    old.total <- totalN
  }

  if (mce+mct==2){
    stop("At least one of eventN/totalN must be provided")
  }else if (mce+mct==0){
    message("Both number of events and number of subjects are specified.
             The asymptotic power is calculated based on number of events.")
  }


  tot_time <- maxfup
  num <- k*tot_time
  # create the subintervals
  x <- seq(1/k,tot_time,by=1/k)
  ctrlRate <- CtrlHaz(x)
  haz_val <- hazR(x)*ctrlRate
  haz_point <- x*k
  ## load the transition matrix
  load <- trans.mat(numN=num,x=x,ctrlRate=ctrlRate,haz_val=haz_val,
                    haz_point=haz_point,ratio=ratio,
                    transP1=transP1,transP0=transP0,k=k,
                    fup=fup,entry=entry,entry_pdf0=entry_pdf0,
                    entry_pdf1=entry_pdf1,hazR=hazR,tot_time=tot_time)

  pdat <- load$pdat
  eprob <- stats::weighted.mean(c(pdat$C_E[num],pdat$E_E[num]),w=c(1,ratio))
  if (mce==1){ eventN <- round(totalN*eprob)}
  if (mct==1|mce+mct==0){ totalN <- round(eventN/eprob)}
  wn <- length(Wlist)
  W <- matrix(NA,nrow=nrow(pdat),ncol=wn)
  ## calculate the variance-covariance matrix
  Vmat <- matrix(NA,nrow=wn,ncol=wn)

  for (j in 1:wn){
    W[,j] <- Wlist[[j]](pdat$S)
  }
  dnum <- eventN
  for (k1 in 1:wn){
    for (k2 in 1:wn){
      Vmat[k1,k2] <-  dnum*t(W[,k1]*W[,k2]) %*%(pdat$rho*pdat$eta)
    }
  }
  if (method == "MaxLR"){
    rho_est <- stats::cov2cor(Vmat)
    mu <- as.vector(dnum*t(W)%*%(pdat$rho*pdat$gamma))/sqrt(diag(Vmat))
    if (alternative=="two.sided"){
      ftwo <- function(i){
        crit <- mvtnorm::qmvnorm(1-alpha,tail="both.tails",
                                 mean=rep(0,wn),sigma = rho_est)$quantile
        power <- 1-mvtnorm::pmvnorm(-crit,crit,mean=mu,sigma = rho_est)
        return(c(crit,power))
      }
      ftwod <- apply(do.call(rbind,sapply(1:nreps,ftwo,simplify=FALSE)) ,2,mean)
      power <- ftwod[2]
    }else if (alternative=="less"){ # l1 <l0
      ftwo <- function(i){
        crit <- mvtnorm::qmvnorm(alpha,tail="lower.tail",mean=rep(0,wn),sigma = rho_est)$quantile
        power <- 1-mvtnorm::pmvnorm(crit,Inf,mean=mu,sigma = rho_est)[1]
        return(c(crit,power))
      }
      ftwod <- apply(do.call(rbind,sapply(1:nreps,ftwo,simplify=FALSE)) ,2,mean)
      crit <- ftwod[1];  power <- ftwod[2]
    }else if (alternative=="greater"){
      ftwo <- function(i){
        crit <- qmvnorm(alpha,tail="upper.tail",mean=rep(0,wn),sigma = rho_est)$quantile
        power <- 1-pmvnorm(-Inf,crit,mean=mu,sigma = rho_est)[1]
        return(c(crit,power))
      }
      ftwod <- apply(do.call(rbind,sapply(1:nreps,ftwo,simplify=FALSE)) ,2,mean)
      crit <- ftwod[1];  power <- ftwod[2]
    }
  }
  else if (method == "Projection"){
    if (alternative!="two.sided"){
      message(c("note: only two-sided is supported for projection test."))
    }
    ## get the rank of the variance matrix
    mu <- as.vector(dnum*t(W)%*%(pdat$rho*pdat$gamma))
    vr <- qr(Vmat)$rank
    crit <- stats::qchisq(1-alpha,df=vr)
    ## get the noncentral parameter
    lmd <- t(mu)%*%MASS::ginv(Vmat)%*%mu
    power <- stats::pchisq(crit,df=vr,ncp=lmd,lower.tail = FALSE)

  }
  inn <- c(old.event, old.total)
  names(inn) <- c("event","total")
  outn <- c(eventN, totalN)
  names(outn) <- c("event","total")

  listall <- list( power = as.numeric(power)
                   ,inN  = inn
                   ,outN = outn
                   ,prob_event =eprob
                   ,L_trans = load$L_trans
                   ,pdat = pdat
                   ,studytime=c(entry,fup)
                   ,RandomizationRatio=ratio

  )

  return(listall)


}
