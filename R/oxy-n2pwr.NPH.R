###input parameters
#' @title Power Calculation with Combination Test
#' @description  \code{n2pwr.NPH} calculates the power given either the
#' number of events or number of subjects using combination test
#' @param method a text specifying the calculation method, either
#' \code{"MaxLR"} or \code{"Projection"}. Maximum weighted
#' logrank test is used if \code{"MaxLR"} is specified; otherwise,
#' projection test is used.
#' @param entry a numeric value indicating the enrollment time, Default: 1
#' @param fup a numeric value indicating the minimum follow-up time for subjects.
#'  , Default: 1
#' @param maxfup maximum follow-up time
#' @param CtrlHaz a function,  specifying the hazard function for control group.
#' @param hazR a function, specifying the hazard ratio function between
#' treatment and control group
#' @param transP1 a numeric vector of length 2, consisting of the transition
#' probability from
#' receiving treatment to drop-out (drop-out rate) and
#' from receiving treatment to receiving control (drop-in rate) per time unit.
#' @param transP0 a numeric vector of length 2, consisting of the transition
#' probability from
#' receiving control to drop-out (drop-out rate) and
#' from receiving control to receiving treatment (drop-in rate) per time unit.
#' @param Wlist a list, consisting of weight functions applied to the test.
#' The element of the list must be functions. Default is a list of one constant
#' function, corresponding to the logrank test.
#' @param entry_pdf0 a function, indicating the probability density function (pdf)
#' of enrollment/entry time for control group. The default assumes a uniform distribution
#' corresponding to the constant enrollment rate.
#' Default: function(x) {
#'    (1/entry) * (x >= 0 & x <= entry)
#'}
#' @param entry_pdf1 a pdf function of enrollment/entry time for treatment
#' @param eventN the number of events
#' @param totalN the number of subjects
#' @param ratio allocation ratio, Default: 1
#' @param alpha type i error, Default: 0.05
#' @param alternative alternative hypothesis - one of c(\code{"two.sided", "less", "greater"}),
#' Default: \code{"two.sided"}
#' @param k an integer, indicating number of sub-intervals per time unit, Default: 100
#' @param nreps number of replicates used for calculating quantitle
#' using multivariate normal
#' @return
#' a list of components:
#' \item{power}{asymptotic power }
#' \item{inN}{a vector consisting of the input of \code{eventN} and \code{totalN}}
#' \item{outN}{a vector including the output of number of events
#' and total sample. See details. }
#' \item{prob_event}{event probability at the end of trial}
#' \item{L_trans}{a list, consisting of transition matrix at each interval}
#' \item{pdat}{ a data frame including all the intermediate variables in the calculation.
#' }
#'  \item{studytime}{a vector of length 2, including the entry and follow-up time as input}
#'  \item{RandomizationRatio}{as input}
#' @details
#' Function \code{npwr.NPH} calculates the asymptotic power given number
#' of events or number of subjects using maximum weighted logrank test or
#' projection type test. If only \code{eventN} is provided, the
#' asymptotic power is based on provided number of events. If
#' only \code{totalN} is given, the pooled event probability (\eqn{eprob}) is calculated
#' according input design parameters including entry time, follow-up time
#' and hazard functions, etc. The number of events is calculated as
#' \code{totalN}*\eqn{eprob}, which is given in returned vector \code{outN}.
#' Similarly, if only \code{eventN} is given, the total sample
#' size is given as \code{eventN}/\eqn{eprob}. However, if both
#' \code{eventN} and \code{totalN} are provided, we only use
#' \code{eventN} for calculation.
#' Check function \code{pwr2n.NPH} for more calculation details.
#' @examples
#' # entry time
#' t_enrl <- 12
#' # follow-up time
#' t_fup <- 18
#' # baseline hazard
#' lmd0 <- -log(0.2)/10
#' # delayed treatment effects
#' f_hr_delay <- function(x){(x<=6)+(x>6)*0.75}
#' # maxcombo test
#' maxc <- gen.wgt(method="Maxcombo")
#' pwr1 <- n2pwr.NPH(entry   = t_enrl
#'                  ,fup      = t_fup
#'                  ,CtrlHaz = function(x){x^0*lmd0}
#'                  ,hazR = f_hr_delay
#'                  ,transP1 = c(0,0)
#'                  ,transP0 = c(0,0)
#'                  ,Wlist = maxc
#'                  ,eventN = 50 # targeted number of events
#')
#' @seealso
#'  \code{\link{pwr2n.NPH}}
#' @rdname n2pwr.NPH
#' @export
#' @importFrom mvtnorm qmvnorm pmvnorm
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

  #old.total <- totalN
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

