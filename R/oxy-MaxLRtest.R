#' @title Maximum Weighted Logrank Test
#' @description \code{MaxLRtest} performs the maximum weighted logrank test if
#' multiple weight functions are provided. It is the regular weighted logrank test,
#' if a single weight function is specified,
#'
#' @param dat a dataframe or matrix. The first three columns of the data
#' set are survival
#' time, event status indicator  and group label. The status indicator, normally
#' 0=alive, 1=dead/event. Other choices are TRUE/FALSE (TRUE=death) or 1/2 (2=death).
#' The group label can be either numeric values like 0=control, 1=treatment or text
#' like C=control, T=treatment.
#'
#' @param Wlist a list with components of weight functions
#' @param base a text must be one of c("\code{KM}", "\code{Combined}", "\code{N}"),
#'  Default: c("KM")
#' @param alpha a number indicating type I error rate, Default: 0.05
#' @param alternative a text must be one of c("\code{two.sided}", "\code{less}",
#'  \code{"greater"}), indicating the alternative hypothesis, Default: c("two.sided")
#' @return
#' a list of components including
#' \item{stat}{a numeric value indicating the test statistic. It is logrank or weighted logrank test
#' statistic if one weight function is specified. Otherwise, it gives
#' the maximum weighted logrank test statistic, which takes the maximum
#' of absolute values of all the statistics. }
#' \item{stat.mat}{a matrix with the first column showing weighted
#' logrank test statistics and other columns displaying the variance and
#' covariance between statistics}
#' \item{critV}{a numeric value indicating the critical value corresponding to the nominal
#' level - \code{alpha}}
#' \item{details}{a dataframe showing the intermediate variables used in the
#' calculation. }
#' \item{p.value}{a numeric value indicating the p-value of the test}
#' @details
#' \code{MaxLRtest} function performs logrank, weighted logrank test such as
#' Fleming-Harrington test and maximum weighted logrank test depending on
#' the type and number of weight functions. Let \eqn{w(x_t)} denote the weight applied
#' at event time point \eqn{t}, where \eqn{x_t} is the base function. There are three options
#' for \code{base}. If \code{KM} is used, \eqn{x_t=1-S_t}, where \eqn{S_t}
#' is pooled Kaplan-Meier estimate of survival rate at time point t. A FH(1,0) test
#' needs a weight function \eqn{1-x_t}. If \code{Combined} base is selected,
#' \eqn{x_t=1-S^*_t}, where \eqn{S^*_t=w_1S^1_t+w_0S^0_t}, the weighted average
#' of KM estimate of survival rate for treatment (\eqn{S^1_t}) and control group
#' (\eqn{S^0_t}). It is considered more robust in case of unbalanced data.
#' For option \code{N}, \eqn{x_t=1-\frac{Y_t}{N}}, where \eqn{Y_t} is the subjects
#' at risk at time t and \eqn{N} is the total number of subjects.The Wilcoxon and
#' tarone test should use this base. The base \eqn{x_t} in all three cases is an
#' increasing function of time t. Function \code{gen.wgt} helps to generate the commonly
#' used weight functions.
#'
#'  Let \eqn{\Lambda_1} and \eqn{\Lambda_0} denote the cumulative hazard for
#'  treatment and control group. The alternative of a two-sided test is
#'  \eqn{H_a: \Lambda_1 \neq \Lambda_0}. The \code{"less"} alternative
#'  corresponds to \eqn{H_a: \Lambda_1 < \Lambda_0} and the \code{"greater"}
#'  alternative is \eqn{H_a: \Lambda_1 > \Lambda_0}.
#'
#' A p-value is obtained from a multivariate normal distribution if multiple weights
#' are provided. The function \code{pmvnorm} from R package \pkg{mvtnorm} is used.
#' Because the algorithm is slightly seed-dependent,the p-value and critical value
#' is the average of 10 runs.
#' @examples
#' data(lung)
#' #Only keep variables for analysis
#' tmpd <- with(lung, data.frame(time=SurvTime,stat=1-censor,grp=Treatment))
#' #logrank test
#' wlr <- gen.wgt(method = "LR")
#' t1 <- MaxLRtest(tmpd, Wlist = wlr, base = c("KM") )
#' t1$stat ;t1$p.value
#'
#'
#' # maxcombo test
#' wmax <- gen.wgt(method="Maxcombo")
#' t2 <- MaxLRtest(tmpd, Wlist = wmax, base = c("KM") )
#' t2$stat ;t2$p.value
#' #visualize the weight functions
#' plot(t2)
#' @seealso
#'  \code{\link{pwr2n.NPH}}, \code{\link{gen.wgt}}
#' @export
#' @importFrom mvtnorm qmvnorm pmvnorm
#' @importFrom survival Surv survfit
MaxLRtest <- function(dat
                      ,Wlist
                      ,base=c("KM")
                      ,alpha=0.05
                      ,alternative=c("two.sided")

){

  if (!missing(base)&!base %in% c("KM","Combined","N")){
    stop("base must be one of 'KM','Combined','N'")
  }
  ## transform the input data to the matrix ready for analysis
  datM <- load.mat(dat)
  ## only keep the timepoints with event
  datM <- datM[datM$event>0,]
  if (base=="Combined"){

    tnew0 <- c(1,
               cumprod(chkV(datM$n.risk.x,1-datM$n.event.x/datM$n.risk.x))*table(dat[,3])[2]/nrow(dat)+
                 cumprod(chkV(datM$n.risk.y,1-datM$n.event.y/datM$n.risk.y))*table(dat[,3])[1]/nrow(dat))
    tnew <- 1-tnew0[1:nrow(datM)]

  }
  if (base=="KM"|missing(base)){
    # based on the pooled survival
    s_fit<- survival::survfit(survival::Surv(dat[,1], dat[,2])~1 , data = dat)
    f_s <- stats::stepfun(s_fit$time,y=c(0,1-s_fit$surv),right = TRUE)
    tnew <- f_s(datM$time)
  }
  if (base=="N"){
    # based on time/person at risk
    # tnew <- datM$time/max(datM$time)
    tnew <- 1-datM$risk/nrow(dat)
  }
  ## number of weight functions
  wn <- length(Wlist)
  W <- matrix(NA,nrow=nrow(datM),ncol=wn)
  ## calculate the variance-covariance matrix
  Vmat <- matrix(NA,nrow=wn,ncol=wn)
  for (j in 1:wn){
    W[,j] <- Wlist[[j]](tnew)
  }
  for (k1 in 1:wn){
    for (k2 in 1:wn){
      Vmat[k1,k2] <-  t(W[,k1]*W[,k2]) %*%datM$V
    }
  }
  Zstat <- apply(W,2,function(x) {CalZ(x,data=datM)})
  ## correlation matrix
  rho_est <- stats::cov2cor(Vmat)
  if (alternative=="two.sided"){
    statistic <- max(abs(Zstat))
    ftwo <- function(i){
      crit <- mvtnorm::qmvnorm(1-alpha,tail="both.tails",mean=rep(0,wn),
                               sigma = rho_est)$quantile
      p.value <- 1-mvtnorm::pmvnorm(-1*statistic,statistic,mean=rep(0,wn),
                                    sigma = rho_est)[1]
      return(c(crit,p.value))
    }
    ftwod <- apply(do.call(rbind,sapply(1:10,ftwo,simplify=FALSE)) ,2,mean)
    crit <- as.numeric(ftwod[1]);  p.value <- as.numeric(ftwod[2])
  }else if (alternative=="greater"){
    statistic <- max(abs(Zstat))*sign(Zstat[1])
    ftwo <- function(i){

      crit <- qmvnorm(1-alpha,tail="lower.tail",mean=rep(0,wn),sigma = rho_est)$quantile
      p.value <- 1-pmvnorm(-Inf,statistic,mean=rep(0,wn),sigma = rho_est)[1]
      return(c(crit,p.value))
    }
    ftwod <- apply(do.call(rbind,sapply(1:10,ftwo,simplify=FALSE)) ,2,mean)
    crit <- as.numeric(ftwod[1]);  p.value <- as.numeric(ftwod[2])

  }else if (alternative=="less"){
    statistic <- max(abs(Zstat))*sign(Zstat[1])
    ftwo <- function(i){

      crit <- qmvnorm(1-alpha,tail="upper.tail",mean=rep(0,wn),sigma = rho_est)$quantile
      p.value <- 1-pmvnorm(statistic,Inf,mean=rep(0,wn),sigma = rho_est)[1]
      return(c(crit,p.value))
    }

    ftwod <- apply(do.call(rbind,sapply(1:10,ftwo,simplify=FALSE)) ,2,mean)
    crit <- as.numeric(ftwod[1]);  p.value <- as.numeric(ftwod[2])
  }


  res <- statistic>=crit
  stat.mat <- data.frame(Zstat,rho_est)
  names(stat.mat) <- c("Zstat",paste0("W",1:wn))
  Rdat <- data.frame(datM[,c("time","event","risk","n.event.x","exp.x","n.risk.x","V")],W)
  names(Rdat)[-c(1:7)] <- paste0("W",1:wn)
  list <-list(stat=statistic
              ,stat.mat=stat.mat
              ,critV=crit
              ,details=Rdat
              ,p.value=p.value)
  class(list) <- 'MaxLR'

  return(list)

}

#' @title Graphical Display of Weight Functions
#' @description Display weight functions used in the function \code{MaxLRtest}
#' @param x object of \code{MaxLRtest} function
#' @param ... additional graphical arguments passed to the plot function
#' @return
#' Plots are produced on the current graphics device
#' @examples
#' # See examples in the help file of function MaxLRtest
#' @rdname plot.MaxLR
#' @seealso
#'  \code{\link{MaxLRtest}}
#' @export
#' @importFrom graphics points legend mtext
plot.MaxLR <- function(x,...) {
  datM <- x$details
  xtime <- datM$time
  wn <- ncol(datM)-7
  plot(xtime,datM$n.event.x,cex=1,col=1,xlab="time",
       ylab="weight",pch=1,ylim=c(-1,1),...)
  for (i in 1:wn){
    graphics::points(xtime,datM[,7+i],cex=0.3,col=1+i,pch=1+i)
  }

  graphics::legend("bottomright",legend=c("orginal data",paste("weight",1:wn)),
                   col=1:(wn+1),pch = 1:(wn+1),cex=0.8)
  graphics::mtext(paste(paste0("Statistic",1:wn,sep=":"),round(x$stat.mat[,1],3),
                        collapse = "  "),side=3,cex=0.7)
}

