#' @title Projection test
#' @description Perform projection test as proposed by Brendel (2014)
#' @param dat a dataframe or matrix, of which the first three columns are survival
#' time, event status indicator  and group label. The status indicator, normally
#' 0=alive, 1=dead/event. Other choices are TRUE/FALSE (TRUE=death) or 1/2 (2=death).
#' The group label can be either numeric values like 0=control, 1=treatment or text
#' like C=control, T=treatment.
#' @param Wlist a list object with components of weight functions
#' @param base a text must be one of c("\code{KM}","\code{Combined}","\code{N}"), Default: c("KM")
#' @param alpha a number indicating type I error rate, Default: 0.05
#' @return
#' a list of components including
#' \item{chisq}{a numeric value indicating the chi-square statistic}
#' \item{df.chis}{a numeric value indicating the degree freedom of the test}
#' \item{pvalue}{a numeric value giving the p-value of the test }
#' \item{details}{a data frame consisting of statistics from multiple
#' weight functions and the variance-covariance matrix}
#' @details
#' The base functions are the same as those described in function
#' \code{MaxLRtest}. The method detail can be found in Brendel (2014)
#' paper. The main idea is to map the multiple weighted logrank statistics
#' into a chi-square distribution. The degree freedom of the chi-square
#' is the rank of the generalized inverse of covariance matrix. Only two-sided
#' test is supported in the current function.
#' @references
#' Brendel, M., Janssen, A., Mayer, C. D., & Pauly, M. (2014). Weighted logrank
#' permutation tests for randomly right censored life science
#' data. Scandinavian Journal of Statistics, 41(3), 742-761.
#' @examples
#' # load and prepare data
#' data(lung)
#' tmpd <- with(lung, data.frame(time=SurvTime,stat=1-censor,grp=Treatment))
#' # two weight functions are defined.
#' # one is constant weight; the other emphasize diverging hazards
#' timef1 <- function(x){1}
#' timef2 <- function(x){(x)}
#' test1 <- projection.test(tmpd,list(timef1,timef2),base="KM")
#' test1$chisq; test1$pvalue; test1$df.chisq
#' @seealso
#'  \code{\link{MaxLRtest}}
#' @rdname projection.test
#' @export
#' @importFrom survival survfit
#' @importFrom stats stepfun qchisq pchisq
#' @importFrom MASS ginv
projection.test <- function(dat
                            ,Wlist
                            ,base
                            ,alpha=0.05
                            ,side=c("two.sided")
){

  ## transform the input data to the matrix ready for analysis
  datM <- load.mat(dat)
  #---------------------create the weight --------------------------------;
  datM <- datM[datM$event>0,] # only keep the timepoints with event
  if (base=="Combined"){

    tnew0 <- c(1,
               cumprod(chkV(datM$n.risk.x,1-datM$n.event.x/datM$n.risk.x))*table(dat$V3)[2]+
                 cumprod(chkV(datM$n.risk.y,1-datM$n.event.y/datM$n.risk.y))*table(dat$V3)[1])
    tnew <- 1-tnew0[1:nrow(datM)]

  }else if (base=="KM"){
    # based on the pooled survival
    s_fit<- survival::survfit(Surv(dat[,1], dat[,2])~1 , data = dat)
    f_s <- stats::stepfun(s_fit$time,y=c(0,1-s_fit$surv),right = TRUE)
    tnew <- f_s(datM$time)
  }else if (base=="N"){
    # based on time/person at risk
    # tnew <- datM$time/max(datM$time)
    tnew <- 1-datM$risk/nrow(dat)
  }
  # number of weight functions
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
  Zstat <- apply(W,2,function(x){CalN(x,data=datM)})
  statistic <- t(Zstat)%*%MASS::ginv(Vmat)%*%Zstat
  crit <- stats::qchisq(1-alpha,qr(Vmat)$rank)

  pvalue <- stats::pchisq(statistic,qr(Vmat)$rank,lower.tail = FALSE)

  list <-list(
    chisq=statistic,
    df.chisq=qr(Vmat)$rank,
    pvalue=pvalue,
    details=cbind(Zstat,Vmat)
  )
  return(list)

}
