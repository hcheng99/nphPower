projection.test <- function(dat
                        ,Wlist
                        ,base
                        ,alpha=0.05

){

  ## transform the input data to the matrix ready for analysis
  datM <- load.mat(dat)
  #---------------------create the weight --------------------------------;
  datM <- datM[datM$event>0,] # only keep the timepoints with event
  if (base=="Combined"){

    tnew0 <- c(1,
               cumprod(chkV(datM$n.risk.x,1-datM$n.event.x/datM$n.risk.x))*table(dat[,3])[2]/nrow(dat)+
                 cumprod(chkV(datM$n.risk.y,1-datM$n.event.y/datM$n.risk.y))*table(dat[,3])[1]/nrow(dat))
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
