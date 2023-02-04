
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

