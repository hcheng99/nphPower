

###input parameters
pwr2n.prj<-function(entry   = 1
                       ,fup      = 1
                       ,CtrlHaz
                       ,hazR
                       ,transP1
                       ,transP0
                       ,Wlist = list(function(x){ x^0})
                       ,entry_pdf0=function(x){(1/entry)*(x>=0&x<=entry)}
                       ,entry_pdf1=entry_pdf0
                       ,ratio    = 1
                       ,alpha    = 0.05
                       ,beta     = 0.1
                       ,alternative = c("two.sided")
                       ,criteria = 500
                       ,k        = 100
){

  if (!alternative %in% c("two.sided","greater","less")){
    stop("The alternative must be one of 'two.sided','greater','less'.")
  }
  tot_time <- entry+fup
  num <- k*tot_time
  # create the subintervals
  x <- seq(0,tot_time,by=1/k)
  ctrlRate <- CtrlHaz(x)
  haz_val <- hazR(x)*ctrlRate
  haz_point <- x*k
  ## load the transition matrix
  load <- trans.mat(numT=num,x=x,ctrlRate=ctrlRate,haz_val=haz_val,
                    haz_point=haz_point,ratio=ratio,
                    transP1=transP1,transP0=transP0,k=k,
                    fup=fup,entry=entry,entry_pdf0=entry_pdf0,
                    entry_pdf1=entry_pdf1)

  pdat <- load$pdat
  wn <- length(Wlist)
  W <- matrix(NA,nrow=nrow(pdat),ncol=wn)
  ## calculate the variance-covariance matrix
  Vmat <- matrix(NA,nrow=wn,ncol=wn)
  event <- c()
  if (alternative=="two.sided"){
    part2=(qnorm(1-alpha/2)+qnorm(1-beta))^2
  }else {
    part2=(qnorm(1-alpha)+qnorm(1-beta))^2
  }
  for (j in 1:wn){
    W[,j] <- Wlist[[j]](pdat$S)
    wgt <- W[,j]
    part1=t(pdat$rho)%*%(pdat$eta*wgt^2)
    part3=t(pdat$rho)%*%(pdat$gamma*wgt)
    part3=part3^2
    event[j] <- part1*part2/part3
  }
  #print(c(part1,part2,part3))
  ## get the min and max sample size
  E_min=min(event)
  E_max=max(event)
  eta=-1
  dnum <- E_min
  count <- 0
  while(eta <0){
    count <- count+1
    for (k1 in 1:wn){
      for (k2 in 1:wn){
        Vmat[k1,k2] <-  dnum*t(W[,k1]*W[,k2]) %*%(pdat$rho*pdat$eta)
      }
    }
    ## get the rank of the variance matrix
    mu <- as.vector(dnum*t(W)%*%(pdat$rho*pdat$gamma))
    vr <- qr(Vmat)$rank
    crit <- stats::qchisq(1-alpha,df=vr)
    ## get the noncentral parameter
    lmd <- t(mu)%*%MASS::ginv(Vmat)%*%mu
    power <- stats::pchisq(crit,df=vr,ncp=lmd,lower.tail = FALSE)
    if (power<1-beta&count<criteria) {dnum<-dnum+1}
    else if (count>=criteria){
      warning(paste0("the algoritm doesn't converge within ",criteria," iterations;
                     the current power is ",power,"; event size: ",dnum,"
                     please consider increase the criteria option."))
      ;break}
    else {break}

  }

  eprob <- stats::weighted.mean(c(pdat$C_E[num],pdat$E_E[num]),w=c(1,ratio))
  Nsize <- dnum/eprob
 print(c(dnum,Nsize))

  listall <- list( eventN  = dnum
                   ,totalN = Nsize
                   ,pwr = as.numeric(power)
                   ,prob_event =eprob
                   ,L_trans = load$L_trans
                   ,pdat = pdat
                   ,studytime=c(entry,fup)
                   ,RandomizationRatio=ratio

  )
  class(listall) <-"MaxLRpwr"
  return(listall)


}

#*********************************************
#*show the survival plot/ hazards plots
#*********************************************
plot.MaxLRpwr<- function(x,type=c("hazard","survival","dropout","event","censor"),...) {
  datM <- x$pdat
  totalN <- x$totalN
  ratio <- x$RandomizationRatio
  tval <-1
  if( missing(type)){ tval <- 0}
  ## draw the survival curves
  if (tval==0|"survival" %in% type ){
    with(datM,{
      plot(ti,S1,cex=0.1,lty=1,ylim=c(0,1),col=1,xlab="Time",
           ylab="Survival Probability",main="Survival Curves",...)
      graphics::lines(ti,S0,col=2,lty=2,...)
      graphics::legend("bottomleft",legend=c("treatment","control"),
                       col=1:2,lty=1:2,cex=0.8)

    })
  }

  ## hazard functions
  if (tval==0|"hazard" %in% type ){

    ymax <- max(c(datM$hazard_C,datM$hazard_E))*1.1
    with(datM,{
      plot(ti,hazard_E,cex=0.1,lty=1,col=1,xlab="Time",ylab="Hazard Rate",
           ylim=c(0,ymax),main="Hazard Curves",...)
      graphics::lines(ti,hazard_C,col=2,lty=2,...)
      graphics::legend("bottomright",legend=c("treatment","control"),
                       col=1:2,lty=1:2,cex=0.8)
    })
    ## hazard ratio
    with(datM,{
      plot(ti,theta,cex=0.1,lty=1,col=1,xlab="Time",
           ylim=c(min(theta)*0.9,max(theta)*1.1),
           ylab="Hazard Ratio (treatment over control)",
           main="Hazard Ratio over Time",...)


    })

  }
  ## drop out
  if (tval==0|"dropout" %in% type){
    ymax <- max(c(datM$E_L,datM$C_L))*1.1
    with(datM,{
      plot(ti,E_L,cex=0.1,lty=1,col=1,xlab="Time",ylab="proporion of drop out",
           ylim=c(0,ymax),main="Drop-out overtime",...)
      graphics::lines(ti,C_L,col=2,lty=2,...)
      graphics::legend("topleft",legend=c("treatment","control"),
                       col=1:2,lty=1:2,cex=0.8)
    })
  }
  ## censor
  if (tval==0|"censor" %in% type){
    ymax <- max(c(datM$E_C,datM$C_C))*1.1
    with(datM,{
      plot(ti,E_C,cex=0.1,lty=1,col=1,xlab="Time",
           ylab="proporion of censor",
           ylim=c(0,ymax),main="Administrative censoring overtime",...)
      graphics::lines(ti,C_C,col=2,lty=2,...)
      graphics::legend("topleft",legend=c("treatment","control"),
                       col=1:2,lty=1:2,cex=0.8)
    })
  }
  ##event number
  if (tval==0|"event" %in% type ){
    with(datM,{
      plot(ti,round(eprob*totalN,digits=0),cex=0.1,lty=1,col=1,xlab="Time",
           ylab="number of events",
           main="Events over time",...)
      graphics::lines(ti,round(E_E*totalN*ratio/(ratio+1),digits=0),col=2,lty=2,...)
      graphics::lines(ti,round(C_E*totalN/(ratio+1),digits=0),col=3,lty=3,...)
      graphics::legend("topleft",legend=c("overall","treatment","control"),
                       col=1:3,lty=1:3,cex=0.8)
    })
  }


  #
  #   graphics::mtext(paste(paste0("Statistic",1:wn,sep=":"),round(x$stat.mat[,1],3),
  #                         collapse = "  "),side=3,cex=0.7)
}
