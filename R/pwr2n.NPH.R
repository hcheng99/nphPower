pwr2n.NPH<- function(method = "MaxLR"
                     ,entry   = 1
                     ,fup      = 1
                     ,maxfup = entry+fup
                     ,CtrlHaz
                     ,hazR
                     ,transP1 = c(0, 0)
                     ,transP0 = c(0, 0)
                     ,Wlist = list(function(x){ x^0})
                     ,entry_pdf0=function(x){(1/entry)*(x>=0&x<=entry)}
                     ,entry_pdf1=entry_pdf0
                     ,ratio    = 1
                     ,alpha    = 0.05
                     ,beta     = 0.1
                     ,alternative = c("two.sided")
                     ,criteria = 500
                     ,k        = 100
                     ,m = 0
                     ,nreps = 10
                     ,summary = TRUE

){

  if (!method %in% c("MaxLR","Projection")){
    stop("The 'method' must be one of 'MaxLR','Projection'.")
  }
  if (!alternative %in% c("two.sided","greater","less")){
    stop("The 'alternative' must be one of 'two.sided','greater','less'.")
  }

  tot_time <- maxfup
  num <- k*tot_time
  # create the subintervals
  if (m >1 |m <0){stop("m must be within 0 and 1")}
  x <- seq(1/k*m,tot_time,by=1/k)
  # if control hazard is zero at x=0, then use a small value for x like
  # 0.0001
  ## 2022-01-19 or control hazard is Inf
  if (m==0){
    if (CtrlHaz(0)==0|CtrlHaz(0)==Inf|hazR(0)==0|hazR(0)==Inf){
      x[1] <- 1/k*0.1}
  }

  ctrlRate <- CtrlHaz(x)
  haz_val <- hazR(x)*ctrlRate
  haz_point <- x*k
  ## load the transition matrix
  load <- trans.mat(numN=num,x=x,ctrlRate=ctrlRate,haz_val=haz_val,
                    haz_point=haz_point,ratio=ratio,
                    transP1=transP1,transP0=transP0,k=k,
                    fup=fup,entry=entry,entry_pdf0=entry_pdf0,
                    entry_pdf1=entry_pdf1,hazR=hazR,tot_time=tot_time)

  old.nrow <- nrow(load$pdat)
  pdat <- na.omit(load$pdat)
  if (nrow(pdat)<old.nrow){
    message(paste("At least",old.nrow-nrow(pdat), "NAs are produced in the transition matrix. Please Check
                  whether initial hazards are zero."))
    num <- nrow(pdat)
  }
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
    W[,j] <- Wlist[[j]](1-pdat$S) # to be consistent with MWLR test
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
    if (method=="MaxLR"){
      rho_est <- stats::cov2cor(Vmat)
      mu <- as.vector(dnum*t(W)%*%(pdat$rho*pdat$gamma))/sqrt(diag(Vmat))
      if (alternative=="two.sided"){
        ftwo <- function(i){
          #set.seed(i)
          crit <- mvtnorm::qmvnorm(1-alpha,tail="both.tails",
                                   mean=rep(0,wn),sigma = rho_est)$quantile
          power <- 1-mvtnorm::pmvnorm(-crit,crit,mean=mu,sigma = rho_est)
          return(c(crit,power))
        }
        ftwod <- apply(do.call(rbind,sapply(1:nreps,ftwo,simplify=FALSE)) ,2,mean)
        power <- ftwod[2]
        crit <- ftwod[1]
      }else if (alternative=="less"){ # l1 <l0
        ftwo <- function(i){
          # set.seed(i)
          crit <- mvtnorm::qmvnorm(1-alpha,tail="upper.tail",mean=rep(0,wn),sigma = rho_est)$quantile
          power <- 1-mvtnorm::pmvnorm(crit,Inf,mean=mu,sigma = rho_est)[1]
          return(c(crit,power))
        }
        ftwod <- apply(do.call(rbind,sapply(1:nreps,ftwo,simplify=FALSE)) ,2,mean)
        crit <- ftwod[1];  power <- ftwod[2]
      }else if (alternative=="greater"){
        ftwo <- function(i){
          #set.seed(i)
          crit <- qmvnorm(1-alpha,tail="lower.tail",mean=rep(0,wn),sigma = rho_est)$quantile
          power <- 1-pmvnorm(-Inf,crit,mean=mu,sigma = rho_est)[1]
          return(c(crit,power))
        }
        ftwod <-  apply(do.call(rbind,sapply(1:nreps,ftwo,simplify=FALSE)) ,2,mean)
        crit <- ftwod[1];  power <- ftwod[2]
      }
    }
    else if (method=="Projection"){
      if (alternative!="two.sided"){
        cat(c("note: only two-sided is supported for projection test."))

      }
      ## get the rank of the variance matrix
      mu <- as.vector(dnum*t(W)%*%(pdat$rho*pdat$gamma))
      vr <- qr(Vmat)$rank
      crit <- stats::qchisq(1-alpha,df=vr)
      ## get the noncentral parameter
      lmd <- t(mu)%*%MASS::ginv(Vmat)%*%mu
      power <- stats::pchisq(crit,df=vr,ncp=lmd,lower.tail = FALSE)
    }

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
  ## 2022-01-19
  eprob1 <- pdat$E_E[num]; eprob0 <- pdat$C_E[num]
  inparam <- c("Method","Entry Time", "Follow-up Time",
               "Allocation Ratio", "Type I Error", "Type II Error",
               "Alternative","Number of Weights")
  inval <- c(method,entry,fup,ratio, alpha, beta,
             alternative,length(Wlist))
  inputdata <- data.frame(parameter=inparam, value=inval)
  outparam <- c("Number of Events", "Number of Total Sampe Size",
                "Asymptotic Power", "Overall Event Rate")
  outval <- c(round(c( dnum, Nsize, as.numeric(power),eprob),digits = 3))
  outputdata <- data.frame(parameter=outparam, value=outval)

  if (summary==TRUE){
    cat("-----Summary of the Input Parameters----- \n")

    print(inputdata, row.names = FALSE)
    cat("-----Summary of the Output Parameters----- \n ")

    print(outputdata, row.names = FALSE)
  }
  summaryoutput <- list(input = inputdata, output = outputdata)

  inputfun <-list(
    alpha = alpha,
    beta = beta,
    transP1 = transP1,
    transP0 = transP0,
    k= k,
    criteria = criteria,
    controalhazard = CtrlHaz,
    hazardratio = hazR,
    entrypdf0  = entry_pdf0,
    entry_pdf1 = entry_pdf1,
    Weightfunctions = Wlist
  )
  names(event) <- paste0("w", 1:wn)
  listall <- list( eventN  = dnum
                   ,totalN = Nsize
                   ,pwr = as.numeric(power)
                   ,prob_event =eprob
                   ,prob1 = eprob1
                   ,prob0 = eprob0
                   ,L_trans = load$L_trans
                   ,pdat = pdat
                   ,studytime=c(entry,fup)
                   ,RandomizationRatio=ratio
                   ,eventList = event
                   ,inputfun = inputfun
                   ,summaryout = summaryoutput

  )
  class(listall) <-"NPHpwr"
  return(listall)


}
summary.NPHpwr <- function(object,...){
  summary <- object$summaryout
  wl  <- lapply(object$inputfun$Weightfunctions,deparse)
  wlc <- paste0(noquote(unlist(noquote(lapply(wl,"[",3)))),collapse=";")
  wlc <- gsub("\\s","",wlc)
  input <- data.frame(Parameter = c(summary$input$parameter,"Weight Functions"),
                      Value = c(summary$input$value, wlc))

  output <- summary$output
  names(input) <- c("__Parameter__", "__Value__")
  names(output) <- c("__Parameter__", "__Value__")
  cat("------------------------------------------ \n ")
  cat("----Summary of the Input Parameters---- \n")
  cat("------------------------------------------ \n ")

  print(input, row.names = FALSE, right=FALSE, justify = "left")
  cat("------------------------------------------ \n ")
  cat("-----Summary of the Output Parameters----- \n ")
  cat("------------------------------------------ \n ")
  print(output, row.names = FALSE, right=FALSE,  justify = "left")
  cat("------------------------------------------ \n ")
  cat("Notes: the base (x) of weight function is \n K-M estimate of CDF ")
}
plot.NPHpwr<- function(x,type=c("hazard","survival","dropout","event","censor"),...) {
  datM <- x$pdat
  totalN <- x$totalN
  ratio <- x$RandomizationRatio
  tval <-1
  if( missing(type)){ tval <- 0}
  ## draw the survival curves
  if (tval==0|"survival" %in% type ){
    with(datM,{
      plot(ti,S1,cex=0.1,lty=1,ylim=c(0,1),col=1,xlab="Analysis Time",
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
      plot(ti,hazard_E,cex=0.1,lty=1,col=1,xlab="Analysis Time",ylab="Hazard Rate",
           ylim=c(0,ymax),main="Hazard Curves",...)
      graphics::lines(ti,hazard_C,col=2,lty=2,...)
      graphics::legend("bottomright",legend=c("treatment","control"),
                       col=1:2,lty=1:2,cex=0.8)
    })
    ## hazard ratio
    with(datM,{
      plot(ti,theta,cex=0.1,lty=1,col=1,xlab="Analysis Time",
           ylim=c(min(theta)*0.9,max(theta)*1.1),
           ylab="Hazard Ratio (treatment over control)",
           main="Hazard Ratio",...)


    })

  }
  ## drop out
  if (tval==0|"dropout" %in% type){
    ymax <- max(c(datM$E_L,datM$C_L))*1.1
    with(datM,{
      plot(ti,E_L,cex=0.1,lty=1,col=1,xlab="Analysis Time",ylab="Proporion of drop out",
           ylim=c(0,ymax),main="Drop-out",...)
      graphics::lines(ti,C_L,col=2,lty=2,...)
      graphics::legend("topleft",legend=c("treatment","control"),
                       col=1:2,lty=1:2,cex=0.8)
    })
  }
  ## censor
  if (tval==0|"censor" %in% type){
    ymax <- max(c(datM$E_C,datM$C_C))*1.1
    with(datM,{
      plot(ti,E_C,cex=0.1,lty=1,col=1,xlab="Analysis Time",
           ylab="proporion of censor",
           ylim=c(0,ymax),main="Administrative Censoring",...)
      graphics::lines(ti,C_C,col=2,lty=2,...)
      graphics::legend("topleft",legend=c("treatment","control"),
                       col=1:2,lty=1:2,cex=0.8)
    })
  }
  ##event number
  if (tval==0|"event" %in% type ){
    with(datM,{
      plot(ti,round(eprob*totalN,digits=0),cex=0.1,lty=1,col=1,xlab="Analysis Time",
           ylab="Number of events",
           main="Events",...)
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
