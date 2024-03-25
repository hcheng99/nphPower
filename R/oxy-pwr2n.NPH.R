

###input parameters
#' @title Sample Size Calculation with Combination Test
#'
#' @description \code{pwr2n.NPH} calculates the number of events and
#' subjects required to achieve pre-specified power in the setup of two groups.
#' The method extends the calculation in the framework of the Markov model by Lakatos, allowing
#' for using the maximum weighted logrank tests or projection test with an arbitrary number of weight
#' functions. For maximum weighted logrank type test, if only one weight function
#' is provided, the test is essentially
#' the classic (weighted) logrank test.
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
#' of enrollment time for control group. The default assumes a uniform distribution
#' corresponding to the constant enrollment rate.
#' Default: function(x) {(1/entry) * (x >= 0 & x <= entry)}
#' @param entry_pdf1 a pdf of enrollment time for treatment
#' group. See \code{entry_pdf0}, Default: assume same pdf as control group.
#' @param ratio an integer, indicating the randomization ratio between treatment
#' and control group, Default: 1
#' @param alpha type I error rate, Default: 0.05
#' @param beta type II error rate, Default: 0.1
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "\code{two.sided}", "\code{greater}","\code{less}". See details.
#' For \code{"Projection"} method, only \code{"two-sided"}
#' alternative is supported. Default: c("\code{two.sided}")
#' @param criteria an integer indicating the maximum iteration allowed in
#' obtaining the number of events. See details , Default: 500
#' @param k an integer, indicating number of sub-intervals per time unit,
#'  Default: 100
#' @param m a value within 0 and 1.
#' @param nreps an integer, indicating number of iterations in calculating
#' the quantile of multivariate normal. See Details.
#' @param varianceType Default: \code{equal}. Indicates different variance assumptions
#'  for the sample size calculation.
#' It is not applicable for the maximum weighted logrank test. See details.
#' @param weightBase A character, either "\code{F}" or "\code{T}". \code{F} indicates
#' a CDF is the base for the weight function used in the weighted logrank or maximum weighted
#' logrank test. \code{T} indicates time is the base for weight function. Default: \code{F}
#' @param summary a logical value, controlling whether to print the
#' summary of calculation, Default: TRUE
#' @return
#' An object of class "\code{NPHpwr}" with corresponding \code{plot} function.
#' The object is a list containing the following components:
#'  \item{eventN}{total number of events}
#'  \item{totalN}{total number of subjects}
#'  \item{pwr}{actual power given the number of events}
#'  \item{prob_event}{event probability at the end of trial}
#'  \item{prob1}{event probability for the treatment group}
#'  \item{prob0}{event probability for the control group}
#'  \item{L_trans}{a list, consisting of transition matrix at each interval}
#'  \item{pdat}{ a dataframe including all the intermediate variables in the calculation.
#'  see Details.}
#'  \item{studytime}{a vector of length 2, including the entry and follow-up time as input}
#'  \item{RandomizationRatio}{as input}
#'  \item{eventlist}{a vector containing the number of events using each
#'  weight function alone}
#'  \item{inputfun}{a list containing all the input functions specified by
#'  users}
#' @details
#' The detailed methods can be found in the reference papers. The number
#' of subjects is determined by several factors, including the control hazard function,
#' hazard ratio function, entry time distribution, follow-up time, etc.  Under proportional
#' hazard assumption, the number of events is mainly determined by the hazard ratio besides
#' type i/ii error rates. However, under nonproportional hazards, all the above design parameters
#' may have an impact on the number of events.
#' The study design assumes \code{entry} time units of
#' enrollment and at least \code{fup} time units of follow-up. If enrollment
#' time \code{entry} is set to zero, all subjects are enrolled simultaneously,
#' so there is no staggered entry. Otherwise, if
#' \code{entry} is greater than 0, administrative censoring is considered. The user-defined
#'enrollment time function, hazard function for the control group and hazard ratio function can be either discrete or continuous.
#'Various non-proportional hazards types are accommodated. See examples below.
#'If multiple weight functions are provided in \code{Wlist}, a maximum weighted logrank
#'test or combination test is implemented. An iterative procedure
#'is used to obtain the event number based on the multivariate normal distribution.  Package
#'\pkg{mvtnorm} is used to calculate the quantiles. Because the algorithm is slightly
#'seed dependent, the quantiles are mean values of ten replicates by default. The
#'number of replicates is controlled by argument \code{ninter}.
#'
#'The "\code{alternative}" option supports both two-sided and one-sided test.
#' Let \eqn{\Lambda_1} and \eqn{\Lambda_0} denote the cumulative hazard of
#' treatment and control group. The \code{less} option tests
#' \eqn{H_0: \Lambda_1 > \Lambda_0} against
#' \eqn{H_a: \Lambda_1 <= \Lambda_0}. The \code{greater} option tests
#'\eqn{H_0: \Lambda_1 < \Lambda_0} against \eqn{H_a: \Lambda_1 >= \Lambda_0}.
#'
#'When \code{varianceType} is \code{equal}, the sample size for a two sided test is
#'\eqn{(z_{1-\alpha/2}+z_{1-\beta})^2\tilde{\sigma}^2/\mu_w^2}, where \eqn{\tilde{\sigma}^2}
#' is the variance estimate under alternative. when \code{varianceType} is not
#' \code{equal}. The formula is \eqn{(z_{1-\alpha/2}\sigma_w+z_{1-\beta}\tilde{\sigma})^2/\mu_w^2}.
#' Please use \code{equal} variance type for the maximum weighted logrank test.
#'
#' @references
#'
#' Brendel, M., Janssen, A., Mayer, C. D., & Pauly, M. (2014). Weighted logrank
#' permutation tests for randomly right censored life science data. Scandinavian
#' Journal of Statistics, 41(3), 742-761.
#'
#' Cheng, H., & He, J. (2021). A Maximum Weighted Logrank Test in Detecting
#' Crossing Hazards. arXiv preprint arXiv:2110.03833.
#'
#'
#' Cheng H, He J. Sample size calculation for the combination test under
#' nonproportional hazards. Biom J. 2023 Apr;65(4):e2100403. doi: 10.1002/bimj.202100403.
#' Epub 2023 Feb 15. PMID: 36789566
#' @examples
#'  #------------------------------------------------------------
#'  ## Delayed treatment effects using maxcombo test
#'  ## generate a list of weight functions for maxcombot test
#'  wmax <- gen.wgt(method = "Maxcombo" )
#'  t_enrl <- 12; t_fup <- 18; lmd0 <- log(2)/12
#'  ## delayed treatment effects
#'  f_hr_delay <- function(x){(x<=6)+(x>6)*0.75}
#'  f_haz0 <- function(x){lmd0*x^0}
#'  ##  The following code takes more than 5 seconds to run
#'  \donttest{
#'  snph1 <- pwr2n.NPH(entry = t_enrl, fup = t_fup, Wlist = wmax,
#'                     k = 100, ratio = 2, CtrlHaz = f_haz0, hazR = f_hr_delay)
#'  }
#'  #-------------------------------------------------------------
#'  # same setting using projection test
#'  snph2 <- pwr2n.NPH(method = "Projection", entry = t_enrl,
#'   fup = t_fup, Wlist = wmax, k = 10, ratio = 2, CtrlHaz = f_haz0,
#'   hazR = f_hr_delay)
#'
#'  #-------------------------------------------------------------
#'  #proportional hazards with weibull survival for control group
#'  #logrank test
#'  wlr <- gen.wgt(method = "LR" )
#'  b0 <- 3
#'  th0 <- 10/(-log(0.2))^(1/b0)
#'  #Weibull hazard function
#'  f_hz_weibull <- function(x){b0/th0^b0*x^(b0-1)}
#'  #hazard ratio function
#'  f_hr <- function(x){0.5*x^0}
#'  # define entry and follow-up time
#'  t_enrl <- 5; t_fup <- 5
#'  exph1 <- pwr2n.NPH(entry = t_enrl, fup = t_fup, k = 100,
#'  Wlist = wlr,  CtrlHaz = f_hz_weibull, hazR = f_hr, summary = FALSE)
#'  summary(exph1)
#' @seealso
#'  \code{\link{pwr2n.LR}}
#'  \code{\link{gen.wgt}}, \code{\link{evalfup}}
#' @rdname pwr2n.NPH
#' @export
#' @importFrom stats cov2cor weighted.mean na.omit
#' @importFrom mvtnorm qmvnorm pmvnorm

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
                     ,varianceType = c("equal") # not applicable for combination test
                     ,weightBase = "C" # C for CDF ； T for time
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
  # print("step1:")
  ## calculate the variance-covariance matrix
  Vmat <- matrix(NA,nrow=wn,ncol=wn)
  event <- c()
  if (alternative=="two.sided"){
    div <- 2
    #part2=(qnorm(1-alpha/2)+qnorm(1-beta))^2
  }else {
    #part2=(qnorm(1-alpha)+qnorm(1-beta))^2
    div <- 1
  }
  part2=(qnorm(1-alpha/div)+qnorm(1-beta))^2
  for (j in 1:wn){
    if (weightBase == "C"){
      W[,j] <- Wlist[[j]](1-pdat$S) # to be consistent with MWLR test
    } else if (weightBase == "T"){
      W[,j] <- Wlist[[j]](pdat$ti)
    }

    wgt <- W[,j]

    part3=t(pdat$rho)%*%(pdat$gamma*wgt)
    part3=part3^2
    if (varianceType=="equal"){
      part1=t(pdat$rho)%*%(pdat$eta*wgt^2)
      event[j] <- part1*part2/part3
    } else{
      # different variance for the sample size
      pdat$eta_new <- ratio/(ratio+1)^2
      part1_1 <- t(pdat$rho)%*%(pdat$eta*wgt^2)
      part1_2 <- t(pdat$rho)%*%(pdat$eta_new*wgt^2)
      event[j] <- (qnorm(1-alpha/div)*sqrt(part1_2)+qnorm(1-beta)*sqrt(part1_1))^2/part3
    }
  }
  #print("step2:")
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
  #print("step3:")
  #print(c(event, power,dnum, wn))
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
  # print("step4:")
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
  #print("step5:")
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
  #print("step6:")
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
#' @title Summary of the \code{pwr2n.NPH} function
#' @description Summarize and print the results of \code{pwr2n.NPH} function
#' @param object object of the \code{pwr2n.NPH} function
#' @param ... additional arguments passed to the summary function
#' @return
#' No return value. Summary results are printed to Console.
#' @seealso
#'  \code{\link{pwr2n.NPH}}
#' @rdname summary.NPHpwr
#' @export
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

#*********************************************
#*show the survival plot/ hazards plots
#*********************************************
#' @title Graphical Display of Design Parameters in Sample Size
#' Calculation
#' @description Displays graphs of survival, hazards, drop-out and censor over time
#' as specified in the calculation.
#' @param x object of the \code{pwr2n.NPH} function
#' @param type a vector of string, specifying the graphs to display. The options
#' include "\code{hazard}", "\code{survival}", "\code{dropout}",
#'  "\code{event}", and
#'  "\code{censor}". If \code{type} is not provided, all the available graphs are
#' generated.
#' @param ... additional graphical arguments passed to the plot function
#' @return
#' plots are produced on the current graphics device
#' @details
#' The \code{type} argument provides five options to visualize the trial in design.
#' Option \code{survival} shows the survival probabilities of treatment and control
#' group over time. Option \code{hazard} provides the hazard rates and hazard ratio
#' over time. Option \code{dropout} shows the proportion of drop-out subjects across
#' the trial duration.
#' Option \code{censor} shows the proportion of censored subjects over time.
#' @examples
#'  # generate weight function
#'  wlr <- gen.wgt(method = "LR" )
#'  t_enrl <- 12; t_fup <- 18; lmd0 <- log(2)/12
#'  # delayed treatment effects, the crossign point is at 6.
#'  f_hr_delay <- function(x){(x<=6)+(x>6)*0.75}
#'  f_haz0 <- function(x){lmd0*x^0}
#'  snph1 <- pwr2n.NPH(entry = t_enrl, fup = t_fup, Wlist = wlr,
#'                     k = 100, ratio = 2, CtrlHaz = f_haz0,
#'                     hazR = f_hr_delay)
#'  # display the hazards plot only
#'  plot(snph1, type="hazard")
#'  # display all plots
#'  plot(snph1)
#' @seealso
#'  \code{\link{pwr2n.NPH}}
#' @rdname plot.NPHpwr
#' @export
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

