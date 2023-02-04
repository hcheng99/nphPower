
#' @title Simulate Survival Trial Data
#' @description \code{simu.trial} simulates survival data allowing flexible input
#' of design parameters. It supports both event-driven design and fixed study duration
#' design.
#' @param type indicates whether event-driven trial ("\code{event}) or fixed study duration
#' trial ("\code{time}"), Option: c("\code{event}", "\code{time}")
#' @param trial_param a vector of length 3 with components for required subject size, enrollment
#' time and required number of events ("\code{event}" type trial)/follow-up time
#' ("\code{time}" type trial)
#' @param bsl_dist indicates the survival distribution for control group, option:
#'  c("\code{weibull}", "\code{loglogistic}")
#' @param bsl_param a vector of length 2 with the shape and rate (scale) parameter for the
#'  survival distribution of control group. See details.
#' @param drop_param0 a vector of length 2 with shape and scale parameter for the
#' weibull distribution of drop-out time for control group
#' @param drop_param1 a vector of length 2 with shape and scale parameter for the
#' weibull distribution of drop-out time for treatment group
#' @param entry_pdf0 a function describing the pdf of the entry time for control. Default: uniform enrollment
#' @param entry_pdf1 a function describing the pdf of the entry time for treatment.
#' @param HR_fun a function describing the hazard ratio function between treatment
#' and control group
#' @param ratio allocation ratio between treatment and control group.
#' For example, \code{ratio}=2 if 2:1 allocation is used.
#' @param upInt a value indicating the upper bound used in the \code{uniroot} function.
#' See details.  Default: 100
#' @param summary a logical indicating whether basic information summary is printed
#' to the console or not, Default: TRUE
#' @return
#' A list containing the following components
#'
#' \code{data}: a dataframe (simulated dataset) with columns:
#' \describe{
#' \item{id}{identifier number from 1:n, n is the total sample size}
#' \item{group}{group variable with 1 indicating treatment and 0 indicating control}
#' \item{aval}{observed survival time, the earliest time among event time,
#' drop-out time and end of study time}
#' \item{cnsr}{censoring indicator with 1 indicating censor and 0 indicating event}
#' \item{cnsr.desc}{description of the \code{cnsr} status, including three options-
#' drop-out, event and end of study. Both drop-out and end of study are considered as
#' censor.}
#' \item{event}{event indicator with 1 indicating event and 0 indicating censor}
#' \item{entry.time}{time when the patient is enrolled in the study}
#'}
#' a list of summary information of the trial including
#'
#' \describe{
#' \item{\code{type}}{a character indicating input design type - \code{event} or \code{time}}
#' \item{\code{entrytime}}{{a number indicating input enrollment period}}
#' \item{\code{maxob}}{a number indicating the maximum study duration. For \code{time}
#' type of design, the value is equal to the enrollment period plus the follow-up
#' period. For \code{event} type of design, the value is the calendar time of the
#' last event.}
#' }
#'
#'@details
#' The loglogistic distribution for the event time has the
#' survival function \eqn{S(x)=b^a/(b^a+x^a)} and hazard function
#' \eqn{\lambda(x)=a/b(x/b)^(a-1)/(1+(t/b)^a)},, where \eqn{a} is the shape parameter
#' and \eqn{b} is the scale parameter. The weibull distribution for event time and drop-out time has the survival function \eqn{f(x)=exp(-(xb)^a)}
#' and hazard function \eqn{\lambda(x)=ab(xb)^a-1}, where \eqn{a} is the shape parameter
#' and \eqn{b} is the rate parameter. The median of weibull
#' distribution is \eqn{(ln(2)^(1/a)/b)}. If drop out or loss to follow-up are
#' do not need to be considered, a very small rate parameter \eqn{b} can be chosen
#' such that the median time is greatly larger than the study duration. The default
#' entry time is uniformly distributed within the enrollment period by default.
#' Other options are allowed through providing the density function.
#'
#' The \code{simu.trial} function simulates survival times for control and
#' treatment groups separately. The control survival times are drawn from standard parametric
#' distribution (Weibull, Loglogistic). Let \eqn{\lambda_1(t)} and \eqn{\lambda_0(t)}
#' denote the hazard function for treatment and control. It is assumed that
#' \eqn{\lambda_1(t)/\lambda_0(t)=g(t)}, where \eqn{g(t)} is the user-defined
#' function, describing the change of hazard ratio over time. In case of proportional
#' hazards, \eqn{g(t)} is a constant. The survival function for treatment group
#' is \eqn{S_1(t)=exp(-\int_0^t g(s)\lambda_0(s)ds)}. The Survival time T is
#' given by \eqn{T=S_1^(-1)(U)}, where U is drawn from uniform (0,1). More details
#' can be found in Bender, et al. (2005). Function \code{uniroot} from
#' \code{stats} package is used to generate the random variable. The argument
#' \code{upInt} in \code{simu.trial} function corresponds to the upper end point
#' of the search interval and it can be adjusted by user if the default value 100
#' is not appropriate. More details can be found in help file of \code{uniroot}
#' function.
#'
#' @references
#' Bender, R., Augustin, T., & Blettner, M. (2005). Generating survival times to simulate Cox proportional
#' hazards models. Statistics in medicine, 24(11), 1713-1723.

#' @examples
#' # total sample size
#' N <- 300
#' # target event
#' E <- 100
#' # allocation ratio
#' RR <- 2
#' # enrollment time
#' entry <- 12
#' # follow-up time
#' fup <- 18
#' # shape and scale parameter of weibull for entry time
#' b_weibull <- c(1,log(2)/18)
#' # shape and scale parameter of weibull for drop-out time
#' drop_weibull <- c(1,log(2)/30)
#' # hazard ratio function (constant)
#' HRf <- function(x){0.8*x^0}
#'
#' ### event-driven trial
#' set.seed(123445)
#' data1 <- simu.trial(type="event",trial_param=c(N,entry,E),bsl_dist="weibull",
#'                     bsl_param=b_weibull,drop_param0=drop_weibull,HR_fun=HRf,
#'                     ratio=RR)
#'
#' ### fixed study duration
#' set.seed(123445)
#' data2 <- simu.trial(type="time",trial_param=c(N,entry,fup),bsl_dist="weibull",
#'                     bsl_param=b_weibull,drop_param0=drop_weibull,HR_fun=HRf,
#'                     ratio=RR)
#' @seealso
#'  \code{\link[stats]{uniroot}}
#' @rdname simu.trial
#' @export
#' @importFrom stats rweibull integrate rlogis uniroot runif
simu.trial <- function(type=c("event","time")
                       ,trial_param # include the total sample size,entry time,
                       # target event (event type)/fup time (time )
                       ,bsl_dist=c("weibull","loglogistic")
                       ,bsl_param   # alpha=1 corresponds to exponential
                       ,drop_param0
                       ,drop_param1=drop_param0
                       ,entry_pdf0=function(x){(1/trial_param[2])*(x>=0&x<=trial_param[2])}
                       ,entry_pdf1=entry_pdf0
                       ,HR_fun #the non proportion hazard function
                       ,ratio # # of trt/# of placebo
                       ,upInt=100
                       ,summary=TRUE
){

  if (length(trial_param) !=3) {stop("The trial parameters must include
                                     total sample size, entry time, targeted events/
                                     follow-up time ")
  }else {
    t_p1 <- trial_param[1]
    t_p2 <- trial_param[2]
    t_p3 <- trial_param[3]
  }
  if (!type %in% c("event","time")){stop("The type must be in ('event','time')")}

  # to note: n1 is the # of subjects in treatment group, which is assumed to
  # have fewer events
  prop <- ratio/(ratio+1)
  n1 <- ceiling(t_p1*prop)
  n0 <- t_p1-n1
  if (ceiling(n1)!=n1|ceiling(n0)!=n0){
    stop("The number of subjects in each group must be interger. Check!")
  }
  trt <- c(rep(0,n0),rep(1,n1))
  a=bsl_param[1]
  b=bsl_param[2]
  #****************************
  #* Simulate Event Time
  #****************************
  if (bsl_dist=="weibull"){
    #set.seed(seed)
    T_0 <- stats::rweibull(n0,a,1/b)
    #-- get the cummulative hazard and survival function----#
    Hf <- function(t){exp(-1* a*b*stats::integrate( function(x){(x*b)^(a-1)*HR_fun(x)},0,t)$value)}
  }else if (bsl_dist=="loglogistic"){
    #set.seed(seed)
    T_0 <- exp(stats::rlogis(n0,log(b),1/a))
    Hf <- function(t){exp(-1* a/b*integrate( function(x){(x/b)^{a-1}/(1+(x/b)^a)*HR_fun(x)},0,t)$value)}
  }
  gen_t <- function(y){stats::uniroot(function(x){Hf(x)-y},interval = c(0,upInt),extendInt="yes")$root}
 # set.seed(seed*10)
  U1 <- stats::runif(n1)
  T_1 <-as.vector(unlist(lapply(U1, gen_t)))
  Tm<-c(T_0,T_1)
  ## in extreme case, time is 0, add 0.1;
  #****************************
  #* Simulate drop-out Time
  #****************************
  if (missing(drop_param0)&missing(drop_param1)){
    cat("Notes: No drop-out parameters are provided.
Drop-out is not considered in the simulation.")

  }else if (missing(drop_param0)){
    drop_param0 <- drop_param1
  }else if (missing(drop_param1)){
    drop_param1 <- drop_param0
  }
  Tm <- ifelse(Tm==0,0.1,Tm)
  if (!missing(drop_param0)&!missing(drop_param1)){
    a0_drop <- drop_param0[1]
    b0_drop <- drop_param0[2]
    a1_drop <- drop_param1[1]
    b1_drop <- drop_param1[2]
    #set.seed(seed*99)
    drop_time0 <- stats::rweibull(n0,a0_drop,1/b0_drop)
    drop_time1 <- stats::rweibull(n1,a1_drop,1/b1_drop)
    drop_time <- c(drop_time0,drop_time1)

  }else{
    drop_time <- rep(Inf,t_p1)
  }
  #****************************
  #* Simulate entry Time
  #****************************
  #- create the CDF;
  ent_cdf0 <- function(t){stats::integrate(entry_pdf0,lower=0,upper=t)$value}
  gen_ent0 <- function(y){stats::uniroot(function(x){ent_cdf0(x)-y},
                                  interval = c(0,upInt),extendInt="yes")$root}
  ent_cdf1 <- function(t){stats::integrate(entry_pdf1,lower=0,upper=t)$value}
  gen_ent1 <- function(y){stats::uniroot(function(x){ent_cdf1(x)-y},
                                  interval = c(0,upInt),extendInt="yes")$root}

  #set.seed(seed+1)
  tu0_0 <- runif(n0)
  tu0_1 <- runif(n1)
  t0_0 <-as.vector(unlist(lapply(tu0_0, gen_ent0)))
  t0_1 <-as.vector(unlist(lapply(tu0_1, gen_ent1)))
  t0 <- c(t0_0,t0_1)
  ot <- t0+Tm

  dat <- data.frame(id=1:t_p1,ent=t0,time=Tm,trt=trt,ot=ot,
                    drop_time=drop_time,ot_drop=t0+drop_time)
  if (type=="event"){
    # find the smallest time between drop-out and event time
    min_ind0 <- apply(cbind(dat$time,dat$drop_time),1,which.min)
    dat$i1 <- min_ind0==1
    dat <- dat[order(dat$ot),]
    dat$c0 <- cumsum(dat$i1)
    if (max(dat$c0)<t_p3){stop(paste("The target event "),t_p3,
                               " cannot achieve. Please check the parameters
                               for event, entry and drop-out parameters")}

    Dur <- min(dat[dat$c0==t_p3,]$ot )
    }else{
    Dur <- t_p2+t_p3 # the length of study
  }
  min_ind <- apply(cbind(dat$ot,Dur,dat$ot_drop),1,which.min)
  status <- c("event","end of study","drop-out")
  tot_len <- cbind(dat$ot,Dur,dat$ot_drop)[cbind(seq_along(min_ind),min_ind)]
  dat$t_val <- tot_len-dat$ent
  dat$cnsr_desc <- status[min_ind]
  dat$cnsr <- 1-(min_ind==1) #cnsr=1 indicates censoring

  final <- with (dat,data.frame(
    id=id
    ,group=trt
    ,aval=t_val
    ,cnsr=cnsr
    ,cnsr.desc=cnsr_desc
    ,event=1-cnsr
    ,entry.time=ent
    ,event.time=time
    ,drop.time=drop_time
    ,obs.time=tot_len
  )
  )
  #table(final$group,final$cnsr.desc)
  if (summary==TRUE){
    ctext <- c("Trial Type:", "Entry Time:", "Maximum Study Duration:",
               "Number of Subjects:", "Number of Events:")
    cval <- c(type,t_p2,round(Dur,digits = 2),t_p1, sum(final$event))
    csum <- data.frame(parameter=ctext, value=cval)
    cat("\n -------- Summary of the Simulation -------- \n")
    print(csum)

  }

  list <- list(data=final,
               type=type,
               entrytime=t_p2,
               maxobs=Dur)
  class(list) <- 'SimuTrial'
  return(list)
}
