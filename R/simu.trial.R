
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
    cat("Notes: Drop-outs are not considered in the simulation.")

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
