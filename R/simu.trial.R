simu.trial <- function(type=c("event","time")
                       ,trial_param # include the total sample size,entry time,
                       # target event (event type)/fup time (time )
                       ,bsl_dist=c("weibull","loglogistic", "mix-weibull")
                       ,bsl_param   # alpha=1 corresponds to exponential
                       ,drop_param0
                       ,drop_param1=drop_param0
                       ,entry_pdf0=function(x){(1/trial_param[2])*(x>=0&x<=trial_param[2])}
                       ,entry_pdf1=entry_pdf0
                       ,enrollmentType = NULL # "piecewise uniform"
                       ,entryP = list(10000, 1) # one person one time unit
                       ,HR_fun #the non proportion hazard function
                       ,ratio # # of trt/# of placebo
                       ,cureModel = NULL # c("PHCM", "PHCRM")
                       ,cureRate1 = NULL # the cure rate for experimental group under PHCM
                       ,HR_data = NULL # an object
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
  if (missing(HR_fun) & is.null(HR_data)){
    stop("Either HR_fun or HR_data must be provided.")
  } else if (!missing(HR_fun) & !is.null(HR_data)){
    warning("HR_fun is used and HR_data is ignored. ")
  }
  if (!is.null(HR_data) & bsl_dist !="weibull"){
    stop("The input of HR_data is only supported for weibull distribution.")
  }
  if (is.null(cureModel)){ cureModel <- "None"}
  if (cureModel == "PHCM" &is.null(cureRate1)){
    stop("cure rate for the experimental group must be provided if PHCM model is selected.")
  }
  if (bsl_dist == "mix-weibull"  &length(bsl_param)!=3){
    stop("bsl_param should must three parameters including shape and scale of weibull and cure rate for the contrl group.")
  }
  # parameters for the distribution
  a=bsl_param[1]
  b=bsl_param[2]

  #****************************
  #* Simulate Event Time
  #****************************
  # if HR_fun is provided
  if (!missing(HR_fun)){
    # to note: n1 is the # of subjects in treatment group, which is assumed to
    # have fewer events
    prop <- ratio/(ratio+1)
    n1 <- ceiling(t_p1*prop)
    n0 <- t_p1-n1
    if (ceiling(n1)!=n1|ceiling(n0)!=n0){
      stop("The number of subjects in each group must be interger. Check!")
    }
    trt <- c(rep(0,n0),rep(1,n1))
    # add pi for mix-weibull
    pi_1 <- -0.1
    if (bsl_dist=="weibull" |(bsl_dist == "mix-weibull" & cureModel == "PHCM")){
      #set.seed(seed)
      T_0 <- stats::rweibull(n0,a,1/b)
      #-- get the cummulative hazard and survival function----#
      Hf <- function(t){exp(-1* a*b*stats::integrate( function(x){(x*b)^(a-1)*HR_fun(x)},0,t)$value)}
    }else if (bsl_dist=="loglogistic"){
      #set.seed(seed)
      T_0 <- exp(stats::rlogis(n0,log(b),1/a))
      Hf <- function(t){exp(-1* a/b*integrate( function(x){(x/b)^{a-1}/(1+(x/b)^a)*HR_fun(x)},0,t)$value)}
    }else if (bsl_dist == "mix-weibull" | cureModel == "PHCRM" ){ # added on DEC 18, 2023
      pi_0 <- bsl_param[3]
      u_n0 <- stats::runif(n0)
      T_0 <- u_n0
      T_0[u_n0>pi_0] <- -((log((u_n0[u_n0>pi_0]-pi_0)/(1-pi_0)))^1/a)/b
      T_0[u_n0<=pi_0] <- 10^10

      # control hazard
      ctrlhr <- function(x) {(1-pi_0)*exp(-(x*b)^a)*a*b*(x*b)^(a-1)/(pi_0+(1-pi_0)*exp(-(x*b)^a))}
      Hf <- function(t){exp(-1* stats::integrate( function(x){ctrlhr(x)*HR_fun(x)},0,t)$value)}
      pi_1 <-ceiling(stats::optimize(Hf,interval=c(0,10^5))$objective*10000)/10000

    }
    gen_t <- function(y){stats::uniroot(function(x){Hf(x)-y},interval = c(0,upInt),extendInt="yes")$root}
    # set.seed(seed*10)
    U1 <- stats::runif(n1)
    U1_sub <- U1[U1>pi_1]
    T1_sub <- as.vector(unlist(lapply(U1_sub, gen_t)))
    T_1 <- U1
    T_1[U1>pi_1] <- T1_sub
    T_1[U1<=pi_1] <- 10^10

    Tm<-c(T_0,T_1)
  }
  # if HR_data is provided for weibull distribution
  if (!missing(HR_data)){
    if (!"grp" %in% colnames(HR_data)){
      stop("column grp must exist in HR_data.")
    }
    if (nrow(HR_data)!=t_p1){
      stop("the number of observations of HR_data must be the same as
            as the specified sample size")
    }
    HR_data <-HR_data[ order(HR_data$grp),]
    trt <- HR_data$grp
    # the hazard for each observation is
    n0 <- sum(HR_data$grp==0)
    n1 <- sum(HR_data$grp==1)
    Tm <- stats::rweibull(t_p1, a,1/(b*HR_data$HR^(1/a)))
  }

  #**************************
  # PHCM model
  #**************************
  # if (!is.null(cure)){
  if (bsl_dist == "mix-weibull" & cureModel == "PHCM"){

    ucure <- runif(n0+n1)
    # the following method is not accurate for small sample size with large cure rate
    #cureThreshold <- c(rep(cure[1],n0), rep(cure[2],n1))
    cureThreshold <- c(rep(bsl_param[3],n0), rep(cureRate1,n1))
    Tm <- ifelse(ucure<cureThreshold, 10^10, Tm)
    # ensure the cure rate is almost exactly the specified proportion
    # Tm[sample(1:n0,round(cure[1]*n0))] <- 10^10
    # Tm[sample((n0+1):(n0+n1),round(cure[2]*n1))] <- 10^10
    #
  }

  #
  # remap to scheduled visits
  # if (!is.null(scheduleVisit)){
  #   #1st element is the schedule matrix
  #   #2nd element is the visit window
  #   Tm <- re_map(scheduleVisit[[1]], visit_window = scheduleVisit[[2]],time=Tm)$timeRemap
  # }
  #****************************
  #* Simulate drop-out Time
  #****************************
  if (missing(drop_param0)&missing(drop_param1)){
    #cat("Notes: Drop-outs are not considered in the simulation.")

  }else if (missing(drop_param0)){
    drop_param0 <- drop_param1
  }else if (missing(drop_param1)){
    drop_param1 <- drop_param0
  }
  # in extreme case, time is 0, add 0.1;
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
  if (is.null(enrollmentType)){
    #- create the CDF;
    ent_cdf0 <- function(t){stats::integrate(entry_pdf0,lower=0,upper=t)$value}
    gen_ent0 <- function(y){stats::uniroot(function(x){ent_cdf0(x)-y},
                                           interval = c(0,trial_param[2]),extendInt="yes")$root}
    ent_cdf1 <- function(t){stats::integrate(entry_pdf1,lower=0,upper=t)$value}
    gen_ent1 <- function(y){stats::uniroot(function(x){ent_cdf1(x)-y},
                                           interval = c(0,trial_param[2]),extendInt="yes")$root}

    #set.seed(seed+1)
    tu0_0 <- runif(n0)
    tu0_1 <- runif(n1)
    t0_0 <-as.vector(unlist(lapply(tu0_0, gen_ent0)))
    t0_1 <-as.vector(unlist(lapply(tu0_1, gen_ent1)))
    t0 <- c(t0_0,t0_1)
    ot <- t0+Tm
  }else{

    entry_pwu <- function(n,entryP){
      Tbreak <- entryP[[1]]
      erate <- entryP[[2]]
      Tbreakdiff <- c(Tbreak[1], diff(Tbreak))
      nsum <- floor(Tbreakdiff*erate)
      # how long to recruit one pt
      rt <- 1/erate
      for (i in 1:length(Tbreak)){
        if (i==1) rtc <- rep(rt[i],nsum[i])
        else rtc <- c(rtc, rep(rt[i],nsum[i]))
      }
      etime <- cumsum(rtc)
      val <- etime[1:n]
      return(val)
    }
    t0_01 <- entry_pwu(n0+n1, entryP)
    t0_012 <- sample(t0_01,n0+n1)
    t0_1 <- t0_012[1:n1]
    t0_0 <- t0_012[(n1+1):(n1+n0)]
    t0 <- c(t0_0,t0_1)
    ot <- t0+Tm
    t_p2 <- max(t0)
  }


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
  ))
  if (sum(final$aval<0)>0){
    # the enrollment has not completed when the target events are obtained
    #final <- final %>% dplyr::filter(aval>0)
    final <- final[final$aval>0,]
    # stop("Check: the aval has negative values.")
    cat("the enrollment has not compelted when reaching target events")
  }
  if (!missing(HR_data)) {
    HR_data1 <- HR_data[,!names(HR_data) %in% c("HR", "grp")]
    nf <- ncol(final)
    final <- data.frame(final, HR_data1)
    colnames(final)[-c(1:nf)] <- paste0("X",1:ncol(data.frame(HR_data1)))
  }

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

