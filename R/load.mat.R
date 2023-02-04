load.mat <- function(dat){
  ## dat: 1-3rd columns  are time, status and treatment
  var1 <- dat[,1]
  var2 <- dat[,2]
  var3 <- dat[,3]
  nsize <- dim(dat)[1]
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package \"survival\" needed for this function to work. Please install it.",
         call. = FALSE)
  }


  fit <-  survival::survfit(survival::Surv(var1, var2) ~var3)
  dat0 <- data.frame(time = fit$time,
                     n.risk = fit$n.risk,
                     n.event = fit$n.event,
                     n.censor = fit$n.censor
  )
  ## get number of events
  en <- fit$strata
  dat1 <- dat0[1:en[1],]
  dat2 <- dat0[-c(1:en[1]),]
  #n.event.x refers to the 2nd strata
  # if the trt is 1,0, then 1 is treated as x
  datM  <- merge(x=dat2,y=dat1,by="time",all = TRUE)

  # datM <- datM %>% dplyr::mutate_at(datM,
  #   dplyr::vars(datM$n.event.x,datM$n.censor.x,datM$n.event.y,datM$n.censor.y),
  #                    ~replace(., is.na(.), 0))

  datM[,c("n.event.x","n.censor.x","n.event.y","n.censor.y")] <-
    apply(cbind(datM$n.event.x,datM$n.censor.x,datM$n.event.y,datM$n.censor.y),2,
          function(x){replace(x,is.na(x),0)})
  # get the risk at event for last time point
  riskly <- dat1[en[1],2]-dat1[en[1],3]-dat1[en[1],4]
  risklx <- dat2[en[2],2]-dat2[en[2],3]-dat2[en[2],4]
  if (is.na(datM[dim(datM)[1],"n.risk.x"])) {datM[dim(datM)[1],"n.risk.x"] <- risklx}
  if (is.na(datM[dim(datM)[1],"n.risk.y"])) {datM[dim(datM)[1],"n.risk.y"] <- riskly}
  datM$n.risk.x <- zoo::na.locf(datM$n.risk.x,fromLast = TRUE)
  rep_y <- zoo::na.locf(datM$n.risk.y,fromLast = TRUE)
  l <- length(rep_y)
  datM$n.risk.y[1:l] <- rep_y

  datM$n.risk.y[-c(1:l)] <-
    datM$n.risk.y[l]-
    (datM$n.event.y[-c(1:l)]+datM$n.censor.y[-c(1:l)])

  datM$event <- (datM$n.event.x+datM$n.event.y)
  datM$risk <- (datM$n.risk.x+datM$n.risk.y)
  datM$exp.x <- (datM$event*datM$n.risk.x/datM$risk)
  datM$V <- ifelse(datM$risk==1,0,
                datM$event*datM$n.risk.x/datM$risk*(datM$risk-datM$event)*
                  (datM$n.risk.y)/datM$risk/(datM$risk-1))

	return(datM)
}
