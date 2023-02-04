#' @title Sample Size Calculation under Proportional Hazards
#' @description \code{pwr2n.LR} calculates the total number of events and total
#' number of subjects required given the provided design parameters based on either
#' schoenfeld or freedman formula.
#' @param method calculation formula, Default: c("schoenfeld", "freedman")
#' @param lambda0 hazard rate for the control group
#' @param lambda1 hazard rate for the treatment group
#' @param ratio randomization ratio between treatment and control. For example,
#'  ratio=2 if randomization ratio is 2:1 to treatment and control group.
#'  Default:1
#' @param entry enrollment time. A constant enrollment rate is assumed,
#'  Default: 0
#' @param fup follow-up time.
#' @param alpha type I error rate, Default: 0.05
#' @param beta type II error rate. For example,if the target power is 80%, beta is 0.2.
#'  Default: 0.1
#' @param alternative a value must be one of ("two.sided", "one.sided"), indicating
#' whether a two-sided or one-sided test to use.  Default: c("two.sided")
#' @param Lparam a vector of shape and scale parameters for the drop-out Weibull distribution,
#'  See Details below. Default: NULL
#' @param summary a logical controlling whether a brief summary is printed or not
#' , Default: TRUE
#' @return
#' a list of components including
#' \item{eventN}{a numeric value giving the total number of events}
#' \item{totalN}{a numeric value giving the total number of subjects}
#' \item{summary}{a list containing the input parameters and output results}
#' @aliases phsize
#' @details
#'  Both Schoenfeld's formula and Freedman's formula are included in the
#'  function \code{pwr2n.LR}.
#'  The total event number is determined by \eqn{\alpha, \beta} and
#'  hazard ratio, i.e., \eqn{\lambda_1/\lambda_0}. Other design parameters such as
#'  enrollment period affects the event probability and thus the total sample size.
#'  A fixed duration design is assumed in the calculation. All patients are enrolled
#'  at a constant rate within \code{entry} time and have at least \code{fup}
#'  time of follow-up. So the total study duration is \code{entry}+\code{fup}.
#'  If drop-out is expected, a Weibull distribution with shape parameter -\eqn{\alpha}
#'  and scale parameter - \eqn{\beta} is considered. The CDF of Weibull is
#'  \eqn{F(x)=1-exp(-(x/\beta)^\alpha)}, where \eqn{\alpha} is the shape
#'  parameter and \eqn{\beta} is the scale parameter. The event rate
#'  is calculated through numeric integration. See more details in
#'  \code{\link{cal_event}}.
#'
#' @references
#' Schoenfeld, D. (1981) The asymptotic properties of nonparametric
#'  tests for comparing survival distributions. Biometrika, 68,
#' 316–319.
#'
#' Freedman, L. S. (1982) Tables of the number of patients required
#' in clinical trials using the logrank test. Statistics in medicine, 1, 121–129.

#' @examples
#' # define design parameters
#' l0 <- log(2)/14; HR <- 0.8; RR <- 2; entry <- 12; fup <- 12;
#' eg1 <- pwr2n.LR( method    = c("schoenfeld")
#'                  ,l0
#'                  ,l0*HR
#'                  ,ratio=RR
#'                  ,entry
#'                  ,fup
#'                  ,alpha     = 0.05
#'                  ,beta      = 0.1
#' )
#' # event number, total subjects, event probability
#' c(eg1$eventN,eg1$totalN,eg1$eventN/eg1$totalN)
#'
#' # example 2: drop-out from an exponential with median time is 30
#' eg2 <- pwr2n.LR( method    = c("schoenfeld")
#'                  ,l0
#'                  ,l0*HR
#'                  ,ratio=RR
#'                  ,entry
#'                  ,fup
#'                  ,alpha     = 0.05
#'                  ,beta      = 0.1
#'                  ,Lparam = c(1,30/log(2))
#' )
#' # event number, total subjects, event probability
#' c(eg2$eventN,eg2$totalN,eg2$eventN/eg2$totalN)
#' @seealso
#'  \code{\link{pwr2n.NPH}},
#'  \code{\link{evalfup}}, \code{\link{cal_event}}
#' @rdname pwr2n.LR
#' @export
#' @importFrom stats integrate weighted.mean qnorm
pwr2n.LR <- function( method    = c("schoenfeld","freedman")
                          ,lambda0
                          ,lambda1
                          ,ratio  =  1
                          ,entry  = 0
                          ,fup
                          ,alpha = 0.05
                          ,beta  = 0.1
                          ,alternative = c("two.sided")
                          ,Lparam=NULL
                          ,summary = TRUE
)
  # the calculation assumes the exponential distribution for control and
  # treatment group
  # ratio: randomization ratio
  #
{
  if (missing(method)) {stop("method must be specified")}
  if (!(match.arg(method)  %in% (c("schoenfeld","freedman")))){
    stop("Error: method must be schoenfeld or freedman")
  }
  if (!is.null(Lparam)){
    l_shape=Lparam[1]; l_scale=Lparam[2]
    calp <- function(lambda1){
      f1 <- function(x){lambda1*exp(-x*lambda1-(x/l_scale)^l_shape)}
      F1 <- Vectorize(function(u){integrate(f1,lower=0,upper=u)$value})
      if (entry==0){
        e1 <-F1(fup)
      }else {
        e1 <- stats::integrate(F1,lower=fup,upper=entry+fup)$value/entry
      }

      return(e1)
    }
    }
  HR <- lambda1/lambda0
  if (entry==0){
    if (is.null(Lparam)){
      e0 <- 1-exp(-lambda0*fup)
      e1 <- 1-exp(-lambda1*fup)
    }else {
      e0 <- calp(lambda0)
      e1 <- calp(lambda1)
    }


  }
  else {
    if (is.null(Lparam)){
      intef <- function(x,l){(1-exp(-l*x))}
      e0 <-stats::integrate(function(x){intef(x,l=lambda0)},lower=fup,upper=entry+fup)
      e1 <-stats::integrate(function(x){intef(x,l=lambda1)},lower=fup,upper=entry+fup)
      e0 <- e0$value/entry
      e1 <- e1$value/entry
    }else {
      e0 <- calp(lambda0)
      e1 <- calp(lambda1)
    }



  }
  ## calculate the event rate
  erate <- stats::weighted.mean(c(e0,e1),
                         w=c(1/(1+ratio),ratio/(1+ratio)))
  ## calculate the number of events
 if (alternative == "two.sided"){
   numerator=(stats::qnorm(1-alpha/2)+stats::qnorm(1-beta))^2
 } else if (alternative == "one.sided"){
   numerator=(stats::qnorm(1-alpha)+stats::qnorm(1-beta))^2
 } else {
   stop ("alternative must be either 'two-sided' or 'one-sided'!")
 }
  if (method == "schoenfeld"){
    Dnum=numerator*(1+ratio)^2/ratio/log(HR)^2

  }else if (method=="freedman"){
    Dnum=numerator*(1+HR*ratio)^2/ratio/(1-HR)^2
  }

  N=Dnum /erate
  inparam <- c("Method", "Lambda1/Lambda0/HR","Entry Time", "Follow-up Time",
               "Allocation Ratio", "Type I Error", "Type II Error",
               "Alternative","Drop-out Parameter")
  if (is.null(Lparam)) {Lparam <- "Not Provided"
  }else {Lparam <- round(Lparam, digits = 3)}

  inval <- c(method, paste0(round(lambda1,digits=3),"/",round(lambda0,digits=3), "/",
                            round(lambda1/lambda0,digits=3)),
             entry, fup,ratio, alpha, beta,
             alternative,paste0(Lparam,collapse = ","))
  inputdata <- data.frame(parameter=inparam, value=inval)
  outparam <- c("Number of Events", "Number of Total Sampe Size",
                "Overall Event Rate")
  outval <- round(c(Dnum, N, Dnum/N),digits=3)
  outputdata <- data.frame(parameter=outparam, value=outval)
  summaryout <- list(input=inputdata,output=outputdata)
  if(summary ==TRUE){
    cat("------------------------------------------ \n ")
    cat("-----Summary of the Input Parameters----- \n")
    cat("------------------------------------------ \n ")
    names(inputdata) <- c("__Parameter__", "__Value__")
    print(inputdata, row.names = FALSE, right=FALSE)
    cat("------------------------------------------ \n ")
    cat("-----Summary of the Output Parameters----- \n ")
    cat("------------------------------------------ \n ")
    names(outputdata) <- c("__Parameter__", "__Value__")
    print(outputdata, row.names = FALSE, right=FALSE)

  }
  return(list(eventN=Dnum,totalN=N,summary=summaryout))
}

#' @title Event Rate Calculation
#' @description Calculate the event rate given the hazards and drop-out
#' distribution parameters
#' @param ratio allocation ratio
#' @param lambda1 hazard rate for treatment group
#' @param lambda0 hazard rate for control group
#' @param entry enrollment period time
#' @param fup follow-up period time
#' @param l_shape shape parameter of weibull distribution for drop-out
#' @param l_scale scale parameter of weibull distribution for drop-out
#' @return
#' a list of components:
#' \item{ep1}{event rate for treatment group}
#' \item{ep0}{event rate for control group}
#' \item{ep}{mean event rate weighted by the randomization ratio}
#' @details
#' The event rate is calculated based on the following assumptions: 1)
#' patients are uniformly enrolled within \code{entry} time; 2) survival
#' times for treatment and control are from exponential distribution; 3)
#' the drop-out times for treatment and control follow the weibull distribution.
#' The final rate is obtained via numeric integration:
#'
#' \deqn{
#'P=\int_{t_{fup}}^{t_{enrl}+t_{fup}} \Big \{
#' \int_0^{t}r(u)exp\big [-\int_0^u[r(x)+l(x)]dx \big]d(u) \Big \}
#' \frac{1}{t_{enrl}} dt
#' }
#' where \eqn{r(x)} is the hazard of event and \eqn{l(x)} is the hazard
#' of drop-out; \eqn{t_{enrl}} is the entry time and \eqn{t_{fup}} is the
#' follow-up duration.
#' @examples
#' # median survival time for treatment and control: 16 months vs 12  months
#' # entry time: 12 months ; follow-up time: 18 months
#' # the shape parameter for weibull drop-out : 0.5
#' # median time for drop-out : 48 =>
#' # scale parameter: 48/log(2)^(1/0.5)=100
#'   RR <- 1; l1 <- log(2)/16; l0 <- log(2)/12
#'   t_enrl <- 12; t_fup <- 18
#'
#'   cal_event(1,l1,l0,t_enrl,t_fup,0.5,100)
#' @rdname cal_event
#' @export
cal_event <- function(
  ratio
  ,lambda1
  ,lambda0
  ,entry
  ,fup
  ,l_shape
  ,l_scale
){
  calp <- function(lambda1){
    f1 <- function(x){lambda1*exp(-x*lambda1-(x/l_scale)^l_shape)}
    F1 <- Vectorize(function(u){integrate(f1,lower=0,upper=u)$value})
    e1 <- integrate(F1,lower=fup,upper=entry+fup)$value/entry
    return(e1)
  }

  ## treatment group
  e1 <- calp(lambda1)
  ## placebo
  e0 <- calp(lambda0)
  e <- ratio/(1+ratio)*e1+e0/(1+ratio)
  return(list(ep1 = e1,
              ep0 = e0,
              ep = e))
}

