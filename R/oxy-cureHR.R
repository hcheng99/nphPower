#' @title Control hazard and hazard ratio generation function
#' @description Generate control hazard and hazard ratio function used for sample size calculation for cure model
#' @param pi0 cure rate for the control group
#' @param pi1 cure rate for the treatment group, Default: NULL
#' @param k0 shape parameter of the Weibull distribution for the control group
#' @param lmd0 rate parameter of the Weibull distribution for the control group
#' @param theta hazard ratio function
#' @param HRType hazard ratio function type. \code{susceptible} indicates the hazard
#' ratio function applies to the susceptible only; \code{overall} indicates the hazard
#' ratio function applies to the overall population; \code{delayed} indicates a cure model
#' with delayed treatment effects. See details.
#' @param tchg delayed timepoint for \code{HRType} = \code{delayed}, Default: NULL
#' @return
#' a list of components including
#' \item{ctrl_hr}{a hazard function for the control group}
#' \item{hr}{a hazard ratio function}
#' @details DETAILS
#' The control group has a survival function of \eqn{S_o0=\pi_0+(1-\pi_0)S_0},where
#' \eqn{\pi_0} is the cure rate and \eqn{S_0} is the survival function for the susceptible
#' population. For \code{HRType} = \code{susceptible}, the user also needs to provide
#' the cure rate for the experimental group. The provided hazard ratio applies to the susceptible
#' population only. The returned hazard ratio function is the overall one. For \code{HRType}=\code{delayed},
#' the returned hazard ratio is derived based on the paper of Wei and Wu (2020) .
#' @references
#' Wei, J. and Wu, J., 2020. Cancer immunotherapy trial design with
#' cure rate and delayed treatment effect. Statistics in medicine, 39(6), pp.698-708.
#' @examples
#' p0 <- 0.2; p1 <- 0.3; param <- c(1, log(2)/12);
#' theta_eg <-function(t){t^0*0.7}
#' fit <- cureHR(p0, p1, param[1], param[2],theta_eg, HRType="susceptible")
#' # with delayed effects
#' theta_eg2 <- function(t){(t<=9)+(t>9)*0.7}
#' fit2 <- cureHR(p0, p1, param[1], param[2],theta_eg2, HRType="delayed", tchg=9)
#' @seealso
#'  \code{\link[stats]{integrate}}
#' @rdname cureHR
#' @export
#' @importFrom stats integrate
cureHR<- function(pi0, pi1=NULL, k0, lmd0, theta, HRType,tchg=NULL){
  #model 1: PHCM
  #model 2: PHCRM
  lmd0 <-lmd0^k0
  # functions for susceptible patients
  s0 <- function ( x ) {exp ( - lmd0 *x ^ k0 ) }
  f0 <- function(x){exp ( - lmd0 *x ^ k0 )*lmd0*k0*x^(k0-1)}
  h0 <- function ( x ) { lmd0 *k0*x ^( k0 -1) }
  # control hazard function for all patients
  o0h <- function ( x ) {(1 - pi0 )*s0 ( x )*h0 ( x )/( pi0 +(1 - pi0 )*s0 ( x ) ) }
  if (HRType=="susceptible"){
    # survival function for susceptible patients from experimental arm
    ss1 <- function(x){exp(-as.numeric(stats::integrate(function(s){h0(s)*theta(s)},0,x)$value))}
    s1 <- Vectorize(ss1)
    f1 <- function(x){s1(x)*theta(x)*h0(x)}
    # hazard function susceptible patients from experimental arm
    o1h <- function(x){(1-pi1)*f1(x)/(pi1+(1-pi1)*s1(x))}

    # define the overall hazard ratio function
    HRf <- function(x){
      num = o1h(x)
      den = o0h(x)
      ifelse(den>0, num/den,Inf)
    }
    hr <- HRf
  }else if (HRType=="overall"){
    hr <- theta
  }else if (HRType=="delayed"){
    # due to the definition of s, no need to add s_0^1-hr.
    St0 <-  s0(tchg)
    pi1.tilde=1/(1+((pi0+(1-pi0)*St0)/pi1-1)/St0)
    c=(pi0+(1-pi0)*St0)/(pi1.tilde+(1-pi1.tilde)*St0)
    c.tilde=c*(1-pi1.tilde)/(1-c*pi1.tilde)
    # survival function for susceptible patients from experimental arm
    ss1 <- function(x){exp(-as.numeric(stats::integrate(function(s){h0(s)*theta(s)},0,x)$value))}
    s1 <- Vectorize(ss1)
    f1 <- function(x){s1(x)*theta(x)*h0(x)}
    o1h <- function(x){((1-pi1)*f1(x)*c.tilde/(pi1+(1-pi1)*s1(x)*c.tilde))*(x>tchg)+
        ((1-pi0)*f1(x)/(pi0+(1-pi0)*s1(x)))*(x<=tchg)
    }
    # define the overall hazard ratio function
    HRf <- function(x){
      num = o1h(x)
      den = o0h(x)
      ifelse(den>0, num/den,Inf)
    }
    hr <- HRf
  }

  return(list(ctrl_hr = o0h,
              hr= hr)
  )

}
