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
