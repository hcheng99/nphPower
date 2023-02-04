gen.wgt <- function(method=c("LR")
                    , param
                    , theta=0.5 )
{

  if (missing(method)){stop("method must be specified.")}
  if (sum(method=="FH")>1){stop("one FH method a time.")}

  lst <- list()
  nlen <- length(method)
  for (i in 1:nlen){
    if (!method[i] %in% c("LR","FH","Wilcoxon","Tarone","Maxcombo",
                          "Maxcross")){
      stop("method must be one of 'LR','FH','Wilcoxon','Tarone',
           'Maxcombo','Maxcross'")
    }
    if ("LR" %in% method[i]){
      lst <- c(lst,LR=function(x){x^0})
    }
    if ("FH" %in% method[i]){
      if (missing(param)){stop("param must be specified for FH method")}
      if (length(param)!=2){stop("both rho and gamma must be
                              specified for FH method")}
      lst <- c(lst, FH=function(x){(1-x)^param[1]*x^(param[2])})
    }
    if ("Wilcoxon" %in% method[i]){
      lst <- c(lst,Wil=function(x){1-x})
    }
    if ("Tarone" %in% method[i]){
      lst <- c(lst,tar= function(x){(1-x)^0.5})
    }
    if ("Maxcombo" %in% method[i]){
      lst <- c(lst,list(
        conf=function(x){x^0},
        f01=function(x){x},
        f10=function(x){1-x},
        f11=function(x){(1-x)*x}
      ))
    }
    if ("Maxcross" %in% method[i]){
      if (missing(theta)){stop(
        "theta must be specified for Maxcross test."
      )
      }
      lst <- c(lst,list(
        conf=function(x){x^0},
        f01=function(x){x},
        f10=function(x){1-x},
        fcross=function(x,pp=theta){
          (x<=pp)*(1/pp*x-1)+(x>pp)*(1/(1-pp)*(x-pp))}
      ))

    }



  }

  return(lst)

}
