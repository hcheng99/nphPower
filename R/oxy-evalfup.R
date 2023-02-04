##################
## evaluate the relationship between follow-up time and sample size
###################
#' @title Visualization of the Relationship between Follow-up and Sample Size
#' @description \code{evalfup} function displays the graph showing the relationship
#' between the follow-up time and the total sample size/event number required to
#' achieve the the same power
#' @param object returned object by function \code{pwr2n.NPH}
#' @param lower.time a numeric value specifying the shortest duration time
#' @param upper.time a numeric value specifying the longest duration time
#' @param size an integer specifying the planned total sample size
#' @param increment a numeric value specifying an increment number used for creating
#' a sequence of duration times in plotting, Default: 0.5
#' @param xlabel a text for labeling the x axis in the plot, Default: 'Follow-up Time'
#' @param ylabel a text for labeling the y axis in the plot, Default: 'Total Sample Size'
#' @param title a text for title in the plot: 'Relationship between Follow-up and
#' Total Sample Size'
#' @return
#' a graph showing the relationship and a list of components:
#' \item{approx.time}{approximate follow-up time corresponding to specified sample size
#' to reach the same target power }
#' \item{original}{a list with elements of \code{x} and \code{y}. Vector \code{x} contains
#' the follow-up duration and vector \code{y} contains the corresponding sample
#' size}
#' \item{interp}{a list containing the interpolated \code{x} and \code{y} included
#' in \code{original}}
#' \item{Esize}{a vector of events number corresponding to \code{x} in \code{original}}
#' @details
#' The \code{evalfun} function helps to evaluate the relationship between
#' sample size/event number and follow-up duration. It retrieves the trial
#' design information from the \code{object} returned by \code{pwr2n.NPH}
#' function. A sequence of follow-up times starting from \code{lower.time}
#' and ending with \code{upper.time} are generated. The number of subjects and
#' number of events required for achieving the specified power in \code{object}
#' are calculated at each time point. An interpolation function \code{approx}
#' from \pkg{stats} is applied to smooth the curves. In case of
#' proportional hazards, the follow-up duration has little impact on the
#' event number except for variations from numeric approximations, while in
#' case of nonproportional hazards, the follow-up time imposes an important impact
#' on both the total sample size and event number.
#' @examples
#' # The following code takes more than 5 seconds to run.
#' \donttest{
#' # define design parameters
#'   t_enrl <- 12; t_fup <- 18; lmd0 <- log(2)/12
#' # define hazard ratio function
#'   f_hr_delay <- function(x){(x<=6)+(x>6)*0.75}
#' # define control hazard
#'   f_haz0 <- function(x){lmd0*x^0}
#' # perform sample size calculation using logrank test
#' # generate weight for test
#'   wlr <- gen.wgt(method="LR")
#'   snph1 <- pwr2n.NPH(entry = t_enrl, fup = t_fup, Wlist = wlr,
#'                     k = 100, ratio = 2, CtrlHaz = f_haz0, hazR = f_hr_delay)
#'
#' # suppose the follow-up duration that are taken into consideration ranges
#' # from 12 to 24. The planned number of patients to recruit 2200.
#' # draw the graph
#'   efun <- evalfup(snph1,lower.time = 12, upper.time = 24, size = 2200,
#'                title = NULL)
#' }
#' @rdname evalfup
#' @export
#' @importFrom stats approx
#' @importFrom graphics text points segments
evalfup <- function(object, lower.time, upper.time, size,
                    increment=0.5, xlabel = "Follow-up Time",
                    ylabel = "Total Sample Size/Event Number",
                    title = "Relationship between Follow-up and \n Total Sample Size"
){
  fupseq <- seq(from=lower.time, to=upper.time, by=increment)
  N <- c()
  D <- c()
  for (i in 1: length(fupseq)){
    #print(fupseq[i])
    tmp <- pwr2n.NPH(entry     = object$studytime[1]
                       ,fup      = fupseq[i]
                       ,k        = object$inputfun$k
                       ,ratio    = object$RandomizationRatio
                       ,Wlist  = object$inputfun$Weightfunctions
                       ,CtrlHaz=object$inputfun$controalhazard
                       ,transP1=object$inputfun$transP1
                       ,transP0=object$inputfun$transP0
                       ,hazR     = object$inputfun$hazardratio
                       ,alpha    = object$inputfun$alpha
                       ,beta     = object$inputfun$beta
                       ,entry_pdf0=object$inputfun$entrypdf0
                       ,entry_pdf1 =object$inputfun$entry_pdf1
                       ,summary=FALSE
    )
    N <-c(N,tmp$totalN)
    D <- c(D,tmp$eventN)
  }
  Nint <- ceiling(N/(object$RandomizationRatio+1))*(object$RandomizationRatio+1)
  Dint <- ceiling(D)
  ymax <- max(size,Nint)
  ymin <- min(size,Nint)
  if (max(Nint)<size){
    message("The largest required sample size ( ",max(Nint)," ) is less than the target  sampe size ", size, ". Consider reduce the follow-up time.")

  }else if (min(Nint) >size){
    message("The smallest required sample size ( ",min(Nint)," ) is greater than the target  sampe size ", size, ". Consider increase the follow-up time.")
  }
  flag <-  (max(Nint)<size)|(min(Nint) >size)
  interp <- stats::approx(fupseq,y=Nint,method="linear")
  interp_E <- stats::approx(fupseq,y=Dint,method="linear")
  if (flag != TRUE){

    y1 <- interp$y
    x1 <- interp$x
    p1 <- max(which(y1 >= size))
    p2 <- min(which(y1  <= size))
    xsize <- mean(c(x1[p1:p2]))
    Ey <- ceiling(mean(c(interp_E$y[p1:p2])))
    move <- (upper.time-lower.time)/7
  }else {
    xsize <- NA
  }

  plot(x=fupseq,y=Nint,ylim=c(min(Dint),ymax),cex=0.8,col="red",
       xlab= xlabel, ylab=ylabel,
       main = title)
  graphics::points(interp$x,interp$y,pch=16,cex=0.3)
  graphics::points(x=fupseq, y = Dint, cex = 0.8, col = "blue", pch = 0)
  graphics::points(x=interp_E$x, y = interp_E$y, cex = 0.25, col = "black",
                   pch = 3)
  if (flag !=TRUE){
    graphics::segments(x0=xsize,x1=xsize,y0=0,y1=size,lty=2)
    graphics::segments(x0=0,x1=xsize,y0=size,y1=size,lty=2)
    graphics::text(paste0("(",round(xsize,digits = 1),",",Ey,",",size,")"),x=xsize+move,
         y=size*1.05,cex=0.7)
  }

  legend("topright",legend=c("N","Interpolated N", "E","Interpolated E"),
         pch=c(1, 16, 0,3),
         col=c("red","black","blue","black"),cex=c(0.5), horiz = TRUE)
  # return values
  original <- list(x=fupseq,y=N)
  return(list(
    approx.time = xsize,
    original = original,
    interp = interp,
    Esize = D

  ))
}



