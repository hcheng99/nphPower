##################
## evaluate the relationship between follow-up time and sample size
###################
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




