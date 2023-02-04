
## calculate the Z statistic based on the weight
CalZ <- function(w,data){
    num_w1 <-t(w) %*%(data$n.event.x-data$exp.x)
    v_w1 <- t(w^2)%*%data$V
    z_w1 <- num_w1/sqrt(v_w1)
    return(z_w1)
}
## calculate the numerator based on the weight
CalN <- function(w,data){
  num_w1 <-t(w) %*%(data$n.event.x-data$exp.x)
  return(num_w1)
}
## conditional assignment
chkV <- function(x,y=x,chkV=0,assignV=1){
     ifelse(x==chkV,assignV,y)
}

########## the following two functions are for sample size calculation#########
## without staggered entry
CrtTM <- function(Plist,K=10){
  Plist[2] <- 1-exp(-Plist[2])
  Plist[5] <- 1-exp(-Plist[5])
  ylist <- unlist(lapply(Plist,function(x){1-(1-x)^{1/K}}) )
  tranM <- matrix(c(1,0,0,0,0,1,0,0,ylist[1:2],1-sum(ylist[1:3]),ylist[3],
                    ylist[4:6],1-sum(ylist[4:6])),ncol=4)
  return(tranM)
}
#*********************************************
# CrtTM_C ----
## function to create transition matrix with administered censoring
CrtTM_C0 <- function(Plist,K=10,a,b,i=i){
  # a: entry time
  # b: follow-up time
  # i: current unit
  # K: subinterval per unit time
  Plist[2] <- 1-exp(-Plist[2])
  Plist[5] <- 1-exp(-Plist[5])
  ylist <- unlist(lapply(Plist,function(x){1-(1-x)^{1/K}}) )
  tmp <- (a+b)*K-(i-1)

  tranM <- matrix(c(1,0,0,0,0,
                    0,1,0,0,0,
                    ylist[1:2],1-sum(ylist[1:3],1/tmp),ylist[3],1/tmp,
                    ylist[4:6],1-sum(ylist[4:6],1/tmp),1/tmp,
                    0,0,0,0,1), ncol=5)
  if (tmp==1){
    tranM[3:4,3:4] <- 0
    tranM[5,3] <- 1-sum(tranM[1:2,3])
    tranM[5,4] <- 1-sum(tranM[1:2,4])
  }
  return(tranM)
}

#*********************************************
# CrtTM_C ----
## function to create transition matrix with specified enrollment
CrtTM_C <- function(Plist,K=10,tott,epdf0,epdf1,i=i){
  # a: entry time
  # b: follow-up time
  # i: current unit
  # K: subinterval per unit time
  Plist[2] <- 1-exp(-Plist[2])
  Plist[5] <- 1-exp(-Plist[5])
  ylist <- unlist(lapply(Plist,function(x){1-(1-x)^{1/K}}) )
  Np <- tott*K
  epoint <- Np-(i-1)


  tmp0 <- stats::integrate( epdf0,(epoint-1)/Np,epoint/Np)$value/
    stats::integrate( epdf0,0/Np,epoint/Np)$value

  tmp1 <- stats::integrate( epdf1,(epoint-1)/Np,epoint/Np)$value/
    stats::integrate( epdf1,0/Np,epoint/Np)$value


  tranM <- matrix(c(1,0,0,0,0,
                    0,1,0,0,0,
                    ylist[1:2],1-sum(ylist[1:3],tmp0),ylist[3],tmp0,
                    ylist[4:6],1-sum(ylist[4:6],tmp1),tmp1,
                    0,0,0,0,1), ncol=5)

  if (sum(tranM[c(1,2,5),3])>=1){
    tranM[3:4,3] <- 0
    tranM[5,3] <- 1-sum(tranM[1:2,3])

  }
  if (sum(tranM[c(1,2,5),4])>=1){
    tranM[3:4,4] <- 0
    tranM[5,4] <- 1-sum(tranM[1:2,4])

  }

  if (sum(tranM<0) >0){
    print(i)
    print(tranM)

    stop("The transition matrix has negative values. Please check!")}

  return(tranM)
}
