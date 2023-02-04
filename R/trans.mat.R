trans.mat <- function(...){
 fup <- numN <-k <- transP1 <- haz_val <- transP0 <- tot_time<-NULL
 ctrlRate <- entry <- entry_pdf0 <- entry_pdf1 <- ratio <- hazR<- NULL
 getlist <- list(...)
 listname <- names(getlist)
 for (i in 1:length(listname)){
    assign(listname[i],getlist[[i]])
  }

  ##initial states for treatment
  D1 <- c(0,0,1,0,0)
  ##initial states for control
  D0 <- c(0,0,0,1,0)
  temp <- diag(rep(1,4))
  pdat <- matrix(NA,nrow = round(numN,digits = 0),ncol=21)
  ti <-0
  l <- 1
  L_trans <- list()
  ctrl_prope <- c()
  ctrl_prop0 <- c() ## proportion of control in temp0
  ctrl_prop1 <- c() ## proportion of treatment in temp0
  trt_stay <- c()    ## probability of staying at treatment
  trt_evt <- c() ## probability of event for treatment
  #interm <- list()
  #matrix(0,nrow = round(numN,digits = 0),
  #  ncol = round(numN,digits = 0))

  for (i in 1:numN){
    #j <- ceiling(i/k)
    j <- i/k
    trans.prob <- c(transP1[1],haz_val[i],transP1[2],transP0[1],ctrlRate[i],
                    transP0[2])
    # before a, there is no censor
      if (j <=fup){
        L_trans[[i]] <- cbind(rbind(CrtTM(Plist = trans.prob,K=k),rep(0,4)),
          c(0,0,0,0,1))
      }else if (j > fup){
        L_trans[[i]] <- CrtTM_C(Plist = trans.prob,tott=tot_time,
          epdf0=entry_pdf0,epdf1=entry_pdf1,K=k,i=i)
      }

    if (i==1){temp1 <- L_trans[[i]]%*%D1 ; temp0 <- L_trans[[i]]%*%D0
    } else {temp1 <- L_trans[[i]]%*%temp1; temp0 <- L_trans[[i]]%*%temp0;}

     ### to check whether there is drop-in rate (a constant rate)
    if (transP0[2]>0){
      trt_stay[i] <- L_trans[[i]][3,3]
      ctrl_prop0[i] <- temp0[4]
      ctrl_prop1[i] <- temp0[3]
      ctrl_prope[i] <- temp0[2]
      trt_evt[i] <- L_trans[[i]][2,3]
      p_t <- L_trans[[i]][3,4]
      tmp_in <- c()
      ## if NPH exists, there is a lag effect for drop-in effects
      if (i==2){
        tmp_in <- D0[4]*p_t
        e_tmp <- ctrl_prope[i-1]+D0[4]*p_t*trt_evt[1]+ctrl_prop0[i-1]*L_trans[[i]][2,4]
      } else if (i>2){
        for (j in 1:(i-1)){
          p_t <- L_trans[[i-j]][3,4]
          if (i-(j+1)<(i-1)&(i-j-1)>0){
            tseq <- trt_stay[(i-(j+1)):(i-1)]
            cprop <- ctrl_prop0[i-(j+1)]
          }else{
            tseq <- 1
            cprop <- D0[4]
          }
          tmp_in[j] <- p_t*cprop*prod(tseq, na.rm = TRUE)*trt_evt[j]

          itmp <- sum(tmp_in, na.rm = TRUE)+
            ctrl_prop0[i-1]*L_trans[[i]][2,4]
        }
        e_tmp <- ctrl_prope[i-1]+itmp
      }
      tmp_in <- NULL

      if (i>1){
        temp0[2] <- e_tmp;
        if (i <numN){
         temp0[3] <- 1-sum(temp0[-3])
        }else if (i==numN){
         temp0[5] <- 1-sum(temp0[c(1,2)])
        }

        ctrl_prop0[i] <- temp0[4]
        ctrl_prop1[i] <- temp0[3]
        ctrl_prope[i] <- temp0[2]
      }
    }


     #if(numN-i<4){print(i); print(L_trans[[i]])}
    pdat[i,1:11] <- c(ti,t(temp1),t(temp0))
    #print(c(j,fup, ti,t(temp1),t(temp0)))

    if (i==1) {
      #phi
      pdat[i,12] <- 1
      ## hazard of dying for trt=1
      pdat[i,13] <- pdat[i,3]
      # for trt=0
      pdat[i,14] <- pdat[i,8]

    }
    else {
      #phi
      pdat[i,12] <- sum(pdat[i-1,5],pdat[i-1,4])/sum(pdat[i-1,10],pdat[i-1,9])
      #hazard1
      pdat[i,13] <- (pdat[i,3]-pdat[i-1,3])/(pdat[i-1,5]+pdat[i-1,4])
      #hazard 0
      pdat[i,14] <- (pdat[i,8]-pdat[i-1,8])/(pdat[i-1,9]+pdat[i-1,10])
    }
    #theta
    ## in case the initial values are zero or Inf

    if (i==1){
      if (hazR(0)==0|hazR(0)==Inf) {
        pdat[i,15] <- hazR(0.0000001)
      }else{pdat[i,15] <- hazR(ti)}
    }
    else{ pdat[i,15] <- pdat[i,13]/pdat[i,14]}
    ti <- ti+1/k
    # pdat[i,13] <- 2
    #gamma
    pdat[i,16] <- pdat[i,12]*pdat[i,15]*ratio/(1+pdat[i,12]*pdat[i,15]*ratio)-
      pdat[i,12]*ratio/(1+pdat[i,12]*ratio)
    #eta
    pdat[i,17] <- pdat[i,12]*ratio/(1+pdat[i,12]*ratio)^2
    # print(pdat[i,])
    }
  #save(interm,file="interm.Rdata")
  #pdat
  rho <- (pdat[numN,3]*ratio+pdat[numN,8])/(ratio+1)

  ## sample size
  for (i in 1:numN){
    ## S1
    pdat[i,20] <- prod(1-pdat[1:i,13])
    ## S0
    pdat[i,21] <- prod(1-pdat[1:i,14])
    ## mean survival
    pdat[i,19] <- stats::weighted.mean(c(pdat[i,20],pdat[i,21]),w=c(ratio,1))
    ## rho
    if (i==1){
      pdat[i,18] <- (pdat[i,3]*ratio+pdat[i,8])/rho/(ratio+1)
    }else
    {
      pdat[i,18] <- ((pdat[i,3]-pdat[i-1,3])*ratio+pdat[i,8]-pdat[i-1,8])/rho/(ratio+1)

    }
  }

  pdat <- as.data.frame(pdat)
  names(pdat) <-
    c("ti","E_L","E_E","E_Ae","E_Ac","E_C","C_L","C_E","C_Ae",
      "C_Ac","C_C","phi","hazard_E","hazard_C","theta","gamma","eta",
      "rho","S","S1","S0")
  ## numNber of weight functions
  pdat$eprob <- apply(cbind(pdat$C_E,pdat$E_E),1,
                      function(x) {stats::weighted.mean(c(x[1],x[2]),
                                                        w=c(1,ratio))})
  return(list(pdat=pdat,L_trans=L_trans))
}
