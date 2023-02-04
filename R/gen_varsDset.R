gen_varsDset <- function( nsub
                        , distlist
                        , parameterlist
                        , betas ){
  ndist <- length(distlist)
  nbetas <- length(betas)
  if ((ndist+1)!=nbetas){
    stop("betas must have one more element than distlist")
  }
  tmpD <- list()
  for (i in 1:ndist){
    if (distlist[i] == "n"){
      tmpD[[i]] <- stats::rnorm(nsub, parameterlist[[i]][1], parameterlist[[i]][2])
    } else if (distlist[i] == "b"){
      tmpD[[i]] <- stats::rbinom(nsub, 1, parameterlist[[i]][1])
    }
  }
  tmpD1 <- tmpD %>% unlist %>%
    matrix(nrow = nsub, ncol= ndist)
  lgHr <- betas[1] + tmpD1 %*% betas[-1]
  Hr <- exp(lgHr)
  out <- as.data.frame(tmpD1) %>% mutate(lgHr, Hr)
  return(out)
}

nsub <- 100
distlist <- list("n", "b")
parameterlist <- list( c(50,10), c(0.35))
betas <- c(1, -0.03, 1.5)

gen_varsDset(nsub, distlist, parameterlist, betas)
