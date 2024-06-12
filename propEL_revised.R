Gt <- function(z){
  0 * I(z<=-1) + 1* I(z>=1) +(15/16)*(z-(2/3)*z^3+(1/5)*z^5+(8/15))*I(abs(z)<1)
}
####proposed method
fnone <- function(x, y, t0,theta0){
  n0 <- length(x)
  n1 <- length(y)
  h0=n0^(-1/3);h1=n1^(-1/3)
  ## parameter const
  lo <- max(min(x), min(y))
  up <- min(max(x), max(y))
  
  ###########object function
  fn <- function(beta){ ## beta <- c(eta, theta)
    
    Gn <- mean(Gt((beta - x) / h0))
    Fn <- mean(Gt((beta - y) / h1))
    
    part1 <- 2 * n1 * (  Fn * log( Fn / (1 - theta0) )   +  (1 - Fn)  * log( (1 - Fn) / theta0 ) ) #(2.11)
    part2 <- 2 * n0 * (  Gn * log( Gn / (1 - t0    ) )   +  (1 - Gn)  * log( (1 - Gn) / t0     ) ) #(2.11)
    a <- part1 + part2
    return(a)
  }
  
  eta.int <- quantile(x, probs = 1- t0)
  ######################################
  if(eta.int < lo[1]){
    eta.int <- lo[1] + 0.01
  }else if(eta.int > up[1]){
    eta.int <- up[1] - 0.01
  }else{
    eta.int <- eta.int
  }
  #############???Ż?????  ######
  ret <-optim(par=eta.int, fn=fn, lower=lo, upper=up,  method = "L-BFGS-B")
  #BB::BBoptim(par=eta.int, fn=fn, lower=lo, upper=up, control = list(maximize=F, trace=FALSE, maxit=1e4))
  bb <- ret$value  #fn??ֵ
  return(list(beta=ret$par, `-2LLR` = bb , Pval = 1 - pchisq(bb, df = 1)))
  
}
 prop_el=function(x,y,seqq,theta,q90,q95){
   m=length(seqq)
   el_seq=NULL
   for(k in 1:m){
     el_seq[k]=fnone(x, y,seqq[k],theta[k])$`-2LLR`
   }
   prop_val=max(el_seq)
   prop90=ifelse(prop_val<q90,1,0)
   prop95=ifelse(prop_val<q95,1,0)
   return(list(prop90=prop90,prop95=prop95,el_val=prop_val))
   
 }
