findul <- function(step = 0.01, initStep = 0, fun, MLE, level = 3.84146, ...)
{
  value <- 0
  step1 <- step
  
  Lbeta <- MLE - initStep
  while (value < level&&Lbeta>0.002) {
    Lbeta <- Lbeta - step1
    value <-  tryCatch(fun(Lbeta, ...)$"-2LLR", error=function(err) 0)
  }
  Lbeta0 <- Lbeta
  Lbeta1 <- Lbeta + step1
  tempfun <- function(beta){
    return( level - fun(beta, ...)$"-2LLR" )
  }
  if(Lbeta>0){
    temp1 <- tryCatch(uniroot(tempfun, lower=Lbeta0, upper=Lbeta1), error=function(err) list(root=(Lbeta0+Lbeta1)/2))
    Lbeta <- temp1$root
    value1 <- tryCatch(level - temp1$f.root, error=function(err) 0)}
  else{Lbeta=0}
  value <- 0
  Ubeta <- MLE + initStep
  while (value < level&&Ubeta<0.998) {
    Ubeta <- Ubeta + step
    value <- tryCatch(fun(Ubeta, ...)$"-2LLR", error=function(err) 1)
  }
  Ubeta0 <- Ubeta
  Ubeta1 <- Ubeta - step
  if(Ubeta<1){
    temp2 <- tryCatch(uniroot(tempfun, lower=Ubeta1, upper=Ubeta0), error=function(err) list(root=(Ubeta0+Ubeta1)/2))
    Ubeta <- temp2$root
    value <- tryCatch(level - temp2$f.root, error=function(err) 0)}
  else{Ubeta=1}
  return(list(Low = Lbeta, Up = Ubeta))
}
