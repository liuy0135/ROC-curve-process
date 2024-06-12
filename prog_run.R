getwd()
setwd("/Users/lym/Desktop/nutstore/papers in revision/process of roc curve/revised sim code")
source("propEL_revised.R")
source("QZ1_revised.R")
source("QZ2_revised.R")
source("hhz_revised.R")
source("dat_cov_gen.R")
library(mvtnorm)
library(doParallel)
program_run=function(n0,n1,typ=1){
  library(mvtnorm)
  m=100
  B=500
  dat=dat_gen(n0,n1,typ)
  x=dat$x###x1
  y=dat$y###x2
  
  G=ecdf(x)
  from=1-G(quantile(y,0.995))+0.05
  to=1-G(quantile(y,0.005))-0.05
  seqq=seq(from,to,length.out=m)
  theta=find_true(seqq,typ)$theta
  
  #####proposed LD
  Sig1=gaussprocess(from, to, n1,n0,type=typ,m)
  apprx=rmvnorm(100000, mean = rep(0, m), Sig1)
  q95=quantile(apply(apprx^2, 1, max),0.95)
  q90=quantile(apply(apprx^2, 1, max),0.9)
  #####run
  tmp_prop=prop_el(x,y,seqq,theta,q90,q95)
  tmp_hhz=hhz_boot(x,y,seqq,theta,B)
  tmp_qz1=qz_quant1(x,y,seqq,theta,B)
  tmp_qz2=qz_quant2(x,y,seqq,theta,B)
  ####
  
  prop_val=tmp_prop$el_val
  prop90=tmp_prop$prop90
  prop95=tmp_prop$prop95
  
  hhz90=tmp_hhz$hhz90
  hhz95=tmp_hhz$hhz95
  
  qz1_90=tmp_qz1$c90
  qz1_95=tmp_qz1$c95
  
  qz2_90=tmp_qz2$c90
  qz2_95=tmp_qz2$c95
  
  cp90=c(prop90,hhz90,qz1_90,qz2_90)
  cp95=c(prop95,hhz95,qz1_95,qz2_95)
  
  return(list(cp=c(cp90,cp95),elval=prop_val))
}

multiResultClass <- function(cp1=NULL,cp2=NULL,cp3=NULL,val1=NULL,val2=NULL,val3=NULL)
{
  me <- list(
    cp1= cp1, cp2= cp2, cp3= cp3,
    val1 = val1, val2 = val2, val3 = val3
  )
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}


cores=detectCores()

cl <- makeCluster((cores-2))
registerDoParallel(cl)

asim=500

result1=foreach(icount(asim)) %dopar% {
  rec1=program_run(50,50,1)
  rec2=program_run(50,50,2)
  rec3=program_run(50,50,3)
  CP= multiResultClass()
  CP$cp1=rec1$cp
  CP$val1=rec1$elval
  CP$cp2=rec2$cp
  CP$val2=rec2$elval
  CP$cp3=rec3$cp
  CP$val3=rec3$elval
  return(CP)
}
stopCluster(cl)

ma1=ma2=ma3=matrix(0,asim,8)

for(k in 1:asim){
  ma1[k,]=result1[[k]]$cp1
  ma2[k,]=result1[[k]]$cp2
  ma3[k,]=result1[[k]]$cp3
}
M1=apply(ma1,2,mean)
M2=apply(ma2,2,mean)
M3=apply(ma3,2,mean)

MM=rbind(M1,M2,M3)

M=format(round(MM, 3), nsmall = 3)

library(xtable)
xtable(M)

