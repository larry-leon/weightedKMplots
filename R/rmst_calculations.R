# Mean for general S (truncated at L)
# \int_{0}^{L}S(t)dt
# fatx=S, x=t
LS.int<-function(fatx,x,L=Inf){
  f <- fatx[1:(length(fatx)-1)]
  d.x <- diff(x)
  int.f <- x[1]+sum((d.x*f)[which(x[-1]<=L)])
  return(int.f)
}


mu.L<-function(L,S,t){
  mu <- LS.int(fatx=S,x=t,L=L)
  return(mu)
}

# Need to check!
LS.int.del<-function(fx,x,dx,L=Inf){
  int.f < -sum((dx*fx)[which(x<=L)])
  return(int.f)
}


#temp<-NA.CHR.Weighted(time=time0,delta=(delta0==1))

# S0.KM<-temp$S.KM
# se0.KM<-temp$se.KM
# t0.points<-temp$at.points
# R0<-temp$n.risk
# dN0<-temp$n.events
# # May need modification (R0-1 term) for ties
# dC0<-ifelse(R0>1,dN0/(R0*(R0-1)),0.0)
###################################
# Restricted mean (point estimate)
###################################
# m0.L<-LS.int(S0.KM,t0.points,L)
# m0.t<-unlist(lapply(as.list(t0.points),mu.L,S=S0.KM,t=t0.points))
# # Gill notation
# mbar.0<-c(m0.L-m0.t)
# var.m0L<-LS.int.del(fx=mbar.0^2,dx=dC0,x=t0.points,L=L)

