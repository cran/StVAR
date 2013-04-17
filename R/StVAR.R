StVAR <-
function(Data,lag=1,v=1,maxiter=1000,meth="BFGS",hes="FALSE",init="na",tre="C")
{

if (lag-round(lag)!=0) return("lag must be an integer")
if (lag<0) return("lag number must be positive")

Z<-t(embed(Data,lag+1)) ; T<-ncol(Z)+lag ; l<-ncol(Data)  ## No. of variables in VAR

if (ncol(cbind(Data))<2) return("no. of variables must be greater than 1")
if (v<0) return("degrees of freedom must be greater or equal to 0")
if (maxiter<10) return("Iteration must be at least 10")
if (2*T<3*((lag+3)*(l+1)*l/2)+3*l) return(list(((lag+3)*(l+1)*l/2)+3*l,"Too many parameters for given sample size. Reduce the number of lags."))

##Likelihood Function
L<-function(a)
{
 S<-BlockTop(a[((l*(l+1)/2)*(1-1)+1):((lag+3)*(l+1)*l/2)],l,lag)$S ##Var-Cov Matrix
 F<-solve(S)/v ## The Inverse

 M0<-rep(c(a[(((lag+3)*(l+1)*l/2)+1):(((lag+3)*(l+1)*l/2)+l)]),(lag+1)) 
 M1<-c(a[(((lag+3)*(l+1)*l/2)+l+1):(((lag+3)*(l+1)*l/2)+2*l)]) ##Linear
 M2<-c(a[(((lag+3)*(l+1)*l/2)+2*l+1):(((lag+3)*(l+1)*l/2)+3*l)]) ##Quadratic 

##Linear and Quadratic trend
 trenz<-trend(M1,M2,T,lag,l) 
 m1<-trenz$m1 ; m2<-trenz$m2 
 
 if(tre=="Q") M<-(M0+c(m1[,,])+c(m2[,,]))
 if(tre=="L") M<-(M0+c(m1[,,]))
 if(tre=="C") M<-(M0)

 D <- 1 + diag(t(Z-M)%*%F%*%(Z-M)) ##Quadratic form

## Likelihood function
 LLn <- (T-lag)*const - 0.5*(T-lag)*log(det(S)) - 0.5*(v+(lag+1)*l)*sum(log(D))
 neg.LLn <- -LLn
 neg.LLn
}

if(init[1]=="na" & tre=="Q") int<-c(runif(((lag+3)*(l+1)*l/2)+3*l,0,0.1))  ##Initialization
if(init[1]=="na" & tre=="L") int<-c(runif(((lag+3)*(l+1)*l/2)+2*l,0,0.1))  ##Initialization
if(init[1]=="na" & tre=="C") int<-c(runif((((lag+3)*(l+1)*l/2)+l),0,0.1))  ##Initialization

if(init[1]!="na") int<-init
const<-log(gamma((v+(lag+1)*l)/2))-log(gamma(v/2))-0.5*l*(lag+1)*log(pi*v)

op<-optim(int,L,hessian=hes,control=list(trace=1,maxit=maxiter,reltol=1e-12),method=meth) 
a<-op$par ; Like<-op$value 

###################################Parameters##################

S<-BlockTop(a[((l*(l+1)/2)*(1-1)+1):((lag+3)*(l+1)*l/2)],l,lag)$S
PP<-Par(a,l,lag,v,T,tre)
M0<-PP$M0 ; M1<-PP$M1 ; M2<-PP$M2 
Delta0<-PP$Delta0 ; Delta1<-PP$Delta1 ; Delta2<-PP$Delta2
B1<-PP$B1 ; s2<-PP$s2 ; Q<-PP$Q
m1<-PP$m1 ; m2<-PP$m2


#############Fitted values/Residuals/Con. Covariance###############

q<-v/(v+l*lag-2)
Ct<-Ctt<-vector(length=(T-lag)) 
U<-muy<-muyy<-muy1<-matrix(nrow=l,ncol=(T-lag))
for(i in 1:(T-lag)) 
{
Z11<-matrix(nrow=l,ncol=lag)
for(s in 1:lag)
{
Z11[,s]<-c(M0[1:l]+m1[,s+1,i]+m2[,s+1,i])
}
  Z1<-c(Z11)
if(tre=="Q")  M<-M0+c(m1[,,i])+c(m2[,,i])
if(tre=="L")  M<-M0+c(m1[,,i])
if(tre=="C")  M<-M0

  Ct[i] <- 1 + t(Z[(l+1):((lag+1)*l),i]-M[(l+1):((lag+1)*l)])%*%Q%*%(Z[(l+1):((lag+1)*l),i]-M[(l+1):((lag+1)*l)])
if(tre!="C")  Ctt[i] <- 1 + t(Z[(l+1):((lag+1)*l),i])%*%Q%*%(Z[(l+1):((lag+1)*l),i]) #-M0[(l+1):((lag+1)*l)]
if(tre=="C")  Ctt[i] <- Ct[i]

if(tre=="Q")  B0<-M0[1:l]+m1[,1,i]+m2[,1,i]-t(B1)%*%Z1
if(tre=="L")  B0<-M0[1:l]+m1[,1,i]-t(B1)%*%Z1
if(tre=="C")  B0<-M0[1:l]-t(B1)%*%Z1

  muy[,i] <- B0 + t(B1)%*%Z[(l+1):(l*(lag+1)),i] 
  #if(tre!="L")
  muyy[,i] <- t(B1)%*%Z[(l+1):(l*(lag+1)),i] 
  #if(tre=="C") muyy[,i]<-muy[,i]
  U[,i] <- t(Z[1:l,i])-muy[,i]   
}


#############################################################
if(hes=="TRUE") VARth<-solve(op$hessian[])
############################Coefficients/SEs/P-values########

################
##Jacobian##SE##
################

if(hes=="TRUE" & tre=="Q") 
{Jc<-Jacob(a,lag,l,v,tre)$J #jacobian(J,a)
SE<-sqrt(diag(Jc%*%VARth%*%t(Jc)))
p_value<-2*(1-pt(abs(c(Delta0,Delta1,Delta2,B1,vech(s2)))/SE,(T-lag)))
COEF<-round(cbind(c(Delta0,Delta1,Delta2,B1,vech(s2)),c(SE[1:length(SE)]),c(p_value[1:length(SE)])),8)
}
if(hes=="FALSE" & tre=="Q") COEF<-round(cbind(c(Delta0,Delta1,Delta2,B1,vech(s2))),8)

if(hes=="TRUE" & tre=="L") 
{Jc<-Jacob(a,lag,l,v,tre)$J #jacobian(J,a)
SE<-sqrt(diag(Jc%*%VARth%*%t(Jc)))
p_value<-2*(1-pt(abs(c(Delta0,Delta1,B1,vech(s2)))/SE,(T-lag)))
COEF<-round(cbind(c(Delta0,Delta1,B1,vech(s2)),c(SE[1:length(SE)]),c(p_value[1:length(SE)])),8)
}
if(hes=="FALSE" & tre=="L") COEF<-round(cbind(c(Delta0,Delta1,B1,vech(s2))),8)


if(hes=="TRUE" & tre=="C") 
{Jc<-Jacob(a,lag,l,v,tre)$J #jacobian(J,a)
SE<-sqrt(diag(Jc%*%VARth%*%t(Jc)))
p_value<-2*(1-pt(abs(c(Delta0,B1,vech(s2)))/SE,(T-lag)))
COEF<-round(cbind(c(Delta0,B1,vech(s2)),c(SE[1:length(SE)]),c(p_value[1:length(SE)])),8)
}
if(hes=="FALSE" & tre=="C") COEF<-round(cbind(c(Delta0,B1,vech(s2))),8)


if (hes=="TRUE") COEF<-cbind("estimate"=COEF[,1],"Std. Error"=COEF[,2],"p-value"=COEF[,3])
if (hes=="FALSE") COEF<-cbind("estimate"=COEF)

#########################MS-tests#############################

MS<-MS(U,lag,tre,s2,Ct,Ctt,muy,muyy,T,q,v)
Dist<-MS$Dist
MS<-MS$MS


if(hes=="TRUE" & tre=="Q") result<-list(beta=cbind(Delta0,Delta1,Delta2,t(B1)),coef=COEF,like=-Like,sigma=s2,cvar=Ct,res=t(U),fitted=t(muy),ms=MS,dist=Dist,init=a,hes=op$hes,S=S)
if(hes=="FALSE" & tre=="Q") result<-list(beta=cbind(Delta0,Delta1,Delta2,t(B1)),coef=COEF,like=-Like,sigma=s2,cvar=Ct,res=t(U),fitted=t(muy),ms=MS,dist=Dist,init=a,S=S)

if(hes=="TRUE" & tre=="L") result<-list(beta=cbind(Delta0,Delta1,t(B1)),coef=COEF,like=-Like,sigma=s2,cvar=Ct,res=t(U),fitted=t(muy),ms=MS,dist=Dist,init=a,hes=op$hes,S=S)
if(hes=="FALSE" & tre=="L") result<-list(beta=cbind(Delta0,Delta1,t(B1)),coef=COEF,like=-Like,sigma=s2,cvar=Ct,res=t(U),fitted=t(muy),ms=MS,dist=Dist,init=a,S=S)

if(hes=="TRUE" & tre=="C") result<-list(beta=cbind(Delta0,t(B1)),coef=COEF,like=-Like,sigma=s2,cvar=Ct,res=t(U),fitted=t(muy),ms=MS,dist=Dist,init=a,hes=op$hes,S=S)
if(hes=="FALSE" & tre=="C") result<-list(beta=cbind(Delta0,t(B1)),coef=COEF,like=-Like,sigma=s2,cvar=Ct,res=t(U),fitted=t(muy),ms=MS,dist=Dist,init=a,S=S)

return(result)
}
