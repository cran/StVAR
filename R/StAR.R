StAR <-
function(Data,Trend=1,lag=1,v=1,maxiter=1000,meth="BFGS",hes="FALSE",init="na")
{
Data <- cbind(Data)
if (ncol(Data)!=1) return("data must be in a vector form")
if (lag-round(lag)!=0) return("lag must be an integer")
if (lag<0) return("lag number must be positive")

Z <- t(embed(Data,lag+1)) ; T <- nrow(Data) ; l <- 1  ## No. of variables in VAR
if (Trend[1]==1 & length(Trend)==1) {Trend <- matrix(1,T,1); TREND <- embed(Trend,lag+1) ; c <- 1}

if (Trend[1]==0 & length(Trend)==1) {Trend <- matrix(0,T,1); TREND <- embed(Trend,lag+1) ; c <- 0}

if(ncol(cbind(Trend))>1) {TREND <- embed(Trend,lag+1) ; c <- ncol(Trend)}
if(Trend[1]!=1) {if (is.null(colnames(Trend))=="TRUE") {colnames(Trend) <- paste("t",".",seq(1,ncol(Trend),1),sep="")}}

if (Trend[1]!=1) {if (nrow(Data)!= nrow(Trend)) return("Data and Trend matrix must have same number of rows")}
if (ncol(Data)>2) return("no. of variables must equal to 1")
if (Trend[1]!=1) {if (ncol(Trend)<1) return("The Trend matrix must have at least one column")}
if (v<0) return ("degrees of freedom must be greater or equal to 0")
if (maxiter<10) return("Iteration must be at least 10")
if (2*T<3*((lag+3)*(l+1)*l/2)+3*l+c) return(list(((lag+3)*(l+1)*l/2)+3*l+c,"Too many parameters for given sample size. Reduce the number of lags."))

##Likelihood Function
L<-function(a)
{
 S <- BlockTop(a[((l*(l+1)/2)*(1-1)+1):((lag+3)*(l+1)*l/2)],l,lag)$S ##Var-Cov Matrix
 F <- solve(S)/v ## VarCov and its Inverse

if(c!=0) {MU <- matrix(c(a[(((lag+3)*(l+1)*l/2)+1):(((lag+3)*(l+1)*l/2)+l*c)]),l,c)}

if(c!=0) {M0 <- kronecker(diag(1,lag+1),MU)}
if(c==0) {M0 <- kronecker(diag(1,lag+1),0)}

if(c!=0) {M <- M0%*%t(TREND)} ; if (c==0) {M <- 0}

D <- 1 + diag(t(Z-M)%*%F%*%(Z-M)) ##Quadratic form

## Likelihood function
 LLn <- (T-lag)*const - 0.5*(T-lag)*log(det(S)) - 0.5*(v+(lag+1)*l)*sum(log(D))
 neg.LLn <- -LLn
 neg.LLn

}

if(init[1]=="na") int <- c(runif(((lag+3)*(l+1)*l/2)+l*c,0,.1))  ##Initialization

if(init[1]!="na") int<-init
const <- log(gamma((v+(lag+1)*l)/2)) - log(gamma(v/2)) - 0.5*l*(lag+1)*log(pi*v)

op <- optim(int,L,hessian=hes,control=list(trace=1,maxit=maxiter,reltol=1e-14),method=meth) 
a <- op$par ; Like <- op$value ; hess <- op$hes

###################################Parameters##################
###################################Parameters##################
###################################Parameters##################
S <- BlockTop(a[((l*(l+1)/2)*(1-1)+1):((lag+3)*(l+1)*l/2)],l,lag)$S
if(c!=0) {MU <- matrix(c(a[(((lag+3)*(l+1)*l/2)+1):(((lag+3)*(l+1)*l/2)+l*c)]),l,c)}

if(c!=0) {M0 <- kronecker(diag(1,lag+1),MU)}
if(c==0) {M0 <- kronecker(diag(1,lag+1),0)}

if(c!=0) {M <- M0%*%t(TREND)} ; if (c==0) {M <- 0}

PP <- Par.star(a,S,MU,l,lag,v,T,M,c,TREND)
Delta <- PP$Delta ; Delta0 <- PP$Delta0
B1 <- PP$B1 ; s2 <- PP$s2 ; Q <- PP$Q
var.coef <- PP$var.coef

beta <- cbind(Delta0, t(B1))

###################################################################
#############Fitted values/Residuals/Con. Covariance###############
###################################################################
q <- v/(v+l*lag-2)
Ct <- vector(length=(T-lag)) 
U <- muy <- matrix(nrow=l,ncol=(T-lag))
trend <- vector(length=ncol(cbind(Data))-1)

for(i in 1:(T-lag)) 
{
  if(c!=0) {Ct[i] <- 1 + t(Z[(l+1):((lag+1)*l),i]-M[(l+1):((lag+1)*l)])%*%Q%*%(Z[(l+1):((lag+1)*l),i]-M[(l+1):((lag+1)*l)])}
  if(c==0) {Ct[i] <- 1 + t(Z[(l+1):((lag+1)*l),i])%*%Q%*%(Z[(l+1):((lag+1)*l),i])}
  if(c!=0) {trend[i] <- Delta0%*%(Trend[i,])} ; if(c==0) {trend[i] <- 0%*%(Trend[i,])}
  muy[,i] <- trend[i] + t(B1)%*%Z[(l+1):(l*(lag+1)),i] 
  U[,i] <- t(Z[1:l,i])- muy[,i]   
}

res <- t(U)
#############################################################
if(hes=="TRUE") VARth <- solve(op$hessian)
############################Coefficients/SEs/P-values########
################
##Jacobian##SE##
################

if(hes=="TRUE") 
{
Jc <- Jacob.star(a,lag,l,v,c)$J 	#jacobian(J,a)
SE <- sqrt(diag(Jc%*%VARth%*%t(Jc)))
p_value <- 2*(1-pt(abs(c(Delta0,B1,vech(s2)))/SE,(T-lag)))
COEF <- round(cbind(c(Delta0,B1,vech(s2)),c(SE[1:length(SE)]),c(p_value[1:length(SE)])),8)
colnames(COEF) <- c("coef.","std.err.","p-value")

Jv <- ConJacob.star(a,lag,l,v,c)$Jv #jacobian(J,a)
SEv <- sqrt(diag(Jv%*%VARth%*%t(Jv)))
p_valuev <- 2*(1-pt(abs(c(var.coef))/SEv,(T-lag)))
VAR.COEF <- round(cbind(c(var.coef),c(SEv[1:length(SEv)]),c(p_valuev[1:length(SEv)])),8)
colnames(VAR.COEF) <- c("var.coef","std.err.","p-value")
}

if(hes=="FALSE") { 
COEF <- round(cbind(c(Delta0,B1,vech(s2))),8) ; colnames(COEF) <- c("coef.")
VAR.COEF <- round(cbind(var.coef),8) ; colnames(VAR.COEF) <- c("var.coef")
}


head.var.coef <- matrix(nrow=l,ncol=lag) 
for(i in 1:lag)
{
head.var.coef[,i] <- paste(colnames(Data),".",i,sep="")
}

var.head <-matrix(nrow=lag*l,ncol=lag*l)
for(i in 1:(lag*l))
{
for(j in 1:(lag*l))
{
var.head[i,j] <- paste(c(head.var.coef)[j],".",c(head.var.coef)[i],sep="")
}
}

rownames(VAR.COEF) <- c("const.",vech(var.head))

if(Trend[1]!= 0) beta <- cbind(Delta0,t(B1))
if(Trend[1] == 0) beta <- cbind(t(B1))

rownames(beta) <- paste(colnames(Data))
head.x <- matrix(nrow=lag,ncol=ncol(Data)) ; for(i in 1:lag){head.x[i,] <- paste(colnames(Data),".",i,sep="")}

ifelse(mean(Trend[,1])==0 & ncol(Trend)==1, colnames(beta) <- c(t(head.x)), ifelse(mean(Trend[,1])==1 & ncol(Trend)==1,colnames(beta) <- c("const.",t(head.x)),colnames(beta) <- c("const.",colnames(Trend)[2:ncol(Trend)],t(head.x))))

fitted <- t(muy)
res <- cbind(U)

ad <- Student(res,Data,Ct,s2,v,lag,B1)$ad

if(hes=="TRUE") result <- list(beta=beta,coef=COEF,var.coef=VAR.COEF,like=-Like,sigma=s2,cvar=Ct,trend=trend,res=res,fitted=fitted,init=a,hes=hess,S=S,ad=ad)
if(hes=="FALSE") result <- list(beta=beta,coef=COEF,var.coef=VAR.COEF,like=-Like,sigma=s2,cvar=Ct,trend=trend,res=res,fitted=fitted,init=a,S=S,ad=ad)
return(result)
}
