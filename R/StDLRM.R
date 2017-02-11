StDLRM <-
function(y,X,Trend=1,lag=1,v=1,maxiter=1000,meth="BFGS",hes="FALSE",init="na")
{
Data <- cbind(y,X)
if (lag-round(lag)!=0) return("lag must be an integer")
if (lag<0) return("lag number must be positive")

if (length(y)!= nrow(cbind(X))) return("y and X matrix must have same number of rows")

Z <- t(embed(cbind(y,X),lag+1)) ; T <- nrow(cbind(X)) ; l <- ncol(cbind(y,X))  ## No. of variables in VAR
if (mean(Trend)==1 & length(Trend)==1) {Trend <- matrix(1,T,1); TREND <- embed(Trend,lag+1) ; c <- 1}

if (mean(Trend)==0 & length(Trend)==1) {Trend <- matrix(0,T,1); TREND <- embed(Trend,lag+1) ; c <- 0}

if(ncol(cbind(Trend))>1) {TREND <- embed(Trend,lag+1) ; c <- ncol(Trend)}
if(Trend[1]!=1) {if (is.null(colnames(Trend))=="TRUE") {colnames(Trend) <- paste("t",".",seq(1,ncol(Trend),1),sep="")}}

if (Trend[1]!=1) {if (nrow(cbind(X))!= nrow((Trend))) return("Data and Trend matrix must have same number of rows")}
if (ncol(cbind(X))<1) return("no. of variables must be at least 1")
if (Trend[1]!=1) {if (ncol(Trend)<1) return("The Trend matrix must have at least one column")}
if (v<0) return("degrees of freedom must be greater or equal to 0")
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

if(init[1]!="na") int <- init
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

PP <- Par.dlrm(a,S,MU,l,lag,v,T,M,c,TREND)
Delta <- PP$Delta ; Delta0 <- PP$Delta0
B1 <- PP$B1 ; s2 <- PP$s2 ; Q <- PP$Q
var.coef <- PP$var.coef

beta <- cbind(Delta0, t(B1))

###################################################################
#############Fitted values/Residuals/Con. Covariance###############
###################################################################
q <- v/(v+l*lag+l-2)
Ct <- vector(length=(nrow(cbind(X))-lag)) 
U <- muy <- vector(length=(nrow(cbind(X))-lag))
trend <- vector(length=ncol(cbind(X))-1)
for(i in 1:(nrow(cbind(X))-lag)) 
{
  if(c!=0) {Ct[i] <- 1 + t(Z[2:ncol(S),i]-M[2:ncol(S)])%*%Q%*%(Z[2:ncol(S),i]-M[2:ncol(S)])}
  if(c==0) {Ct[i] <- 1 + t(Z[2:ncol(S),i])%*%Q%*%(Z[2:ncol(S),i])}
  if(c!=0) {trend[i] <- Delta0%*%(Trend[i,])} ; if(c==0) {trend[i] <- 0%*%(Trend[i,])}
  muy[i] <- trend[i] + t(B1)%*%Z[c(2:l,rep(1:l,lag)),i] 
  U[i] <- t(Z[1,i])- muy[i]   
}

#############################################################
if(hes=="TRUE") VARth <- solve(op$hessian)
############################Coefficients/SEs/P-values########
################
##Jacobian##SE##
################

if(hes=="TRUE") 
{
Jc <- Jacob.dlrm(a,lag,l,v,c)$J #jacobian(J,a)
SE <- sqrt(diag(Jc%*%VARth%*%t(Jc)))
p_value <- 2*(1-pt(abs(c(Delta0,B1,vech(s2)))/SE,(T-lag)))
COEF <- round(cbind(c(Delta0,B1,vech(s2)),c(SE[1:length(SE)]),c(p_value[1:length(SE)])),8)
colnames(COEF) <- c("coef.","std.err.","p-value")

Jv <- ConJacob.dlrm(a,lag,l,v,c)$Jv #jacobian(J,a)
SEv <- sqrt(diag(Jv%*%VARth%*%t(Jv)))
p_valuev <- 2*(1-pt(abs(c(var.coef))/SEv,(T-lag)))
VAR.COEF <- round(cbind(c(var.coef),c(SEv[1:length(SEv)]),c(p_valuev[1:length(SEv)])),8)
colnames(VAR.COEF) <- c("var.coef","std.err.","p-value")
}

if(hes=="FALSE") { 
COEF <- round(cbind(c(Delta0,B1,vech(s2))),8) ; colnames(COEF) <- c("coef.")
VAR.COEF <- round(cbind(var.coef),8) ; colnames(VAR.COEF) <- c("var.coef")

}

if(Trend[1]!= 0) beta <- cbind(Delta0,t(B1))
if(Trend[1] == 0) beta <- cbind(t(B1))

rownames(beta) <- paste(colnames(Data)[1])

head.x <- matrix(nrow=lag,ncol=ncol(Data)) ; for(i in 1:lag){head.x[i,] <- paste(colnames(Data),".",i,sep="")}

ifelse(mean(Trend[,1])==0 & ncol(Trend)==1, colnames(beta) <- c(colnames(X),t(head.x)), ifelse(mean(Trend)==1 & ncol(Trend)==1,colnames(beta) <- c("const.",colnames(X),t(head.x)),colnames(beta) <- c("const.",colnames(Trend)[2:ncol(Trend)],colnames(X),t(head.x))) )

head.var.coef <- matrix(nrow=(ncol(X)+1),ncol=lag) 
for(i in 1:lag)
{
head.var.coef[,i] <- c(paste(colnames(Data),".",i,sep=""))
}

coef.label <- c(colnames(X),head.var.coef)

var.head <- matrix(nrow=length(coef.label),ncol=length(coef.label))
for(i in 1:length(coef.label))
{
for(j in 1:length(coef.label))
{
var.head[i,j] <- paste(coef.label[j],".",coef.label[i],sep="")
}
}
res <- cbind(U)

rownames(VAR.COEF) <- c("const.",vech(var.head))

ad <- Student(res,Data,Ct,s2,v,lag,B1)$ad

if(hes=="TRUE") result <- list(beta=beta,coef=COEF,var.coef=VAR.COEF,like=-Like,sigma=s2,cvar=Ct,trend=trend,res=res,fitted=muy,init=a,hes=hess,S=S,ad=ad)
if(hes=="FALSE") result <- list(beta=beta,coef=COEF,var.coef=VAR.COEF,like=-Like,sigma=s2,cvar=Ct,trend=trend,res=res,fitted=muy,init=a,S=S,ad=ad)

return(result)
}