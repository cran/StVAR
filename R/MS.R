MS <-
function(U,lag,tre,s2,Ct,Ctt,muy,muyy,T,q,v)
{
l<-nrow(U)
#st<-rmvst((T-lag),l,rep(0,l),Omega=diag(1,l),alpha=rep(0,l),df=v+lag*l)

res<-matrix(nrow=(T-lag),ncol=l)
	for(i in 1:(T-lag))
		{
	      ##Spanos residual 
		res[i,]<-U[,i]/sqrt((q*diag(s2)*Ct[i])) #-sqrt((q*diag(s2)*Ct[i]))*st[i,]
		}
t<-seq(1,T)
##MS tests for StVAR###
	hom<-tinv2<-Ind2<-Lin<-tinv<-Ind<-vector(length=l)
	for(i in 1:l)
		{
		polt<-poly(t,4,raw=FALSE)[,1:4]
		yhatt<-muyy[i,] ; y.hatt<-poly(yhatt,2,raw=FALSE)[,1:2]
		ms.g<-lm(res[(3:(T-(lag))),i]~polt[3:(T-(lag)),]+y.hatt[(3:(T-(lag))),]+res[(2:(T-(lag+1))),i]+res[(1:(T-(lag+2))),i])
		cf<-names(coef(ms.g))
		Lin[i]<-linearHypothesis(ms.g,cf[7])[2,6] ##Linearity
if(tre=="Q") tinv[i]<-linearHypothesis(ms.g,cf[4:5])[2,6] ##t-invariance
if(tre=="L") tinv[i]<-linearHypothesis(ms.g,cf[3:5])[2,6] ##t-invariance
if(tre=="C") tinv[i]<-linearHypothesis(ms.g,cf[2:5])[2,6] ##t-invariance
		
		Ind[i]<-linearHypothesis(ms.g,cf[8:9])[2,6]  ## Independence

		var.h<-q*s2[i,i]*Ct[1:(T-lag)]
		var.hatt<-q*s2[i,i]*Ctt[1:(T-lag)] ; var.hh<-poly(var.hatt,2,raw=FALSE)[,1:2]

if(tre=="Q") tr<-6
if(tre=="L") tr<-4
if(tre=="C") tr<-2

		poltt<-poly(t,tr,raw=FALSE)[,1:tr]

#summary(ms.g2)

		ms.g2<-lm(I(res[3:(T-lag),i]^2)~var.h[3:(T-lag)]+var.hh[3:(T-lag),2]+var.h[2:(T-lag-1)]+var.h[1:(T-lag-2)]+poltt[3:(T-lag),(tr-1):tr])
		cf<-names(coef(ms.g2))
if(tre=="Q") tinv2[i]<-linearHypothesis(ms.g2,cf[6:7])[2,6] ##t-inv
if(tre=="L") tinv2[i]<-linearHypothesis(ms.g2,cf[6:7])[2,6] ##t-inv
if(tre=="C") tinv2[i]<-linearHypothesis(ms.g2,cf[6:7])[2,6] ##t-inv

		hom[i]<-linearHypothesis(ms.g2,cf[3])[2,6] ##homo
		Ind2[i]<-linearHypothesis(ms.g2,cf[4:5])[2,6] ##Indep

		}

##Skewness-Kurtosis test
res1<-matrix(nrow=T-lag-2,ncol=l)
for(i in 1:l)
{
res1[,i]<-res[res[,i]!=max(res[,i]) & res[,i]!=min(res[,i]),i]
}

sk<-skewness(res1) ; kt<-kurtosis(res1)
N<-nrow(res1) ; kap<-6/(v+lag*l-4)
SK2<-(N/(6*(1+kap)))*sk^2+(N/(24*(1+kap)))*((kt-3-kap)^2)
p_chi<-pchisq(SK2,2,lower.tail=FALSE)
Students<-rbind(sk,kt,SK2,p_chi)

kolmo<-matrix(nrow=2,ncol=l)
for(i in 1:l)
{
ks<-ks.test(res[,i],"pt",df=v+lag*l,alternative="greater")
kolmo[,i]<-rbind(ks$statistic,ks$p.value)
}

if(v+lag*l<=4) Dist<-round(matrix(rbind(kolmo[1:2,]),2,l,dimnames=list(c("kolmo","p_kolmo"),c(rownames(U)))),2)

if(v+lag*l>4) Dist<-round(matrix(rbind(kolmo[1:2,],Students[3:4,]),4,l,dimnames=list(c("kolmo","p_kolmo","SK","p_SK"),c(rownames(U)))),4)



MS<-round(matrix(rbind(Ind,Lin,tinv,hom,Ind2,tinv2),6,l,dimnames=list(c("Independence","Linearity","t-Invariance","Homosk","2nd Independ","2nd t-Inv"),c(rownames(U)))),4)

return(list(MS=MS,Dist=Dist))
}
