Jacob <-
function(a,lag,l,v,tre)
{
J<-function(a)
{
	S<-BlockTop(a[((l*(l+1)/2)*(1-1)+1):((lag+3)*(l+1)*l/2)],l,lag)$S

	S11<-S[1:l,1:l] ; S12<-S[1:l,(l+1):(l*(lag+1))] ; S22<-S[(l+1):(l*(lag+1)),(l+1):(l*(lag+1))]
	SS<-solve(S22)  ; Q<-SS/v ; B1<-SS%*%t(S12) ; s2<-S11-S12%*%SS%*%t(S12) 

	M0<-rep(c(a[(((lag+3)*(l+1)*l/2)+1):(((lag+3)*(l+1)*l/2)+l)]),(lag+1)) ## Mean
	M1<-c(a[(((lag+3)*(l+1)*l/2)+l+1):(((lag+3)*(l+1)*l/2)+2*l)])
	M2<-c(a[(((lag+3)*(l+1)*l/2)+2*l+1):(((lag+3)*(l+1)*l/2)+3*l)]) 

	if(tre=="C")  M1<-M2<-rep(0,l)
	if(tre=="L")   M2<-rep(0,l)

	m12<-m22<-matrix(nrow=l,ncol=lag)
	for(s in 1:lag)
	{
		m22[,s]<-2*s*M2/1000
		m12[,s]<-(s*M1/1000+((s/1000)^(2))*M2)
	}

	m222<-c(m22) ; m122<-c(m12)
	Delta0<-M0[1:l]-t(B1)%*%c(rep(M0[1:l],lag)-m122)
	Delta1<-M1-t(B1)%*%c(rep(M1,lag)-m222) ; Delta2<-M2-t(B1)%*%rep(M2,lag)

if(tre=="Q")	Cc<-c(Delta0,Delta1,Delta2,B1,vech(s2))
if(tre=="L")	Cc<-c(Delta0,Delta1,B1,vech(s2))
if(tre=="C")	Cc<-c(Delta0,B1,vech(s2))

			Cc
}

return(list(J=jacobian(J,a)))

}
