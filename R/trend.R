trend <-
function(M1,M2,T,lag,l)
{
m2<-m1<-array(dim=c(l,lag+1,T-lag))
t<-t(seq((1-lag)/1000,(T-lag)/1000,1/1000))
M11<-M1%*%t ; M22<-M2%*%(t^2)


for(i in (T-lag):(1))
{
m1[,,(i)]<-M11[,(i+lag):(i)]
m2[,,(i)]<-M22[,(i+lag):(i)]
}
return(list(m1=m1,m2=m2))
}
