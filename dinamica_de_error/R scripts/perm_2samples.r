perm.2samples<-function(data,B=10000,alt){
	n=table(data[,1])   
	C=length(n)  
	contr = rep(1/n,n)

	contr[-c(1:n[1])]=-contr[-c(1:n[1])]
	
	contr<-round(contr,digits=4)
	T<-array(0,dim=c((B+1),1))
	T[1]=data[,2]%*%contr
	for(bb in 2:(B+1)){
		X.star<-sample(data[,2])
		T[bb]=X.star%*%contr
		}
	#P.dist<-t2p(T)
	if(alt=="greater") {
		P<-t2p(T)
	}
	if(alt=="less"){
		P<-t2p(-T)
	}
	if(alt=="two.sided"){
		P<-t2p(abs(T))
		#P=apply(p,1,function(x) 2*min(x,1-x))
		#P<-array(P,dim=c((B+1),1)) 
	}
	list(p.value=P[1],dist=P)
}

