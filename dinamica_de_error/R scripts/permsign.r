perm.sign<-function(x,y=FALSE,paired=FALSE,alternative=FALSE,B=1000,fun=FALSE){
	if(paired==FALSE){#one sample problem
		n<-dim(x)[1]
  		p<-dim(x)[2]
  	    ##########################################################
 	    if(p==1) cat("Permutation Test for Symmetry \n")
         if(p>1) cat("Multivariate Permutation Test for Symmetry \n")
         ##########################################################
         T<-array(0,dim=c((B+1),p))
         T[1,]<-apply(x,2,sum)
         for(bb in 2:(B+1)){
         	random.sign<-1-2*rbinom(n,1,.5)
     	    T[bb,] = t(x)%*%random.sign
     	 }
     	 T<-abs(T)
     	 #from statistic to p.value
     	 P<-t2p(T)
     	 if(p==1) return(list(p.value=P[1]))
     	 if(p>1){
     	 #Combination
          T2<-comb(P,fcomb=fun)
          P2<-t2p(T2)
          return(list(p.value=P2[1]))
          }
 	}#end one sample problem
 	if(paired==TRUE){#paired samples
 		  cat("Permutation Test for Paired Samples \n")
 		  n<-dim(x)[1]
 		  T<-array(0,dim=c((B+1),1))
         if(alternative=="greater"|alternative=="two.sided"){diff<-x-y}
     	  #if(alternative=="two.sided"){diff=abs(diff)}
     	  if(alternative=="less"){diff<-y-x}
          T[1,]<-apply(diff,2,sum)
          for(bb in 2:(B+1)){
          	random.sign<-1-2*rbinom(n,1,.5)
          	T[bb,] = t(diff)%*%random.sign
          }
          if(alternative=="two.sided"){T=abs(T)}
          P<-t2p(T)
          return(list(p.value=P[1]))
     }#end paired samples
 	   
 }
