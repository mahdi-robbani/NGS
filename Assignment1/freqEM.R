
getLikes<-function(x,d=5,e=0.01,norm=FALSE){
  n<-length(x)
  dep<-rpois(n,d)
  nA<-rbinom(n,dep,c(e,0.5,1-e)[x+1])
  res<-rbind(dbinom(nA,dep,e),dbinom(nA,dep,0.5),dbinom(nA,dep,1-e))
  if(norm)
    res<-t(t(res)/colSums(res))
  res
}

estPem<-function(like){ #numeric optimazition EM
    pk<-1 #start guess
    pkTemp<-0.2
    while(abs(pk-pkTemp)>0.0001){
        pk<-pkTemp
        w0<-like[1,]*(1-pk)^2
        w1<-like[2,]*2*pk*(1-pk)
        w2<-like[3,]*pk^2
        pkTemp<-mean((w1+2*w2)/(2*(w0+w1+w2)))
    }
    pk
}

emLog <- function(p=0.1,loglike,iter=30,accu=0.000001,accu2=0){
    numInds = ncol(loglike)
    temp_p = p
  for(it in 1:iter){
    sum=0
    for(i in 1:numInds){
      W0=exp(loglike[i*3+0-2])*(1-p)^2
      W1=exp(loglike[i*3+1-2])*2*p*(1-p)
      W2=exp(loglike[i*3+2-2])*p^2
      sum = sum + (W1+2*W2)/(2*(W0+W1+W2))
    }

    p=sum/numInds
    if((p-temp_p<accu & temp_p-p<accu) | (p/temp_p<1+accu2 & p/temp_p>1-accu2))
      break
    temp_p=p;
}
    p
}

emLog(0.1,t(r))

ind<-1000
dep<-2
freqPop <- 0.1
genoTrue <- rbinom(ind,2,freqPop)
freqTrue <- mean(genoTrue)/2
l1<-getLikes(genoTrue,dep,norm=TRUE)

c(popFreq=freqPop,trueFreq=freqTrue,estFreq=estPem(l1))
















double **angsd::get3likesRMlow(funkyPars *pars){

  double **loglike = NULL;
  loglike = new double*[pars->numSites]; 
  for(int s=0;s<pars->numSites;s++)
    loglike[s] = new double[3*pars->nInd];
  
  for(int s=0;s<pars->numSites;s++){

    if(pars->keepSites[s]==0)
      continue;
    for(int i=0;i<pars->nInd;i++){
      
      //fprintf(stderr,"mm: %d\t%d\n",pars->major[s],pars->major[s]);
      //fprintf(stderr,"%s\t%d\t%c\t%c\t",pars->sites[s].chromo,pars->sites[s].position+1,intToRef[pars->major[s]],intToRef[pars->minor[s]]);

      loglike[s][i*3+0]=pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->major[s]]];
      loglike[s][i*3+1]=pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->minor[s]]];
      loglike[s][i*3+2]=pars->likes[s][i*10+angsd::majorminor[pars->minor[s]][pars->minor[s]]];
      if(loglike[s][i*3+0] < -20 && loglike[s][i*3+1] < -20 && loglike[s][i*3+2] < -20){
	loglike[s][i*3+0] = 0;
	loglike[s][i*3+1] = 0;
	loglike[s][i*3+2] = 0;
      }
    }
  }
  return loglike;

}
