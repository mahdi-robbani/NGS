## single round (E+M step) of the EM algorithm for the coin example
EMcoinStep<-function(pk,data){ #numeric optimazition EM
    lik <- cbind(#likelihood of data assuming 
        dbinom(data,10,pk[1]), #coin A
        dbinom(data,10,pk[2]) #coin B
    )
    ##probabilty of coin given data 
    wA<-lik[,1]*0.5 / (lik[,1]*0.5 + lik[,2]*0.5) #random coin/uniform prior
    wB<-lik[,2]*0.5 / (lik[,1]*0.5 + lik[,2]*0.5)

    c(#n+1 - parameter
        sum(data*wA)/sum(data*wA + (10-data)*wA), #10 tooses with each coin
        sum(data*wB)/sum(data*wB + (10-data)*wB)
    )
}

## run until converged of the EM algorithm for the coin example
EMcoin<-function(data){ #numeric optimazition EM
    pkTemp<-c(0.6,0.5)  #start guess
    while(any(abs(pk-pkTemp)>0.0001)){#run until next parameter guess is similar
        pk<-pkTemp
        lik <- cbind(#likelihood of data assuming 
            dbinom(data,10,pk[1]), #coin A
            dbinom(data,10,pk[2]) #coin B
        )
        #probabilty of coin given data 
        wA<-lik[,1]*0.5 / (lik[,1]*0.5 + lik[,2]*0.5) #random coin/uniform prior
        wB<-lik[,2]*0.5 / (lik[,1]*0.5 + lik[,2]*0.5)

        pkTemp<-c(
            sum(data*wA)/sum(data*wA + (10-data)*wA), #10 tooses with each coin
            sum(data*wB)/sum(data*wB + (10-data)*wB)
        )
        print(pkTemp)
    }
    pk
}

data <- c(5,9,8,4,7)
EMcoin(data)


### Removing the assumption of the coin prior of 50%

EMcoin2<-function(data){ #numeric optimazition EM
    pk<-0
    pkTemp<-c(0.6,0.5,0.1)  #start guess
    while(any(abs(pk-pkTemp)>0.0001)){#run until next parameter guess is similar
        pk<-pkTemp
        lik <- cbind(#likelihood of data assuming 
            dbinom(data,10,pk[1]), #coin A
            dbinom(data,10,pk[2]) #coin B
        )
        #probabilty of coin given data 
        wA<-lik[,1]*pk[3] / (lik[,1]*pk[3] + lik[,2]*(1-pk[3])) #random coin/uniform prior
        wB<-lik[,2]*(1-pk[3]) / (lik[,1]*pk[3] + lik[,2]*(1-pk[3]))

        pkTemp<-c(
            sum(data*wA)/sum(data*wA + (10-data)*wA), #10-data: #T
            sum(data*wB)/sum(data*wB + (10-data)*wB),
            sum(wA)/sum(wA+wB)
        )
        print(pkTemp)
    }
    pk
}

data <- 10 - c(5,1,2,6,3)
EMcoin2(data)
