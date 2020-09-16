library(R.oo) #package with charToInt function to convert ascii to a value

#load file
dat <- read.delim("MTnice.pileup", as.is=T, comment.char="", head=F)
names(dat) <- c("CHR","POS","REF", c("depth","bases","Qscore"))

## bases for individual 1 as a list
bases <- strsplit(dat$bases,"")

## ascii qualities for individual 1
asciiQ <- strsplit(dat$Qscore,"") #ascii code
Q <- lapply(asciiQ, function(x) R.oo::charToInt(x) - 33) # quality score
E <- lapply(Q, function(x) 10^(-x/10)) # error probability


## row 2
bases[[1]]
## row 2 col 3
b_sub <- bases[[1]][1:10]
e_sub <- E[[1]][1:10]


#test

EMcoinStep<-function(pk,data){
  lik <- cbind(#likelihood of data assuming
    dbinom(data,10,pk[1]), #coin A
    dbinom(data,10,pk[2]) #coin B
  )
  ##probabilty of coin given data
  wA<-lik[,1]*0.5 / (lik[,1]*0.5 + lik[,2]*0.5)
  wB<-lik[,2]*0.5 / (lik[,1]*0.5 + lik[,2]*0.5)
  c(#n+1 parameter guess
    sum(data*wA)/sum(data*wA + (10-data)*wA),
    sum(data*wB)/sum(data*wB + (10-data)*wB)
  )
}
data <- c(5,9,8,4,7)
EMcoinStep(c(0.6,0.5),data)