library(tidyverse)
# load files
info <- read.table("pop.info") #sample information
likelihoods <- read.table("input.gz", header=T) #likelihood for columns
genotypes <- read.table("input.geno") # true genotypes

#1.2.2
#1
#subset dataframe

#head(likelihoods[,c(1:3, col:(col+2))])
get_columns <- function(data, ind){
  col <- 1+3*(ind)
  cols <- col:(col+2)
  return(data[,cols])
}

ind <- 61
lik_61 <- get_columns(likelihoods, ind)
head(lik_61)

#get posterior
get_posterior <- function(likelihood, prior){
  numerator <- likelihood * prior
  denominator <- apply(numerator, 1, sum)
  posterior <- numerator/denominator
  return(posterior)
}

posterior_1 <- get_posterior(lik_61, 1/3)
called_geno_1 <- apply(posterior, 1, max)
hist(called_geno_1)

#2
freq <- read.table("assign3.fopt.gz")

info[info$V2 == "NA12750",]
#Since individual is european, pick the second freq column

get_genotype_freq <- function(q){
  #q <- 1 - p
  p <- 1 - q
  return(data.frame("RR" = p^2, "RA" = 2*p*q, "AA" = q^2))
}

prior_2 <- get_genotype_freq(freq$V2)
posterior_2 <- get_posterior(lik_61, prior_2)
called_geno_2 <- apply(posterior_2, 1, max)
hist(called_geno_2)

#3
beagle_probs <- read.table("imputation.input.gz.gprobs.gz", header = T)
posterior_3 <- get_columns(beagle_probs, ind)
called_geno_3 <- apply(posterior_3, 1, max)
hist(called_geno_3)

#1.2.3
plotAccuracy<-function(x,p,...){
  p<-p[!is.na(x)]
  x<-x[!is.na(x)]
  ord<-order(p,decreasing=T)
  lines(1:length(x)/length(x),cumsum(!x[ord])/1:length(x),...)
}


get_predictions <- function(post, max, true){
  pred <- sweep((post == max), 2, c(1,2,3), `*`)
  pred <- apply(pred, 1, sum) - 1
  return(pred == true)
}

true <- genotypes$NA12750
pred_1 <- get_predictions(posterior_1, called_geno_1, true)
pred_2 <- get_predictions(posterior_2, called_geno_2, true)
pred_3 <- get_predictions(posterior_1, called_geno_3, true)


plot(1,xlim=0:1,ylim=c(0,0.40),col="transparent",xlab="callrate",ylab="error rate")
plotAccuracy(pred_1,called_geno_1,lwd=3,col="hotpink")
plotAccuracy(pred_2,called_geno_2,lwd=3,col="red")
plotAccuracy(pred_3,called_geno_3,lwd=3,col="blue")


###########################
# plot for ind2

ind2 <- 3
lik_3 <- get_columns(likelihoods, ind2)
true2 <- genotypes$NA19663

#get predictions and posteriors for 3 pops
get_pred_maxpost <- function(pop, lik, true){
  prior <- get_genotype_freq(pop)
  post <- get_posterior(lik, prior)
  max_post <- apply(post, 1, max)
  pred <- get_predictions(post, max_post, true)
  return(list(pred, max_post))
}

list_3_pops <- lapply(freq, get_pred_maxpost, lik_3, true2)

#get predictions and posteriors for combined pop
admix <- read.table("assign3.qopt")
prior_combined <- (admix$V1 * get_genotype_freq(freq$V1)) +
  (admix$V2 * get_genotype_freq(freq$V2)) +
  (admix$V3 * get_genotype_freq(freq$V3))

posterior_combined <- get_posterior(lik_3, prior_combined)
called_geno_comb <- apply(posterior_combined, 1, max)
pred_comb <- get_predictions(posterior_combined, called_geno_comb, true2)

# plot
plot(1,xlim=0:1,ylim=c(0,0.40),col="transparent",xlab="callrate",ylab="error rate")
for (i in 1:length(list_3_pops)){
  plotAccuracy(list_3_pops[[i]][[1]],list_3_pops[[i]][[2]],lwd=3,col=i)
}
plotAccuracy(pred_comb,called_geno_comb,lwd=3,col=4)
legend(0, 0.4, legend = c("African", "European", "Asian", "Combined"), 
       col= c(1,2,3,4), lty=1)

