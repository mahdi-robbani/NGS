library(R.oo) #package with charToInt function to convert ascii to a value

#####functions
get_base_probability <- function(genotype, base, error){
  #for a single genotype and a vector of bases and error,
  #return quality
  if (base == genotype | base == tolower(genotype)){
    return(1 - error)
  } 
  else {
    return(error/3)
  }
}

get_new_thetas <- function(bases_vector, errors_vector){
  #set initial frequencies and genotypes and varibales
  thetas_new <- c(0.25, 0.25, 0.25, 0.25)
  genotypes <- c("A", "C", "G", "T")
  difference <- 1
  iteration <- 0
  
  #base_prob <- list(A, T, C, G), A = [b1..bN], T = [b1..bN]...
  base_probabilities <- lapply(genotypes, 
                               get_base_probability, bases_vector, errors_vector)
  
  #run EM algo
  while (difference > 0.0001 | iteration < 1000){
    thetas <- thetas_new
    
    #Q Step #calculate posterior
    numerators <- mapply(function(x, y) x*y, base_probabilities, as.list(thetas))
    denominator <- apply(numerators, 1, sum)
    posterior <- numerators/denominator
    
    #M Step #calculate sum of posterior
    posterior_sum <- apply(posterior, 2, sum)
    denominator <- sum(posterior_sum)
    
    #calculate new thetas
    thetas_new <- posterior_sum/denominator
    difference <- max(abs(thetas_new - thetas))
    iteration <- iteration + 1
  }
  
  return(thetas_new)
}

#test
#get_new_thetas(bases[[2]][1:10], errors[[1]][1:10])
#get_new_thetas(bases[[302]], errors[[302]])

#load file
dat <- read.delim("data/MTnice.pileup", as.is=T, comment.char="", head=F, quote="")
names(dat) <- c("CHR","POS","REF", c("depth","bases","Qscore"))

## bases for individual 1 as a list
bases <- strsplit(dat$bases,"")

## ascii qualities for individual 1
asciiQ <- strsplit(dat$Qscore,"") #ascii code
Q <- lapply(asciiQ, function(x) charToInt(x) - 33) # quality score
errors <- lapply(Q, function(x) 10^(-x/10)) # error probability

#delte useless files
rm(dat, asciiQ, Q)

#calculate allele frequencies
allele_freqs <- mapply(get_new_thetas, bases, errors)
allele_freqs <- t(allele_freqs)
head(allele_freqs)

#find maximum allele frequncies
max_freq <- apply(allele_freqs, 1, max)
#find max allele frequncies less than 0.9
uncertain_sites <- allele_freqs[max_freq < 0.9,]
dim(uncertain_sites)[1]
#errors_sites <- allele_freqs[max_freq == 0,]