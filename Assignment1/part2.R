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
  while (difference > 0.0001 | iteration > 10000){
    thetas <- thetas_new
    
    #E Step calculate posterior
    #calculate expected values
    numerators <- mapply(base_probabilities, as.list(thetas), FUN= function(x, y) x*y)
    denominator <- apply(numerators, 1, sum)
    expected_values <- numerators/denominator
    
    #M Step
    #calculate new thetas
    thetas_new <- apply(expected_values, 2, mean)
    difference <- max(abs(thetas_new - thetas))
    iteration <- iteration + 1
    
    #print
    # print("iteration")
    # print(iteration)
    # print("theta old, theta new")
    # print(thetas)
    # print(thetas_new)
    # print("difference")
    # print(difference)
  }
  
  return(thetas_new)
}


#load file
dat <- read.delim("MTnice.pileup", as.is=T, comment.char="", head=F, quote="")
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
colnames(allele_freqs) <- c("A", "C", "G", "T")
head(allele_freqs)


allele_freqs[apply(allele_freqs, 1, max) < 0.9,]



