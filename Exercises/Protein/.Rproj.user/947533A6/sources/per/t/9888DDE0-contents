library(Biostrings)
library(tidyverse)

dna <- readDNAStringSet("native_DNA.fa")
protein <- translate(dna)
protein

mutate <- function(dna){
  pos <- sample(length(dna), 1)
  base <- sample(c("A", "C", "G", "T"), 1)
  dna[pos] <- base
  return (dna)
}

n_mutates <- function(dna, n){
  for(i in 1:n){
    dna <- mutate(dna)
  }
  return(dna)
}



mutated_dna <- n_mutates(dna$seq, 100)
mutated_protein <- translate(mutated_dna)
mutated_protein
protein$seq

hist(mutated_dna)
dna_freq <- alphabetFrequency(mutated_dna)[1:4]

df <- data.frame(Freq=dna_freq) %>% rownames_to_column("Base")
df %>%
  ggplot(aes(x=Base, y=Freq)) + geom_histogram(stat="identity")
