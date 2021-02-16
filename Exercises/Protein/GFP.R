library(Biostrings)
library(tidyverse)

#matrix
native <- readDNAStringSet("native_DNA.fa")
N <- str_length(native$seq)


#load  and trim datadata
dim_counts <- read.table("dim_GFP_beads.counts", stringsAsFactors = F)
bright_counts <- read.table("bright_GFP_beads.counts", stringsAsFactors = F)


#filter data

filter_df_1 <- function(df, n){
  df <- df[df$V1 > 4,]
  df <- df[!grepl("[^ATCG]", df$V2),] #remove anything without the 4 bases
  df <- df[str_length(df$V2) == n,]
  return(df)
}

conv_and_translate <- function(string){
  string <- DNAString(string)
  string <- translate(string)
  string <- toString(string)
  return(string)
}

dim_filter <- filter_df(dim_counts, N)

x <- dim_filter[1:100,]
x$V3 <- as.character(lapply(x$V2, conv_and_translate))
x <- x[!grepl(".*[*].", x$V3),]


#translate

native_aa <- translate(native$seq)
native_aa_split <- str_split(native_aa, "")[[1]]

get_aa_diff <- function(mutated_aa, native_split){
  mutated_aa <- str_split(mutated_aa, "")[[1]]
  diff_vector <- native_split != mutated_aa
  if(sum(diff_vector) > 0){
    native_base <- native_split[diff_vector]
    mutated_base <- mutated_aa[diff_vector]
    pos <- (1:length(native_split))[diff_vector]
    #return(paste0(native_base, mutated_base))
    return(data.frame(native_base, mutated_base, pos))
   }
}

aa_df <- do.call("rbind", lapply(x$V3, get_aa_diff, native_aa_split))
heatmap(table(aa_df[,1:2]), Rowv = NA, Colv = NA, col=heat.colors(12), legend=c("none", "col"))

ggplot(aa_df, aes(x=native_base, y=mutated_base)) + geom_tile()


