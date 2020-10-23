# AS Oct 2020
# Advanced Topics in Bioinformatics
# homework assignment
library(tidyverse)
library(Biostrings)

base_dir <- "data"

# ordering of amino acid by biochemical properties. You don't need to use this. It makes it easier to spot clusters of similar effects that are due to biochemistry.
colnames.by.AA <- c("G", "A", "V", "L", "I", "M", "F", "W", "P", "S", "T", "C", "Y", "N", "Q", "D", "E", "K", "R", "H")

Sarkisyan.file <- paste(base_dir, "nt_sequences_to_brightness.csv", sep="/")
Sarkisyan.data <- read.csv(Sarkisyan.file, stringsAsFactors = F)
colnames(Sarkisyan.data)
# 

nativeDNA <- "AGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGTCGTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACACTAGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCACGGCATGGACGAGCTGTACAAGTGA"
nativeAA <- translate(DNAString(nativeDNA))

# Task 1: quality control and translation
# ---------------------------------------

# quality control step 1: remove sequences that are too long, too short, or have gaps
# useful function: nchar()

filter_sark <- function(df, n){
  df <- df[!grepl("[^ATCG]", df$sequence),] #remove anything without the 4 bases
  df <- df[str_length(df$sequence) == n,]
  return(df)
}

N <- str_length(nativeDNA)
Sarkisyan.data <- filter_sark(Sarkisyan.data, N)

# translate to protein
conv_and_translate <- function(string){
  string <- DNAString(string)
  string <- translate(string)
  string <- toString(string)
  return(string)
}

Sarkisyan.data$AAsequence <- as.character(lapply(Sarkisyan.data$sequence, 
                                                conv_and_translate))


# quality control step 2: remove sequences with premature STOP codons
data_no_stop <- Sarkisyan.data[!grepl(".*[*].", Sarkisyan.data$AAsequence),]

#filtered_data <- rownames_to_column(filtered_data, "SeqID")
data_no_stop %>% 
  group_by(AAsequence) %>%
  summarize(Barcode = sum(uniqueBarcodes),
            Brightness = mean(medianBrightness), 
            Error = sqrt(sum(stdErr^2, na.rm = T)) * 1/length(stdErr[!is.na(stdErr)]),
            DNA = sequence[1]) -> data_combined_AA

data_combined_AA$ID <- 1:length(data_combined_AA$AAsequence)

# Here, unlike in the GFP dataset we discussed in class, the uniqueBarcodes column tells us how often a particular sequence has been observed, while the medianBrightness tells us how bright fluorescence in those cells was (with stdErr estimating the measurement error, if available).
# How many unique barcodes (=DNA variant sequences) are found? How many unique protein sequences after cleanup? What is the most common protein sequence that is not wild-type? Include your answers in your hand-in.
paste("Unique Barcodes:", length(unique(data_no_stop$sequence)))
paste("Unique Proteins:", length(unique(data_no_stop$AAsequence)))

data_combined_AA %>% 
  filter(AAsequence != as.character(nativeAA)) %>%
  arrange(desc(Barcode)) %>%
  select(AAsequence) %>%
  head(1) %>% 
  as.character()


# Task 2: protein-level variants
# ------------------------------
# Next, determine differences to the native protein sequence. For simplicity we will consider each position independently, i.e., regardless of whether this mutation was only observed in context of other mutations. So, if a protein with A23P and S84Q is reported to have a medianBrightness of 3.5, then the brightness of A23P is 3.5 and the brightness of S84Q is 3.5 (in that protein). We also want to average the brightness across all contexts, so if i.e. there were 3 unique(!) proteins containing the A23P mutation, the averaged brightness of A23P is mean(brightness(protein1),brightness(protein2),brightness(protein3)).
#It may be convenient to generate a dataframe spanning all 20 possible amino acids at all positions, like we did in Exercise 2.
get_aa_diff <- function(ID, df, native){
  mutated <- df$AAsequence[df$ID == ID]
  mutated <- str_split(mutated, "")[[1]]
  diff_vector <- native != mutated
  diff_num <- sum(diff_vector)
  if(sum(diff_vector) > 0){
    native_base <- native[diff_vector]
    mutated_base <- mutated[diff_vector]
    pos <- (1:length(native))[diff_vector]
  } else{
    native_base <- "X"
    mutated_base <- "X"
    pos <- 0
  }
  mut <- paste0(native_base, pos, mutated_base)
  return(data.frame("Ref" = native_base, 
                    "Alt" = mutated_base, 
                    "Pos" = pos,
                    "Mutation" = mut,
                    "NumMuts" = diff_num,
                    "ID" = ID))
}

#remove old data
#rm(Sarkisyan.data)
#rm(data_no_stop)

#create native and find mutations
nativeAA <- translate(DNAString(nativeDNA))
nativeAA <- str_split(nativeAA, "")[[1]]
AA_diff_df <- bind_rows(lapply(data_combined_AA$ID, get_aa_diff, 
                               data_combined_AA, nativeAA))

#Merge
data_AA <- merge(x = AA_diff_df, y = data_combined_AA, by = "ID", all.x = TRUE)

data_AA %>%
  filter(Ref != "X") %>%
  group_by(Ref, Pos, Alt) %>%
  summarize(AvgBrightness = mean(Brightness), 
            AvgStdErr = sqrt(sum(Error^2, na.rm = T)) * 1/length(Error[!is.na(Error)])) %>%
  mutate(Mutation = paste0(Ref, Pos, Alt)) -> AA_avg_brightness

# - see discussion in class on Oct 6 for useful functions to determine differences between wild-type and variant sequence

# As a control for the averaged data across different sequences, create a subset of the dataset where only single-mutation sequences are considered. 
# useful functions: subset() and
substitution_distance <- function(s1,s2) { mapply(function(c1,c2) sum(c1!=c2), strsplit(s1,''), strsplit(s2,'')) }

data_AA %>%
  filter(Ref != "X") %>%
  filter(NumMuts == 1) %>%
  group_by(Mutation) %>%
  summarize(SingleBrightness = mean(Brightness), 
            StdErr = sqrt(sum(Error^2, na.rm = T)) * 1/length(Error[!is.na(Error)])) -> AA_brightness_control


# Then compare the medianBrightness of those single-mutant sequences to the averaged data you created above. In our example above that would mean comparing the brightness of the A23P single mutant protein to the average brightness over all proteins that contain an A23P mutation. You should get scatter plots analogous to those we discussed in class on Oct 9th, see also  https://www.biorxiv.org/content/biorxiv/early/2020/05/26/2020.05.26.116756/F7.large.jpg?width=800&height=600&carousel=1. Include the stdErr in the plot, using geom_errorbar() and geom_errorbarh(). Are the deviations you observe beyond what you expect based on the experimental error? Submit plot and discussion as part of your hand-in.
AA_brightness_combined <- inner_join(AA_avg_brightness, 
                                     AA_brightness_control, by="Mutation")
AA_brightness_combined %>%
  ggplot(aes(x=SingleBrightness, y=AvgBrightness)) + geom_point(col="blue") +
    geom_errorbar(aes(ymin=AvgBrightness-AvgStdErr, ymax=AvgBrightness+AvgStdErr), alpha=0.4) +
    geom_errorbarh(aes(xmin=SingleBrightness-StdErr, xmax=SingleBrightness+StdErr), alpha=0.4)
  

# Next, pick 2 amino acids from your first and last name, respectively -> AA1, AA1
# Visualise the distributions of mutations from AA1 and from AA2 to all other amino acids across all positions in the sequence. Then do the same for mutations from any amino acid into AA1 and AA2 - do they differ? What would you expect based on biochemistry vs. what do you observe?
# - for synonymous mutations?
# - for missense mutations?
# This is analogous to Exercise 3.
# Submit the plots for all 4 distributions as part of your handin. Of course you can plot all 4 distributions in the same figure, so long it is clear which line corresponds to which dataset.

compare_letters <- function(A, B, data){
  data %>%
    filter(Ref == A) %>%
    select(AvgBrightness) %>%
    mutate(Group = "To All", AA = A) -> dist1
  
  data %>%
    filter(Ref == B) %>%
    select(AvgBrightness) %>%
    mutate(Group = "To All", AA = B) -> dist2
  
  data %>%
    filter(Alt == A) %>%
    select(AvgBrightness) %>%
    mutate(Group = "To Res", AA = A) -> dist3
  
  data %>%
    filter(Alt == B) %>%
    select(AvgBrightness) %>%
    mutate(Group = "To Res", AA = B) -> dist4
  
  dist <- rbind(dist1,dist2,dist3,dist4)
  return(dist)
}

compare_letters("S", "P", AA_avg_brightness) %>%
  ggplot(aes(x=AvgBrightness, fill=Group)) + 
  geom_density(alpha=0.5) +
  facet_wrap(~AA)

# Task 3: summary matrix
# ----------------------
# summarise the results across all variants in a 20x20 matrix showing the wild-type and target amino acids, as we did in the exercises in class. Submit a plot of the matrix (see e.g. ex. 3) as part of your homework assignment. You can use the colnames.by.AA to sort the amino acids:
data_AA %>%
  select(Ref, Alt, Brightness) %>%
  add_row(Ref = colnames.by.AA, Alt = "X", Brightness = 0) %>%
  add_row(Alt = colnames.by.AA, Ref = "X", Brightness = 0) %>%
  group_by(Ref, Alt) %>%
  summarize(Brightness = mean(Brightness)) %>%
  filter(Ref != "X" & Alt != "X") -> mut_data

mut_data$Ref <- factor(mut_data$Ref, levels = colnames.by.AA) 
mut_data$Alt <- factor(mut_data$Alt, levels = colnames.by.AA)
mut_data %>%
  ggplot(aes(x=Ref, y=Alt, fill= Brightness)) + geom_tile() +
  scale_fill_gradientn(colours = c("white", "yellow", "red"), 
                       values = c(0,0.5,1))

 
# useful function: ddply(mut.data, ..., summarise, ...)
# https://colorbrewer2.org/ for colour schemes



# Task 4: compare to the other GFP mutagenesis dataset
# ----------------------------------------------------

# Task 4.1
# load in the native DNA from exercise 1
# compare it to the nativeDNA included above (e.g. by pairwise sequence alignment), then translate both sequences to protein and compare those. 
# Write a short paragraph describing what you observe.
nativeDNA_1 <- readDNAStringSet("data/native_DNA.fa")
pairwiseAlignment(nativeDNA, nativeDNA_1)

nativeAA_1 <- translate(nativeDNA_1$seq)
AA_alignment <- pairwiseAlignment(nativeAA, nativeAA_1)
AA_alignment


# Task 4.2
# Load in the GFP dataset we parsed in exercise 2, including the cleanup steps, translation to protein and identification of differences to the wt sequence. 
# How many variants (wt, position, mut.aa) are
#  - observed in both datasets? 
#  - only observed in the Sarkisyan dataset?
#  - only observed in the dataset we worked with in class?
filter_gfp <- function(df, n){
  df <- df[df$V1 > 4,]
  df <- df[!grepl("[^ATCG]", df$V2),] #remove anything without the 4 bases
  df <- df[str_length(df$V2) == n,]
  return(df)
}

translate_gfp <- function(df){
  df$V3 <- as.character(lapply(df$V2, conv_and_translate))
  df <- df[!grepl(".*[*].", df$V3),]
  colnames(df) <- c("Count", "DNA", "AAsequence")
  return(df)
}

get_diff_df_gfp <- function(seq_df, nativeDNA){
  #calculate length of DNA and create AA list
  N <- str_length(nativeDNA$seq)
  nativeAA <- translate(nativeDNA$seq)
  nativeAA <- str_split(nativeAA, "")[[1]]
  #filter, trnslate and add seqID
  seq_df <- filter_gfp(seq_df, N)
  seq_df <- translate_gfp(seq_df)
  seq_df <- rownames_to_column(seq_df, "ID")
  #calculate differences between native and seq
  diff_df <-  bind_rows(lapply(seq_df$ID, get_aa_diff, 
                               seq_df, nativeAA))
  diff_df$Pos <- diff_df$Pos + 132 #alignment
  diff_df$Mutation <- paste0(diff_df$Ref, diff_df$Pos, diff_df$Alt) #fix mutations
  # merge the dataframes
  seq_df <- merge(x = diff_df, y = seq_df, by = "ID", all.x = TRUE)
  return(seq_df)
}

# get dim freq
dim_df <- read.table("data/dim_GFP_beads.counts", stringsAsFactors = F)
dim_df <- get_diff_df_gfp(dim_df, nativeDNA_1)

dim_df %>%
  filter(Ref != "X") %>%
  group_by(Mutation) %>%
  summarize(FreqDim = n()) -> dim_freq


# get bright freq
bright_df <- read.table("data/bright_GFP_beads.counts", stringsAsFactors = F)
bright_df <- get_diff_df_gfp(bright_df, nativeDNA_1)


bright_df %>%
  filter(Ref != "X") %>%
  group_by(Mutation) %>%
  summarize(FreqBright = n()) -> bright_freq


#combined dim and bright to gfp
bright_dim_combined <- inner_join(bright_freq, dim_freq, by="Mutation")

#combine sark with gfp
combined_df <- inner_join(AA_avg_brightness, bright_dim_combined, by="Mutation")
paste("Combined data size:", dim(combined_df)[1])

paste("Sarkisyan data size:", sum(!AA_avg_brightness$Mutation %in% combined_df$Mutation))
paste("Exercise data size:", sum(!bright_dim_combined$Mutation %in% combined_df$Mutation))


# Task 4.3
# For the variants found in both datasets, create a scatterplot to compare their averaged medianBrightness (see task 2) vs. log(bright/dim) ratio. Briefly describe what trends you observe, and whether those are what you would expect.
# Submit the scatter plot and discussion as part of your hand-in.
combined_df %>%
  mutate(LogOdds = log10(FreqBright/FreqDim)) %>%
  ggplot(aes(x=AvgBrightness, y=LogOdds)) + geom_point() +
  geom_smooth(method = "lm")
