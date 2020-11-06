# Section 4: Identification of caste differentially expressed genes
normal_ant = which(sampleTable$species %in% c("Aech",'Mpha',"Lhum",'Sinv',"Lnig"))
ortholog_exp.ant = ortholog_counts[,normal_ant]

ortholog_counts.ant = ortholog_counts[,normal_ant]
ortholog_counts.ant = ortholog_counts.ant[!apply(ortholog_counts.ant, 1, anyNA),]
ortholog_counts.ant.norm = matrix(as.integer(ortholog_counts.ant), ncol = dim(ortholog_counts.ant)[2],
                                  dimnames = list(rownames(ortholog_counts.ant),colnames(ortholog_counts.ant)))



get_dif_exp_genes <- function(target_species, model=1){
  if(model == 1){
    dds <- DESeqDataSetFromMatrix(
      ortholog_counts.ant.norm[,which(sampleTable.ant$species %in% target_species)], sampleTable.ant[which(sampleTable.ant$species %in% target_species),],
      design = ~caste)
  } else{
    dds <- DESeqDataSetFromMatrix(
      ortholog_counts.ant.norm[,which(sampleTable.ant$species %in% target_species)], sampleTable.ant[which(sampleTable.ant$species %in% target_species),],
      design = ~caste +species+caste:species)
  }
  dds = DESeq(dds)
  res.aech = results(dds, contrast = c("caste",c("Gyne",'Worker')),alpha = 0.05)
  return(res.aech)
}





get_overlapping_genecount <- function(model_type){
  ants_names <- c("Mpha", "Aech", "Lhum", "Lnig", "Sinv")
  dseq_results <- list()
  filtered_genes <- list()
  
  for(name in ants_names){
    dseq_results[[name]] <- get_dif_exp_genes(name, model_type)
    filtered_genes[[name]] <- dseq_results[[name]] %>% 
      as.data.frame %>% 
      rownames_to_column("Gene") %>% 
      filter(padj < 0.05) %>% 
      select(Gene)
  }
  
  combinations <- list()
  gene_count <- list()
  
  for(i in model_type:5){
    combinations[[i]] <- combn(ants_names, i) # list of species names combinations
    columns <- dim(combinations[[i]])[2]
    for(c in 1:columns){
      species <- combinations[[i]][,c] #vector of species names
      group <- paste0(species, collapse = "-")
      gene_list <- list() #list of vectors of gene names
      for(s in species){
        gene_list[[s]] <- filtered_genes[[s]]$Gene
      }
      gene_count[[group]] <- length(Reduce(intersect, gene_list))
    }
  }
  
  gene_count_df <- gene_count %>% as.data.frame() %>% t
  colnames(gene_count_df) <- paste0("Model", model_type)
  
  return(gene_count_df)
}



model_1 <- get_overlapping_genecount(1)
model_2 <- get_overlapping_genecount(2)

model_2 <- get_dif_exp_genes(c("Mpha", "Aech", "Lhum", "Lnig", "Sinv"), 2)

model_2 %>% 
  as.data.frame() %>%
  filter(padj < 0.05) %>%
  dim

Reduce(intersect, list(filtered_genes$Mpha$Gene, filtered_genes$Aech$Gene, filtered_genes$Lnig$Gene))


inner_join(filtered_genes$Mpha, filtered_genes$Aech, filtered_genes$Lnig)


for(name in names(filtered_genes)){
  print(paste("Number of overlapping genes in", name, ":", length(filtered_genes[[name]]$Gene) ))
}


summary(dseq_results$Mpha)
filtered_genes$Mpha %>% length




x$rowname
str_count(x$rowname, "[.]")


####################old
# Section 4: Identification of caste differentially expressed genes
normal_ant = which(sampleTable$species %in% c("Aech",'Mpha',"Lhum",'Sinv',"Lnig"))
ortholog_exp.ant = ortholog_counts[,normal_ant]

ortholog_counts.ant = ortholog_counts[,normal_ant]
ortholog_counts.ant = ortholog_counts.ant[!apply(ortholog_counts.ant, 1, anyNA),]
ortholog_counts.ant.norm = matrix(as.integer(ortholog_counts.ant), ncol = dim(ortholog_counts.ant)[2],
                                  dimnames = list(rownames(ortholog_counts.ant),colnames(ortholog_counts.ant)))



get_dif_exp_genes <- function(target_species){
  dds <- DESeqDataSetFromMatrix(ortholog_counts.ant.norm[,which(sampleTable.ant$species %in% target_species)], sampleTable.ant[which(sampleTable.ant$species %in% target_species),],design = ~caste)
  dds = DESeq(dds)
  res.aech = results(dds, contrast = c("caste",c("Gyne",'Worker')),alpha = 0.05)
  return(res.aech)
}




ants_names <- c("Mpha", "Aech", "Lhum", "Lnig", "Sinv")
combinations <- list()
dseq_results <- list()
filtered_genes <- list()
for(i in 1:5){
  combinations[[i]] <- combn(ants_names, i)
  columns <- dim(combinations[[i]])[2]
  for(c in 1:columns){
    species <- combinations[[i]][,c]
    group <- paste0(species, collapse = "-")
    dseq_results[[group]] <- get_dif_exp_genes(species)
    filtered_genes[[group]] <- dseq_results[[group]] %>% 
      as.data.frame %>% 
      rownames_to_column("Gene") %>% 
      filter(padj < 0.05) %>% 
      select(Gene)
  }
}


for(name in names(filtered_genes)){
  print(paste("Number of overlapping genes in", name, ":", length(filtered_genes[[name]]$Gene) ))
}



