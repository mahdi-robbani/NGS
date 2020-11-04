load('inputData.Rdata')
source('shared_functions.R')
# Note 1. ortholog_exp and ortholog_counts are the expression matrix (Transcripts Per Kilobase Million, TPM) and the read count matrix for the 6672 one-to-one orthologous genes across 7 ant species.
# Note 2. sampleTable is the sample information for ortholog_exp, rownames are sample ID.
# Species: Aech, Mpha, Lhum, and Lnig are ants with caste system (life-long reproductive division of labour, queen (reproductive) and worker (sterile)). Dqua and Cbir are queenless ants, they have lost the queen caste and secondarily evolved reproductive worker caste.
# Lab: Samples of these seven ant species were from different labs.
# Note 3. Phylogeny of the seven ant species:
ant_tree = read.newick(text = ant_tree.data)
plot(ant_tree)
table(sampleTable$lab,sampleTable$species)
# Section 1: Normalization, Hierarchy clustering and PCA.

# Task 1.1: Try different normalization methods with ortholog_exp
ortholog_exp.norm = ortholog_counts #
ortholog_exp.norm = log2(ortholog_counts + 1) 
ortholog_exp.norm = log2(normalize.quantiles(ortholog_counts)+1)
#e.g. ortholog_exp.norm = log10(ortholog_counts + 1) Question: Why adding 1?
#e.g. ortholog_exp.norm = log10(normalize.quantiles(ortholog_exp)+1)

colnames(ortholog_exp.norm) = colnames(ortholog_counts)
rownames(ortholog_exp.norm) = rownames(ortholog_counts)
ortholog_exp.norm = ortholog_exp.norm[!apply(ortholog_exp.norm, 1, anyNA),] #Removed genes showing NA (e.g. without expression)
plot(density(ortholog_exp.norm[,21],na.rm=T),col = sp_color[sampleTable$species[1]],lwd=3, ylim = c(0,.3))
for(i in c(2:62)){
  lines(density(ortholog_exp.norm[,i],na.rm=T),col = sp_color[sampleTable$species[i]])}


# Task 1.2: Try different distance measuring methods, 
sampleDists = dist(t(ortholog_exp.norm))
sampleDists = as.dist(1 - cor(ortholog_exp.norm, method = 's', use = 'c'))
#e.g. sampleDists = dist(t(ortholog_exp.norm))
#e.g. sampleDists = as.dist(1 - cor(ortholog_exp.norm,method = 'p', use = 'c')) 

# A simple way to look at the data
sample_cluster = hclust(sampleDists)
plot(sample_cluster)
# For better visualization, we plot the hclust result in heatmap:
pheatmap(sampleDists,annotation_col = sampleTable[,c(1:3)], 
         annotation_colors = ann_colors, 
         # Take a look at the ann_colors in shared_function.R.
         color = colors) 

# Question 1: How does the normalization method influence the clustering result?
# Question 2: How does the distance measuring method influence the result?
# Try:
# ortholog_exp.norm = ortholog_exp
# sampleDists = dist(t(ortholog_exp.norm))

# Question 3:
# What are the dominant factors in transcriptomes?
# Mpha and Sinv belong to the same subfamily Myrmicinae (see ant_tree), why is Mpha more similar to Lhum, which belong to subfamily Dolichoderinae?
# Try with different normalization methods, how does it change the pattern? Why?

# We can further examine the transcriptome patterns among samples using PCA:
n = 10
ortholog_exp.norm.pca <- PCA(t(ortholog_exp.norm),ncp = n, graph = FALSE)

# Take a look at the amount of variations explained by each PC.
fviz_eig(ortholog_exp.norm.pca, addlabels = TRUE,main = 'Explained variance for each PC')

pca.var = ortholog_exp.norm.pca$eig
pca.data = cbind(ortholog_exp.norm.pca$ind$coord[,c(1:n)],sampleTable)

ggplot(pca.data, aes(x = Dim.6, y = Dim.7, color = species, shape = caste)) +
  geom_point(size=3) +
  coord_fixed()+
  theme_bw()

# Question 4:
# How much percentage of variance are being explained by the first two PCs? (Hint: Take a look at pca.var)
# What are the dominant factors in the first two PCs of transcriptomes? (How can we test this association?)
# Try with different normalization methods, how does it change the pattern? 

# Try with ICA approach, does it give different result? (As homework?)
# e.g. 
# sample_genes = sample(dim(ortholog_exp.norm)[1],size = 1000)
# sample_genes = order(apply(ortholog_exp.norm, 1, var),decreasing = T)[1:1000]
# ortholog_exp.norm.ica =fastICA(t(ortholog_exp.norm[sample_genes,]),n.comp = n)
