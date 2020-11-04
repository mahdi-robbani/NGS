# Prepare input: Read in count matrix produced by Salmon, matching orthologous genes.
# Salmon: https://combine-lab.github.io/salmon/getting_started/
# Proteinortho: https://www.bioinf.uni-leipzig.de/Software/proteinortho/manual.html

source('shared_functions.R')
input_dir = "./"
#Read in the ortholog table (based on reciprocal blastp and synteny information)
gene_ortholog_table = read.table(paste0(input_dir,'input/gene_table.poff'),
                                 col.names = c('Aech','Sinv','Mpha','Lnig','Lhum','Cbir','Dqua')) 

#####
# Read in expression data for each species
#Importing small worker and gyne brain samples of A.echinatior (we removed samples 1:3 because ERCC suggested poor quality of these samples)
Aech_exp = read_input('Aech',col_name = c(4,6,7,9,10,12,13,15), input_dir) 
#Samples of L.humile (We removed sample 1 and 2 because this colony was collected from different area)
Lhum_exp = read_input('Lhum',col_name = c(3:10), input_dir) 
Mpha_exp = read_input('Mpha',col_name = c(1:6,'7x',8,'9x',10), input_dir) #Samples of M. pharaonis
Sinv_exp = read_input('Sinv',col_name = c(3,4,7,8,11,12,15,16), input_dir) #Samples of S.invicta
Lnig_exp = read_input('Lnig',col_name = c(1:8), input_dir) #Samples of L.niger
Cbir_exp = read_input('Cbir',col_name = c(1:8), input_dir) # Samples of O.biroi
Dqua_exp = read_input('Dqua',col_name = c(1:5,7:13), input_dir) # D.quadriceps, removed 2 samples without colony information
#####
#Expression matrix for all 1-to-1 orthologous genes
ortholog_exp = cbind(Aech_exp[match(gene_ortholog_table$Aech, rownames(Aech_exp)),],
                     Lhum_exp[match(gene_ortholog_table$Lhum, rownames(Lhum_exp)),],
                     Mpha_exp[match(gene_ortholog_table$Mpha, rownames(Mpha_exp)),],
                     Sinv_exp[match(gene_ortholog_table$Sinv, rownames(Sinv_exp)),],
                     Lnig_exp[match(gene_ortholog_table$Lnig, rownames(Lnig_exp)),],
                     Cbir_exp[match(gene_ortholog_table$Cbir, rownames(Cbir_exp)),],
                     Dqua_exp[match(gene_ortholog_table$Dqua, rownames(Dqua_exp)),])

#Phenotype information: Gyne, Minor worker, Worker, Reproductive worker, and Non-reproductive worker
caste = factor(c(rep(c('Gyne','Minor_worker'),4),
                 rep(c('Gyne','Worker'),4),
                 rep(c('Gyne','Worker'),5),
                 rep(c('Gyne','Worker'),4),
                 rep(c('Gyne','Worker'),4),
                 rep(c('Reproductive','Non_reproductive'),4),
                 rep(c('Reproductive','Non_reproductive'), each = 6)))

#Colony information.
colony = factor(c(rep(c(1:21), each = 2),
                  rep(c(22:25), each = 2), rep(c(26:31),2)))

#Species information.
species_info =  factor(c(rep("Aech",8),rep("Lhum",8),rep("Mpha",10),rep("Sinv",8), rep("Lnig",8),
                         rep('Cbir',8),rep('Dqua',12)))
# Lab_information:
lab_info =  factor(c(rep("CSE_dk",26), rep("RBC_tw",8), rep("DEE_CH",8),
                     rep('LSE_us',8),rep('TBI_uk',12)))

#Construct sample information table.
sampleTable <- data.frame(caste = caste, species = species_info, lab = lab_info, colony = colony)
rownames(sampleTable) <- colnames(ortholog_exp)

ant_tree.data = "(((((((Mpha, Sinv), Aech), Lnig), Lhum), Cbir), Dqua));"
ant_tree.info = data.frame(row.names = c('Aech','Mpha','Sinv','Lnig','Lhum','Cbir','Dqua'),
                           type = c(rep("Q+",5),rep('Q-',2)),
                           lab = c("CSE",'CSE','RBC','DEE','CSE','LSE','TBI'))

# Import count data, for differentially expressed gene analysis
Aech_count = read_input('Aech',col_name = c(4,6,7,9,10,12,13,15), input_dir,data_type = 'counts') 
Lhum_count = read_input('Lhum',col_name = c(3:10), input_dir, data_type = 'counts') 
Mpha_count = read_input('Mpha',col_name = c(1:6,'7x',8,'9x',10), input_dir, data_type = 'counts') 
Sinv_count = read_input('Sinv',col_name = c(3,4,7,8,11,12,15,16), input_dir, data_type = 'counts') 
Lnig_count = read_input('Lnig',col_name = c(1:8), input_dir, data_type = 'counts') 
Cbir_count = read_input('Cbir',col_name = c(1:8), input_dir, data_type = 'counts') 
Dqua_count = read_input('Dqua',col_name = c(1:5,7:13), input_dir, data_type = 'counts') 
#####
ortholog_counts = cbind(Aech_count[match(gene_ortholog_table$Aech, rownames(Aech_count)),],
                        Lhum_count[match(gene_ortholog_table$Lhum, rownames(Lhum_count)),],
                        Mpha_count[match(gene_ortholog_table$Mpha, rownames(Mpha_count)),],
                        Sinv_count[match(gene_ortholog_table$Sinv, rownames(Sinv_count)),],
                        Lnig_count[match(gene_ortholog_table$Lnig, rownames(Lnig_count)),],
                        Cbir_count[match(gene_ortholog_table$Cbir, rownames(Cbir_count)),],
                        Dqua_count[match(gene_ortholog_table$Dqua, rownames(Dqua_count)),])


save(sampleTable, ant_tree.data, ant_tree.info, ortholog_exp, ortholog_counts, file = 'inputData.Rdata')
