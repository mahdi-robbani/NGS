library('tximport')
library("RColorBrewer")
library("pheatmap")
library('preprocessCore')
library("FactoMineR")
library('factoextra')
library('tidyverse')
library('phytools')
library('devtools')
library('Biobase')
library('DESeq2')
library('rjson')
library('fastICA')
library('sva')
library('ggtree')
read_input = function(species_header, col_name, input_dir, data_type = 'abundance'){
  quant_files = c(paste(species_header,col_name,sep ='_'))
  files <- file.path(input_dir, "input/deg_salmon_gemoma/", species_header,
                     'quants', quant_files, "quant.sf") # Output of salmon, transcript level quantification
  names(files) <- quant_files
  tx2gene <- read.csv(paste(input_dir, '/input/deg_salmon_gemoma/', species_header,'/',species_header,
                            '_gemoma_t2g.txt',sep = ''),header = F, sep = '\t') # Transcript to gene information, for gene level quantification
  tx2gene$V1 = toupper(tx2gene$V1)
  txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,
                         countsFromAbundance = 'lengthScaledTPM')
  abundance_data = txi.salmon[[data_type]]
  rownames(abundance_data)  = paste(species_header,  rownames(abundance_data), sep = '_')
  return(abundance_data)
}

colors <- rev(colorRampPalette(brewer.pal(9, "Blues"))(255))
colors_gene <- colorRampPalette(brewer.pal(9, "RdYlBu"))(255)

sp_color = rainbow(7)

ann_colors = list(species = c('Aech' = sp_color[1],'Sinv' = sp_color[2],'Mpha' = sp_color[3],'Lnig' = sp_color[4],'Lhum' = sp_color[5],
                              'Dqua' = sp_color[6],'Cbir' = sp_color[7]),
                  caste = c(Gyne = rgb(1,0,0,0.8), Worker = rgb(0,0,1,0.8), Minor_worker = 'green',
                            Reproductive = 'purple', `Non_reproductive` = 'pink'))

