library(synapseClient)
library("dplyr")
library("org.Hs.eg.db")

options(stringsAsFactors = F)

synapseLogin()


#heatmap code
source("heatmap_code.R")

#get the MsigDB object
cat('Reading the MSIGDB object from synapse...')
MSIGDB_syn<-synGet("syn2227979")
load(MSIGDB_syn@filePath) #available as MSigDB R object
pathways_list <- c(MSigDB$C2.CP.BIOCARTA, MSigDB$C2.CP.KEGG, MSigDB$C2.CP.REACTOME)
cat('..Done\n\n')

#gene annotation
cat('Preparing the hg19 annotation df.....')
k <- keys(org.Hs.eg.db,keytype="SYMBOL")
hg19_gene_annot <- select(org.Hs.eg.db, keys=k, columns=c("GENENAME","ALIAS", "ENSEMBL", "ENSEMBLTRANS", "ENTREZID"), keytype="SYMBOL")
#saveRDS(hg19_gene_annot,"precomputed_hg19_gene_annot.RDS")
hg19_grpd <- hg19_gene_annot %>%
                group_by(ENSEMBL) %>%
                summarise(ALIAS = paste(unique(ALIAS),collapse=", "),
                          SYMBOL = paste(unique(SYMBOL),collapse=", "),
                          GENENAME = paste(unique(GENENAME),collapse=", "),
                          ENTREZID = paste(unique(ENTREZID),collapse=", ")
                          )
hg19_grpd <- as.data.frame(hg19_grpd)
cat('Done \n\n')



#get the expression matrix
cat('Reading expression data from synapse...')
expMat_synId <- 'syn2582946'
expMat <- synGet(expMat_synId)
expMat <- read.table(expMat@filePath, header=T)
rownames(expMat) <- expMat$Gene_ID
expMat$Gene_ID <- NULL
#log2 transform
expMat <- log2(expMat + 1)
boxplot(expMat)
cat('Done \n\n')


#get the enriched geneList
cat('Reading the enriched geneList from synapse...')
enriched_geneLists_synId <- "syn2588320"
enriched_geneLists <- synGet(enriched_geneLists_synId)
enriched_geneLists <- read.table(enriched_geneLists@filePath, skip=1, sep="\t")
colnames(enriched_geneLists) <- c('symbol', 'entrezId', 'ensemblId', 'source', 'url', 
                                  'name', 'logFold', 'BH-adjPval')
enriched_geneLists['geneListName'] <- gsub('-fold2.0_adjp0.05', '', enriched_geneLists$name)
geneLists <- unique(enriched_geneLists['geneListName'])
cat("Done\n\n")



#get the metadata
cat('Reading the metadata from synapse...')
metadata_synId <- "syn2583490"
metadata <- synGet(metadata_synId)
metadata <- read.table(metadata@filePath, header=T, sep="\t")
cat("Done\n\n")


#preselected_gene_list
cat('Reading a sample genelist...')
df <- read.table("precomputed_data/pre_selected_genelist.txt",sep="\t")
sample_gene_list <- paste(as.character(unique(df$V5), collapse=", "), collapse=", ")
cat('Done \n\n')