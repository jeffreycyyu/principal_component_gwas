load('~/scratch/principal_component_GWAS/sladek_gip1/gip1_positions.RData')
load('~/scratch/principal_component_GWAS/sladek_gip1/gip1_positions_filtered.RData')





significant_gip1_positions = gip1_positions[which(gip1_positions$p_val < 0.05/nrow(gip1_positions)),]


significant_gip1_positions = significant_gip1_positions[order(significant_gip1_positions$p_val),]


library(biomaRt)

rsids = significant_gip1_positions$rsid
#Select mart
ensembl <- useEnsembl(biomart = "snp", dataset="hsapiens_snp")

#get genomic position
SNPs <- getBM(attributes=c("refsnp_id",
                           "ensembl_gene_stable_id"),
              filters ="snp_filter", values =rsids, mart = ensembl, uniqueRows=TRUE)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <-  SNPs$ensembl_gene_stable_id
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = genes, mart= mart)
colnames(gene_IDs) = c("ensembl_gene_stable_id", "hgnc_symbol")

gene_IDs_SNPs = merge(SNPs, gene_IDs,by="ensembl_gene_stable_id")
colnames(gene_IDs_SNPs) = c("ensembl_gene_stable_id", "rsid", "hgnc_symbol")

significant_gip1_positions_genes = merge(gene_IDs_SNPs, significant_gip1_positions,by="rsid")


significant_gip1_positions_genes = significant_gip1_positions_genes[order(significant_gip1_positions_genes$p_val),]

significant_gip1_positions_genes = significant_gip1_positions_genes[,c("rsid", "hgnc_symbol", "chr", "start", "p_val", "beta", "se")]


