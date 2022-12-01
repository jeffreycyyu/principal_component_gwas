


#===================
# ORDERING SNPS POSITIONS (in file: positions_dataframe)
#===================
load('~/scratch/principal_component_GWAS/sladek_gip1/gip1_dataframe.RData')


print(1)

library(biomaRt)

rsids = gip1_dataframe$rsid
#Select mart
ensembl <- useEnsembl(biomart = "snp", dataset="hsapiens_snp")

print(2)

#get genomic position
SNPs <- getBM(attributes=c("refsnp_id",
                           "chr_name",
                           "chrom_start",
                           "chrom_end"),
              filters ="snp_filter", values =rsids, mart = ensembl, uniqueRows=TRUE)

print(3)

save(rsids, SNPs, file = '~/scratch/principal_component_GWAS/sladek_gip1/gip1_positions.RData')

print(4)
