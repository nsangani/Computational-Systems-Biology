data = read.csv(file.choose(), header = T)

library(biomaRt)
mart.snp <- useMart("ENSEMBL_MART_SNP", "hsapiens_snp")

getENSG <- function(rs = "rs3043732", mart = mart.snp) {
  results <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id"),
                   filters    = "snp_filter", values = rs, mart = mart)
  return(results)
}

getENSG(rs = data[,1])
