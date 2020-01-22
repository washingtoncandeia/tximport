library(biomaRt)

# Criação de IDs Ensembl
mart <- biomaRt::useMart(biomart = "ensembl", 
                         dataset = "hsapiens_gene_ensembl", 
                         host="www.ensembl.org")


t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id_version", 
                                     "ensembl_gene_id_version"), 
                      mart = mart)


t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id_version,
                     ens_gene = ensembl_gene_id_version)

head(t2g, 12)
nrow(t2g)