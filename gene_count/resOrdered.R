##---------------------------
# Organizar Data Frame resOrdered
# Data: 16/12/2019
# Washington Candeia
# GBS (recuperados) x ZIKV
# Tibble com genes organizados 
# do menor ao maior padj
##---------------------------
library(tidyverse)
library(org.Hs.eg.db)

# Arquivo com ENSEMBL IDs.
## 1. Criar coluna SYMBOL contendo símbolos dos genes 
#  associados aos IDs Ensembl.
res <- read_csv('resOrdered/resOrdered_zika_vs_control.csv')

my_data <- as_tibble(res)
my_data

res <- my_data %>% rename(X1 = "ensembl_id")
res

# A. Anotações dos símbolos a partir do Ensembl.
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res$ensembl_id, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")

# B. Criar a coluna de símbolos a partir dos IDs Ensembl, da primeira coluna.
ens2symbol <- as_tibble(ens2symbol)
ens2symbol

# C. Confirmar:
head(ens2symbol, 10)

# D. Unir a coluna SYMBOL ao data frame:
res <- inner_join(res, ens2symbol, by = c('ensembl_id'='ENSEMBL'))

# E. Confirmar:
head(res, 10)

# F. Organizar data frame com genes Ensembl dos menores p-values aos maiores
# Retirar NAs
resOrdered <- res %>% na.omit()
resOrdered <- resOrdered[order(resOrdered$pvalue), ]
head(resOrdered, 12)

# Gravar um data frame .CSV para observar os nomes dos genes:
write.csv(as.data.frame(resOrdered), './genesOrdered_pvalue/zika_vs_control_orderedPvalue.csv')

