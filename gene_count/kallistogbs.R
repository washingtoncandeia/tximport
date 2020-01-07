##---------------------------------------------
# Análise em Nível de Transcrito Kallisto
# Utilizando Octuplicatas
# alpha = 0.05 (FDR, padj < 0.05)
# Wald test p-value: condition gbs vs control
# Data: 07/01/2020
##---------------------------------------------
library(tximport)
library(DESeq2)
library(apeglm)
library(biomaRt)
library(org.Hs.eg.db)
library(readr)
library(dplyr)
library(rhdf5)
library(IHW)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(genefilter)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Hsapiens.v79)

## Parte 1 - Preparação de dados das amostras de kallisto.
# Caminho dos arquivos (fele path)
dir <- './results'
list.files(dir)

# Nomes de populações
ZIKA <- 'ZIKA'
CHIKV <- 'CHIKV'
CHIKV_REC <- 'CHIKV_REC'
GBS <- 'GBS'
GBS_REC <- 'GBS_REC'
CONTROL <- 'CONTROL'


# Vetor com nomes populações para fazer coluna pop:
# Obs.: Ao fazer o vetor, melhor deixar organizado combinando amostra + nomearquivo
# Coluna 1
pop <- c(rep(GBS, 8),      # GBS        1
         rep(CONTROL, 8),  # Control   10
         rep(CONTROL, 8),  # Control   11
         rep(ZIKA, 8),     # Zika      12
         rep(GBS, 8),      # GBS       13
         rep(GBS, 8),      # GBS       14
         rep(GBS_REC, 8),  # GBSrec    16 
         rep(CHIKV_REC, 8),# Chikv rec 18
         rep(CONTROL, 8),  # Control   19
         rep(GBS, 8),      # GBS        2
         rep(CONTROL, 8),  # Control   20
         rep(CHIKV, 8),    # Chikv     22
         rep(ZIKA, 8),     # Zika      24
         rep(GBS, 8),      # GBS       25
         rep(GBS_REC, 8),  # GBSrec    27
         rep(GBS_REC, 8),  # GBSrec    28
         rep(GBS_REC, 8),  # GBSrec    30
         rep(CHIKV, 8),    # Chikv     31
         rep(CHIKV, 8),    # Chikv     32
         rep(CONTROL, 8),  # Control   33  
         rep(CONTROL, 8),  # Control   34
         rep(CHIKV_REC, 8),# Chikv rec 35
         rep(ZIKA, 8),     # Zika      36
         rep(GBS_REC, 8),  # GBSrec    38
         rep(GBS_REC, 8),  # GBSrec    39
         rep(GBS, 8),      # GBS        4
         rep(CHIKV, 8),    # Chikv     40
         rep(CHIKV_REC, 8),# Chikv rec 41
         rep(CHIKV_REC, 8),# Chikv rec 42
         rep(CHIKV_REC, 8),# Chikv rec 46
         rep(CHIKV_REC, 8),# Chikv rec 47
         rep(ZIKA, 8),     # Zika      48
         rep(GBS_REC, 8),  # GBSrec     7
         rep(GBS_REC, 8),  # GBSrec     8
         rep(GBS_REC, 8))  # GBSrec     9

head(pop, 12)    # ZIKA, CHIKV, CHIKV_REC, GBS, GBS_REC, CONTROL
length(pop)      # 280


# Coluna 2
# Nome de centro de pesquisa para fazer coluna center
center <- rep('IMT-UFRN', 280)
head(center)
length(center)   # 280

# Coluna 3
# Nomes de amostras analisadas para fazer coluna run
run <- list.files(dir)
head(run)
mode(run)
length(run)      # 280


# Coluna 4
condition <- c(rep('gbs', 8),       # GBS        1
               rep('control', 8),   # Control   10
               rep('control', 8),   # Control   11
               rep('zika', 8),      # Zika      12
               rep('gbs', 8),       # GBS       13
               rep('gbs', 8),       # GBS       14
               rep('gbs_rec', 8),   # GBSrec    16 
               rep('chikv_rec', 8), # Chikv rec 18
               rep('control', 8),   # Control   19
               rep('gbs', 8),       # GBS        2
               rep('control', 8),   # Control   20
               rep('chikv', 8),     # Chikv     22
               rep('zika', 8),      # Zika      24
               rep('gbs', 8),       # GBS       25
               rep('gbs_rec', 8),   # GBSrec    27
               rep('gbs_rec', 8),   # GBSrec    28
               rep('gbs_rec', 8),   # GBSrec    30
               rep('chikv', 8),     # Chikv     31
               rep('chikv', 8),     # Chikv     32
               rep('control', 8),   # Control   33  
               rep('control', 8),   # Control   34
               rep('chikv_rec', 8), # Chikv rec 35
               rep('zika', 8),      # Zika      36
               rep('gbs_rec', 8),   # GBSrec    38
               rep('gbs_rec', 8),   # GBSrec    39
               rep('gbs', 8),       # GBS        4
               rep('chikv', 8),     # Chikv     40
               rep('chikv_rec', 8), # Chikv rec 41
               rep('chikv_rec', 8), # Chikv rec 42
               rep('chikv_rec', 8), # Chikv rec 46
               rep('chikv_rec', 8), # Chikv rec 47
               rep('zika', 8),      # Zika      48
               rep('gbs_rec', 8),   # GBSrec     7
               rep('gbs_rec', 8),   # GBSrec     8
               rep('gbs_rec', 8))   # GBSrec     9

length(condition)


# Coluna 5
# Pares
pares <- c(rep('não', 8),  # GBS        1
           rep('não', 8),  # Control   10
           rep('não', 8),  # Control   11
           rep('não', 8),  # Zika      12
           rep('não', 8),  # GBS       13
           rep('não', 8),  # GBS       14
           rep('par', 8),  # GBSrec    16 
           rep('não', 8),  # Chikv rec 18
           rep('não', 8),  # Control   19
           rep('não', 8),  # GBS        2
           rep('não', 8),  # Control   20
           rep('não', 8),  # Chikv     22
           rep('não', 8),  # Zika      24
           rep('não', 8),  # GBS       25
           rep('não', 8),  # GBSrec    27
           rep('não', 8),  # GBSrec    28
           rep('par', 8),  # GBSrec    30
           rep('não', 8),  # Chikv     31
           rep('não', 8),  # Chikv     32
           rep('não', 8),  # Control   33  
           rep('não', 8),  # Control   34
           rep('não', 8),  # Chikv rec 35
           rep('não', 8),  # Zika      36
           rep('não', 8),  # GBSrec    38
           rep('não', 8),  # GBSrec    39
           rep('não', 8),  # GBS        4
           rep('não', 8),  # Chikv     40
           rep('não', 8),  # Chikv rec 41
           rep('não', 8),  # Chikv rec 42
           rep('não', 8),  # Chikv rec 46
           rep('não', 8),  # Chikv rec 47
           rep('não', 8),  # Zika      48
           rep('par', 8),  # GBSrec     7
           rep('não', 8),  # GBSrec     8
           rep('par', 8))  # GBSrec     9

length(pares)

# Coluna 6
# Tipo GBS: desmielinizante x não desmielinizante
# Pares de condições
diagn <- c(rep('Desmielin', 8),   # GBS        1
           rep('não', 8),         # Control   10
           rep('não', 8),         # Control   11
           rep('não', 8),         # Zika      12
           rep('Desmielin', 8),   # GBS       13
           rep('Desmielin', 8),   # GBS       14
           rep('Desmielin', 8),   # GBSrec    16 
           rep('não', 8),         # Chikv rec 18
           rep('não', 8),         # Control   19
           rep('Desmielin', 8),   # GBS        2
           rep('não', 8),         # Control   20
           rep('não', 8),         # Chikv     22
           rep('não', 8),         # Zika      24
           rep('Desmielin', 8),   # GBS       25
           rep('Desmielin', 8),   # GBSrec    27
           rep('Desmielin', 8),   # GBSrec    28
           rep('Inconclus', 8),   # GBSrec    80
           rep('não', 8),         # Chikv     31
           rep('não', 8),         # Chikv     32
           rep('não', 8),         # Control   33  
           rep('não', 8),         # Control   34
           rep('não', 8),         # Chikv rec 35
           rep('não', 8),         # Zika      36
           rep('Desmielin', 8),   # GBSrec    38
           rep('Axonal', 8),      # GBSrec    39
           rep('Desmielin', 8),   # GBS        4
           rep('não', 8),         # Chikv     40
           rep('não', 8),         # Chikv rec 41
           rep('não', 8),         # Chikv rec 42
           rep('não', 8),         # Chikv rec 46
           rep('não', 8),         # Chikv rec 47
           rep('não', 8),         # Zika      48
           rep('Desmielin', 8),   # GBSrec     7
           rep('Desmielin', 8),   # GBSrec     8
           rep('Desmielin', 8))   # GBSrec     9


length(diagn)

# Coluna 7
# Sexo
sex <- c(rep('Fem', 8),  # GBS        1
         rep('Fem', 8),  # Control   10
         rep('Mas', 8),  # Control   11
         rep('Fem', 8),  # Zika      12
         rep('Mas', 8),  # GBS       13
         rep('Mas', 8),  # GBS       14
         rep('Mas', 8),  # GBSrec    16 
         rep('Fem', 8),  # Chikv rec 18
         rep('Mas', 8),  # Control   19
         rep('Fem', 8),  # GBS        2
         rep('Mas', 8),  # Control   20
         rep('Fem', 8),  # Chikv     22
         rep('Fem', 8),  # Zika      24
         rep('Fem', 8),  # GBS       25
         rep('Fem', 8),  # GBSrec    27
         rep('Mas', 8),  # GBSrec    28
         rep('Fem', 8),  # GBSrec    30
         rep('Mas', 8),  # Chikv     31
         rep('Fem', 8),  # Chikv     32
         rep('Fem', 8),  # Control   33  
         rep('Mas', 8),  # Control   34
         rep('Fem', 8),  # Chikv rec 35
         rep('Fem', 8),  # Zika      36
         rep('Mas', 8),  # GBSrec    38
         rep('Fem', 8),  # GBSrec    39
         rep('Fem', 8),  # GBS        4
         rep('Mas', 8),  # Chikv     40
         rep('Fem', 8),  # Chikv rec 41
         rep('Fem', 8),  # Chikv rec 42
         rep('Fem', 8),  # Chikv rec 46
         rep('Fem', 8),  # Chikv rec 47
         rep('Fem', 8),  # Zika      48
         rep('Fem', 8),  # GBSrec     7
         rep('Mas', 8),  # GBSrec     8
         rep('Fem', 8))  # GBSrec     9

length(sex)


# Coluna 8
# Replicatas
replicates <- c('rep01', 'rep02', 'rep03 ', 'rep04',
                'rep05', 'rep06 ', 'rep07', 'rep08')

## Parte II
# Aqui inicia-se a construção do data frame com todas as informações de amostras.
## Parte 2 - Unir cada vetor formando colunas de um data frame:
samples_info <- data.frame(pop = pop,
                           center = center,
                           run = run,
                           condition = condition,
                           diagnostico = diagn,
                           pair = pares,
                           sex = sex,
                           replicate = rep(replicates, 35))   

# Obs.: Replicate: Não colocar 105 pois são de 3 em 3.
# Logo: 35 x 3 = 105. Repete-se a tríade 35 vezes, o que geram 105 replicatas.


# Observações gerais:
head(samples_info, 10)
str(samples_info)
names(samples_info)


## -------------------- Eliminando Linhas do Data Frame por Nomes  -------------------- ##
## Quando necessário o uso de um data frame menor, com apenas algumas variáveis.

# Usando dplyr, função filter e negando com regex (função grepl)

samples_info <- samples_info %>% 
  dplyr::filter(!grepl(CHIKV_REC, pop))

samples_info <- samples_info %>% 
  dplyr::filter(!grepl(CHIKV, pop))

samples_info <- samples_info %>% 
  dplyr::filter(!grepl(ZIKA, pop))

samples_info

## ------------------------------------------------------------------------------------ ##
length(samples_info$condition)
# Salvar a tabela no formato .txt (tsv)
write.table(samples_info, './tables/gbs/condition_gbs_vs_control.txt', sep = '\t')

# Criar um vetor nomeado apontando os arquivos de quantificação.
# Estes arquivos têm seus nomes anotados em uma tabela (samples.txt).
samples <- read.table('./tables/gbs/condition_gbs_vs_control.txt', header = TRUE, row.names = 1)
head(samples, 9)
samples$condition 

# Nomeando as linhas com nome de cada arquivo de amostra:
rownames(samples) <- samples$run

samples

# Relevel: ajustando a condição referência para análise
samples$condition <- relevel(samples$condition, ref = 'control')

# Obtendo cada arquivo de replicata das amostras usadas em kallisto:
files <- file.path(dir, samples$run, 'abundance.h5')
files
names(files) <- samples$run
files
length(files)

# Usando EnsDb.Hsapiens.v86 ou 79 / ensembldb
# EnsDB.Hsapiens.v86
# EX: https://support.bioconductor.org/p/81012/
edb <- EnsDb.Hsapiens.v86
edb79 <- EnsDb.Hsapiens.v79

edb
edb79

txdf <- transcripts(edb, 
                    columns = c(listColumns(edb , "tx"), "gene_name"), 
                    return.type="DataFrame")

tx2gene <- as.data.frame(txdf[,c("tx_id","gene_id")])
head(tx2gene, 10)
nrow(tx2gene)

write.csv(as.data.frame(tx2gene), file = './transcritos/zika/ensembl.csv')

# Versão EnsDb.Hsapiens.v79:
txdf79 <- transcripts(EnsDb.Hsapiens.v79, return.type="DataFrame")

tx2gene79 <- as.data.frame(txdf[,c("tx_id", "gene_id")])
head(tx2gene79, 10)
nrow(txdf79)

## Usando biomaRt para nomear transcritos e genes
mart <- biomaRt::useMart(biomart = "ensembl", 
                         dataset = "hsapiens_gene_ensembl", 
                         host="www.ensembl.org")


t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", 
                                     "ensembl_gene_id"), 
                      mart = mart)


t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id)


head(t2g, 9)
nrow(t2g)

## Obs.: A coluna de transcripts IDs não possui versão.
## Ao utilizar o tximport prestar atenção na opção ignoreTxVersion.

### Parte III - Tximport Utlizando Arquivos kallisto
## Quantificação de Abundâncias de Transcritos com kallisto
## Análise de Expressão Diferencial (DE) com DESeq2.
# Estimativa de contagens a partir de kallisto,
# Usar ignoreTxVersion e ignoreAfterBar para que o data frame de IDs de transcritos
# e genes do Ensembl tenham ignorados as versões e barras |, respectivamente.
txi.kallisto <- tximport(files, 
                         type = 'kallisto',
                         tx2gene = t2g,
                         ignoreTxVersion = TRUE)

# ignoreTxVersion = TRUE deve ser utilizado quando em gene counts analysis.
# Observações gerais
names(txi.kallisto)
head(txi.kallisto$abundance)
head(txi.kallisto$counts)
head(txi.kallisto$length)
head(txi.kallisto$infReps)

# Salvamento de um objeto R para uso posterior:
#dir.create(path = "./count_estimates/")
#save(txi.kallisto, file = "./count_estimates/gbs/gbs.Rdata")


#### Parte IV - DESeq2
## Design com formula simples:
dds.txi <- DESeqDataSetFromTximport(txi = txi.kallisto,
                                    colData = samples,
                                    design = ~condition)

# Relevel: ajustando a condição referência para análise
dds.txi$condition <- relevel(dds.txi$condition, ref = 'control')
str(dds.txi$condition)

# Observar as amostras:
head(as.data.frame(colData(dds.txi)), 12)


# Agora, o objeto dds.Txi pode ser usado como aquele dds nos
# passos subsequentes de DESeq2.
head(dds.txi$condition)

## Pre-filtering
# Filtrar por counts insignificantes.
# Remover genes sem contagens significantes para as amostras. Neste caso, no mínimo 10 contagens.
keep <- rowSums(counts(dds.txi)) >= 10

# Renomear dds.txi para dds:
dds <- dds.txi[keep, ]

head(dds, 9)
# Outra forma:
#dds <- dds.txi[rowSums(counts(dds.txi)) >= 10, ]

# Observar o design
design(dds)


###----------------- Análise de Expressão Diferencial (DE) ----------------- ###
# Objeto dds por DESeq2
# Modelando as contagens com efeitos de condition vs control
dds <- DESeq(dds)

# A função results gera tabelas de resultados.
# condition chikv rec vs control 
res <- results(dds, contrast = c('condition', 'gbs', 'control'), alpha = 0.05)

# Visualizar
res
# Summary
summary(res)

# Informações de metadados do objeto res: colunas.
mcols(res, use.names=TRUE)

# Salvar .csv Wald test p-value: condition gbs/gbs_rec vs control em .csv para fgsea
write.csv(as.data.frame(res), file = './GSEA/chikv/fgsea_condition_gbs_vs_control.csv')

## Reordenando Resultados com p-values e adjusted p-values
# Ordenar os resultados da tabela por menor p value:
resOrdered <- res[order(res$pvalue), ]
resOrdered
write.csv(as.data.frame(resOrdered), file = './tables/gbs/resOrdered/resOrdered_gbs_vs_control.csv')

# Examinando as contagens e contagens normalizadas para os genes com menor p-value:
idx <- which.min(res$pvalue)
idx
counts(dds)[idx, ]

#Normalization
norm.Counts <- counts(dds, normalized=TRUE)[ idx, ]
norm.Counts

## Continuação:
## Log fold change shrinkage for visualization and ranking
# Contração log fod change para visualização e ranqueamento.
# Shrinkage of effect size (LFC estimates)
resultsNames(dds)

# coef = 2 ("condition_chikv_vs_control")
# coef = 3 ("condition_chikv_rec_vs_control") 

### LFC - Log2 Fold Changes Shrinkage
## Visualizando para log2 fold changes shrinkage, LFC (Shrinkage of Effect Size)
# Para contrair (shrink) LFC passar objeto dds para função lfcShrink:
resLFC <- lfcShrink(dds, coef = 3, type = 'apeglm', res = res)

# Observar
resLFC
# Summary
summary(resLFC)

# FDR cutoff, alpha = 0.05.
res001 <- results(dds, alpha=0.01)
summary(res025)

### Independent Hypothesis Weighting
## Ponderação de Hipóteses Independentes
# Filtragem de p value: ponderar (weight) hipóteses para otimizar o poder.
# Está disponível no Bioconductor sob nome IHW.
resIHW <- results(dds, filterFun = ihw, alpha = 0.05)
resIHW
# Summary
summary(res)
# Summary
summary(res025)
# Summary
summary(resIHW)
# Summary
summary(resLFC)

## Continuação:
## Log fold change shrinkage for visualization and ranking
# Contração log fod change para visualização e ranqueamento.
# Shrinkage of effect size (LFC estimates)
resultsNames(dds)

## coeficientes nesta análise: 
# coef = 1 ("condition_gbs_vs_control")     
# coef = 2 ("condition_gbs_rec_vs_control")

##### Parte V - Exploração de Resultados
## MA-plot

# A função plotMA mostra os log2 fold change atribuível a uma dada variável
# sobre a média de contagens normalizadas para todas as amostras no DESeqDataSet.
plotMA(res , ylim = c(-2, 2))

# Objeto com alpha < 0.05 (adjusted p-value < 0.1)
plotMA(res025, ylim = c(-2, 2))

# Objeto resLFC
plotMA(resLFC, ylim = c(-2, 2))

# Objeto resIHW
plotMA(resIHW, ylim = c(-2, 2))

# Agora, observar os plots juntos
# Versão 1: ma_all_v1.jpeg
par(mfrow=c(2,3), mar=c(4,4,4,2))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(res, xlim=xlim, ylim=ylim, main="Síndrome de Guillain-Barré vs Controles \n(pdaj < 0.05)")
plotMA(resLFC, xlim=xlim, ylim=ylim, main="Síndrome de Guillain-Barré vs Controles \n(LFC)")
plotMA(resIHW, xlim=xlim, ylim=ylim, main="Síndrome de Guillain-Barré vs Controles \n(IHW)")
plotMA(res025, xlim=xlim, ylim=ylim, main="Síndrome de Guillain-Barré vs Controles \n(padj < 0.01)")
# Pontos em vermelho: se o adjusted p value for menor que 0.1.

par(mfrow=c(2,2), mar=c(4,4,4,2))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(res, xlim=xlim, ylim=ylim, main="Síndrome de Guillain-Barré vs Controles \n(pdaj < 0.05)")
plotMA(resLFC, xlim=xlim, ylim=ylim, main="Síndrome de Guillain-Barré vs Controles \n(LFC)")
plotMA(res025, xlim=xlim, ylim=ylim, main="Síndrome de Guillain-Barré vs Controles \n(padj < 0.01)")
plotMA(resIHW, xlim=xlim, ylim=ylim, main="Síndrome de Guillain-Barré vs Controles \n(IHW)")

# Pontos em vermelho: se o adjusted p value for menor que 0.1.

# Versão 2
par(mfrow=c(2,2), mar=c(4,4,4,2))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(res, xlim=xlim, ylim=ylim, main="Síndrome de Guillain-Barré vs Control \n(pdaj < 0.1)")
plotMA(resOrdered, xlim=xlim, ylim=ylim, main="Síndrome deGuillain-Barré vs Control \n(Reordenado por menor pvalue)")
plotMA(resLFC, xlim=xlim, ylim=ylim, main="Síndrome de Guillain-Barré vs Control \n(LFC)")
plotMA(res05, xlim=xlim, ylim=ylim, main="Síndrome de Guillain-Barré vs Control \n(padj < 0.01)")
plotMA(resIHW, xlim=xlim, ylim=ylim, main="Síndrome de Guillain-Barré vs Control \n(IHW)")

# Pontos em vermelho: se o adjusted p value for menor que 0.1.

### Alternative Shrinkage Estimators
## Lembrar de objeto LFC e Shrinkage

# Especificar o coeficiente pela ordem em que aparece em results(dds)
# O coeficiente usado em lfcShrink anterior (resNorm) foi "condition chikv vs control"
# Porém, é possível especificar o coeficiente pela ordem em que aparece quando se usa resultsnames(dds):

# Em argumento type, pode-se utilar os parâmetros: ashr, apeglm e normal.
# a. ashr - adaptive shrinkage estimator from the ashr package (Stephens 2016)
# b. apeglm - the adaptive t prior shrinkage estimator from the apeglm package (Zhu, Ibrahim, and Love 2018).
# c. normal - estimador shrinkage original de DESeq2 (an adaptive Normal distribution as prior).

# Usaremos o coeficiente como 2, pois é o que indica condition_chikv_vs_control.
resultsNames(dds)

resNorm <- lfcShrink(dds, coef="condition_gbs_vs_control", type="normal")
resAsh <- lfcShrink(dds, coef="condition_gbs_vs_control", type="ashr")
resLFC <- lfcShrink(dds, coef="condition_gbs_vs_control", type='apeglm')

## Agora, observar os plots juntos para coef = 3
par(mfrow=c(1,3), mar=c(4,4,4,2))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="Síndrome de Guillain-Barré vs Controles \nLFC (apeglm)")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="Síndrome de Guillain-Barré vs Controles \n(Normalizados)")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="Síndrome de Guillain-Barré vs Controles \n(Adaptative Shrinkage Estimator)")


###### Parte VI
## Transformação e Visualização de Dados
# Transformações de contagem

# VST - Variance Stabilizing Transformation 
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)
# Os valores transformados não são contagens, sendo armazenados no slot assay.
# colData está ligado a dds e é acessível:
colData(vsd)

# RLD - Regularized log Transformation 
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 6)
# Os valores transformados não são contagens, sendo armazenados no slot assay.
# colData está ligado a dds e é acessível:
colData(rld)

## Efeitos das Transformações na Variância
# Plot de desvio padrão dos dados transformados através das amostras,
# contra a média, usando shifting logarithm transformation
# fornece log2(n + 1)
ntd <- normTransform(dds)

## library(vsn)
# 1. Objeto ntd
meanSdPlot(assay(ntd))

# 2. Objeto vsd
meanSdPlot(assay(vsd))

# 3. Objeto rld
meanSdPlot(assay(rld))

### PCA - Principal component plot das amostras
# O plot PCA está relacionado à matriz de distância e evidencia as amostras no plano de 2D
# abrangidas por seus primeiros componentes principais.
# Esse gráfico é útil para visualizar o efeito geral de covariantes experimentais (e batch effects).

## Usando VST
# vsd object
plotPCA(vsd, intgroup=c("condition", "run"))

# vsd object
plotPCA(vsd, intgroup=c("condition", "replicate"))

## Usando RLT
# rld object
plotPCA(rld, intgroup=c("condition", "run"))

# rld
plotPCA(rld, intgroup=c("condition", "replicate"))

plotPCA(rld, 
        intgroup=c("condition", "replicate"),
        returnData = TRUE)


# Outras formas:
pcaVSD <- plotPCA(vsd, 
                  ntop = nrow(counts(dds)),
                  returnData=FALSE)

pcaVSD2 <- plotPCA(vsd, 
                   ntop = nrow(counts(dds)),
                   intgroup=c("condition", "replicate"),
                   returnData=FALSE)

pcaVSD
pcaVSD2

pcaRLD <- plotPCA(rld, 
                  ntop = nrow(counts(dds)),
                  returnData=FALSE)

pcaRLD2 <- plotPCA(rld, 
                   ntop = nrow(counts(dds)),
                   intgroup=c("condition", "replicate"),
                   returnData=FALSE)

pcaRLD
pcaRLD2

# PCA sob TRUE
plotPCA(rld, 
        ntop = nrow(counts(dds)),
        intgroup=c("condition", "replicate"),
        returnData=TRUE)


## Utilizando ggplot2 com os dados
## ggplot2
# rsd object
pcaData <- plotPCA(vsd, intgroup=c("condition", "replicate"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# Outra forma , de acordo com pcaRLD, pcaVSD:
pcaData <- plotPCA(vsd,  ntop = nrow(counts(dds)), intgroup=c("condition", "replicate"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()



## Qualidade de Dados por Clusterização e Visualização
# Heatmap da matriz de contagem 
select <- order(rowMeans(counts(dds,normalized = TRUE)),
                decreasing = TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","replicate")])  # Todas as linhas (genes) e variáveis (colunas condition e replicate)
pheatmap(assay(ntd)[select,], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = FALSE, annotation_col = df)


# Utilizando variável condition e pop
df2 <- as.data.frame(colData(dds)[,c("condition","pop")])
pheatmap(assay(ntd)[select,], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = FALSE, annotation_col = df2)


## VST - Variance Stabilizing Transformation 
# vsd - condition, replicate (df)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# vsd - condition, pop (df2)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df2)

## RLD - Regularized log Transformation 
# rld - condition, replicate (df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# rld - condition, pop (df2)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df2)


### Heatmap das Distâncias Amostra-Amostra
# O utro uso de dados transformados: sample clustering.
# Usando a função dist para transposição de matriz de contagem transformada.
sampleDists <- dist(t(assay(vsd)))

# library(RColorBrewer)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$run, sep=" - ")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


## Plot counts
# É útil examinar a contagem de reads para um único gene entre os grupos (control e zika).
# Existe a função plotCounts que pode fazer isso, a qual normaliza as contagens por profundidade
# de sequenciamento (sequencing depth) e adiciona uma pseudocontagem de 1/2 para permitir a plotagem
# em escala de log.
# Pode-se selecionar o gene de interesse a ser plotado por rowname ou por índice numérico.
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

# Customização com ggplot2
a <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(a, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


# Outros genes: adicionar on ID (nome):
plotCounts(dds, gene='ENSG00000135845', intgroup="condition")

# Outros genes
b <- plotCounts(dds, gene='ENSG00000135845', intgroup="condition", returnData = T)
library("ggplot2")
ggplot(b, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

# Outros genes
c <- plotCounts(dds, gene='ENSG00000168300', intgroup="condition", returnData = T)
library("ggplot2")
ggplot(c, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

d <- plotCounts(dds, gene='ENSG00000254708', intgroup="condition", returnData = T)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


## Reordenando Resultados com p-values e adjusted p-values
# Ordenar os resultados da tabela por menor p value:
resOrdered <- res[order(res$pvalue), ]

# Criar csv (pode ser usado em fgsea)
write.csv(as.data.frame(resOrdered), file = './GSEA/gbs/condition_gbs_rec_vs_control_resOrdered_GSEA.csv')
# Outros genes: adicionar on ID (nome): TYMP
# Gene 2 da lista organizada por pvalue (resOrdered)
plotCounts(dds, gene='ENSG00000025708', intgroup="condition")
# Gene3
plotCounts(dds, gene='ENSG00000135845', intgroup="condition")
# Gene 4
plotCounts(dds, gene='ENSG00000104765', intgroup="condition")
# Gene 5
plotCounts(dds, gene='ENSG00000153574', intgroup="condition")
# Gene 6
plotCounts(dds, gene='ENSG00000137267', intgroup="condition")
# Gene 7
plotCounts(dds, gene='ENSG00000198692', intgroup="condition")
# Gene 8
plotCounts(dds, gene='ENSG00000141574', intgroup="condition")
# Gene 9
plotCounts(dds, gene='ENSG00000013561', intgroup="condition")
# Gene 10
plotCounts(dds, gene='ENSG00000155749', intgroup="condition")


## Mais informações na coluna Results
mcols(res)$description

### Lembrar dos objetos:
# res
# res05 (padj < 0.05)
# resOrdered
# resLFC
# resNorm
# resAsh

########################################### SELEÇÃO GENES UP e DOWN ###########################################
contr_gbs <- as.data.frame(res)

# Criar uma nova coluna com os nomes (SYMBOLS) dos genes.
contr_gbs$genes <- rownames(contr_gbs)

# Remoção de NAs na coluna de padj.
contr_gbs$padj[is.na(contr_gbs$padj)] <- 1
DEG_gbs <- subset(contr_gbs, padj <= 0.05 & abs(log2FoldChange) > 1)
head(DEG_gbs, 9)
nrow(DEG_gbs)

# Seleção dos genes DEs
gbs.DE <- subset(DEG_gbs, padj <= 0.05 & abs(log2FoldChange) > 0.5)
nrow(gbs.DE)
head(gbs.DE, 15)
# Abaixo um arquivo que pode servir como teste para GSEA/fgsea:
write.csv(as.data.frame(gbs.DE), file = './DE/gbs/gbs_DE_subset_gbs_vs_control.csv')


# Filtrando por log2FC
DE_subset <- subset(gbs.DE, log2FoldChange >= 1.5 | log2FoldChange <= -1.5)
head(DE_subset, 9)
nrow(DE_subset)

# Filtrando por padj
DE_subset <- subset(DE_subset, padj <= 0.01)  # ou padj <= 0.05
nrow(DE_subset)

# Pode ser usado em fgsea
write.csv(DE_subset, file = './DE/gbs/DE_subset_gbs_vs_control.csv')

# Separando em genes up e down regulados
DE_down <- subset(DE_subset, log2FoldChange <= -1.5)  # Downregulated genes
length(DE_down$genes)
head(DE_down, 16)

# Tabela de contagem de genes 'downregulados'
#write.table(as.data.frame(DE_down), sep = ",", "./DE/gbs/transcript_count_gbs_rec_DOWN.txt", row.names = FALSE)
write.csv(as.data.frame(DE_down), file = './DE/gbs/transcript_count_gbs_DOWN.csv')

DE_up <- subset(DE_subset, log2FoldChange >= 1.5)     # Upregulated genes
length(DE_up$genes)
head(DE_up, 16)

# Tabela de contagem de genes 'upregulados'
#write.table(as.data.frame(DE_up), sep = ",", "./DE/gbs/transcript_count_gbs_rec_UP.txt", row.names = FALSE)
write.csv(as.data.frame(DE_up), file = './DE/gbs/transcript_count_gbs_rec_UP.csv')

###################################################################################################################