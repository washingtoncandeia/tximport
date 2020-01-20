##---------------------------------------------
# Análise em Nível de Transcrito Kallisto
# Utilizando Octuplicatas
# alpha = 0.05 (FDR, padj < 0.05)
# Wald test p-value: condition gbs_rec vs control
# Data: 19/01/2020
##---------------------------------------------
library(tximport)
library(DESeq2)
library(apeglm)
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
library(ggplot2)
library(sva)

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

#
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

# Relevel: ajustando a condição referência para análise
samples$condition <- relevel(samples$condition, ref = 'control')


# Obtendo cada arquivo de replicata das amostras usadas em kallisto:
files <- file.path(dir, samples$run, 'abundance.h5')
files
names(files) <- samples$run
head(files, 9)

# Transcript to Gene (SYMBOLS)
# Criação de IDs Ensembl
# Arquivo t2g data frame
t2g <- read_csv('./t2sg.csv')
head(t2g, 12)
nrow(t2g)

head(t2g, 12)
nrow(t2g)
### Parte III - Tximport Utlizando Arquivos kallisto
txi.kallisto <- tximport(files, 
                         type = 'kallisto',
                         tx2gene = t2g)

# Observações gerais
names(txi.kallisto)
head(txi.kallisto$abundance)

#### Parte IV - DESeq2
## Design com formula simples:
dds.txi <- DESeqDataSetFromTximport(txi = txi.kallisto,
                                    colData = samples,
                                    design = ~condition)

head(dds.txi$condition)

## Pre-filtering
# Remover genes sem contagens significantes para as amostras. Neste caso, no mínimo 10 contagens.
keep <- rowSums(counts(dds.txi)) >= 10
dds <- dds.txi[keep, ]
head(dds$condition, 9)

### Análise de Expressão Diferencial (DE)
# Objeto dds por DESeq2
# Modelando as contagens com efeitos de condition vs control
dds <- DESeq(dds)
res <- results(dds, contrast = c('condition', 'gbs', 'control'), filterFun=ihw, alpha = 0.05)
res
summary(res)

write.csv(as.data.frame(res), file = './GSEA/fgsea_gbs_vs_control_SYMBOL.csv')

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

# Filtrando por log2FC
DE_subset <- subset(gbs.DE, log2FoldChange >= 1.5 | log2FoldChange <= -1.5)
head(DE_subset, 9)
nrow(DE_subset)

# Filtrando por padj
DE_subset <- subset(DE_subset, padj <= 0.05)  # ou padj <= 0.05
nrow(DE_subset)

# Separando em genes up e down regulados
DE_down <- subset(DE_subset, log2FoldChange <= -1.5)  # Downregulated genes
length(DE_down$genes)
head(DE_down, 16)

# Tabela de contagem de genes 'downregulados'
#write.table(as.data.frame(DE_down), sep = ",", "./DE/gbs/transcript_count_gbs_rec_DOWN.txt", row.names = FALSE)
write.csv(as.data.frame(DE_down), file = './genes_DE/gbs/transcript_count_gbs_DOWN.csv')

DE_up <- subset(DE_subset, log2FoldChange >= 1.5)     # Upregulated genes
length(DE_up$genes)
head(DE_up, 16)

# Tabela de contagem de genes 'upregulados'
#write.table(as.data.frame(DE_up), sep = ",", "./DE/gbs/transcript_count_gbs_rec_UP.txt", row.names = FALSE)
write.csv(as.data.frame(DE_up), file = './genes_DE/gbs/transcript_count_gbs_UP.csv')

###################################################################################################################

# Dispersão dos dados
#http://www.sthda.com/english/wiki/rna-seq-differential-expression-work-flow-using-deseq2
plotDispEsts(dds, ylim = c(1e-6, 1e1))

# Histogram of pvalue
#http://www.sthda.com/english/wiki/rna-seq-differential-expression-work-flow-using-deseq2
hist(res$pvalue, breaks=20, col="grey")

##### Parte V - Exploração de Resultados
## MA-plot

# LFC - Log2 Fold Changes Shrinkage

# A função plotMA mostra os log2 fold change atribuível a uma dada variável
# sobre a média de contagens normalizadas para todas as amostras no DESeqDataSet.
plotMA(res , ylim = c(-2, 2))

resLFC <- lfcShrink(dds, coef = 'condition_gbs_vs_control', type = 'apeglm', res = res)
res001 <- results(dds, filterFun=ihw, alpha=0.01)

# Agora, observar os plots juntos
# Versão 1: ma_all_v1.jpeg
par(mfrow=c(2,3), mar=c(4,4,4,2))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(res, xlim=xlim, ylim=ylim, main="Guillain-Barré Doentes vs Controles \n(pdaj < 0.05 + IHW)")
plotMA(resLFC, xlim=xlim, ylim=ylim, main="Guillain-Barré Doentes vs Controles \n(LFC)")
plotMA(res001, xlim=xlim, ylim=ylim, main="Guillain-Barré Doentes vs Controles \n(padj < 0.01)")
# Pontos em vermelho: se o adjusted p value for menor que 0.1.

### Alternative Shrinkage Estimators
# Em argumento type, pode-se utilar os parâmetros: ashr, apeglm e normal.
# a. ashr - adaptive shrinkage estimator from the ashr package (Stephens 2016)
# b. apeglm - the adaptive t prior shrinkage estimator from the apeglm package (Zhu, Ibrahim, and Love 2018).
# c. normal - estimador shrinkage original de DESeq2 (an adaptive Normal distribution as prior).

# Usaremos o coeficiente como 2, pois é o que indica condition_zika_vs_control.
resultsNames(dds)
resNorm <- lfcShrink(dds, coef="condition_gbs_vs_control", type="normal")
resAsh <- lfcShrink(dds, coef="condition_gbs_vs_control", type="ashr")
resLFC <- lfcShrink(dds, coef="condition_gbs_vs_control", type='apeglm')

## Agora, observar os plots juntos para coef = 2
par(mfrow=c(1,3), mar=c(4,4,4,2))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="Guillain-Barré Doentes vs Control \nLFC (apeglm)")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="Guillain-Barré Doentes vs Control \n(Normalizados)")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="Guillain-Barré Doentes vs Control \n(Adaptative Shrinkage Estimator-ASHR)")

###### Parte VI
## Transformação e Visualização de Dados
# VST - Variance Stabilizing Transformation 
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)

# RLD - Regularized log Transformation 
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 6)
# Os valores transformados não são contagens, sendo armazenados no slot assay.

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
pcaVSD <- plotPCA(vsd, 
                  ntop = nrow(counts(dds)),
                  returnData=FALSE)

pcaVSD

pcaRLD <- plotPCA(rld, 
                  ntop = nrow(counts(dds)),
                  returnData=FALSE)

pcaRLD



## Qualidade de Dados por Clusterização e Visualização
# Heatmap da matriz de contagem 
select <- pheatorder(rowMeans(counts(dds,normalized = TRUE)),
                     decreasing = TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","replicate")])
pheatmap(assay(ntd)[select,], cluster_rows = TRUE, show_rownames = TRUE,
         cluster_cols = TRUE, annotation_col = df)


# Utilizando variável condition e pop
df2 <- as.data.frame(colData(dds)[,c("condition","pop")])
pheatmap(assay(ntd)[select,], cluster_rows = T, show_rownames = T,
         cluster_cols = T, annotation_col = df2)


## VST - Variance Stabilizing Transformation 
# vsd - condition, replicate (df)
select <- pheatorder(rowMeans(counts(dds,normalized = TRUE)),
                     decreasing = TRUE)[1:20]
df.vsd <- as.data.frame(colData(dds)[,c("condition","replicate")])
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df.vsd)


## RLD - Regularized log Transformation 
# rld - condition, replicate (df)
pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)

# rld - condition, pop (df2)
pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df2)


## Os 20 genes Com Mais Alta Variância Através das Amostras
# Um subconjunto dos genes mais altamente veriáveis:
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)

mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("condition", "replicate")])
pheatmap(mat, annotation_col = anno)
# Heatmap of relative VST-transformed values across samples.

# RLD
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)

mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("condition", "pop")])
pheatmap(mat, annotation_col = anno)

# NTD
topVarGenes <- head(order(rowVars(assay(ntd)), decreasing = TRUE), 20)

mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(ntd)[, c("condition", "pop")])
pheatmap(mat, annotation_col = anno)

# Sample Dists-------------------------------------------------------------------------------------
# http://www.sthda.com/english/wiki/rna-seq-differential-expression-work-flow-using-deseq2

# SampleDists
sampleDists <- dist(t(assay(rld)))
as.matrix(sampleDists)[ 1:3, 1:3 ]
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$run, sep="-")
colnames(sampleDistMatrix) <- NULL   
library(gplots)
colours = colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
heatmap.2(sampleDistMatrix, trace="none", col=colours)

# http://www.sthda.com/english/wiki/rna-seq-differential-expression-work-flow-using-deseq2

###---------------------------------------------------------------------------------------

# RLD
sampleDists <- dist(t(assay(rld)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$run)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Reds")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# NTD
sampleDists <- dist(t(assay(ntd)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(ntd$condition, ntd$run, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# VSD
sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$run, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Greens")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

## MDS plot
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = condition)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS (VST) - \n(Guillain-Barré Vs Controles)")



## MDS plot
mds <- as.data.frame(colData(rld))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = condition)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS (RLD) - \n(Guillain-Barré Vs Controles)")




# Using SVA with DESeq2
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)

# Surrogate variables 1 and 2 plotted over cell line.
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds$condition, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + condition
ddssva
res2 <- results(ddssva, contrast = c('condition', 'gbs', 'control'), filterFun=ihw, alpha = 0.05)

res2
summary(res2)



# VST - Variance Stabilizing Transformation 
vsd.2 <- vst(ddssva, blind=FALSE)


# Um subconjunto dos genes mais altamente veriáveis:
topVarGenes <- head(order(rowVars(assay(vsd.2)), decreasing = TRUE), 20)

mat  <- assay(vsd.2)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("condition", "replicate")])
pheatmap(mat, annotation_col = anno)

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup="condition")

plotCounts(ddssva, gene=which.min(res$padj), intgroup="condition")
topGene <- rownames(res2)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup="condition")


########################################### SELEÇÃO GENES UP e DOWN - SVA correction ###########################################
contr_gbs <- as.data.frame(res2)

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

# Filtrando por log2FC
DE_subset <- subset(gbs.DE, log2FoldChange >= 1.5 | log2FoldChange <= -1.5)
head(DE_subset, 9)
nrow(DE_subset)

# Filtrando por padj
DE_subset <- subset(DE_subset, padj <= 0.05)  # ou padj <= 0.05
nrow(DE_subset)

# Separando em genes up e down regulados
DE_down <- subset(DE_subset, log2FoldChange <= -1.5)  # Downregulated genes
length(DE_down$genes)
head(DE_down, 16)

# Tabela de contagem de genes 'downregulados'
#write.table(as.data.frame(DE_down), sep = ",", "./DE/gbs/transcript_count_gbs_rec_DOWN.txt", row.names = FALSE)
write.csv(as.data.frame(DE_down), file = './genes_DE/gbs/sva_correction/transcript_SVA_count_gbs_DOWN.csv')

DE_up <- subset(DE_subset, log2FoldChange >= 1.5)     # Upregulated genes
length(DE_up$genes)
head(DE_up, 16)

# Tabela de contagem de genes 'upregulados'
#write.table(as.data.frame(DE_up), sep = ",", "./DE/gbs/transcript_count_gbs_rec_UP.txt", row.names = FALSE)
write.csv(as.data.frame(DE_up), file = './genes_DE/gbs/sva_correction/transcript_SVA_count_gbs_UP.csv')

###################################################################################################################
