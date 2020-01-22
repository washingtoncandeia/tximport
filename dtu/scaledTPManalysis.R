##---------------------------------------------
# Análise em Nível de Transcrito Kallisto
# Utilizando Octuplicatas
# alpha = 0.05 (FDR, padj < 0.05) + IWH
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
library(iCOBRA)
library(reshape2)
library(Hmisc)
library(DEXSeq)
#BiocManager::install("DRIMSeq")
#BiocManager::install("DEXSeq")
#BiocManager::install("iCOBRA")
#install.packages("reshape2")
#install.packages("Hmisc")


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


# Coluna 3
# Nomes de amostras analisadas para fazer coluna run
run <- list.files(dir)
length(run)

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

# Coluna 8
# Replicatas
replicates <- c('rep01', 'rep02', 'rep03 ', 'rep04',
                'rep05', 'rep06 ', 'rep07', 'rep08')

## Parte II
# Aqui inicia-se a construção do data frame com todas as informações de amostras.
## Parte 2 - Unir cada vetor formando colunas de um data frame:
samples_info <- data.frame(pop = pop,
                           run = run,
                           condition = condition,
                           replicate = rep(replicates, 35))   

head(samples_info)
length(samples_info$condition)

## -------------------- Eliminando Linhas do Data Frame por Nomes  -------------------- ##
## Quando necessário o uso de um data frame menor, com apenas algumas variáveis.

# Usando dplyr, função filter e negando com regex (função grepl)

samples_info <- samples_info %>% 
  filter(!grepl(CHIKV_REC, pop))

samples_info <- samples_info %>% 
  filter(!grepl(CHIKV, pop))

samples_info <- samples_info %>% 
  filter(!grepl(GBS, pop))

samples_info <- samples_info %>% 
  filter(!grepl(GBS_REC, pop))

samples_info

## ------------------------------------------------------------------------------------ ##
length(samples_info$condition)
# Salvar a tabela no formato .txt (tsv)
write.table(samples_info, './tables/dtu/zika_gbs_vs_control.txt', sep = '\t')

# Criar um vetor nomeado apontando os arquivos de quantificação.
# Estes arquivos têm seus nomes anotados em uma tabela (samples.txt).
samples <- read.table('./tables/dtu/zika_gbs_vs_control.txt', header = TRUE, row.names = 1)
head(samples, 9)
samples$condition 

# Relevel: ajustando a condição referência para análise
samples$condition <- relevel(samples$condition, ref = 'control')

# Nomeando as linhas com nome de cada arquivo de amostra:
rownames(samples) <- samples$run
head(samples, 25)

# Obtendo cada arquivo de replicata das amostras usadas em kallisto:
files <- file.path(dir, samples$run, 'abundance.h5')
files
names(files) <- samples$run
head(files, 9)

# Criação de IDs Ensembl
# Arquivo t2g data frame
t2g <- read_csv('./t2sg.csv')
head(t2g, 12)
nrow(t2g)


# tximport:
# 1. simple_sum 
txi.kallisto <- tximport(files, 
                         type = 'kallisto',
                         tx2gene = t2g)
# 2. scaledTPM
txi.kallisto.scaledTPM <- tximport(files, 
                         type = 'kallisto',
                         tx2gene = t2g, 
                         txOut = FALSE,
                         countsFromAbundance = "scaledTPM")

# 3. Offset Matrix para uso com simple_sum
txi.kallisto.tx <- tximport(files,
                            type = "kallisto",
                            tx2gene = t2g, 
                            txOut = TRUE,
                            countsFromAbundance = "no")

# 4. DTU
txi.kallisto.dtu <- tximport(files, 
                             type = 'kallisto',
                             tx2gene = t2g,
                             txOut = TRUE,
                             countsFromAbundance = "dtuScaledTPM")


# Análise baseada em Sonensens et al (2016)
# https://f1000researchdata.s3.amazonaws.com/datasets/7563/e22fa771-ee2b-486a-8cf4-31b645163d6b_gse64570_quantification.html#deseq2

kallisto_quant <- list(geneCOUNT_kall_simplesum = txi.kallisto$counts,
                       geneCOUNT_kall_scaledTPM = txi.kallisto.scaledTPM$counts,
                       avetxlength = txi.kallisto$length,
                       geneTPM_kall = txi.kallisto$abundance,
                       txTPM_kall = txi.kallisto.tx$abundance,
                       txCOUNT_kall = txi.kallisto.tx$counts,
                       txi_kallistosimplesum = txi.kallisto,
                       txi_kallistoscaledTPM = txi.kallisto.scaledTPM,
                       txi_kallistotx = txi.kallisto.tx)



res_kall_simplesum_avetxl_deseq2 <-  diff_expression_DESeq2(txi = kallisto_quant$txi_kallistosimplesum,
                                                            counts = NULL,
                                                            meta = samples, cond_name = "condition",
                                                            sample_name = "run",
                                                            level1 = "control",
                                                            level2 = "gbs")

res_kall_scaledTPM_deseq2 <- diff_expression_DESeq2(txi = NULL,
                                                    counts = kallisto_quant$geneCOUNT_kall_scaledTPM,
                                                    meta = samples, cond_name = "condition",
                                                    sample_name = "run",
                                                    level1 = "control",
                                                    level2 = "gbs")
  
res_kall_simplesum_deseq2 <- diff_expression_DESeq2(txi = NULL,
                                                    counts = kallisto_quant$geneCOUNT_kall_simplesum, 
                                                    meta = samples, cond_name = "condition",
                                                    sample_name = "run",
                                                    level1 = "control",
                                                    level2 = "gbs")



dfh <- data.frame(pvalue = c(res_kall_scaledTPM_deseq2$res$pvalue,
                             res_kall_simplesum_avetxl_deseq2$res$pvalue,
                             res_kall_simplesum_deseq2$res$pvalue),
                  mth = c(rep("A. scaledTPM_kallisto, Guillain-Barré", nrow(res_kall_scaledTPM_deseq2$res)),
                          rep("B. simplesum_kallisto_avetxl, Guillain-Barré", nrow(res_kall_simplesum_avetxl_deseq2$res)),
                          rep("C. simplesum_kallisto, Guillain-Barré)", nrow(res_kall_simplesum_deseq2$res))))
ggplot(dfh, aes(x = pvalue)) + geom_histogram() + facet_wrap(~mth) + 
  plot_theme() + 
  xlab("p-value") + ylab("count") 

# PlotDispEsts DESeq2
par(mfrow = c(1, 3))
plotDispEsts(res_kall_scaledTPM_deseq2$dsd, main = "A. scaledTPM_kallisto, Guillain-Barré")
plotDispEsts(res_kall_simplesum_avetxl_deseq2$dsd, main = "B. simplesum_kallisto_avetxl, Guillain-Barré")
plotDispEsts(res_kall_simplesum_deseq2$dsd, main = "C. simplesum_kallisto, Guillain-Barré")


# Ma plots
# avetxl = average transcripts lenght
par(mfrow = c(1, 3))
DESeq2::plotMA(res_kall_scaledTPM_deseq2$dsd, main = "A. scaledTPM_kallisto, Guillain-Barré")
DESeq2::plotMA(res_kall_simplesum_avetxl_deseq2$dsd, main = "B. simplesum_kallisto_avetxl, Guillain-Barré") 
DESeq2::plotMA(res_kall_simplesum_deseq2$dsd, main = "C. simplesum_kallisto, Guillain-Barré")

# Comparação de genes significantes encontrados em diferentes matrizes
cobra_deseq2 <- COBRAData(padj = data.frame(simplesum_kallisto = res_kall_simplesum_deseq2$res$padj, 
                                            row.names = rownames(res_kall_simplesum_deseq2$res)))

cobra_deseq2 <- COBRAData(padj = data.frame(scaledTPM_kallisto = res_kall_scaledTPM_deseq2$res$padj, 
                                            row.names = rownames(res_kall_scaledTPM_deseq2$res)),
                          object_to_extend = cobra_deseq2)

cobra_deseq2 <- COBRAData(padj = data.frame(simplesum_kallisto_avetxl = res_kall_simplesum_avetxl_deseq2$res$padj, 
                                            row.names = rownames(res_kall_simplesum_avetxl_deseq2$res)),
                          object_to_extend = cobra_deseq2)

cobraperf_deseq2 <- calculate_performance(cobra_deseq2, aspects = "overlap", thr_venn = 0.05)


cobraplot1_deseq2 <- prepare_data_for_plot(cobraperf_deseq2, incltruth = FALSE, 
                                           colorscheme = c("blue", "red", "black"))
plot_overlap(cobraplot1_deseq2, cex = c(1, 0.7, 0.7))
title("Febre Guillain-Barré Rec Vs Controles", line = 0)


# Comparação de estimativas logFC todos os metodos
df1 <- Reduce(function(...) merge(..., by = "gene", all = TRUE), 
              list(data.frame(gene = rownames(res_kall_scaledTPM_deseq2$res),
                              scaledTPM_kallisto = res_kall_scaledTPM_deseq2$res$log2FoldChange,
                              stringsAsFactors = FALSE),
                   data.frame(gene = rownames(res_kall_simplesum_deseq2$res),
                              simplesum_kallisto = res_kall_simplesum_deseq2$res$log2FoldChange, 
                              stringsAsFactors = FALSE),
                   data.frame(gene = rownames(res_kall_simplesum_avetxl_deseq2$res),
                              simplesum_kallisto_avetxl =
                                res_kall_simplesum_avetxl_deseq2$res$log2FoldChange,
                              stringsAsFactors = FALSE)))
rownames(df1) <- df1$gene
df1$gene <- NULL
pairs(df1, upper.panel = panel_smooth, lower.panel = panel_cor)

# Comparison of logFC estimates - simplesum vs scaledTPM
df2 <- Reduce(function(...) merge(..., by = "gene", all = TRUE), 
              list(data.frame(gene = rownames(kallisto_quant$geneTPM_kall),
                              kallisto_quant$geneTPM_kall,
                              stringsAsFactors = FALSE),
                   data.frame(gene = rownames(res_kall_scaledTPM_deseq2$res),
                              scaledTPM_kallisto_logFC = res_kall_scaledTPM_deseq2$res$log2FoldChange,
                              scaledTPM_kallisto_basemean = res_kall_scaledTPM_deseq2$res$baseMean, 
                              stringsAsFactors = FALSE),
                   data.frame(gene = rownames(res_kall_simplesum_deseq2$res),
                              simplesum_kallisto_logFC = res_kall_simplesum_deseq2$res$log2FoldChange, 
                              simplesum_kallisto_basemean = res_kall_simplesum_deseq2$res$baseMean, 
                              stringsAsFactors = FALSE),
                   data.frame(gene = rownames(res_kall_simplesum_avetxl_deseq2$res),
                              simplesum_kallisto_avetxl_logFC =
                                res_kall_simplesum_avetxl_deseq2$res$log2FoldChange,
                              simplesum_kallisto_avetxl_basemean = 
                                res_kall_simplesum_avetxl_deseq2$res$baseMean, 
                              stringsAsFactors = FALSE)))

meta <- samples
rownames(df2) <- df2$gene
df2$gene <- NULL
df2$scaledTPM_kallisto_basemeanbinary <- Hmisc::cut2(df2$scaledTPM_kallisto_basemean, g = 2)
df2$simplesum_kallisto_basemeanbinary <- Hmisc::cut2(df2$simplesum_kallisto_basemean, g = 2)
df2$simplesum_kallisto_avetxl_basemeanbinary <- Hmisc::cut2(df2$simplesum_kallisto_avetxl_basemean, g = 2)
df2$sumA <- rowSums(df2[, meta$run[meta$condition == "control"]])
df2$sumB <- rowSums(df2[, meta$run[meta$condition == "gbs_rec"]])
df2$allzero_onecond <- "Expresso em ambos os grupos"
df2$allzero_onecond[union(which(df2$sumA == 0), which(df2$sumB == 0))] <- "Expresso em um grupo"
df2$onecol <- rep("", nrow(df2))

ggplot(df2, aes(x = simplesum_kallisto_logFC, y = scaledTPM_kallisto_logFC, col = onecol)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_point(size = 2, alpha = 0.5) + 
  plot_theme() + ggtitle("Estimativas logFC, Guillain-Barré rec") + theme(legend.position = "bottom") + 
  scale_color_manual(values = c("blue"), name = "") + 
  theme(legend.background = element_rect(fill = "white"), legend.key = element_blank()) + 
  xlab("simplesum_kallisto, logFC") + ylab("scaledTPM_kallisto, logFC") + 
  guides(colour = guide_legend(override.aes = list(size = 0)))


# Em ambos e em um só grupo
ggplot(df2, aes(x = simplesum_kallisto_logFC, y = scaledTPM_kallisto_logFC, col = allzero_onecond)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_point(size = 2, alpha = 0.5) + 
  plot_theme() + ggtitle("Estimativas logFC, Guillain-Barré rec") + theme(legend.position = "bottom") + 
  scale_color_manual(values = c("blue", "red"), name = "") + 
  xlab("simplesum_kallisto, logFC") + ylab("scaledTPM_kallisto, logFC") + 
  guides(colour = guide_legend(override.aes = list(size = 7)))


# base mean
ggplot(subset(df2, !is.na(scaledTPM_kallisto_basemeanbinary)), 
       aes(x = simplesum_kallisto_logFC, y = scaledTPM_kallisto_logFC, 
           col = scaledTPM_kallisto_basemeanbinary)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_point(size = 2, alpha = 0.5) + 
  facet_wrap(~scaledTPM_kallisto_basemeanbinary) + 
  plot_theme() + ggtitle("Estimativas logFC, Guillain-Barré rec") + theme(legend.position = "bottom") + 
  xlab("simplesum_kallisto, logFC") + ylab("scaledTPM_kallisto, logFC") + 
  scale_color_manual(values = c("red", "blue"), name = "scaledTPM_kallisto, base mean") + 
  guides(colour = guide_legend(override.aes = list(size = 7)))


# Criação de IDs Ensembl
# Arquivo t2g data frame
t2g <- read_csv('./t2sg.csv')
head(t2g, 12)
nrow(t2g)


# DTU analysis com DEXSeq
BPPARAM = MulticoreParam(6)
stopifnot(all(colnames(kallisto_quant$txCOUNT_kall) == rownames(meta)))
dxd <- DEXSeqDataSet(countData = round(kallisto_quant$txCOUNT_kall), sampleData = meta, 
                     design = ~sample + exon + condition:exon,
                     featureID = rownames(kallisto_quant$txCOUNT_kall),
                     groupID = t2g$SYMBOL[match(rownames(kallisto_quant$txCOUNT_kall), 
                                                  t2g$target_id)])
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)
