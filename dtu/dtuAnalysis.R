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
t2g <- read_csv('./t2gs.csv')
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
                       geneCOUNT_kall_dtu = txi.kallisto.dtu$counts,
                       geneCOUNT_kall_scaledTPM = txi.kallisto.scaledTPM$counts,
                       avetxlength = txi.kallisto$length,
                       geneTPM_kall = txi.kallisto$abundance,
                       txTPM_kall = txi.kallisto.tx$abundance,
                       txCOUNT_kall = txi.kallisto.tx$counts,
                       txi_kallistosimplesum = txi.kallisto,
                       txi_kallistosdtu = txi.kallisto.dtu,
                       txi_kallistoscaledTPM = txi.kallisto.scaledTPM,
                       txi_kallistotx = txi.kallisto.tx)


#-----------Função: dff_expression_DESeq2 

diff_expression_DESeq2 <- function(txi = NULL, counts, meta, cond_name, 
                                   level1, level2, sample_name) {
  ## Differential expression analysis with DESeq2
  
  suppressPackageStartupMessages(library(DESeq2))
  
  ## If tximport object provided, generate DESeqDataSet from it. Otherwise, 
  ## use the provided count matrix.
  if (!is.null(txi)) {
    txi$counts <- round(txi$counts)
    keep_feat <- rownames(txi$counts[rowSums(is.na(txi$counts)) == 0 & rowSums(txi$counts) != 0, ])
    txi <- lapply(txi, function(w) {
      if (!is.null(dim(w))) w[match(keep_feat, rownames(w)), ]
      else w
    })
    dsd <- DESeqDataSetFromTximport(txi, 
                                    colData = meta[match(colnames(txi$counts), 
                                                         meta[, sample_name]), ],
                                    design = as.formula(paste0("~", cond_name)))
  } else {
    counts <- round(counts)
    cts = counts[rowSums(is.na(counts)) == 0, ]
    cts <- cts[rowSums(cts) != 0, ]
    dsd <- DESeqDataSetFromMatrix(countData = round(cts), 
                                  colData = meta[match(colnames(cts), 
                                                       meta[, sample_name]), ],
                                  design = as.formula(paste0("~", cond_name)))
  }
  
  ## Estimate dispersions and fit model
  dsd <- DESeq(dsd, test = "Wald", fitType = "local", betaPrior = TRUE)
  res <- as.data.frame(results(dsd, contrast = c(cond_name, level2, level1),
                               cooksCutoff = FALSE, independentFiltering = FALSE))
  return(list(dsd = dsd, res = res))
}
#---


###-- Função: plot_theme 

plot_theme <- function() {
  ## ggplot2 plotting theme
  theme_grey() +
    theme(legend.position = "right",
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text = element_text(size = 10),
          strip.background = element_rect(fill = NA, colour = "black"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(colour = "black", size = 20))
}

##--

###--- Função: panel_cor 
panel_cor <- function(x, y, digits = 3, cex.cor) {
  ## Panel function to print Pearson and Spearman correlations
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1 <- abs(cor(x, y, method = "pearson", use = "complete"))
  txt1 <- format(c(r1, 0.123456789), digits = digits)[1]
  r2 <- abs(cor(x, y, method = "spearman", use = "complete"))
  txt2 <- format(c(r2, 0.123456789), digits = digits)[1]
  text(0.5, 0.35, paste("pearson =", txt1), cex = 1.1)
  text(0.5, 0.65, paste("spearman =", txt2), cex = 1.1)
}


panel_smooth <-function (x, y, col = "blue", bg = NA, pch = ".", 
                        cex = 0.8, ...) {
  ## Panel function to plot points
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
}
  
  
  
# Estimativas de transcritos: DGE e DTU analysis  
  
res_kall_simplesum_avetxl_deseq2 <-  diff_expression_DESeq2(txi = kallisto_quant$txi_kallistosimplesum,
                                                              counts = NULL,
                                                              meta = samples, cond_name = "condition",
                                                              sample_name = "run",
                                                              level1 = "control",
                                                              level2 = "zika")

res_kall_dtu_deseq2 <- diff_expression_DESeq2(txi = NULL,
                                                    counts = kallisto_quant$geneCOUNT_kall_dtu,
                                                    meta = samples, cond_name = "condition",
                                                    sample_name = "run",
                                                    level1 = "control",
                                                    level2 = "zika")

res_kall_simplesum_deseq2 <- diff_expression_DESeq2(txi = NULL,
                                                    counts = kallisto_quant$geneCOUNT_kall_simplesum, 
                                                    meta = samples, cond_name = "condition",
                                                    sample_name = "run",
                                                    level1 = "control",
                                                    level2 = "zika")


res_kall_scaledTPM_deseq2 <- diff_expression_DESeq2(txi = NULL,
                                                    counts = kallisto_quant$geneCOUNT_kall_scaledTPM,
                                                    meta = samples, cond_name = "condition",
                                                    sample_name = "run",
                                                    level1 = "control",
                                                    level2 = "zika")

dfh <- data.frame(pvalue = c(res_kall_dtu_deseq2$res$pvalue,
                             res_kall_scaledTPM_deseq2$res$pvalue,
                             res_kall_simplesum_avetxl_deseq2$res$pvalue,
                             res_kall_simplesum_deseq2$res$pvalue),
                  mth = c(rep("dtu_kallisto, Zika", nrow(res_kall_dtu_deseq2$res)),
                          rep("scaledTPM_kallisto, Zika", nrow(res_kall_scaledTPM_deseq2$res)),
                          rep("simplesum_kallisto_avetxl, Zika", nrow(res_kall_simplesum_avetxl_deseq2$res)),
                          rep("simplesum_kallisto, Zika", nrow(res_kall_simplesum_deseq2$res))))
ggplot(dfh, aes(x = pvalue)) + geom_histogram() + facet_wrap(~mth) + 
  #plot_theme() + 
  xlab("p-value") + ylab("count") 

# PlotDispEsts DESeq2
par(mfrow = c(2, 2))
plotDispEsts(res_kall_dtu_deseq2$dsd, main = "dtu_kallisto, Zika")
plotDispEsts(res_kall_scaledTPM_deseq2$dsd, main = "scaledTPM_kallisto, Zika")
plotDispEsts(res_kall_simplesum_avetxl_deseq2$dsd, main = "simplesum_kallisto_avetxl, Zika")
plotDispEsts(res_kall_simplesum_deseq2$dsd, main = "simplesum_kallisto, Zika")


# Ma plots
# avetxl = average transcripts lenght
par(mfrow = c(2, 2))
DESeq2::plotMA(res_kall_dtu_deseq2$dsd, main = "dtu_kallisto, Zika")
DESeq2::plotMA(res_kall_scaledTPM_deseq2$dsd, main = "scaledTPM_kallisto, Zika")
DESeq2::plotMA(res_kall_simplesum_avetxl_deseq2$dsd, main = "simplesum_kallisto_avetxl, Zika") 
DESeq2::plotMA(res_kall_simplesum_deseq2$dsd, main = "simplesum_kallisto, Zika")


par(mfrow = c(1, 1))


# Comparação de genes significantes encontrados em diferentes matrizes
cobra_deseq2 <- COBRAData(padj = data.frame(simplesum_kallisto = res_kall_simplesum_deseq2$res$padj, 
                                            row.names = rownames(res_kall_simplesum_deseq2$res)))

cobra_deseq2 <- COBRAData(padj = data.frame(dtu_kallisto = res_kall_dtu_deseq2$res$padj, 
                                            row.names = rownames(res_kall_dtu_deseq2$res)),
                          object_to_extend = cobra_deseq2)

cobra_deseq2 <- COBRAData(padj = data.frame(scaledTPM_kallisto = res_kall_scaledTPM_deseq2$res$padj, 
                                            row.names = rownames(res_kall_scaledTPM_deseq2$res)),
                          object_to_extend = cobra_deseq2)

cobra_deseq2 <- COBRAData(padj = data.frame(simplesum_kallisto_avetxl = res_kall_simplesum_avetxl_deseq2$res$padj, 
                                            row.names = rownames(res_kall_simplesum_avetxl_deseq2$res)),
                          object_to_extend = cobra_deseq2)

cobraperf_deseq2 <- calculate_performance(cobra_deseq2, aspects = "overlap", thr_venn = 0.05)


cobraplot1_deseq2 <- prepare_data_for_plot(cobraperf_deseq2, incltruth = FALSE, 
                                           colorscheme = c("blue", "black", "green", "red"))
plot_overlap(cobraplot1_deseq2, cex = c(1, 0.7, 0.7))
title("Febre Zika Vs Controles, DTU", line = 0)



df1 <- Reduce(function(...) merge(..., by = "gene", all = TRUE), 
              list(data.frame(gene = rownames(res_kall_scaledTPM_deseq2$res),
                              scaledTPM_kallisto = res_kall_scaledTPM_deseq2$res$log2FoldChange,
                              stringsAsFactors = FALSE),
                   data.frame(gene = rownames(res_kall_dtu_deseq2$res),
                              dtu_kallisto = res_kall_dtu_deseq2$res$log2FoldChange,
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


res_kall_dtu_deseq2$res
res_kall_dtu_deseq2$res$log2FoldChange
head(df1)
