# Função: dff_expression_DESeq2 

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



# Função: plot_theme 
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


# Função: panel_cor 
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


# Função: panel_smooth
panel_smooth <-function (x, y, col = "blue", bg = NA, pch = ".", 
                         cex = 0.8, ...) {
  ## Panel function to plot points
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
}




# Função: dff_expression_DESeq2

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



# Função: plot_theme 

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


# Função: panel_cor
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


panel_smooth<-function (x, y, col = "blue", bg = NA, pch = ".", 
                        cex = 0.8, ...) {
  ## Panel function to plot points
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
}

