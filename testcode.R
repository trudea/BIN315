library(DESeq2)

# Differential expression analysis of Colon and Rectum data

SampleTable <- data.frame(classes2)
expr <- round(expr2)
dds <- DESeqDataSetFromMatrix(expr, colData = SampleTable, design = ~classes2)


dds <- DESeq(dds)
# Sets sample name according to where it was retrieved
colnames(expr) <- classes2

res <- results(dds)

# With a FDR cutoff of 0.05
res0.05 <- subset(res, padj < 0.05)
res0.05[order(res0.05$padj),]


res_new <- na.omit(res)
res_1 <-  subset(res_new, res_new$log2FoldChange >= 1 & padj < 0.05)
res_2 <- subset(res_new, res_new$log2FoldChange <= -1 & padj < 0.05)

# Skjøter sammen tabellene
res_n <- rbind(res_1, res_2)


DESeq2::plotMA(res, alpha = 0.05)

id_table <- read.table("Human_ensembl_ids_to_symbols.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Example vector of gene names
gene.names <- c("ENSG00000143842", "ENSG00000144810", "ENSG00000146966")

# Merge gene names 
gene.symbols <- merge(data.frame(ensembl_id = gene.names), id_table)$gene_symbol

write(gene.symbols, file="genes.txt", sep = "\n")

up_regulated <- attributes(res_1)$rowname
up_regulated_symbol <- merge(data.frame(ensembl_id = up_regulated), id_table)$gene_symbol
write(up_regulated_symbol, file="upregulatedgenes.txt", sep="\n")

down_regulated <- attributes(res_2)$rowname
down_regulated_symbol <- merge(data.frame(ensembl_id = down_regulated), id_table)$gene_symbol
write(down_regulated, file="downregulatedgenes.txt", sep="\n")






modules <- blockwiseModules(datExpr = t(expr2), power = 10, minModuleSize = 40, minAbsSplitHeight = 0.1)

# Number of genes per module
mod_genes = c()
for (i in 1:ncol(modules$MEs)) {
  mod_genes[i] <- sum(paste("ME", modules$colors, sep = "") == colnames(modules$MEs)[i])
}
names(mod_genes) <- colnames(modules$MEs)
xx <- barplot(mod_genes, col = gsub("ME", "", colnames(modules$MEs)),
              main = "Number of genes per module",
              las = 2, cex.names = 0.65,)
text(x = xx, y = mod_genes, label = mod_genes, pos = 3, cex = 0.8)


# Compute correlation between all genes and all module eigengenes
kME <- cor(exp, modules$MEs)

# Correlation within modules
intra_cor <- c()
for (m in colnames(kME)) {
  if (m != "MEgrey") {
    intra_cor <- c(intra_cor, kME[paste("ME", modules$colors, sep = "") == m, m])
  }
}
hist(intra_cor, xlim = c(-1,1), breaks = seq(-1,1,0.1),
     main = "Correlations with module eigengene (within module correlation)",
     xlab = "Correlation")


# Correlation between modules
idx <- which(colnames(modules$MEs) == "MEgrey")
MEs_R <- cor(modules$MEs[,-idx], modules$MEs[,-idx])
hist(MEs_R[upper.tri(MEs_R)], xlim = c(-1,1), breaks = seq(-1,1,0.1),
     main = "Correlations of module eigengenes (between module correlation)",
     xlab = "Correlation")


# Find samples labled as colon
colon <- dds$classes2 == "Colon"

# Run for each ME in modules: a t.test between autism and control samples, saving the pvalues from the test
MEpval <- unlist(lapply(modules$MEs, function(ME){ t.test(ME[colon], ME[!colon])$p.value })) 

# Calculate pvalues adjusted for multiple testing
MEpval.adj <- p.adjust(MEpval, method="fdr")

# Turn adjusted pvalue into a positive, scaled "significance" value by taking the inverse log10
MEsignificance <- -log10(MEpval.adj)

# Plot ME significance
barplot(MEsignificance, col = substr(names(MEpval), 3, 20), las = 2, cex.names = 0.65,
        ylab = "module eigengene significance")
abline(h=-log10(0.05), lty = 3)

# ME midnight blue and dark orange seem to be significant.... 

df <- data.frame(modules$MEs[c("MEmidnightblue", "MEdarkorange2")], colon)


colors = c("red", "grey")
# then we plot with inedxing, appears to be more upregulated in midnight blue in the colon and dark orange....

barplot(df$MEmidnightblue, main= "Model Eigengene Midnight Blue" ,col = colors[as.factor(df$colon)])
legend("bottomright", legend = c("Colon", "Rectum"), fill = c("Red", "Grey"))


barplot(df$MEdarkorange2, main= "Model Eigengene Dark Orange 2" ,col = colors[as.factor(df$colon)])
legend("bottomright", legend = c("Colon", "Rectum"), fill = c("Red", "Grey"))





# Measure module membership by calculating the correlation between genes and module eigengenes
kME <- cor(t(expr2), modules$MEs)
kME <- as.data.frame(kME) # Convert to a data.frame

symbols <- data.frame(rownames(expr2))
Purple.genes <- data.frame(kME$MEmidnightblue, symbols)
# Sort by correlation
Purple.genes <- Purple.genes[order(-Purple.genes$kME.MEmidnightblue), ]
Central.genes <- Purple.genes[1:10, ]
