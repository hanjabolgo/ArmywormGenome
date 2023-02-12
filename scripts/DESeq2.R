library("DESeq2")
# first day pheromone glands ~ third day pheromone glands
countdataraw <- read.csv("TPM.tab", header=T, sep=',',row.names=1)
countdataraw <- as.matrix(countdataraw)
head(countdataraw)
countdata <- countdataraw[, c(1:6)]
countData <- countdata
condition <- factor(c("A","A","A","B","B","B"))
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
summary(res)
table(res$padj<0.05)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_deseq2_1PGvs3PG_005LFC1 <- row.names(diff_gene_deseq2)
length(diff_gene_deseq2_1PGvs3PG_005LFC1)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file= "1PG_vs_3PG.cvs",row.names = F)

# third day pheromone glands ~ mated pheromone glands
countdata <- countdataraw[, c(4:9)]
head(countdata)
countData <- countdata
condition <- factor(c("A","A","A","B","B","B"))
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
head(dds)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
summary(res)
table(res$padj<0.05)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_deseq2_3PGvs4PG_005LFC1 <- row.names(diff_gene_deseq2)
length(diff_gene_deseq2_3PGvs4PG_005LFC1)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file= "3PG_vs_4PG.cvs",row.names = F)

# Male antenna ~ female antenna
countdata <- countdataraw[, c(10:15)]
head(countdata)
countData <- countdata
condition <- factor(c("A","A","A","B","B","B"))
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
head(dds)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
summary(res)
table(res$padj<0.05)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_deseq2_MAvsFA_005LFC1 <- row.names(diff_gene_deseq2)
length(diff_gene_deseq2_MAvsFA_005LFC1)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file= "MA_vs_FA.cvs",row.names = F)
