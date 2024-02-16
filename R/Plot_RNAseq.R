
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Dm.eg.db)





# Load your count matrix and sample information
# Assuming 'countData' is your matrix of read counts (genes x samples)
# and 'colData' is a DataFrame describing the condition of each sample

RNA_to_prot <- function(






  )


cts <- as.matrix(read.csv("data/pasilla_gene_counts.tsv", row.names="gene_id", sep ="\t"))

coldata <- read.csv("data/pasilla_sample_annotation.csv", sep = ",", row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)
rownames(coldata) <- sub("fb", "", rownames(coldata))
cts <- cts[, rownames(coldata)]

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

# Run differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# Export gene expression results
res_df <- as.data.frame(res)
#write.csv(as.data.frame(res), "Output/gene_expression.csv")


# Functional annotation and pathway analysis
# Identify significantly differentially expressed genes
sig_genes <- rownames(subset(res, padj < 0.05))

# Convert gene symbols to Entrez IDs (example for human genes)
entrez_ids <- mapIds(org.Dm.eg.db, keys = sig_genes, column = "ENTREZID",
                     keytype = "FLYBASE",
                     multiVals = "first")

# GO enrichment analysis
go_enrichment <- enrichGO(gene = entrez_ids, OrgDb = org.Dm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)

# KEGG pathway analysis
kegg_pathway <- enrichKEGG(gene = entrez_ids, organism = 'dme', keyType = 'entrez', pAdjustMethod = "BH", qvalueCutoff = 0.05)

# Export pathway analysis results
write.csv(as.data.frame(go_enrichment), "Output/GO_enrichment_results.csv")
write.csv(as.data.frame(kegg_pathway), "Output/KEGG_pathway_results.csv")














