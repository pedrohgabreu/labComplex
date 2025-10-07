library(DESeq2)
library(edgeR)
library(openxlsx)
library(SummarizedExperiment) 
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("affy")
library("ggfortify")
library(clusterProfiler)
library(fgsea)

#criou a matriz e colocou o geneid como nome da linha - DADOS BAIXADOS DO GEO(Desktop/Pedro/GSE272957_CountsTable_v2.txt')
count_matrix <- read.table('Desktop/Pedro/GSE272957_CountsTable_v2.txt', row.names = 1)
rownames(count_matrix) <- count_matrix$Geneid
#excluiu a coluna geneid
count_matrix$Geneid <- NULL

#pegar s metadados direto do site geo omnibus(GSE272957 CODIGO DE ACESSO AOS DADDOS)
meta <- getGEO("GSE272957", GSEMatrix = TRUE)
metadata_GSE272957 <- pData(meta[[1]])


# UTLIZAR FILTRO AUTOMATICO
y <- DGEList(counts=GSE224004_filt, group=metadata_GSE224004_filt$condition)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
GSE224004_filt <- y$counts

# FILTRO MANUAL
# Filter rows with row sum higher than 10
keep <- rowSums(count_matrix) > 10
count_matrix_filt <- count_matrix[keep,]

# Filter rows with values higher than 10 and rowsums higher than 3
keep <- rowSums(count_matrix_filt >= 10) >= 3
count_matrix_filt <- count_matrix_filt[keep,]


boxplot(log2(count_matrix_filt+1), main = "Distribuição de Contagens", 
        #  xlab = "Amostras", 
        ylab = "Counts", 
        #  col = c("skyblue", "pink", "lightgreen"),
        las = 2) 




hist(as.numeric(unlist(log2(count_matrix_filt+1))), breaks = 100,
     col = "tomato",
     main = "Distribuição log10 das contagens (RNA-Seq)",
     xlab = "log10(Contagem + 1)",
     ylab = "Frequência")


# Transpor matriz de counts transformados
counts_genes_filt_t <- t(count_matrix_filt)

# Setar nomes dos grupos
rownames(counts_genes_filt_t) <- metadata_GSE272957$`genotype:ch1`

# Clusterização hierárquica
d <- dist(as.matrix(counts_genes_filt_t), method = "euclidean")
clusters <- hclust(d, method = "complete")

# Plotar
par(mar = c(1, 1, 1, 1))
plot(clusters)

# Normalization
library('preprocessCore')
set.seed(262)


#normalizacao da matriz
count_matrix_filt_norm <- normalize.quantiles(as.matrix(count_matrix_filt))
#colocando o nome das colunas e linhas na matriz normalizada
colnames(count_matrix_filt_norm) <- colnames(count_matrix_filt)
rownames(count_matrix_filt_norm) <- rownames(count_matrix_filt)
count_matrix_filt_norm <- as.data.frame(count_matrix_filt_norm)


boxplot(log2(count_matrix_filt_norm+1), main = "Distribuição de Contagens", 
        #  xlab = "Amostras", 
        ylab = "Counts", 
        #  col = c("skyblue", "pink", "lightgreen"),
        las = 2) 

hist(as.numeric(unlist(log2(count_matrix_filt_norm+1))), breaks = 100,
     col = "tomato",
     main = "Distribuição log10 das contagens (RNA-Seq)",
     xlab = "log10(Contagem + 1)",
     ylab = "Frequência")


# copiado de RNAseq_keyfiles

# "second normalization"
d0 <- DGEList(count_matrix_filt_norm)

d0 <- calcNormFactors(d0)
d0


cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

# Limma

group <- as.factor(metadata_GSE272957$`genotype:ch1`)
group

# MDS
plotMDS(d, col = as.numeric(group), labels=metadata_GSE272957$`genotype:ch1`)
abline(h=0,col="grey")
abline(v=0,col="grey")


mm <- model.matrix(~0 + group)
colnames(mm) <- gsub('-', '_', levels(group))
colnames(mm) <-  gsub(' ', '_', colnames(mm))
mm

y <- voom(d, mm, plot = T)

fit <- lmFit(y, mm)
head(coef(fit))


contr <- makeContrasts('EWS_FLI1-Control_lentivirus', 'EWS_FLI1-CTRL', levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr)

tmp <- eBayes(tmp)

summary(decideTests(tmp, method="global"))

EWS_FLI1_Control_lentivirus_sf <- topTable(tmp, sort.by ="P", n = Inf, coef = 1)
EWS_FLI1_CTRL_sf <- topTable(tmp, sort.by ="P", n = Inf, coef = 1)

EWS_FLI1_Control_lentivirus <- topTable(tmp, sort.by ="P", n = Inf, coef = 1, p.value = 0.05, lfc = 0.58)
EWS_FLI1_CTRL <- topTable(tmp, sort.by ="P", n = Inf, coef = 2, p.value = 0.05, lfc = 0.58)

