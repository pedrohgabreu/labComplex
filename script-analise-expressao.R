# instalar packages do CRAN e Bioconductor
cran_pkgs <- c("tidyverse", "dplyr")
bio_pkgs <- c("oligo", "GEOquery")

install.packages(cran_pkgs)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(bio_pkgs)

# carregar libs
library(tidyverse)
library(dplyr)
library(oligo)
library(GEOquery)

# primeiramente, baixar e extrair os dados do site
# codigo para extrair os metadados e dados do chip (IDs para simbolos de genes) 
meta <- getGEO("GSE162785", GSEMatrix = TRUE)
metadata_GSE162785 <- pData(meta[[1]])
platforms <- lapply(meta, annotation)
platform_id <- platforms[[1]]  
platform_info <- getGEO(platform_id, destdir = ".")
platform_GSE162785 <- Table(platform_info)

# Formatando metadados
metadata_GSE162785$condition <- gsub(".*\\+\\s*", "", metadata_GSE162785$title)

# fitros para deixar controle 1 e controle 2 como controle
metadata_GSE162785[metadata_GSE162785$condition == 'control 1' | metadata_GSE162785$condition == 'control 2','condition'] <- 'control'

 # colocar caminho dos dados baixados e extraidos
celFiles <- list.files('Desktop/Pedro/Amostras/', full = TRUE, pattern = '\\.CEL.gz', ignore.case = TRUE)

#armazenar dados em rawdata
rawData <- read.celfiles(celFiles)

#Normalizar com RMA
normalized.data_GSE162785 <- oligo::rma(rawData)

#salvar dados matrix de expressao em normalized.expr
GSE162785 <- as.data.frame(exprs(normalized.data_GSE162785))

# agrupar os metadados dentro da matriz de expressao
metadata_GSE162785$title
pData(normalized.data_GSE162785)

# Comparar odernamento de nomes das amostras entre pData e matrix de expressao
all(colnames(exprs(normalized.data_GSE162785)) == rownames(pData(normalized.data_GSE162785)))

# Se TRUE
pData(normalized.data_GSE162785)$condition <- metadata_GSE162785$condition

head(pData(normalized.data_GSE162785))

#formatar nome das linhas do pdata para ficar igual o nome das linhas do metadata
#rownames(pData(normalized.data_GSE162785)) <- gsub("_.*$", "", rownames(pData(normalized.data_GSE162785)))
#head(pData(normalized.data_GSE162785))

#colnames(exprs(normalized.data_GSE162785)) <- rownames(pData(normalized.data_GSE162785))

# Controle de qualidade
library(arrayQualityMetrics)

arrayQualityMetrics(expressionset = normalized.data_GSE162785,
                    outdir = 'Desktop/Pedro/relatorio/',
                    force = TRUE,
                    intgroup = 'condition')

#daqui para baixo códigos para plotPCA
# PCA
library(PCAtools)

p <- pca(as.matrix(GSE162785), metadata = pData(normalized.data_GSE162785), removeVar = 0.1)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p, lab = NULL, colby = "condition",  legendPosition = 'right')


## 3d PCA
library(plotly)
library(stats)
library(RColorBrewer)

prin_comp <- prcomp(t(as.matrix(GSE162785)), rank. = 3)


components <- prin_comp[["x"]]

components <- data.frame(components)

components$PC2 <- -components$PC2

components$PC3 <- -components$PC3

components = cbind(components, pData(normalized.data_GSE162785)$condition)


tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]

tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio)


#tit = 'Total Explained Variance = 99.48'


fig <- plot_ly(components,
               x = ~PC1, y = ~PC2, z = ~PC3,
               color = ~pData(normalized.data_GSE162785)$condition,
               #symbol = ~pData(normalized.data_GSE176190)$genotype,
               colors = brewer.pal(n = 8, name = "Set1"),  # personalize as cores
               symbols = c('circle', 'square', 'diamond', 'x', 'cross', 'triangle-up')) %>%  # personalize os símbolos
  add_markers(marker = list(size = 10)) 



fig <- fig %>%
  
  layout(
    
    # title = tit,
    
    scene = list(bgcolor = "#e5ecf6")
    
  )


fig

#Boxplot
boxplot(normalized.data_GSE162785)





# Filtro de intensidade 

media_matrizcontagem <- rowMeans(GSE162785)

hist_matriz_contagem <- hist(media_matrizcontagem, 100, col = "azure3", freq = FALSE, 
                             main = "Histogram of the mean intensities", 
                             border = "antiquewhite4",
                             xlab = "Mean intensities")

library(genefilter)
limite_histo <- shorth(media_matrizcontagem) 


hist_res_PC <- hist(media_matrizcontagem, 100, col = "azure3", freq = FALSE, 
                    main = "Histogram of the mean intensities with cutoff line",
                    border = "antiquewhite4",
                    xlab = "Mean intensities")

abline(v = limite_histo, col = "coral4", lwd = 2)


filtered.set <- GSE162785[media_matrizcontagem >=limite_histo,] # número de sondas que passaram pelo limite de corte mais curto
print(paste('Amostras antes do Filtro : ', dim(GSE162785)[1]))
print(paste('Amostras depois do Filtro : ', dim(filtered.set)[1]))


no_of_samples_PC <- table(paste0(pData(normalized.data_GSE162785)$condition)) # number of samples
no_of_samples_PC

samples_cutoff_PC <- min(no_of_samples_PC) # select the group with minimum number of samples (54 in our case) 

idx_man_threshold_PC <- apply(GSE162785, 1,function(x){sum(x > limite_histo) >=samples_cutoff_PC})
table(idx_man_threshold_PC) # número de sondas que passaram pelo corte curto e em pelo menos 50% das amostras

rmaData_filtrado <- subset(normalized.data_GSE176190, idx_man_threshold_PC) # filtrar as sondas que não atendem aos critérios
print(paste('Amostras depois do filtro de intensidade : ', dim(rmaData_filtrado)[1]))


# checking 

media_matrizcontagem2 <- rowMeans(Biobase::exprs(rmaData_filtrado))

hist_matriz_contagem <- hist(media_matrizcontagem2, 100, col = "azure3", freq = FALSE, 
                             main = "Histogram of the mean intensities", 
                             border = "antiquewhite4",
                             xlab = "Mean intensities")

hist(as.numeric(unlist(Biobase::exprs(rmaData_filtrado))),
     breaks = 100,
     col = "steelblue",
     main = "Distribuição das intensidades normalizadas (RMA)",
     xlab = "Intensidade log2",
     ylab = "Frequência")



boxplot(rmaData_filtrado)

GSE162785_filt <- as.data.frame(exprs(rmaData_filtrado))

# trancrever o codigo da matrix d exxpressao para nome do gene

library(stringr)

platform_GSE162785 <- platform_GSE162785 %>%
  mutate(
    Gene_symbol = sapply(str_split(gene_assignment, " // "), function(x) ifelse(length(x) > 1, x[2], NA))
  )


# Verificando  
table(is.na(platform_GSE162785$Gene_symbol) & platform_GSE162785$gene_assignment == '---')
View(platform_GSE162785[platform_GSE162785$gene_assignment == '---',])


# juntar na matriz de expressao o nome do gene e o id (se ja tiver filtrado, usar o filtrado)

ds_norm <- merge(
  GSE162785,                       # A matriz de contagem como dataframe
  platform_GSE162785[, c("ID", "Gene_symbol")], # Selecionar apenas as colunas relevantes do dataframe de símbolos
  by.x = "row.names",
  by.y = 'ID', # Unir pelo "gene_id"
  all.x = TRUE                           # Preservar todas as linhas da matriz de contagem
)

#organizacao das colunas
ds_norm <- ds_norm[, c(1,44,2:43)]
#vizualizacao
head(ds_norm)

# transformar os dados de expressao em numerico, para nao ler como caracter
ds_norm <- ds_norm %>%
  mutate(across(-c(1,2), ~ as.numeric(as.character(.))))

#eliminar os NA
ds_norm <- ds_norm[!is.na(ds_norm$Gene_symbol),]
#eliminar os vazios
ds_norm <- ds_norm[ds_norm$Gene_symbol != '',]

# eliminar dados repitidos, selecionnar o gene com maior variancia

# Calcular a variância de cada linha (ignorando a primeira coluna)
ds_norm$variancia <- apply(ds_norm[,-c(1,2)], 1, var, na.rm = TRUE)

# Ordenar do maior para o menor valor de variância
ds_norm <- ds_norm[order(-ds_norm$variancia), ]

#Dessa forma, pode-se eliminar as duplicatas e pegar aquelas com a maior variância através da funcao duplicated, pois os dados já foram ordenados do maior para o menor
ds_norm <- ds_norm[!duplicated(ds_norm$Gene_symbol ), ]

# Verificar se ainda existe duplicata
table(duplicated(ds_norm$Gene_symbol))

# colocamos os nomes dos genes como linhas da matriz
rownames(ds_norm) <- ds_norm$Gene_symbol

#excluir as colunas ja tratadas
ds_norm[, c(1,2)] <- NULL

ds_norm$variancia <- NULL

# salvar a matrix de expressao formatada
write.csv(ds_norm, '/media/daner/DATA/Arquivos leituras/Amostras/GSE27524/Counts_GSE27524_filt.csv', quote = FALSE, row.names = TRUE)



###############################################################################################################################################################

library(limma)

# separamos os grupos
Groups <- unlist(metadata_GSE162785$condition)
#criar matriz relacionada a separacao dos grupos
design <- model.matrix(~factor(Groups) -1)
design

#salvar os grupos como objeto fatores
fatores <- levels(as.factor(Groups))

#substituir / e - por _
fatores <- gsub("[-/]", "_", fatores)


#trocar nome das colunas
colnames(design) <- fatores
design
#aplicar o modelo linear baseado na matriz de expressao e design
fit <- lmFit(ds_norm, design)

# comparacoes entre os grupos das duas matrrize
colnames(design) 
constrast.matrix <- makeContrasts('Doxorubicin-control','FK228-control',
                                  levels = design)
fitc <- contrasts.fit(fit, constrast.matrix)
fitc <- eBayes(fitc)
#mostra o resuno dos resultados
summary(decideTests(fitc, method="global"))

#mostra os resultados com os filtros dos coeficientes
Doxorubicin_control <- topTable(fitc, number=Inf, p.value=0.05, lfc = 0.58, coef=1)
FK228_control <- topTable(fitc, number=Inf, p.value=0.05, lfc = 0.58, coef=2)

# mostrar os resultados de todos os genes 
Doxorubicin_control_sf <- topTable(fitc, number=Inf, coef=1)
FK228_control_sf <- topTable(fitc, number=Inf, coef=2)

# salvar dados
write.csv(Doxorubicin_control, '/media/daner/DATA/Arquivos leituras/Amostras/GSE27524/Counts_GSE27524_filt.csv', quote = FALSE, row.names = TRUE)
write.csv(FK228_control, '/media/daner/DATA/Arquivos leituras/Amostras/GSE27524/Counts_GSE27524_filt.csv', quote = FALSE, row.names = TRUE)
write.csv(Doxorubicin_control_sf, '/media/daner/DATA/Arquivos leituras/Amostras/GSE27524/Counts_GSE27524_filt.csv', quote = FALSE, row.names = TRUE)
write.csv(FK228_control_sf, '/media/daner/DATA/Arquivos leituras/Amostras/GSE27524/Counts_GSE27524_filt.csv', quote = FALSE, row.names = TRUE)

###############################################################################################################################################################################################

# Enriquecimento funcional


# utilzando dados sem filtro

library('fgsea')
gene_sets <- gmtPathways('h.all.v2024.1.Hs.symbols.gmt')


stats <- Doxorubicin_control_sf$logFC
names(stats) <- rownames(Doxorubicin_control_sf)
stats <- stats[is.finite(stats)]
stats <- sort(stats, decreasing = TRUE)

# Executar fgsea
fgsea_res <- fgsea(pathways = gene_sets, 
                   stats = stats, 
                   eps = 0.0,
                   minSize = 15,
                   maxSize = 500)


