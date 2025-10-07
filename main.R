# Script principal para análise de RNA-seq

# instalar packages do CRAN e Bioconductor
cran_pkgs <- c("tidyverse", "pheatmap", "PCAtools")
bio_pkgs <- c("limma", "BiocManager", "DESeq2", "edgeR", "fgsea")

install.packages(cran_pkgs)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(bio_pkgs)


# Carregar bibliotecas necessárias
library(tidyverse)
library(limma)
library(edgeR)
library(fgsea)
library(DESeq2)
library(pheatmap)
library(PCAtools)

# Carregar configurações e funções
source("config.R")
source("functions.R")

# Criar diretórios de saída
for (dir in config$output_dirs) {
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
    }
}

# 1. Carregar e pré-processar dados
cat("Carregando dados...\n")
data <- load_and_preprocess("GSE162785", config)

# 2. Realizar controle de qualidade
cat("Realizando controle de qualidade...\n")
qc_results <- perform_qc(data$expr_data, data$metadata, config)

# 3. Filtrar dados
cat("Filtrando dados...\n")
filtered_data <- filter_data(data$expr_data, config)

# 4. Análise de expressão diferencial
cat("Realizando análise de expressão diferencial...\n")
de_results <- differential_expression(filtered_data, data$metadata, config)

# 5. Análise de enriquecimento
cat("Realizando análise de enriquecimento...\n")
enrichment_results <- perform_enrichment(de_results, config)

# 6. Salvar resultados
cat("Salvando resultados...\n")
write.csv(de_results, 
          file.path(config$output_dirs$results, "differential_expression.csv"))
write.csv(enrichment_results, 
          file.path(config$output_dirs$results, "enrichment_results.csv"))

cat("Análise completa!\n")