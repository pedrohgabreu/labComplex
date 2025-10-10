# Script principal para análise de RNA-seq

#=============================================================================
# INSTALAÇÃO DE DEPENDÊNCIAS DO SISTEMA (executar no terminal antes do R)
#=============================================================================
# No AnduinOS/Ubuntu, execute estes comandos no terminal primeiro:
# sudo apt update
# sudo apt install build-essential r-base-dev
# sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev
# sudo apt install libgit2-dev libharfbuzz-dev libfribidi-dev
# sudo apt install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
# sudo apt install libfontconfig1-dev libcairo2-dev
# sudo apt install cmake libgsl-dev

#=============================================================================
# CONFIGURAÇÃO DO BIOCONDUCTOR E INSTALAÇÃO DE PACOTES
#=============================================================================

# Instalar BiocManager com versão correta para R 4.4
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.20")

# Pacotes do CRAN necessários
cran_pkgs <- c(
    "tidyverse",        # Manipulação e visualização de dados
    "pheatmap",         # Heatmaps
    "PCAtools",         # Análise de PCA
    "RColorBrewer",     # Paletas de cores
    "VennDiagram",      # Diagramas de Venn
    "corrplot",         # Gráficos de correlação
    "ggrepel",          # Labels sem sobreposição
    "plotly",           # Gráficos interativos
    "DT",              # Tabelas interativas
    "openxlsx",        # Leitura/escrita Excel
    "readr",           # Leitura de arquivos
    "dplyr",           # Manipulação de dados
    "ggplot2",         # Visualização
    "stringr",         # Manipulação de strings
    "forcats",         # Fatores
    "purrr",           # Programação funcional
    "tibble",          # Data frames modernos
    "tidyr"            # Organização de dados
)

# Pacotes do Bioconductor necessários
bio_pkgs <- c(
    "limma",                # Análise de expressão diferencial
    "edgeR",               # RNA-seq analysis
    "DESeq2",              # RNA-seq differential expression
    "fgsea",               # Gene set enrichment analysis
    "clusterProfiler",     # Análise de enriquecimento
    "org.Hs.eg.db",       # Anotação genômica humana
    "DOSE",                # Disease ontology
    "pathview",            # Visualização de pathways
    "EnhancedVolcano",     # Volcano plots
    "ComplexHeatmap",      # Heatmaps complexos
    "GenomicFeatures",     # Features genômicas
    "AnnotationDbi",       # Interface para anotações
    "GO.db",               # Gene Ontology
    "KEGG.db",             # KEGG pathways
    "ReactomePA",          # Reactome pathway analysis
    "msigdbr",             # MSigDB gene sets
    "biomaRt",             # Interface Ensembl/BioMart
    "ArrayExpress",        # Acesso ao ArrayExpress
    "GEOquery",            # Acesso ao GEO
    "Biobase",             # Classes básicas Bioconductor
    "SummarizedExperiment", # Estruturas de dados
    "GenomicRanges",       # Ranges genômicos
    "IRanges",             # Integer ranges
    "S4Vectors",           # Estruturas S4
    "BiocGenerics",        # Funções genéricas
    "preprocessCore",      # Normalização
    "affy",                # Análise Affymetrix (se necessário)
    "oligo"                # Análise oligonucleotide arrays
)

# Instalar pacotes CRAN
cat("Instalando pacotes CRAN...\n")
install.packages(cran_pkgs, dependencies = TRUE)

# Instalar pacotes Bioconductor
cat("Instalando pacotes Bioconductor...\n")
BiocManager::install(bio_pkgs, update = TRUE, ask = FALSE)

# Verificar instalação
cat("Verificando instalação...\n")
BiocManager::valid()

#=============================================================================
# CARREGAR BIBLIOTECAS NECESSÁRIAS
#=============================================================================

# Manipulação e visualização de dados
library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(forcats)
library(purrr)
library(tibble)
library(tidyr)

# Análise de expressão diferencial
library(limma)
library(edgeR)
library(DESeq2)

# Análise de enriquecimento
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(ReactomePA)
library(msigdbr)

# Visualização
library(pheatmap)
library(PCAtools)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(RColorBrewer)
library(VennDiagram)
library(corrplot)
library(ggrepel)
library(plotly)

# Acesso a dados
library(ArrayExpress)
library(GEOquery)
library(biomaRt)

# Estruturas de dados Bioconductor
library(Biobase)
library(SummarizedExperiment)
library(GenomicRanges)
library(GenomicFeatures)
library(AnnotationDbi)

# Utilitários
library(DT)
library(openxlsx)

#=============================================================================
# CONFIGURAÇÕES GERAIS
#=============================================================================

# Configurar diretório de trabalho
setwd("~/analise_rnaseq")

# Criar diretórios se não existirem
dirs_to_create <- c("data", "results", "plots", "reports")
for(dir in dirs_to_create) {
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
        cat("Diretório criado:", dir, "\n")
    }
}

# Configurações de plot
theme_set(theme_bw())
options(ggplot2.discrete.colour = scale_colour_viridis_d)
options(ggplot2.discrete.fill = scale_fill_viridis_d)

#=============================================================================
# CARREGAR CONFIGURAÇÕES E FUNÇÕES PERSONALIZADAS
#=============================================================================

# Carregar configurações e funções (se existirem)
if(file.exists("config.R")) {
    source("config.R")
    cat("Configurações carregadas de config.R\n")
}

if(file.exists("functions.R")) {
    source("functions.R") 
    cat("Funções personalizadas carregadas de functions.R\n")
}

# Criar diretórios de saída se config.R existir
if(exists("config") && "output_dirs" %in% names(config)) {
    for (dir in config$output_dirs) {
        if (!dir.exists(dir)) {
            dir.create(dir, recursive = TRUE)
        }
    }
}

#=============================================================================
# PIPELINE PRINCIPAL DE ANÁLISE
#=============================================================================

cat("=== INICIANDO ANÁLISE DE RNA-SEQ ===\n")
cat("R versão:", R.version.string, "\n")
cat("Bioconductor versão:", BiocManager::version(), "\n\n")

# 1. Carregar e pré-processar dados
cat("1. Carregando dados...\n")
if(exists("load_and_preprocess") && file.exists("config.R")) {
    data <- load_and_preprocess("GSE162785", config)
} else {
    cat("   Função load_and_preprocess não encontrada. Defina em functions.R\n")
}

# 2. Realizar controle de qualidade
cat("2. Realizando controle de qualidade...\n")
if(exists("perform_qc") && exists("data")) {
    qc_results <- perform_qc(data$expr_data, data$metadata, config)
} else {
    cat("   Função perform_qc não encontrada. Defina em functions.R\n")
}

# 3. Filtrar dados
cat("3. Filtrando dados...\n")
if(exists("filter_data") && exists("data")) {
    filtered_data <- filter_data(data$expr_data, config)
} else {
    cat("   Função filter_data não encontrada. Defina em functions.R\n")
}

# 4. Análise de expressão diferencial
cat("4. Realizando análise de expressão diferencial...\n")
if(exists("differential_expression") && exists("filtered_data")) {
    de_results <- differential_expression(filtered_data, data$metadata, config)
} else {
    cat("   Função differential_expression não encontrada. Defina em functions.R\n")
}

# 5. Análise de enriquecimento
cat("5. Realizando análise de enriquecimento...\n")
if(exists("perform_enrichment") && exists("de_results")) {
    enrichment_results <- perform_enrichment(de_results, config)
} else {
    cat("   Função perform_enrichment não encontrada. Defina em functions.R\n")
}

# 6. Salvar resultados
cat("6. Salvando resultados...\n")
if(exists("de_results")) {
    output_dir <- if(exists("config")) config$output_dirs$results else "results"
    write.csv(de_results, 
              file.path(output_dir, "differential_expression.csv"),
              row.names = FALSE)
    cat("   Resultados DE salvos em:", file.path(output_dir, "differential_expression.csv"), "\n")
}

if(exists("enrichment_results")) {
    output_dir <- if(exists("config")) config$output_dirs$results else "results"
    write.csv(enrichment_results, 
              file.path(output_dir, "enrichment_results.csv"),
              row.names = FALSE)
    cat("   Resultados de enriquecimento salvos em:", file.path(output_dir, "enrichment_results.csv"), "\n")
}

cat("\n=== ANÁLISE COMPLETA! ===\n")

# Informações da sessão para reprodutibilidade
cat("\nInformações da sessão R:\n")
sessionInfo()