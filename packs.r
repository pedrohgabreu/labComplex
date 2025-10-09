if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Passo 2: Definir a lista de pacotes a serem instalados

# Lista de pacotes do CRAN (repositório padrão do R)
cran_packages <- c(
    "tidyverse",         # Inclui dplyr, ggplot2, stringr, etc.
    "openxlsx",          # Para ler e escrever arquivos Excel
    "pheatmap",          # Para criar heatmaps bonitos
    "RColorBrewer",      # Paletas de cores
    "plotly",            # Gráficos 3D e interativos
    "ggfortify",         # Facilita o uso do ggplot2 com outros objetos
    "Matrix"             # Pacote essencial para matrizes, que deu problema antes
)

# Lista de pacotes do Bioconductor (repositório de bioinformática)
bioc_packages <- c(
    "affy",              # A alternativa que encontramos para o 'oligo' para seus dados
    "GEOquery",          # Para baixar dados do repositório GEO
    "limma",             # Essencial para análise de expressão diferencial
    "arrayQualityMetrics", # Para controle de qualidade de microarrays
    "PCAtools",          # Ferramentas para análise de PCA
    "genefilter",        # Filtragem de genes
    "fgsea",             # Análise de enriquecimento de conjuntos de genes
    "DESeq2",            # Para análise de RNA-Seq
    "edgeR",             # Alternativa para análise de RNA-Seq
    "SummarizedExperiment", # Estrutura de dados central no Bioconductor
    "clusterProfiler",   # Para análise de enriquecimento funcional
    "preprocessCore"     # Dependência importante para normalização
    # "oligo"            # O pacote original que causou o erro de compilação. Deixei comentado.
)


# Passo 3: Instalar os pacotes

# Instala os pacotes do CRAN
cat("\n--- Instalando pacotes do CRAN ---\n")
install.packages(cran_packages)

# Instala os pacotes do Bioconductor
cat("\n--- Instalando pacotes do Bioconductor ---\n")
BiocManager::install(bioc_packages)