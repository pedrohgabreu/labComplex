# Funções utilitárias para análise de RNA-seq

#' Função para carregar e pré-processar dados
#' @param geo_id ID do GEO dataset
#' @param config Lista de configurações
#' @return Lista com dados processados e metadados
load_and_preprocess <- function(geo_id, config) {
    # Carregar dados do GEO
    meta <- getGEO(geo_id, GSEMatrix = TRUE)
    metadata <- pData(meta[[1]])
    
    # Carregar plataforma
    platforms <- lapply(meta, annotation)
    platform_id <- platforms[[1]]
    platform_info <- getGEO(platform_id)
    platform_data <- Table(platform_info)
    
    return(list(
        metadata = metadata,
        platform = platform_data
    ))
}

#' Função para controle de qualidade
#' @param expr_data Matriz de expressão
#' @param metadata Metadados
#' @param config Configurações
#' @return Lista com resultados do QC
perform_qc <- function(expr_data, metadata, config) {
    # Boxplot
    png(file.path(config$output_dirs$figures, "boxplot.png"))
    boxplot(expr_data, main="Distribution of Expression Values")
    dev.off()
    
    # Histograma
    png(file.path(config$output_dirs$figures, "histogram.png"))
    hist(as.numeric(unlist(expr_data)), 
         breaks=100,
         main="Distribution of Expression Values",
         xlab="Expression")
    dev.off()
    
    # PCA
    pca_res <- prcomp(t(expr_data))
    png(file.path(config$output_dirs$figures, "pca.png"))
    plot(pca_res$x[,1:2], 
         col=as.factor(metadata$condition),
         main="PCA Plot")
    dev.off()
    
    return(pca_res)
}

#' Função para filtrar dados
#' @param expr_data Matriz de expressão
#' @param config Configurações
#' @return Matriz filtrada
filter_data <- function(expr_data, config) {
    # Filtro por contagem mínima
    keep <- rowSums(expr_data) > config$filter$min_count
    expr_filtered <- expr_data[keep,]
    
    # Filtro por número mínimo de amostras
    keep <- rowSums(expr_filtered >= config$filter$min_count) >= config$filter$min_samples
    expr_filtered <- expr_filtered[keep,]
    
    return(expr_filtered)
}

#' Função para análise de expressão diferencial
#' @param expr_data Matriz de expressão
#' @param metadata Metadados
#' @param config Configurações
#' @return Resultados da análise diferencial
differential_expression <- function(expr_data, metadata, config) {
    # Criar design
    design <- model.matrix(~metadata$condition)
    
    # Fit linear model
    fit <- lmFit(expr_data, design)
    fit <- eBayes(fit)
    
    # Get results
    results <- topTable(fit, 
                       number=Inf,
                       p.value=config$diff_expression$p_value_cutoff,
                       lfc=config$diff_expression$log2fc_cutoff)
    
    return(results)
}

#' Função para análise de enriquecimento
#' @param de_results Resultados da expressão diferencial
#' @param config Configurações
#' @return Resultados do enriquecimento
perform_enrichment <- function(de_results, config) {
    # Preparar estatísticas
    stats <- de_results$logFC
    names(stats) <- rownames(de_results)
    stats <- stats[is.finite(stats)]
    stats <- sort(stats, decreasing = TRUE)
    
    # Carregar gene sets
    gene_sets <- gmtPathways('h.all.v2024.1.Hs.symbols.gmt')
    
    # Executar fgsea
    fgsea_res <- fgsea(pathways = gene_sets,
                       stats = stats,
                       eps = 0.0,
                       minSize = config$enrichment$min_size,
                       maxSize = config$enrichment$max_size)
    
    return(fgsea_res)
}