# Configuração global para análise de RNA-seq
config <- list(
    # Parâmetros de filtragem
    filter = list(
        min_count = 10,
        min_samples = 3,
        variance_cutoff = 0.1
    ),
    
    # Parâmetros de normalização
    normalization = list(
        method = "rma",  # ou "quantile"
        log2_transform = TRUE
    ),
    
    # Parâmetros para análise diferencial
    diff_expression = list(
        p_value_cutoff = 0.05,
        log2fc_cutoff = 0.58
    ),
    
    # Parâmetros para enriquecimento
    enrichment = list(
        min_size = 15,
        max_size = 500,
        p_value_cutoff = 0.05
    ),
    
    # Diretórios de saída
    output_dirs = list(
        figures = "output/figures",
        results = "output/results",
        qc = "output/qc"
    )
)