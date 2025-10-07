hmap <- function(set, tset, col, scale, clustering, hcexRow, hcexCol, colsidecolors,
                 scex = .2, ggex = .5, main){
  require(gplots)
  require(dplyr)
  require(dendextend)
  # browser()
  if (clustering) {
    # Distancia, clustering y dendograma para los genes
    genes.distance <- as.dist(1 - cor(tset))
    genes.clusters <- hclust(genes.distance, method = "ward.D2")
    genes.dendrogr <- as.dendrogram(genes.clusters, h = 5)
    # Distancia, clustering y dendograma para las muestras
    samples.distance <- as.dist(1 - cor(set))
    samples.clusters <- hclust(samples.distance, method = "ward.D2")
    samples.dendrogr <- as.dendrogram(samples.clusters)
    # Change cex values
    samples.dendrogr <- samples.dendrogr %>% set("labels_cex", scex)
    genes.dendrogr   <- genes.dendrogr   %>% set("labels_cex", ggex) 
    # Representamos las muestras
    plot(samples.dendrogr)
    par(mar = c(2,10,2,10))
    # plot(samples.dendrogr, horiz = TRUE, cex = 0.2)
    # Representamos los genes
    # plot(genes.clusters, labels = F, cex = 0.4)
    # par(mar = c(2,10,2,10), cex = .2)
    # browser()
    plot(genes.dendrogr, horiz = T)
    # Reordenamos los genes
    fils <- nrow(set)
    cols <- ncol(set)
    genes.dendrogr <- reorder(genes.dendrogr, 1:fils)
    # Reordenamos las muestras
    samples.dendrogr <- reorder(samples.dendrogr, 1:cols)
    # Generamos el heatmap
    heatmap.2(set, Rowv = genes.dendrogr, 
              ColSideColors = colsidecolors,
              Colv = samples.dendrogr, 
              scale = scale, #srtCol = 0,
              cexRow = hcexRow,
              cexCol = hcexCol,
              col = col, lwd = 5,
              trace = 'none',
              keysize = 1,
              main = main,
              density.info = 'none')
  }
  else{
    # Generamos el heatmap
    heatmap.2(set, Rowv = NULL, 
              ColSideColors = colsidecolors,
              Colv = NULL, scale = 'row', #srtCol = 0,
              cexRow = hcexRow,
              cexCol = hcexCol,
              col = col, lwd = 5,
              trace = 'none',
              keysize = 1,
              main = main,
              density.info = 'none')
  }
}