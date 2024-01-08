## Author : Ilango Guy
## Contact : guy.ilango@univ-tours.fr
## Doc : usefull function to source for R analysis


annotate_de <- function(df) {
  df <- mutate(df, annot = case_when(
    log2FoldChange > 0 ~ 'up',
    log2FoldChange < 0 ~ "down"
     
)) 
}


basic_pipe <- function(seurat_obj){
    seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj , resolution = 0.2)
return(seurat_obj) 
}


doublet_finder <- function(seurat_obj){ 
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Find significant PCs
  stdv <- seurat_obj[["pca"]]@stdev
  sum.stdv <- sum(seurat_obj[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  # finish pre-processing
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:min.pc)
  seurat_obj <- FindNeighbors(object = seurat_obj, dims = 1:min.pc)              
  seurat_obj <- FindClusters(object = seurat_obj, resolution = 0.2)
sweep.list <- paramSweep(seurat_obj, PCs = 1:min.pc, num.cores = detectCores() - 1)
sweep.stats <- summarizeSweep(sweep.list)
bcmvn <- find.pK(sweep.stats)
  
# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
## Homotypic doublet proportion estimate
annotations <- seurat_obj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp.poi <- round(optimal.pk * nrow(seurat_obj@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
# run DoubletFinder
seurat_obj <- doubletFinder_v3(seu = seurat_obj, 
                                   PCs = 1:min.pc, 
                                   pK = optimal.pk,
                                   nExp = nExp.poi.adj)
colnames(seurat_obj@meta.data)[8] <- "Doublet_Finder"
DimPlot(seurat_obj , group.by = "Doublet_Finder")    
seurat_obj.singlets <- subset(seurat_obj, Doublet_Finder == "Singlet") 
seurat_obj.singlets
seurat_obj.singlets    
return(seurat_obj.singlets)
}

plot_cinetic <- function(table1, cluster, table2){
  p<-ggplot(table1 %>% filter(clusters == cluster) , aes(x = timepoint  , y = freq   , group = 1) ) +geom_line( linewidth = 1)
  
  p<- p + annotate("text", x = "J4", y = as.numeric(table1 %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) -5, label = head(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[5] , color ="blue")
  p<- p + annotate("text", x = "J4", y = as.numeric(table1 %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) -4, label = head(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[4] , color ="blue")
  p <- p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) -3, label = head(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[3] , color ="blue")
  p<-p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) -2, label = head(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[2] , color ="blue")
  p<- p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) -1, label =head(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[1] , color ="blue")
  p<- p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) +5, label = tail(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[1] , color ="red")
  p<- p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters ==cluster , timepoint == "J4") %>% select(freq)) +4, label = tail(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[2] , color ="red")
  p<- p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) +3, label = tail(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[3] , color ="red")
  p<- p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) +2, label = tail(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[4] , color ="red")
  p<- p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) +1, label = tail(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[5] , color ="red")
  
  p<- p + annotate("text", x = "J7", y = as.numeric(table1 %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) -5, label = head(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[5] , color ="blue")
  p<- p + annotate("text", x = "J7", y = as.numeric(table1 %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) -4, label = head(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[4] , color ="blue")
  p <- p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) -3, label = head(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[3] , color ="blue")
  p<-p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) -2, label = head(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[2] , color ="blue")
  p<- p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) -1, label =head(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[1] , color ="blue")
  p<- p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) +5, label = tail(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[1] , color ="red")
  p<- p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters ==cluster , timepoint == "J7") %>% select(freq)) +4, label = tail(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[2] , color ="red")
  p<- p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) +3, label = tail(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[3] , color ="red")
  p<- p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) +2, label = tail(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[4] , color ="red")
  p<- p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) +1, label = tail(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[5] , color ="red")
  
  p<- p + annotate("text", x = "J12", y = as.numeric(table1 %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) -5, label = head(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[5] , color ="blue")
  p<- p + annotate("text", x = "J12", y = as.numeric(table1 %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) -4, label = head(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[4] , color ="blue")
  p <- p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) -3, label = head(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[3] , color ="blue")
  p<-p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) -2, label = head(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[2] , color ="blue")
  p<- p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) -1, label =head(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[1] , color ="blue")
  p<- p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) +5, label = tail(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[1] , color ="red")
  p<- p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters ==cluster , timepoint == "J12") %>% select(freq)) +4, label = tail(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[2] , color ="red")
  p<- p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) +3, label = tail(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[3] , color ="red")
  p<- p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) +2, label = tail(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[4] , color ="red")
  p<- p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) +1, label = tail(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[5] , color ="red")
  return(p)
  
}


dotplot_go_spe_mouse <- function(table , cluster){ 
cluster_down <- enrichGO(gene         = rownames(table %>% filter(avg_log2FC < 0)),
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

cluster_up <- enrichGO(gene         = rownames(table %>% filter( avg_log2FC > 0)),
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

a<-dotplot(cluster_down) + ggtitle(paste0("Cluster ", as.character(cluster) ," downregulated gene"))

b <- dotplot(cluster_up) + ggtitle(paste0("Cluster ",as.character(cluster) ," upregulated gene"))
return(list(a , b))

                                   }

dotplot_go_mouse <- function(table , cluster){ 
cluster_down <- enrichGO(gene         = rownames(table %>% filter(avg_log2FC < 0 , cluster == cluster) ),
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

cluster_up <- enrichGO(gene         = rownames(table %>% filter( avg_log2FC > 0, cluster == cluster)),
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

a<-dotplot(cluster_down) + ggtitle(paste0("Cluster ", as.character(cluster) ," downregulated gene"))

b <- dotplot(cluster_up) + ggtitle(paste0("Cluster ",as.character(cluster) ," upregulated gene"))
return(list(a , b))

                                   }

