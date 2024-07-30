## Author : Ilango Guy
## Contact : guy.ilango@univ-tours.fr
## Doc : usefull function to source for R analysis


library(dplyr)
library(stringr)

prep_GO <- function(df){
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- rownames(df)

# omit any NA values 
gene_list<- original_gene_list

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
return(gse)
}

baranal <- function(tab , frequence , sample , cluster , curvetype){
    tab %>% ggplot(aes(y = frequence, x = sample, fill = as.character(cluster))) +
  geom_flow(aes(alluvium = cluster), alpha= .5, color = "white",
            curve_type = curvetype, 
            width = .5) +
  geom_col(width = .5, color = "white") +
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(panel.grid.major = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())}

heatmap_bulk_zscore <- function(dds , title){
library(ClassDiscovery)
ntd <- normTransform(dds)
t <- assay(ntd)
zscores <- apply(t, 1, function(gene_expr) {
  (gene_expr - mean(gene_expr)) / sd(gene_expr)
})
zscores <- zscores[ , colSums(is.na(zscores)) == 0]
hc.features <- t(zscores) %>% t() %>%  
  distanceMatrix(metric="pearson") %>%
  hclust(method="average")               
 
hc.samples <- t(zscores) %>% 
  distanceMatrix(metric="pearson") %>%
  hclust(method="average")
 
return(pheatmap::pheatmap(t(zscores), 
                   scale = "row",
                   
                   show_colnames  = T
                   ,main = title, 
                   cluster_rows   = hc.features,
                   cluster_cols   = hc.samples) )
}


save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}





sc_annotate_de <- function(df) {
    df$gene <- rownames(df)
  df <- mutate(df, annot = case_when(
      str_detect(gene, '^Gm') & avg_log2FC > 0 & p_val_adj <= 0.05 ~ 'pas_interessant_up_sig',
      str_detect(gene, '^Gm') & avg_log2FC < 0 & p_val_adj <= 0.05 ~ 'pas_interessant_down_sig',
      str_detect(gene, '^Gm') &   p_val_adj > 0.05 ~ 'pas_interessant_unsig',
      str_detect(gene, '^mt') & avg_log2FC > 0 & p_val_adj <= 0.05 ~ 'pas_interessant_up_sig',
      str_detect(gene, '^mt') & avg_log2FC < 0 & p_val_adj <= 0.05 ~ 'pas_interessant_down_sig',
      str_detect(gene, '^mt') &   p_val_adj > 0.05 ~ 'pas_interessant_unsig',
      str_detect(gene, '^Rp[l|s]') & avg_log2FC > 0 & p_val_adj <= 0.05 ~ 'pas_interessant_up_sig',
      str_detect(gene, '^Rp[l|s]') & avg_log2FC < 0 & p_val_adj <= 0.05 ~ 'pas_interessant_down_sig',
      str_detect(gene, '^Rp[l|s]') &  p_val_adj > 0.05 ~ 'pas_interessant_unsig',
      str_detect(gene, '^Mir') & avg_log2FC > 0 & p_val_adj <= 0.05 ~ 'pas_interessant_up_sig',
      str_detect(gene, '^Mir') & avg_log2FC < 0 & p_val_adj <= 0.05 ~ 'pas_interessant_down_sig',
      str_detect(gene, '^Mir') &   p_val_adj > 0.05 ~ 'pas_interessant_unsig',
     avg_log2FC > 0 & p_val_adj <= 0.05 ~ 'up_sig',
      avg_log2FC < 0 & p_val_adj <= 0.05 ~ "down_sig",
      p_val_adj > 0.05 ~ "unsig",
      
      
         
     
)) 
}

sc_volcano_plotly <- function(df,title) {
  plot_ly(data = df, x = df$avg_log2FC, y = -log10(df$p_val_adj), text = rownames(df), mode = "markers", color = df$annot) %>%   layout(title = title)

 
}
sc_volcano <- function(df){
ggplot(df , aes( x = df$avg_log2FC , y = -log10(df$p_val_adj) , color = df$annot)) + geom_point()}



bulk_volcano_plotly <- function(df,title) {
  plot_ly(data = df, x = df$log2FoldChange, y = -log10(df$pvalue), text = rownames(df), mode = "markers", color = df$annot) %>%   layout(title = title)

 
}
bulk_volcano <- function(df){
ggplot(df , aes( x = df$log2FoldChange , y = -log10(df$pvalue) , color = df$annot)) + geom_point()}


make_it_bulk <- function(seurat_obj , sampleID){
    seurat_obj@meta.data$bulk <- seurat_obj@meta.data$sampleID
    Idents(seurat_obj) <- seurat_obj@meta.data$bulk
    markers <- FindAllMarkers(seurat_obj)
    return(markers)
    }

make_readable <- function(df){
df <- dplyr::filter(df , !grepl('Rik', rownames(df)))
df <- dplyr::filter(df , !grepl('^mt', rownames(df)))
df <- dplyr::filter(df , !grepl('^Rp[l|s]', rownames(df)))
df <- dplyr::filter(df , !grepl('^Mir', rownames(df)))
    df <- dplyr::filter(df , !grepl('^Gm[0-9]+', rownames(df)))
}

bulk_annotate_de <- function(df) {
  df <- mutate(df, annot = case_when(
    log2FoldChange > 0 ~ 'up',
    log2FoldChange < 0 ~ "down"
     
)) 
}


basic_pipe <- function(seurat_obj , dim , res ){
    seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

seurat_obj <- RunUMAP(seurat_obj, dims = 1:dim)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dim)
seurat_obj <- FindClusters(seurat_obj , resolution = res)
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



##### signatures


th17 <- c("Mir873a","Ccr6","Il2","Il4","Il6","Il6ra","Irf4","Loxl3","Ly9","Smad7","Ascl2","Rora","Rorc","Ccl20","Foxp3","Stat3","Nlrp3","Tgfb1","Zbtb7b","Zc3h12a","Malt1","Tnfsf18","Nfkbid","Slamf6","Rc3h2","Rc3h1","Batf","Otud5","Tbx21","Mir301","Mir326","Nfkbiz","Il23a","Entpd7")

Tissue_residency <- (c("2900026A02Rik","AA467197","Abcb1b","Adam19","Atf3","B4galnt4","Bag3","Bhlhe40","Btg2","Btg3","Ccl4","Cd244","Cd69","Cdh1","Chd7","Chn2","Cish","Coq10b","Cpd","Crem","Ctla4","Dennd4a","Dgat1","Dnaja1","Dnaja4","Dnajb1","Dnajb4","Dnajb6","Dusp1","Dusp4","Dusp6","Egr1","Ehd1","Ell2","Fasl","Fbxo30","Fgl2","Fmnl3","Fos","Fosl2","Gabarapl1","Gadd45b","Gem","Glrx","Gpr171","Gzmb","H3f3b","Havcr2","Hip1","Hpgds","Hsp90aa1","Hspa5","Icos","Ifng","Ifrd1","Il21r","Il4ra","Inpp4b","Irf4","Isg20","Isy1","Itga1","Jun","Junb","Kdm6b","Litaf","Lmnb1","Ly6g5b","Mxd1","Mxi1","Neurl3","Nfil3","Nfkbie","Nr4a1","Nr4a2","Nr4a3","Ntan1","Odc1","P2ry10","P4hb","Pdcd1","Per1","Plk3","Plscr1","Pmepa1","Ppp1r15a","Prdx6","Ptp4a1","Qpct","Rbpj","Rgs1","Rgs10","Rgs16","Rhob","Rnf149","Rrad","Sik1","Skil","Slc16a6","Slc7a5","Spty2d1","Ssbp2","Stk17b","Tgif1","Tigit","Tiparp","Tjp1","Tnf","Tnfaip3","Tnfrsf1b","Tnfrsf9","Tra2a","Trib1","Trp53inp2","Ube2s","Vps37b","Wsb1","Xcl1","Zfand5","Zfp36l1"))

Tissue_repair1 <- (c("Flg","Flg2","Furin","Areg","Bmp2","Bmp4","Bmp7","Ngf","Tgfa","Tgfb1","Tgfb2","Tgfb3","Pdgfa","Pdgfb","Csf2","Fgf1","Fgf7","Fgf10","Fgf22","Igf1","Igf2","Vegfa","Ctgf","Mmp3","Mmp10","Mmp13","Mmp15","Mmp25","Mmp28","Defb1","Defb6","Epgn","Hbegf","Egf","Mmp2","Tnf","Il1a","Il1b","Csf1","Ccl2","Ccl3","Cxcl1","Cxcl2","Cxcl10","Cxcl15","Il6","Pgf","Angpt1","Angpt2","Angpt4","Hgf","Cyr61","Inhba","Inhbb","Btc","Ereg","Nrg1","Nrg2","Nrg4","Nrg3","Lep","Vegfb","Vegfc","Vegfd","Pdgfc","Hmgb1","Ptges2","Cxcl12","Shh","Dhh","Ihh","Disp1","Dll1","Dll3","Dll4","Jag1","Jag2","Wnt1","Wnt2","Wnt2b","Wnt3","Wnt3a","Wnt4","Wnt5a","Wnt5b","Wnt6","Wnt7a","Wnt7b","Wnt8a","Wnt8b","Wnt9a","Wnt9b","Wnt10a","Wnt10b","Wnt11","Wnt16","Chat","Thbs1","Hif1a"))

Tissue_repair2 <- (c("Pcyt1a","Pcsk1","Ccr1","Areg","Frmd5","Havcr2","Slc15a3","Ccr3","Il10rb","Lyn","Npnt","Ctsh","Arnt2","Il23r","Tubb6","Neb","Snx9","Sik1","Rgs2","Dgat2","Fam46a","Trf","Kcna4","Arl5a","Tnfrsf10b","Padi2","Gpr55","Sdc4","Cd200r1","Pparg","Camk2n1","Ehd4","Tmbim1","Ifrd1","Serpinb6a","Plxnd1","Rrad","Rasgef1b","Plin2","Msrb3","Per1","Sgms1","Dusp14","Spty2d1","Myo1d","Coq10b","Il9r","Evi5","Cables1","Cd80","Cpm","Ehd1","Nr1d1","Alcam","Epcam","Arl5b","Ets2","Acot11","Nebl","Klf4","Klf8","Slc25a19","Cnnm2","L1cam","Gprc5a","Sytl3","Tagln2","Axl","Zdhhc23","Ky","Ttbk2","Fam129b","Trim46","Gla","Gm2a","Cd74","Clic4","Slc16a3","Dusp3","Pdgfb","Itgav","Tdpoz4","Rnf125","Lgmn","Crem","Tnfrsf9","Rassf6","Bmpr1a","N4bp1","Rab11fip1","Zc3h12c","Got1","Hpse","Dot1l","Zfp36","Icam1","Rab4a","Plxdc1","Rorc","Asns","Irak3","Kdm6b","Ero1l","Fos"))

Tissue_repair3 <- (c("Il18r1","Il17re","Itgb1","Ramp3","Tnfsf11","Cd40lg","Plcb4","B3gnt5","Cep290","Tmem176b","Plxnc1","Tmem176a","Lif","Nrp1","Itga7","Colq","Sema4f"))

Tissue_repair4 <- (c("Adrb2","Aimp1","Aqp3","Arhgap1","Arnt","B4galt1","Bcl11b","Bub3","Cap2","Casp1","Cd151","Cd9","Cdkn1b","Cebpg","Csf2","Ctnnb1","Epb41","Fam129b","Fgfr1op2","Fhl2","Flna","Hmox2","Icam1","Icos","Il22","Il33","Ilk","Itga3","Itga9","Itgb1","Jun","Junb","Klf10","Klf4","Klk8","Lgals1","Map3k1","Map3k5","Mfge8","Mgat5","Mif","Myc","Myd88","Ncoa3","Nf1","Nfe2l2","Nfic","Nr4a1","Panx1","Pkd1","Ppard","Prkca","Prkce","Prkch","Ptk2","Ptk2b","Rac1","Rac2","Raf1","Ramp1","Ripk3","Sdc1","Sdc4","Siah2","Smad4","Smad7","Sod1","Stat3","Tcf3","Tgfb1","Thra","Thy1","Tnf","Tnfsf14","Trp53","Vav1","Vav3","Vegfa","Vim","Slc6a6","Itgb2","Hif1a","Clic4","Akt1","Anxa1","Dst","Egr1","Kit","Plaur","Pten","Il4","Ncoa6","Cd44"))

Circulating <- (c("1700025G04Rik","A430078G23Rik","Acss1","Aff3","Apobec2","Arhgap26","Arhgef18","Armc7","As3mt","Aven","B3gnt5","Bin2","Car5b","Cdc14b","Cmah","D1Ertd622e","Dapl1","Dennd2d","Dgka","Dock5","Dusp7","Ehd3","Elovl7","Eomes","F2rl2","Fam49a","Fam53b","Fgf13","Flna","Gab3","Gramd4","Gzmm","Haao","Icam2","Il17ra","Il18rap","Kbtbd11","Kcnj8","Klf2","Klhdc1","Klrb1c","Klrg1","Lair1","Lfng","Limd1","Lpin1","Ly9","Mboat1","Ms4a4c","Nck2","Nfatc3","Nod1","Nsg2","Pde2a","Prss12","Prune1","Pycard","Qprt","Rap1gap2","Rasa3","Rasgrp2","Rgs19","S1pr1","S1pr4","Samd3","Sell","Sfxn3","Sh2d1a","Sh3bp5","Sidt1","Sike1","Slamf6","Smpdl3b","Spn","St3gal1","Stard10","Tcf7","Tlr1","Tmem71","Tnfrsf26","Traf3ip3","Ubxn2b","Usp33","Ypel1","Zeb2","Zfp595","Zfp760"))

iT_17 <- (c("Rorc","Itgae","Il23r","Serpinb1a","Mmp25","1700040D17Rik","Pxdc1","Il1r1","Lingo4","Ccr8"))

iT_cycling_S <- (c("Actb","Ptma","Slc25a5","Ppia","Hmgn2","Lgals1","Cfl1","Gapdh","Ran","H2afz","Ldha","Eif5a","Hmgb1","Npm1","Actg1","Ybx1","Hmgb1","Ranbp1","Set","Mif"))

iT_cycling_G2M <- (c("Tuba1b","Stmn1","Tubb5","Ptma","Pclaf","H2afz","Hmgb2","Hmgn2","Hmgb1","Cks1b","Top2a","Dut","Hnrnpab","Lmnb1","Dek","Ran","Ppia","Smsc2","Hist1h2ap","Ranbp1"))

TCR_signaling <- str_to_title(c("AKT1","AKT2","AKT3","ARAF","BRAF","CD247","CD28","CD3D","CD3E","CD3G","CD80","CD86","CDC42","CHUK","CSK","FOS","GRAP2","HRAS","IKBKB","ITPR1","JUN","LAT","LCK","LCP2","MAP2K1","MAP2K2","MAP3K1","MAPK1","MAPK3","MAPK8","MAPK9","NCK1","NCK2","NFKBIA","NRAS","PAK1","PAK2","PAK3","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PLCG1","PPP3CA","PPP3CB","PPP3CC","PRKCQ","PTPRC","RAC1","RAF1","SOS1","SOS2","VAV1","VAV2","VAV3","WAS","ZAP70","WAS","ZAP70"))

Fas_pathway <- (c("Arhgdib","Casp3","Casp6","Casp7","Casp8","Cflar","Daxx","Dffa","Dffb","Fadd","Faf1","Fas","Fasl","Jun","Lmna","Lmnb1","Lmnb2","Map2k4","Map3k1","Map3k7","Mapk8","Pak1","Pak2","Parp1","Prkdc","Ptpn13","Rb1","Sptan1"))

Il1R_signaling <- (c("Chuk","Ecsit","Ifna1","Ifnb1","Ikbkb","Il1a","Il1b","Il1r1","Il1rap","Il1rn","Il6","Irak1","Irak2","Irak3","Jun","Map2k3","Map2k6","Map3k1","Map3k14","Map3k7","Mapk14","Mapk8","Myd88","Nfkb1","Nfkbia","Rela","Tgfb1","Tgfb2","Tgfb3","Tnf","Traf6"))

ISG <- (c("Gbp5","I830012O16Rik","Ifit3","Ifi205","Ifit2","Cxcl10","Ifit1","Rsad2","Ifitm3","E030037K03Rik","Ifi204","Gbp3","E430029J22Rik","Iigp1","Tgtp1","Mx1","Slfn5","Tgtp1","Tnfsf10","Cmpk2","Slfn9","Mx2","Usp18","BC023105","Gbp2","Ifi44","Herc6","Ifi47","Isg20","AI607873","Ppm1k","Parp12","Gbp4","Fcgr1","Cd69","Gbp6","Pnp","Slfn1","Nlrc5","Mmp13","Nt5c3","Pnpt1","Irg1","Ms4a4c","Pyhin1","Oas3","Ifih1","Setdb2","BC094916","Oasl2","Nlrc5","Oas2","Oasl1","Pnp2","Eif2ak2","2010002M12Rik","D14Ertd668e","Ddx60","Tmem67","Irf7","Tagap","Ifi203","Rtp4","Nlrc5","Acer3","A530064D06Rik","Gpr141","Zufsp","Ms4a4a","Gdap10","Ly6i","Il18"))
