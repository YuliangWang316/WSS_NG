library(Seurat)
library(patchwork)
library(dplyr)
library(cowplot)

pbmc.data<-Read10X(data.dir = "c:/Users/xjmik/Downloads/BRCA_NG/Data/")
pbmc.metadata<-read.csv("c:/Users/xjmik/Downloads/BRCA_NG/Whole_miniatlas_meta.csv",sep = ",",header = TRUE,row.names = 1)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", meta.data = pbmc.metadata)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 1.2)
pbmc <- RunUMAP(pbmc, dims = 1:30)
pbmc <- RunTSNE(pbmc, dims = 1:30)
Idents(pbmc)<-pbmc$celltype_major
remove(all.genes,pbmc.metadata,pbmc.data)
DimPlot(pbmc,raster = FALSE)
pbmc <- RunLDA(pbmc,labels = pbmc$celltype_major,reduction.name = "lda_1")
pbmc <- RunLDA(pbmc,labels = pbmc$celltype_minor,reduction.name = "lda_2")
pbmc <- RunLDA(pbmc,labels = pbmc$celltype_subset,reduction.name = "lda_3")
pbmc <- RunUMAP(pbmc, dims = 1:8,reduction = "lda_1",reduction.name = "lda1_umap")
pbmc <- RunTSNE(pbmc, dims = 1:8,reduction = "lda_1",reduction.name = "lda1_tsne")
pbmc <- RunUMAP(pbmc, dims = 1:28,reduction = "lda_2",reduction.name = "lda2_umap")
pbmc <- RunTSNE(pbmc, dims = 1:28,reduction = "lda_2",reduction.name = "lda2_tsne")
pbmc <- RunUMAP(pbmc, dims = 1:48,reduction = "lda_3",reduction.name = "lda3_umap")
pbmc <- RunTSNE(pbmc, dims = 1:48,reduction = "lda_3",reduction.name = "lda3_tsne")
DimPlot(pbmc,raster = FALSE,reduction = "tsne")
DimPlot(pbmc,raster = FALSE,reduction = "lda1_umap")
DimPlot(pbmc,raster = FALSE,reduction = "lda1_tsne")
Idents(pbmc)<-pbmc$celltype_minor
DimPlot(pbmc,raster = FALSE,reduction = "umap")
DimPlot(pbmc,raster = FALSE,reduction = "tsne")
DimPlot(pbmc,raster = FALSE,reduction = "lda2_umap")
DimPlot(pbmc,raster = FALSE,reduction = "lda2_tsne")
Idents(pbmc)<-pbmc$celltype_subset
DimPlot(pbmc,raster = FALSE,reduction = "umap")
DimPlot(pbmc,raster = FALSE,reduction = "tsne")
DimPlot(pbmc,raster = FALSE,reduction = "lda3_umap")
DimPlot(pbmc,raster = FALSE,reduction = "lda3_tsne")
Idents(pbmc)<-pbmc$celltype_major
new.cluster.ids <- c("Endothelial", "CAFs", "PVL", "B-cells", "T-cells", "Myeloid",
                     "Normal Epithelial", "B-cells", "Cancer Epithelial")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$cluster_new<-Idents(pbmc)
DimPlot(pbmc,raster = FALSE,reduction = "umap")
DimPlot(pbmc,raster = FALSE,reduction = "tsne")
pbmc <- RunLDA(pbmc,labels = pbmc$cluster_new,reduction.name = "lda_4")
pbmc <- RunUMAP(pbmc, dims = 1:7,reduction = "lda_4",reduction.name = "lda4_umap")
pbmc <- RunTSNE(pbmc, dims = 1:7,reduction = "lda_4",reduction.name = "lda4_tsne")

library(ggplot2)
library(dplyr)
library(ggalluvial)
library(scales)
Ratio <- pbmc@meta.data %>%
  group_by(subtype, cluster_new) %>% 
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))


mycolor = hue_pal()(8)
Ratio$subtype<-factor(Ratio$subtype,levels = c("TNBC","HER2+","ER+"))
ggplot(Ratio, aes(x =subtype, y= relative_freq, fill = cluster_new,
                  stratum=cluster_new, alluvium=cluster_new)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.8, knot.pos=0.6)+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  # coord_flip()+
  scale_fill_manual(values = mycolor)

Cellratio <- prop.table(table(Idents(pbmc), pbmc$subtype), margin = 2)
Cellratio <- data.frame(Cellratio)

Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  # coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

Cellratio <- prop.table(table(Idents(pbmc), pbmc$Patient), margin = 2)
Cellratio <- data.frame(Cellratio)

Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  # coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))


library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]



sample <- rownames(cellper)
a<-pbmc@meta.data
a<-a[!duplicated(a$Patient),]
rownames(a)<-a$Patient
a<-a[sample,]
group <- a$subtype
remove(a)
samples <- data.frame(sample, group)

rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']
cellper$group <- samples[rownames(cellper),'group']


pplist = list()
sce_groups = unique(levels(Idents(pbmc)))
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group',group_)))
  colnames(cellper_) = c('sample','group','percent')
  cellper_$percent = as.numeric(cellper_$percent)
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + 
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list( c("TNBC", "HER2+"), c("TNBC", "ER+"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,label.y = c(90, 90, 90)) +stat_compare_means(label.y = 100)
  
  pplist[[group_]] = pp1
  
}

DimPlot(pbmc,raster = FALSE,reduction = "umap")
DimPlot(pbmc,raster = FALSE,reduction = "tsne")
DimPlot(pbmc,raster = FALSE,reduction = "lda4_umap")
DimPlot(pbmc,raster = FALSE,reduction = "lda4_tsne")
Bcells<-subset(pbmc,idents = "B-cells")
Idents(Bcells)<-Bcells$subtype
Bcells.markers<-FindAllMarkers(Bcells,only.pos = TRUE,min.pct = 0)
write.table(Bcells.markers,file = "Bcell.markers.txt",sep = "\t")
Idents(Bcells)<-Bcells$Patient
Bsub<-subset(Bcells,idents = c("CID4495","CID4515","CID44041","CID4465","CID4513"))
new.cluster.ids <- c("T2", "T2", "T1", "T3", "T1")
names(new.cluster.ids) <- levels(Bsub)
Bsub <- RenameIdents(Bsub, new.cluster.ids)
Bsub$Tumor<-Idents(Bsub)
Bsub.markers<-FindAllMarkers(Bsub,only.pos = TRUE,min.pct = 0)
write.table(Bsub.markers,file = "Bsub.markers.txt",sep = "\t")
Idents(pbmc)<-pbmc$Patient
pbmcsub<-subset(pbmc,idents = c("CID4495","CID4515","CID44041","CID4465","CID4513"))
new.cluster.ids <- c("T2", "T2", "T1", "T3", "T1")
names(new.cluster.ids) <- levels(pbmcsub)
pbmcsub <- RenameIdents(pbmcsub, new.cluster.ids)
pbmcsub$Tumor<-Idents(pbmcsub)
Ratio <- pbmcsub@meta.data %>%
  group_by(Tumor, cluster_new) %>% 
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))


mycolor = hue_pal()(8)
Ratio$Tumor<-factor(Ratio$Tumor,levels = c("T1","T2","T3"))
ggplot(Ratio, aes(x =Tumor, y= relative_freq, fill = cluster_new,
                  stratum=cluster_new, alluvium=cluster_new)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.8, knot.pos=0.6)+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  # coord_flip()+
  scale_fill_manual(values = mycolor)

Cellratio <- prop.table(table(pbmcsub$cluster_new, pbmcsub$Tumor), margin = 2)
Cellratio <- data.frame(Cellratio)

Cellratio <- as.data.frame(Cellratio)
Cellratio$Var2<-factor(Cellratio$Var2,levels = c("T1","T2","T3"))

colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  # coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

Cellratio <- prop.table(table(pbmcsub$cluster_new, pbmcsub$Patient), margin = 2)
Cellratio <- data.frame(Cellratio)

Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  # coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))


library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]



sample <- rownames(cellper)
a<-pbmcsub@meta.data
a<-a[!duplicated(a$Patient),]
rownames(a)<-a$Patient
a<-a[sample,]
group <- a$Tumor
remove(a)
samples <- data.frame(sample, group)

rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']
cellper$group <- samples[rownames(cellper),'group']


pplist = list()
Idents(pbmcsub)<-pbmcsub$cluster_new
sce_groups = unique(levels(Idents(pbmcsub)))
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group',group_)))
  colnames(cellper_) = c('sample','group','percent')
  cellper_$percent = as.numeric(cellper_$percent)
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + 
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list( c("T1", "T2"), c("T1", "T3"),c("T2","T3"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,label.y = c(90, 80, 70)) +stat_compare_means(label.y = 100)
  
  pplist[[group_]] = pp1
  
}
DimPlot(pbmcsub)
DimPlot(pbmcsub,reduction = "tsne")
DimPlot(pbmcsub,reduction = "lda4_umap")
DimPlot(pbmcsub,reduction = "lda4_tsne")
library(clusterProfiler)
library(org.Hs.eg.db)
ids<-bitr(Bsub.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
Bsub.markers_new<-merge(Bsub.markers,ids,by.x="gene",by.y="SYMBOL")
gcSample<-split(Bsub.markers_new$ENTREZID,Bsub.markers_new$cluster)
S1<-compareCluster(geneClusters = gcSample,fun = "enrichKEGG",organism="hsa",pvalueCutoff=0.05,pAdjustMethod = "BH") #"enrichKEGG" "enrichGO"
S2<-compareCluster(geneClusters = gcSample,fun = "enrichGO",OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05) #"enrichKEGG" "enrichGO"
S3<-compareCluster(geneClusters = gcSample,fun = "enrichMKEGG",organism="hsa",pvalueCutoff=0.05,pAdjustMethod = "BH") #"enrichKEGG" "enrichGO"
S4<-compareCluster(geneClusters = gcSample,fun = "enrichWP",organism = "Homo sapiens",pAdjustMethod = "BH",pvalueCutoff = 0.05) #"enrichKEGG" "enrichGO"
# S5<-compareCluster(geneClusters = gcSample,fun = "enrichDAVID",idType = "ENTREZ_GENE_ID",pAdjustMethod = "BH",pvalueCutoff = 0.05,david.user = "clusterProfiler@hku.hk") #"enrichKEGG" "enrichGO"
write.table(S1,file = "KEGG.txt",sep = "\t")
write.table(S2,file = "GO.txt",sep = "\t")
write.table(S3,file = "MKEGG.txt",sep = "\t")
write.table(S4,file = "WP.txt",sep = "\t")
# write.table(S5,file = "DAVID.txt",sep = "\t")
Bsub.markers$Group<-Bsub.markers$pct.1 - Bsub.markers$pct.2
Bsub.markers_sub<-Bsub.markers[which(Bsub.markers$p_val_adj <0.05 & Bsub.markers$Group > 0.5),]
Bsub.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(Bsub, features = top10$gene) 
DotPlot(Bsub,features = top10$gene)+coord_flip()
# saveRDS(Bcells,file = "d:/Bcells.rds")
# saveRDS(Bsub,file = "d:/Bsub.rds")
# saveRDS(pbmcsub,file = "d:/pbmcsub.rds")
Bsub <- NormalizeData(Bsub)
Bsub <- FindVariableFeatures(Bsub, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(Bsub)
Bsub <- ScaleData(Bsub, features = all.genes)
Bsub <- RunPCA(Bsub, features = VariableFeatures(object = Bsub))
Bsub <- FindNeighbors(Bsub, dims = 1:30)
Bsub <- FindClusters(Bsub, resolution = 1.2)
Bsub <- RunUMAP(Bsub, dims = 1:30)
Bsub <- RunTSNE(Bsub, dims = 1:30)
Bsub.markers<-FindAllMarkers(Bsub,only.pos = TRUE,min.pct = 0)
Bsub.markers$Group<-Bsub.markers$pct.1 - Bsub.markers$pct.2
write.table(Bsub.markers,file = "Bsub.markers.txt",sep = "\t")
ids<-bitr(Bsub.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
Bsub.markers_new<-merge(Bsub.markers,ids,by.x="gene",by.y="SYMBOL")
gcSample<-split(Bsub.markers_new$ENTREZID,Bsub.markers_new$cluster)
S1<-compareCluster(geneClusters = gcSample,fun = "enrichKEGG",organism="hsa",pvalueCutoff=0.05,pAdjustMethod = "BH") #"enrichKEGG" "enrichGO"
S2<-compareCluster(geneClusters = gcSample,fun = "enrichGO",OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05) #"enrichKEGG" "enrichGO"
S3<-compareCluster(geneClusters = gcSample,fun = "enrichMKEGG",organism="hsa",pvalueCutoff=0.05,pAdjustMethod = "BH") #"enrichKEGG" "enrichGO"
S4<-compareCluster(geneClusters = gcSample,fun = "enrichWP",organism = "Homo sapiens",pAdjustMethod = "BH",pvalueCutoff = 0.05) #"enrichKEGG" "enrichGO"
S5<-compareCluster(geneClusters = gcSample,fun = "enrichDAVID",idType = "ENTREZ_GENE_ID",pAdjustMethod = "BH",pvalueCutoff = 0.05,david.user = "clusterProfiler@hku.hk") #"enrichKEGG" "enrichGO"
write.table(S1,file = "KEGG.txt",sep = "\t")
write.table(S2,file = "GO.txt",sep = "\t")
write.table(S3,file = "MKEGG.txt",sep = "\t")
write.table(S4,file = "WP.txt",sep = "\t")
write.table(S5,file = "DAVID.txt",sep = "\t")
DimPlot(Bsub,split.by = "celltype_minor")
library(tidyverse)
Bsub@meta.data <- unite(Bsub@meta.data, "cluster_new2",celltype_minor ,seurat_clusters, sep = "_", remove = FALSE)
Idents(Bsub)<-Bsub$cluster_new2
DimPlot(Bsub)
new.cluster.ids <- c("B cells Memory_1", 
                     "B cells Memory_2", 
                     "B cells Memory_3", 
                     "B cells Memory_Macrophage like B cells",
                     "Plasmablasts_Macrophage like B cells", 
                     "B cells Memory_4",
                     "B cells Memory_5",
                     "B cells Memory__Macrophage like B cells",
                     "B cells Memory_6",
                     "B cells Memory_7",
                     "B cells Memory_Macrophage like B cells",
                     "B cells Naive_1",
                     "B cells Naive_2",
                     "B cells Naive_3",
                     "Plasmablasts_1",
                     "Plasmablasts_2",
                     "Plasmablasts_3",
                     "Plasmablasts_4",
                     "Plasmablasts_5",
                     "Plasmablasts_6",
                     "Plasmablasts_Macrophage like B cells",
                     "Plasmablasts_Macrophage like B cells",
                     "B cells Naive_4",
                     "B cells Naive_5",
                     "B cells Naive_6",
                     "Plasmablasts_7",
                     "Plasmablasts_8"
                    
                      )
names(new.cluster.ids) <- levels(Bsub)
Bsub <- RenameIdents(Bsub, new.cluster.ids)
Bsub$cluster_new3<-Idents(Bsub)
new.cluster.ids <- c("B cells Memory_1", 
                     "B cells Memory_2", 
                     "B cells Memory_3", 
                     "Macrophage like B cells",
                     "Macrophage like B cells", 
                     "B cells Memory_4",
                     "B cells Memory_5",
                     "Macrophage like B cells",
                     "B cells Memory_6",
                     "B cells Memory_7",
                     
                     "B cells Naive_1",
                     "B cells Naive_2",
                     "B cells Naive_3",
                     "Plasmablasts_1",
                     "Plasmablasts_2",
                     "Plasmablasts_3",
                     "Plasmablasts_4",
                     "Plasmablasts_5",
                     "Plasmablasts_6",
                     
                     "B cells Naive_4",
                     "B cells Naive_5",
                     "B cells Naive_6",
                     "Plasmablasts_7",
                     "Plasmablasts_8"
                     
)
names(new.cluster.ids) <- levels(Bsub)
Bsub <- RenameIdents(Bsub, new.cluster.ids)
Bsub$cluster_new4<-Idents(Bsub)
new.cluster.ids <- c("B cells Memory", 
                     "B cells Memory", 
                     "B cells Memory", 
                     "Macrophage like B cells",
                     
                     "B cells Memory",
                     "B cells Memory",
                     
                     "B cells Memory",
                     "B cells Memory",
                     
                     "B cells Naive",
                     "B cells Naive",
                     "B cells Naive",
                     "Plasmablasts",
                     "Plasmablasts",
                     "Plasmablasts",
                     "Plasmablasts",
                     "Plasmablasts",
                     "Plasmablasts",
                     
                     "B cells Naive",
                     "B cells Naive",
                     "B cells Naive",
                     "Plasmablasts",
                     "Plasmablasts"
                     
)
names(new.cluster.ids) <- levels(Bsub)
Bsub <- RenameIdents(Bsub, new.cluster.ids)
Bsub$cluster_new5<-Idents(Bsub)
Ratio <- Bsub@meta.data %>%
  group_by(Tumor, cluster_new5) %>% 
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))

library(scales)
library(ggalluvial)
mycolor = hue_pal()(8)
Ratio$Tumor<-factor(Ratio$Tumor,levels = c("T1","T2","T3"))
ggplot(Ratio, aes(x =Tumor, y= relative_freq, fill = cluster_new5,
                  stratum=cluster_new5, alluvium=cluster_new5)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.8, knot.pos=0.6)+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  # coord_flip()+
  scale_fill_manual(values = mycolor)
Bsub.markers_new_new<-FindAllMarkers(Bsub,only.pos = TRUE,min.pct = 0)
write.table(Bsub.markers_new_new,file = "Bsub.markers_new_new.txt",sep = "\t")
Macro_NB<-FindMarkers(Bsub,ident.1 = "Macrophage like B cells",ident.2 = "B cells Naive",min.pct = 0)
write.table(Macro_NB,file = "Macro_NB.txt",sep = "\t")
DimPlot(Bsub)
DimPlot(Bsub,reduction = "tsne")
Bsub<-RunLDA(Bsub,labels = Bsub$cluster_new5,reduction.name = "lda_1")
Bsub <- RunUMAP(Bsub, dims = 1:3,reduction = "lda_1",reduction.name = "lda1_umap")
Bsub <- RunTSNE(Bsub, dims = 1:3,reduction = "lda_1",reduction.name = "lda1_tsne")
DimPlot(Bsub,reduction = "lda1_umap")
DimPlot(Bsub,reduction = "lda1_tsne")
FeaturePlot(Bsub,features = "CD163",reduction = "lda1_tsne",order = TRUE)
FeaturePlot(Bsub,features = "MS4A6A",reduction = "lda1_tsne",order = TRUE)
FeaturePlot(Bsub,features = "CD163L1",reduction = "lda1_tsne",order = TRUE)
FeaturePlot(Bsub,features = "IL10",reduction = "lda1_tsne",order = TRUE)
FeaturePlot(Bsub,features = "IL12A",reduction = "lda1_tsne",order = TRUE)
FeaturePlot(Bsub,features = "EBI3",reduction = "lda1_tsne",order = TRUE)
FeaturePlot(Bsub,features = "CD14",reduction = "lda1_tsne",order = TRUE)
FeaturePlot(Bsub,features = "CCL18",reduction = "lda1_tsne",order = TRUE)
FeaturePlot(Bsub,features = "TGFB1",reduction = "lda1_tsne",order = TRUE)
FeaturePlot(Bsub,features = "LILRB2",reduction = "lda1_tsne",order = TRUE)
FeaturePlot(Bsub,features = "CCL22",reduction = "lda1_tsne",order = TRUE)
FeaturePlot(Bsub,features = "MRC1",reduction = "lda1_tsne",order = TRUE)
DotPlot(Bsub,features = c("CD14","CCL18","CCL22","LILRB2","CD163","CD163L1","IL10","MS4A6A","MRC1"))

Macrophage_B<-subset(Bsub,idents = "Macrophage like B cells")
NB<-subset(Bsub,idents = "B cells Naive")
write.table(cbind(Macrophage_B@assays$RNA@data,NB@assays$RNA@data),file = "data.txt",sep = "\t")
FeaturePlot(Bsub,features = "MRC1",reduction = "lda1_tsne",order = TRUE)
library(monocle)
trace('project2MST',edit = T,where = asNamespace("monocle"))
data<-as.sparse(Bsub@assays$RNA@counts)
pd <-Bsub@meta.data
pd <- new('AnnotatedDataFrame', data = pd)  
fData<-data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data, phenoData = pd,featureData = fd,lowerDetectionLimit = 0.1,
                              expressionFamily = VGAM::negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
# 
# 
# 
# pbmc.marker<-FindAllMarkers(pbmc_new_new,only.pos = TRUE,min.pct = 0)
diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = "~cluster_new5",cores = 20)

# pbmcmarkers_new<-pbmc.marker[which(pbmc.marker$p_val_adj < 1e-10 ),]
# pbmcmarkers_new$filter<-pbmcmarkers_new$pct.1 - pbmcmarkers_new$pct.2
# ordering_genes <- pbmcmarkers_new$gene
ordering_genes <- row.names (subset(diff_test_res, qval < 1e-2))
monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes = ordering_genes)
monocle_cds <-reduceDimension(monocle_cds,reduction_method = "DDRTree",max_components = 2,norm_method = "log")

monocle_cds <-orderCells(monocle_cds)
plot_cell_trajectory(monocle_cds, color_by = "State",cell_size = 0.75) 
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size = 0.75)
plot_cell_trajectory(monocle_cds, color_by = "cluster_new5",cell_size = 0.75,)
plot_cell_trajectory(monocle_cds, color_by = "cluster_new5",cell_size = 0.75,)+ facet_wrap(~cluster_new5, nrow = 1)

#preparefor slingshot
Bsub_new<-as.SingleCellExperiment(Bsub)

library(slingshot)

sce_6<-slingshot(Bsub_new,clusterLabels = "cluster_new5",reducedDim = 'LDA1_TSNE',start.clus = c("B cells Memory","B cells Naive"),end.clus = c("Macrophage like B cells","Plasmablasts"),dist.method= "slingshot") 
library(TrajectoryUtils)
a<-data.frame(sce_6$slingPseudotime_1,sce_6$slingPseudotime_2)

b<-averagePseudotime(a)
remove(a)
sce_6$slingPseudotime_3<-b
remove(b)
library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_6$slingPseudotime_3, breaks=100)]
plot(reducedDims(sce_6)$LDA1_TSNE, col = plotcol, pch=16, asp = 1) 
lines(SlingshotDataSet(sce_6), lwd=2, col='black')
data<-sce_6@int_colData$reducedDims@listData$LDA1_TSNE
color<-data.frame(sce_6$slingPseudotime_3,plotcol)
rownames(color)<-rownames(data)
data<-cbind(data,color)
colnames(data)[3]<-"Pseudotime"
colnames(data)[4]<-"color"
remove(color)
library(ggplot2)
library(scales)
library(ggsn)
ggplot(data, aes(x = lda1tsne_1, y = lda1tsne_2,fill=Pseudotime)) +
  geom_point(col = plotcol ) + 
  theme_classic() +
  scale_fill_gradientn(colors = colors)
DimPlot(pbmc_new_new,reduction = "lda_tsne")
pbmcsub@meta.data$group<-rep("NO hits",length(rownames(pbmcsub@meta.data)))
pbmc@meta.data$group<-rep("NO hits",length(rownames(pbmc@meta.data)))
a<-colnames(Macrophage_B)
for (i in 1:length(a)) {
  for (j in 1:length(colnames(pbmcsub))) {
    if (a[i] == rownames(pbmcsub@meta.data)[j]){
      pbmcsub@meta.data$group[j] <- "Hits"
    }
  }
  
}
for (i in 1:length(a)) {
  for (j in 1:length(colnames(pbmc))) {
    if (a[i] == rownames(pbmc@meta.data)[j]){
      pbmc@meta.data$group[j] <- "Hits"
    }
  }
  
}
Idents(pbmcsub)<-pbmcsub$group
Macrophage_B_pbmcsub.markers<-FindAllMarkers(pbmcsub,only.pos = TRUE,min.pct = 0)
Idents(pbmc)<-pbmc$group
Macrophage_B_pbmc.markers<-FindAllMarkers(pbmc,only.pos = TRUE,min.pct = 0)
write.table(Macrophage_B_pbmcsub.markers,file = "Macrophage_B_pbmcsub.markers.txt",sep = "\t")
write.table(Macrophage_B_pbmc.markers,file = "Macrophage_B_pbmc.markers.txt",sep = "\t")
