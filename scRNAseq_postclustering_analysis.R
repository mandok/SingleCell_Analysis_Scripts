##########################################
#    This function runs:
#    a) UMAPs for input of seurat object (already clustered and potentially annotated)
#    b) Proportion of celltypes per condition
#    c) DEG (wilcox-test) per cluster
#    d) Geneset enrichment of DEG per cluster using REACTOME, BIOCARTA, KEGG, WIKIPATHWAYS, GO:BP and HALLMARKS
#    e) DEG (wilcox-test) treatment vs control per cluster
#    f) Geneset enrichment of DEG (treatment vs control) per cluster using REACTOME, BIOCARTA, KEGG, WIKIPATHWAYS, GO:BP and HALLMARKS
#    g) Geneset enrichment of Hypoxia/Angiogenesis signatures (from HALLMARKS) from genes in e)
#  
##########################################

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Use the arguments
results_folder_path <- args[1] #path where folder with resuls will be created e.g. "/PATH/test"
clust <- args[2] #seurat object e.g. ( file = "/PATH/Monocyticcells.Rds")
celltype_name <- args[3] # Monocyticcells,  try to not leave spaces
enrichment_tcell_signatures <- args[4] # yes or no,  try to not leave spaces


#0. Load libraries
library(openxlsx)
library(ggpubr)
library(rstatix)
library(dplyr)
library(Seurat)
library(patchwork)
library(scater)
library(biomaRt)
library(scater)
library(robustbase)
library(cowplot)
library(ggplot2)
library(reshape2)
#library(ArchR)
library(gdata)
library(scCustomize)
library(pals)
library(clusterProfiler)
library(org.Mm.eg.db)
library(purrr)
library(msigdbr)
library(viridis)
library(SeuratObject)
library(hypeR)
library(readxl)

colorpalettes<-c("alphabet", "alphabet2", "cols25", "glasbey", "kelly", "okabe", "polychrome","tableau20", "tol", "watlington")

set.seed(23)


#Read seurat object
clust_so<-readRDS(clust)

#Create results folder
results_folder <- paste(celltype_name, length(levels(clust_so$seurat_clusters)), "clusters", sep = "_")
results_folder <- paste(results_folder_path, results_folder, sep = "/")
dir.create(path = results_folder)


#1. Generate UMAPs
#----
cols2<-c("MycP53"="#D51F26","MycP53_Treatment"="#272E6A")
cols4<-sample(colorpalettes, size = 1, replace = F)
cols4<-eval(rlang::parse_expr(paste(cols4, "(", length(levels(clust_so$seurat_clusters)),")", sep="")))
names(cols4)<-names(table(Idents(clust_so)))

p1 <- DimPlot(clust_so, cols = cols2,reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(clust_so, reduction = "umap",cols = cols4[!is.na(names(cols4))],label = TRUE, repel = TRUE)
p1+p2
ggsave(paste(results_folder, "/UMAP_",  celltype_name, "_clusters.pdf", sep=""), width = 12, height = 5)
#----

#2. Calculate DEG per cluster
#----
DefaultAssay(clust_so) <- "RNA"


clust_so.markers <- FindAllMarkers(clust_so, logfc.threshold = 0.1,  min.pct =  0.1, test.use = "wilcox",
                                                 min.diff.pct = 0.1,assay = "RNA",slot="data")

clust_so.markers_top50DEG_Upregulated<-clust_so.markers %>%
  group_by(cluster) %>% filter(p_val_adj<=0.5 & avg_log2FC>0) %>%
  slice_max(n = 100, order_by = avg_log2FC) %>% dplyr::select(cluster, gene) %>% mutate(index = row_number()) %>% tidyr::pivot_wider(names_from = "cluster", values_from = "gene")

#Export as excel
SourceData <- createWorkbook()
addWorksheet(SourceData, paste("Zoomin_", celltype_name, sep = ""))#Add sheets

# Write the data to the sheets
writeData(SourceData, rowNames = T, sheet = paste("Zoomin_", celltype_name, sep = ""), x = clust_so.markers_top50DEG_Upregulated )

# Export the file
saveWorkbook(SourceData, paste(results_folder, "/DEG_", celltype_name,"_clusters.xlsx", sep = ""))
#----
             
#3. Proportions
#----
clust_so$Idents<-Idents(clust_so)
proportions<-clust_so@meta.data %>% group_by(orig.ident,Idents, replicate) %>% summarise(orig.ident,Idents, replicate) %>% table()
#R1
t<-apply(proportions[,,replicate = "R1"], 1, sum)
t2<-sapply(1:2, function(x){proportions[x,,replicate = "R1"]/t[x]})
colnames(t2)<-rownames(proportions)

#R2
t3<-apply(proportions[,,replicate = "R2"], 1, sum)
t4<-sapply(1:2, function(x){proportions[x,,replicate = "R2"]/t3[x]})
colnames(t4)<-rownames(proportions)

#R3
t5<-apply(proportions[,,replicate = "R3"], 1, sum)
t6<-sapply(1:2, function(x){proportions[x,,replicate = "R3"]/t5[x]})
colnames(t6)<-rownames(proportions)

#R4
t7<-apply(proportions[,,replicate = "R4"], 1, sum)
t8<-sapply(1:2, function(x){proportions[x,,replicate = "R4"]/t7[x]})
colnames(t8)<-rownames(proportions)




t<-rbind(data.frame("rep"="R1",reshape2::melt(t2)),
         data.frame("rep"="R2",reshape2::melt(t4)),
         data.frame("rep"="R3",reshape2::melt(t6)),
         data.frame("rep"="R4",reshape2::melt(t8)))
colnames(t)<-c("Rep","cluster","Model","value")


t$cluster<-as.character(t$cluster)

p <- ggbarplot(t[complete.cases(t),], x = "Model", y = "value", add = "mean_se",#facet.by ="Rep",
               fill =  "cluster",palette = cols4[!is.na(names(cols4))])
p
ggsave(paste(results_folder, "/Proportions_", celltype_name, "_clusters.pdf", sep = ""), width = 4, height = 6)
#----

#4. Calculate DEG per cluster per condition (Treatment - Vec)
#----
#Option 1: wilcox test, comparison of individual cells

DefaultAssay(clust_so) <- "RNA"
XLvsVeh_wilcox<-list()
for(cluster in levels(clust_so$seurat_clusters)){
  name<-paste("cluster", cluster, sep="")
  tryCatch({ XLvsVeh_wilcox[[name]]<-FindMarkers(clust_so, logfc.threshold = 0.5,  min.pct =  0.25, ident.1 = "MycP53_Treatment", ident.2 = "MycP53", subset.ident = cluster,group.by = "orig.ident", 
                                                 min.diff.pct = 0.25,assay = "RNA",slot="data")
  },  error=function(e){})
   
}

#Export as excel
SourceData <- createWorkbook()
for(sheet in names(XLvsVeh_wilcox)){
  addWorksheet(SourceData,paste("XLvsVeh_", sheet, sep=""))
  writeData(SourceData, rowNames = T, sheet = paste("XLvsVeh_", sheet, sep=""), x = XLvsVeh_wilcox[[sheet]])
}

#Create subfolder for TreatmentvsVec_Analysis
results_folder_sub<- "/TreatmentvsVec_Analysis"
dir.create(path = paste(results_folder, results_folder_sub, sep = ""))
saveWorkbook(SourceData, paste(results_folder,results_folder_sub,"/DEG_TreatmentvsVec_", celltype_name, "_clusters.xlsx", sep = ""))

#Option2: Pseudo-bulk PENDING
#----

#5. Enrichment for DEGs
#----
#Select multiple databases
m_t2g_Reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g_Biocarta <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:BIOCARTA") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g_Kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g_Wikipathways <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g_gobp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g_hallmarks <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

lm_t2g_dbs<-list(m_t2g_Reactome,m_t2g_Biocarta, m_t2g_Kegg, m_t2g_Wikipathways, m_t2g_gobp, m_t2g_hallmarks)
names(lm_t2g_dbs)<-c("Reactome", "Biocarta", "KEGG", "Wikipathways", "GOBP","Hallmarks")

m_t2g_hallmarks_hypoxia <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol) 
HALLMARKS_Hypoxia<-m_t2g_hallmarks_hypoxia[grep(pattern = "ANGIOGENESIS|VEGF|HYPOXIA", m_t2g_hallmarks_hypoxia$gs_name),] %>% group_by(gs_name) %>% summarise(gene_symbol = list(gene_symbol))
names(HALLMARKS_Hypoxia$gene_symbol)<-HALLMARKS_Hypoxia$gs_name

###---DEG per cluster
results_folder_sub_Enrichment<-"/Enrichment_clusters"
dir.create(path = paste(results_folder, results_folder_sub_Enrichment, sep = ""))

geneid.ls.DEGpercluster <- clust_so.markers_top50DEG_Upregulated[,-1] %>% map(~{
  
  gene.df <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = .x,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
  
  gene <- gene.df$ENTREZID
  gene <- gene[which(!is.na(gene))]
  gene <- unique(gene)
  print(gene)
  return(gene)
})
print("map correct")
#Run enrichment + Plots + Excel export 
SourceData <- createWorkbook()
for(pathway in names(lm_t2g_dbs)){
  tryCatch({ 
    ck <- compareCluster(geneCluster = geneid.ls.DEGpercluster, 
                         fun = enricher,
                         TERM2GENE=lm_t2g_dbs[[pathway]], 
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05,
                         pAdjustMethod = "BH")
    
    ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
    
    #Plot results
    dotplot(ck) + theme_bw()+scale_fill_viridis(option = "D")+
      ggtitle(pathway)+
      theme(
        legend.justification=c(1,0), 
        axis.text.x  = element_text(angle=45, hjust = 1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=8, color="black"))
    ggsave(paste(results_folder,results_folder_sub_Enrichment, "/Enrichment_", pathway, "_DEG_", celltype_name, "_clusters.pdf", sep = ""), width = 6, height = 12)
    
    #Export as excel the deg in the pathways
    addWorksheet(SourceData,pathway)
    writeData(SourceData, rowNames = T, sheet = pathway, x = ck@compareClusterResult)
    },  error=function(e){})

  
}
saveWorkbook(SourceData, paste(results_folder,results_folder_sub_Enrichment,"/Enrichment_results_DEG", celltype_name, "_clusters.xlsx", sep=""))

#Enrichment for Hypoxia signatures from HALLMARKS
hyp_obj <- hypeR(signature = as.list(clust_so.markers_top50DEG_Upregulated[,-1]), genesets = HALLMARKS_Hypoxia$gene_symbol, test = "hypergeometric")
hyp.combined<-hyp_obj$as.list()
hyp.combined2<-do.call(rbind.data.frame, lapply(hyp.combined, function(x){x[c(1:3,6)]}))
#rownames(hyp.combined2)<-gsub(pattern = "\\.HALLMARK_ANGIOGENESIS|\\.HALLMARK_HYPOXIA", replacement = "",rownames(hyp.combined2) )
ggplot(data = hyp.combined2[hyp.combined2$label=="HALLMARK_HYPOXIA",], aes(x=reorder(rownames(hyp.combined2[hyp.combined2$label=="HALLMARK_HYPOXIA",]), -1*log10(fdr)), y=-1*log10(fdr),fill=overlap))+geom_bar(stat="identity", col="black")+coord_flip()+theme_bw()+
  ggtitle("HALLMARK_HYPOXIA")+
  geom_hline(yintercept = 1.3)+
  scale_fill_distiller( palette = "Blues", direction = 1)+
  xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.justification=c(1,0), 
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
ggsave(paste(results_folder,results_folder_sub_Enrichment, "/Enrichment_HALLMARK_HYPOXIA", "_DEG_", celltype_name, "_clusters.pdf", sep = ""), width = 8, height = 5)


ggplot(data = hyp.combined2[hyp.combined2$label=="HALLMARK_ANGIOGENESIS",], aes(x=reorder(rownames(hyp.combined2[hyp.combined2$label=="HALLMARK_ANGIOGENESIS",]), -1*log10(fdr)), y=-1*log10(fdr),fill=overlap))+geom_bar(stat="identity", col="black")+coord_flip()+theme_bw()+
  ggtitle("HALLMARK_ANGIOGENESIS")+
  geom_hline(yintercept = 1.3)+
  scale_fill_distiller( palette = "Greens", direction = 1)+
  xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.justification=c(1,0), 
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
ggsave(paste(results_folder,results_folder_sub_Enrichment, "/Enrichment_HALLMARK_ANGIOGENESIS", "_DEG_", celltype_name, "_clusters.pdf", sep = ""), width = 8, height = 5)

#Export as excel
SourceData <- createWorkbook()
addWorksheet(SourceData,"HALLMARKS_HYPOXIA_ANG")
writeData(SourceData, rowNames = T, sheet = "HALLMARKS_HYPOXIA_ANG", x = do.call(rbind.data.frame, lapply(hyp.combined, function(x){x[c(1:3,6,8)]})))
#Save as excel the deg in the pathways
saveWorkbook(SourceData, paste(results_folder,results_folder_sub_Enrichment,"/Enrichment_results_DEG_TreatmentvsVec_HALLMARKS_HYPOXIA", celltype_name, "_clusters.xlsx", sep=""))


if(enrichment_tcell_signatures=="yes"){
  
  Tcellsignatures_public<-read_excel("/Users/m.ando/surfdrive - Masami Ando Kuri@surfdrive.surf.nl/Documents/HCC/Public/Tcellsignatures.xlsx")
  Tcellsignatures_Zheng_public<-read_excel("/Users/m.ando/surfdrive - Masami Ando Kuri@surfdrive.surf.nl/Documents/HCC/Public/Tcellsignatures_Zheng_Zhang_summarised.xlsx", sheet = 1)
  Tcellsignatures_Zheng_public<-reshape2::melt(as.data.frame(Tcellsignatures_Zheng_public), measure.vars = colnames(Tcellsignatures_Zheng_public))
  
  Tcellsignatures_Zheng_public_mouse<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values =Tcellsignatures_Zheng_public$value, mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 105), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version = 105), uniqueRows=T)
  Tcellsignatures_Zheng_public_mouse$MGI.symbol[match(x = Tcellsignatures_Zheng_public$value, table = Tcellsignatures_Zheng_public_mouse$HGNC.symbol)]
  Tcellsignatures_Zheng_public$mousegene<-Tcellsignatures_Zheng_public_mouse$MGI.symbol[match(x = Tcellsignatures_Zheng_public$value, table = Tcellsignatures_Zheng_public_mouse$HGNC.symbol)]
  
  Tcellsignatures_public_Exhaustion<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$Exhaustion[!is.na(Tcellsignatures_public$Exhaustion)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 105,verbose = TRUE ), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", verbose = TRUE, dataset="mmusculus_gene_ensembl",  version = 105), uniqueRows=T)
  Tcellsignatures_public_G1S<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$`G1/S`[!is.na(Tcellsignatures_public$`G1/S`)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 105), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version = 105), uniqueRows=T)
  Tcellsignatures_public_G2M<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$`G2/M`[!is.na(Tcellsignatures_public$`G2/M`)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 105), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version = 105), uniqueRows=T)
  Tcellsignatures_public_Exhaustion_CD8<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$Exhaustion_CD8[!is.na(Tcellsignatures_public$Exhaustion_CD8)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 105), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version = 105), uniqueRows=T)
  Tcellsignatures_public_Exhaustion_CD4<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$Exhaustion_CD4[!is.na(Tcellsignatures_public$Exhaustion_CD4)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 105), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version = 105), uniqueRows=T)
  
  
  Tcellsignatures_public_list<-list("Memory"=Tcellsignatures_public$Memory,"Effector"=Tcellsignatures_public$Effector,
                                    "Exhaustion"=Tcellsignatures_public_Exhaustion$MGI.symbol,
                                    "G1S"=Tcellsignatures_public_G1S$MGI.symbol,  "G2M"=Tcellsignatures_public_G2M$MGI.symbol,
                                    "Exhaustion_CD8"=Tcellsignatures_public_Exhaustion_CD8$MGI.symbol,"Exhaustion_CD4"=Tcellsignatures_public_Exhaustion_CD4$MGI.symbol,
                                    "Progenitor_Exh"=Tcellsignatures_public$Progenitor_Exh,"Effector_like"=Tcellsignatures_public$Effector_like,
                                    "Terminally_Exh"=Tcellsignatures_public$Terminally_Exh, "Proliferating"=Tcellsignatures_public$Proliferating)
  
  Tcellsignatures_Zheng_list<-list("C1_CD8−LEF1"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C1_CD8−LEF1"],
                                   "C2_CD8-CX3CR1"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C2_CD8-CX3CR1"],
                                   "C3_CD8-SLC4A10"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C3_CD8-SLC4A10"],
                                   "C4_CD8-LAYN"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C4_CD8-LAYN"],
                                   "C5_CD8-GZMK"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C5_CD8-GZMK"],
                                   "C6_CD4-CCR7"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C6_CD4-CCR7"],
                                   "C7_CD4-FOXP3"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C7_CD4-FOXP3"],
                                   "C8_CD4-CTLA4"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C8_CD4-CTLA4"],
                                   "C9_CD4-GZMK"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C9_CD4-GZMK"],
                                   "C10_CD4-CXCL13"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C10_CD4-CXCL13"],
                                   "C11_CD4-GNLY"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C11_CD4-GNLY"])
  
  
  hyp_obj <- hypeR(signature = as.list(clust_so.markers_top50DEG_Upregulated[,-1]), genesets = Tcellsignatures_Zheng_list, test = "hypergeometric")
  hyp_dots(hyp_obj,merge=TRUE,fdr = 0.05)+
    theme(text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
  
  hyp.combined<-hyp_obj$as.list()
  hyp.combined2<-do.call(rbind.data.frame, lapply(hyp.combined, function(x){x[c(1:3,5:6)]}))
  hyp.combined2$cluster<-gsub(pattern = "\\.C[0-9]+_CD[0-9]+-[A-Z0-9]+|.C1_CD8−LEF1",replacement = "", rownames(hyp.combined2))
  hyp.combined2$group<-unlist(lapply(strsplit(hyp.combined2$label, split = "_"), function(x){gsub(pattern = "-[A-Z0-9]+|−LEF1", replacement = "", x[2])}))
  ggplot(data = hyp.combined2[hyp.combined2$fdr<=0.05, ], aes(x=label, y= cluster,col=fdr))+geom_point(aes(size = geneset))+
    facet_wrap(~group, scales="free_x")+theme_bw()+
    theme(#panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      legend.justification=c(1,0), 
      #legend.title = element_text("Clusters"),  
      axis.text.x  = element_text(angle=45, hjust = 1),
      legend.background = element_blank(),
      panel.background = element_blank(),
      legend.key = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))+
    scale_color_viridis(option = "B", direction = 1)
  ggsave(paste(results_folder,results_folder_sub_Enrichment, "/Enrichment_Zheng_Tcellsignatures_DEG_", celltype_name, "_clusters.pdf", sep = ""), width = 10, height = 5)
  
  
  hyp_obj <- hypeR(signature = as.list(clust_so.markers_top50DEG_Upregulated[,-1]), genesets = Tcellsignatures_public_list, test = "hypergeometric")
  hyp_dots(hyp_obj,merge=TRUE,fdr = 0.05)+
    theme(text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
  
  hyp.combined<-hyp_obj$as.list()
  hyp.combined2<-do.call(rbind.data.frame, lapply(hyp.combined, function(x){x[c(1:3,5:6)]}))
  hyp.combined2$cluster<-gsub(pattern = "\\.[A-Za-z0-9_]+",replacement = "", rownames(hyp.combined2))
  ggplot(data = hyp.combined2[hyp.combined2$fdr<=0.05, ], aes(x=label, y= cluster,col=fdr))+geom_point(aes(size = geneset))+
    theme_bw()+
    theme(#panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      legend.justification=c(1,0), 
      #legend.title = element_text("Clusters"),  
      axis.text.x  = element_text(angle=45, hjust = 1),
      legend.background = element_blank(),
      panel.background = element_blank(),
      legend.key = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))+
    scale_color_viridis(option = "B", direction = 1)
  ggsave(paste(results_folder,results_folder_sub_Enrichment, "/Enrichment_Public_Tcellsignatures_DEG_", celltype_name, "_clusters.pdf", sep = ""), width = 8, height = 5)
  
}

###---DEG per cluster per condition (Treatment - Vec)
results_folder_sub_sub<- "/Enrichment"
dir.create(path = paste(results_folder, results_folder_sub,results_folder_sub_sub, sep = ""))

#Select only the Gene names
UP_TreatmentvsVec<-lapply(XLvsVeh_wilcox, function(x){
  x %>% filter(p_val_adj<=0.05) ->filt
  filt$gene<-rownames(filt)
  up<-filt$gene[filt$avg_log2FC>0 ]
})

DOWN_TreatmentvsVec<-lapply(XLvsVeh_wilcox, function(x){
  x %>% filter(p_val_adj<=0.05) ->filt
  filt$gene<-rownames(filt)
  up<-filt$gene[filt$avg_log2FC<0 ]
})


#Convert genename to entrez id
geneid.ls.UP_TreatmentvsVec <- UP_TreatmentvsVec %>% map(~{
  
  
  gene.df <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = .x,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
  
  gene <- gene.df$ENTREZID
  gene <- gene[which(!is.na(gene))]
  gene <- unique(gene)
  
  return(gene)
})

geneid.ls.DOWN_TreatmentvsVec <- DOWN_TreatmentvsVec %>% map(~{
  gene.df <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = .x,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
  
  gene <- gene.df$ENTREZID
  gene <- gene[which(!is.na(gene))]
  gene <- unique(gene)
  
  return(gene)
})

#Run enrichment + Plots + Excel export 
SourceData <- createWorkbook()
for(pathway in names(lm_t2g_dbs)){
  
  #Upregulated genes
  ck_up <- compareCluster(geneCluster = geneid.ls.UP_TreatmentvsVec, 
                       fun = enricher,
                       TERM2GENE=lm_t2g_dbs[[pathway]], 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
  
  ck_up <- setReadable(ck_up, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  
  #Plot results
  dotplot(ck_up) + theme_bw()+scale_fill_viridis(option = "D")+
    ggtitle(pathway)+
    theme(
      legend.justification=c(1,0), 
      axis.text.x  = element_text(angle=45, hjust = 1),
      legend.background = element_blank(),
      legend.key = element_blank(),
      panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=8, color="black"))
  ggsave(paste(results_folder,results_folder_sub,results_folder_sub_sub, "/Enrichment_", pathway, "_Upregulated_TreatmentvsVec_",celltype_name, "_clusters.pdf", sep = ""), width = 6, height = 12)
  
  #Save as excel the deg in the pathways
  addWorksheet(SourceData,paste("Up_TreatmentvsVec_", pathway, sep = ""))
  writeData(SourceData, rowNames = T, sheet = paste("Up_TreatmentvsVec_", pathway, sep = ""), x = ck_up@compareClusterResult)
  
  #Downregulated genes
  ck_down <- compareCluster(geneCluster = geneid.ls.DOWN_TreatmentvsVec, 
                          fun = enricher,
                          TERM2GENE=lm_t2g_dbs[[pathway]], 
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05,
                          pAdjustMethod = "BH")
  
  ck_down <- setReadable(ck_down, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  
  #Plot results
  dotplot(ck_down) + theme_bw()+scale_fill_viridis(option = "D")+
    ggtitle(pathway)+
    theme(
      legend.justification=c(1,0), 
      axis.text.x  = element_text(angle=45, hjust = 1),
      legend.background = element_blank(),
      legend.key = element_blank(),
      panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=8, color="black"))
  ggsave(paste(results_folder,results_folder_sub, results_folder_sub_sub,"/Enrichment_", pathway, "_Downregulated_TreatmentvsVec_",celltype_name, "clusters.pdf", sep = ""), width = 6)
  
  #Export as excel the deg in the pathways
  addWorksheet(SourceData,paste("Down_TreatmentvsVec_", pathway, sep = ""))
  writeData(SourceData, rowNames = T, sheet = paste("Down_TreatmentvsVec_", pathway, sep = ""), x = ck_down@compareClusterResult)
  
  
  
  saveWorkbook(SourceData, paste(results_folder,results_folder_sub,results_folder_sub_sub,"/Enrichment_results_DEG_TreatmentvsVec_", pathway, "_",celltype_name,"_clusters.xlsx", sep=""))
  
}

#Enrichment for Hypoxia signatures from HALLMARKS
SourceData <- createWorkbook()
hyp_obj <- hypeR(signature = UP_TreatmentvsVec[!unlist(lapply(UP_TreatmentvsVec, isEmpty))], genesets = HALLMARKS_Hypoxia$gene_symbol, test = "hypergeometric")
hyp.combined<-hyp_obj$as.list()
hyp.combined2<-do.call(rbind.data.frame, lapply(hyp.combined, function(x){x[c(1:3,6)]}))
#rownames(hyp.combined2)<-gsub(pattern = "\\.HALLMARK_ANGIOGENESIS|\\.HALLMARK_HYPOXIA", replacement = "",rownames(hyp.combined2) )
ggplot(data = hyp.combined2[hyp.combined2$label=="HALLMARK_HYPOXIA",], aes(x=reorder(rownames(hyp.combined2[hyp.combined2$label=="HALLMARK_HYPOXIA",]), -1*log10(fdr)), y=-1*log10(fdr),fill=overlap))+geom_bar(stat="identity", col="black")+coord_flip()+theme_bw()+
  ggtitle("HALLMARK_HYPOXIA")+
  geom_hline(yintercept = 1.3)+
  scale_fill_distiller( palette = "Blues", direction = 1)+
  xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.justification=c(1,0), 
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
ggsave(paste(results_folder,results_folder_sub, results_folder_sub_sub,"/Enrichment_HALLMARK_HYPOXIA", "_UPregulated_TreatmentvsVec_",celltype_name, "clusters.pdf", sep = ""),width = 8, height = 5)

ggplot(data = hyp.combined2[hyp.combined2$label=="HALLMARK_ANGIOGENESIS",], aes(x=reorder(rownames(hyp.combined2[hyp.combined2$label=="HALLMARK_ANGIOGENESIS",]), -1*log10(fdr)), y=-1*log10(fdr),fill=overlap))+geom_bar(stat="identity", col="black")+coord_flip()+theme_bw()+
  ggtitle("HALLMARK_ANGIOGENESIS")+
  geom_hline(yintercept = 1.3)+
  scale_fill_distiller( palette = "Greens", direction = 1)+
  xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.justification=c(1,0), 
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
ggsave(paste(results_folder,results_folder_sub, results_folder_sub_sub,"/Enrichment_HALLMARK_ANGIOGENESIS", "_UPregulated_TreatmentvsVec_",celltype_name, "clusters.pdf", sep = ""),width = 8, height = 5)

#Save as excel the deg in the pathways
addWorksheet(SourceData,"Up_TreatmentvsVec_HYPOXIA")
writeData(SourceData, rowNames = T, sheet = "Up_TreatmentvsVec_HYPOXIA", x = do.call(rbind.data.frame, lapply(hyp.combined, function(x){x[c(1:3,6,8)]})))



#Enrichment for Hypoxia signatures from HALLMARKS
hyp_obj <- hypeR(signature = DOWN_TreatmentvsVec[!unlist(lapply(DOWN_TreatmentvsVec, isEmpty))], genesets = HALLMARKS_Hypoxia$gene_symbol, test = "hypergeometric")
hyp.combined<-hyp_obj$as.list()
hyp.combined2<-do.call(rbind.data.frame, lapply(hyp.combined, function(x){x[c(1:3,6)]}))
#rownames(hyp.combined2)<-gsub(pattern = "\\.HALLMARK_ANGIOGENESIS|\\.HALLMARK_HYPOXIA", replacement = "",rownames(hyp.combined2) )
ggplot(data = hyp.combined2[hyp.combined2$label=="HALLMARK_HYPOXIA",], aes(x=reorder(rownames(hyp.combined2[hyp.combined2$label=="HALLMARK_HYPOXIA",]), -1*log10(fdr)), y=-1*log10(fdr),fill=overlap))+geom_bar(stat="identity", col="black")+coord_flip()+theme_bw()+
  ggtitle("HALLMARK_HYPOXIA")+
  geom_hline(yintercept = 1.3)+
  scale_fill_distiller( palette = "Blues", direction = 1)+
  xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.justification=c(1,0), 
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
ggsave(paste(results_folder,results_folder_sub, results_folder_sub_sub,"/Enrichment_HALLMARK_HYPOXIA", "_Downregulated_TreatmentvsVec_",celltype_name, "clusters.pdf", sep = ""),width = 8, height = 5)


ggplot(data = hyp.combined2[hyp.combined2$label=="HALLMARK_ANGIOGENESIS",], aes(x=reorder(rownames(hyp.combined2[hyp.combined2$label=="HALLMARK_ANGIOGENESIS",]), -1*log10(fdr)), y=-1*log10(fdr),fill=overlap))+geom_bar(stat="identity", col="black")+coord_flip()+theme_bw()+
  ggtitle("HALLMARK_ANGIOGENESIS")+
  geom_hline(yintercept = 1.3)+
  scale_fill_distiller( palette = "Greens", direction = 1)+
  xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.justification=c(1,0), 
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
ggsave(paste(results_folder,results_folder_sub, results_folder_sub_sub,"/Enrichment_HALLMARK_ANGIOGENESIS", "_Downregulated_TreatmentvsVec_",celltype_name, "clusters.pdf", sep = ""),width = 8, height = 5)

#Export results in excel file

addWorksheet(SourceData,"Down_TreatmentvsVec_HYPOXIA")
writeData(SourceData, rowNames = T, sheet = "Down_TreatmentvsVec_HYPOXIA", x = do.call(rbind.data.frame, lapply(hyp.combined, function(x){x[c(1:3,6,8)]})))
#Save as excel the deg in the pathways

saveWorkbook(SourceData, paste(results_folder,results_folder_sub,results_folder_sub_sub,"/Enrichment_results_DEG_TreatmentvsVec_HALLMARKS_HYPOXIA" ,celltype_name,"_clusters.xlsx", sep=""))

#----



