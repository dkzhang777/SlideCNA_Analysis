### Ggplot Parameters
#
text_size = 16
title_size = 18
legend_size_pt = 4
legend_height_bar = 1.5
#

### Get counts from adata object
get_adata_counts <- function(adata) {
    counts <- adata$obsm$counts %>% 
        t() %>% # transpose
        as.data.frame() %>%
        `rownames<-`(adata$uns$counts_var) %>%
        `colnames<-`(adata$obs_names)
    return(counts)
   }

### Make a seurat object
make_seurat_annot <- function(cb, md, ft=TRUE, seed = 1){
    so <- CreateSeuratObject(counts = cb,min.features = 0, min.cells = 3)
    so <- PercentageFeatureSet(so,pattern = "^MT-",col.name = "percent.mito")
    so <- NormalizeData(object = so)
    so <- FindVariableFeatures(object = so)
    so <- ScaleData(object = so,vars.to.regress = c("nCount_RNA","percent.mito"))
    so <- RunPCA(object = so)
    so <- FindNeighbors(object = so)
    so <- FindClusters(object = so,algorithm = 1,random.seed = seed)
    so <- RunTSNE(object = so,dims = 1:10, check_duplicates = FALSE, seed.use = seed)
    so <- RunUMAP(object = so, dims = 1:10, seed.use = seed)
    so <- AddMetaData(so, metadata = md)
        
    return(so)   
}

### Set cluster type labels (malignant or non-malignant status)
set_cluster_type <- function(md, normal_seurat_clusters=None) {
    # use seurat cluster-based annotations of malignant vs non-malignant cell types if not already provided
    if(!('cluster_type' %in% colnames(md))) {
        md[,cluster_type:=ifelse(seurat_clusters%in%normal_seurat_clusters, 'Normal', "Malignant")]
    }
    return(md)
}

### Run enrichR
run_enrichr=function(genes,ngenes){
    res=as.data.table(enrichr(genes[1:ngenes],databases = c("GO_Biological_Process_2018")))
    return(res)
}