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

### Set cluster type labels (malignant or non-malignant status)
set_cluster_type <- function(md, normal_seurat_clusters = None) {
    # use seurat cluster-based annotations of malignant vs non-malignant cell types if not already provided
    if(!('cluster_type' %in% colnames(md))) {
        md[,cluster_type:=ifelse(seurat_clusters%in%normal_seurat_clusters, 'Normal', "Malignant")]
    }
    return(md)
}