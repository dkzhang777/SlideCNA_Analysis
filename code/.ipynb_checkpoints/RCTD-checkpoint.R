### Functions for RCTD implementation on Slide-seq data
source(glue("/ahg/regevdata/users/dzhang/projects/HTAPP_MBC/src/slide_CNV_official_v2/code/utils.R"))

### Make a reference object
makeRef <- function(sn_adata, RCTD_dir) {    
    counts <- get_adata_counts(sn_adata)
    
    cell_types <- sn_adata$obs$cell_type
    names(cell_types) <- rownames(sn_adata$obs) # create cell_types named list
    cell_types <- factor(cell_types, ordered = FALSE) # convert to factor data type
    
    head(counts)
    
    nUMI <- sn_adata$obs$n_counts; names(nUMI) <- rownames(sn_adata$obs) # create nUMI named list

    ### Create the Reference object
    reference <- Reference(counts, cell_types, nUMI)
    #> Warning in Reference(counts, cell_types, nUMI): Reference: nUMI does not match
    #> colSums of counts. If this is unintended, please correct this discrepancy. If
    #> this is intended, there is no problem.
    saveRDS(reference, glue("{RCTD_dir}/RCTD_SCRef.rds"))
    return(reference)
}

### Make a puck object
makePuck <- function(so_adata, RCTD_dir) {
    counts <- so_adata$obsm$counts %>% 
            t() %>% # transpose
            as.data.frame() %>%
            `rownames<-`(so_adata$uns$counts_var) %>%
            `colnames<-`(so_adata$obs_names)
    
    coords <- so_adata$obs %>%
                select(c("x","y")) %>% 
                rename(xcoord=x, ycoord=y)
    
    nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
    
    #> Warning in Reference(counts, cell_types, nUMI): Reference: nUMI does not match
    #> colSums of counts. If this is unintended, please correct this discrepancy. If
    #> this is intended, there is no problem.
    ### Create SpatialRNA object
    puck <- SpatialRNA(coords, counts, nUMI)
    
    ## Examine SpatialRNA object (optional)
    print(dim(puck@counts)) # observe Digital Gene Expression matrix
    hist(log(puck@nUMI,2)) # histogram of log_2 nUMI
    
    saveRDS(puck, glue("{RCTD_dir}/RCTD_puck.rds"))
    return(puck)
}

### Run RCTD
runRCTD <- function(puck, reference, RCTD_dir, max_cores, doublet_mode) {
    myRCTD <- create.RCTD(puck, reference, max_cores=max_cores, CELL_MIN_INSTANCE=1)
    myRCTD <- run.RCTD(myRCTD, doublet_mode=doublet_mode)
    saveRDS(myRCTD, glue("{RCTD_dir}/myRCTD_{doublet_mode}.rds"))
    return(myRCTD)
}