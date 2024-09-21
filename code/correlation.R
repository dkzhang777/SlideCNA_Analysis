### Compute the Spearman correlation between the average CNA score per chromosome arm for each pair of CNA analyses

### Libraries
library(tidyr)
library(stringr)
library(biomaRt)
library(ggpubr)
library(corrplot)
library(GenomicRanges)
library(rtracklayer)

### Map shorthand of analysis to its full label
analysis_to_label <- list("infercnv_slide"= "Infercnv and Slide-seq",
                            "infercnv_sc"= "Infercnv and Single Nucleus",
                            "slidecna_slide"= "SlideCNA and Slide-seq",
                            "slidecna_sc"= "SlideCNA and Single Nucleus",
                              "copykat_slide"= "CopyKAT and Slide-seq",
                             "copykat_sc"="CopyKAT and Single Nucleus",
                          "absolute_wes"="ABSOLUTE and WES",
                           "in_silico_sep1"="In Silico Non-Malignant Separation Cluster 1",
                            "in_silico_sep2"="In Silico Non-Malignant Separation Cluster 2",
                           "in_silico_sep3"="In Silico Non-Malignant Separation Cluster 3",
                          "in_silico_mix1"="In Silico Non-Malignant Mixing Cluster 1",
                            "in_silico_mix2"="In Silico Non-Malignant Mixing Cluster 2",
                           "in_silico_mix3"="In Silico Non-Malignant Mixing Cluster 3",
                          "HTAPP_895_infercnv_sc"="HTAPP-895 Single Nucleus and Infercnv",
                             "HTAPP_944_infercnv_sc"="HTAPP-944 Single Nucleus and Infercnv",
                         "slidecna_rep1"= "SlideCNA Replicate 1",
                         "slidecna_rep2"= "SlideCNA Replicate 2",
                         "slidecna_joint"= "SlideCNA Combined")

# abbreviated analysis_to_label dictionary
analysis_to_label_abr <- list("infercnv_slide"= "Infercnv and Slide-seq",
                            "infercnv_sc"= "Infercnv and Single Nucleus",
                            "slidecna_slide"= "SlideCNA and Slide-seq",
                            "slidecna_sc"= "SlideCNA and Single Nucleus",
                              "copykat_slide"= "CopyKAT and Slide-seq",
                             "copykat_sc"="CopyKAT and Single Nucleus",
                              "absolute_wes"="ABSOLUTE and WES",
                           "in_silico_sep1"="In Silico Cluster 1",
                            "in_silico_sep2"="In Silico Cluster 2",
                           "in_silico_sep3"="In Silico Cluster 3",
                          "in_silico_mix1"="In Silico Cluster 1",
                            "in_silico_mix2"="In Silico Cluster 2",
                           "in_silico_mix3"="In Silico Cluster 3",
                          "HTAPP_895_infercnv_sc"="HTAPP-895 SN and Infercnv",
                             "HTAPP_944_infercnv_sc"="HTAPP-944 SN and Infercnv",
                             "slidecna_rep1"= "SlideCNA Replicate 1",
                         "slidecna_rep2"= "SlideCNA Replicate 2",
                         "slidecna_joint"= "SlideCNA Combined")

### Get dataframe of average CNA score per chromosome arm given Infercnv Object
infercnv_bands <- function(infercnv_obj, gene_pos_complete, analysis, malig_filt=TRUE) {
    
    # Annotate average CNA score per cell/bead per gene with chromosome info and Ensembl ID
    # filtering for cell indices with status 'Malignant'
    if (malig_filt) { 
        malig_i <- infercnv_obj@observation_grouped_cell_indices$Malignant
        cna_score <- merge(as.data.frame(gene_pos_complete),
                    as.data.frame(infercnv_obj@expr.data[,c(malig_i)]),
                    by.x="GENE", by.y="row.names")
        }
    else {
        cna_score <- merge(as.data.frame(gene_pos_complete),
                    as.data.frame(infercnv_obj@expr.data),
                    by.x="GENE", by.y="row.names")
    }
    print(paste0("N genes: ", nrow(cna_score)))
    
    # Get avg cna score per chromosome arm using Ensembl biomart
    avg_cna_score <- get_band_avg(cna_score, analysis)

    return(avg_cna_score)
   }

### Get dataframe of average CNA score per chromosome arm given SlideCNA cna scores and dat_bin
slidecna_bands <- function(sub_wide, dat_bin, gene_pos_complete, analysis) {
    
    # Reformat CNA score dataframe to get gene labels
    sub_wide <- sub_wide %>% t()
    genes <- dat_bin %>% 
        as.data.frame() %>%
        dplyr::select(GENE, chr, start, end) %>%
        distinct()
    
    slide_expr_data <- cbind(genes, sub_wide)
    
    # Annotate average CNA score per cell/bead per gene with chromosome info and Ensembl ID
    cna_score <- merge(as.data.frame(gene_pos_complete),
                    slide_expr_data,
                    by=c("GENE", "chr", "start", "end"))
    print(paste0("N genes: ", nrow(cna_score)))

    # Get avg cna score per chromosome arm using Ensembl biomart
    avg_cna_score <- get_band_avg(cna_score, analysis)

    return(avg_cna_score)
   }

### Get dataframe of average CNA score per chromosome arm given CopyKAT cna scores
copykat_bands <- function(ck_results, 
                          md, 
                          sample, 
                          analysis, 
                          malig_filt=TRUE) {
    
    # filter to only include malignant beads
    if(malig_filt) {
        # get row names of CopyKAT results that are malignant beads

        md_malig <- md %>% dplyr::select(c(bc, cluster_type)) %>% 
        mutate(bc_alt=str_replace_all(bc, "-", ".")) %>%
        filter(cluster_type=="Malignant")
        
        # if(analysis=="copykat_slide") {
        #     md_malig <- md %>% dplyr::select(c(bc, cluster_type)) %>% 
        #     mutate(bc_alt=str_replace_all(bc, "-", ".")) %>%
        #     filter(cluster_type=="Malignant")
        # }
        # else if(analysis == "copykat_sc"){
        #     md_malig <- md %>% dplyr::select(c(bc, cluster_type)) %>% 
        #     # mutate(bc_alt=paste0(sample,".TST.channel1_",bc)) %>%
        #     mutate(bc_alt=str_replace_all(bc, "-", ".")) %>%
        #     filter(cluster_type=="Malignant")
        # }

        # # get rid of sample.TST.channel... prefix in bead names           
        # colnames(ck_results)[8:ncol(ck_results)] <- sub('.*_', '', colnames(ck_results)[8:ncol(ck_results)])
        
         # select only malignant beads in CopyKAT CNA results
        malig_i <- which(colnames(ck_results) %in% md_malig$bc_alt)
        
        ck_results <- ck_results[,c(1:7, malig_i)]
        
    }

    # get chromosome arms
    ck_results <- ck_results %>% 
        mutate(arm = paste0(chromosome_name,
                            substr(band, 1, 1))) %>%
        relocate(arm, .before =1)

    # collect average CNA score per gene across all cells/beads
    # add 1 to CNA score to put on same scale as other analyses
    ck_results$cna_avg <- rowMeans(ck_results[,c(9:ncol(ck_results))]) + 1

    # collect average CNA score per chrom arm across all beads and genes within that arm
    avg_cna_score <- aggregate(cna_avg~arm, ck_results, FUN=mean)
    avg_cna_score <- avg_cna_score[(str_order(avg_cna_score$arm, numeric = TRUE)),]
        
    # rename the avg cna score to be the type of analysis (method + data type)
    colnames(avg_cna_score) <- c("arm", analysis)
    
    # rename chromosome 23 as X
    if ("23p" %in% avg_cna_score$arm) {
        avg_cna_score[avg_cna_score$arm == "23p",]$arm <- "Xp"
    }
    if ("23q" %in% avg_cna_score$arm) {
        avg_cna_score[avg_cna_score$arm == "23q",]$arm <- "Xq"
    }
    
    return(avg_cna_score)
}

### Get dataframe of average CNA score per chromosome arm given ABSOLUTE cna values per segment
absolute_bands <- function(wes_data, 
                          analysis) {
    
    # get centromere locations
    
    # Connect to UCSC and set the genome to hg19
    ucsc_session <- browserSession("UCSC")
    genome(ucsc_session) <- "hg19"

    # get cytoband talbe
    query <- ucscTableQuery(ucsc_session, table = "gap") 
    cent_tbl <- getTable(query)

    # filter for centromeres
    centromere <- subset(cent_tbl, grepl("centromere", type)) %>%
        mutate_at('chrom', ~ sub(paste0(".*", "chr"), "", .)) 
    
    # initialize WES data frame with arm column
    wes_arm <- data.frame(Sample=character(),
                      Chromosome=character(),
                      Start=integer(),
                      End=integer(),
                      Num_Probes=integer(),
                      Segment_Mean=double(),
                      arm=character()) 
    
    # loop over each segment
    for(i in 1:nrow(wes_data)) {
    
        Sample <- wes_data$Sample[i]
        Chromosome <- wes_data$Chromosome[i]
        Start <- wes_data$Start[i]
        End <- wes_data$End[i]
        Num_Probes <- wes_data$Num_Probes[i]
        Segment_Mean <- wes_data$Segment_Mean[i]

        # get centromere start and end positions
        cent_start <- centromere[centromere$chrom==as.character(Chromosome),]$chromStart
        cent_end <- centromere[centromere$chrom==as.character(Chromosome),]$chromEnd

        # segment starts on p arm
        if (Start < cent_end) { 
            if (End < cent_end) {
                # segment ends on p arm
                new_arm <- glue("{Chromosome}p")
                # Return the original row with the arm label
                new_row <- data.frame(Sample = Sample, 
                                   Chromosome = Chromosome, 
                                   Start = Start, 
                                   End = End,
                                   Num_Probes = Num_Probes,
                                   Segment_Mean = Segment_Mean,
                                   arm = new_arm)

                wes_arm <- wes_arm %>% rbind(new_row)
            }
            # segment starts on p and ends on q
            else {
                # separate segment into 2 arms
                p_row <- data.frame(Sample = Sample, 
                                    Chromosome = Chromosome, 
                                    Start = Start, 
                                    End = cent_end,
                                    Num_Probes = Num_Probes / 2,
                                    Segment_Mean = Segment_Mean,
                                    arm = glue("{Chromosome}p"))
                q_row <- data.frame(Sample = Sample, 
                                   Chromosome = Chromosome, 
                                   Start = cent_end, 
                                   End = End,
                                   Num_Probes = Num_Probes / 2,
                                   Segment_Mean = Segment_Mean,
                                   arm = glue("{Chromosome}q"))

                wes_arm <- wes_arm %>% rbind(p_row) %>% rbind(q_row)
            }

        }
        # segment starts on q arm
        else {
            new_arm <- glue("{Chromosome}q")
            # Return the original row with the arm label
            new_row <- data.frame(Sample = Sample, 
                               Chromosome = Chromosome, 
                               Start = Start, 
                               End = End,
                               Num_Probes = Num_Probes,
                               Segment_Mean = Segment_Mean,
                               arm = new_arm)

            wes_arm <- wes_arm %>% rbind(new_row)
        }
    }
    
    # Normalize by length of segment
    wes_arm <- wes_arm %>% 
        mutate(arm = as.character(arm)) %>% 
        mutate(Length = End - Start)
    
    arm_lengths <- aggregate(Length~arm, wes_arm, FUN=sum) %>% 
        rename(Length = "Arm_Length")
    
    wes_arm <- wes_arm %>% 
        merge(arm_lengths, by = "arm")
    
        # normalize and center at 1
    wes_arm <- wes_arm %>% 
mutate(Segment_Mean_Norm = (Segment_Mean * Length / Arm_Length)+1)
    
    avg_cna_score <- aggregate(Segment_Mean_Norm~arm, wes_arm, FUN=mean) %>%
rename(Segment_Mean_Norm = analysis)
    
    avg_cna_score <- avg_cna_score[(str_order(avg_cna_score$arm, numeric = TRUE)),]
    
    return(avg_cna_score)
}

### Given a dataframe (genes x (cells/beads/bins + gene metadata)) of cna scores
### calculate the avg cna score per chromosome arm using Ensembl biomart to get karyotype bands
get_band_avg <- function(cna_score, analysis) {
    
    # get ensembl biomart for humans
    ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

    # Only use normal chromosomes in humans
    normal_chroms <- c(1:22, "X", "Y", "MT")
    
    # for every ensembl ID, annotate with karyotype band and chromosome arm
    my_genes <- getBM(c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", 
                          "start_position", "end_position", "band"),
                        filters = c("ensembl_gene_id", "chromosome_name"),
                        values = list(ensembl_gene_id=cna_score$ENSEMBL_ID,
                                     chromosome_name=normal_chroms),
                                      mart = ensembl) %>%
                 mutate(arm = paste0(chromosome_name, substr(band, 1, 1))) 
    
    print(paste0("N genes with karyotype band annotations: ", length(unique(my_genes$hgnc_symbol))))
          
    # annotate CNA scores with karyotype band and chromosome arm
    cna_band <- merge(my_genes, cna_score, by.x="ensembl_gene_id", by.y="ENSEMBL_ID")
    
    # collect average CNA score per gene across all cells/beads
    cna_band$cna_avg <- rowMeans(cna_band[,c(12:ncol(cna_band))])
    
    # collect average CNA score per chrom arm across all beads and genes within that arm
    avg_cna_score <- aggregate(cna_avg~arm, cna_band, FUN=mean)
    avg_cna_score <- avg_cna_score[(str_order(avg_cna_score$arm, numeric = TRUE)),]
    
    # rename the avg cna score to be the type of analysis (method + data type)
    colnames(avg_cna_score) <- c("arm", analysis)
    
    return(avg_cna_score)
}

### Make heat map of average CNAs across chromosome arms for an analysis or analyses
avg_cna_heatmap <- function(all_cna_avgs, modality, out_dir) {
    
    # set which analyses we will plot
    if (modality=='infercnv') {
        analyses_to_plot <- c('Infercnv and Slide-seq', 'Infercnv and Single Nucleus')
    }
    else if (modality == 'slidecna') {
        analyses_to_plot <- c('SlideCNA and Slide-seq', 'SlideCNA and Single Nucleus')
    }
    else if (modality == 'copykat') {
        analyses_to_plot <- c('CopyKAT and Slide-seq', 'CopyKAT and Single Nucleus')
    }
    else if (modality == 'absolute') {
        analyses_to_plot <- c('ABSOLUTE and WES')
    }
    else if (modality == 'HTAPP_infercnv') {
        analyses_to_plot <- c('HTAPP-895 Single Nucleus and Infercnv',
                              'HTAPP-944 Single Nucleus and Infercnv')
    }
    else if (modality == 'HTAPP_in-silico_sep') {
        analyses_to_plot <- c('In Silico Non-Malignant Separation Cluster 1', 
                              'In Silico Non-Malignant Separation Cluster 2', 
                              'In Silico Non-Malignant Separation Cluster 3')
    }
    else if (modality == 'HTAPP_in-silico_mix') {
        analyses_to_plot <- c('In Silico Non-Malignant Mixing Cluster 1', 
                              'In Silico Non-Malignant Mixing Cluster 2', 
                              'In Silico Non-Malignant Mixing Cluster 3')
    }
    else if (modality == 'replicate') {
        analyses_to_plot <- c('SlideCNA Replicate 1',
                              'SlideCNA Replicate 2',
                              'SlideCNA Combined')
    }

    # convert data to long format
    cna_avgs_long <- reshape2::melt(all_cna_avgs,
                                    id.vars = c("arm")) %>% 
        `colnames<-`(c("arm", "analysis", "avg_cna_score")) %>% 
        mutate(analysis=as.character(analysis)) %>%
        mutate(analysis=as.character(analysis_to_label[analysis])) 
    
    # set order of chromosome arms and analyses
    cna_avgs_long$arm <- factor(cna_avgs_long$arm, levels=unique(cna_avgs_long$arm)) 
    cna_avgs_long$analysis <- factor(cna_avgs_long$analysis, 
                                     levels=rev(analyses_to_plot))

    # Make heat map
    options(repr.plot.width = 25, repr.plot.height = 3)

    gg <- ggplot(cna_avgs_long[cna_avgs_long$analysis %in% analyses_to_plot,], 
                 aes(x = arm, 
                     y = analysis, 
                     fill = avg_cna_score)) +
    geom_tile(color = "black",
              lwd = 0.5,
              linetype = 1)  +
    scale_fill_gradient2(midpoint=1, 
                         low="navy",
                         mid="white", 
                         high="firebrick3") +
    theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(),
             axis.text=element_text(size=text_size),
            axis.title=element_text(size=text_size),
             legend.title=element_text(size=title_size),
            legend.text=element_text(size=text_size)) +
        xlab("Chromosome Arm") +
        ylab("Analysis") +
        labs(fill = "Average CNA Score")
    
    # shorten height of plot if only one analysis
    if (length(analyses_to_plot) == 1) {
        h = 2
    }
    else {
        h = 3
    }
    
    pdf(file = glue("{out_dir}/avg_cna_heatmap_{modality}.pdf"), width = 25, height = h)
    print(gg)
    dev.off()
    print(gg)
}

### Calculate and plot heat map of  pairwise spearman correlation of average CNA score per arm for all analyses
all_corr_plot <- function(all_cna_avgs, out_dir, modality="HTAPP") {
    
    # Get matrix of spearman correlation for all analysis-analysis pairs
    if (modality == "HTAPP") {
        all_corr <- all_cna_avgs %>% 
                        dplyr::select(-arm) %>%
                        dplyr::select(slidecna_slide, slidecna_sc, 
                                      infercnv_slide, infercnv_sc,
                                      copykat_slide, copykat_sc) %>% # choose order of analyses
                        cor(method="spearman")
    }
    else if (modality == "HTAPP-895-SMP-7359") {
        all_corr <- all_cna_avgs %>% 
                        dplyr::select(-arm) %>%
                        dplyr::select(slidecna_slide, slidecna_sc, 
                                      infercnv_slide, infercnv_sc,
                                      copykat_slide, copykat_sc,
                                     absolute_wes) %>% # choose order of analyses
                        cor(method="spearman")
    }
    else if (modality == "HTAPP_in-silico_sep") {
        all_corr <- all_cna_avgs %>% 
                        dplyr::select(-arm) %>%
                        dplyr::select(HTAPP_895_infercnv_sc,
                                      HTAPP_944_infercnv_sc,
                                      in_silico_sep1,
                                      in_silico_sep2, 
                                      in_silico_sep3) %>% # choose order of analyses
                        cor(method="spearman")
    }
    else if (modality == "HTAPP_in-silico_mix") {
        all_corr <- all_cna_avgs %>% 
                        dplyr::select(-arm) %>%
                        dplyr::select(HTAPP_895_infercnv_sc,
                                      HTAPP_944_infercnv_sc,
                                      in_silico_mix1,
                                      in_silico_mix2, 
                                      in_silico_mix3) %>% # choose order of analyses
                        cor(method="spearman")
    }
    else if (modality == "replicate") {
        all_corr <- all_cna_avgs %>% 
                        dplyr::select(-arm) %>%
                        dplyr::select(slidecna_rep1,
                                      slidecna_rep2,
                                      slidecna_joint) %>% # choose order of analyses
                        cor(method="spearman")
    }

    # Rename rows and columns
    all_corr <- all_corr %>% 
        `colnames<-`(analysis_to_label_abr[colnames(all_corr)]) %>%
        `rownames<-`(analysis_to_label_abr[colnames(all_corr)])
    
    options(repr.plot.width = 8, repr.plot.height = 8)
    ## text labels rotated 45 degrees and  wider color legend with numbers right aligned
    pdf(file = glue("{out_dir}/all_corr_plot_{modality}.pdf"), width = 8, height = 8)
    corrplot(all_corr, 
                     type = 'lower', 
                     order = 'original', 
                     tl.col = 'black',
                     cl.ratio = 0.2, 
                     addCoef.col = 'black',
                     tl.srt = 45,  
                     col = COL2( "PiYG"),
                     diag = FALSE)
    dev.off()
    # print out again to console
    corrplot(all_corr, 
            type = 'lower', 
             order = 'original', 
            tl.col = 'black',
             cl.ratio = 0.2, 
             addCoef.col = 'black',
              tl.srt = 45,  
                col = COL2( "PiYG"),
            diag = FALSE)

}



### Calculate and plot spearman rank correlation of average CNA score per arm between 2 analyses 
cna_correlation <- function(all_cna_avgs, analysis1, analysis2, out_dir) {
    
    options(repr.plot.width = 8, repr.plot.height =8)

    gg_corr <- ggscatter(all_cna_avgs, x = analysis1, y = analysis2, 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                    cor.coef.size=8,
                  xlab = analysis_to_label[analysis1], ylab = analysis_to_label[analysis2]) +
                theme(axis.text=element_text(size=text_size),
                      axis.title=element_text(size=text_size))
    
    pdf(file = glue("{out_dir}/{analysis1}_{analysis2}_corr.pdf"), width = 6, height = 6)
    print(gg_corr)
    dev.off()
    print(gg_corr)
    
    test <- cor.test(all_cna_avgs[,analysis1], all_cna_avgs[,analysis2], method=c("spearman"))
    print(test)
    return(test)
}

