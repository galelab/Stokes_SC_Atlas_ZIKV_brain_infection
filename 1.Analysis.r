#Author Leanne Whitmore
#Description This will be a code comprising the analysis for both
# Caleb Stokes single cell projects iNPC_Zika project (original project ZIKV_07_iNPC) and 
# human fetal brain tissues (SC_FetalBrain)
# as of right now for me this has to be run on Frankline Rv3.6
#################
# Libraries - Analysis run on R v3.6.0
##################
library(dplyr)
library(Seurat) # v3.2.3
library(ggplot2)
# devtools::install_github("cole-trapnell-lab/garnett", ref = "monocle3")
library(monocle3)
library(garnett)
library(stringr)
library(org.Hs.eg.db)
library(SingleCellExperiment)
library(data.table)
library(UpSetR)
library(VennDiagram)
library(ggrepel)
source('SC_DE_Fig.r')

######################
# Functions - general
######################
generate_folder <- function(foldername) {
    workDir <- getwd()
    subDir <- foldername
    results_path <- file.path(workDir, subDir)
    if (file.exists(subDir)) {
    } else {
        dir.create(results_path)
    }
    return(results_path)
}

#####################
#Functions --iNPC
#####################
run_normalizationinpc <- function(combined, SD = TRUE, SCTRANSFORM = FALSE) {
    message("STATUS: normalizing data")
    all.genes <- rownames(combined)
    if (SD == TRUE) {
        combined <- ScaleData(combined,
            vars.to.regress = "percent.mt",
            features = all.genes, verbose = FALSE
        )
    } else if (SCTRANSFORM == TRUE) {
        combined <- SCTransform(combined,
            vars.to.regress = "percent.mt",
            verbose = FALSE
        )
    }
    return(combined)
}


feature_reductioninpc <- function(combined, result_folder) {
    message("STATUS: performing PCA and UMAP")

    combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
    combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)

    DimPlot(combined,
        reduction = "umap",
        group.by = "orig.ident",
        label = FALSE
    )
    ggsave(file.path(result_folder, "UMAP_samples.png"), height = 5, width = 8, dpi = 500)

    DimPlot(combined,
        reduction = "umap",
        group.by = "ZIKA",
        cols = c("#ddd2d2", "#2d0446"),
        label = FALSE
    )
    ggsave(file.path(result_folder, "UMAP_ZIKA.png"), height = 5, width = 8, dpi = 500)


    return(combined)
}

initial_SC_Seurat_inpc <- function(data_dir, sampleID, mock = FALSE, results_path) {
    pbmc.data <- Read10X(data.dir = data_dir)

    pbmc <- CreateSeuratObject(
        counts = pbmc.data,
        project = sampleID,
        min.cells = 0,
        min.features = 200
    )
    # print("STATUS: number of cells before filtering for ", sampleID, " is ", length(pbmc$orig.ident))
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    ggsave(file.path(results_path, paste0(sampleID, "_VnPlotMt.png")), dpi = 300)

    VlnPlot(pbmc, features = c("zika-positive"))
    ggsave(file.path(results_path, paste0(sampleID, "_brazilzika_VnPlotMt.png")), dpi = 300)

    pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 10)
    VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    ggsave(file.path(results_path, paste0(sampleID, "_VnPlotMt_Filtered.png")), dpi = 300)
    # print("STATUS: number of cells after filtering for ", sampleID, " is ", length(pbmc$orig.ident))
    z <- pbmc["zika-positive", ]
    x <- z@meta.data[z@meta.data$nCount_RNA > 0, ]
    x1 <- z@meta.data[z@meta.data$nFeature_RNA == 1, ]
    if (length(rownames(x1)) > 0) {
        pbmczika <- subset(pbmc, cells = rownames(x1))

        if (length(rownames(x1)) > 1) {
            write(paste0("STATUS: FOR SAMPLE ", sampleID, " read count distribution"),
                file.path(results_path, "2.ZikaReadCountdist.txt"),
                append = TRUE
            )
            write(quantile(as.matrix(x1$nCount_RNA)),
                file.path(results_path, "2.ZikaReadCountdist.txt"),
                append = TRUE
            )
        }
        VlnPlot(pbmczika, features = c("zika-positive"))
        ggsave(file.path(results_path, paste0(sampleID, "_brazilzika_VnPlotMt_onlyzikacells.png")), height = 5, width = 5, dpi = 300)

        VlnPlot(pbmczika, slot = "counts", features = c("zika-positive"))
        ggsave(file.path(results_path, paste0(sampleID, "_brazilzika_VnPlotMt_onlyzikacells_counts.png")), height = 5, width = 5, dpi = 300)
    } else {
        message("STATUS: no zika in any of the cells for sample ", sampleID)
    }


    ## CLASSIFY CELLS AS ZIKA OR NO ZIKA HAVE TO HAVE GREATER THAN ONE TRANSCRIPT MAPPING TO THE ZIKA GENOME
    zikadata <- c()
    for (cell in colnames(pbmc)) {
        if (cell %in% rownames(x)) {
            zikadata <- c(zikadata, "zika")
        } else {
            zikadata <- c(zikadata, "nozika")
        }
    }
    pbmc$ZIKA <- zikadata


    message("STATUS: Percentage of cells with any Zika transcript ", dim(x1)[1] / length(pbmc$orig.ident))
    message("STATUS: Percentage of cells Zika transcript at least 10 read counts ", dim(x)[1] / length(pbmc$orig.ident))
    write(paste0("STATUS FOR SAMPLE ", sampleID, ": Percentage of cells with any Zika transcript ", dim(x1)[1] / length(pbmc$orig.ident), " Cells with Zika ", dim(x1)[1], " Total cells ", length(pbmc$orig.ident)),
        file.path(results_path, "2.ZikaExpressingCellStats.txt"),
        append = TRUE
    )
    write(paste0("STATUS FOR SAMPLE ", sampleID, ": Percentage of cells with any Zika transcript at least 10 read counts ", dim(x)[1] / length(pbmc$orig.ident), " Cells with Zika ", dim(x)[1], " Total cells ", length(pbmc$orig.ident)),
        file.path(results_path, "2.ZikaExpressingCellStats.txt"),
        append = TRUE
    )

    if (isTRUE(mock)) {
        message("STATUS: removing Zika infected cells from analysis")
        cells <- pbmc@meta.data[pbmc@meta.data$ZIKA == "zika", ]
        message("STATUS: number of HIV cells being removed ", length(rownames(cells)))
        if (length(cells) > 0) {
            nohivcells <- setdiff(rownames(pbmc@meta.data), rownames(cells))
            pbmc <- subset(pbmc, cells = nohivcells)
        }
    }

    # Generate Monocle OBject
    cds <- generate_monocle_cds(str_remove(data_dir, "outs/filtered_feature_bc_matrix"))

    cds <- cds[, colnames(pbmc)]
    pbmc <- NormalizeData(pbmc,
        normalization.method = "LogNormalize",
        scale.factor = 10000
    )


    # Find features (genes) that are highly variable from cell-to-cell
    pbmc <- FindVariableFeatures(pbmc,
        selection.method = "vst",
        nfeatures = 2000
    )

    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(pbmc), 10)

    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(pbmc)
    LabelPoints(plot = plot1, points = top10, repel = TRUE)
    ggsave(file.path(results_path, paste0(sampleID, "_HighlyVariableGenes.png")), dpi = 300)

    pbmc <- ScaleData(pbmc)
    pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
    pbmc <- RunUMAP(pbmc, dims = 1:10)

    z <- pbmc["IFNB1", ]
    x <- z@meta.data[z@meta.data$nCount_RNA > 0, ]
    x1 <- z@meta.data[z@meta.data$nFeature_RNA == 1, ]
    if (length(rownames(x1)) > 0) {
        pbmcifnb1 <- subset(pbmc, cells = rownames(x1))
        VlnPlot(pbmcifnb1, features = c("IFNB1"))
        ggsave(file.path(results_path, paste0(sampleID, "_IFNB1_VnPlotMt_onlyIFNB1cells.png")), height = 5, width = 5, dpi = 300)

        VlnPlot(pbmcifnb1, slot = "counts", features = c("IFNB1"))
        ggsave(file.path(results_path, paste0(sampleID, "_IFNB1_VnPlotMt_onlyIFNB1cells_counts.png")), height = 5, width = 5, dpi = 300)
    }
    VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
    ggsave(file.path(results_path, paste0(sampleID, "_PCA_Loadings.png")), dpi = 300)

    DimPlot(pbmc, reduction = "pca")
    ggsave(file.path(results_path, paste0(sampleID, "_PCA.png")), dpi = 300)

    DimPlot(pbmc, reduction = "umap")
    ggsave(file.path(results_path, paste0(sampleID, "_UMAP.png")), dpi = 300)

    FeaturePlot(pbmc, reduction = "umap", features = c("zika-positive"))
    ggsave(file.path(results_path, paste0(sampleID, "_brazilzika.png")), dpi = 300)

    FeaturePlot(pbmc, reduction = "umap", features = c("IFNB1"))
    ggsave(file.path(results_path, paste0(sampleID, "_IFNB1.png")), dpi = 300)

    VlnPlot(pbmc, features = c("IFNB1"))
    ggsave(file.path(results_path, paste0(sampleID, "_IFNB1Vnplot.png")), dpi = 300)

    z <- pbmc["IFNB1", ]
    x <- z@meta.data[z@meta.data$nCount_RNA > 1, ]
    x1 <- z@meta.data[z@meta.data$nFeature_RNA == 1, ]
    if (length(rownames(x1)) > 0) {
        pbmcifnb <- subset(pbmc, cells = rownames(x1))

        VlnPlot(pbmcifnb, features = c("IFNB1", "zika-positive", "AQP4", "OLIG1", "SNAP25", "BIRC5"), group.by = "orig.ident", same.y.lims = TRUE)
        ggsave(file.path(results_path, paste0(sampleID, "1_VnPlotMt_onlyIFNB1cells.png")), dpi = 300)

        DoHeatmap(pbmcifnb, features = c("IFNB1", "zika-positive", "AQP4", "OLIG1", "SNAP25", "BIRC5"), group.by = "orig.ident", draw.lines = FALSE) + scale_fill_gradientn(colors = c("black", "#eb3306")) +
            theme(panel.grid = element_line(color = "white"))
        ggsave(file.path(results_path, paste0(sampleID, "_Heatmap_onlyIFNB1cells.png")), width = 4, height = 4, units = "in", dpi = 300)
    }
    DoHeatmap(pbmc, features = c("IFNB1", "zika-positive", "AQP4", "OLIG1", "SNAP25", "BIRC5"), group.by = "orig.ident", draw.lines = FALSE) + scale_fill_gradientn(colors = c("black", "#eb3306")) +
        theme(panel.grid = element_line(color = "white"))
    ggsave(file.path(results_path, paste0(sampleID, "_Heatmap_IFNB1cells.png")), width = 4, height = 4, units = "in", dpi = 300)


    # FeaturePlot(pbmc, reduction = "umap", features = c("PAX6", "EOMES", "SOX2", "SLC1A3", "zika-positive"), cols = c("#e6d8d8", "#eb3306"))
    # ggsave(file.path(results_path, paste0(sampleID, "_markersNPC.png")), dpi = 300)

    # FeaturePlot(pbmc, reduction = "umap", features = c("BIRC5", "CCNB1", "zika-positive"), cols = c("#e6d8d8", "#eb3306"))
    # ggsave(file.path(results_path, paste0(sampleID, "_markersNPCmaybe.png")), dpi = 300)

    # FeaturePlot(pbmc, reduction = "umap", features = c("SCG2", "RELN", "SYNPR", "VIP", "SNAP25", "SCN2A", "zika-positive"), cols = c("#e6d8d8", "#eb3306"))
    # ggsave(file.path(results_path,  paste0(sampleID, "_markersNeurons.png")), dpi = 300)

    # FeaturePlot(pbmc, reduction = "umap", features = c("AQP4", "SLCA4", "GJB6", "GJA1", "S100B", "GFAP", "zika-positive"), cols = c("#e6d8d8", "#eb3306"))
    # ggsave(file.path(results_path,  paste0(sampleID, "_markersAstrocyte.png")), dpi = 300)

    # FeaturePlot(pbmc, reduction = "umap", features = c("NEUROD6", "zika-positive"), cols = c("#e6d8d8", "#eb3306"))
    # ggsave(file.path(results_path, paste0(sampleID, "_markersExNeurons.png")), width = 8, height = 4, units = "in", dpi = 300)

    # FeaturePlot(pbmc, reduction = "umap", features = c("SLC17A6", "zika-positive"), cols = c("#e6d8d8", "#eb3306"))
    # ggsave(file.path(results_path, paste0(sampleID, "_markersProjectNeurons.png")), width = 8, height = 4, units = "in", dpi = 300)

    # FeaturePlot(pbmc, reduction = "umap", features = c("GAD1", "GAD2", "zika-positive"), cols = c("#e6d8d8", "#eb3306"))
    # ggsave(file.path(results_path, paste0(sampleID, "_markersInNeurons.png")), dpi = 300)

    # FeaturePlot(pbmc, reduction = "umap", features = c(Oligodendrocytemarkers, c("zika-positive")), cols = c("#e6d8d8", "#eb3306"))
    # ggsave(file.path(results_path, paste0(sampleID, "_markersOligodendrocyte.png")), dpi = 300)

    # FeaturePlot(pbmc, reduction = "umap", features = c("PADI2", "PAD2", "zika-positive"), cols = c("#e6d8d8", "#eb3306"))
    # ggsave(file.path(results_path, paste0(sampleID, "_PADI2.png")), dpi = 300)


    pbmc <- FindNeighbors(pbmc, dims = 1:10)
    pbmc <- FindClusters(pbmc, resolution = 0.4) # Default values ued
    DimPlot(pbmc, reduction = "umap")
    ggsave(file.path(results_path, paste0(sampleID, "_clusters.png")), dpi = 300)

    # unique_clusters <- unique(Idents(pbmc))
    # for (cluster in unique_clusters) {
    #     cluster.markers <- FindMarkers(pbmc, ident.1 = cluster, min.pct = 0.25)
    #     write.csv(head(cluster.markers, n = 5), paste0(results_path, sampleID, "cluster", cluster, "markers.csv"))
    # }

    return(list("seurat" = pbmc, "monocle" = cds))
}

get_percentages <- function(tabletemp_orig) {
    tabletemp <- data.frame(table(tabletemp_orig))
    tabletemp$Freq <- tabletemp$Freq / length(tabletemp_orig)
    return(tabletemp)
}
get_zika_percentages <- function(tabletemp_orig, combined, sample) {
    zika_final <- data.frame()
    for (cell in unique(tabletemp_orig)) {
        x <- which(tabletemp_orig == cell)
        total_per <- length(x) / length(tabletemp_orig)
        z <- combined$ZIKA[names(x)]
        zika <- as.data.frame(table(z))
        zika$Freq <- (zika$Freq / length(x)) * total_per
        zika$cell <- rep(cell, length(zika$Freq))
        zika$sample <- rep(sample, length(zika$Freq))
        zika_final <- rbind(zika_final, zika)
    }
    return(zika_final)
}

get_zika_percentages_total <- function(tabletemp_orig, sample, combined) {
    zika_final <- data.frame()
    for (cell in unique(tabletemp_orig)) {
        x <- which(tabletemp_orig == cell)
        total_per <- length(x) / length(tabletemp_orig)
        z <- combined$ZIKA[names(x)]
        zika <- as.data.frame(table(z))
        zika$Freq <- (zika$Freq / length(x))
        zika$cell <- rep(cell, length(zika$Freq))
        zika$sample <- rep(sample, length(zika$Freq))
        zika_final <- rbind(zika_final, zika)
    }
    return(zika_final)
}

num_cells_expressing_gene <- function(gene, combined, results_folder) {
    z <- combined[gene, ]
    x1 <- z@meta.data[z@meta.data$nFeature_RNA == 1, ]
    t <- table(x1$garnett_cluster_extend_lw)
    # z <- combined_sd_sub@assays$RNA["IFNB1", ]
    T <- x1 %>%
        group_by(orig.ident, garnett_cluster_extend_lw, ZIKA) %>%
        summarise(Count = n())
    T <- data.frame(T)
    write.csv(T, file.path(results_folder, paste0(gene, "_table.csv")))
    ggplot(data = T, aes(
        x = factor(garnett_cluster_extend_lw, levels = c("NPC", "Early Neurons", "Neurons", "Astrocytes", "Oligodendrocyte")),
        y = Count, fill = factor(ZIKA, levels = c("nozika", "zika"))
    )) +
        geom_bar(stat = "identity") +
        theme_Publication() +
        facet_wrap(~ factor(orig.ident), ncol = 2) +
        theme(strip.background = element_rect(fill = c("white"))) +
        theme(legend.text = element_text(size = 8)) +
        scale_color_manual(values = c("#ddd2d2", "black")) +
        scale_fill_manual(values = c("#ddd2d2", "black")) +
        labs(x = "", fill = "Zika", title = gene) +
        # ylim(c(0, 1)) +
        ylab("# of Cells") +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
            axis.text.y = element_text(size = 12)
        )
    ggsave(file.path(results_folder, paste0("CellCount_", gene, "_zika.png")), height = 5.5, width = 6, dpi = 500)
    ggsave(file.path(results_folder, paste0("CellCount_", gene, "_zika.svg")), height = 5.5, width = 6, dpi = 500)
    ggsave(file.path(results_folder, paste0("CellCount_", gene, "_zika.pdf")), height = 5.5, width = 6, dpi = 500)

}

get_DE_between_conditions <- function(ident_1, ident_2, compare,
                                      combined,
                                      result_folder, fontsize = 9, height = 8,
                                      foldchange = 0.26,
                                      pval = 0.05, percent_cells = NULL) {
    message("STATUS: getting DEs...")
    ### Parameters for FindMarkers
    #### test.use: mast (default is wilcox)
    #### min.pct: 0.1(default: 0.1) filter out genes (features) that are detected at less than 10 percent frequency in cells in ident_1 or ident_2
    #### logfc: 0 (default is .25) logfc must be higher than 0 (I set this at 0 because I filter this out at a later stage - this was primarily done to
    ####                          understand how the tool (Seurat) works
    #### min.cells.feature: 3 (default is 3) minimum number of cells expressing the feature in at least one of the two groups (similar to min.pct)
    #### min.cells.group: 3 (defualt is 3) minimum number of cells in the group
    #### max.cells.per.ident: Inf (default Inf-means no down sampling) Down sample each identity class (cluster of cells) to a max number
    #### min.dif.pct: Inf (default Inf) Only test genes that show minimum difference in the fraction of detection between the two identities (cluster)
    if (file.exists(file.path(result_folder, paste0("DEgenes_full_", ident_1, "_", ident_2, "_", compare, ".txt")))) {
        message("STATUS: DE has already been done filtering from full file ")
        DEgenes <- read.table(file.path(result_folder, paste0("DEgenes_full_", ident_1, "_", ident_2, "_", compare, ".txt")), header = TRUE, row.names = 1)
    } else {
        DEgenes <- FindMarkers(combined,
            ident.1 = ident_1,
            ident.2 = ident_2,
            test.use = "MAST",
            logfc.threshold = 0
        )
        write.table(data.frame(DEgenes),
            file.path(
                result_folder,
                paste0(
                    "DEgenes_full_",
                    ident_1,
                    "_", ident_2, "_", compare,
                    ".txt"
                )
            ),
            sep = "\t", quote = FALSE
        )
    }

    DEgenes[["gene_name"]] <- rownames(DEgenes)

    DE_sig_final <- DEgenes %>%
        filter(avg_logFC >= foldchange | avg_logFC <= -foldchange) %>%
        dplyr::select(gene_name, p_val, avg_logFC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
    DE_sig_final <- DE_sig_final %>%
        filter(p_val_adj <= pval) %>%
        dplyr::select(gene_name, p_val, avg_logFC, pct.1, pct.2, p_val_adj)
    if (!is.null(percent_cells)) {
        DE_sig_final <- DE_sig_final %>%
            filter(pct.1 >= percent_cells | pct.2 >= percent_cells) %>%
            dplyr::select(gene_name, p_val, avg_logFC, pct.1, pct.2, p_val_adj)
    }
    rownames(DE_sig_final) <- DE_sig_final$gene_name

    DE_sig_final$gene_name <- NULL
    if (is.null(percent_cells)) {
        write.table(data.frame(DE_sig_final),
            file.path(
                result_folder,
                paste0(
                    "DEgenes_sig_",
                    ident_1, "_", ident_2, "_", compare,
                    ".txt"
                )
            ),
            sep = "\t", quote = FALSE
        )
    } else {
        write.table(data.frame(DE_sig_final),
            file.path(
                result_folder,
                paste0(
                    "DEgenes_sig_",
                    ident_1, "_", ident_2, "_", compare,
                    "_", percent_cells, ".txt"
                )
            ),
            sep = "\t", quote = FALSE
        )
    }
    if (length(rownames(DE_sig_final)) > 0) {
        DotPlot(combined,
            idents = c(ident_1, ident_2), features = rownames(DE_sig_final),
            cols = c("blue", "#eb3306"), dot.scale = 3
        ) + coord_flip() +
            theme(axis.text.y = element_text(size = fontsize))
        if (is.null(percent_cells)) {
            ggsave(file.path(
                result_folder, paste0(
                    "DEgenes_sig_", ident_1, "_",
                    ident_2, "_", compare, ".png"
                )
            ), height = 8, width = 6, units = "in", dpi = 350)
        } else {
            ggsave(file.path(
                result_folder, paste0(
                    "DEgenes_sig_", ident_1, "_",
                    ident_2, "_", compare, "_", percent_cells, ".png"
                )
            ), height = 8, width = 6, units = "in", dpi = 350)
        }
    }
    if (length(rownames(DE_sig_final)) > 1) {
        message("STATUS: Making BarPlot")
        pl <- SC_DE_barplot(DE_sig_final, ident_1, ident_2, fontsize = fontsize)
        if (is.null(percent_cells)) {
            ggsave(file.path(
                result_folder, paste0(
                    "DEgenes_sig_", ident_1, "_",
                    ident_2, "_", compare, "barPlot.png"
                )
            ),
            width = 5, height = height, units = "in", dpi = 300
            )
        } else {
            ggsave(file.path(
                result_folder, paste0(
                    "DEgenes_sig_", ident_1, "_",
                    ident_2, "_", compare, "_", percent_cells, "barPlot.png"
                )
            ),
            width = 5, height = height, units = "in", dpi = 300
            )
        }
    }
    return(DE_sig_final)
}

generate_monocle_cds <- function(sample_path) {
    SAMPLE <- load_cellranger_data(sample_path, umi_cutoff = 200)

    SAMPLE.counts <- SAMPLE@assays@data@listData$counts
    gene_meta_data <- rowData(SAMPLE)
    cell_meta_data <- colData(SAMPLE)

    cds <- new_cell_data_set(SAMPLE.counts,
        cell_metadata = cell_meta_data,
        gene_metadata = gene_meta_data
    )
    return(cds)
}
####################
#iNPC analsyis
####################

# --output folder for iNPC -- #
inpc_results <- "inpc_resultsv2"
generate_folder(inpc_results)

# path to aligned files
SAMPLE_PATH <- "/vol08/ngs/P51/ZIKV/ZIKV_07_iNPC/WhitmoreAnalysis/NWGC_GaleLab-09/"

##### ---Run individual samples of inpc data--#####

results_folder <- generate_folder(file.path(inpc_results, "single_analysis_results"))

G001 <- initial_SC_Seurat_inpc(file.path(SAMPLE_PATH, "G001_CS132_iNPCs_GM_Mock_IFN_cDNA_Lib_a_GL09_wbz/outs/filtered_feature_bc_matrix"),
    sampleID = "GM_NP_Mock", mock=TRUE, results_folder
)
G002 <- initial_SC_Seurat_inpc(file.path(SAMPLE_PATH, "/G002_CS132_iNPCs_GM_IFNb_IFN_cDNA_Lib_a_GL09_wbz/outs/filtered_feature_bc_matrix"),
    sampleID = "GM_NP_IFNb", mock=TRUE, results_folder
)
G003 <- initial_SC_Seurat_inpc(file.path(SAMPLE_PATH, "/G003_CS132_iNPCs_DD_Mock_IFN_cDNA_Lib_a_GL09_wbz/outs/filtered_feature_bc_matrix"),
    sampleID = "DD_NP_Mock", mock=TRUE, results_folder
)
G004 <- initial_SC_Seurat_inpc(file.path(SAMPLE_PATH, "/G004_CS132_iNPCs_DD_IFNb_IFN_cDNA_Lib_a_GL09_wbz/outs/filtered_feature_bc_matrix"),
    sampleID = "DD_NP_IFNb", mock=TRUE, results_folder
)
G005 <- initial_SC_Seurat_inpc(file.path(SAMPLE_PATH, "/G005_CS132_iNPCs_GM_Mock_Inf_cDNA_Lib_a_GL09_wbz/outs/filtered_feature_bc_matrix"),
    sampleID = "GM_NP_Mock2zika", mock=TRUE, results_folder
)
G006 <- initial_SC_Seurat_inpc(file.path(SAMPLE_PATH, "/G006_CS132_iNPCs_GM_ZIKV_Inf_cDNA_Lib_a_GL09_wbz/outs/filtered_feature_bc_matrix"),
    sampleID = "GM_NP_ZIKV", mock=FALSE, results_folder
)
G007 <- initial_SC_Seurat_inpc(file.path(SAMPLE_PATH, "/G007_CS132_iNPCs_DD_Mock_Inf_cDNA_Lib_a_GL09_wbz/outs/filtered_feature_bc_matrix"),
    sampleID = "DD_NP_Mock2zika", mock=TRUE, results_folder
)
G008 <- initial_SC_Seurat_inpc(file.path(SAMPLE_PATH, "/G008_CS132_iNPCs_DD_ZIKV_Inf_cDNA_Lib_a_GL09_wbz/outs/filtered_feature_bc_matrix"),
    sampleID = "DD_NP_ZIKV", mock=FALSE, results_folder)


##### ---USE GARNETT TO CLASSIFY CELL TYPES---#####

results_folder <- generate_folder(file.path(inpc_results, "garnett_classification_inpc_results"))

#------build monocle CDS list
batch.name <- c("G001", "G002", "G003", "G004", "G005", "G006", "G007", "G008")
cds.list <- list(
    G001$monocle, G002$monocle, G003$monocle, G004$monocle,
    G005$monocle, G006$monocle, G007$monocle, G008$monocle
)
names(cds.list) <- batch.name
big_cds <- combine_cds(cds.list)
saveRDS(big_cds, file.path(results_folder, "2.big_cds.rds"))
# big_cds <- readRDS(file.path(results_folder, "2.big_cds.rds"))

## ------Align/Normalize CDS object
big_cds <- preprocess_cds(big_cds, num_dim = 100, preprocess_method = 100)
big_cds <- align_cds(big_cds, alignment_group = "sample")

big_cds <- reduce_dimension(big_cds,
    preprocess_method = "PCA",
    reduction_method = "UMAP"
)
## ------Build training data set using markers
marker_file_path <- "iNPCmarkers.txt"
npc_classifier <- train_cell_classifier(
    cds = big_cds,
    marker_file = marker_file_path,
    db = org.Hs.eg.db,
    cds_gene_id_type = "ENSEMBL",
    num_unknown = 500,
    marker_file_gene_id_type = "SYMBOL"
)
saveRDS(npc_classifier,  file.path(results_folder, "2.npc_classifier.rds"))
npc_classifier <- readRDS( file.path(results_folder, "2.npc_classifier.rds"))

## ------Classify cells
big_cds1 <- classify_cells(big_cds, npc_classifier,
    db = org.Hs.eg.db,
    cluster_extend = TRUE,
    cds_gene_id_type = "ENSEMBL"
)
saveRDS(big_cds1,  file.path(results_folder, "2.big_cds_clusterext.rds"))
# big_cds1 <- readRDS( file.path(results_folder, "2.big_cds_clusterext.rds"))

## ------Generate Marker files
marker_check <- check_markers(big_cds1, marker_file_path,
    db = org.Hs.eg.db,
    cds_gene_id_type = "ENSEMBL",
    marker_file_gene_id_type = "SYMBOL"
)
plot_markers(marker_check)
ggsave(file.path(results_folder,"2.Marker_Fig.png"), dpi = 500)

##### ---Integrate inpc data ---#####
results_folder <- generate_folder(file.path(inpc_results, "integrated_results"))

# --Compile a list of all samples 
sample.list <- list(
    "GM_NP_Mock" = G001$seurat, "GM_NP_IFNb" = G002$seurat,
    "DD_NP_Mock" = G003$seurat, "DD_NP_IFNb" = G004$seurat,
    "GM_NP_Mock2zika" = G005$seurat, "GM_NP_ZIKV" = G006$seurat,
    "DD_NP_Mock2Zika" = G007$seurat, "DD_NP_ZIKV" = G008$seurat
)
inpc.anchors <- FindIntegrationAnchors(object.list = sample.list, dims = 1:20)
inpc.combined <- IntegrateData(anchorset = inpc.anchors, dims = 1:20)
saveRDS(inpc.combined, file.path(results_folder, "inpc.combinedobject.rds"))
# inpc.combined <- readRDS(file.path(results_folder, "inpc.combinedobject.rds"))

DefaultAssay(inpc.combined) <- "integrated"
inpc.combined_sd <- run_normalizationinpc(inpc.combined)

# DefaultAssay(inpc.combined_sd) <- "RNA"
# inpc.combined_sd <- run_normalization(inpc.combined_sd)
saveRDS(inpc.combined_sd, file.path(results_folder, "inpc.combinedobject_norm.rds"))
# inpc.combined_sd <- readRDS(file.path(results_folder, "inpc.combinedobject_norm.rds"))

## --perform feature reduction on seurat object
DefaultAssay(inpc.combined_sd) <- "integrated"
inpc.combined_sd <- feature_reductioninpc(inpc.combined_sd, results_folder)

# saveRDS(inpc.combined_sd, file.path(results_folder, "inpc.combinedobject_norm_fr.rds"))
inpc.combined_sd <- readRDS(file.path(results_folder, "inpc.combinedobject_norm_fr.rds"))

##### ---Extend clusters from garnett  inpc data ---#####
results_folder <- generate_folder(file.path(inpc_results, "classification_extend_results"))

# ------add cell type classifications to integrated seurat object
g_cluster_ext <- big_cds1$cluster_ext_type
inpc.combined_sd[["garnett_cluster_extend"]] <- g_cluster_ext
Idents(inpc.combined_sd) <- inpc.combined_sd$garnett_cluster_extend
DimPlot(inpc.combined_sd, reduction = "umap")
ggsave(file.path(results_folder, "garnett_clusters.png"), width = 6.5, height = 4, dpi = 500)

## ------find clusters in integrated seurat object
inpc.combined_sd <- FindNeighbors(inpc.combined_sd, dims = 1:20)
inpc.combined_sd <- FindClusters(inpc.combined_sd, resolution = 0.4) # Default values ued
Idents(inpc.combined_sd) <- inpc.combined_sd$seurat_clusters
DimPlot(inpc.combined_sd, reduction = "umap")
ggsave(file.path(results_folder, "seruat_clusters.png"), width = 6.5, height = 4, dpi = 500)

## ------for cells with the unknown classification assign them the cell type that is most prevelant in a seurat cluster
inpc.combined_sd[["garnett_cluster_extend_lw"]] <- inpc.combined_sd$garnett_cluster_extend
for (cluster in unique(inpc.combined_sd$seurat_clusters)) {
    celltype_count <- c()
    clust_list <- which(inpc.combined_sd$seurat_clusters == cluster)
    tempdf <- inpc.combined_sd@meta.data[clust_list, ]
    celltypes <- unique(tempdf$garnett_cluster_extend)
    for (celltype in celltypes) {
        if (celltype != "Unknown") {
            celltype_count[celltype] <- length(which(tempdf$garnett_cluster_extend == celltype))
        }
    }
    maxvalue <- max(celltype_count)
    celltypefinal <- names(which(celltype_count == maxvalue))
    message("STATUS: max cell count in cluster ", cluster, " is ", celltypefinal)
    tempdfunknown <- inpc.combined_sd@meta.data[inpc.combined_sd@meta.data$garnett_cluster_extend == "Unknown", ]
    tempdfunknown <- tempdfunknown[tempdfunknown$seurat_clusters == cluster, ]
    message("STATUS number of unknowns in cluster ", cluster, " is ", dim(tempdfunknown)[1])
    for (row in rownames(tempdfunknown)) {
        inpc.combined_sd@meta.data[row, "garnett_cluster_extend_lw"] <- celltypefinal
    }
}
Idents(inpc.combined_sd) <- inpc.combined_sd$garnett_cluster_extend_lw
DimPlot(inpc.combined_sd, reduction = "umap")
ggsave(file.path(results_folder, "garnett_clusters_lw.png"), width = 6.5, height = 4, dpi = 500)
saveRDS(inpc.combined_sd, file.path(results_folder, "inpc.combinedobject_norm_fr_celltype.rds"))
# inpc.combined_sd <- readRDS(file.path(results_folder, "inpc.combinedobject_norm_fr_celltype.rds"))

##--Generate UMAPs 4 Caleb publication --##
cells <- inpc.combined_sd@meta.data[inpc.combined_sd@meta.data$orig.ident == "GM_NP_Mock2zika" |
    inpc.combined_sd@meta.data$orig.ident == "GM_NP_ZIKV" | inpc.combined_sd@meta.data$orig.ident == "DD_NP_Mock2zika" |
    inpc.combined_sd@meta.data$orig.ident == "DD_NP_ZIKV", ]
inpc.combined_sd_sub <- inpc.combined_sd[, rownames(cells)]

inpc.combined_sd_sub$orig.ident <- factor(x = inpc.combined_sd_sub$orig.ident, levels = c("GM_NP_Mock2zika", "GM_NP_ZIKV", "DD_NP_Mock2zika", "DD_NP_ZIKV"))

DimPlot(inpc.combined_sd_sub,
    reduction = "umap",
    group.by = "garnett_cluster_extend_lw",
    split.by = "orig.ident",
    ncol = 2,
    # cells.highlight = rownames(cells),
    # sizes.highlight = 0.01,
    # cols = c("#ddd2d2", "#2d0446"),
    label = FALSE
) + scale_color_manual(values = c("NPC" = "#CCBB44", "Early Neurons" = "#66CCEE", "Neurons" = "#4477AA", "Astrocytes" = "#228833", "Oligodendrocyte" = "#AA3377")) +
    theme(legend.text = element_text(size = 7)) +
    labs(x = "UMAP 1", y = "UMAP 2", title = "") +
    theme(axis.title.x = element_text(size = 13)) +
    theme(axis.title.y = element_text(size = 13)) +
    theme(axis.text = element_text(size = 10)) +
    # theme(legend.key.size = unit(0.2, "cm")) +
    theme(legend.position = "none") +
    theme(strip.background = element_blank(), strip.text.x = element_blank())
# theme(legend.direction = "horizontal")
# guides(colour = guide_legend(override.aes = list(size = 1)))
ggsave(file.path(results_folder, "Celltype_UMAP4pub.png"), width = 4, height = 4, dpi = 500)
ggsave(file.path(results_folder, "Celltype_UMAP4pub.svg"), width = 4, height = 4, dpi = 500)
ggsave(file.path(results_folder, "Celltype_UMAP4pub.pdf"), width = 4, height = 4, dpi = 500)

zika_relabeled <- str_replace_all(inpc.combined_sd@meta.data$ZIKA, "nozika", "ZIKA RNA -")
zika_relabeled <- str_replace_all(zika_relabeled, "^zika$", "ZIKA RNA +")
inpc.combined_sd@meta.data$zika_relabeled <- zika_relabeled
cells <- inpc.combined_sd@meta.data[
    inpc.combined_sd@meta.data$ZIKA == "zika",
]
inpc.combined_sd_sub <- subset(inpc.combined_sd, subset = orig.ident == "GM_NP_ZIKV" |
    orig.ident == "DD_NP_ZIKV")
DimPlot(inpc.combined_sd_sub,
    reduction = "umap",
    split.by = "orig.ident",
    cells.highlight = rownames(cells),
    sizes.highlight = 0.01, ncol = 1,
    # cols = c("#ddd2d2", "#2d0446"),
    label = FALSE
) + scale_color_manual(labels = c("ZIKA RNA -", "ZIKA RNA +"), values = c("#ddd2d2", "black")) +
    theme(
        legend.text = element_text(size = 7)
    ) +
    labs(x = "UMAP 1", y = "UMAP 2", title = "") +
    theme(axis.title.x = element_text(size = 13)) +
    theme(axis.title.y = element_text(size = 13)) +
    theme(axis.text = element_text(size = 10)) +
    theme(legend.key.size = unit(0.2, "cm")) +
    theme(legend.position = "none")
# theme(legend.direction = "horizontal") +
# guides(colour = guide_legend(override.aes = list(size = 1)))
ggsave(file.path(results_folder, "ZIKA_UMAP4pub.png"), width = 4.5, height = 8, dpi = 500)
ggsave(file.path(results_folder, "ZIKA_UMAP4pub.svg"), width = 4.5, height = 8, dpi = 500)
ggsave(file.path(results_folder, "ZIKA_UMAP4pub.pdf"), width = 4.5, height = 8, dpi = 500)

inpc.combined_sd_sub_gm <- subset(inpc.combined_sd, subset = orig.ident == "GM_NP_ZIKV")
inpc.combined_sd_sub_dd <- subset(inpc.combined_sd, subset = orig.ident == "DD_NP_ZIKV")
VlnPlot(inpc.combined_sd_sub_gm,
    features = c("zika-positive"), pt.size=0.5,
    group.by = "garnett_cluster_extend_lw",
    assay = "RNA", slot = "counts", ncol = 1
) + scale_fill_manual(values = c(
    "NPC" = "#CCBB44", "Early Neurons" = "#66CCEE", "Neurons" = "#4477AA",
    "Astrocytes" = "#228833", "Oligodendrocyte" = "#AA3377"
)) + labs(title="")+ ylab("Zika UMI count")
ggsave(file.path(results_folder, "ZIKA_UMIcounts4pubGM.png"), width = 6.5, height = 4.5, dpi = 500)
ggsave(file.path(results_folder, "ZIKA_UMIcounts4pubGM.svg"), width = 6.5, height = 4.5, dpi = 500)
ggsave(file.path(results_folder, "ZIKA_UMIcounts4pubGM.pdf"), width = 6.5, height = 4.5, dpi = 500)

VlnPlot(inpc.combined_sd_sub_dd,
    features = c("zika-positive"), pt.size = 0.5,
    group.by = "garnett_cluster_extend_lw",
    assay = "RNA", slot = "counts", ncol = 1
) + scale_fill_manual(values = c(
    "NPC" = "#CCBB44", "Early Neurons" = "#66CCEE", "Neurons" = "#4477AA",
    "Astrocytes" = "#228833", "Oligodendrocyte" = "#AA3377"
)) + labs(title = "") + ylab("Zika UMI count") +ylim(0, 5000)
ggsave(file.path(results_folder, "ZIKA_UMIcounts4pubDD.png"), width = 6.5, height = 4.5, dpi = 500)
ggsave(file.path(results_folder, "ZIKA_UMIcounts4pubDD.svg"), width = 6.5, height = 4.5, dpi = 500)
ggsave(file.path(results_folder, "ZIKA_UMIcounts4pubDD.pdf"), width = 6.5, height = 4.5, dpi = 500)
##### --- generate plots of cell type percentages ---#####

## pull out cells for each sample
Idents(inpc.combined_sd) <- inpc.combined_sd$orig.ident
GM_NP_Mock <- WhichCells(object = inpc.combined_sd, idents = "GM_NP_Mock")
GM_NP_IFNb <- WhichCells(object = inpc.combined_sd, idents = "GM_NP_IFNb")
GM_NP_Mock2zika <- WhichCells(object = inpc.combined_sd, idents = "GM_NP_Mock2zika")
GM_NP_ZIKV <- WhichCells(object = inpc.combined_sd, idents = "GM_NP_ZIKV")
DD_NP_Mock <- WhichCells(object = inpc.combined_sd, idents = "DD_NP_Mock")
DD_NP_IFNb <- WhichCells(object = inpc.combined_sd, idents = "DD_NP_IFNb")
DD_NP_Mock2zika <- WhichCells(object = inpc.combined_sd, idents = "DD_NP_Mock2zika")
DD_NP_ZIKV <- WhichCells(object = inpc.combined_sd, idents = "DD_NP_ZIKV")

## pull out cell type informatin for each sample
table1 <- inpc.combined_sd$garnett_cluster_extend_lw[names(inpc.combined_sd$garnett_cluster_extend_lw) %in% GM_NP_Mock]
table2 <- inpc.combined_sd$garnett_cluster_extend_lw[names(inpc.combined_sd$garnett_cluster_extend_lw) %in% GM_NP_IFNb]
table3 <- inpc.combined_sd$garnett_cluster_extend_lw[names(inpc.combined_sd$garnett_cluster_extend_lw) %in% GM_NP_Mock2zika]
table4 <- inpc.combined_sd$garnett_cluster_extend_lw[names(inpc.combined_sd$garnett_cluster_extend_lw) %in% GM_NP_ZIKV]
table5 <- inpc.combined_sd$garnett_cluster_extend_lw[names(inpc.combined_sd$garnett_cluster_extend_lw) %in% DD_NP_Mock]
table6 <- inpc.combined_sd$garnett_cluster_extend_lw[names(inpc.combined_sd$garnett_cluster_extend_lw) %in% DD_NP_IFNb]
table7 <- inpc.combined_sd$garnett_cluster_extend_lw[names(inpc.combined_sd$garnett_cluster_extend_lw) %in% DD_NP_Mock2zika]
table8 <- inpc.combined_sd$garnett_cluster_extend_lw[names(inpc.combined_sd$garnett_cluster_extend_lw) %in% DD_NP_ZIKV]

table_gm <- inpc.combined_sd$garnett_cluster_extend_lw[names(inpc.combined_sd$garnett_cluster_extend_lw) %in% c(GM_NP_Mock2zika, GM_NP_ZIKV)]
table_dd <- inpc.combined_sd$garnett_cluster_extend_lw[names(inpc.combined_sd$garnett_cluster_extend_lw) %in% c(DD_NP_Mock2zika, DD_NP_ZIKV)]

## ------get zika percentages for each cell type
zika_GM <- get_zika_percentages(table4,inpc.combined_sd,  "GM_NP_ZIKV")
zika_DD <- get_zika_percentages(table8, inpc.combined_sd, "DD_NP_ZIKV")
total_zika <- rbind(zika_GM, zika_DD)

zika_GM <- get_zika_percentages_total(table4, "GM_NP_ZIKV", inpc.combined_sd)
zika_DD <- get_zika_percentages_total(table8, "DD_NP_ZIKV", inpc.combined_sd)
total_zika_100 <- rbind(zika_GM, zika_DD)

table4zika <- inpc.combined_sd$ZIKA[names(inpc.combined_sd$ZIKA) %in% GM_NP_ZIKV]
gm_zika <- length(which(table4zika == "zika"))
gm_nozika <- length(which(table4zika == "nozika"))
table8zika <- inpc.combined_sd$ZIKA[names(inpc.combined_sd$ZIKA) %in% DD_NP_ZIKV]
dd_zika <- length(which(table8zika == "zika"))
dd_nozika <- length(which(table8zika == "nozika"))

# get percentages for each cell type (no zika info)
table1 <- get_percentages(table1)
table2 <- get_percentages(table2)
table3 <- get_percentages(table3)
table4 <- get_percentages(table4)
table5 <- get_percentages(table5)
table6 <- get_percentages(table6)
table7 <- get_percentages(table7)
table8 <- get_percentages(table8)

table_gm <- get_percentages(table_gm)
table_dd <- get_percentages(table_dd)

# integrate cell type information
countstable <- merge(table1, table2, by.x = "tabletemp_orig", by.y = "tabletemp_orig", all = TRUE)
colnames(countstable) <- c("cell", "GM_NP_Mock", "GM_NP_IFNb")
countstable <- merge(countstable, table3, by.x = "cell", by.y = "tabletemp_orig", all = TRUE)
colnames(countstable)[4] <- "GM_NP_Mock2zika"
countstable <- merge(countstable, table4, by.x = "cell", by.y = "tabletemp_orig", all = TRUE)
colnames(countstable)[5] <- "GM_NP_ZIKV"
countstable <- merge(countstable, table5, by.x = "cell", by.y = "tabletemp_orig", all = TRUE)
colnames(countstable)[6] <- "DD_NP_Mock"
countstable <- merge(countstable, table6, by.x = "cell", by.y = "tabletemp_orig", all = TRUE)
colnames(countstable)[7] <- "DD_NP_IFNb"
countstable <- merge(countstable, table7, by.x = "cell", by.y = "tabletemp_orig", all = TRUE)
colnames(countstable)[8] <- "DD_NP_Mock2zika"
countstable <- merge(countstable, table8, by.x = "cell", by.y = "tabletemp_orig", all = TRUE)
colnames(countstable)[9] <- "DD_NP_ZIKV"
countstable[is.na(countstable)] <- 0
counts_melt <- melt(countstable)
countstablesample <- merge(table_gm, table_dd, by.x = "tabletemp_orig", by.y = "tabletemp_orig", all = TRUE)
colnames(countstablesample) <- c("cell", "GM", "DD")
countstablesample[is.na(countstablesample)] <- 0
countstablesample_melt <- melt(countstablesample)

# Generate output tables
write.csv(counts_melt, file.path(results_folder, "Count_table.csv"))
write.csv(countstablesample_melt, file.path(results_folder, "Countsample_table.csv"))
write.csv(total_zika_100, file.path(results_folder, "total_zika_100_table.csv"))

# Generate bar and violin plot for publication
ggplot(data = countstablesample_melt, aes(x = factor(cell, levels = c("NPC", "Early Neurons", "Neurons",
             "Astrocytes", "Oligodendrocyte")), y = value, 
             fill = factor(cell, levels = c("NPC", "Early Neurons", "Neurons", "Astrocytes", "Oligodendrocyte")))) +
    geom_bar(stat = "identity") +
    theme_Publication() +
    facet_wrap(~ factor(variable, levels = c("GM", "DD")), ncol = 1) +
    theme(legend.text = element_text(size = 8)) +
    labs(fill = "", x = "") +
    theme(strip.background = element_blank(), strip.text.x = element_blank()) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_color_manual(values = c("NPC" = "#CCBB44", "Early Neurons" = "#66CCEE", 
                "Neurons" = "#4477AA", "Astrocytes" = "#228833", "Oligodendrocyte" = "#AA3377")) +
    scale_fill_manual(values = c("NPC" = "#CCBB44", "Early Neurons" = "#66CCEE", 
                "Neurons" = "#4477AA", "Astrocytes" = "#228833", "Oligodendrocyte" = "#AA3377")) +
    xlab("") +
    # ylim(c(0, 1)) +
    ylab("% of Cells") +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 12)
    ) +
    theme(legend.position = "none") +
    theme(panel.spacing = unit(1, "lines"))
ggsave(file.path(results_folder, "Cell_percentages4pub.png"), height = 6, width = 3, dpi = 500)
ggsave(file.path(results_folder, "Cell_percentages4pub.svg"), height = 6, width = 3, dpi = 500)
ggsave(file.path(results_folder, "Cell_percentages4pub.pdf"), height = 6, width = 3, dpi = 500)



ggplot(data = total_zika_100, aes(x = factor(cell, levels = c(
    "NPC", "Early Neurons", "Neurons", "Astrocytes", "Oligodendrocyte"
)), y = Freq, fill = z)) +
    geom_bar(stat = "identity") +
    theme_Publication() +
    facet_wrap(~ factor(sample, level = c("GM_NP_ZIKV", "DD_NP_ZIKV"))) +
    theme(strip.background = element_rect(fill = c("white"))) +
    theme(legend.text = element_text(size = 8)) +
    scale_color_manual(values = c("#ddd2d2", "black")) +
    scale_fill_manual(values = c("#ddd2d2", "black")) +
    xlab("") +
    ylim(c(0, 1)) +
    ylab("% of Cells") +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 12)
    )
ggsave(file.path(results_folder, "CellCount_percentage_zika_1004pub.png"), height = 4, width = 5.5, dpi = 500)
ggsave(file.path(results_folder, "CellCount_percentage_zika_1004pub.svg"), height = 4, width = 5.5, dpi = 500)
ggsave(file.path(results_folder, "CellCount_percentage_zika_1004pub.pdf"), height = 4, width = 5.5, dpi = 500)

##### --- generate ISG counts ---#####
results_folder <- generate_folder(file.path(inpc_results, "ISG_results"))
DefaultAssay(inpc.combined_sd) <- "RNA"
num_cells_expressing_gene("IL6", inpc.combined_sd, results_folder)
num_cells_expressing_gene("IL10", inpc.combined_sd, results_folder)
num_cells_expressing_gene("TNF", inpc.combined_sd, results_folder) # TNFa
num_cells_expressing_gene("IFNA1", inpc.combined_sd, results_folder)
num_cells_expressing_gene("IFNA2", inpc.combined_sd, results_folder)
num_cells_expressing_gene("IFNB1", inpc.combined_sd, results_folder)

##### --- DE Analysis ---#####
results_folder <- generate_folder(file.path(inpc_results, "DE_results"))

#--make sure active assay is RNA
DefaultAssay(inpc.combined_sd) <- "RNA"

#-- Zika positive vs Zika negative cells with in an infected sample
results_folder <- generate_folder(file.path(inpc_results, "DE_results", "zika+vszika-"))

cells <- inpc.combined_sd@meta.data[
    inpc.combined_sd@meta.data == "GM_NP_ZIKV",
]
inpc_subset <- inpc.combined_sd[, rownames(cells)]
inpc_subset$cell_zika <- paste(inpc_subset@meta.data$garnett_cluster_extend_lw, inpc_subset@meta.data$ZIKA, sep = "_")
Idents(inpc_subset) <- inpc_subset$cell_zika

DEgeneszikagmNPC <- get_DE_between_conditions("NPC_zika", "NPC_nozika", "GM_ZIKA_NPC", inpc_subset, results_folder, fontsize = 12, height = 4.5)
DEgeneszikagmEarlyNeurons <- get_DE_between_conditions("Early Neurons_zika", "Early Neurons_nozika", "GM_ZIKA_Neurons", inpc_subset, results_folder, fontsize = 10)
DEgeneszikagmNeurons <- get_DE_between_conditions("Neurons_zika", "Neurons_nozika", "GM_ZIKA_Neurons", inpc_subset, results_folder, fontsize = 10)
DEgeneszikagmAstrocytes <- get_DE_between_conditions("Astrocytes_zika", "Astrocytes_nozika", "GM_ZIKA_Astrocytes", inpc_subset, results_folder, fontsize = 10, height = 4.5)
DEgeneszikagmoligo <- get_DE_between_conditions("Oligodendrocyte_zika", "Oligodendrocyte_nozika", "GM_ZIKA_Oligodendrocyte", inpc_subset, results_folder, fontsize = 12, height = 4.5)

## -DD Zika vs nozika
cells <- inpc.combined_sd@meta.data[
    inpc.combined_sd@meta.data == "DD_NP_ZIKV",
]
inpc_subset <- inpc.combined_sd[, rownames(cells)]
inpc_subset$cell_zika <- paste(inpc_subset@meta.data$garnett_cluster_extend_lw, inpc_subset@meta.data$ZIKA, sep = "_")
Idents(inpc_subset) <- inpc_subset$cell_zika

DEgeneszikaddNPC <- get_DE_between_conditions("NPC_zika", "NPC_nozika", "DD_ZIKA_NPC", inpc_subset, results_folder, fontsize = 12, height = 4.5)
DEgeneszikaddEarlyNeurons <- get_DE_between_conditions("Early Neurons_zika", "Early Neurons_nozika", "DD_ZIKA_Neurons", inpc_subset, results_folder, fontsize = 11, height = 6)
DEgeneszikaddNeurons <- get_DE_between_conditions("Neurons_zika", "Neurons_nozika", "DD_ZIKA_Neurons", inpc_subset, results_folder, fontsize = 11, height = 6)
DEgeneszikaddAstrocytes <- get_DE_between_conditions("Astrocytes_zika", "Astrocytes_nozika", "DD_ZIKA_Astrocytes", inpc_subset, results_folder, fontsize = 12, height = 5)
DEgeneszikaddoligo <- get_DE_between_conditions("Oligodendrocyte_zika", "Oligodendrocyte_nozika", "DD_ZIKA_Oligodendrocyte", inpc_subset, results_folder, fontsize = 12, height = 4.5)

#-- Zika positive vs Zika negative cells with in an infected sample
results_folder <- generate_folder(file.path(inpc_results, "DE_results", "IFNvsCRL"))

## -GM IFNb vs Mock
inpc.combined_sd$cell_condition <- paste(inpc.combined_sd@meta.data$garnett_cluster_extend_lw, inpc.combined_sd@meta.data$orig.ident, sep = "_")

Idents(inpc.combined_sd) <- inpc.combined_sd$cell_condition

DEgenesifnbgmNPC <- get_DE_between_conditions("NPC_GM_NP_IFNb", "NPC_GM_NP_Mock", "GM_IFNb_NPC", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)
DEgeneifnbagmEarlyNeurons <- get_DE_between_conditions("Early Neurons_GM_NP_IFNb", "Early Neurons_GM_NP_Mock", "GM_IFNb_Neurons", inpc.combined_sd, results_folder, fontsize = 10)
DEgeneifnbgmNeurons <- get_DE_between_conditions("Neurons_GM_NP_IFNb", "Neurons_GM_NP_Mock", "GM_IFNb_Neurons", inpc.combined_sd, results_folder, fontsize = 10)
DEgenesifnbgmAstrocytes <- get_DE_between_conditions("Astrocytes_GM_NP_IFNb", "Astrocytes_GM_NP_Mock", "GM_IFNb_Astrocytes", inpc.combined_sd, results_folder, fontsize = 10, height = 4.5)
# Oligo's have to few cells to run 
# DEgenesifnbagmoligo <- get_DE_between_conditions("Oligodendrocyte_GM_NP_IFNb", "Oligodendrocyte_GM_NP_Mock", "GM_IFNb_Oligodendrocyte", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)

## -DD IFNb vs Mock
Idents(inpc.combined_sd) <- inpc.combined_sd$cell_condition

DEgenesifnbddNPC <- get_DE_between_conditions("NPC_DD_NP_IFNb", "NPC_DD_NP_Mock", "DD_IFNb_NPC", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)
DEgeneifnbaddEarlyNeurons <- get_DE_between_conditions("Early Neurons_DD_NP_IFNb", "Early Neurons_DD_NP_Mock", "GM_IFNb_Neurons", inpc.combined_sd, results_folder, fontsize = 10)
DEgeneifnbddNeurons <- get_DE_between_conditions("Neurons_DD_NP_IFNb", "Neurons_DD_NP_Mock", "DD_IFNb_Neurons", inpc.combined_sd, results_folder, fontsize = 10)
DEgenesifnbddAstrocytes <- get_DE_between_conditions("Astrocytes_DD_NP_IFNb", "Astrocytes_DD_NP_Mock", "DD_IFNb_Astrocytes", inpc.combined_sd, results_folder, fontsize = 10, height = 4.5)
DEgenesifnbaddoligo <- get_DE_between_conditions("Oligodendrocyte_DD_NP_IFNb", "Oligodendrocyte_DD_NP_Mock", "DD_IFNb_Oligodendrocyte", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)

### --- Generate an upset plot 4 publication 
genelisttotal <- list(
    "DE Astros DD Zika" = rownames(DEgeneszikaddAstrocytes), "DE Neurons DD Zika" = rownames(DEgeneszikaddNeurons),
    "DE Early Neurons DD Zika" = rownames(DEgeneszikaddEarlyNeurons), "DE NPC DD Zika" = rownames(DEgenesifnbddNPC),
    "DE NPC GM Zika" = rownames(DEgeneszikagmNPC),
    "DE Astros DD IFN" = rownames(DEgenesifnbddAstrocytes), "DE Neurons DD IFN" = rownames(DEgeneifnbddNeurons),
    "DE Early Neurons DD IFN" = rownames(DEgeneifnbaddEarlyNeurons), "DE NPC DD IFN" = rownames(DEgenesifnbddNPC),
    "DE NPC GM IFN" = rownames(DEgenesifnbgmNPC)
)

svg(file.path(inpc_results, "DE_results", "Upsetgraph_DEgenesZika+VsZika-andIFN.svg"))
# png(file.path(inpc_results, "DE_results", "Upsetgraph_DEgenesZika+VsZika-andIFN.png"), res = 250, units = "in", width = 8, height = 6)
upset(fromList(genelisttotal), sets = rev(c(
    "DE NPC GM IFN", "DE NPC DD IFN", "DE Early Neurons DD IFN", "DE Neurons DD IFN",
    "DE Astros DD IFN", "DE NPC GM Zika", "DE NPC DD Zika", "DE Early Neurons DD Zika", "DE Neurons DD Zika", "DE Astros DD Zika"
)), keep.order = TRUE, order.by = "freq")
dev.off()
svg(file.path(inpc_results, "DE_results", "Upsetgraph_DEgenes_degreeZika+VsZika-andIFN.svg"))
# png(file.path(inpc_results, "DE_results", "Upsetgraph_DEgenes_degreeZika+VsZika-andIFN.png"), res = 250, units = "in", width = 8, height = 5)
upset(fromList(genelisttotal), sets = rev(c(
    "DE NPC GM IFN", "DE NPC DD IFN", "DE Early Neurons DD IFN", "DE Neurons DD IFN",
    "DE Astros DD IFN", "DE NPC GM Zika", "DE NPC DD Zika", "DE Early Neurons DD Zika", "DE Neurons DD Zika", "DE Astros DD Zika"
)), keep.order = TRUE, order.by = "degree")
dev.off()


#-- Zika positive vs Zika negative cells with in an infected sample
results_folder <- generate_folder(file.path(inpc_results, "DE_results", "ZikavsMock"))

# -GM ZIKA vs Mock
DEgenesZIKAgmNPC <- get_DE_between_conditions("NPC_GM_NP_ZIKV", "NPC_GM_NP_Mock2zika", "GM_ZIKA_NPC2mock", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)
DEgeneZIKAgmEarlyNeurons <- get_DE_between_conditions("Early Neurons_GM_NP_ZIKV", "Early Neurons_GM_NP_Mock2zika", "GM_ZIKA_Neurons2mock", inpc.combined_sd, results_folder, fontsize = 10)
DEgeneZIKAgmNeurons <- get_DE_between_conditions("Neurons_GM_NP_ZIKV", "Neurons_GM_NP_Mock2zika", "GM_ZIKA_Neurons2mock", inpc.combined_sd, results_folder, fontsize = 10)
DEgenesZIKAgmAstrocytes <- get_DE_between_conditions("Astrocytes_GM_NP_ZIKV", "Astrocytes_GM_NP_Mock2zika", "GM_ZIKA_Astrocyte2mocks", inpc.combined_sd, results_folder, fontsize = 10, height = 4.5)
DEgenesZIKAgmoligo <- get_DE_between_conditions("Oligodendrocyte_GM_NP_ZIKV", "Oligodendrocyte_GM_NP_Mock2zika", "GM_ZIKA_Oligodendrocyte2mock", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)


# -DD ZIKA vs Mock
DEgenesZIKAddNPC <- get_DE_between_conditions("NPC_DD_NP_ZIKV", "NPC_DD_NP_Mock2zika", "DD_ZIKA_NPC2mock", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)
DEgeneZIKAddEarlyNeurons <- get_DE_between_conditions("Early Neurons_DD_NP_ZIKV", "Early Neurons_DD_NP_Mock2zika", "GM_ZIKA_Neurons2mock", inpc.combined_sd, results_folder, fontsize = 10)
DEgeneZIKAddNeurons <- get_DE_between_conditions("Neurons_DD_NP_ZIKV", "Neurons_DD_NP_Mock2zika", "GM_ZIKA_Neurons2mock", inpc.combined_sd, results_folder, fontsize = 10)
DEgenesZIKAddAstrocytes <- get_DE_between_conditions("Astrocytes_DD_NP_ZIKV", "Astrocytes_DD_NP_Mock2zika", "GM_ZIKA_Astrocytes2mock", inpc.combined_sd, results_folder, fontsize = 10, height = 4.5)
DEgenesZIKAddoligo <- get_DE_between_conditions("Oligodendrocyte_DD_NP_ZIKV", "Oligodendrocyte_DD_NP_Mock2zika", "GM_ZIKA_Oligodendrocyte2mock", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)


## zika positive cells vs mock and Zika- cells vs mock
results_folder <- generate_folder(file.path(inpc_results, "DE_results", "Zika+&Zika-vsMock"))

inpc.combined_sd$cell_zika_condition <- paste(inpc.combined_sd@meta.data$garnett_cluster_extend_lw, 
    inpc.combined_sd@meta.data$ZIKA, inpc.combined_sd@meta.data$orig.ident, sep = "_")
Idents(inpc.combined_sd) <- inpc.combined_sd$cell_zika_condition

# -GM ZIKA vs Mock
DEgenesZIKAposgmNPC <- get_DE_between_conditions("NPC_zika_GM_NP_ZIKV", "NPC_nozika_GM_NP_Mock2zika", "GM_ZIKApos_NPC2mock", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)
DEgeneZIKAposgmEarlyNeurons <- get_DE_between_conditions("Early Neurons_zika_GM_NP_ZIKV", "Early Neurons_nozika_GM_NP_Mock2zika", "GM_ZIKApos_Neurons2mock", inpc.combined_sd, results_folder, fontsize = 10)
DEgeneZIKAposgmNeurons <- get_DE_between_conditions("Neurons_zika_GM_NP_ZIKV", "Neurons_nozika_GM_NP_Mock2zika", "GM_ZIKApos_Neurons2mock", inpc.combined_sd, results_folder, fontsize = 10)
DEgenesZIKAposgmAstrocytes <- get_DE_between_conditions("Astrocytes_zika_GM_NP_ZIKV", "Astrocytes_nozika_GM_NP_Mock2zika", "GM_ZIKApos_Astrocyte2mocks", inpc.combined_sd, results_folder, fontsize = 10, height = 4.5)
DEgenesZIKAposgmoligo <- get_DE_between_conditions("Oligodendrocyte_zika_GM_NP_ZIKV", "Oligodendrocyte_nozika_GM_NP_Mock2zika", "GM_ZIKApos_Oligodendrocyte2mock", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)

DEgenesZIKAneggmNPC <- get_DE_between_conditions("NPC_nozika_GM_NP_ZIKV", "NPC_nozika_GM_NP_Mock2zika", "GM_ZIKAneg_NPC2mock", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)
DEgeneZIKAneggmEarlyNeurons <- get_DE_between_conditions("Early Neurons_nozika_GM_NP_ZIKV", "Early Neurons_nozika_GM_NP_Mock2zika", "GM_ZIKAneg_Neurons2mock", inpc.combined_sd, results_folder, fontsize = 10)
DEgeneZIKAneggmNeurons <- get_DE_between_conditions("Neurons_nozika_GM_NP_ZIKV", "Neurons_nozika_GM_NP_Mock2zika", "GM_ZIKAneg_Neurons2mock", inpc.combined_sd, results_folder, fontsize = 10)
DEgenesZIKAneggmAstrocytes <- get_DE_between_conditions("Astrocytes_nozika_GM_NP_ZIKV", "Astrocytes_nozika_GM_NP_Mock2zika", "GM_ZIKAneg_Astrocyte2mocks", inpc.combined_sd, results_folder, fontsize = 10, height = 4.5)
DEgenesZIKAneggmoligo <- get_DE_between_conditions("Oligodendrocyte_nozika_GM_NP_ZIKV", "Oligodendrocyte_nozika_GM_NP_Mock2zika", "GM_ZIKAneg_Oligodendrocyte2mock", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)

# --DD ZIKApos vs Mock
DEgenesZIKAposddNPC <- get_DE_between_conditions("NPC_zika_DD_NP_ZIKV", "NPC_nozika_DD_NP_Mock2zika", "DD_ZIKApos_NPC2mock", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)
DEgeneZIKAposddEarlyNeurons <- get_DE_between_conditions("Early Neurons_zika_DD_NP_ZIKV", "Early Neurons_nozika_DD_NP_Mock2zika", "GM_ZIKApos_Neurons2mock", inpc.combined_sd, results_folder, fontsize = 10)
DEgeneZIKAposddNeurons <- get_DE_between_conditions("Neurons_zika_DD_NP_ZIKV", "Neurons_nozika_DD_NP_Mock2zika", "DD_ZIKApos_Neurons2mock", inpc.combined_sd, results_folder, fontsize = 10)
DEgenesZIKAposddAstrocytes <- get_DE_between_conditions("Astrocytes_zika_DD_NP_ZIKV", "Astrocytes_nozika_DD_NP_Mock2zika", "DD_ZIKApos_Astrocytes2mock", inpc.combined_sd, results_folder, fontsize = 10, height = 4.5)
DEgenesZIKAposdoligo <- get_DE_between_conditions("Oligodendrocyte_zika_DD_NP_ZIKV", "Oligodendrocyte_nozika_DD_NP_Mock2zika", "DD_ZIKApos_Oligodendrocyte2mock", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)

DEgenesZIKAnegddNPC <- get_DE_between_conditions("NPC_nozika_DD_NP_ZIKV", "NPC_nozika_DD_NP_Mock2zika", "DD_ZIKAneg_NPC2mock", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)
DEgeneZIKAnegddEarlyNeurons <- get_DE_between_conditions("Early Neurons_nozika_DD_NP_ZIKV", "Early Neurons_nozika_DD_NP_Mock2zika", "DD_ZIKAneg_Neurons2mock", inpc.combined_sd, results_folder, fontsize = 10)
DEgeneZIKAnegddNeurons <- get_DE_between_conditions("Neurons_nozika_DD_NP_ZIKV", "Neurons_nozika_DD_NP_Mock2zika", "DD_ZIKAneg_Neurons2mock", inpc.combined_sd, results_folder, fontsize = 10)
DEgenesZIKAnegddAstrocytes <- get_DE_between_conditions("Astrocytes_nozika_DD_NP_ZIKV", "Astrocytes_nozika_DD_NP_Mock2zika", "DD_ZIKAneg_Astrocytes2mock", inpc.combined_sd, results_folder, fontsize = 10, height = 4.5)
DEgenesZIKAnegddoligo <- get_DE_between_conditions("Oligodendrocyte_nozika_DD_NP_ZIKV", "Oligodendrocyte_nozika_DD_NP_Mock2zika", "DD_ZIKAneg_Oligodendrocyte2mock", inpc.combined_sd, results_folder, fontsize = 12, height = 4.5)

##### --- Scaterplot of Zika + & Zika - vs Mock ---#####
results_folder <- generate_folder(file.path(inpc_results, "DE_results", "scatterplots"))

#### --NPCs-- ####

# genes sig in both
both <- intersect(rownames(DEgenesZIKAposgmNPC), rownames(DEgenesZIKAneggmNPC))
# genes sig in only zika postive cells
pos <-  setdiff(rownames(DEgenesZIKAposgmNPC), rownames(DEgenesZIKAneggmNPC))
# genes sig in only zika negative cells 
neg <- setdiff(rownames(DEgenesZIKAneggmNPC), rownames(DEgenesZIKAposgmNPC))

allgenes <- union(rownames(DEgenesZIKAposgmNPC), rownames(DEgenesZIKAneggmNPC))

fullNPC_Zikapos <- read.table(file.path(
    inpc_results, "DE_results", "Zika+&Zika-vsMock",
    "DEgenes_full_NPC_zika_GM_NP_ZIKV_NPC_nozika_GM_NP_Mock2zika_GM_ZIKApos_NPC2mock.txt"
), row.names = 1)

fullNPC_Zikaneg <- read.table(file.path(
    inpc_results, "DE_results", "Zika+&Zika-vsMock",
    "DEgenes_full_NPC_nozika_GM_NP_ZIKV_NPC_nozika_GM_NP_Mock2zika_GM_ZIKAneg_NPC2mock.txt"
), row.names = 1)


tab <- merge(fullNPC_Zikapos[allgenes, ], fullNPC_Zikaneg[allgenes,], by = "row.names")
lmfitdata <- lm(avg_logFC.y ~ avg_logFC.x, data = tab)
lmfitdata.res <- resid(lmfitdata)
tab$res <- lmfitdata.res
tab$res_abs <- abs(lmfitdata.res)
genes_median <- c()

for (i in tab$Row.names) {
    if (tab[tab$Row.names == i, "res_abs"] > median(tab$res_abs)) {
        genes_median <- c(genes_median, i)
    } else {
        genes_median <- c(genes_median, "")
    }
}
tab$geneswithresgtmedian <- genes_median

significance <- c()
for (g in tab$Row.names) {
    if (g %in% both) {
        significance <- c(significance, "both")
    } else if (g %in% pos) {
        significance <- c(significance, "Sig. Zika+ vs Mock")
    } else if (g %in% neg) {
        significance <- c(significance, "Sig. Zika- vs Mock")
    }
}
tab$significance <- significance
write.csv(tab, file.path(results_folder, "NPC_table.csv"))

ggplot(tab, aes(x = avg_logFC.x, y = avg_logFC.y)) + # , fill=significance)) +
    geom_point(aes(colour = significance)) +
    labs(x = "Zika+ vs Mock LFCs", y = "Zika- vs Mock LFCs") +
    geom_smooth(method = lm) +
    scale_color_manual(name = "significance", values = c("both" = "black",
     "Sig. Zika+ vs Mock" = "purple", "Sig. Zika- vs Mock" = "#861212")) +
    # geom_text(
    #     data = subset(tab, res_abs > median(tab$res_abs)),
    #     aes(avg_log2FC.x, avg_log2FC.y, label = Row.names), nudge_x = 0.25, size = 2
    # ) +
      geom_text_repel(
        data = subset(tab, res_abs >  median(tab$res_abs)),# quantile(tab$res_abs, 0.95)),
        mapping = aes(avg_logFC.x, avg_logFC.y, label = geneswithresgtmedian),
        size = 3, size = 3, box.padding = unit(0.05, "lines")
      ) +
    #  scale_y_continuous(limits = c(-2, 3), breaks = seq(-2, 3, by = 0.5)) +
    #  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 0.5)) +
    theme_Publication()
ggsave(file.path(results_folder, "GM_NPC_Zika+&Zika+vsMock.png"), units = "in", width = 6, height = 4, dpi = 300)
ggsave(file.path(results_folder, "GM_NPC_Zika+&Zika+vsMock.svg"), units = "in", width = 6, height = 4, dpi = 300)
ggsave(file.path(results_folder, "GM_NPC_Zika+&Zika+vsMock.pdf"), units = "in", width = 6, height = 4, dpi = 300)


#### --Astrocytes -- ####

# genes sig in both
both <- intersect(rownames(DEgenesZIKAposddAstrocytes), rownames(DEgenesZIKAnegddAstrocytes))
# genes sig in only zika postive cells
pos <- setdiff(rownames(DEgenesZIKAposddAstrocytes), rownames(DEgenesZIKAnegddAstrocytes))
# genes sig in only zika negative cells
neg <- setdiff(rownames(DEgenesZIKAnegddAstrocytes), rownames(DEgenesZIKAposddAstrocytes))

allgenes <- union(rownames(DEgenesZIKAposddAstrocytes), rownames(DEgenesZIKAnegddAstrocytes))

fullNPC_Zikapos <- read.table(file.path(
    inpc_results, "DE_results", "Zika+&Zika-vsMock",
    "DEgenes_full_Astrocytes_zika_DD_NP_ZIKV_Astrocytes_nozika_DD_NP_Mock2zika_DD_ZIKApos_Astrocytes2mock.txt"
), row.names = 1)

fullNPC_Zikaneg <- read.table(file.path(
    inpc_results, "DE_results", "Zika+&Zika-vsMock",
    "DEgenes_full_Astrocytes_nozika_DD_NP_ZIKV_Astrocytes_nozika_DD_NP_Mock2zika_DD_ZIKAneg_Astrocytes2mock.txt"
), row.names = 1)

tab <- merge(fullNPC_Zikapos[allgenes, ], fullNPC_Zikaneg[allgenes, ], by = "row.names")
lmfitdata <- lm(avg_logFC.y ~ avg_logFC.x, data = tab)
lmfitdata.res <- resid(lmfitdata)
tab$res <- lmfitdata.res
tab$res_abs <- abs(lmfitdata.res)
genes_median <- c()

for (i in tab$Row.names) {
    if (tab[tab$Row.names == i, "res_abs"] > median(tab$res_abs)) {
        genes_median <- c(genes_median, i)
    } else {
        genes_median <- c(genes_median, "")
    }
}
tab$geneswithresgtmedian <- genes_median

significance <- c()
for (g in tab$Row.names) {
    if (g %in% both) {
        significance <- c(significance, "both")
    } else if (g %in% pos) {
        significance <- c(significance, "Sig. Zika+ vs Mock")
    } else if (g %in% neg) {
        significance <- c(significance, "Sig. Zika- vs Mock")
    }
}
tab$significance <- significance
write.csv(tab, file.path(results_folder, "Astrocytes_table.csv"))


ggplot(tab, aes(x = avg_logFC.x, y = avg_logFC.y)) + # , fill=significance)) +
    geom_point(aes(colour = significance)) +
    labs(x = "Zika+ vs Mock LFCs", y = "Zika- vs Mock LFCs") +
    geom_smooth(method = lm) +
    scale_color_manual(name = "significance", values = c(
        "both" = "black",
        "Sig. Zika+ vs Mock" = "purple", "Sig. Zika- vs Mock" = "#861212"
    )) +
    # geom_text(
    #     data = subset(tab, res_abs > median(tab$res_abs)),
    #     aes(avg_log2FC.x, avg_log2FC.y, label = Row.names), nudge_x = 0.25, size = 2
    # ) +
    geom_text_repel(
        data = subset(tab, res_abs > median(tab$res_abs)), # quantile(tab$res_abs, 0.95)),
        mapping = aes(avg_logFC.x, avg_logFC.y, label = geneswithresgtmedian),
        size = 3, size = 3, box.padding = unit(0.05, "lines")
    ) +
    #  scale_y_continuous(limits = c(-2, 3), breaks = seq(-2, 3, by = 0.5)) +
    #  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 0.5)) +
    theme_Publication()
ggsave(file.path(results_folder, "DD_Astro_Zika+&Zika+vsMock.png"), units = "in", width = 6, height = 4, dpi = 300)

#### --Neurons -- ####

# genes sig in both
both <- intersect(rownames(DEgeneZIKAposddNeurons), rownames(DEgeneZIKAnegddNeurons))
# genes sig in only zika postive cells
pos <- setdiff(rownames(DEgeneZIKAposddNeurons), rownames(DEgeneZIKAnegddNeurons))
# genes sig in only zika negative cells
neg <- setdiff(rownames(DEgeneZIKAnegddNeurons), rownames(DEgeneZIKAposddNeurons))

allgenes <- union(rownames(DEgeneZIKAposddNeurons), rownames(DEgeneZIKAnegddNeurons))

fullNPC_Zikapos <- read.table(file.path(
    inpc_results, "DE_results", "Zika+&Zika-vsMock",
    "DEgenes_full_Neurons_zika_DD_NP_ZIKV_Neurons_nozika_DD_NP_Mock2zika_DD_ZIKApos_Neurons2mock.txt"
), row.names = 1)

fullNPC_Zikaneg <- read.table(file.path(
    inpc_results, "DE_results", "Zika+&Zika-vsMock",
    "DEgenes_full_Neurons_nozika_DD_NP_ZIKV_Neurons_nozika_DD_NP_Mock2zika_DD_ZIKAneg_Neurons2mock.txt"
), row.names = 1)

significance <- c()
for (g in allgenes) {
    if (g %in% both) {
        significance <- c(significance, "both")
    } else if (g %in% pos) {
        significance <- c(significance, "Sig. Zika+ vs Mock")
    } else if (g %in% neg) {
        significance <- c(significance, "Sig. Zika- vs Mock")
    }
}
tab <- merge(fullNPC_Zikapos[allgenes, ], fullNPC_Zikaneg[allgenes, ], by = "row.names")
lmfitdata <- lm(avg_logFC.y ~ avg_logFC.x, data = tab)
lmfitdata.res <- resid(lmfitdata)
tab$res <- lmfitdata.res
tab$res_abs <- abs(lmfitdata.res)
genes_median <- c()

for (i in tab$Row.names) {
    if (tab[tab$Row.names == i, "res_abs"] > median(tab$res_abs)) {
        genes_median <- c(genes_median, i)
    } else {
        genes_median <- c(genes_median, "")
    }
}
tab$geneswithresgtmedian <- genes_median

significance <- c()
for (g in tab$Row.names) {
    if (g %in% both) {
        significance <- c(significance, "both")
    } else if (g %in% pos) {
        significance <- c(significance, "Sig. Zika+ vs Mock")
    } else if (g %in% neg) {
        significance <- c(significance, "Sig. Zika- vs Mock")
    }
}
tab$significance <- significance
write.csv(tab, file.path(results_folder, "Neurons_table.csv"))


ggplot(tab, aes(x = avg_logFC.x, y = avg_logFC.y)) + # , fill=significance)) +
    geom_point(aes(colour = significance)) +
    labs(x = "Zika+ vs Mock LFCs", y = "Zika- vs Mock LFCs") +
    geom_smooth(method = lm) +
    scale_color_manual(name = "significance", values = c(
        "both" = "black",
        "Sig. Zika+ vs Mock" = "purple", "Sig. Zika- vs Mock" = "#861212"
    )) +
    # geom_text(
    #     data = subset(tab, res_abs > median(tab$res_abs)),
    #     aes(avg_log2FC.x, avg_log2FC.y, label = Row.names), nudge_x = 0.25, size = 2
    # ) +
    geom_text_repel(
        data = subset(tab, res_abs > median(tab$res_abs)), # quantile(tab$res_abs, 0.95)),
        mapping = aes(avg_logFC.x, avg_logFC.y, label = geneswithresgtmedian),
        size = 3, size = 3, box.padding = unit(0.05, "lines")
    ) +
    #  scale_y_continuous(limits = c(-2, 3), breaks = seq(-2, 3, by = 0.5)) +
    #  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 0.5)) +
    theme_Publication()
ggsave(file.path(results_folder, "DD_Neurons_Zika+&Zika+vsMock.png"), units = "in", width = 6, height = 4, dpi = 300)

#### --Early Neurons -- ####

# genes sig in both
both <- intersect(rownames(DEgeneZIKAposddEarlyNeurons), rownames(DEgeneZIKAnegddEarlyNeurons))
# genes sig in only zika postive cells
pos <- setdiff(rownames(DEgeneZIKAposddEarlyNeurons), rownames(DEgeneZIKAnegddEarlyNeurons))
# genes sig in only zika negative cells
neg <- setdiff(rownames(DEgeneZIKAnegddEarlyNeurons), rownames(DEgeneZIKAposddEarlyNeurons))

allgenes <- union(rownames(DEgeneZIKAposddEarlyNeurons), rownames(DEgeneZIKAnegddEarlyNeurons))

fullNPC_Zikapos <- read.table(file.path(
    inpc_results, "DE_results", "Zika+&Zika-vsMock",
    "DEgenes_full_Early Neurons_zika_DD_NP_ZIKV_Early Neurons_nozika_DD_NP_Mock2zika_GM_ZIKApos_Neurons2mock.txt"
), row.names = 1)

fullNPC_Zikaneg <- read.table(file.path(
    inpc_results, "DE_results", "Zika+&Zika-vsMock",
    "DEgenes_full_Early Neurons_nozika_DD_NP_ZIKV_Early Neurons_nozika_DD_NP_Mock2zika_DD_ZIKAneg_Neurons2mock.txt"
), row.names = 1)

significance <- c()
for (g in allgenes) {
    if (g %in% both) {
        significance <- c(significance, "both")
    } else if (g %in% pos) {
        significance <- c(significance, "Sig. Zika+ vs Mock")
    } else if (g %in% neg) {
        significance <- c(significance, "Sig. Zika- vs Mock")
    }
}
tab <- merge(fullNPC_Zikapos[allgenes, ], fullNPC_Zikaneg[allgenes, ], by = "row.names")
lmfitdata <- lm(avg_logFC.y ~ avg_logFC.x, data = tab)
lmfitdata.res <- resid(lmfitdata)
tab$res <- lmfitdata.res
tab$res_abs <- abs(lmfitdata.res)
genes_median <- c()

for (i in tab$Row.names) {
    if (tab[tab$Row.names == i, "res_abs"] > median(tab$res_abs)) {
        genes_median <- c(genes_median, i)
    } else {
        genes_median <- c(genes_median, "")
    }
}
tab$geneswithresgtmedian <- genes_median

significance <- c()
for (g in tab$Row.names) {
    if (g %in% both) {
        significance <- c(significance, "both")
    } else if (g %in% pos) {
        significance <- c(significance, "Sig. Zika+ vs Mock")
    } else if (g %in% neg) {
        significance <- c(significance, "Sig. Zika- vs Mock")
    }
}
tab$significance <- significance
write.csv(tab, file.path(results_folder, "EarlyNeurons_table.csv"))


ggplot(tab, aes(x = avg_logFC.x, y = avg_logFC.y)) + # , fill=significance)) +
    geom_point(aes(colour = significance)) +
    labs(x = "Zika+ vs Mock LFCs", y = "Zika- vs Mock LFCs") +
    geom_smooth(method = lm) +
    scale_color_manual(name = "significance", values = c(
        "both" = "black",
        "Sig. Zika+ vs Mock" = "purple", "Sig. Zika- vs Mock" = "#861212"
    )) +
    # geom_text(
    #     data = subset(tab, res_abs > median(tab$res_abs)),
    #     aes(avg_log2FC.x, avg_log2FC.y, label = Row.names), nudge_x = 0.25, size = 2
    # ) +
    geom_text_repel(
        data = subset(tab, res_abs > median(tab$res_abs)), # quantile(tab$res_abs, 0.95)),
        mapping = aes(avg_logFC.x, avg_logFC.y, label = geneswithresgtmedian),
        size = 3, size = 3, box.padding = unit(0.05, "lines")
    ) +
    #  scale_y_continuous(limits = c(-2, 3), breaks = seq(-2, 3, by = 0.5)) +
    #  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 0.5)) +
    theme_Publication()
ggsave(file.path(results_folder, "DD_EarlyNeurons_Zika+&Zika+vsMock.png"), units = "in", width = 6, height = 4, dpi = 300)

#### --Oligo -- ####

# genes sig in both
both <- intersect(rownames(DEgenesZIKAposdoligo), rownames(DEgenesZIKAnegddoligo))
# genes sig in only zika postive cells
pos <- setdiff(rownames(DEgenesZIKAposdoligo), rownames(DEgenesZIKAnegddoligo))
# genes sig in only zika negative cells
neg <- setdiff(rownames(DEgenesZIKAnegddoligo), rownames(DEgenesZIKAposdoligo))

allgenes <- union(rownames(DEgenesZIKAposdoligo), rownames(DEgenesZIKAnegddoligo))
allgenes <- allgenes[!(allgenes %in% c("zika-antisense", "zika-positive"))]

fullNPC_Zikapos <- read.table(file.path(
    inpc_results, "DE_results", "Zika+&Zika-vsMock",
    "DEgenes_full_Oligodendrocyte_zika_DD_NP_ZIKV_Oligodendrocyte_nozika_DD_NP_Mock2zika_DD_ZIKApos_Oligodendrocyte2mock.txt"
), row.names = 1)

fullNPC_Zikaneg <- read.table(file.path(
    inpc_results, "DE_results", "Zika+&Zika-vsMock",
    "DEgenes_full_Oligodendrocyte_nozika_DD_NP_ZIKV_Oligodendrocyte_nozika_DD_NP_Mock2zika_DD_ZIKAneg_Oligodendrocyte2mock.txt"
), row.names = 1)


tab <- merge(fullNPC_Zikapos[allgenes, ], fullNPC_Zikaneg[allgenes, ], by = "row.names")
tab<- tab[!is.na(tab), ]
lmfitdata <- lm(avg_logFC.y ~ avg_logFC.x, data = tab)
lmfitdata.res <- resid(lmfitdata)
tab$res <- lmfitdata.res
tab$res_abs <- abs(lmfitdata.res)
genes_median <- c()

for (i in tab$Row.names) {
    if (tab[tab$Row.names == i, "res_abs"] > median(tab$res_abs)) {
        genes_median <- c(genes_median, i)
    } else {
        genes_median <- c(genes_median, "")
    }
}
tab$geneswithresgtmedian <- genes_median

significance <- c()
for (g in tab$Row.names) {
    if (g %in% both) {
        significance <- c(significance, "both")
    } else if (g %in% pos) {
        significance <- c(significance, "Sig. Zika+ vs Mock")
    } else if (g %in% neg) {
        significance <- c(significance, "Sig. Zika- vs Mock")
    }
}
tab$significance <- significance
write.csv(tab, file.path(results_folder, "Oligo_table.csv"))


ggplot(tab, aes(x = avg_logFC.x, y = avg_logFC.y)) + # , fill=significance)) +
    geom_point(aes(colour = significance)) +
    labs(x = "Zika+ vs Mock LFCs", y = "Zika- vs Mock LFCs") +
    geom_smooth(method = lm) +
    scale_color_manual(name = "significance", values = c(
        "both" = "black",
        "Sig. Zika+ vs Mock" = "purple", "Sig. Zika- vs Mock" = "#861212"
    )) +
    # geom_text(
    #     data = subset(tab, res_abs > median(tab$res_abs)),
    #     aes(avg_log2FC.x, avg_log2FC.y, label = Row.names), nudge_x = 0.25, size = 2
    # ) +
    geom_text_repel(
        data = subset(tab, res_abs > median(tab$res_abs)), # quantile(tab$res_abs, 0.95)),
        mapping = aes(avg_logFC.x, avg_logFC.y, label = geneswithresgtmedian),
        size = 3, size = 3, box.padding = unit(0.05, "lines")
    ) +
    #  scale_y_continuous(limits = c(-2, 3), breaks = seq(-2, 3, by = 0.5)) +
    #  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 0.5)) +
    theme_Publication()
ggsave(file.path(results_folder, "DD_Oligo_Zika+&Zika+vsMock.png"), units = "in", width = 6, height = 4, dpi = 300)


##### ---Expression of genes of interest---#####
results_folder <- generate_folder(file.path(inpc_results, "DotPlots"))
genes <- c("DDX58", "IFIH1", "TLR3", "TLR4", "IRF3", "IRF5", "IRF7", "IFNAR1", "IFNAR2", "TLR7", "NFKB1") # IFH1=MDA5, DDX58=RIGI

DotPlot(inpc.combined_sd, features = genes, group.by = "cell_condition", assay = "RNA", scale = TRUE) + coord_flip() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    scale_colour_gradient2(midpoint = 0, mid = "gray", high = "red", low = "blue")
ggsave(file.path(results_folder, "PRRS_dotplot_scaled4pub.png"), width = 9, height = 6, units = "in", dpi = 500)
ggsave(file.path(results_folder, "PRRS_dotplot_scaled4pub.svg"), width = 9, height = 6, units = "in", dpi = 500)
ggsave(file.path(results_folder, "PRRS_dotplot_scaled4pub.pdf"), width = 9, height = 6, units = "in", dpi = 500)


ifngenes <- c("IFIT1", "IFITM1", "IFITM3", "OAS1", "MX1", "ISG15", "USP18")
DotPlot(inpc.combined, features = ifngenes, group.by = "celltype_sample", assay = "RNA", scale = TRUE) + coord_flip() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_colour_gradient2(midpoint = 0, mid = "gray", high = "red", low = "blue")
ggsave(file.path(extra_results_path, "IFIT_dotplot_scaled.png"), width = 9, height = 6, units = "in", dpi = 300)
ggsave(file.path(extra_results_path, "IFIT_dotplot_scaled.svg"), width = 9, height = 6, units = "in", dpi = 300)
ggsave(file.path(extra_results_path, "IFIT_dotplot_scaled.pdf"), width = 9, height = 6, units = "in", dpi = 300)


#####################
# Functions --HFB
#####################

generate_monocle_cds <- function(sample_path) {
    SAMPLE <- load_cellranger_data(sample_path, umi_cutoff = 200)

    SAMPLE.counts <- SAMPLE@assays@data@listData$counts
    gene_meta_data <- rowData(SAMPLE)
    cell_meta_data <- colData(SAMPLE)

    cds <- new_cell_data_set(SAMPLE.counts,
        cell_metadata = cell_meta_data,
        gene_metadata = gene_meta_data
    )
    return(cds)
}

run_normalizationhfb <- function(combined, SD = TRUE, SCTRANSFORM = FALSE) {
    message("STATUS: normalizing data")
    all.genes <- rownames(combined)
    if (SD == TRUE) {
        combined <- ScaleData(combined,
            vars.to.regress = NULL,
            features = all.genes, verbose = FALSE
        )
    } else if (SCTRANSFORM == TRUE) {
        combined <- SCTransform(combined,
            vars.to.regress = NULL,
            verbose = FALSE
        )
    }
    return(combined)
}

feature_reductionhfb <- function(combined, result_folder, ndims = 20) {
    message("STATUS: performing PCA and UMAP")

    combined <- RunPCA(combined, npcs = 60, verbose = FALSE)
    ElbowPlot(combined, ndims = 60)
    ggsave(file.path(result_folder, "ElbowPlot.png"))
    combined <- RunUMAP(combined, reduction = "pca", dims = 1:ndims)
    combined <- RunTSNE(combined, reduction = "pca", dims = 1:ndims)

    DimPlot(combined,
        reduction = "umap",
        group.by = "orig.ident",
        label = FALSE
    )
    ggsave(file.path(result_folder, "UMAP_samples.png"), height = 5, width = 6, dpi = 500)

    cells <- combined@meta.data[
        combined@meta.data$ZIKA == "zika",
    ]

    DimPlot(combined,
        reduction = "umap",
        group.by = "ZIKA",
        cells.highlight = rownames(cells),
        sizes.highlight = 0.04,
        label = FALSE
    ) + scale_color_manual(labels = c("nozika", "zika"), values = c("#ddd2d2", "black")) +
        labs(x = "UMAP 1", y = "UMAP 2", title = "")
    ggsave(file.path(result_folder, "UMAP_ZIKA.png"), height = 5, width = 6, dpi = 500)

    DimPlot(combined,
        reduction = "umap",
        split.by = "condition",
        cells.highlight = rownames(cells),
        sizes.highlight = 0.04,
        label = FALSE
    )
    ggsave(file.path(result_folder, "UMAP_ZIKA_condition.png"), height = 5, width = 8, dpi = 500)
    combined <- FindNeighbors(combined, dims = 1:ndims)
    combined <- FindClusters(combined, resolution = 0.2)
    DimPlot(combined, reduction = "umap")
    ggsave(file.path(result_folder, "UMAP_clusters.png"), dpi = 300)

    DimPlot(combined,
        reduction = "umap",
        split.by = "condition",
        label = FALSE
    )
    ggsave(file.path(result_folder, "UMAP_cluster_condition.png"), height = 5, width = 8, dpi = 500)


    for (cluster in unique(Idents(combined))) {
        cluster.markers <- FindMarkers(combined, ident.1 = cluster, min.pct = 0.25)
        x <- head(cluster.markers, n = 6)
        VlnPlot(combined, pt.size = 0.1, features = rownames(x))
        ggsave(file.path(result_folder, paste0(cluster, "_cluster_markers.png")), height = 6, width = 6, dpi = 300)
        write.csv(cluster.markers, file.path(result_folder, paste0(cluster, "_cluster_markers.csv")))
    }
    return(combined)
}

initial_SC_Seurathfb <- function(data_dir, sampleID, ndims = 20, results_path) {
    message("STATUS: processing sample ", sampleID)

    pbmc.data <- Read10X(data.dir = data_dir)

    pbmc <- CreateSeuratObject(
        counts = pbmc.data,
        project = sampleID,
        min.cells = 0,
        min.features = 200
    )

    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    ggsave(file.path(results_path, paste0(sampleID, "_VnPlotMt.png")), dpi = 300)

    VlnPlot(pbmc, features = c("zika-positive"))
    ggsave(file.path(results_path,  paste0(sampleID, "_brazilzika_VnPlotMt.png")), dpi = 300)

    pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 10)
    VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    ggsave(file.path(results_path, paste0(sampleID, "_VnPlotMt_Filtered.png")), dpi = 300)

    message("STATUS: Number of cells ", length(colnames(pbmc)))
    if (length(colnames(pbmc)) > 100) {
        z <- pbmc["zika-positive", ]
        x <- z@meta.data[z@meta.data$nCount_RNA > 0, ]
        x1 <- z@meta.data[z@meta.data$nFeature_RNA == 1, ]
        if (length(rownames(x1)) > 0) {
            pbmczika <- subset(pbmc, cells = rownames(x1))

            VlnPlot(pbmczika, features = c("zika-positive"))
            ggsave(file.path(results_path, paste0(sampleID, "_brazilzika_VnPlotMt_onlyzikacells.png")), height = 5, width = 5, dpi = 300)

            VlnPlot(pbmczika, slot = "counts", features = c("zika-positive"))
            ggsave(file.path(results_path, paste0(sampleID, "_brazilzika_VnPlotMt_onlyzikacells_counts.png")), height = 5, width = 5, dpi = 300)
        } else {
            print(paste0("STATUS: no zika in any of the cells for sample ", sampleID))
        }

        if (isTRUE(str_detect(sampleID, "mock_"))) {
            pbmc$condition <- rep("Control", length(colnames(pbmc)))
        } else {
            pbmc$condition <- rep("Zika", length(colnames(pbmc)))
        }
        ## CLASSIFY CELLS AS ZIKA OR NO ZIKA HAVE TO HAVE GREATER THAN ONE TRANSCRIPT MAPPING TO THE ZIKA GENOME
        zikadata <- c()
        for (cell in colnames(pbmc)) {
            if (cell %in% rownames(x)) {
                zikadata <- c(zikadata, "zika")
            } else {
                zikadata <- c(zikadata, "nozika")
            }
        }
        pbmc$ZIKA <- zikadata


        print(paste0("STATUS: Percentage of cells with any Zika transcript ", dim(x1)[1] / length(pbmc$orig.ident)))
        print(paste0("STATUS: Percentage of cells Zika transcript at least 10 read counts ", dim(x)[1] / length(pbmc$orig.ident)))
        write(paste0("STATUS FOR SAMPLE ", sampleID, ": Percentage of cells with any Zika transcript ", dim(x1)[1] / length(pbmc$orig.ident), " Cells with Zika ", dim(x1)[1], " Total cells ", length(pbmc$orig.ident)),
            file.path(results_path, "2.ZikaExpressingCellStats.txt"),
            append = TRUE
        )
        write(paste0("STATUS FOR SAMPLE ", sampleID, ": Percentage of cells with any Zika transcript at least 10 read counts ", dim(x)[1] / length(pbmc$orig.ident), " Cells with Zika ", dim(x)[1], " Total cells ", length(pbmc$orig.ident)),
            file.path(results_path, "2.ZikaExpressingCellStats.txt"),
            append = TRUE
        )

        # Generate Monocle OBject
        cds <- generate_monocle_cds(str_remove(data_dir, "outs/filtered_feature_bc_matrix"))
        cds <- cds[, colnames(pbmc)]

        pbmc <- NormalizeData(pbmc,
            normalization.method = "LogNormalize",
            scale.factor = 10000
        )

        # Find features (genes) that are highly variable from cell-to-cell
        pbmc <- FindVariableFeatures(pbmc,
            selection.method = "vst",
            nfeatures = 2000
        )

        # Identify the 10 most highly variable genes
        top10 <- head(VariableFeatures(pbmc), 10)

        # plot variable features with and without labels
        plot1 <- VariableFeaturePlot(pbmc)
        LabelPoints(plot = plot1, points = top10, repel = TRUE)
        ggsave(file.path(results_path, paste0(sampleID, "_HighlyVariableGenes.png")), dpi = 300)

        pbmc <- ScaleData(pbmc)
        pbmc <- RunPCA(pbmc, npcs = 60)
        ElbowPlot(pbmc, ndims = 60)
        ggsave(file.path(results_path, paste0(sampleID, "_ElbowPlot.png")))
        pbmc <- RunUMAP(pbmc, dims = 1:ndims)

        VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
        ggsave(file.path(results_path, paste0(sampleID, "_PCA_Loadings.png")), dpi = 300)

        DimPlot(pbmc, reduction = "pca")
        ggsave(file.path(results_path, paste0(sampleID, "_PCA.png")), dpi = 300)

        DimPlot(pbmc, reduction = "umap")
        ggsave(file.path(results_path, paste0(sampleID, "_UMAP.png")), dpi = 300)

        FeaturePlot(pbmc, reduction = "umap", features = c("zika-positive"))
        ggsave(file.path(results_path, paste0(sampleID, "_brazilzika.png")), dpi = 300)

        pbmc <- FindNeighbors(pbmc, dims = 1:ndims)
        pbmc <- FindClusters(pbmc, resolution = 0.4) # Default values ued
        DimPlot(pbmc, reduction = "umap")
        ggsave(file.path(results_path, paste0(sampleID, "_clusters.png")), dpi = 300)

        for (cluster in unique(Idents(pbmc))) {
            cluster.markers <- FindMarkers(pbmc, ident.1 = cluster, min.pct = 0.25)
            x <- head(cluster.markers, n = 6)
            VlnPlot(pbmc, pt.size = 0.1, features = rownames(x))
            ggsave(file.path(results_path, paste0(cluster, "_cluster_markers.png")), height = 6, width = 6, dpi = 300)
            write.csv(cluster.markers, file.path(results_path, paste0(cluster, "_cluster_markers.csv")))
        }


        return(list("seurat" = pbmc, "monocle" = cds))
    } else {
        message("STATUS: Not enough cells for sample ", sampleID)
        return(NULL)
    }
}

integrate_and_process_data <- function(data, results_folder, ndims = 20) {
    anchors <- FindIntegrationAnchors(object.list = data, dims = 1:ndims)
    combined <- IntegrateData(anchorset = anchors, dims = 1:ndims)
    saveRDS(combined, file.path(results_folder, "HFB_combinedobject.rds"))

    DefaultAssay(combined) <- "integrated"
    combined <- run_normalizationhfb(combined)
    saveRDS(combined, file.path(results_folder, "HFB_combinedobject_norm.rds"))
    DefaultAssay(combined) <- "integrated"
    combined <- feature_reductionhfb(combined, results_folder)
    saveRDS(combined, file.path(results_folder, "HFB_combinedobject_norm_fr.rds"))
    return(combined)
}

get_DE_between_conditions <- function(ident_1, ident_2, compare,
                                      immune.combined,
                                      result_folder, fontsize = 9, height = 8,
                                      foldchange = 0.26,
                                      pval = 0.05, percent_cells = NULL) {
    print("STATUS: getting DEs...")
    ### Parameters for FindMarkers
    #### test.use: mast (default is wilcox)
    #### min.pct: 0.1(default: 0.1) filter out genes (features) that are detected at less than 10 percent frequency in cells in ident_1 or ident_2
    #### logfc: 0 (default is .25) logfc must be higher than 0 (I set this at 0 because I filter this out at a later stage - this was primarily done to
    ####                          understand how the tool (Seurat) works
    #### min.cells.feature: 3 (default is 3) minimum number of cells expressing the feature in at least one of the two groups (similar to min.pct)
    #### min.cells.group: 3 (defualt is 3) minimum number of cells in the group
    #### max.cells.per.ident: Inf (default Inf-means no down sampling) Down sample each identity class (cluster of cells) to a max number
    #### min.dif.pct: Inf (default Inf) Only test genes that show minimum difference in the fraction of detection between the two identities (cluster)
    if (file.exists(file.path(result_folder, paste0("DEgenes_full_", ident_1, "_", ident_2, "_", compare, ".txt")))) {
        message("STATUS: DE has already been done filtering from full file ")
        DEgenes <- read.table(file.path(result_folder, paste0("DEgenes_full_", ident_1, "_", ident_2, "_", compare, ".txt")), header = TRUE, row.names = 1)
    } else {
        DEgenes <- FindMarkers(immune.combined,
            ident.1 = ident_1,
            ident.2 = ident_2,
            test.use = "MAST",
            logfc.threshold = 0
        )
        write.table(data.frame(DEgenes),
            file.path(
                result_folder,
                paste0(
                    "DEgenes_full_",
                    ident_1,
                    "_", ident_2, "_", compare,
                    ".txt"
                )
            ),
            sep = "\t", quote = FALSE
        )
    }

    DEgenes[["gene_name"]] <- rownames(DEgenes)

    DE_sig_final <- DEgenes %>%
        filter(avg_logFC >= foldchange | avg_logFC <= -foldchange) %>%
        dplyr::select(gene_name, p_val, avg_logFC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
    DE_sig_final <- DE_sig_final %>%
        filter(p_val_adj <= pval) %>%
        dplyr::select(gene_name, p_val, avg_logFC, pct.1, pct.2, p_val_adj)
    if (!is.null(percent_cells)) {
        DE_sig_final <- DE_sig_final %>%
            filter(pct.1 >= percent_cells | pct.2 >= percent_cells) %>%
            dplyr::select(gene_name, p_val, avg_logFC, pct.1, pct.2, p_val_adj)
    }
    rownames(DE_sig_final) <- DE_sig_final$gene_name

    DE_sig_final$gene_name <- NULL
    if (is.null(percent_cells)) {
        write.table(data.frame(DE_sig_final),
            file.path(
                result_folder,
                paste0(
                    "DEgenes_sig_",
                    ident_1, "_", ident_2, "_", compare,
                    ".txt"
                )
            ),
            sep = "\t", quote = FALSE
        )
    } else {
        write.table(data.frame(DE_sig_final),
            file.path(
                result_folder,
                paste0(
                    "DEgenes_sig_",
                    ident_1, "_", ident_2, "_", compare,
                    "_", percent_cells, ".txt"
                )
            ),
            sep = "\t", quote = FALSE
        )
    }
    if (length(rownames(DE_sig_final)) > 0) {
        DotPlot(immune.combined,
            idents = c(ident_1, ident_2), features = rownames(DE_sig_final),
            cols = c("blue", "#eb3306"), dot.scale = 3
        ) + coord_flip() +
            theme(axis.text.y = element_text(size = fontsize))
        if (is.null(percent_cells)) {
            ggsave(file.path(
                result_folder, paste0(
                    "DEgenes_sig_", ident_1, "_",
                    ident_2, "_", compare, ".png"
                )
            ), height = 8, width = 6, units = "in", dpi = 350)
        } else {
            ggsave(file.path(
                result_folder, paste0(
                    "DEgenes_sig_", ident_1, "_",
                    ident_2, "_", compare, "_", percent_cells, ".png"
                )
            ), height = 8, width = 6, units = "in", dpi = 350)
        }
    }
    if (length(rownames(DE_sig_final)) > 1) {
        message("STATUS: Making BarPlot")
        pl <- SC_DE_barplot(DE_sig_final, ident_1, ident_2, fontsize = fontsize, horizontal = TRUE)
        if (is.null(percent_cells)) {
            ggsave(file.path(
                result_folder, paste0(
                    "DEgenes_sig_", ident_1, "_",
                    ident_2, "_", compare, "barPlot.png"
                )
            ),
            width = height, height = 4, units = "in", dpi = 300
            )
        } else {
            ggsave(file.path(
                result_folder, paste0(
                    "DEgenes_sig_", ident_1, "_",
                    ident_2, "_", compare, "_", percent_cells, "barPlot.png"
                )
            ),
            width = height, height = 4, units = "in", dpi = 300
            )
        }
    }
    return(DE_sig_final)
}

super_extend_garnett_clusters <- function(combined, result_folder) {
    combined[["garnett_cluster_extend_lw"]] <- combined$garnett_cluster_extend
    for (cluster in unique(combined$seurat_clusters)) {
        celltype_count <- c()
        clust_list <- which(combined$seurat_clusters == cluster)
        tempdf <- combined@meta.data[clust_list, ]
        celltypes <- unique(tempdf$garnett_cluster_extend)
        for (celltype in celltypes) {
            if (celltype != "Unknown") {
                celltype_count[celltype] <- length(which(tempdf$garnett_cluster_extend == celltype))
            }
        }
        maxvalue <- max(celltype_count)
        celltypefinal <- names(which(celltype_count == maxvalue))
        if (length(celltypefinal) > 1) {
            print(celltypefinal[1])
            celltypefinal <- celltypefinal[1]
        }
        print(paste0("STATUS: max cell count in cluster ", cluster, " is ", celltypefinal))
        tempdfunknown <- combined@meta.data[combined@meta.data$garnett_cluster_extend == "Unknown", ]
        tempdfunknown <- tempdfunknown[tempdfunknown$seurat_clusters == cluster, ]
        print(paste0("STATUS number of unknowns in cluster ", cluster, " is ", dim(tempdfunknown)[1]))
        for (row in rownames(tempdfunknown)) {
            combined@meta.data[row, "garnett_cluster_extend_lw"] <- celltypefinal
        }
    }
    Idents(combined) <- combined$garnett_cluster_extend_lw
    DimPlot(combined, reduction = "umap")
    ggsave(file.path(result_folder, "2.garnett_clusters_lw.png"), width = 6.5, height = 4, dpi = 500)

    return(combined)
}

num_cells_expressing_genehfb <- function(gene, combined, results_folder) {
    z <- combined[gene, ]
    x1 <- z@meta.data[z@meta.data$nFeature_RNA == 1, ]
    t <- table(x1$garnett_cluster_extend_lw)
    # z <- immune.combined_sd_sub@assays$RNA["IFNB1", ]
    T <- x1 %>%
        group_by(orig.ident, garnett_cluster_extend_lw, ZIKA) %>%
        summarise(Count = n())
    T <- data.frame(T)
    write.csv(T, file.path(results_folder, paste0(gene, "_table.csv")))
    ggplot(data = T, aes(
        x = factor(garnett_cluster_extend_lw, levels = c("NPC", "Early Neurons", "Neurons",
             "Neurons excitatory", "Neurons inhibitory", "Astrocytes", "Oligodendrocyte")),
        y = Count, fill = factor(ZIKA, levels = c("nozika", "zika"))
    )) +
        geom_bar(stat = "identity") +
        theme_Publication() +
        facet_wrap(~ factor(orig.ident), ncol=2) +
        theme(strip.background = element_rect(fill = c("white"))) +
        theme(legend.text = element_text(size = 8)) +
        scale_color_manual(values = c("#ddd2d2", "black")) +
        scale_fill_manual(values = c("#ddd2d2", "black")) +
        labs(x = "", fill = "Zika", title=gene) +
        # ylim(c(0, 1)) +
        ylab("# of Cells") +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
            axis.text.y = element_text(size = 12)
        )
    ggsave(file.path(results_folder, paste0("2.CellCount_", gene, "_zika.png")), height = 5.5, width = 6, dpi = 500)

}
####################
# HFB analsyis
####################

hfb_results <- "hfb_results"
generate_folder(hfb_results)

# path to aligned files 
SAMPLE_PATH <- "/vol08/ngs/Globus/GaleLab-19_done/WhitmoreAnalysis/SCGodata/"
sample_files <- list.dirs(SAMPLE_PATH, full.names = FALSE, recursive = FALSE)

##### ---Run individual samples of hfb data--#####
results_folder <- generate_folder(file.path(hfb_results, "single_analysis_results"))

allSeuratobjects <- list()
allSeuratobjectscds <- list()
for (sample in sample_files) {
    if (isTRUE(str_detect(sample, "FetalBrain"))) {
        sample_name <- str_remove(sample, "G\\d+_ZIKV_09_SC_FetalBrain_")
        sample_name <- str_remove(sample_name, "_Lib1_GL19$")
        print(sample_name)
        sobject <- initial_SC_Seurathfb(file.path(SAMPLE_PATH, sample, "outs", "filtered_feature_bc_matrix"),
                sampleID = sample_name, ndims=30, results_folder
        )
        allSeuratobjects[[sample_name]]=sobject$seurat
        allSeuratobjectscds[[sample_name]] = sobject$monocle
    }
}
saveRDS(allSeuratobjects, file.path(hfb_results, "single_analysis_results", "2.ListOfSeuratObjects.rds"))
saveRDS(allSeuratobjectscds, file.path(hfb_results, "single_analysis_results", "2.ListOfSeuratObjectscds.rds"))

allSeuratobjects <- readRDS(file.path(hfb_results, "single_analysis_results", "2.ListOfSeuratObjects.rds"))
allSeuratobjectscds <- readRDS(file.path(hfb_results, "single_analysis_results", "2.ListOfSeuratObjectscds.rds"))


##### ---USE GARNETT TO CLASSIFY CELL TYPES---#####

results_folder <- generate_folder(file.path(hfb_results, "garnett_classification_hfb_results"))

message("STATUS: Starting cell classification")
big_cds <- combine_cds(allSeuratobjectscds)
saveRDS(big_cds, file.path(hfb_results, "garnett_classification_hfb_results", "big_cds_svz.rds"))
# big_cds <- readRDS(file.path(hfb_results, "garnett_classification_hfb_results", "big_cds.rds"))

# Align/Normalize CDS object
message("STATUS: Cell classification alignment")
big_cds <- preprocess_cds(big_cds, num_dim = 100, preprocess_method = 100)
big_cds <- align_cds(big_cds, alignment_group = "sample")

big_cds <- reduce_dimension(big_cds,
    preprocess_method = "PCA",
    reduction_method = "UMAP"
)
# Build training data set using markers
message("STATUS: Cell classification build training")
marker_file_path <- "hfbmarkers.txt"
hfb_classifier <- train_cell_classifier(
    cds = big_cds,
    marker_file = marker_file_path,
    db = org.Hs.eg.db,
    cds_gene_id_type = "ENSEMBL",
    num_unknown = 500,
    marker_file_gene_id_type = "SYMBOL",
    cores=8
)
saveRDS(hfb_classifier, file.path(hfb_results, "garnett_classification_hfb_results", "hfb_classifier.rds"))
hfb_classifier_fpbmc <- readRDS(file.path(hfb_results, "garnett_classification_hfb_results", "hfb_classifier.rds"))

# Classify cells
message("STATUS: Cell classification")
big_cds1 <- classify_cells(big_cds, hfb_classifier,
    db = org.Hs.eg.db,
    cds_gene_id_type = "ENSEMBL",
    cluster_extend = TRUE
)
saveRDS(big_cds1, file.path(hfb_results, "garnett_classification_hfb_results", "big_cds_clusterext.rds"))
# big_cds1 <- readRDS(file.path(hfb_results, "garnett_classification_hfb_results", "big_cds_clusterext.rds"))

# Generate Marker files
marker_check <- check_markers(big_cds1, marker_file_path,
    db = org.Hs.eg.db,
    cds_gene_id_type = "ENSEMBL",
    marker_file_gene_id_type = "SYMBOL")
plot_markers(marker_check)
ggsave(file.path(hfb_results, "garnett_classification_hfb_results", "Marker_Fig.png"), dpi = 500)

##### ---Integrate inpc data ---#####
results_folder <- generate_folder(file.path(hfb_results, "integrated_results"))

message("STATUS: Integrate data")
brain.combined <- integrate_and_process_data(allSeuratobjects, file.path(hfb_results, "integrated_results"))

##### ---Extend clusters from garnett hfb data ---#####
results_folder <- generate_folder(file.path(hfb_results, "classification_extend_results"))

# add cell type classifications to integrated seurat object
g_cluster_ext <- big_cds1$cluster_ext_type
brain.combined[["garnett_cluster_extend"]] <- g_cluster_ext
Idents(brain.combined) <- brain.combined$garnett_cluster_extend
DimPlot(brain.combined, reduction = "umap")
ggsave(file.path(results_folder, "garnett_clusters.png"), width = 6.5, height = 4, dpi = 500)

brain.combined <- super_extend_garnett_clusters(brain.combined, file.path(hfb_results, "classification_extend_results"))
saveRDS(brain.combined, file.path(results_folder, "HFB_combinedobject_norm_fr_celltype.rds"))
# brain.combined <- readRDS(file.path(hfb_results, "classification_extend_results", "HFB_combinedobject_norm_fr_celltype.rds"))

Idents(brain.combined) <- brain.combined$garnett_cluster_extend_lw
DimPlot(brain.combined, reduction = "umap") + scale_color_manual(values = c("NPC" = "#CCBB44", "Neurons" = "#4477AA", "Neurons inhibitory" = "#454546", "Neurons excitatory" = "#4889ca", "Astrocytes" = "#228833", "Oligodendrocyte" = "#AA3377")) +
    theme(legend.text = element_text(size = 8))
ggsave(file.path(results_folder, "garnett_clusters_recolroed_lw.png"), width = 6, height = 4, dpi = 500)

DimPlot(brain.combined, ncol = 2, split.by = "orig.ident", reduction = "umap") + scale_color_manual(values = c("NPC" = "#CCBB44", "Neurons" = "#4477AA", "Neurons inhibitory" = "#454546", "Neurons excitatory" = "#4889ca", "Astrocytes" = "#228833", "Oligodendrocyte" = "#AA3377")) +
    theme(legend.text = element_text(size = 8))
ggsave(file.path(results_folder, "garnett_clusters_recolroed_sample_lw.png"), width = 8, height = 8, dpi = 500)
ggsave(file.path(results_folder, "garnett_clusters_recolroed_sample_lw.svg"), width = 8, height = 8, dpi = 500)
ggsave(file.path(results_folder, "garnett_clusters_recolroed_sample_lw.pdf"), width = 8, height = 8, dpi = 500)

## percentage table 
Idents(brain.combined) <- brain.combined$orig.ident
mock_48hr <- WhichCells(object = brain.combined, idents = "mock_48hr")
IFNb_48hr <- WhichCells(object = brain.combined, idents = "IFNb_48hr")
BRZ_48hr <- WhichCells(object = brain.combined, idents = "BRZ_48hr")
FSS_48hr <- WhichCells(object = brain.combined, idents = "FSS_48hr")

table1 <- brain.combined$garnett_cluster_extend_lw[names(brain.combined$garnett_cluster_extend_lw) %in% mock_48hr]
table2 <- brain.combined$garnett_cluster_extend_lw[names(brain.combined$garnett_cluster_extend_lw) %in% IFNb_48hr]
table3 <- brain.combined$garnett_cluster_extend_lw[names(brain.combined$garnett_cluster_extend_lw) %in% BRZ_48hr]
table4 <- brain.combined$garnett_cluster_extend_lw[names(brain.combined$garnett_cluster_extend_lw) %in% FSS_48hr]

zika_BRZ <- get_zika_percentages(table3, "BRZ_48hr")
zika_FSS <- get_zika_percentages(table4, "FSS_48hr")
total_zika <- rbind(zika_BRZ, zika_FSS)
total_zika$Freq <- total_zika$Freq * 100

zika_BRZ <- get_zika_percentages_total(table3, "BRZ_48hr")
zika_FSS <- get_zika_percentages_total(table4, "FSS_48hr")
total_zika_100 <- rbind(zika_BRZ, zika_FSS)
total_zika_100$Freq <- total_zika_100$Freq * 100

table1 <- get_percentages(table1)
table2 <- get_percentages(table2)
table3 <- get_percentages(table3)
table4 <- get_percentages(table4)

countstable <- merge(table1, table2, by.x = "tabletemp_orig", by.y = "tabletemp_orig", all = TRUE)
colnames(countstable) <- c("cell", "Mock", "IFNb")
countstable <- merge(countstable, table3, by.x = "cell", by.y = "tabletemp_orig", all = TRUE)
colnames(countstable)[4] <- "BRZ"
countstable <- merge(countstable, table4, by.x = "cell", by.y = "tabletemp_orig", all = TRUE)
colnames(countstable)[5] <- "FSS"
countstable[is.na(countstable)] <- 0
countstable$cell <- as.character(countstable$cell)
# countstable$cell[7] <- "Erythroid\nprogenitor\ncells"
counts_melt <- melt(countstable)
write.csv(counts_melt, file.path(cellclassification_results, "Count_table.csv"))

counts_melt$value <- counts_melt$value * 100
counts_melt$variable <- as.character(counts_melt$variable)
# counts_melt$variable <- factor(counts_melt$variable, levels = c("Control", "Zika"))

p <- ggplot(data = counts_melt, aes(x = factor(cell, levels = c("NPC", "Early Neurons", "Neurons",
             "Neurons excitatory", "Neurons inhibitory", "Astrocytes", "Oligodendrocyte")), 
             y = value, fill = cell)) +
    geom_bar(stat = "identity", position = position_dodge(width = .9)) +
    facet_wrap(~ factor(variable), ncol=4) +
    theme_classic() +
    theme(legend.text = element_text(size = 6)) +
    scale_fill_manual(values = c("NPC" = "#CCBB44", "Neurons" = "#4477AA", "Neurons inhibitory" = "#8f9092", 
    "Neurons excitatory" = "#4477AA", "Astrocytes" = "#228833", "Oligodendrocyte" = "#AA3377")) +
    labs(fill = "") +
    xlab("") +
    ylim(c(0, 70)) +
    ylab("% of cells") +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "top",
        legend.direction = "horizontal"
    )
ggsave(file.path(results_folder, "CellCount_percentage.png"), height = 4.5, width = 8, dpi = 500)
ggsave(file.path(results_folder, "CellCount_percentage.svg"), height = 4.5, width = 8, dpi = 500)
ggsave(file.path(results_folder, "CellCount_percentage.pdf"), height = 4.5, width = 8, dpi = 500)

##### --- generate ISG counts ---#####
results_folder <- generate_folder(file.path(hfb_results, "ISG_results"))
DefaultAssay(brain.combined_sd) <- "RNA"
num_cells_expressing_genehfb("IL6", brain.combined, results_folder)
num_cells_expressing_genehfb("IL10", brain.combined, results_folder)
num_cells_expressing_genehfb("TNF", brain.combined, results_folder) # TNFa
num_cells_expressing_genehfb("IFNB1", brain.combined, results_folder)
# - no counts for these two genes
# num_cells_expressing_genehfb("IFNA1", brain.combined, results_folder)
# num_cells_expressing_genehfb("IFNA2", brain.combined, results_folder)

##### ---DE Analysis ---#####
results_folder <- generate_folder(file.path(hfb_results, "DE_results"))

DefaultAssay(brain.combined) <- "RNA"
Idents(brain.combined) <- brain.combined$orig.ident

results_folder <- generate_folder(file.path(hfb_results, "DE_results", "ZikaVsMck"))

DEIFNb <- get_DE_between_conditions("IFNb_48hr", "mock_48hr", "IFNb", brain.combined, results_folder, fontsize = 8, height = 10)
DEBRZ <- get_DE_between_conditions("BRZ_48hr", "mock_48hr", "BRZ", brain.combined, results_folder, fontsize = 8, height = 10)
DEFSS <- get_DE_between_conditions("FSS_48hr", "mock_48hr", "FSS", brain.combined, results_folder, fontsize = 7, height = 10)

commongenes <- intersect(intersect(rownames(DEIFNb), rownames(DEBRZ)), rownames(DEFSS))
commonviral <- setdiff(intersect(rownames(DEFSS), rownames(DEBRZ)), rownames(DEIFNb))
brz_uniq <- setdiff(rownames(DEBRZ), unique(c(commonviral, commongenes)))
brz_fss <- setdiff(rownames(DEFSS), unique(c(commonviral, commongenes)))
ifnb_uniq <- setdiff(rownames(DEIFNb), c(rownames(DEBRZ), rownames(DEFSS)))


DoHeatmap(brain.combined,
    features = c(commongenes, commonviral, brz_fss, ifnb_uniq), assay = "RNA", slot = "data", disp.max = 3,
    label = FALSE, size = 5, group.by = "orig.ident"
) + scale_fill_gradientn(colors = c("black", "red")) + theme(axis.text.y = element_text(size = 6))
ggsave(file.path(results_folder, "DEheatmap_all_bulk_de.png"), width = 8, height = 8, dpi = 500)

 ## --cell type
 brain.combined$celltype_idents <- paste(brain.combined@meta.data$garnett_cluster_extend_lw, brain.combined@meta.data$orig.ident, sep = "_")
 Idents(brain.combined) <- brain.combined$celltype_idents

 DoHeatmap(brain.combined,
     features = c(commongenes, commonviral, brz_fss, ifnb_uniq), assay = "RNA", slot = "data", disp.max = 3,
     label = FALSE, size = 5, group.by = "celltype_idents"
 ) + scale_fill_gradientn(colors = c("black", "red")) + theme(axis.text.y = element_text(size = 6))
 ggsave(file.path(results_folder, "DEheatmap_all_bulk_de_celltypes.png"), width = 16, height = 8, dpi = 500)
 ggsave(file.path(results_folder, "DEheatmap_all_bulk_de_celltypes.svg"), width = 16, height = 8, dpi = 500)
 ggsave(file.path(results_folder, "DEheatmap_all_bulk_de_celltypes.pdf"), width = 16, height = 8, dpi = 500)


brain.combined@active.ident <- factor(brain.combined@active.ident, levels = rev(c(
    "NPC_IFNb_48hr", "NPC_BRZ_48hr", "NPC_FSS_48hr", "NPC_mock_48hr",
    "Early Neurons_mock_48hr",
    "Neurons_IFNb_48hr", "Neurons_BRZ_48hr", "Neurons_FSS_48hr", "Neurons_mock_48hr",
    "Neurons inhibitory_IFNb_48hr", "Neurons inhibitory_BRZ_48hr", "Neurons inhibitory_FSS_48hr", "Neurons inhibitory_mock_48hr",
    "Neurons excitatory_IFNb_48hr", "Neurons excitatory_BRZ_48hr", "Neurons excitatory_FSS_48hr", "Neurons excitatory_mock_48hr",
    "Astrocytes_IFNb_48hr", "Astrocytes_BRZ_48hr", "Astrocytes_FSS_48hr", "Astrocytes_mock_48hr",
    "Oligodendrocyte_IFNb_48hr", "Oligodendrocyte_BRZ_48hr", "Oligodendrocyte_FSS_48hr", "Oligodendrocyte_mock_48hr"
)))

# -- remove Early neurons
brain.combined_sub <- subset(x = brain.combined, idents = c(
    "NPC_IFNb_48hr", "NPC_BRZ_48hr", "NPC_FSS_48hr", "NPC_mock_48hr",
    "Neurons_IFNb_48hr", "Neurons_BRZ_48hr", "Neurons_FSS_48hr", "Neurons_mock_48hr",
    "Neurons inhibitory_IFNb_48hr", "Neurons inhibitory_BRZ_48hr", "Neurons inhibitory_FSS_48hr", "Neurons inhibitory_mock_48hr",
    "Neurons excitatory_IFNb_48hr", "Neurons excitatory_BRZ_48hr", "Neurons excitatory_FSS_48hr", "Neurons excitatory_mock_48hr",
    "Astrocytes_IFNb_48hr", "Astrocytes_BRZ_48hr", "Astrocytes_FSS_48hr", "Astrocytes_mock_48hr",
    "Oligodendrocyte_IFNb_48hr", "Oligodendrocyte_BRZ_48hr", "Oligodendrocyte_FSS_48hr", "Oligodendrocyte_mock_48hr"
))

DotPlot(brain.combined_sub, features = c(commongenes, commonviral, brz_fss, ifnb_uniq), assay = "RNA", cols = c("black", "#eb3306"), col.min = 0, scale = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9)) + theme(legend.position = "top") +
    theme(legend.text = element_text(size = 8))
ggsave(file.path(results_folder, "DEdotplot_all_bulk_de_celltypes.png"), width = 15, height = 6.5, units = "in", dpi = 300)


results_folder <- generate_folder(file.path(hfb_results, "DE_results", "Celltype_ZikaVsMck"))

DENPCsIFNbvsMock <- get_DE_between_conditions("NPC_IFNb_48hr", "NPC_mock_48hr", "NPC_IFNbvsMock", brain.combined, results_folder, fontsize = 9)
DENPCsBRZvsMock <- get_DE_between_conditions("NPC_BRZ_48hr", "NPC_mock_48hr", "NPC_BRZvsMock", brain.combined, results_folder, fontsize = 9)
DENPCsFSSvsMock <- get_DE_between_conditions("NPC_FSS_48hr", "NPC_mock_48hr", "NPC_FSSvsMock", brain.combined, results_folder, fontsize = 9)

DEENIFNbvsMock <- get_DE_between_conditions("Early Neurons_IFNb_48hr", "Early Neurons_mock_48hr", "Early Neurons_IFNbvsMock", brain.combined, results_folder, fontsize = 9)
DEENBRZvsMock <- get_DE_between_conditions("Early Neurons_BRZ_48hr", "Early Neurons_mock_48hr", "Early Neurons_BRZvsMock", brain.combined, results_folder, fontsize = 9)
DEENFSSvsMock <- get_DE_between_conditions("Early Neurons_FSS_48hr", "Early Neurons_mock_48hr", "Early Neurons_FSSvsMock", brain.combined, results_folder, fontsize = 9)

DENIFNbvsMock <- get_DE_between_conditions("Neurons_IFNb_48hr", "Neurons_mock_48hr", "Neurons_IFNbvsMock", brain.combined, results_folder, fontsize = 9)
DENBRZvsMock <- get_DE_between_conditions("Neurons_BRZ_48hr", "Neurons_mock_48hr", "Neurons_BRZvsMock", brain.combined, results_folder, fontsize = 9)
DENFSSvsMock <- get_DE_between_conditions("Neurons_FSS_48hr", "Neurons_mock_48hr", "Neurons_FSSvsMock", brain.combined, results_folder, fontsize = 9)

DENEIFNbvsMock <- get_DE_between_conditions("Neurons excitatory_IFNb_48hr", "Neurons excitatory_mock_48hr", "Neurons excitatory_IFNbvsMock", brain.combined, results_folder, fontsize = 9)
DENEBRZvsMock <- get_DE_between_conditions("Neurons excitatory_BRZ_48hr", "Neurons excitatory_mock_48hr", "Neurons excitatory_BRZvsMock", brain.combined, results_folder, fontsize = 9)
DENEFSSvsMock <- get_DE_between_conditions("Neurons excitatory_FSS_48hr", "Neurons excitatory_mock_48hr", "Neurons excitatory_FSSvsMock", brain.combined, results_folder, fontsize = 9)

DENIIFNbvsMock <- get_DE_between_conditions("Neurons inhibitory_IFNb_48hr", "Neurons inhibitory_mock_48hr", "Neurons inhibitory_IFNbvsMock", brain.combined, results_folder, fontsize = 9)
DENIBRZvsMock <- get_DE_between_conditions("Neurons inhibitory_BRZ_48hr", "Neurons inhibitory_mock_48hr", "Neurons inhibitory_BRZvsMock", brain.combined, results_folder, fontsize = 9)
DENIFSSvsMock <- get_DE_between_conditions("Neurons inhibitory_FSS_48hr", "Neurons inhibitory_mock_48hr", "Neurons inhibitory_FSSvsMock", brain.combined, results_folder, fontsize = 9)

DEAIFbvsMock <- get_DE_between_conditions("Astrocytes_IFNb_48hr", "Astrocytes_mock_48hr", "Astrocytes_IFNbvsMock", brain.combined, results_folder, fontsize = 9)
DEABRZvsMock <- get_DE_between_conditions("Astrocytes_BRZ_48hr", "Astrocytes_mock_48hr", "Astrocytes_BRZvsMock", brain.combined, results_folder,  fontsize = 9)
DEAFSSvsMock <- get_DE_between_conditions("Astrocytes_FSS_48hr", "Astrocytes_mock_48hr", "Astrocytes_FSSvsMock", brain.combined, results_folder, fontsize = 9)

DEONIFbvsMock <- get_DE_between_conditions("Oligodendrocyte_IFNb_48hr", "Oligodendrocyte_mock_48hr", "Oligodendrocyte_IFNbvsMock", brain.combined, results_folder, fontsize = 9)
DEOBRZvsMock <- get_DE_between_conditions("Oligodendrocyte_BRZ_48hr", "Oligodendrocyte_mock_48hr", "Oligodendrocyte_BRZvsMock", brain.combined, results_folder, fontsize = 9)
DEOFSSvsMock <- get_DE_between_conditions("Oligodendrocyte_FSS_48hr", "Oligodendrocyte_mock_48hr", "Oligodendrocyte_FSSvsMock", brain.combined, results_folder, fontsize = 9)

genelisttotal <- list(
    "DE Astros BRZ" = rownames(DEABRZvsMock), "DE Astros FSS" = rownames(DEAFSSvsMock),
    "DE Astros IFN" = rownames(DEAIFbvsMock),
    "DE NeuronsIn BRZ" = rownames(DENIBRZvsMock), "DE NeuronsIn FSS" = rownames(DENIFSSvsMock),
    "DE NeuronsIn IFN" = rownames(DENIIFNbvsMock),
    "DE NeuronsEx BRZ" = rownames(DENEBRZvsMock), "DE NeuronsEx FSS" = rownames(DENEFSSvsMock),
    "DE NeuronsEx IFN" = rownames(DENEIFNbvsMock),
    "DE NPC BRZ" = rownames(DENPCsBRZvsMock), "DE NPC FSS" = rownames(DENPCsFSSvsMock),
    "DE NPC IFN" = rownames(DENPCsIFNbvsMock)
)

png(file.path(hfb_results, "DE_results", "Upsetgraph_DEgenes_degree.png"), res = 250, units = "in", width = 8, height = 5)
upset(fromList(genelisttotal), sets = rev(c(
    "DE NPC IFN", "DE NPC BRZ", "DE NPC FSS", "DE NeuronsIn IFN", "DE NeuronsIn BRZ", "DE NeuronsIn FSS", "DE NeuronsEx IFN",
    "DE NeuronsEx BRZ", "DE NeuronsEx FSS", "DE Astros IFN",
    "DE Astros BRZ", "DE Astros FSS"
)), keep.order = TRUE, order.by = "degree")
dev.off()
svg(file.path(hfb_results, "DE_results", "Upsetgraph_DEgenes.svg"))
# png(file.path(hfb_results, "DE_results", "Upsetgraph_DEgenes.png"), res = 250, units = "in", width = 8, height = 5)
upset(fromList(genelisttotal), sets = rev(c(
    "DE NPC IFN", "DE NPC BRZ", "DE NPC FSS", "DE NeuronsIn IFN", "DE NeuronsIn BRZ", "DE NeuronsIn FSS", "DE NeuronsEx IFN",
    "DE NeuronsEx BRZ", "DE NeuronsEx FSS", "DE Astros IFN",
    "DE Astros BRZ", "DE Astros FSS"
)), keep.order = TRUE,  order.by = "freq")
dev.off()


##### ---Expression of genes of interest---#####

results_folder <- generate_folder(file.path(hfb_results, "DotPlots"))
genes <- c("DDX58", "IFIH1", "TLR3", "TLR4", "IRF3", "IRF5", "IRF7", "IFNAR1", "IFNAR2", "TLR7", "NFKB1") # IFH1=MDA5, DDX58=RIGI

DotPlot(brain.combined, features = genes, group.by = "celltype_idents", assay = "RNA", scale = TRUE) + coord_flip() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    scale_colour_gradient2(midpoint = 0, mid = "gray", high = "red", low = "blue")

ggsave(file.path(results_folder, "PRRS_dotplot_scaled4pub.png"), width = 9, height = 6, units = "in", dpi = 500)
ggsave(file.path(results_folder, "PRRS_dotplot_scaled4pub.svg"), width = 9, height = 6, units = "in", dpi = 500)
ggsave(file.path(results_folder, "PRRS_dotplot_scaled4pub.pdf"), width = 9, height = 6, units = "in", dpi = 500)


ifngenes <- c("IFIT1", "IFITM1", "IFITM3", "OAS1", "MX1", "ISG15", "USP18")
DotPlot(brain.combined, features = ifngenes, group.by = "celltype_idents", assay = "RNA", scale = TRUE) + coord_flip() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_colour_gradient2(midpoint = 0, mid = "gray", high = "red", low = "blue")
ggsave(file.path(results_folder, "IFIT_dotplot_scaled.png"), width = 9, height = 6, units = "in", dpi = 300)
ggsave(file.path(results_folder, "IFIT_dotplot_scaled.svg"), width = 9, height = 6, units = "in", dpi = 300)
ggsave(file.path(results_folder, "IFIT_dotplot_scaled.pdf"), width = 9, height = 6, units = "in", dpi = 300)
