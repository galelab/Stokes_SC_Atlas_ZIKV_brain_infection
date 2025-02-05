######## --Libraries--########
library(ggplot2)
library(stringr)
library(data.table)
library(Seurat)
####### --Functions --########
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

theme_Publication <- function(base_size = 14, base_family = "arial") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour = "#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.4, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(size = 10, face = "bold"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

theme_sctour <- function(base_size = 14, base_family = "helvetica") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2, size=10),
            axis.title.x = element_text(vjust = -0.2, size=10),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.4, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(size = 10, face = "bold"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "white", fill = "white"),
            strip.text = element_text(face = "bold")
        ))
}

####### --Main Code --########
input_folder <- "2.TrajectoryAnalysisSCtouriNPC_NPC"
results_folder <- file.path("2.TrajectoryAnalysisSCtouriNPC_NPC", "figs")
generate_folder(results_folder)
inpc <- readRDS("/share/lwhitmo/projects/Stokes_HFB_iNPC_Analysis/zika_scViralQuant_aligned/inpc_results/Convert4SCtour/inpc.combinedsub.RDS")
df <- read.csv(file.path(input_folder, "umapcoordinates.csv"), row.names=1)
df$ptime <- as.numeric(df$ptime)

results_folder=file.path(input_folder,"figs")
generate_folder(results_folder)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color=ptime)) +
    geom_point(size=0.5) +
    scale_color_viridis_c(option = "inferno") +
    theme_sctour()
ggsave(file.path(results_folder, "ptime.png"),width=4.5, height=4, units = "in", dpi=400)
ggsave(file.path(results_folder, "ptime.pdf"),width=4.5, height=4, units = "in", dpi=400)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, color=ptime)) +
    geom_point(size=0.5) + facet_wrap(~orig.ident, ncol=2) +
    scale_color_viridis_c(option = "inferno") +
    theme_sctour()
ggsave(file.path(results_folder, "ptimesplit.png"),width=8, height=8, units = "in", dpi=400)
ggsave(file.path(results_folder, "ptimesplit.pdf"),width=8, height=8, units = "in", dpi=400)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, color=ZIKA)) +
    geom_point(size=0.5) + facet_wrap(~orig.ident, ncol=2) +
    theme_sctour()
ggsave(file.path(results_folder, "ptimesplitzika.png"),width=8, height=8, units = "in", dpi=400)
ggsave(file.path(results_folder, "ptimesplitzika.pdf"),width=8, height=8, units = "in", dpi=400)

####### --genes of interest --########
genes <- c("RFC3", "PCNA", "MCM4")
g<- inpc@assays$SCT$data[genes,]
g <- as.data.frame(g)
g <- t(g)

if (all.equal(rownames(g), rownames(df)==FALSE)) {
    message("MASSIVE WARNING ROWNAMES NOT EQUAL TRYING TO FIX")
    g <- g[rownames(df), ]
    if (isTRUE(all.equal(rownames(g), rownames(df)))){
        message("FIXED")

    }
}
df <- cbind(df, g)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, color=RFC3)) +
    geom_point(size=0.5) + facet_wrap(~orig.ident, ncol=2) +
    scale_color_viridis_c(option = "mako") +
    theme_sctour()
ggsave(file.path(results_folder, "RFC3ptimesplit.png"),width=8, height=8, units = "in", dpi=400)
ggsave(file.path(results_folder, "RFC3ptimesplit.pdf"),width=8, height=8, units = "in", dpi=400)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color=PCNA)) +
    geom_point(size=0.5) + facet_wrap(~orig.ident, ncol=2) +
    scale_color_viridis_c(option = "mako") +
    theme_sctour()
ggsave(file.path(results_folder, "PCNAptimesplit.png"),width=8, height=8, units = "in", dpi=400)
ggsave(file.path(results_folder, "PCNAptimesplit.pdf"),width=8, height=8, units = "in", dpi=400)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color=MCM4)) +
    geom_point(size=0.5) + facet_wrap(~orig.ident, ncol=2) +
    scale_color_viridis_c(option = "mako") +
    theme_sctour()
ggsave(file.path(results_folder, "MCM4ptimesplit.png"),width=8, height=8, units = "in", dpi=400)
ggsave(file.path(results_folder, "MCM4ptimesplit.pdf"),width=8, height=8, units = "in", dpi=400)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, color=RFC3)) +
    geom_point(size=0.5) +
    scale_color_viridis_c(option = "mako") +
    theme_sctour()
ggsave(file.path(results_folder, "RFC3ptime.png"),width=4.5, height=4, units = "in", dpi=400)
ggsave(file.path(results_folder, "RFC3ptime.pdf"),width=4.5, height=4, units = "in", dpi=400)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, color=PCNA)) +
    geom_point(size=0.5) +
    scale_color_viridis_c(option = "mako") +
    theme_sctour()
ggsave(file.path(results_folder, "PCNAptime.png"),width=4.5, height=4, units = "in", dpi=400)
ggsave(file.path(results_folder, "PCNAptime.pdf"),width=4.5, height=4, units = "in", dpi=400)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, color=MCM4)) +
    geom_point(size=0.5) +
    scale_color_viridis_c(option = "mako") +
    theme_sctour()
ggsave(file.path(results_folder, "MCM4ptime.png"),width=4.5, height=4, units = "in", dpi=400)
ggsave(file.path(results_folder, "MCM4ptime.pdf"),width=4.5, height=4, units = "in", dpi=400)
