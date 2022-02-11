library(ggplot2)

theme_Publication <- function(base_size = 14) {
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
            axis.title = element_text(size = rel(1)),
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
            legend.key.size = unit(0.5, "cm"),
            legend.margin = unit(0.1, "cm"),
            legend.title = element_text(size = 10, face = "bold"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(size = 10, face = "bold")
        ))
}

SC_DE_barplot <- function(DE, ident_1, ident_2, fontsize=9,
    high.color="darkgreen", middle.color="white",
    low.color="darkorange3", horizontal=FALSE, sort=TRUE) {
    DE <- DE
    DE$genes <- rownames(DE)
    DE$pct.diff <- DE$pct.1 - DE$pct.2
    DE$pct.diff.abs <- abs(DE$pct.1 - DE$pct.2)
    if (isTRUE(horizontal)) {
        if (isTRUE(sort)) {
            DE <- DE[order(-DE$avg_logFC), ]
        }
        pl <- ggplot( data = DE, aes(y = avg_logFC,
                x = factor(genes, levels = rownames(DE)),
                fill = pct.diff)) +
            geom_bar(stat = "identity", color = "black") +
            theme_Publication() +
            scale_fill_gradient2(low = low.color, mid = middle.color, high = high.color, midpoint = 0, limits = c(-1, 1)) +
            labs(fill = paste0("diff in % of\ncells exp genes:\n", ident_1, " -\n", ident_2),
                y = "lfc", x = "genes") + theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size = fontsize))
    } else {
        if (isTRUE(sort)) {
            DE <- DE[order(DE$avg_logFC), ]
        }
        pl <- ggplot(data = DE,
                    aes(x = avg_logFC,
                    y = factor(genes, levels=rownames(DE)),
                    fill=pct.diff)) +
            geom_bar(stat="identity", color="black") +
            theme_Publication() +
        scale_fill_gradient2(low = low.color, mid = middle.color, high =high.color, midpoint = 0, limits = c(-1, 1)) +
        labs(fill=paste0("diff in % of\ncells exp genes:\n", ident_1, " -\n", ident_2),
            x = "lfc", y = "genes") +
        theme(axis.text.y = element_text(size = fontsize))
    }
    return(pl)
}

SC_AE_barplot <- function(DE, AE, fill_val, fontsize = 9,
                          high.color = "firebrick4", 
                          low.color = "white", horizontal = FALSE, sort = TRUE) {
    DE$genes <- rownames(DE)
    if (isTRUE(horizontal)) {
        if (isTRUE(sort)) {
            DE <- DE[order(-DE$avg_logFC), ]
        }
        pl <- ggplot(data = DE, aes(
            y = DE[[AE]],
            x = factor(genes, levels = rownames(DE)),
            fill = DE[[fill_val]]
        )) +
            geom_bar(stat = "identity", color = "black") +
            theme_Publication() +
            scale_fill_gradient(low = low.color, high = high.color, limits=c(0,1)) +
            labs(
                fill = paste0("% of\ncells exp genes:\n", AE),
                y = "Average Expression", x = "genes"
            ) +
            theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2, size = fontsize))
    } else {
        if (isTRUE(sort)) {
            DE <- DE[order(DE$avg_logFC), ]
        }
        pl <- ggplot(
            data = DE,
            aes(
                x = DE[[AE]],
                y = factor(genes, levels = rownames(DE)),
                fill = DE[[fill_val]]
            )
        ) +
            geom_bar(stat = "identity", color = "black") +
            theme_Publication() +
            scale_fill_gradient(low = low.color, high = high.color, limits=c(0,1)) +
            labs(
                fill = paste0("% of\ncells exp genes:\n", AE),
                x = "Average Expression", y = "genes"
            ) +
            theme(axis.text.y = element_text(size = fontsize))
    }
    return(pl)
}

SC_DE_dotplot <- function(DE, ident_1, ident_2, compare, fontsize=9,
    high.color="red", middle.color="white",
    low.color="blue") {
        DE$genes <- rownames(DE)
        DE$pct.diff <- DE$pct.1 - DE$pct.2
        DE$pct.diff.abs <- abs(DE$pct.1 - DE$pct.2)
        DE$Exp <- rep(compare, length(rownames(DE)))
        exp_max <- c()
        for (i in DE$pct.diff) {
            if (i > 0) {
                exp_max <- c(exp_max, ident_1)
            } else {
                exp_max <- c(exp_max, ident_2)
            }
        }
        DE$Exp_max <- exp_max
        DE <- DE[order(DE$avg_logFC), ]
        pl <- ggplot(data = DE, aes(
            x = Exp, y = factor(genes, levels = rownames(DE)),
            size = pct.diff.abs, fill = avg_logFC, colour = factor(Exp_max, levels = c(ident_1, ident_2))
        )) +
        geom_point(shape = 21, width = 1, stroke = 1) +
        theme_Publication() +
        scale_colour_manual(values = c("firebrick", "steelblue4")) +
        scale_fill_gradient2(low = low.color, mid =middle.color, high =high.color, midpoint = 0) +
        labs(
            fill = "log fold change", colour = "conditions",
            size = "diff in % of cells", y = "genes", x = "comparison"
        ) +
        scale_size_continuous(range = c(2, dotscale + 1)) +
        theme(axis.text.y = element_text(size = fontsize))
    return(pl)
}
