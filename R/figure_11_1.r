library(TDbook)
taxa <- df_alltax_info
dt <- df_difftax
p <- ggdiffclade(obj = taxa,
                nodedf = dt,
                factorName = "DIAGNOSIS",
                layout = "radial",
                skpointsize = 0.6,
                cladetext = 2,
                linewd = 0.2,
                taxlevel = 3,
                # This argument is to remove the branch of unknown taxonomy.
                reduce = TRUE) +
    scale_fill_manual(values = c("#00AED7", "#009E73")) +
    guides(color = guide_legend(
        keywidth = 0.1,
        keyheight = 0.6,
        order = 3,
        ncol = 1
    )) +
    theme(
        panel.background = element_rect(fill = NA),
        legend.position = "right",
        plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.spacing.y = unit(0.5, "cm"),
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 5.5),
        legend.box.spacing = unit(0, "cm")
    )

print(p)

ggsave("reports/figures/figure_11_1.png", plot = p, width = 12, height = 12, dpi = 300)


