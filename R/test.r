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
        plot.margin = margin(0, 0, 0, 0),
        legend.spacing.y = unit(0.02, "cm"),
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 5.5),
        legend.box.spacing = unit(0.02, "cm")
    )

print(p)


diffclad_obj <- diff_analysis(
                                subset_physeq,
                                classgroup = "Group",
                                mlfun = "lda",
                                ldascore = 3,
                                normalization = 1000000,
                                clmin=1
                              )

df_diffclad_obj <- as.data.frame(diffclad_obj)



diffclade_obj <- grt_ggdiffclade_obj(subset_physeq)
df_diffclade_obj <- as.data.frame(diffclade_obj)
tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")


