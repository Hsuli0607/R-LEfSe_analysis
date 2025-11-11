source("R/check_and_load_packages.r")
source("R/generate_physeq_obj.r")
source("R/subset_physeq_obj.r")
source("R/ggdiffclade_obj.r")
source("R/ggdiffclade_nodedf.r")

)


physeq_obj <- generate_physeq_obj(raw_data = "data/processed/ALL-update.csv",
                                  META_cols = c("index", "Group"),
                                  sample_id_col = "index")

subset_physeq <- subset_physeq_obj(physeq_obj, sample_col = "Group",
                                   subset_vals = c("EHI", "EHI.abt"),
                                   abundance_filter_threshold = 1000,
                                   prevalence_filter_threshold = 3)


