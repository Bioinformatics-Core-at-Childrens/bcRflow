######################################################################################################
#                                                                                                    #
# ░▒▓███████▓▒░ ░▒▓██████▓▒░░▒▓███████▓▒░░▒▓████████▓▒░▒▓█▓▒░      ░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░ #
# ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░ #
# ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░ #
# ░▒▓███████▓▒░░▒▓█▓▒░      ░▒▓███████▓▒░░▒▓██████▓▒░ ░▒▓█▓▒░     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░ #
# ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░ #
# ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░ #
# ░▒▓███████▓▒░ ░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓████████▓▒░▒▓██████▓▒░ ░▒▓█████████████▓▒░  #
#                                                                                                    #
######################################################################################################

# load packages, load custom functions to environment:
source("./utils.R")

# read in user input:
option_list <- list(
  make_option(
    c("-r", "--reports"),
    type = "character",
    default = "../MiXCR",
    help = "Path to MiXCR output (IGH) [default %default]"
  ),
  make_option(
    c("-s", "--species"),
    type = "character",
    default = "hsa",
    help = "Species of samples [default %default]"
  ),
  make_option(
    c("-m", "--metadata"),
    type = "character",
    default = file.path("./samplesList.csv"),
    # points to symLink of metadata table in working dir, file must be named samplesList.csv!
    help = "Path to sample metadata table [default %default]"
  ),
  make_option(
    c("-d", "--downsample"),
    type = "character",
    default = "TRUE",
    help = "Whether or not to downsample the data to the size of the smallest repertoire [default %default]"
  ),
  make_option(
    c("-t", "--threshold"),
    type = "double",
    default = 0.7,
    help = "Sequence similarity threshold for clonal clustering [default %default]"
  ),
  make_option(
    c("-o", "--output_dir"),
    type = "character",
    default = file.path("./immunarch"),
    help = "Output directory [default %default]"
  )
)
arguments <- parse_args(OptionParser(option_list = option_list))
arguments$downsample <- as.logical(arguments$downsample)

# argument handling:
species <- arguments$species
valid_species <- c("hsa", "mmu")
stopifnot(species %in% valid_species)

outdir <- file.path(arguments$output_dir)
metadata <- read.csv(arguments$metadata)
stopifnot(nrow(metadata) > 0)

# check if there's group metadata:
if (is.null(metadata$Group) |
    length(metadata$Group) != nrow(metadata)) {
  stop("Group metadata is missing or incomplete! Please check, and try again.")
}
rownames(metadata) <- metadata$SampleID

# check and create output directory:
if (!dir.exists(outdir)) {
  dir.create(outdir, showWarnings = T)
} else {
  print("Output directory already exists...")
}
setwd(outdir)

# diversity dir for diversity metric figures:
if (!dir.exists(file.path("./diversity"))){
  dir.create(file.path("./diversity"), showWarnings = T)
} else {
  print("Diversity directory already exists...")
}

# chord dir for V-J usage chord plots:
if (!dir.exists(file.path("./chords"))){
  dir.create(file.path("./chords"), showWarnings = T)
} else {
  print("Chords directory already exists...")
}

# network dir for convergent network plots:
if (!dir.exists(file.path("./network"))){
  dir.create(file.path("./network"), showWarnings = T)
} else {
  print("Network directory already exists...")
}

# list MiXCR output reports:
IGH_reports <-
  list.files(
    file.path(arguments$reports),
    full.names = T,
    pattern = ".tsv",
    recursive = T
  )

# load data using Immunarch:
immdata <- repLoad(file.path(IGH_reports))
imm.meta <- metadata[, c("SampleID", "Group")]
colnames(imm.meta) <-  c("Sample", "Group")
rownames(imm.meta) <- imm.meta$Sample
imm.meta <- imm.meta[names(immdata$data),]
immdata$meta <- tibble(imm.meta)

if(arguments$downsample){
  immdata$data <- repSample(immdata$data,
                       .method = "downsample")
}

# #sample order for plots (custom for manuscript):
# plot_order <- c("Healthy","Exposed","Mild","Severe")
# immdata$meta$Group <- factor(immdata$meta$Group)

# Preprocessing MiXCR output:
# extract the top hit for each V, D, J and C alignment:
for (i in 1:length(immdata$data)) {
  immdata$data[[i]]$sample_id <- names(immdata$data)[i]
  immdata$data[[i]]$V.name <- sapply(immdata$data[[i]]$V.name, function(x) ifelse(length(unlist(strsplit(x, ","))) > 1, extract_first_hit(x), x))
  immdata$data[[i]]$D.name <- sapply(immdata$data[[i]]$D.name, function(x) ifelse(length(unlist(strsplit(x, ","))) > 1, extract_first_hit(x), x))
  immdata$data[[i]]$J.name <- sapply(immdata$data[[i]]$J.name, function(x) ifelse(length(unlist(strsplit(x, ","))) > 1, extract_first_hit(x), x))
  immdata$data[[i]]$C.name <- sapply(immdata$data[[i]]$C.name, function(x) ifelse(length(unlist(strsplit(x, ","))) > 1, extract_first_hit(x), x))
}

# bind samples into one data frame:
data <- do.call(rbind, immdata$data)
data$group_id <- imm.meta[data$sample_id, "Group"]
data$group_id <- factor(data$group_id)

# genes w/o mutations, just convert to NA:
data[data == "region_not_covered"] <- NA

# reshape data to include a column of Vgene name and allele for CSR calculation:
data$CloneID <- data$Clone.ID
data$CloneID <- paste(data$sample_id, data$CloneID, sep = "_")
data$vgene_allele <- str_extract(data$V.name, ".+?(?=\\x28)")
data$cgene <- str_extract(data$C.name, ".+?(?=\\*)")
data$vgene <- str_extract(data$V.name, ".+?(?=\\*)")
data$dgene <- str_extract(data$D.name, ".+?(?=\\*)")
data$jgene <- str_extract(data$J.name, ".+?(?=\\*)")

#split data by group:
data_split         <- group_split(data, group_id, .keep = TRUE)
names(data_split)  <- data$group_id %>% unique() %>% sort()

#####################################################
#IGH-V Gene Usage:
species <- str_sub(arguments$species, end = -2)
ighv_gu <- geneUsage(immdata$data, paste0(species, ".ighv"))
ighv_gu <- ighv_gu %>% replace(is.na(.), 0) %>% data.frame()

# set rownames, drop clones with no aligned V-Gene ("NA")
rownames(ighv_gu) <- ighv_gu$Names
ighv_gu$Names <- NULL
ighv_gu <- ighv_gu[rownames(ighv_gu)!= "NA",]

## ComplexHeatmap:
# create column annotation df from Group metadata:
ighv_annotation_df <-
  data.frame("Sample" = colnames(ighv_gu), "Group" = imm.meta[colnames(ighv_gu), "Group"])
ighv_annotation_df$Group <- factor(ighv_annotation_df$Group)
ighv_annotation_df <-
  ighv_annotation_df[order(ighv_annotation_df$Group),]

rownames(ighv_annotation_df) <- ighv_annotation_df$Sample
immdata$data <- immdata$data[ighv_annotation_df$Sample]
immdata$meta <- immdata$meta %>% arrange(match(Sample, rownames(ighv_annotation_df)))

# sample colors:
sample_cols <-
  custom_colors$discrete[1:length(unique(ighv_annotation_df$Sample))]
names(sample_cols) <- ighv_annotation_df$Sample %>% unique()

# group colors:
group_cols <-
  colors_dutch[1:length(unique(ighv_annotation_df$Group))]
names(group_cols) <- ighv_annotation_df$Group %>% unique()

# drop rows with no mapped IGH-V gene:
ighv_gu <- ighv_gu[row.names(ighv_gu) != "NA", , drop = FALSE]

# scale & center gene util counts:
ighv_gu <- scale(ighv_gu, center = T, scale = T)
col_fun = colorRamp2(c(min(ighv_gu), 0 , max(ighv_gu)), c("blue", "white", "red"))

ighv_gu <- ighv_gu[, ighv_annotation_df$Sample]

# plot with ComplexHeatmap:
ighv_gu_heatmap <- Heatmap(
  column_title = "IGH-V Gene Usage",
  column_title_gp = gpar(fontsize = 29, fontface = "bold"),
  cluster_columns =  F,
  cluster_rows = T,
  ighv_gu[1:50,],
  show_row_dend = T,
  show_column_dend = F,
  name = "Gene Usage",
  col = col_fun,
  top_annotation = columnAnnotation(
    Group = ighv_annotation_df$Group,
    col = list(Group = group_cols),
    show_annotation_name = TRUE,
    show_legend = c(TRUE),
    annotation_name_gp = gpar(fontsize = 20),
    annotation_legend_param = list(fontsize = 20, labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20, fontface = "bold"))
  ),
  column_labels = NULL,
  row_names_gp = gpar(fontsize = 20, lwd = 2),
  column_names_gp = gpar(fontsize = 20, hjust = 0.5),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 20, fontface = "bold"),
    labels_gp = gpar(fontsize = 20)
  ),
  show_column_names = F,
  column_dend_side = "bottom",
  row_names_side = "left",
  row_names_centered = F
)
ighv_gu_heatmap

# save the gene usage heatmap:
tiff(
  "IGHV_gene_usage_heatmap.tiff",
  width = 16,
  height = 14,
  units = "in",
  compression = "lzw",
  res = 300,
  bg = "white"
)
ighv_gu_heatmap
dev.off()

## TODO: Users should be able to set the reference level of Group meta for Odds Ratio analyses:
# #Odds ratio vs Healthy:
# vgene_usage <-  geneUsage(immdata$data, paste0(species, ".ighv"))
# vgene_usage <- vgene_usage %>% replace(is.na(.), 0) %>% data.frame()
# rownames(vgene_usage) <- vgene_usage$Names
# vgene_usage$Names <- NULL
# vgene_usage <- vgene_usage[rownames(vgene_usage) != "NA", ]
# colnames(vgene_usage) <- metadata[colnames(vgene_usage), "Group"]
# vgene_usage <-
#   cbind(sapply(unique(colnames(vgene_usage)[duplicated(colnames(vgene_usage))]), function(x)
#     rowSums(vgene_usage[, grepl(paste(x, "$", sep = ""), colnames(vgene_usage))])),
#     vgene_usage[,!duplicated(colnames(vgene_usage)) &
#                   !duplicated(colnames(vgene_usage), fromLast = TRUE)]) %>% data.frame()
# 
# #Define a function to compute odds ratio for each gene
# calculate_odds_ratio <- function(response_counts, control_counts, genes) {
#   # Create a data frame to store results
#   results <- data.frame(
#     Gene = genes,
#     Odds_Ratio = rep(NA, length(response_counts)),
#     Lower_CI =   rep(NA, length(response_counts)),
#     Upper_CI =   rep(NA, length(response_counts)),
#     P_Value =    rep(NA, length(response_counts)),
#     Adjusted_P_Value = rep(NA, length(response_counts))
#   )
#   
#   for (i in 1:length(response_counts)) {
#     # Create a contingency table
#     contingency_table <- matrix(
#       c(response_counts[i], sum(response_counts) - response_counts[i],
#         control_counts[i], sum(control_counts) - control_counts[i]),
#       nrow = 2, byrow = TRUE
#     )
#     
#     # Add a small constant value to each cell to avoid division by zero
#     contingency_table <- contingency_table + 1
#     
#     # Perform Fisher's exact test
#     fisher_test <- fisher.test(contingency_table)
#     
#     # Get the p-value from the test
#     p_value <- fisher_test$p.value
#     
#     # Calculate odds ratio
#     odds_ratio <- fisher_test$estimate
#     
#     # Calculate standard error of log odds ratio
#     se_log_odds_ratio <- sqrt(1 / contingency_table[1, 1] + 
#                                 1 / contingency_table[1, 2] + 
#                                 1 / contingency_table[2, 1] + 
#                                 1 / contingency_table[2, 2])
#     
#     # Calculate the log odds ratio
#     log_odds_ratio <- log(odds_ratio)
#     
#     # Calculate confidence intervals for log odds ratio
#     z_value <- qnorm(1 - 0.05 / 2)
#     lower_ci <- log_odds_ratio - z_value * se_log_odds_ratio
#     upper_ci <- log_odds_ratio + z_value * se_log_odds_ratio
#     
#     # Convert confidence intervals back to odds ratio scale
#     lower_ci_odds <- exp(lower_ci)
#     upper_ci_odds <- exp(upper_ci)
#     
#     # Store the results in the data frame
#     results$Odds_Ratio[i] <- odds_ratio
#     results$P_Value[i] <- p_value
#     results$Lower_CI[i] <- lower_ci_odds
#     results$Upper_CI[i] <- upper_ci_odds
#   }
#   
#   # Perform Bonferroni correction
#   results$Adjusted_P_Value <- pmin(1, results$P_Value * length(response_counts))
#   
#   return(results)
# }
# 
# # Calculate odds ratio for V-gene usage:
# odds_severe  <- calculate_odds_ratio(vgene_usage$Severe, vgene_usage$Healthy, genes = rownames(vgene_usage))
# odds_mild    <- calculate_odds_ratio(vgene_usage$Mild, vgene_usage$Healthy, genes = rownames(vgene_usage))
# odds_exposed <- calculate_odds_ratio(vgene_usage$Exposed, vgene_usage$Healthy, genes = rownames(vgene_usage))
# 
# 
# odds_severe_filtered <- odds_severe[(odds_severe$Upper_CI - odds_severe$Lower_CI < 3),]
# odds_severe_filtered$significant <- odds_severe_filtered$P_Value < 0.05
# odds_severe_filtered <- odds_severe_filtered[order(odds_severe_filtered$significant),]
# 
# odds_mild_filtered <- odds_mild[(odds_mild$Upper_CI - odds_mild$Lower_CI < 3),]
# odds_mild_filtered$significant <- odds_mild_filtered$P_Value < 0.05
# odds_mild_filtered <- odds_mild_filtered[order(odds_mild_filtered$significant),]
# 
# odds_exposed_filtered <- odds_exposed[(odds_exposed$Upper_CI - odds_exposed$Lower_CI < 3),]
# odds_exposed_filtered$significant <- odds_exposed_filtered$P_Value < 0.05
# odds_exposed_filtered <- odds_exposed_filtered[order(odds_exposed_filtered$significant),]
# 
# odds_plot_severe <- ggplot(odds_severe_filtered, aes(x = Odds_Ratio, y = Gene, color = significant)) +
#   geom_point() +  # Add points for the odds ratio
#   geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0) +  # Add horizontal error bars for confidence intervals
#   geom_text(aes(label = sprintf("p = %.4f", P_Value)), hjust = -1, vjust = 0.5) +  # Add annotations for p-values
#   labs(x = "Odds Ratio", y = "Gene") +  # Label the axes
#   theme_minimal() +  # Use minimal theme
#   theme(axis.text.y = element_text(size = 8),  # Adjust size of y-axis text
#         axis.text.x = element_text(size = 8),  # Adjust size of x-axis text
#         axis.title = element_text(size = 10))+  # Adjust size of axis titles
#   ggtitle("Severe versus Healthy") + theme_bw()
# 
# 
# odds_plot_severe <- ggplot(odds_severe_filtered, aes(x = Odds_Ratio, y = Gene, color = significant)) +
#   geom_point(size = 3) +
#   geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0, color = "black", position = position_dodge(width = 0.5)) +
#   #geom_text(aes(label = sprintf("p = %.5f", P_Value)), hjust = -0.1, vjust = 0.25,color = "black", size = 3,  position = position_dodge(width = 0.5)) +
#   labs(x = "Odds Ratio", y = "V-Gene", title = "V-Gene Usage Analysis", subtitle = "Severe versus Healthy") +
#   scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"), name = "P-value < 0.05") +
#   theme_bw() +
#   theme(
#     plot.title = element_text(family = "Arial", size = 18, face = "bold"),
#     axis.title = element_text(family = "Arial", size = 14),
#     axis.text = element_text(family = "Arial", size = 12),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor = element_blank()
#   )
# 
# odds_plot_mild <- ggplot(odds_mild_filtered, aes(x = Odds_Ratio, y = Gene, color = significant)) +
#   geom_point(size = 3) +
#   geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0, color = "black", position = position_dodge(width = 0.5)) +
#   #geom_text(aes(label = sprintf("p = %.5f", P_Value)), hjust = -0.1, vjust = 0.25,color = "black", size = 3,  position = position_dodge(width = 0.5)) +
#   labs(x = "Odds Ratio", y = "V-Gene", title = "V-Gene Usage Analysis", subtitle = "Mild versus Healthy") +
#   scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"), name = "P-value < 0.05") +
#   theme_bw() +
#   theme(
#     plot.title = element_text(family = "Arial", size = 18, face = "bold"),
#     axis.title = element_text(family = "Arial", size = 14),
#     axis.text = element_text(family = "Arial", size = 12),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor = element_blank()
#   )
# 
# odds_plot_exposed <- ggplot(odds_exposed_filtered, aes(x = Odds_Ratio, y = Gene, color = significant)) +
#   geom_point(size = 3) +
#   geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0, color = "black", position = position_dodge(width = 0.5)) +
#   #geom_text(aes(label = sprintf("p = %.5f", P_Value)), hjust = -0.1, vjust = 0.25,color = "black", size = 3,  position = position_dodge(width = 0.5)) +
#   labs(x = "Odds Ratio", y = "V-Gene", title = "V-Gene Usage Analysis", subtitle = "Exposed versus Healthy") +
#   scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"), name = "P-value < 0.05") +
#   theme_bw() +
#   theme(
#     plot.title = element_text(family = "Arial", size = 18, face = "bold"),
#     axis.title = element_text(family = "Arial", size = 14),
#     axis.text = element_text(family = "Arial", size = 12),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor = element_blank()
#   )
# 
# tiff("odds_exposed_healthy.tiff", height = 8, width = 7, units = "in", compression = "lzw", res = 300)
# odds_plot_exposed
# dev.off()
# 
# tiff("odds_severe_healthy.tiff", height = 8.5, width = 7, units = "in", compression = "lzw", res = 300)
# odds_plot_severe
# dev.off()
# 
# tiff("odds_mild_healthy.tiff", height = 8, width = 7, units = "in", compression = "lzw", res = 300)
# odds_plot_mild
# dev.off()

#####################################################
# V-J paired gene usage chord plots:
## custom fx, see utils.R. Threshold is minimum number of hits required to plot gene, simplifies the chords.
v_j_counts <- lapply(immdata$data, function(x){v_j_matrix(x, threshold = 5)})

#plot the chords:
for(i in 1:length(v_j_counts)){
  tiff(paste0("./chords/",names(v_j_counts)[[i]],"_vj_chord.tiff"), width = 6, height = 6, bg = "white", units = "in", res =300, compression = "lzw")
  suppressMessages(vj_circos(v_j_counts[[i]][[1]])) # custom circos plots, see utils.R
  dev.off()
}

save.image("downstream.RData")
#####################################################
#Somatic Hyper Mutation (SHM):
## Calculate the number of mutations compared to Germline using Immunarch:
## Note: repGermline can only run on one thread - there's a known bug with the multithreading.
shm <- immdata$data %>%
  seqCluster(seqDist(immdata$data), .fixed_threshold = 3) %>%
  repGermline(.threads = 1) %>%
  repAlignLineage(.min_lineage_sequences = 2, .align_threads = 4, .nofail = TRUE) %>%
  repClonalFamily(.threads = 4, .nofail = TRUE) %>%
  repSomaticHypermutation(.threads = 4, .nofail = TRUE)

# Create an empty list to store the merged data frames of shm counts per sample:
merged_list <- list()

# Iterate over each pair of data frames in the lists and merge them:
for (key in names(shm)) {
  merged_data <- merge(shm[[key]], immdata$data[[key]], by = "Clone.ID", all.x = TRUE)
  merged_list[[key]] <- merged_data
}
merged_list <- do.call(rbind, merged_list)

# Estimate mutation rate over clonal nucleotide seq w/o CDR3
merged_list <- merged_list %>%
  mutate(Mutation.Rate = Mutations / (nchar(Sequence.x) - nchar(CDR3.nt.x)))

merged_list$Group  <- imm.meta[merged_list$sample_id, "Group"]
merged_list$C.name <- str_extract(merged_list$C.name, ".+?(?=\\*)")

shm_rate <-
  data.table(
    "Group"     = merged_list$Group,
    "Sample"    = merged_list$sample_id,
    "C.name"    = merged_list$C.name,
    "CDR3.nt"   = merged_list$CDR3.nt.x,
    "SHM.rate"  = merged_list$Mutation.Rate,
    "SHM.count" = merged_list$Mutations,
    "Sequence.length" =  str_length(merged_list$Sequence.x)
)

# drop NA values (div by 0?)
shm_rate <- shm_rate[complete.cases(shm_rate),]
shm_rate <- shm_rate %>% group_by(Sample, C.name) %>% summarize(Mutation.Rate = mean(SHM.rate))
shm_rate$Group <- imm.meta[shm_rate$Sample, "Group"]

# plot the SHM rates across groups, facet by C gene name:
shm.plot <- shm_rate %>% ggplot(aes(x = Group, y = Mutation.Rate, color = Group)) +
  geom_hline(yintercept=0.01, linetype='dotted', col = 'red') +
  ylab("Rate") +
  xlab(NULL) +
  geom_boxplot() + theme_bw(base_size = 20) + theme(panel.grid.major = element_blank(),
                                                    panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = group_cols, name = "Group") +
  facet_wrap(~ C.name)+
ggtitle("Somatic Hypermutation Rates", subtitle = "Immunarch repSomaticHypermutation")

# save the SHM plot:
tiff(
  "SHM_rates.tiff",
  width = 13,
  height = 10,
  bg = "white",
  units = "in",
  compression = "lzw",
  res = 300
)
shm.plot + theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) &
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

#####################################################
# CDR3 Length Distribution:
## extract CDR3 lengths across all groups
data %>% ggplot(aes(
  x = str_length(CDR3.aa),
  y = ..density..,
  fill = group_id
)) + geom_histogram(bins = 30) + ggtitle("CDR3 Length Distribution") +xlab("CDR3aa Length") + ylab("Proportion") + facet_wrap(~group_id) &
scale_fill_manual(values = group_cols, name = "Group")&
theme_classic() & scale_y_continuous(labels = scales::percent) & scale_x_binned(n.breaks = 30)

# violin plot of CDR3 AA length, split by Group meta:
cdr3.violin <-
  data %>% ggplot(aes(
    x = group_id,
    y = str_length(CDR3.aa),
    fill = group_id
  )) + geom_violin() +
  scale_fill_manual(values = group_cols, name = "Group") + ylab("CDR3 Length") + xlab(NULL)  +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("CDR3 Length Distribution")
cdr3.violin

# save the CDR3 length plot:
tiff(
  "CDR3_length_distribution.tiff",
  width = 8,
  height = 6,
  units = "in",
  bg = "white",
  compression = "lzw",
  res = 300
)
cdr3.violin + ggtitle("CDR3 Length Distribution") + theme_bw(base_size = 20) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + # add quants to vln plot
  scale_fill_manual(values = group_cols, name = "Group") + ylab("CDR3 Length") + xlab(NULL)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(panel.grid.major = element_blank(),
                                                                   panel.grid.minor = element_blank())
dev.off()

## Logo plot of CDR3 region across groups:
# split data df (all clones) by group:
data_split <- group_split(data, group_id, .keep = TRUE)
names(data_split)  <- data$group_id %>% unique() %>% sort() # immunarch sorts the group names alphabetically for some reason...
logoplots <- lapply(data_split, logoplots_fun)

# Save logo plots as a multi-page PDF:
pdf(
  "CDR3_logoplots.pdf",
  width = 8,
  height = 4,
  bg = "white"
)
for (i in 1:length(logoplots)) {
  print(
    logoplots[[i]] + ggtitle(paste0(
      "CDR3 Amino Acid Composition: ", names(logoplots)[[i]]
    )) +
      theme_bw(base_size = 20) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      scale_y_continuous(labels = scales::percent)
  )
}
dev.off()

#####################################################
#C-gene Isotype Distribution:
if_2 <- data[!is.na(data$cgene), ]  %>%
  mutate(cgene = factor(cgene),
         group_id = factor(group_id)) %>%
  group_by(cgene, group_id) |>
  summarise(n = sum(round(Clones), na.rm = TRUE)) %>%
  group_by(group_id) %>%
  mutate(pct = prop.table(n)) %>%
  ggplot(aes(pct, cgene, fill = group_id)) +
  geom_col(position = "dodge") +
  ggtitle("Isotype Frequency: Per Group") +
  xlab("Percent") +
  ylab(NULL) +
  scale_fill_manual(values = group_cols, name = "Group") +
  scale_x_continuous(labels = scales::percent)

tiff(
  "isotype_frequency.tiff",
  width = 6,
  height = 6,
  units = "in",
  compression = "lzw",
  res = 300,
  bg = "white"
)
if_2 + theme_bw(base_size = 20) + theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank())
dev.off()

#####################################################
#Clonal proportions - immunarch:
imm_top <-
  repClonality(immdata$data,
               .method = "top",
               .head = c(10, 100, 1000, 3000, 10000))

imm_top.plt <- vis(imm_top, .by = "Group", .meta = immdata$meta) +
  scale_fill_manual(values = group_cols) + xlab(NULL) &
  theme_bw(base_size = 20) +
  theme(axis.text.x =
          element_text(
            angle = 45,
            vjust = 1,
            hjust = 1
          )) &
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
imm_top.data <- imm_top.plt$data
imm_top.data$Group <- factor(imm_top.data$Group)
imm_top.plt$data %>% ggplot(aes(x = Grouping.var, y = Value, fill = Group)) + theme_classic(base_size = 20)+ 
  scale_fill_manual(values = group_cols, name = "Group") +
geom_bar(stat="identity", position = "dodge")

tiff(
  "rep_clonality.tiff",
  width = 8,
  height = 6,
  bg = "white",
  compression = "lzw",
  res = 300,
  units = "in"
)
imm_top.plt$data %>% ggplot(aes(x = Grouping.var, y = Value, fill = Group)) + theme_classic(base_size = 20)+ 
  scale_fill_manual(values = group_cols, name = "Group") +
  geom_bar(stat="identity", position = "dodge")
dev.off()

#####################################################
# Diversity Metrics - immunarch and vegan:
# Chao1 richness measure:
div_chao <- repDiversity(immdata$data, "chao1")

# Hill numbers:
div_hill <- repDiversity(immdata$data, "hill")

# D50:
div_d50 <- repDiversity(immdata$data, "d50")

# Ecological diversity measure, ie "true diversity":
div_div <- repDiversity(immdata$data, "div")

# Gini coefficient:
ginic = repDiversity(immdata$data, "gini")
ginic = data.frame(Gini = ginic[,1], Sample = row.names(ginic))
ginic$Group <- imm.meta[ginic$Sample, "Group"]
comparisons <- combn(unique(ginic$Group), 2)
p_df <- compare_means(Gini ~ Group, ginic, comparisons = comparisons, p.adjust.method = "holm")
y_max <- max(ginic$Gini)
p.value.y.coord <- rep(y_max, nrow(p_df))
step.increase <- (1:nrow(p_df)) * (y_max / 10)
p.value.y.coord <- p.value.y.coord + step.increase

p_df <- p_df %>%
  mutate(
    y.coord = p.value.y.coord,
    p.adj = format.pval(p.adj, digits = 1)
  )
ginic_plot <-  ginic %>% ggplot(aes(x = Group, y = Gini, fill = Group)) +
  geom_boxplot() +
  theme_bw(base_size = 20) + scale_fill_manual(values = group_cols) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) +
  ggtitle("Gini Coefficient")
ginic_plot <- ginic_plot + geom_signif(
  data = p_df,
  aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y.coord),
  manual = TRUE, tip_length = 0.03, size = .5, inherit.aes = FALSE
)
ginic_plot

# Gini-Simpson index:
gini_simp_plot <- repDiversity(
  immdata$data,
  .method = "gini.simp",
  .col = "aa") %>% vis(.by = c("Group"), .meta = immdata$meta) +
  scale_fill_manual(values = group_cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
gini_simp_plot 

# Inverse Simpson index:
inv_simpson_plt <- repDiversity(
  immdata$data,
  .method = "inv.simp",
  .col = "aa") %>% vis(.by = c("Group"), .meta = immdata$meta) +
  scale_fill_manual(values = group_cols) + 
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
inv_simpson_plt 

## Repertoire Overlap:
# Morisita-Horn:
mhorn_plt <- repOverlap(
  data_split,
  .method = "morisita",
  .col = "aa") %>% vis() + ggtitle("Morisita-Horn Similarity") + xlab(NULL) + ylab(NULL) +
  theme_bw(base_size = 20) + xlab(NULL) + ylab(NULL) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
mhorn_plt

# Jaccard:
jaccard_plt <- repOverlap(
  data_split,
  .method = "jaccard",
  .col = "aa") %>% vis() + ggtitle("Jaccard Similarity") +
  theme_bw(base_size = 20) + xlab(NULL) + ylab(NULL) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
jaccard_plt

# Cosine Similarity:
cosine_plt <- repOverlap(
  data_split,
  .method = "cosine",
  .col = "aa") %>% vis() + ggtitle("Cosine Similarity") +
  theme_bw(base_size = 20) + xlab(NULL) + ylab(NULL) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
cosine_plt

# Repertoire Clonality - immunarch:
imm_hom <- repClonality(immdata$data,
                        .method = "homeo",
                        .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1) 
)  %>% vis(.by = c("Group"), .meta = immdata$meta) + 
  scale_fill_manual(values = group_cols) + 
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
imm_hom

# Clones per Kiloread (CPK):
cpk_stats <- lapply(immdata$data, clones_per_kilo) %>% unlist() %>% data.frame()
cpk_stats$Group <- imm.meta[rownames(cpk_stats),"Group"] %>% factor(.)
cpk_stats$Sample <- rownames(cpk_stats)
colnames(cpk_stats) <- c("CPK", "Group", "Sample")

comparisons <- combn(unique(cpk_stats$Group), 2)

p_df <- compare_means(CPK ~ Group, cpk_stats, comparisons = comparisons, p.adjust.method = "holm")

y_max <- max(cpk_stats$CPK)
p.value.y.coord <- rep(y_max, nrow(p_df))
step.increase <- (1:nrow(p_df)) * (y_max / 10)
p.value.y.coord <- p.value.y.coord + step.increase

p_df <- p_df %>%
  mutate(
    y.coord = p.value.y.coord,
    p.adj = format.pval(p.adj, digits = 1)
  )
cpk_plt <- cpk_stats %>% ggplot(aes(x = Group, y = CPK, fill = Group)) +
  geom_boxplot() +
  theme_bw(base_size = 20) + scale_fill_manual(values = group_cols) +
   
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) +
  ggtitle("Clones per Kiloread")
cpk_plt
cpk_plt <- cpk_plt + geom_signif(
  data = p_df,
  aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y.coord),
  manual = TRUE, tip_length = 0.03, size = .5, inherit.aes = FALSE
)
cpk_plt

# Pielou evenness- vegan package
Pielou <- function(.data){
  mtx <- .data[,c("Clone.ID","Clones")] %>% pivot_wider(names_from = "Clone.ID", values_from = "Clones", values_fn = sum)
  H <- diversity(mtx)
  J <- H/log(specnumber(mtx))
}
Pielou_index <- lapply(immdata$data, Pielou) %>% unlist %>% data.frame()
Pielou_index$Sample <- rownames(Pielou_index)
Pielou_index$Group  <- imm.meta[rownames(Pielou_index), "Group"] 
colnames(Pielou_index) <- c("Pielou","Sample","Group")
comparisons <- combn(unique(Pielou_index$Group), 2)
p_df <- compare_means(Pielou ~ Group, Pielou_index, comparisons = comparisons, p.adjust.method = "holm")
y_max <- max(Pielou_index$Pielou)
p.value.y.coord <- rep(y_max, nrow(p_df))
step.increase <- (1:nrow(p_df)) * (y_max / 10)
p.value.y.coord <- p.value.y.coord + step.increase

p_df <- p_df %>%
  mutate(
    y.coord = p.value.y.coord,
    p.adj = format.pval(p.adj, digits = 1)
  )

Pielou_plt <- Pielou_index %>% ggplot(aes(x = Group, y = Pielou, fill = Group)) +
  geom_boxplot() +
  theme_bw(base_size = 20) + scale_fill_manual(values = group_cols) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Pielou Evenness") + geom_signif(
  data = p_df,
  aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y.coord),
  manual = TRUE, tip_length = 0.03, size = .5, inherit.aes = FALSE
)
Pielou_plt

# ACE index - vegan package:
ACE <- function(.data) {
  mtx <-
    .data[, c("Clone.ID", "Clones")] %>% pivot_wider(names_from = "Clone.ID", values_from = "Clones", values_fn = sum)
  ACE.idx <- estimateR(round(mtx)) %>% t() %>% data.frame()
  return(ACE.idx)
}
ACE_index <- lapply(immdata$data, ACE) %>% do.call(rbind, .)
ACE_index$Sample <- rownames(ACE_index)
ACE_index$Group  <- imm.meta[rownames(ACE_index), "Group"]

ACE_plt <- ACE_index %>% ggplot(aes(x = Group, y = S.ACE, fill = Group)) +
  geom_boxplot() +
  theme_bw(base_size = 20) + scale_fill_manual(values = group_cols) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) +
  ylab("ACE") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("ACE Index") 

comparisons <- combn(unique(ACE_index$Group), 2)
p_df <- compare_means(S.ACE ~ Group, ACE_index, comparisons = comparisons, p.adjust.method = "holm")
y_max <- max(ACE_index$S.ACE)
p.value.y.coord <- rep(y_max, nrow(p_df))
step.increase <- (1:nrow(p_df)) * (y_max / 10)
p.value.y.coord <- p.value.y.coord + step.increase

p_df <- p_df %>%
  mutate(
    y.coord = p.value.y.coord,
    p.adj = format.pval(p.adj, digits = 1)
  )

ACE_plt <- ACE_plt + geom_signif(
  data = p_df,
  aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y.coord),
  manual = TRUE, tip_length = 0.03, size = .25, inherit.aes = FALSE
)
ACE_plt

# Shannon index - vegan package:
Shannon <- function(.data){
  mtx <- .data[,c("Clone.ID","Clones")] %>% pivot_wider(names_from = "Clone.ID", values_from = "Clones", values_fn = sum)
  H <- diversity(mtx)
  return(H)
}
Shannon_index <- lapply(immdata$data, Shannon) %>% unlist %>% data.frame

Shannon_index$Sample <- rownames(Shannon_index)
Shannon_index$Group  <- imm.meta[rownames(Shannon_index), "Group"] 
colnames(Shannon_index) <- c("Shannon","Sample","Group")
Shannon_plt <- Shannon_index %>% ggplot(aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot() +
  theme_bw(base_size = 20) + scale_fill_manual(values = group_cols) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Shannon Index")
Shannon_plt

comparisons <- combn(unique(Shannon_index$Group), 2)
p_df <- compare_means(Shannon ~ Group, Shannon_index, comparisons = comparisons, p.adjust.method = "holm")
y_max <- max(Shannon_index$Shannon)
p.value.y.coord <- rep(y_max, nrow(p_df))
step.increase <- (1:nrow(p_df)) * (y_max / 10)
p.value.y.coord <- p.value.y.coord + step.increase

p_df <- p_df %>%
  mutate(
    y.coord = p.value.y.coord,
    p.adj = format.pval(p.adj, digits = 1)
  )

Shannon_plt <- Shannon_plt + geom_signif(
  data = p_df,
  aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y.coord),
  manual = TRUE, tip_length = 0.03, size = .5, inherit.aes = FALSE
)
Shannon_plt

# Plotting calculated Immunarch diversity metrics:
chao_plt <-
  vis(div_chao, .by = "Group", .meta = immdata$meta) + scale_fill_manual(values = group_cols) + 
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
chao_plt

hill_plt <-
  vis(div_hill, .by = "Group", .meta = immdata$meta) + scale_color_manual(values = group_cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL)  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
hill_plt

d50_plt <-
  vis(div_d50, .by = "Group", .meta = immdata$meta) + scale_fill_manual(values = group_cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
d50_plt

eco_div_plt <-
  vis(div_div, .by = "Group", .meta = immdata$meta) + scale_fill_manual(values = group_cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
eco_div_plt

#Save all of the diversity metric plots:
tiff("./diversity/chao.tiff", width = 8, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
print(chao_plt)
dev.off()

tiff("./diversity/hill.tiff", width = 8, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
hill_plt
dev.off()

tiff("./diversity/d50.tiff", width = 8, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
d50_plt
dev.off()

tiff("./diversity/true_div.tiff", width = 8, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
eco_div_plt
dev.off()

tiff("./diversity/gini_coef.tiff", width = 8, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
ginic_plot 
dev.off()

tiff("./diversity/shannon.tiff", width = 8, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
Shannon_plt
dev.off()

tiff("./diversity/ACE.tiff", width = 8, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
ACE_plt
dev.off()

tiff("./diversity/pielou.tiff", width = 8, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
Pielou_plt
dev.off()

tiff("./diversity/cpk.tiff", width = 8, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
cpk_plt
dev.off()

tiff("./diversity/cosine.tiff", width = 8, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
cosine_plt
dev.off()

tiff("./diversity/jaccard.tiff", width = 8, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
jaccard_plt
dev.off()

tiff("./diversity/morisita_horn.tiff", width = 8, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
mhorn_plt
dev.off()

tiff("./diversity/inverse_simpson.tiff", width = 8.5, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
inv_simpson_plt
dev.off()

tiff("./diversity/gini_simpson.tiff", width = 8, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
gini_simp_plot
dev.off()

tiff("./diversity/clonal_hom.tiff", width = 9, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
imm_hom
dev.off()

#####################################################
# Convergent clustering / network analysis of CDR3 sequences:
## Step 1: Group clones (regardless of Group) based on V/J genes, CDR3.aa length
## Step 2: Calculate Levenshtein distance between grouped CDR3 sequences
network_dist <-
  seqDist(
    list("network"=data),
    .col = 'CDR3.aa',
    .trim_genes = T,
    .group_by_seqLength = T,
    .method = 'lv',
    .group_by = c('V.name', 'J.name')
  )

# Step 3: Cluster the grouped CDR3 sequences
network_clust <- seqCluster(list("network"=data), network_dist, .perc_similarity = arguments$threshold)
network <- network_clust$network
clusters <- table(network$Cluster)

# Step 4: select all of the large clusters (>= 10 clones):
clusters <- clusters[clusters >= 10]
clusters <- clusters[order(clusters,decreasing = T)]
network <- network[network$Cluster %in% names(clusters),]

# Calculate cluster sizes and get the top 100 clusters
top_clusters <- network %>%
  group_by(Cluster) %>%
  tally(sort = TRUE) %>%
  top_n(50, n) %>%
  pull(Cluster)

# Filter the original data to include only the top 100 clusters
network_top50 <- network %>%
  filter(Cluster %in% top_clusters)

# Plot the top 50 clusters (bar plot)
network_clusters <- network_top50 %>%
  ggplot(aes(
    x = reorder(Cluster, Cluster, function(x) -length(x)),
    fill = group_id
  )) + 
  geom_bar() +
  scale_fill_manual(values = group_cols, name = "Group") + 
  xlab("Cluster") + 
  ylab("# Clones") + 
  theme_classic2(base_size = 20) +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  ggtitle("Global Clustering of Clonotypes", subtitle = "CDR3aa Levenshtein Distance, Sequence Similarity >= 0.7")

tiff(
  "./network_clusters.tiff",
  width = 12,
  height = 10,
  units = "in",
  bg = "white",
  compression = "lzw",
  res = 300
)
network_clusters
dev.off()

# Step 5: Networks analysis of large clusters:
# Split the data frame into groups based on the "Cluster" column
cluster_groups <- network_top50 %>%
  group_split(Cluster)

# Initialize an empty list to store network plots
network_plots <- list()

# Create a network plot for each cluster
for (i in seq_along(cluster_groups)) {
  # Get the current cluster data
  cluster_data <- cluster_groups[[i]]
  
  # Check the structure of cluster_data
  print(str(cluster_data))
  
  # Ensure that CDR3.aa and group_id columns are present and non-empty
  if ("CDR3.aa" %in% colnames(cluster_data) && "group_id" %in% colnames(cluster_data) &&
      nrow(cluster_data) > 0) {
    
    # Calculate the Levenshtein distance matrix for the current cluster
    dist_matrix <- calculate_distance_matrix(cluster_data) 
    
    # Create a graph from the distance matrix
    g <- graph_from_adjacency_matrix(dist_matrix, mode = "min", weighted = TRUE)
    
    # Assign colors to nodes based on group membership
    node_colors <- group_cols[cluster_data$group_id]
    
    # Plot k-mer motif (assuming CDR3.aa contains sequences and getKmers function works correctly)
    tryCatch({
      # Code that might produce an error
      kmer_plot <-
        getKmers(cluster_data, .k = nchar(cluster_data$CDR3.aa[1])) %>%
        kmer_profile(.method = "prob") %>%
        vis(.plot = "seq") +
        theme_bw(base_size = 20) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_y_continuous(labels = scales::percent) + 
        ggtitle(unique(cluster_data$Cluster))
    }, error = function(e) {
      # Do nothing on error
    })
    
    # Plot the graph with colored nodes
    tiff(paste0("./network/Cluster_", gsub("[[:punct:]]", "_", unique(cluster_data$Cluster)), ".tiff"), 
         width = 8, height = 8, units = "in", res = 300, compression = "lzw")
    plot(g, main = paste("Cluster", unique(cluster_data$Cluster)), vertex.color = node_colors, 
         edge.width = E(g)$weight, vertex.label = NA, cex.main = 20)
    dev.off()
    
    # Plot the motif for the cluster
    tiff(paste0("./network/Cluster_", gsub("[[:punct:]]", "_", unique(cluster_data$Cluster)), "_kmers.tiff"), 
         width = 10, height = 6, units = "in", res = 300, compression = "lzw")
    if (exists("kmer_plot")) print(kmer_plot)
    dev.off()
    
    # Add the graph to the list, named by the cluster
    network_plots[[paste("Cluster", unique(cluster_data$Cluster))]] <- g
  } else {
    warning(paste("Cluster data for", unique(cluster_data$Cluster), "is missing 'CDR3.aa' or 'group_id' columns, or the cluster is empty."))
  }
}

#####################################################

#Class-Switch Recombination (CSR) w/ BrepPhylo:
# bind samples into one data frame:
data <- do.call(rbind, immdata$data)
data$group_id <- imm.meta[data$sample_id, "Group"]
data$group_id <- factor(data$group_id)

# genes w/o mutations, just convert to NA:
data[data == "region_not_covered"] <- NA

# reshape data to include a column of Vgene name and allele for CSR calculation:
data$CloneID <- data$Clone.ID
data$CloneID <- paste(data$sample_id, data$CloneID, sep = "_")
data$vgene_allele <- str_extract(data$V.name, ".+?(?=\\x28)")
data$cgene <- str_extract(data$C.name, ".+?(?=\\*)")
data$vgene <- str_extract(data$V.name, ".+?(?=\\*)")
data$dgene <- str_extract(data$D.name, ".+?(?=\\*)")
data$jgene <- str_extract(data$J.name, ".+?(?=\\*)")

#split data by group:
data_split <- group_split(data, sample_id, .keep = TRUE)
names(data_split)  <- data$sample_id %>% unique() %>% sort()

distBCR <-
  seqDist(
    data_split,
    .col = 'CDR3.nt',
    group_by_seqLength = FALSE,
    .method = 'lv',
    .group_by = c('vgene_allele')
  )

#clustering BCR by CDR3 regions (to summarize across groups for CSR calculation)
clustBCR <- seqCluster(data_split, distBCR, .perc_similarity = 0.7)

clustBCR <- lapply(clustBCR, function(x){ 
  #get only the aligned V-gene sequence for CSR calc:
  x$V.seq <- substr(x$Sequence, 1, x$V.end)
  return(x)})

# Run CSR calculation
batch_results_all <- lapply(clustBCR, CSR_calculation)
batch_results_all <- batch_results_all[unlist(lapply(batch_results_all, function(x) length(unlist(x)) != 0))]
batch_summaries_all <- lapply(batch_results_all, function(x){summary <- getSummaryFromBatch(x)
return(summary$csr_events)})
batch_summaries_all <- do.call(rbind, batch_summaries_all)
batch_summaries_all$Group <- str_extract(batch_summaries_all$CloneID, ".+?(?=_)")
csr_summary <- summariseCSR( batch_summaries_all, 
                             dist_column = "distFromGermline", 
                             cloneID_column = "CloneID",
                             summarise_variables = "Group")

# next we can plot a graphical summary
# first set the order of the isotypes so that they reflect
# the actual physical order on the genome (hence the 
# order possible in CSR)
csr_summary$startIsotype <- factor(csr_summary$startIsotype,
                                   levels = c("M", "D", "G3", "G1", "A1", "G2",
                                              "G4", "E", "A2"),
                                   labels = c("M/D", "M/D", "G3", "G1", "A1", "G2",
                                              "G4", "E", "A2"))
csr_summary$endIsotype <- factor(csr_summary$endIsotype,
                                 levels = c("M", "D", "G3", "G1", "A1", "G2",
                                            "G4", "E", "A2"),
                                 labels = c("M/D", "M/D", "G3", "G1", "A1", "G2",
                                            "G4", "E", "A2"))

# plot it; it uses ggplot2 so you can extend this with ggplot2 functions
# here separate into different panels using the column 'Group'
library(ggplot2)
CSR_plot_per_pt <- plotCSRsummary(csr_summary) + facet_wrap(~ Group) & scale_fill_viridis_b()

# save the CSR plot:
tiff(
  "CSR_events_per_pt.tiff",
  width = 12,
  height = 8,
  units = "in",
  bg = "white",
  res = 300,
  compression = "lzw"
)
CSR_plot_per_pt + theme_bw(base_size = 20) + theme(panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank())
dev.off()

#####################################################
# save the environment!
save.image("downstream.RData")

