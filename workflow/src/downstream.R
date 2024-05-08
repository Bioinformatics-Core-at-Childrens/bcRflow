library(ComplexHeatmap)
library(viridisLite)
library(colorRamp2)
library(immunarch)
library(BrepPhylo)
library(Platypus)
library(optparse)
library(circlize)
library(ggplot2)
library(stringr)
library(ggpubr)
library(dplyr)
library(tidyr)
library(vegan)
library(NAIR)
source("utils.R")

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

#read in user input:
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
    #points to symLink of metadata table, file must be named samplesList.csv!
    help = "Path to sample metadata table [default %default]"
  ),
  make_option(
    c("-o", "--output_dir"),
    type = "character",
    default = file.path("./immunarch"),
    help = "Output directory [default %default]"
  )
)
arguments <- parse_args(OptionParser(option_list = option_list))

#argument handling:
species <- arguments$species
valid_species <- c("hsa", "mmu")
stopifnot(species %in% valid_species)

outdir <- file.path(arguments$output_dir)
metadata <- read.csv(arguments$metadata)
stopifnot(nrow(metadata) > 0)

#check if there's group metadata:
if (is.null(metadata$Group) |
    length(metadata$Group) != nrow(metadata)) {
  stop("Group metadata is missing or incomplete! Please check, and try again.")
}
rownames(metadata) <- metadata$SampleID

#check and create output directory:
if (!dir.exists(outdir)) {
  dir.create(outdir, showWarnings = T)
} else {
  print("Output directory already exists...")
}
setwd(outdir)

#reformatted dir for reshaped MiXCR output:
if (!dir.exists(file.path("./reformatted"))){
  dir.create(file.path("./reformatted"), showWarnings = T)
} else {
  print("Reformatted directory already exists...")
}

#diversity dir for diversity metric figures:
if (!dir.exists(file.path("./diversity"))){
  dir.create(file.path("./diversity"), showWarnings = T)
} else {
  print("Diversity directory already exists...")
}

#chord dir for V-J usage chord plots:
if (!dir.exists(file.path("./chords"))){
  dir.create(file.path("./chords"), showWarnings = T)
} else {
  print("Chords directory already exists...")
}

#list MiXCR output reports:
IGH_reports <-
  list.files(
    file.path(arguments$reports),
    full.names = T,
    pattern = ".tsv",
    recursive = T
  )

#load data using Immunarch:
immdata <- repLoad(file.path(IGH_reports))
imm.meta <- metadata[, c("SampleID", "Group")]
colnames(imm.meta) <-  c("Sample", "Group")
immdata$meta <- tibble(imm.meta)

#Preprocessing MiXCR output:
# 1) add sample ID to data
# 2) get only the best gene annotation for V, D, J and C(first in ranked list)
# 3) add Somatic Hypermutations for V, D, and J genes

for (i in 1:length(immdata$data)) {
  immdata$data[[i]]$sample_id <- names(immdata$data)[i]
  immdata$data[[i]]$V.name <-
    str_extract(immdata$data[[i]]$V.name, ".+?(?=\\,)")
  immdata$data[[i]]$D.name <-
    str_extract(immdata$data[[i]]$D.name, ".+?(?=\\,)")
  immdata$data[[i]]$J.name <-
    str_extract(immdata$data[[i]]$J.name, ".+?(?=\\,)")
  immdata$data[[i]]$C.name <-
    str_extract(immdata$data[[i]]$C.name, ".+?(?=\\,)")
  
  #add mutations from MiXCR output (for SHM)
  tmp <-
    fread(
      IGH_reports[[i]],
      select = c(
        "cloneId",
        "nMutationsVRegion",
        "nMutationsJRegion",
        "nMutationsDRegion"
      )
    )
  immdata$data[[i]] <-
    merge(immdata$data[[i]], tmp, by.x = "Clone.ID", by.y = "cloneId")
}

#bind samples into one data frame:
data <- do.call(rbind, immdata$data)
data$group_id <- imm.meta[data$sample_id, "Group"]

#genes w/o mutations, just convert to NA:
data[data == "region_not_covered"] <- NA

#reshape for CSR calc:
data$CloneID <- data$Clone.ID
data$CloneID <- paste(data$sample_id, data$CloneID, sep = "_")
data$vgene_allele <- str_extract(data$V.name, ".+?(?=\\x28)")
data$cgene <- str_extract(data$C.name, ".+?(?=\\*)")
data$vgene <- str_extract(data$V.name, ".+?(?=\\*)")
data$dgene <- str_extract(data$D.name, ".+?(?=\\*)")
data$jgene <- str_extract(data$J.name, ".+?(?=\\*)")

#split data by group:
data_split <- group_split(data, group_id, .keep = TRUE)
names(data_split)  <- data$group_id %>% unique() %>% sort()

#####################################################
#IGH-V Gene Usage:
species <- str_sub(arguments$species, end = -2)
ighv_gu <- geneUsage(immdata$data, paste0(species, ".ighv"))
ighv_gu <- ighv_gu %>% replace(is.na(.), 0) %>% data.frame()
rownames(ighv_gu) <- ighv_gu$Names
ighv_gu$Names <- NULL

##ComplexHeatmap:
#annotation df:
ighv_annotation_df <-
  data.frame("Sample" = colnames(ighv_gu), "Group" = imm.meta[colnames(ighv_gu), "Group"])
ighv_annotation_df <-
  ighv_annotation_df[order(ighv_annotation_df$Group),]

#sample colors:
sample_cols <-
  custom_colors$discrete[1:length(unique(ighv_annotation_df$Sample))]
names(sample_cols) <- ighv_annotation_df$Sample %>% unique()

#group colors:
group_cols <-
  colors_dutch[1:length(unique(ighv_annotation_df$Group))]
names(group_cols) <- ighv_annotation_df$Group %>% unique()

#drop rows with no mapped IGH-V gene
ighv_gu <- ighv_gu[row.names(ighv_gu) != "NA", , drop = FALSE]

#scale gene util counts:
ighv_gu <- scale(ighv_gu, center = T, scale = T)
col_fun = colorRamp2(c(min(ighv_gu), 0 , max(ighv_gu)), c("blue", "white", "red"))

ighv_annotation_df <-
  ighv_annotation_df[order(ighv_annotation_df$Group), ]
ighv_gu <- ighv_gu[, ighv_annotation_df$Sample]

#plot with ComplexHeatmap
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
    df = ighv_annotation_df,
    col = list(Sample = sample_cols, Group = group_cols),
    show_annotation_name = TRUE,
    show_legend = c(FALSE, TRUE),
    annotation_name_gp = gpar(fontsize = 20),
    annotation_legend_param = list(fontsize = 20, labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20, fontface = "bold"))
  ),
  column_labels = NULL,
  #legend_label_dp = list(title_gp = gpar(fontsize = 20, fontface = "bold"),labels_gp = gpar(fontsize = 20)),
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

#####################################################
#V-J Gene Usage:
#generate the v-j counts matrix for each sample
v_j_matrix <- function(.data, threshold = 3){
  #extract clones with both V and J-gene alignments
  v_j_data <- .data[!is.na(.data$vgene) & !is.na(.data$jgene),]
  
  #tally the unique v-j pairs for each sample
  v_j_counts <- v_j_data %>% count(sample_id, vgene, jgene)
  
  #apply the lower threshold
  v_j_counts <- v_j_counts[v_j_counts$n >= threshold,]
  v_j_samples <- unique(v_j_counts$sample_id)
  v_j_counts <- v_j_counts %>% group_split(sample_id, .keep = F) 
  names(v_j_counts) <- v_j_samples %>% unique()
  
  v_j_counts <- lapply(v_j_counts, function(x){
    x <- x %>% pivot_wider(names_from = vgene, values_from = n,values_fill = 0) %>% data.frame
    rownames(x) <- x$jgene
    x$jgene <- NULL
    x <- as.matrix(x)
    colnames(x) <- colnames(x) %>% str_replace_all("\\.",replacement = "\\-")
    return(x)
  })
  
  #return named list of V-J counts >= thresh
  return(v_j_counts)
}
v_j_counts <- v_j_matrix(data)

#plot the chords:
for(i in 1:length(v_j_counts)){
  tiff(paste0("./chords/",names(v_j_counts)[[i]],"_vj_chord.tiff"), width = 4.5, height = 4.5, bg = "white", units = "in", res =300, compression = "lzw")
  suppressMessages(vj_circos(v_j_counts[[i]]))
  dev.off()
}

# v_j_counts$sample_id <-
#   paste(v_j_counts$sample_id, imm.meta[v_j_counts$sample_id, "Group"], sep = ".")
# 
# v_j_counts <-
#   v_j_counts %>% pivot_wider(names_from = c(sample_id), values_from = n) %>% data.frame()
# 
# v_j_counts[is.na(v_j_counts)] <- 0
# rownames(v_j_counts) <- v_j_counts$v_j_genes
# 
# v_j_counts$v_j_genes <- NULL
# 
# colnames(v_j_counts) <-
#   sub(".*\\.", "", colnames(v_j_counts)) %>% make.unique
# 
# v_j_counts <-
#   v_j_counts[, str_sort(colnames(v_j_counts), numeric = TRUE)]
# 
# vj_corrplot <- v_j_counts %>% cor() %>% ggcorrplot::ggcorrplot(
#   method = "square",
#   show.diag = F,
#   hc.order = T,
#   ggtheme = theme_bw
# ) + ggtitle("Paired IGHV/J Usage Correlation", subtitle = "Samples ordered via hierarchical clustering")
# 
# tiff(
#   "IGHV_IGHJ_usage.tiff",
#   width = 0.25 * ncol(v_j_counts) + 1.5,
#   height = 0.25 * ncol(v_j_counts) + 1.5,
#   units = "in",
#   compression = "lzw",
#   res = 300,
#   bg = "white"
# )
# vj_corrplot  + theme_bw(base_size = 22) + xlab(NULL) + ylab(NULL) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# dev.off()

#####################################################
#Somatic Hyper Mutation (SHM):
#mutation count per clone for SHM - sum all listed mutations for V,D,J regions output by MiXCR:
data$SHM <-
  sapply(
    strsplit(data$nMutationsVRegion, ","),
    FUN = function(x) {
      length(x[!is.na(x)])
    }
  ) + sapply(
    strsplit(data$nMutationsDRegion, ","),
    FUN = function(x) {
      length(x[!is.na(x)])
    }
  ) + sapply(
    strsplit(data$nMutationsJRegion, ","),
    FUN = function(x) {
      length(x[!is.na(x)])
    }
  )
data$SHM.rate <- data$SHM / str_length(data$CDR3.nt)

shm_rate <-
  data.table(
    "Group"     = data$group_id,
    "Sample"    = data$sample_id,
    "C.name"    = data$C.name,
    "CDR3.nt"   = data$CDR3.nt,
    "SHM.rate"  = data$SHM.rate,
    "SHM.count" = data$SHM,
    "CDR3.length" = str_length(data$CDR3.nt)
  )
shm_rate$C.name <- str_extract(shm_rate$C.name, ".+?(?=\\*)")
shm_rate$C.name[is.na(shm_rate$C.name)] <- "Unknown"
shm_rate <- shm_rate[complete.cases(shm_rate),]

#plot SHM rate across all groups & C-gene isotypes:
shm.plot <-
  shm_rate %>% ggplot(aes(x = Group, y = SHM.rate, color = Group)) +
  ylab("Rate") +
  xlab(NULL) +
  geom_boxplot() + theme_bw(base_size = 20) + theme(panel.grid.major = element_blank(),
                                                    panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = group_cols, name = "Group") +
  facet_wrap(~ C.name, scales = "free") +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("Somatic Hypermutation Rates")
shm.plot

tiff(
  "SHM_rates.tiff",
  width = 12,
  height = 9,
  bg = "white",
  units = "in",
  compression = "lzw",
  res = 300
)
shm.plot + theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) &
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

#####################################################
#CDR3 Length Distribution:
cdr3.violin <-
  data %>% ggplot(aes(
    x = group_id,
    y = str_length(.$CDR3.aa),
    fill = group_id
  )) +
  geom_violin() + scale_fill_manual(values = group_cols, name = "Group") + ylab("CDR3 Length") + xlab(NULL)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("CDR3 Length Distribution")
cdr3.violin

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
  geom_violin() + scale_fill_manual(values = group_cols, name = "Group") + ylab("CDR3 Length") + xlab(NULL)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(panel.grid.major = element_blank(),
                                                                   panel.grid.minor = element_blank())
dev.off()

#logo plot of CDR3 region across groups:
#split data by group:
data_split <- group_split(data, group_id, .keep = TRUE)
names(data_split)  <- data$group_id %>% unique() %>% sort()

#calculate kmer probability based on average CDR3 length:
logoplots_fun <- function(x) {
  avg_cdr3_length <- mean(str_length(x$CDR3.aa)) %>% round
  x_avg <- x[str_length(x$CDR3.aa) == avg_cdr3_length,]
  plt <-
    getKmers(x_avg, .k = avg_cdr3_length) %>% kmer_profile(.method = "prob") %>%  vis(.plot = "seq")
  return(plt)
}
logoplots <- lapply(data_split, logoplots_fun)

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
if_2 <- data[!is.na(data$cgene), ]  |>
  mutate(cgene = factor(cgene),
         group_id = factor(group_id)) |>
  group_by(cgene, group_id) |>
  summarise(n = sum(round(Clones), na.rm = TRUE)) |>
  group_by(group_id) |>
  mutate(pct = prop.table(n)) |>
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

tiff(
  "rep_clonality.tiff",
  width = 8,
  height = 6,
  bg = "white",
  compression = "lzw",
  res = 300,
  units = "in"
)
vis(imm_top, .by = "Group", .meta = immdata$meta) +
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
dev.off()

#####################################################
# Diversity Metrics - immunarch and vegan:
# Chao1 diversity measure
div_chao <- repDiversity(immdata$data, "chao1")

# Hill numbers
div_hill <- repDiversity(immdata$data, "hill")

# D50
div_d50 <- repDiversity(immdata$data, "d50")

# Ecological diversity measure
div_div <- repDiversity(immdata$data, "div")

#Gini-Simpson index:
gini_simp_plot <- repDiversity(
  immdata$data,
  .method = "gini.simp",
  .col = "aa") %>% vis(.by = c("Group"), .meta = immdata$meta) +
  scale_fill_manual(values = group_cols) + 
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Inverse Simpson index:
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

#Repertoire Overlap:
#Morisita-Horn:
mhorn_plt <- repOverlap(
  data_split,
  .method = "morisita",
  .col = "aa") %>% vis() + ggtitle("Morisita-Horn Similarity") + xlab(NULL) + ylab(NULL) +
  theme_bw(base_size = 20) + xlab(NULL) + ylab(NULL) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Jaccard:
jaccard_plt <- repOverlap(
  data_split,
  .method = "jaccard",
  .col = "aa") %>% vis() + ggtitle("Jaccard Similarity") +
  theme_bw(base_size = 20) + xlab(NULL) + ylab(NULL) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Cosine Similarity:
cosine_plt <- repOverlap(
  data_split,
  .method = "cosine",
  .col = "aa") %>% vis() + ggtitle("Cosine Similarity") +
  theme_bw(base_size = 20) + xlab(NULL) + ylab(NULL) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Clonal Expansion Index (CEI):
# clonal_expansion_index <- function(.data){
#   raw <- .data[order(.data$Clones, decreasing = T),]
#   total_clones <- sum(.data$Clones)
#   r <- raw$Clones / total_clones
#   cei = 0
#   for(i in 1:length(r)){
#     cei = cei + ((r[i] - (i/total_clones)) / total_clones)
#   }
#   return(cei)
# }
# immdata$data %>% length()
# cei_stats <- lapply(immdata$data, clonal_expansion_index)

#Clonality - immunarch:
imm_hom <- repClonality(immdata$data,
                        .method = "homeo",
                        .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1) 
)  %>% vis(.by = c("Group"), .meta = immdata$meta) + 
  scale_fill_manual(values = group_cols) + 
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Clones per Kiloread (CPK) - Expansion Index:
clones_per_kilo <- function(.data){
  num_CDR3 <- nrow(.data)
  sum_of_counts = sum(.data$Clones)
  normalized_data = (num_CDR3 / sum_of_counts)
  return(normalized_data)
}
cpk_stats <- lapply(immdata$data, clones_per_kilo) %>% unlist() %>% data.frame()
cpk_stats$Group <- imm.meta[rownames(cpk_stats),"Group"]
cpk_stats$Sample <- rownames(cpk_stats)
colnames(cpk_stats) <- c("CPK", "Group", "Sample")

cpk_plt <- cpk_stats %>% ggplot(aes(x = Group, y = CPK, fill = Group)) +
  geom_boxplot() +
  theme_bw(base_size = 20) + scale_fill_manual(values = group_cols) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) +
  ggtitle("Clones per Kiloread")

#Pielou evenness- vegan package
Pielou <- function(.data){
  mtx <- .data[,c("Clone.ID","Clones")] %>% pivot_wider(names_from = "Clone.ID", values_from = "Clones")
  H <- diversity(mtx)
  J <- H/log(specnumber(mtx))
}
Pielou_index <- lapply(immdata$data, Pielou) %>% unlist %>% data.frame()
Pielou_index$Sample <- rownames(Pielou_index)
Pielou_index$Group  <- imm.meta[rownames(Pielou_index), "Group"] 
colnames(Pielou_index) <- c("Pielou","Sample","Group")
Pielou_plt <- Pielou_index %>% ggplot(aes(x = Group, y = Pielou, fill = Group)) +
  geom_boxplot() +
  theme_bw(base_size = 20) + scale_fill_manual(values = group_cols) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Pielou Evenness")
Pielou_plt

#ACE index - vegan package:
ACE <- function(.data) {
  mtx <-
    .data[, c("Clone.ID", "Clones")] %>% pivot_wider(names_from = "Clone.ID", values_from = "Clones")
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
ACE_plt

#Shannon index - vegan package:
Shannon <- function(.data){
  mtx <- .data[,c("Clone.ID","Clones")] %>% pivot_wider(names_from = "Clone.ID", values_from = "Clones")
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

#Gini coefficient - vegan package:
Gini <- function(.data){
  mtx <- .data[,c("Clone.ID","Clones")] %>% pivot_wider(names_from = "Clone.ID", values_from = "Clones")
  H <- diversity(mtx, index="simpson")
  return(H)
}
Gini_index <- lapply(immdata$data, Gini) %>% unlist %>% data.frame

Gini_index$Sample <- rownames(Gini_index)
Gini_index$Group  <- imm.meta[rownames(Gini_index), "Group"] 
colnames(Gini_index) <- c("Gini","Sample","Group")
Gini_plt <- Gini_index %>% ggplot(aes(x = Group, y = Gini, fill = Group)) +
  geom_boxplot() +
  theme_bw(base_size = 20) + scale_fill_manual(values = group_cols) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Gini Coefficient")
Gini_plt

#Immunarch diversity calculations:
chao_plt <-
  vis(div_chao, .by = "Group", .meta = immdata$meta) + scale_fill_manual(values = group_cols) + 
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

hill_plt <-
  vis(div_hill, .by = "Group", .meta = immdata$meta) + scale_color_manual(values = group_cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

d50_plt <-
  vis(div_d50, .by = "Group", .meta = immdata$meta) + scale_fill_manual(values = group_cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

eco_div_plt <-
  vis(div_div, .by = "Group", .meta = immdata$meta) + scale_fill_manual(values = group_cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Save the diversity metric plots:
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
Gini_plt 
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

tiff("./diversity/rep_clonality.tiff", width = 9, height = 6, units = "in", bg = "white", compression = "lzw", res = 300)
imm_hom
dev.off()

#####################################################
#NAIR: Clonal clustering/network analysis:
reports <-
  IGH_reports[which(lapply(IGH_reports, function(x) {
    length(count.fields(x, skip = 1)) %>% c()
  }) > 0)]
ids     <- basename(reports) %>% str_remove_all("\\.tsv")
group_labels   <- metadata[ids, "Group"]

# Identify clones in a neighborhood around each associated sequence
findPublicClusters(
  file_list = reports,
  input_type = "tsv",
  seq_col = "aaSeqImputedCDR3",
  sample_ids = ids,
  count_col = "readCount",
  top_n_clusters = 10,
  output_dir = "./NAIR",
  dist_cutoff = 0.5
)

# Directory of node metadata from step 1
dir_filtered_samples_node <-
  file.path("./NAIR", "node_meta_data")

# Vector of file paths to node metadata from step 1
files_filtered_samples_node <-
  list.files(dir_filtered_samples_node, full.names = TRUE)

#Colors for network plot
color_list <-
  custom_colors$discrete[1:length(ighv_annotation_df$Sample)]
names(color_list) <- ighv_annotation_df$Sample

#Build the network of public clonotypes based on the CDR3 AA seq:
public_clusters <-
  buildPublicClusterNetwork(
    files_filtered_samples_node,
    seq_col = "aaSeqImputedCDR3",
    count_col = "readCount",
    size_nodes_by = 1,
    print_plots = TRUE,
    color_legend = T
  )
#Label the top 30 clusters:
public_clusters <-
  labelClusters(
    public_clusters,
    top_n_clusters = 30,
    cluster_id_col = "ClusterIDPublic",
    size = 4
  )

network_cols <- group_cols[metadata[names(sample_cols), "Group"]]
names(network_cols) <- names(sample_cols)

#Plot the clusters:
network_plot <- public_clusters$plots$SampleID +
  scale_color_manual(values = network_cols, name = "Group") +
  ggtitle(
    "Global Network of Public Clusters",
    subtitle = "Each node denotes a single public BCR clone
Edges denote a maximum hamming distance of 1 between
receptor sequences

Nodes colored by Group"
  ) + theme_bw(base_size = 15) + theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
  ) + xlab(NULL) + ylab(NULL)

#Fetch the NAIR clusters:
public_clusters$node_data$cdr3aa <-
  public_clusters$node_data$aaSeqImputedCDR3

#save plot
tiff(
  "NAIR_global_clusters.tiff",
  width = 8,
  height = 8,
  units = "in",
  bg = "white",
  compression = "lzw",
  res = 300
)
network_plot
dev.off()

NAIR_clusters <-
  public_clusters$node_data %>% group_split(ClusterIDPublic)

#reshape top clusters to be read in to create Motif logos:
for (i in 1:length(NAIR_clusters)) {
  table = NAIR_clusters[[i]]
  write.table(
    table,
    file =
      paste0("./reformatted/NAIR_cluster_", i, ".tsv"),
    quote = F,
    col.names = T,
    row.names = F,
    sep = "\t"
  )
}
NAIR_imm <- repLoad("./reformatted", .coding = FALSE)

cluster1_plt <- getKmers(NAIR_imm$data$NAIR_cluster_1, .k = 13) %>% 
                  kmer_profile(.method = "prob") %>%  
                  vis(.plot = "seq") + 
                  theme_bw(base_size = 20) + 
                  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
                  scale_y_continuous(labels = scales::percent) + ggtitle("NAIR: Cluster 1") 

cluster2_plt <- getKmers(NAIR_imm$data$NAIR_cluster_2, .k = 13) %>% kmer_profile(.method = "prob")%>%  vis(.plot = "seq") +
                  theme_bw(base_size = 20) + 
                  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
                  scale_y_continuous(labels = scales::percent) + ggtitle("NAIR: Cluster 2")

cluster3_plt <- getKmers(NAIR_imm$data$NAIR_cluster_3, .k = 12)  %>% kmer_profile(.method = "prob")%>%  vis(.plot = "seq") +
                  theme_bw(base_size = 20) + 
                  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
                  scale_y_continuous(labels = scales::percent) + ggtitle("NAIR: Cluster 3")

cluster4_plt <- getKmers(NAIR_imm$data$NAIR_cluster_4, .k = 16) %>% kmer_profile(.method = "prob")%>%  vis(.plot = "seq") +
                  theme_bw(base_size = 20) + 
                  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
                  scale_y_continuous(labels = scales::percent) + ggtitle("NAIR: Cluster 4")

cluster5_plt <- getKmers(NAIR_imm$data$NAIR_cluster_5, .k = 14) %>% kmer_profile(.method = "prob")%>%  vis(.plot = "seq") +
                  theme_bw(base_size = 20) + 
                  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
                  scale_y_continuous(labels = scales::percent) + ggtitle("NAIR: Cluster 5")

pdf("NAIR_top5_cluster_motif.pdf", width = 8, height = 4, bg = "white")
cluster1_plt
cluster2_plt
cluster3_plt
cluster4_plt
cluster5_plt
dev.off()

#####################################################
#Class-Switch Recombination (CSR) w/ BrepPhylo:
data$CloneID <- data$Clone.ID
data$CloneID <- paste(data$sample_id, data$CloneID, sep = "_")
data$vgene_allele <- str_extract(data$V.name, ".+?(?=\\x28)")
data$cgene <- str_extract(data$C.name, ".+?(?=\\*)")
data$vgene <- str_extract(data$V.name, ".+?(?=\\*)")
data$dgene <- str_extract(data$D.name, ".+?(?=\\*)")
data$jgene <- str_extract(data$J.name, ".+?(?=\\*)")

#split data by group:
data_split <- group_split(data, group_id, .keep = TRUE)
names(data_split)  <- data$group_id %>% unique() %>% sort()

##Clustering BCR clonotypes within groups based on CDR3 AA sequence using levenstein distance:
#calculate distance
distBCR <-
  seqDist(
    data_split,
    .col = 'CDR3.aa',
    group_by_seqLength = FALSE,
    .trim_genes = T,
    .perc_similarity = 0.6,
    .method = 'lv',
    .group_by = c('V.name')
  )

#clustering BCR by CDR3 regions (to summarize across groups for CSR calculation)
clustBCR <- seqCluster(data_split, distBCR, .perc_similarity = 0.6)

#merge the clustered clonotypes into one dataframe:
clones <- do.call(rbind, clustBCR)

#add unique ID to prevent overlap of cluster ids:
clones[!is.na(clones$Cluster), "Cluster"] <-
  paste0(clones[!is.na(clones$Cluster), "group_id"], "_", clones[!is.na(clones$Cluster), "Cluster"])
clones <- clones[!is.na(clones$Cluster),]

#subset to desired columns only:
clones <-
  clones[, c(
    "Cluster",
    "CloneID",
    "sample_id",
    "group_id",
    "cgene",
    "vgene",
    "vgene_allele",
    "CDR3.nt"
  )]

#note: dnapars is available w/ the BrepPhylo package, but only works with Linux!
dnapars_executable <-
  system.file("exe/dnapars", package = "BrepPhylo")
outputFolder <- path.expand("./CSR_batchAnalysis")
dir.create(outputFolder, showWarnings = TRUE)

csr_species <- "Homo_sapiens"
if (species == "mmu") {
  csr_species <- "Mus_musculus"
}
clones <- clones[complete.cases(clones),]

clones$Subclass <- str_replace(clones$cgene, "IGH", "")
clones$Class <- clones$Subclass %>% gsub('[0-9.]', '', .)
clones <- clones %>% data.frame()

clones$Cluster <- clones$Cluster %>% str_replace_all("-", "_")

#clonal lineage analysis:
clones$Clone.ID <- clones$CloneID
clones$CloneID <- clones$Cluster
batch_results <- doBatchCloneAnalysis(
  clones,
  outputFolder = outputFolder,
  species = csr_species,
  sequence_column = "CDR3.nt",
  cloneID_column = "CloneID",
  IGHVgeneandallele_column = "vgene_allele",
  plotFormat = "pdf",
  label_column = "cgene",
  phyloTreeType = "dnapars",
  phyloTreeOptions = list("executable" = dnapars_executable),
  useTempDir = FALSE,
  minCloneSize = 3
)

#stop if BrepPhylo fails:
if (length(unlist(batch_results)) == 0) {
  save.image("downstream.RData")
  stop(
    "BrepPhylo was unable to perform clonal lineage analysis, likely due to a lack of clones.\nPlease see the saved .RData file in the output directory (downstream.RData)"
  )
}

#group CSR events by Group metadata:
batch_summary <- getSummaryFromBatch(batch_results)
batch_summary.csr <- batch_summary$csr_events
batch_summary.csr$Group <-
  str_extract(batch_summary.csr$CloneID, ".+?(?=\\_)")

#summarize CSR events:
csr_summary <- summariseCSR(
  batch_summary.csr,
  dist_column = "distFromGermline",
  cloneID_column = "CloneID",
  summarise_variables = "Group"
)

#plot CSR events:
CSR_plot <- plotCSRsummary(csr_summary) + facet_wrap( ~ Group) +
  scale_fill_viridis_c(name = "mean distance\nfrom germline") +
  ggtitle(label = "Class Switch Recombination", subtitle = "BrepPhylo")
CSR_plot + theme_bw(base_size = 20) + theme(panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank())

tiff(
  "CSR_events.tiff",
  width = 9,
  height = 8,
  units = "in",
  bg = "white",
  res = 300,
  compression = "lzw"
)
CSR_plot + theme_bw(base_size = 20) + theme(panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank())
dev.off()

#####################################################
#save the environment!
save.image("downstream.RData")
