# load required packages:
library(ComplexHeatmap)
library(viridisLite)
library(colorRamp2)
library(immunarch)
library(BrepPhylo)
library(Platypus)
library(magrittr)
library(purrr)
library(optparse)
library(circlize)
library(epitools)
library(ggplot2)
library(stringr)
library(igraph)
library(ggpubr)
library(dplyr)
library(tidyr)
library(vegan)
library(fpc)
library(stringdist)
library(NAIR)

# create custom discrete color palette:
custom_colors <- list()

colors_dutch <- c(
  '#FFC312',
  '#BF3EFF',
  '#12CBC4',
  '#FDA7DF',
  '#ED4C67',
  '#F79F1F',
  '#A3CB38',
  '#1289A7',
  '#D980FA',
  '#B53471',
  '#EE5A24',
  '#009432',
  '#0652DD',
  '#9980FA',
  '#833471',
  '#EA2027',
  '#006266',
  '#1B1464',
  '#5758BB',
  '#6F1E51'
)

colors_spanish <- c(
  '#40407a',
  '#706fd3',
  '#f7f1e3',
  '#34ace0',
  '#33d9b2',
  '#2c2c54',
  '#474787',
  '#aaa69d',
  '#227093',
  '#218c74',
  '#ff5252',
  '#ff793f',
  '#d1ccc0',
  '#ffb142',
  '#ffda79',
  '#b33939',
  '#cd6133',
  '#84817a',
  '#cc8e35',
  '#ccae62'
)
colors_blind <- c(
  "#FFBF7FFF",
  "#FF7F00FF",
  "#FFFF99FF",
  "#FFFF32FF",
  "#B2FF8CFF",
  "#32FF00FF",
  "#A5EDFFFF",
  "#19B2FFFF",
  "#CCBFFFFF",
  "#654CFFFF",
  "#FF99BFFF",
  "#E51932FF"
)
custom_colors$discrete <- c(colors_dutch, colors_spanish, colors_blind)


# Create a function to extract the first element of a comma-separated list
extract_first_hit <- function(string) {
  elements <- unlist(strsplit(string, ","))
  first_element <- trimws(elements[1]) # Trim whitespace
  return(first_element)
}

#### Draw Circos-Plot ####
vj_circos <- function(.data) {
  group <- c(rep("IGHJ", nrow(.data)), rep("IGHV", ncol(.data)))
  names(group) <- c(rownames(.data), colnames(.data))
  
  grid.col <-
    stats::setNames(grDevices::rainbow(length(union(
      rownames(.data), colnames(.data)
    ))), sample(union(rownames(.data), colnames(.data))))
  gene.label.length <-
    max(nchar(gsub("^\\D+", "", gsub(
      "TRAV", "", gsub("TRBV", "",  gsub("TRAJ", "", gsub(
        "TRBJ", "", gsub("IGKV", "", gsub("IGLV", "", gsub(
          "IGHV", "", gsub("IGKJ", "", gsub("IGLJ", "", gsub(
            "IGHJ", "", append(dimnames(.data)[[1]], dimnames(.data)[[2]])
          )))
        )))
      )))
    )))) # Determine length of longest gene label
  gene.label <- T
  label.threshold <- 0
  c.count.label = T
  c.count.label.size = 0.8
  gene.label.size = 0.8
  axis = "basic"
  
  circlize::circos.clear()
  circlize::circos.par(points.overflow.warning = FALSE)
  
  circlize::chordDiagram(
    .data,
    link.sort = T,
    link.decreasing = T,
    directional = 0,
    direction.type = c("arrows"),
    link.arr.length = 0.2,
    annotationTrack = c("grid"),
    annotationTrackHeight = c(0.02),
    preAllocateTracks = list(
      list(
        track.height = circlize::mm_h(0.2),
        track.margin = c(0.01, 0)
      ),
      # track 1: margin around plot
      list(
        track.height = circlize::mm_h(0.25),
        track.margin = c(0.01, 0)
      ),
      # track 2: line indicating groups
      list(
        track.height = circlize::mm_h(0.8 + gene.label.length * 0.6),
        track.margin = c(0.01, 0)
      ),
      # track 3: gene label
      list(
        track.height = circlize::mm_h(3.5),
        track.margin = c(0.01, 0)
      )
    ),
    # track 4: clone count
    grid.col = grid.col
  )
  # circlize::circos.info()
  
  
  #### Add Gene/Cluster labels to circos plot ####
  if (gene.label) {
    circlize::circos.track(
      track.index = 4,
      panel.fun = function(x, y) {
        if (circlize::get.cell.meta.data("xrange") > label.threshold) {
          ycenter <- circlize::get.cell.meta.data("ycenter")
          xcenter <- circlize::get.cell.meta.data("xcenter")
          sector.index <- circlize::get.cell.meta.data("sector.index")
          
          sector.index <- gsub("TRAV", "", sector.index)
          sector.index <- gsub("TRBV", "", sector.index)
          sector.index <- gsub("TRAJ", "", sector.index)
          sector.index <- gsub("TRBJ", "", sector.index)
          sector.index <- gsub("IGKV", "", sector.index)
          sector.index <- gsub("IGLV", "", sector.index)
          sector.index <- gsub("IGHV", "", sector.index)
          sector.index <- gsub("IGKJ", "", sector.index)
          sector.index <- gsub("IGLJ", "", sector.index)
          sector.index <- gsub("IGHJ", "", sector.index)
          sector.index <- gsub("^\\D+", "", sector.index)
          
          if (gene.label.size == "undef") {
            if (nchar(sector.index) > 2) {
              label.cex <- 0.6
            }
            if (nchar(sector.index) > 6) {
              label.cex <- 0.4
            }
            if (nchar(sector.index) <= 2) {
              label.cex <- 0.7
            }
          } else{
            label.cex <- gene.label.size
          }
          
          circlize::circos.text(
            xcenter,
            ycenter - circlize::mm_h(7 + gene.label.length * 0.1),
            sector.index,
            facing = "bending.inside",
            niceFacing = TRUE,
            col = "black",
            font = 2,
            cex = label.cex
          )
        } else{
          # ylim <- circlize::get.cell.meta.data("ylim")
          # xcenter <- circlize::get.cell.meta.data("xcenter")
          # sector.index <- circlize::get.cell.meta.data("sector.index")
          # circlize::circos.text(xcenter, ylim[1], sector.index,
          #                       facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.6), cex=0.5, col="white")
        }
      },
      bg.border = NA
    ) # here set bg.border to NA is important
    
  }
  
  #### Add Group labels to circos plot ####
  
  for (p in unique(unname(group))) {
    sub_group = names(group[group == p])
    print(sub_group)
  }
  
  
  for (p in unique(unname(group))) {
    sub_group = names(group[group == p])
    sector_indices <- circlize::get.all.sector.index()
    sub_group <- sub_group[which(sub_group %in% sector_indices)]
    
    circlize::highlight.sector(
      sector.index = sub_group,
      track.index = 2,
      col = "black",
      text = p,
      font = 2,
      cex = 1,
      text.vjust = -0.5,
      facing = "bending.inside",
      niceFacing = TRUE
    )
  }
  return(plot)
}

#### generate the v-j counts matrix for each sample ####
v_j_matrix <- function(.data, threshold = 3){
  
  .data$vgene_allele <- str_extract(.data$V.name, ".+?(?=\\x28)")
  .data$cgene <- str_extract(.data$C.name, ".+?(?=\\*)")
  .data$vgene <- str_extract(.data$V.name, ".+?(?=\\*)")
  .data$dgene <- str_extract(.data$D.name, ".+?(?=\\*)")
  .data$jgene <- str_extract(.data$J.name, ".+?(?=\\*)")
  
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

# calculate kmer probability based on average CDR3 length in each group:
## use only the seqs where CDR3.aa length == avg CDR3.aa length
## (this method is similar to Platypus' approach w/ candidate sequences)
logoplots_fun <- function(x) {
  avg_cdr3_length <- mean(str_length(x$CDR3.aa)) %>% round
  x_avg <- x[str_length(x$CDR3.aa) == avg_cdr3_length,]
  plt <-
    getKmers(x_avg, .k = avg_cdr3_length) %>% kmer_profile(.method = "prob") %>%  vis(.plot = "seq")
  return(plt)
}

# Clones per Kiloread (CPK) - Expansion Index:
clones_per_kilo <- function(.data){
  num_CDR3 <- nrow(.data)
  sum_of_counts = sum(.data$Clones)
  normalized_data = (num_CDR3 / sum_of_counts)
}

# function to calculate Levenshtein distance matrix within each cluster
calculate_distance_matrix <- function(cluster_data) {
  # Convert CDR3.aa sequences to matrix
  aa_seqs <- as.matrix(as.AAbin(AAStringSet(cluster_data$CDR3.aa)))
  # Calculate pairwise Levenshtein distance matrix
  dist_matrix <- stringdistmatrix(cluster_data$CDR3.aa, cluster_data$CDR3.aa, method = "lv",useNames = "string")
  return(dist_matrix)
}

