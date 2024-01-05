# Functions for cell divisions in polyG project -------------------------
library(data.table)
library(ape)
library(ggtree)
library(phytools)
library(patchwork)
library(tidyverse)

theme_martin  <- function () {
    theme_classic() + theme(axis.line = element_line(linewidth = 1),
        axis.title = element_text(size = 32, color = "black"),
        axis.text = element_text(size = 32, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.ticks.length = unit(0.2, "cm"), plot.title = element_text(size=22))
}

get_markerlengths  <- function(dir) { 

  subject <- str_split(dir, "/") %>%
    purrr::map(tail, n = 1) %>%
    unlist() %>%
    str_split("_") %>%
    purrr::map(1) %>%
    unlist()

  message(subject)

  ## path to poly-G raw data directory (marker length files)
  marker_dir <- list.files(dir, pattern = "repre_repli", full.names = TRUE)

  ## load marker lengths and get the average length of each marker in each sample
  get_markers <- function(marker_dir) {

  marker <- suppressMessages(read_tsv(marker_dir))
  
  # getting the marker name
  int_file <- tail(str_split(marker_dir, "/")[[1]], 1) 
  marker_name <- head(str_split(int_file, "_")[[1]], 1)

  cols <- NCOL(marker) +1 

  marker  %>% 
    rowid_to_column("length")  %>% 
    pivot_longer(cols=c(2:cols), names_to="sample", values_to="Frequency")  %>% 
    group_by(sample)  %>% 
    mutate(Frequency=Frequency/sum(Frequency),
            sample=str_remove(sample, "_[1-3]$")) %>%
    summarize(length = sum(length * Frequency)) %>%
    mutate(marker=marker_name) 
    
  }
  
  marker_files <- list.files(marker_dir, full.names = TRUE)

  markers <- lapply(marker_files, get_markers)
  markers <- bind_rows(markers)

    # subtract marker length of normal from each marker   
    root_sample <- str_remove(str_subset(unique(markers$sample), "N"), "_[1-3]$")[1] 
    
    # in the science cohort, the repre_repli data doesn't contain the final sample names
    # the samples won't be name "N", so the subtraction of the normal length must be 
    # performed after renaming 
    # thus this step needs to be skipped

if (str_detect(dir[1], "science|isc|hnscc")) {markers <- markers %>% 
      ungroup %>% 
      mutate(sample=str_remove(sample, "_[1-3]$"))  %>% 
      as.data.table } else {
    
    markers <- markers %>% 
      mutate(sample=str_remove(sample, "_[1-3]$"))  %>% 
      group_by(marker) %>% 
      filter(any(sample == root_sample)) %>% 
      mutate(length = length - length[sample == root_sample]) %>% 
      ungroup %>% 
      as.data.table
      }  

    markers$subject <- subject

 return(markers)
}


# Creating Heatmap from Mean Lengths --------------------------------------

create_heatmap  <- function(subject_i, mean_lengths) {

  # rotating mean length table for heatmap plotting
  filtered_matrix <- mean_lengths %>%
    filter(subject == subject_i) %>%
    mutate(sample=str_remove_all(sample, subject_i)) %>% 
    ungroup() %>% 
    select(-subject)  %>% 
    pivot_wider(names_from = marker, values_from=length)  %>% 
    column_to_rownames("sample")  %>% 
    as.matrix()

  # find first normal sample
  non_root <- str_subset(rownames(filtered_matrix), "N")[-1]

  filtered_matrix <- filtered_matrix[!rownames(filtered_matrix) %in% non_root,]

  # setting minmax to center the heatmap color at 0
  minmax <- max(abs(filtered_matrix), na.rm = TRUE)
  
  hm <- pheatmap::pheatmap(filtered_matrix,
    clustering_distance_rows = "manhattan",
    clustering_distance_cols = "manhattan",
    display_numbers = TRUE,
    number_color = "black",
    fontsize = 16,
    color = colorRampPalette(c("seagreen", "white", "purple3"))(42),
    na_col = "yellow",
    breaks = seq(-minmax, minmax, length.out = 43),
    silent = TRUE
  )
  hm
}


# Function to calculate L1 and cor all sample combinations ----------------

# function to get L1 and cor for a specific combination of two samples
get_l1_r_for_combination <- function(i, combos, markerlengths) {

  sample_a <- combos$a[i]
  sample_b <- combos$b[i]

  markerlengths_a  <- markerlengths %>% 
    dplyr::filter(sample==sample_a) %>% 
    arrange(marker)  %>% 
    pull(length)

  markerlengths_b  <- markerlengths %>% 
    dplyr::filter(sample==sample_b) %>% 
    arrange(marker)  %>% 
    pull(length)
  
  n_markers_a  <- length(markerlengths_a)
  n_markers_b  <- length(markerlengths_b)
  
  if (n_markers_a != n_markers_b) (error)
  
  l1  <- sum(abs(markerlengths_a - markerlengths_b))/n_markers_a
  
  r <- suppressWarnings(cor(markerlengths_a, markerlengths_b))
  
  list(a=sample_a, b=sample_b, l1=l1, r=r, marker=n_markers_a)
}


# Plotting MRCA trees -----------------------------------------------------

plot_tree <- function(subject_i, cell_divs_tbl) { 

  subject_i <- paste0(subject_i, "(?=[:alpha:])")

  sample_tbl <- cell_divs_tbl %>%
    filter(str_detect(a, subject_i), str_detect(b, subject_i)) %>%
    select(a, b, divs)

  # only select the first normal sample
  root_sample <- str_subset(unique(c(sample_tbl$a, sample_tbl$b)), "N")[1]

  sample_tbl <- sample_tbl  %>% 
    filter((a==root_sample)|str_detect(a, "N", negate = TRUE),
           (b==root_sample)|str_detect(b, "N", negate = TRUE)) 

  zero_tbl <- tibble("a"=unique(sample_tbl$a), "b"=unique(sample_tbl$a), "divs"=0) %>% 
    bind_rows(tibble("a"=unique(sample_tbl$b), "b"=unique(sample_tbl$b), "divs"=0)) %>%  
    distinct()
  
  dist_mat_int <- sample_tbl %>% 
    bind_rows(zero_tbl) %>% 
    pivot_wider(names_from=a, values_from=divs)

  # reordering the matrix to be in the correct order for the tree
  dist_mat <- dist_mat_int[match(colnames(dist_mat_int), dist_mat_int$b) %>% na.omit, ] %>% 
    column_to_rownames("b")

  # building an nj tree
  tree <- nj(as.dist(dist_mat))

  # rooting the tree
  phytools::reroot(tree, node.number=which(tree$tip.label==root_sample))
  
}

extract_mrca <- function(subject_i, tree, mean_lengths, width, height, save){
  
  # set mrca timing to be measured as NA in case a subject doesn't have the given sample type
  ad_ad_divs <- NA

  subject_ggtree <- ggtree(tree, size = 1.5) +
    theme_tree2() +
    labs(
      title = paste0(subject_i, " Cell Division Tree"),
      x = "Cell divisions from zygote"
    ) +
    theme(axis.line.x = element_line(linewidth = 1.5), 
        axis.title.x = element_text(size = 35, color = "black"), 
        axis.text.x = element_text(size = 35, color = "black"), 
        axis.ticks.x = element_line(color = "black", linewidth = 1), 
        axis.ticks.length.x = unit(0.2, "cm"), plot.title = element_text(size=35))

   # maximum y dimension for plotting
  max.y <- max(subject_ggtree$data$y)+1
  
  # maximum x dimension for plotting
  max.x <- max(subject_ggtree$data$x)*1.1
  
  # annotating the sample type
  subject_ggtree$data <- subject_ggtree$data %>%
    mutate(
      type = str_extract(label, "(?<=C[0-9]{1,3})[:alpha:]"),
      label = str_remove(label, subject_i)
    )
  
  # finding the mrca of all tumor samples
  all_mrca <- getMRCA(tree, c(str_subset(tree$tip.label, "N|Ad", negate = TRUE)))
  
  # finding the cell divisions from zygote to mrca 
  all_mrca_divs <- subject_ggtree$data[subject_ggtree$data$node==all_mrca & subject_ggtree$data$isTip==F,]$x

  p <- subject_ggtree +
    geom_tiplab(aes(color = type), size = 11) +
    coord_cartesian(xlim = c(0, max.x)) +
    scale_color_manual(
      values = c(A = "purple2", L = "darkgreen", M = "orange2", N = "grey23", P = "firebrick2"),
      guide = "none"
    )
  
  # add adenoma information if samples are present
  # adenomas with multiple samples:
  if (subject_i=="C90") {ad="Ad2"} else {ad="Ad1"}

  if (sum(str_detect(tree$tip.label, ad))>1){  
    
    ad_ad <- getMRCA(tree, c(str_subset(tree$tip.label, ad), str_subset(tree$tip.label, ad)))
    ad_ad_divs <- subject_ggtree$data[subject_ggtree$data$node==ad_ad & subject_ggtree$data$isTip==F,]$x
    
  }
  
  
if (save) {
  # create heatmap
  heatmap <- create_heatmap(subject_i, mean_lengths)
  p_hm <- p + ggtitle(subject_i)
  pdf(paste0("../plots/Supplementary_figures/crc_trees_hm/", subject_i, "_tree_hm.pdf"), width = width * 2, height = height)
  grid.arrange(p_hm, heatmap[[4]], widths = 1:2, vp = grid::viewport(width = 1, height = 0.7))
  dev.off()

  ggsave(paste0("../plots/Supplementary_figures/crc_trees/", subject_i, "_mrca_tree.pdf"), p,
    height = height,
    width = width
  )
}
list(subject=subject_i, all_mrca=all_mrca_divs, ad_ad_mrca=ad_ad_divs)   
}

plot_tree_hm_extract_mrcas  <- function(subject_i, cell_divs_tbl, mean_lengths, save) { 
  
  message(subject_i)
  
  tree <-  plot_tree(subject_i, cell_divs_tbl)
  
  l  <- extract_mrca(subject_i, tree, mean_lengths, 14, 12, save)
  
  return(l)
  
}


# Divergence Timing -------------------------------------------------------

# functions to get the divergence time of metastases from the primary tumor
# and to get the divergence times within the primary tumor


get_lesion_primary_mrca_dist <- function(lesion, tree) {
    # finding all the mrca nodes for this met lesion and all primary regions 
    get_primary_mrcas <- function(P,tree,lesion){
      mrca <- getMRCA(tree, c(str_subset(tree$tip.label, lesion), P))
      list(P=P, mrca=mrca)
    }
    nodes <- lapply(str_subset(tree$tip.label, "P"), get_primary_mrcas, tree, lesion)
    
    # extracting the node number
    node_list <- nodes %>% purrr::map(2) 
    
    # getting the distance of this node from the root
    get_node_height <- function(node, tree) {nodeheight(tree, node)}
    height_list <- lapply(node_list, get_node_height, tree)
    
    # the largest distance is the divergence time of this met lesion and the primary tumor
    met_primary_mrca <- max(unlist(height_list))
    
    list(samples=lesion, mrca = met_primary_mrca, comparison="met_primary")
}


  
get_primary_primary_mrca_dist <- function(i, cbm, tree) {
  s1 <- cbm[1, i]
  s2 <- cbm[2, i]

  height <- fastHeight(tree, s1, s2)
  samples <- paste0(s1, "_", s2)
  list(samples = samples, mrca = height, comparison = "primary_primary")
}

timing_pt_divergence <- function(subject_i, cell_divs_tbl) {

message(subject_i)    

  tree <- plot_tree(subject_i, cell_divs_tbl)
  
  # find unique lesions 
  unique_lesions <- tree$tip.label %>% 
    str_remove("[:alpha:]$") %>% 
    str_subset("M|L") %>% 
    unique()
  
  met_primary_mrca <- lapply(unique_lesions, get_lesion_primary_mrca_dist, tree)
  bind_rows(met_primary_mrca)
  
  # all primary tumor samples
  p_samples <- tree$tip.label %>% 
    str_subset("P")
  
  # all unique combinations of primary tumor samples
  cbm <- combn(p_samples, m=2) %>% as.data.frame() 
  
  primary_primary_mrca <- lapply(1:ncol(cbm), get_primary_primary_mrca_dist, cbm, tree) 

  ### Adenoma divergence 
  get_adenoma_adenoma_mrca_dist <- function(i, cbm, tree){
    
    s1 <- cbm[1,i]   
    s2 <- cbm[2,i]    
    if (str_extract(s1, "(?<=C[0-9]{1,3})Ad[0-9]") != str_extract(s2, "(?<=C[0-9]{1,3})Ad[0-9]")) return()
    height <- fastHeight(tree, s1, s2)
    samples <- paste0(s1, "_", s2)
    list(samples=samples, mrca=height, comparison="adenoma_adenoma")
    
  }
  
  # all primary tumor samples
  ad_samples <- tree$tip.label %>% 
    str_subset("Ad")
  
  # only calculate if more than one Adenoma sample is present 
  if ((length(ad_samples)) > 1) {

    # all unique combinations of adenoma tumor samples
    cbm <- combn(ad_samples, m = 2) %>% as.data.frame()
    adenoma_adenoma_mrca <- lapply(1:ncol(cbm), get_adenoma_adenoma_mrca_dist, cbm, tree)
  } else {adenoma_adenoma_mrca <- NULL}
  
  return(bind_rows(primary_primary_mrca, met_primary_mrca,adenoma_adenoma_mrca) %>% 
           mutate(subject=subject_i))
}



# Time to MRCA of a lesion ------------------------------------------------


get_all_mrcas <- function(i, tree, cbm) {
  
  s1 <- cbm[1,i]   
  s2 <- cbm[2,i]    
  
  mrca <- getMRCA(tree, c(str_subset(tree$tip.label, s1), str_subset(tree$tip.label, s2)))
  mrca_dist <- nodeheight(tree, mrca)

  list(samples=paste0(s1, "_", s2), mrca_dist=mrca_dist)
}

time_mrca <- function(subject_i, cell_divs_tbl) {
  
  message(subject_i)    
  
  sample_tbl  <- cell_divs_tbl  %>% 
    filter(str_detect(a, subject_i), str_detect(b, subject_i))  %>% 
    select(a,b,divs) 
  
  root_sample <- str_subset(unique(c(sample_tbl$a, sample_tbl$b)), "N")[1]
  
  # only select the first normal sample
  sample_tbl <- sample_tbl  %>% 
    filter((a==root_sample)|str_detect(a, "N", negate = TRUE),
           (b==root_sample)|str_detect(b, "N", negate = TRUE)) 
  
  zero_tbl <- tibble("a"=unique(sample_tbl$a), "b"=unique(sample_tbl$a), "divs"=0) %>% 
    bind_rows(tibble("a"=unique(sample_tbl$b), "b"=unique(sample_tbl$b), "divs"=0)) %>% 
    distinct()
  
  dist_mat_int <- sample_tbl %>% 
    bind_rows(zero_tbl) %>% 
    pivot_wider(names_from=a, values_from=divs)
  
  # reordering the matrix to be in the correct order for the tree
  dist_mat <- dist_mat_int[match(colnames(dist_mat_int), dist_mat_int$b) %>% na.omit, ] %>% 
    column_to_rownames("b")
  
  # building an nj tree
  tree <- nj(as.dist(dist_mat))
  
  # rooting the tree
  tree <- phytools::reroot(tree, node.number=which(tree$tip.label==root_sample))
  
  # find unique lesions 
  unique_lesions <- tree$tip.label %>% 
    str_remove("[:alpha:]$") %>% 
    str_subset("M|L") %>% 
    unique()

  mets <- if (sum(str_detect(unique_lesions, "M"))>1) {str_subset(unique_lesions, "M")} else {NA}
  ln <- if (sum(str_detect(unique_lesions, "L"))>1) {str_subset(unique_lesions, "L")} else {NA}
  
  mets_mrcas <- if (sum(is.na(mets))<1) {
    
    cbm <- combn(mets, m=2) %>% as.data.frame() 
    
    lapply(1:ncol(cbm), get_all_mrcas, tree, cbm)
  } else {list(samples=NA, mrca_dist=NA)}
  
  ln_mrcas <- if (sum(is.na(ln))<1) {
    
    cbm <- combn(ln, m=2) %>% as.data.frame() 
    
    lapply(1:ncol(cbm), get_all_mrcas, tree, cbm)
  } else {list(samples=NA, mrca_dist=NA)}
  
  return(bind_rows(mets_mrcas, ln_mrcas) %>% 
           mutate(subject=subject_i))
}


