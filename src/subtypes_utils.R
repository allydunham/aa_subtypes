#!/usr/bin/env Rscript
# Utility functions for Subtypes project

# Consistent filename from gene
gene_to_filename <- function(x){
  return(str_replace_all(str_to_lower(x), ' ', ''))
}

# Function to generate a pretty study name from the id using YAML config data
format_study <- function(x, max_width=60, study_dir='data/studies', mark_filtered=FALSE){
  yaml <- read_yaml(str_c(study_dir, '/', x, '/', x, '.yaml'))
  
  study <- str_c(yaml$authour, ' ', yaml$year, ' (', yaml$gene, ')') %>%
    str_wrap(width = max_width)
  
  if (mark_filtered & yaml$qc$filter){
    study <- str_c(study, '*')
  }
  
  return(study)
}

# Make pretty p value categories for plots
pretty_p_values <- function(p, breaks = c(0.001, 0.01, 0.05)){
  p_out <- rep(str_c('> ', breaks[length(breaks)]), length(p))
  for (cutoff in sort(breaks, decreasing = TRUE)){
    p_out[p < cutoff] <- str_c('< ', cutoff)
  }
  p_out <- factor(p_out, levels = c(str_c('> ', breaks[length(breaks)]),
                                    sapply(sort(breaks, decreasing = TRUE), function(x){str_c('< ', x)})))
  return(p_out)
}

# Split numeric vector into limits/break points/labels inteligently for display on legends etc.
# Using rough_n is approximate
# Sym allows you to make the range symetrical around a given value (usually 0)
pretty_break <- function(x, step=NULL, rough_n=NULL, sig_figs=4, sym=NULL){
  if (is.null(step) & is.null(rough_n)){
    rough_n <- 3
  }
  
  limits <- range(x, na.rm = TRUE)
  if (!is.null(sym)){
    m <- max(abs(limits - sym))
    limits <- m * c(-1, 1) + sym
  }
  
  if (is.null(step) & !is.null(rough_n)){
    step = signif((limits[2] - limits[1])/rough_n, 1)
  }
  
  limits_trunc <- trunc(limits/step) * step
  breaks <- seq(limits_trunc[1], limits_trunc[2], step)
  if (abs(breaks[1]) < abs(limits[1])){
    breaks <- c(limits[1], breaks)
  }
  if (abs(breaks[length(breaks)]) < abs(limits[2])){
    breaks <- c(breaks, limits[2])
  }
  
  labels <- str_remove_all(signif(breaks, sig_figs), "\\.0*$")
  
  if (abs(breaks[1] - breaks[2]) < 0.5 * step){
    labels[2] <- ''
  }
  if (abs(breaks[length(breaks)] - breaks[length(breaks) - 1]) < 0.5 * step){
    labels[length(breaks) - 1] <- ''
  }
  
  return(list(limits=limits, breaks=breaks, labels=labels))
}

# Clamp a value between two limts
clamp <- function(x, upper=Inf, lower=-Inf){
  x[x > upper] <- upper
  x[x < lower] <- lower
  return(x)
}

# Offset Uniprot positions to PDB positions
offset_uniprot_position <- function(position, sections){
  for (sec in sections){
    if ('region' %in% names(sec)){
      reg <- sec$region
    } else {
      reg <- c(-Inf, Inf)
    }
    pos <- position - sec$offset
    if (pos > reg[1] & pos < reg[2]){
      return(list(sec$chain, pos))
    }
  }
  return(list(NA, NA))
}

# Get PDB positions for all positions in a deep scanning tbl
dms_pdb_positions <- function(tbl, sections){
  x <- map2(tbl$position, tbl$gene, ~offset_uniprot_position(.x, sections[[gene_to_filename(.y)]]))
  chn <- map_chr(x, extract2, 1)
  pos <- as.integer(map_dbl(x, extract2, 2))
  return(list(chain=chn, position=pos))
}

# Calculate silhouette data for a clustering
# expects a tibble with columns cluster and cols, with cols giving the position
cluster_silhouette <- function(tbl, cols, distance_method = 'manhattan'){
  cols <- enquo(cols)
  
  if (n_distinct(tbl$cluster) == 1){
    warning('Silhouette undefined for a single cluster, returning NA')
    return(rep(NA, nrow(tbl)))
  }
  
  # Calculate distance between all points
  mat <- tibble_to_matrix(tbl, !!cols)
  if (distance_method == 'cosine'){
    distance <- cosine_distance_matrix(mat)
  } else {
    distance <- as.matrix(dist(mat, method = distance_method))
  }
  diag(distance) <- NA
  
  # Calculate mean distance each point and all clusters
  get_silhouette_distance <- function(x){
    ind <- tbl$cluster == x
    if (sum(ind) == 1){
      return(distance[,ind])
    } else {
      return(rowMeans(distance[,tbl$cluster == x], na.rm = TRUE))
    }
  }
  mean_dists <- sapply(unique(tbl$cluster), get_silhouette_distance)
  
  # Calculate mean dist within cluster
  same_cluster <- cbind(1:nrow(mean_dists), match(tbl$cluster, colnames(mean_dists)))
  a <- mean_dists[same_cluster]
  
  # Calculate mean dist to other clusters
  mean_dists[same_cluster] <- NA
  b <- apply(mean_dists, 1, min, na.rm=TRUE)
  
  # Calculate silhouette score
  s <- (b - a)/(pmax(a, b))
  
  singlton_clusters <- table(tbl$cluster)
  singlton_clusters <- names(singlton_clusters)[singlton_clusters == 1]
  
  s[tbl$cluster %in% singlton_clusters] <- 0
  
  return(s)
}

row_cosine_similarity <- function(x, y){
  rowSums(x*y) / (sqrt(rowSums(x^2) * rowSums(y^2)))
}

cosine_similarity_matrix <- function(mat){
  combs <- expand.grid(1:nrow(mat), 1:nrow(mat))
  cosine <- row_cosine_similarity(mat[combs$Var1,], mat[combs$Var2,])
  m <- matrix(cosine, nrow = nrow(mat), ncol = nrow(mat))
  rownames(m) <- rownames(mat)
  colnames(m) <- rownames(mat)
  return(m)
}

cosine_distance_matrix <- function(mat){
  acos(cosine_similarity_matrix(mat)) / pi
}