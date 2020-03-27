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

# Make colour lookup table for clusters
cluster_colourmap <- function(x){
  x <- as.character(unique(x))
  structure(AA_COLOURS[str_sub(x, end = 1)], names=x)
}

# Add ggtext markdown features to a string or factor (features added as needed) based on lookup maps
# returns a factor sorted by the original sort order (alphabetical or current levels)
# Currently only works with simple 'span' class
add_markdown <- function(x, colour=NULL, order=NULL){
  if (is.null(order)){
    if (is.factor(x)){
      order <- levels(x)
      x <- as.character(x)
    } else {
      order <- sort(unique(x))
    }
  } else {
    x <- as.character(x)
  }
  
  order_codes <- rep('', length(order))
  x_codes <- rep('', length(x))
  
  if (!is.null(colour)){
    order_codes <- str_c(order_codes, 'color:', colour[order], ';')
    x_codes <- str_c(x_codes, 'color:', colour[x], ';')
  }
  
  # Add other features when needed
  
  order <- str_c("<span style = '", str_remove(order_codes, ';$'), "'>", order, "</span>")
  x <- str_c("<span style = '", str_remove(x_codes, ';$'), "'>", x, "</span>")
  
  return(factor(x, levels = order))
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
  
  # Calculate mean distance all points and a cluster (x)
  get_silhouette_distance <- function(x){
    ind <- tbl$cluster == x
    if (sum(ind) == 1){
      return(distance[,ind])
    } else {
      return(rowMeans(distance[,ind], na.rm = TRUE))
    }
  }
  mean_dists <- sapply(as.character(unique(tbl$cluster)), get_silhouette_distance)
  
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

# Plot gene data heatmaps
plot_gene_er <- function(dms, target_gene){
  dms_gene <- filter(dms, gene == target_gene) %>% 
    select(position, wt, A:Y) %>% 
    pivot_longer(A:Y, names_to = 'mut', values_to = 'er') %>% 
    mutate(er = clamp(er, 2, -2),
           mut = add_markdown(mut, colour = AA_COLOURS))
  breaks <- pretty_break(dms_gene$er, rough_n = 5, sig_figs = 3, sym = 0)
  
  ggplot(dms_gene, aes(x=mut, y=position, fill=er)) + 
    geom_tile() +
    geom_tile(data = filter(dms_gene, wt==mut), fill='grey') +
    scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction,
                         limits = breaks$limits, breaks=breaks$breaks, labels=breaks$labels) + 
    guides(fill=guide_colourbar(title = 'ER')) + 
    labs(title = target_gene, caption = 'Nb. |ER| clamped to < 2') + 
    theme(axis.ticks = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_markdown(),
          plot.caption = element_text(hjust = 0.5))
}

plot_gene_pcs <- function(dms, target_gene){
  dms_gene <- filter(dms, gene == target_gene) %>% 
    select(position, wt, PC1:PC20) %>% 
    pivot_longer(PC1:PC20, names_to = 'PC', values_to = 'val', names_prefix = 'PC') %>%
    mutate(PC = as.integer(PC))
  breaks <- pretty_break(dms_gene$val, rough_n = 5, sig_figs = 3, sym = 0)
  
  ggplot(dms_gene, aes(x=PC, y=position, fill=val)) + 
    geom_tile() +
    scale_fill_distiller(type = 'div', palette = 'PuOr', direction = -1,
                         limits = breaks$limits, breaks=breaks$breaks, labels=breaks$labels) + 
    guides(fill=guide_colourbar(title = 'PC Value')) + 
    labs(title = target_gene, x='PC', y='') + 
    theme(axis.ticks = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.background = element_blank())
}