#!/usr/bin/env Rscript
# Functions for dimensionality reduction analysis

### tSNE/UMAP Plots ###
tsne_umap_plots <- function(tbl, x, y, name){
  x <- enquo(x)
  y <- enquo(y)
  p <- list(study = plot_dim_red_study(tbl, !!x, !!y),
            aa = plot_dim_red_aa(tbl, !!x, !!y),
            hydrophobicity = plot_dim_red_hydrophobicity(tbl, !!x, !!y),
            surface_accessibility = plot_dim_red_surface_accessibility(tbl, !!x, !!y),
            mean_er = plot_dim_red_mean_er(tbl, !!x, !!y))
  names(p) <- str_c(name, '_', names(p))
  return(p)
}

plot_dim_red_study <- function(tbl, x, y){
  x <- enquo(x)
  y <- enquo(y)
  (ggplot(tbl, aes(x=!!x, y=!!y, colour=study)) +
      geom_point() +
      facet_wrap(~study, labeller = as_labeller(sapply(unique(dms_wide$study), format_study, max_width=18)), nrow = 4) +
      guides(colour=FALSE)) %>%
    labeled_plot(units='cm', height = 20, width = 32)
}

plot_dim_red_aa <- function(tbl, x, y){
  x <- enquo(x)
  y <- enquo(y)
  (mutate(tbl, aa_class = AA_REDUCED_HASH[wt]) %>%
      ggplot(aes(x=!!x, y=!!y, colour=wt)) +
      geom_point() +
      facet_wrap(~aa_class) +
      scale_colour_manual(values = AA_COLOURS) +
      guides(colour = guide_legend(title='AA')) +
      theme(panel.grid.major.x = element_line(colour = 'gray', linetype = 'dotted'))) %>%
    labeled_plot(units='cm', height = 25, width = 25)
}

plot_dim_red_hydrophobicity <- function(tbl, x, y){
  x <- enquo(x)
  y <- enquo(y)
  (ggplot(tbl, aes(x=!!x, y=!!y, colour=hydrophobicity)) +
      geom_point() +
      scale_colour_gradient2() +
      guides(colour = guide_colourbar(title = 'Hydrophobicity'))) %>%
    labeled_plot(units='cm', height = 10, width = 15)
}

plot_dim_red_surface_accessibility <- function(tbl, x, y){
  x <- enquo(x)
  y <- enquo(y)
  (drop_na(tbl, all_atom_abs) %>%
      ggplot(aes(x=!!x, y=!!y, colour=all_atom_abs)) +
      geom_point() +
      scale_colour_viridis_c() +
      guides(colour = guide_colourbar(title = 'Surface Accessibility\n(All Atom Abs)'))) %>%
    labeled_plot(units='cm', height = 10, width = 15)
}

plot_dim_red_mean_er <- function(tbl, x, y){
  x <- enquo(x)
  y <- enquo(y)
  (drop_na(tbl, all_atom_abs) %>%
      ggplot(aes(x=!!x, y=!!y, colour=mean_score)) +
      geom_point() +
      scale_colour_gradient2() +
      guides(colour = guide_colourbar(title = 'Mean Norm. ER'))) %>%
    labeled_plot(units='cm', height = 10, width = 15)
}

