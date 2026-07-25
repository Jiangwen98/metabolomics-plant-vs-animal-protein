# ---------------------------------------------------------------------------
# ANALYSIS FUNCTIONS
#
# Utility functions used by the analysis scripts in this repository: quality
# control filtering, transformation, scaling, outlier handling, and plotting.
#
# These functions originate in earlier analysis code from the Drum Laboratory,
# inherited from a former lab member together with the B&B project. They are
# not my own work. Tags in the headers below preserve the distinction made in
# the original source: "adapted from earlier lab code" marks the functions I
# modified for the requirements of this analysis; the others are as inherited.
#
# Included so the analysis scripts can be read and run. See ATTRIBUTION.md.
# ---------------------------------------------------------------------------

# You MUST correctly indicate whether the data should be logged. The volcano output follows -lg(P-value) against log2(FC).
#' Plots a Volcano Plot of -lg(P-Value) against log2(Fold Change). -- adapted from earlier lab code
#'
#' @param df Dataframe. Dataframe to be used for plotting.
#' @param pval.threshold Integer. The threshold for the -lg(P-Value) value. Default is 1.3.
#' @param fc.threshold.left Integer. The negative threshold for the log2(Fold Change) value. Default is -1.
#' @param fc.threshold.right Integer. The positive threshold for the log2(Fold Change) value. Default is 1.
#' @param volc.title String. The title of the Volcano Plot. Default is blank ('').
#' @param text_size Integer. The text size of the labelled metabolites fulfilling the mentioned thresholds. Default is 2.5.
#' @param toLog2FC Boolean. Do you want to log2 the FC column?
#' @param toNegLog10Pval Boolean. Do you want to -log10 the Pval column?
#' @param x.title String. The x-axis title of the Volcano Plot. Default is 'log2(Fold Change)'.
#' @param y.title String. The y-axis title of the Volcano Plot. Default is '-lg(P-Value)'.
#' @param x_limit Numeric Vector. The x-limit window to set for the volcano plot. toLog2FC must be F. Default is NULL.
#' @param repel_labels Boolean. Do you want to repel the labels for better organization? geom_text_repel will be used. Default is F.
#' @param margins Numeric Vector of 4. To create margins around the graphs to include labels. c(y.max, x.max, y.min, x.min). Default is a vector of NAs.
#'
#' @return A ggplot2 object of a Volcano Plot of -lg(P-Value) against log2(Fold Change).
#' @export
#'
#' @examples plot_volcano1(df, volc.title = 'Group A vs Group B')
plot_volcano1 <- function (df, pval.threshold = 1.3, fc.threshold.left = NULL, fc.threshold.right = NULL,
                           volc.title = '', x.title = 'log2(Fold Change)', y.title = '-lg(p-value)',
                           text_size = 2.5, toLog2FC = F, toNegLog10Pval= F, x_limit = NULL, y_limit = NULL, repel_labels = F, dot_labels = T,
                           margins = NULL){
  require(tidyverse)
  require(ggrepel)
  # think about the 4 situation, toLog2FC = T/F & toNegLog10Pval = T/F, then you can understand the following -- Jiangwen
  # Checks for foldchange & p-value configuration. Consider all cases. X must be FC. Y must be Pval.
  # Untransformed Pval lies from 0-1. Untransformed FC usually hovers around 1, can be more or less than 1 but always positive.
  # Transformed Pval lies from 0 to infinity. Transformed FC can be negative now and hovers around 0.
  # if(toNegLog10Pval){ # Entirely Raw / Untransformed Pval => can filter for 0 < Pval < 1; let's hope your FC data ain't < 1 only.
  #   if(all((unlist(df[,2]) > 0) & (unlist(df[,2] <= 1)))){df <- df[,c(1,3,2)]}
  # }else if(!toLog2FC){ # All already Transformed => log2(FC) has both +ve and -ve, while -log(p-value) must always be positive
  #   if(any(unlist(df[,3] < 0))){df <- df[,c(1,3,2)]}
  # }else{ # Transformed Pval, non-transformed FC - we'll use distribution mean/median for comparisons. FC should have mean closer to 1. #not so reliable I think...- Jiangwen
  #   absmean1 <- abs(mean(unlist(df[,2]), na.rm = T) - 1)
  #   absmean2 <- abs(mean(unlist(df[,3]), na.rm = T) - 1)
  #   if (absmean2 < absmean1){df <- df[,c(1,3,2)]} #This means the df[,3] contains the FC, which should not be the case.
  # }
  
  ### -- Perform Transformations--- (The columns should be in the right positions)
  if(toLog2FC){df[,2] <- log(unlist(df[,2]), base = 2)}
  if(toNegLog10Pval){df[,3] <- log(unlist(df[,3]), base = 10) * -1}
  
  names(df) <- c('Name', 'X', 'Y') ### X is FC, Y is -log(p-value)... well it should be!
  
  ### -- Initializing color scheme --
  df$color<- 'gray'
  df$color[which(df$Y > pval.threshold)] <- 'black'
  
  if((!is.null(fc.threshold.left)) & (!is.null(fc.threshold.right)))
  {
    df$color[which(df$X > fc.threshold.right & df$Y > pval.threshold)] <- 'darkred'
    df$color[which(df$X < fc.threshold.left & df$Y > pval.threshold)] <- 'deepskyblue4'}
  
  df$color[which(df$Name == "Serum.Glucose.(mmol.per.l)")] <- "darkviolet"
  # df$color[which(df$Name == "isoleucine")] <- "darkgreen"
  # df$color[which(df$Name == "isoleucine_copy")] <- "darkgreen"
  # df$color[which(df$Name == "Leucine/Isoleucine")] <- "darkgreen"
  # df$color[which(df$Name == "Valine")] <- "darkgreen"
  # df$color[which(df$Name == "leucine")] <- "darkgreen"
  
  df <- df[order(df$color),]
  
  ## Only provide labels to those which are colored to avoid over-cluttering, while also need to show glucose in the plot...
  df$labels <- df$Name
  df$labels[which(df$color == 'gray')] <- ''
  df$labels[which(df$color == 'black')] <- ''
  # df$labels[which(df$color == 'darkred')] <- ''
  # df$labels[which(df$color == 'deepskyblue4')] <- ''
  df$labels[which(df$Name == "Serum.Glucose.(mmol.per.l)")] <- "Serum Glucose"
  # df$labels[which(df$Name == "isoleucine")] <- "Isoleucine feature1"
  # df$labels[which(df$Name == "isoleucine_copy")] <- "Isoleucine feature2"
  # df$labels[which(df$Name == "Leucine/Isoleucine")] <- "Leucine or Isoleucine"
  # df$labels[which(df$Name == "Valine")] <- "Valine"
  # df$labels[which(df$Name == "leucine")] <- "Leucine"
  
  df$left_labels <- df$Name
  df$left_labels[which(df$color == 'gray')] <- ''
  df$left_labels[which(df$color == 'black')] <- ''
  df$left_labels[which(df$Name == "Serum.Glucose.(mmol.per.l)")] <- "Serum Glucose"
  # df$left_labels[which(df$Name == "isoleucine")] <- "Isoleucine feature1"
  # df$left_labels[which(df$Name == "isoleucine_copy")] <- "Isoleucine feature2"
  # df$left_labels[which(df$Name == "Leucine/Isoleucine")] <- "Leucine or Isoleucine"
  # df$left_labels[which(df$Name == "Valine")] <- "Valine"
  # df$left_labels[which(df$Name == "leucine")] <- "Leucine"
  df$left_labels[which(df$X > 0)] <- ''
  
  df$right_labels <- df$Name
  df$right_labels[which(df$color == 'gray')] <- ''
  df$right_labels[which(df$color == 'black')] <- ''
  df$right_labels[which(df$Name == "Serum.Glucose.(mmol.per.l)")] <- "Serum Glucose"
  # df$right_labels[which(df$Name == "isoleucine")] <- "Isoleucine feature1"
  # df$right_labels[which(df$Name == "isoleucine_copy")] <- "Isoleucine feature2"
  # df$right_labels[which(df$Name == "Leucine/Isoleucine")] <- "Leucine or Isoleucine"
  # df$right_labels[which(df$Name == "Valine")] <- "Valine"
  # df$right_labels[which(df$Name == "leucine")] <- "Leucine"
  df$right_labels[which(df$X <= 0)] <- ''
  
  #Setting xlim and ylim if they exceed certain thresholds (which skews the graphs badly) - only if non-Logged FC
  # if(!is.null(x_limit)){
  #   if(!any(is.na(x_limit))){x_limit <- sort(x_limit)}
  #   if(!toLog2FC){
  #     if(!is.na(x_limit[1])){
  #       if(min(df$X) > x_limit[1]){x_limit[1] <- NA} # To not set an absurdly large window
  #     }
  #     if(!is.na(x_limit[2])){
  #       if(max(df$X) < x_limit[2]){x_limit[2] <- NA}
  #     }
  #   }
  #   if(all(is.na(x_limit))){rm(x_limit)}
  # }
  
  p1 <- ggplot2::ggplot(df, aes(x = X, y = Y, color = color)) + geom_point(size = 1.5)
  
  if((repel_labels == T) & (dot_labels == T)){ ### Repelling or not repelling
    p1 <- (p1 + coord_cartesian(clip = 'off') ### I'm modifiying this...
           # + geom_text_repel(aes(fontface = 'bold'),
           #                   label = unlist(df$labels),
           #                   size = text_size, alpha = 0.7, vjust = 0.5,
           #                   force = 10, hjust = 0.5,
           #                   segment.size = 0.6, segment.alpha = 0.5, max.overlaps = Inf,
           #                   )
           + geom_text_repel(aes(fontface = "bold"),
                             label = unlist(df$left_labels), # the left side, cold ones
                             size = text_size,
                             alpha = 0.85,
                             force = 25, #default force = 20 or 99
                             nudge_x = -0.3, nudge_y = 0.3,
                             segment.color = "gray56",
                             segment.linetype = 5,
                             segment.curvature = -1e-20,
                             direction = "both", 
                             segment.size = 0.6, segment.alpha = 1, max.overlaps = Inf
                             , seed = 123
           )
           + geom_text_repel(aes(fontface = "bold"),
                             label = unlist(df$right_labels), # the right side, hot ones
                             alpha = 0.85,
                             size = text_size,
                             force = 25,
                             nudge_x = 0.3, nudge_y = 0.3,
                             segment.color = "gray56",
                             segment.linetype = 5,
                             segment.curvature = 1e-20,
                             direction = "both", 
                             segment.size = 0.6, segment.alpha = 1, max.overlaps = Inf
                             ,seed = 123
           )
    ) ## unleash the labels
  }else if((repel_labels == F) & (dot_labels == T)){
    p1 <- p1 + geom_text(aes(fontface = 'bold'),
                         label = unlist(df$labels),
                         size = text_size, alpha = 0.7, vjust = 1.3,)
  }else if((repel_labels == T) & (dot_labels == F))
  {warning("Please define correct value for repel_labels and dot_labels!")
  }
  
  p1 <- (p1 + geom_hline(yintercept = pval.threshold, linetype = 'dotted', color = 'gray12', linewidth = 1) ## pval limit
         + geom_vline(xintercept = fc.threshold.right, linetype = 'dotted', color = 'gray12', linewidth = 1) ## fc limit right
         + geom_vline(xintercept = fc.threshold.left, linetype = 'dotted', color = 'gray12', linewidth = 1) ## fc limit left
         + geom_vline(xintercept = 0, linetype = 'solid', color = 'gray11', linewidth = 1) ## 0 Point 
         + xlab(x.title) + ylab(y.title) + ggtitle(volc.title) ## Labelling
         + theme_minimal() + scale_color_manual(values = c(unique(df$color)))
         + theme(plot.title = element_text(face = 'bold', size = 16, hjust = 0.5),
                 axis.title.x = element_text(face = 'bold', size = 12),
                 axis.title.y = element_text(face = 'bold', size = 12),
                 legend.position = 'none')) ## Font size/type configuration
  
  if(!is.null(margins)){p1 <- p1 + theme(plot.margin = unit(margins, 'in'))}
  if(!is.null(x_limit) & exists('x_limit')){p1 <- p1 + xlim(x_limit)}
  if(!is.null(y_limit) & exists('y_limit')){p1 <- p1 + ylim(y_limit)}
  
  return(p1)
}
