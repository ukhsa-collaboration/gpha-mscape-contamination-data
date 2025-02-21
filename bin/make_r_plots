#!/usr/bin/env Rscript

library(ggplot2)
library(vegan)
library(dplyr)
library(crayon)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

dir.create(output_dir, showWarnings = FALSE)

total_files <- list.files(input_dir, pattern = "*.total.txt", full.names = TRUE)
dna_files <- list.files(input_dir, pattern = "*.dna.txt", full.names = TRUE)
rna_files <- list.files(input_dir, pattern = "*.rna.txt", full.names = TRUE)

#eventually change this so that the number of mscape_labels are carried over from input
mscape_labels = c("BIRM", "GSTT (N2)", "GSTT (NTC)", "", "", "", "", "", "", "", "")
colours = c("pink", "lightblue", "olivedrab3", "gray","gray","gray","gray","gray","gray","gray","gray")

richness <- read.table(file.path(paste0(input_dir, "richness_table.txt")), stringsAsFactors = TRUE, header=TRUE, row.names="index")

taxa_list <- c(1, 2, 3)

for (type_of_taxa in taxa_list) {
  
  #Run code below for All, DNA, and RNA
  if (type_of_taxa == 1) {
    name_list <- lapply(strsplit(total_files, split='/|\\.'), tail, 3) ## Output a list of dataset names ## Output a list of dataset names
    file_list <- total_files
  } else if (type_of_taxa == 2) {
    name_list <- lapply(strsplit(dna_files, split='/|\\.'), tail, 3) ## Output a list of dataset names ## Output a list of dataset names
    file_list <- dna_files
  } else if (type_of_taxa == 3) {
    name_list <- lapply(strsplit(rna_files, split='/|\\.'), tail, 3) ## Output a list of dataset names ## Output a list of dataset names
    file_list <- rna_files
  }
  
  
  loop <- 1
  str_name_list <- c()
  median_shannon <- c()
  for (file in file_list) {
    #read data
    data <- read.table(file, header=TRUE)
    # Calculate Shannon's diversity index
    shannon <- (diversity(data, index = "shannon"))
    median_shannon <- append(median_shannon, median(shannon)) #get the median value in shannon calculations
    
    sub_list <- name_list[loop]
    current_list <- sub_list[[1]]
    current_name <- current_list[1]
    str_name_list <- append(str_name_list, current_name)
    
    if (loop == 1) {
      merge_df <- data.frame(
        Datasets = c(rep(current_name, length(shannon))),
        Shannon_Diversity = shannon
      )
    } else if (loop > 1) {
      new_df <- data.frame(
        Datasets = c(rep(current_name, length(shannon))),
        Shannon_Diversity = shannon
      )
      # merge two dataframes vertically
      merge_df <- rbind(merge_df, new_df)
    }
    loop <- loop + 1
  }
  
  #create a dataframe to sort the shannon boxes by descending order
  median_shannon_df <- data.frame(Names = str_name_list, Median = median_shannon)
  sorted_shannon <- median_shannon_df[order(-median_shannon_df$Median), ]
  
  #colour and label the mscape boxes
  colours <- c()
  mscape_labels <- c()
  for (name in sorted_shannon$Names) {
    if (name == "BIRM") {
      colours <- append(colours, "pink")
      mscape_labels <- append(mscape_labels, "BIRM")
    } else if (name == "GSTT_N2") {
      colours <- append(colours, "lightblue")
      mscape_labels <- append(mscape_labels, "GSTT (N2)")
    } else if (name == "GSTT_NTC") {
      colours <- append(colours, "olivedrab3")
      mscape_labels <- append(mscape_labels, "GSTT (NTC)")
    } else {
      colours <- append(colours, "gray")
      mscape_labels <- append(mscape_labels, "")
    }
  }
  
  # Convert Datasets to a factor with the desired order
  merge_df$Datasets <- factor(merge_df$Datasets, levels = sorted_shannon$Names)
  
  # Create a boxplot of Shannon's diversity for each niche
  p <- ggplot(merge_df, aes(x = Datasets, y = Shannon_Diversity, fill = Datasets)) +
    geom_boxplot() +
    scale_fill_manual(name = "Dataset",
                      labels = sorted_shannon$Names,
                      values = colours) +
    labs(title = "Diversity", x = "", y = "Shannon's Diversity") +
    stat_summary(geom = 'text', label = mscape_labels, fun = max, vjust = -1) +
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    guides(fill = FALSE)
  
  if (type_of_taxa == 1) {
    ggsave(file.path(paste0(output_dir, "total_diversity.png")), plot = p, width = 9, height = 6, dpi = 300) # Save the plot as a PNG
  } else if (type_of_taxa == 2) {
    ggsave(file.path(paste0(output_dir, "dna_diversity.png")), plot = p, width = 9, height = 6, dpi = 300) # Save the plot as a PNG
  } else if (type_of_taxa == 3) {
    ggsave(file.path(paste0(output_dir, "rna_diversity.png")), plot = p, width = 9, height = 6, dpi = 300) # Save the plot as a PNG
  }
  
  if (type_of_taxa == 1) {
    current_richness <- select(richness, All)
  } else if (type_of_taxa == 2) {
    current_richness <- select(richness, DNA)
  } else if (type_of_taxa == 3) {
    current_richness <- select(richness, RNA)
  }
  
  # To ensure not to divide by zero
  current_richness <- current_richness + 1
  transposed <- t(current_richness)
  
  loop <- 1
  str_name_list <- c()
  median_shannon <- c()
  for (file in file_list) {
    #read data
    data <- read.table(file, header=TRUE)
    # Calculate Shannon's diversity index
    shannon <- (diversity(data, index = "shannon"))
    
    sub_list <- name_list[loop]
    current_list <- sub_list[[1]]
    current_name <- current_list[1]
    str_name_list <- append(str_name_list, current_name)
    
    if (loop == 1) {
      merge_df <- data.frame(
        Datasets = c(rep(current_name, length(shannon))),
        Shannon_Evenness = shannon/transposed[loop]
      )
    } else if (loop > 1) {
      new_df <- data.frame(
        Datasets = c(rep(current_name, length(shannon))),
        Shannon_Evenness = shannon/transposed[loop]
      )
      # merge two dataframes vertically
      merge_df <- rbind(merge_df, new_df)
    }
    loop <- loop + 1
  }

  
  #colour and label the mscape boxes
  colours <- c()
  mscape_labels <- c()
  for (name in sorted_shannon$Names) {
    if (name == "BIRM") {
      colours <- append(colours, "pink")
      mscape_labels <- append(mscape_labels, "BIRM")
    } else if (name == "GSTT_N2") {
      colours <- append(colours, "lightblue")
      mscape_labels <- append(mscape_labels, "GSTT_N2")
    } else if (name == "GSTT_NTC") {
      colours <- append(colours, "olivedrab3")
      mscape_labels <- append(mscape_labels, "GSTT_NTC")
    } else {
      colours <- append(colours, "gray")
      mscape_labels <- append(mscape_labels, "")
    }
  }
  
  # Convert Datasets to a factor with the desired order from Diversity measures
  merge_df$Datasets <- factor(merge_df$Datasets, levels = sorted_shannon$Names)
  
  # Create a boxplot of Shannon's evenness for each niche
  p <- ggplot(merge_df, aes(x = Datasets, y = Shannon_Evenness, fill = Datasets)) +
    geom_boxplot() +
    scale_fill_manual(name = "Dataset",
                      labels = sorted_shannon$Names,
                      values = colours) +
    labs(title = "Evenness", x = "", y = "Shannon's Evenness") +
    stat_summary(geom = 'text', label = mscape_labels, fun = max, vjust = -1) +
    theme_bw() +
    theme(axis.text.x = element_blank())
  
  if (type_of_taxa == 1) {
    ggsave(file.path(paste0(output_dir, "total_evenness.png")), plot = p, width = 12, height = 6, dpi = 300) # Save the plot as a PNG
  } else if (type_of_taxa == 2) {
    ggsave(file.path(paste0(output_dir, "dna_evenness.png")), plot = p, width = 12, height = 6, dpi = 300) # Save the plot as a PNG
  } else if (type_of_taxa == 3) {
    ggsave(file.path(paste0(output_dir, "rna_evenness.png")), plot = p, width = 12, height = 6, dpi = 300) # Save the plot as a PNG
  }
  
}
