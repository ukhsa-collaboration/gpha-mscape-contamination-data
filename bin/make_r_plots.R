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
  if (type_of_taxa == 1) {
    #palette from https://sashamaps.net/docs/resources/20-colors/ in "convenient" order
    palette <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075')
    colours <- c()
    mscape_labels <- c()

    #record order of names to attach them to order of colour
    name_order <- c()

    colour_loop <- 1
    for (name in sorted_shannon$Names) {
      name_order <- append(name_order, name)
      if (grepl("public", tolower(name))) {
        colours <- append(colours, "gray")
        mscape_labels <- append(mscape_labels, "")
      } else {
        colours <- append(colours, palette[colour_loop])
        mscape_labels <- append(mscape_labels, name)
        colour_loop <- colour_loop + 1
      }
    }
    # Make and save colour_df in workflow
    colour_df <- cbind(data.frame(name_order), data.frame(colours))
    write.table(colour_df, file.path(paste0(output_dir, "colour_palette.txt")), sep="\t", quote = FALSE, row.names = FALSE) # Write tab delimiter text file

  } else {
    colours <- c()
    mscape_labels <- c()
    for (name in sorted_shannon$Names) {
        # Find corresponding output/colour
        colour <- colour_df$colours[colour_df$name_order == name]
        colours <- append(colours, colour) # This ensures that the same colours are used for datasets across figures
        mscape_labels <- append(mscape_labels, name)
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
    labs(title = "Shannon's Diversity", x = "", y = "") +
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
      # Find corresponding output/colour
      colour <- colour_df$colours[colour_df$name_order == name]
      colours <- append(colours, colour) # This ensures that the same colours are used for datasets across figures
      if (grepl("public", tolower(name))) {
        mscape_labels <- append(mscape_labels, "")
      } else {
        mscape_labels <- append(mscape_labels, name)
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
    labs(title = "Shannon's Evenness", x = "", y = "") +
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
