library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)
library(viridis)

################################################################################
############################## New GWAS Figures ################################
################################################################################

run_all_plots <- function(
    myG,
    myY,
    model,
    trait_pheno,
    y_axis_title
) {
  
  # Call the function to generate Manhattan plots
  generate_manhattan_plots(model = model, trait_pheno = trait_pheno)
  
  # Define the PVE file
  gapit_pve_file <- paste0("GAPIT.Association.PVE.", model, ".", trait_pheno, ".csv")
  
  # Call the function to generate box plots
  generate_box_plots(
    myG = myG,
    myY = myY,
    gapit_pve_file = gapit_pve_file,
    model = model,
    trait_pheno = trait_pheno,
    y_axis_title = y_axis_title
  )
  
  # Call the function to generate violin plots
  generate_violin_plots(
    myG = myG,
    myY = myY,
    gapit_pve_file = gapit_pve_file,
    model = model,
    trait_pheno = trait_pheno,
    y_axis_title = y_axis_title
  )
  
  # Display phenotypic data
  create_plots_pdf(myY = myY, trait = trait_pheno, axis_label = y_axis_title)
  
  # Chromosome map with density
  create_chromosome_map(myG)
  
  # Chromosome map with significant SNPs
  create_chromosome_map_snp(model, trait_pheno, myG)
}


################################################################################
############################### Call function ##################################
################################################################################

#run_all_plots(
#  myG = myG,
#  myY = myY,
#  model = "MLMM",
#  trait_pheno = "IP_Nplus",
#  y_axis_title = "IP_Nplus"  # Y-axis title for box plots and violin plots
#)


################################################################################
#################### Function to generate Manhattan plots ######################
################################################################################

generate_manhattan_plots <- function(model, trait_pheno) {
  fichier <- paste0("GAPIT.Association.GWAS_Results.", model, ".", trait_pheno, ".csv")
  gwas_results <- read.csv(fichier)
  
  gwas_results$logP <- -log10(gwas_results$P.value)
  
  gwas_results <- gwas_results %>%
    group_by(Chr) %>%
    mutate(chr_len = max(Pos)) %>%
    ungroup() %>%
    arrange(Chr, Pos) %>%
    mutate(cum_pos = Pos + cumsum(as.numeric(chr_len)) - as.numeric(chr_len))
  
  gwas_results$adj_pvalue <- p.adjust(gwas_results$P.value, method = "BH")
  
  bonferroni_threshold <- -log10(0.05 / nrow(gwas_results))
  fdr_threshold <- -log10(max(gwas_results$adj_pvalue[gwas_results$adj_pvalue <= 0.05], na.rm = TRUE))
  
  chr_labels <- gwas_results %>%
    group_by(Chr) %>%
    summarize(cum_pos = (max(cum_pos) + min(cum_pos)) / 2)
  
  # Calculate a position slightly below the lines for annotations
  bonferroni_label_position <- bonferroni_threshold - 0.5  # Adjust the value here if necessary
  fdr_label_position <- fdr_threshold - 0.5  # Adjust the value here if necessary
  
  manhattan_plot_bonf_fdr <- ggplot(gwas_results, aes(x = cum_pos, y = logP, color = factor(Chr))) +
    geom_point(alpha = 0.8, size = 1.3) +
    scale_color_manual(values = rep(c("seagreen", "seagreen2"), length(unique(gwas_results$Chr)))) +
    geom_hline(yintercept = bonferroni_threshold, linetype = "dashed", color = "red3") +
    geom_hline(yintercept = fdr_threshold, linetype = "dotted", color = "royalblue") +
    annotate("text", x = max(gwas_results$cum_pos) * 1.05, y = bonferroni_label_position, label = "Bonferroni", color = "red3", hjust = 0, vjust = 0.5, size = 4) +
    annotate("text", x = max(gwas_results$cum_pos) * 1.05, y = fdr_label_position, label = "FDR", color = "royalblue", hjust = 0, vjust = 0.5, size = 4) +
    scale_x_continuous(breaks = chr_labels$cum_pos, labels = as.character(chr_labels$Chr)) +
    labs(title = "Manhattan Plot with Bonferroni and FDR thresholds", x = "Chromosome", y = "-log10(P.value)") +
    theme_minimal() +
    theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5), plot.margin = margin(1, 2, 1, 1, "cm")) +
    coord_cartesian(clip = 'off')
  
  ggsave("remix_gapit_manhattan_plot_bonferronni_fdr.pdf", plot = manhattan_plot_bonf_fdr, width = 12, height = 6, units = "in", dpi = 300)
  
  message("Manhattan plot created")
  
  manhattan_plot_bonf <- ggplot(gwas_results, aes(x = cum_pos, y = logP, color = factor(Chr))) +
    geom_point(alpha = 0.8, size = 1.3) +
    scale_color_manual(values = rep(c("seagreen", "seagreen2"), length(unique(gwas_results$Chr)))) +
    geom_hline(yintercept = bonferroni_threshold, linetype = "dashed", color = "red3") +
    annotate("text", x = max(gwas_results$cum_pos) * 1.05, y = bonferroni_label_position, label = "Bonferroni", color = "red3", hjust = 0, vjust = 0.5, size = 4) +
    scale_x_continuous(breaks = chr_labels$cum_pos, labels = as.character(chr_labels$Chr)) +
    labs(title = "Manhattan Plot with Bonferroni threshold", x = "Chromosome", y = "-log10(P.value)") +
    theme_minimal() +
    theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5), plot.margin = margin(1, 2, 1, 1, "cm")) +
    coord_cartesian(clip = 'off')
  
  ggsave("remix_gapit_manhattan_plot_bonferronni.pdf", plot = manhattan_plot_bonf, width = 12, height = 6, units = "in", dpi = 300)
  
  message("Manhattan plot created")
  
  manhattan_plot_fdr <- ggplot(gwas_results, aes(x = cum_pos, y = logP, color = factor(Chr))) +
    geom_point(alpha = 0.8, size = 1.3) +
    scale_color_manual(values = rep(c("seagreen", "seagreen2"), length(unique(gwas_results$Chr)))) +
    geom_hline(yintercept = fdr_threshold, linetype = "dotted", color = "royalblue") +
    annotate("text", x = max(gwas_results$cum_pos) * 1.05, y = fdr_label_position, label = "FDR", color = "royalblue", hjust = 0, vjust = 0.5, size = 4) +
    scale_x_continuous(breaks = chr_labels$cum_pos, labels = as.character(chr_labels$Chr)) +
    labs(title = "Manhattan Plot with FDR threshold", x = "Chromosome", y = "-log10(P.value)") +
    theme_minimal() +
    theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5), plot.margin = margin(1, 2, 1, 1, "cm")) +
    coord_cartesian(clip = 'off')
  
  ggsave("remix_gapit_manhattan_plot_fdr.pdf", plot = manhattan_plot_fdr, width = 12, height = 6, units = "in", dpi = 300)
  
  message("Manhattan plot created")
}

# Example function call
# generate_manhattan_plots("MLMM", "IP_Nplus")  # correct

################################################################################
################# box plot and spreadsheet of significant snp ##################
################################################################################

generate_box_plots <- function(
    myG,
    myY,
    gapit_pve_file,
    model,
    trait_pheno,
    y_axis_title
) {
  # Load GAPIT data
  gwas_PVE <- read.csv(gapit_pve_file)
  
  # Rename the column "Phenotype_Variance_Explained(%)" to "variance"
  colnames(gwas_PVE)[ncol(gwas_PVE)] <- "variance"
  gwas_PVE$variance <- round(gwas_PVE$variance, 2)
  
  # Extract genotype header
  entete <- myG[1, ]
  
  # Filter SNPs based on those listed in gwas_PVE
  myG_filtered <- myG %>%
    filter(!!sym(names(myG)[1]) %in% gwas_PVE$SNP)
  
  # Reintegrate the header
  myG_filtered <- rbind(entete, myG_filtered)
  
  # Copy the coding of significant SNPs
  codage <- myG_filtered %>% dplyr::select(1:2)
  
  # Transform the first row into a header
  colnames(codage) <- as.character(codage[1, ])
  
  # Remove the first row now used as a header
  codage <- codage[-1, ]
  
  # Remove irrelevant columns (assuming columns 2 to 11 are irrelevant)
  myG_filtered <- myG_filtered[, -c(2:11)]
  
  # Transpose and reformat data
  myG_filtered_transposed <- t(myG_filtered)
  
  # Remove leading/trailing spaces from values
  myG_filtered_transposed <- apply(myG_filtered_transposed, 2, trimws)
  
  # Rename columns with names from the first row
  colnames(myG_filtered_transposed) <- as.character(unlist(myG_filtered_transposed[1, ]))
  myG_filtered_transposed <- myG_filtered_transposed[-1, ]
  
  # Convert to data frame
  myG_filtered_transposed <- as.data.frame(myG_filtered_transposed, stringsAsFactors = FALSE)
  
  # Convert columns to factors
  myG_filtered_transposed[] <- lapply(myG_filtered_transposed, as.factor)
  
  # Rename the first column of myG_filtered_transposed to match myY
  names(myG_filtered_transposed)[1] <- names(myY)[1]
  
  # Merge genotypic and phenotypic data
  merged_data <- merge(myY, myG_filtered_transposed, by = names(myY)[1], all.y = TRUE)
  
  # Create box plots
  num_columns <- ncol(merged_data)
  n_col_myY <- ncol(myY)
  k <- n_col_myY + 1
  
  for (i in k:num_columns) {
    column_name <- colnames(merged_data)[i]
    
    # Check if the column is present in merged_data
    if (!(column_name %in% colnames(merged_data))) {
      next
    }
    
    # Find the row in 'codage' that corresponds to the current marker
    marker_row <- codage[codage$markers == column_name, ]
    
    # If the marker is found and has allele coding
    if (nrow(marker_row) > 0 && !is.na(marker_row$alleles)) {
      alleles <- strsplit(marker_row$alleles, "/")[[1]]
      
      # Replace numeric codes with allele codes only if the value is 0, 1, or 2
      merged_data[[column_name]] <- ifelse(merged_data[[column_name]] == 0, paste0(alleles[1], alleles[1]),
                                           ifelse(merged_data[[column_name]] == 1, paste0(alleles[1], alleles[2]),
                                                  ifelse(merged_data[[column_name]] == 2, paste0(alleles[2], alleles[2]),
                                                         merged_data[[column_name]])))
    }
    
    # Create a box plot
    p <- ggplot(merged_data, aes(x = !!sym(column_name), y = !!sym(trait_pheno))) +
      geom_boxplot(aes(color = !!sym(column_name)), fill = "slategray2", alpha = 0.7, outlier.shape = NA) + 
      geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
      labs(title = paste("Box Plot for SNP", column_name, "\nModel:", model),
           x = column_name,
           y = y_axis_title) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5), # Horizontal rotation of labels
            legend.position = "none")
    
    # Find the information of the current marker
    marker_info <- gwas_PVE[gwas_PVE$SNP == column_name, ]
    chromosome <- ifelse(nrow(marker_info) > 0, marker_info$Chr, NA)
    position <- ifelse(nrow(marker_info) > 0, marker_info$Pos, NA)
    percentage <- ifelse(nrow(marker_info) > 0, paste0(round(marker_info$variance, 2), "%"), NA)
    
    # Add text annotations to the plot
    text_grob <- textGrob(
      label = paste("Chromosome:", chromosome, "\nPosition:", position, "\nPhenotypic variance percentage:", percentage),
      gp = gpar(fontsize = 10),
      hjust = 0,
      x = unit(0.02, "npc")
    )
    
    # Define the PDF file name for each box plot
    pdf_file_name <- paste0("remix_gapit_box_plot_", gsub(" ", "_", column_name), ".pdf")
    pdf(pdf_file_name, width = 8, height = 6)
    
    # Combine plot and text below and add to PDF
    grid.arrange(p, text_grob, ncol = 1, heights = c(4, 1))
    
    # Close the PDF file
    dev.off()
    
    # Display a message for each generated file (optional)
    message(paste("PDF file created:", pdf_file_name))
  }
  
  output_file_2 <- paste0("remix_gapit_list_MUNQ_SNP_", model, "_", trait_pheno, ".csv")
  # Save the CSV file
  write.csv(merged_data, file = output_file_2, row.names = FALSE)
  
  # Confirmation message
  message("The CSV file has been saved as: ", output_file_2)
}

# Call the function
#generate_box_plots(
#  myG = myG,
#  myY = myY,
#  gapit_pve_file = paste0("GAPIT.Association.PVE.", "MLMM", ".", "IP_Nplus", ".csv"),
#  model = "MLMM",
#  trait_pheno = "IP_Nplus",
#  y_axis_title = "IP_Nplus"  # Y-axis title
#)



################################################################################
####################### Violin plot of significant SNPs ########################
################################################################################

generate_violin_plots <- function(
    myG,
    myY,
    gapit_pve_file,
    model,
    trait_pheno,
    y_axis_title
) {
  # Load GAPIT data
  gwas_PVE <- read.csv(gapit_pve_file)
  
  # Rename the column "Phenotype_Variance_Explained(%)" to "variance"
  colnames(gwas_PVE)[ncol(gwas_PVE)] <- "variance"
  gwas_PVE$variance <- round(gwas_PVE$variance, 2)
  
  # Extract the header from the genotypes
  header <- myG[1, ]
  
  # Filter SNPs based on the SNPs listed in gwas_PVE
  myG_filtered <- myG %>%
    filter(!!sym(names(myG)[1]) %in% gwas_PVE$SNP)
  
  # Reattach the header
  myG_filtered <- rbind(header, myG_filtered)
  
  # Copy the coding of significant SNPs
  coding <- myG_filtered %>% dplyr::select(1:2)
  
  # Transform the first row into the header
  colnames(coding) <- as.character(coding[1, ])
  
  # Remove the first row now used as the header
  coding <- coding[-1, ]
  
  # Remove irrelevant columns (assuming columns 2 to 11 are irrelevant)
  myG_filtered <- myG_filtered[, -c(2:11)]
  
  # Transpose and reformat the data
  myG_filtered_transposed <- t(myG_filtered)
  
  # Remove spaces before/after values
  myG_filtered_transposed <- apply(myG_filtered_transposed, 2, trimws)
  
  # Rename columns with the names from the first row
  colnames(myG_filtered_transposed) <- as.character(unlist(myG_filtered_transposed[1, ]))
  myG_filtered_transposed <- myG_filtered_transposed[-1, ]
  
  # Convert to dataframe
  myG_filtered_transposed <- as.data.frame(myG_filtered_transposed, stringsAsFactors = FALSE)
  
  # Convert columns to factors
  myG_filtered_transposed[] <- lapply(myG_filtered_transposed, as.factor)
  
  # Rename the first column of myG_filtered_transposed to match myY
  names(myG_filtered_transposed)[1] <- names(myY)[1]
  
  # Merge genotype and phenotype data
  merged_data <- merge(myY, myG_filtered_transposed, by = names(myY)[1], all.y = TRUE)
  
  # Create the violin plots
  num_columns <- ncol(merged_data)
  n_col_myY <- ncol(myY)
  k <- n_col_myY + 1
  
  for (i in k:num_columns) {
    column_name <- colnames(merged_data)[i]
    
    # Check that the column is present in merged_data
    if (!(column_name %in% colnames(merged_data))) {
      next
    }
    
    # Find the row in 'coding' that corresponds to the current marker
    marker_row <- coding[coding$markers == column_name, ]
    
    # If the marker is found and it has allele coding
    if (nrow(marker_row) > 0 && !is.na(marker_row$alleles)) {
      alleles <- strsplit(marker_row$alleles, "/")[[1]]
      
      # Replace the numeric codes with allele codes only if the value is 0, 1, or 2
      merged_data[[column_name]] <- ifelse(merged_data[[column_name]] == 0, paste0(alleles[1], alleles[1]),
                                           ifelse(merged_data[[column_name]] == 1, paste0(alleles[1], alleles[2]),
                                                  ifelse(merged_data[[column_name]] == 2, paste0(alleles[2], alleles[2]),
                                                         merged_data[[column_name]])))
    }
    
    # If alleles is NA, no changes are made to the column
    
    # Create a violin plot
    p <- ggplot(merged_data, aes(x = !!sym(column_name), y = !!sym(trait_pheno))) +
      geom_violin(fill = "cyan4", color = "black", alpha = 0.5, trim = FALSE) +  # Plot the violins
      geom_boxplot(width = 0.1, alpha = 0.5) +
      labs(title = paste("Violin Plot for SNP", column_name, "\nModel:", model),
           x = column_name,
           y = y_axis_title) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5), # Rotate the labels horizontally
            legend.position = "none") # Hide the legend if not needed
    
    # Find information about the current marker
    marker_info <- gwas_PVE[gwas_PVE$SNP == column_name, ]
    chromosome <- ifelse(nrow(marker_info) > 0, marker_info$Chr, NA)
    position <- ifelse(nrow(marker_info) > 0, marker_info$Pos, NA)
    percentage <- ifelse(nrow(marker_info) > 0, paste0(round(marker_info$variance, 2), "%"), NA)
    
    # Add text annotations to the plot
    text_grob <- textGrob(
      label = paste("Chromosome:", chromosome, "\nPosition:", position, "\nPercentage of phenotypic variance:", percentage),
      gp = gpar(fontsize = 10),
      hjust = 0,
      x = unit(0.02, "npc")
    )
    
    # Define the PDF file name for each violin plot
    pdf_file_name <- paste0("remix_gapit_violin_plot_", gsub(" ", "_", column_name), ".pdf")
    pdf(pdf_file_name, width = 8, height = 6)
    
    # Combine the plot and the text below and add it to the PDF
    grid.arrange(p, text_grob, ncol = 1, heights = c(4, 1))
    
    # Close the PDF file
    dev.off()
    
    # Display a message for each generated file (optional)
    message(paste("PDF file created:", pdf_file_name))
  }
}

# Call the function
#generate_violin_plots(
#  myG = myG,
#  myY = myY,
#  gapit_pve_file = paste0("GAPIT.Association.PVE.", "MLMM", ".", "IP_Nplus", ".csv"),
#  model = "MLMM",
#  trait_pheno = "IP_Nplus",
#  y_axis_title = "IP_Nplus"  # Title for the Y-axis
#)



################################################################################
####################### Study of Phenotypic Data ###############################
################################################################################

# Function to create a PDF with four plots
create_plots_pdf <- function(myY, trait, axis_label) {
  
  # Box plot with points
  box_plot <- ggplot(myY, aes_string(x = "1", y = trait)) +
    geom_boxplot(fill = "lightblue", color = "black", outlier.shape = NA) +  # Box plot without outliers
    geom_jitter(width = 0.2, color = "blue", alpha = 0.5) +  # Add points with a slight horizontal offset
    ylab(axis_label) +
    ggtitle("Box Plot") +
    theme(axis.title.x = element_blank(),  # Remove x-axis label
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = rel(1.3)))
  
  # Simplified violin plot
  violin_plot <- ggplot(myY, aes_string(x = "1", y = trait)) +
    geom_violin(fill = "lightblue", color = "black") +  # Violin plot
    ylab(axis_label) +
    ggtitle("Violin Plot") +
    theme(axis.title.x = element_blank(),  # Remove x-axis label
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = rel(1.3)))
  
  # Histogram with a density line
  hist_plot <- ggplot(myY, aes_string(x = trait)) +
    geom_histogram(aes(y = after_stat(density)), fill = "blue", color = "black", bins = 20, alpha = 0.5) +
    geom_density(color = "red", size = 1) +
    xlab(axis_label) +
    ylab("Density") +
    ggtitle("Histogram with Density") +
    theme(axis.text.x = element_text(size = 10),
          axis.title = element_text(size = rel(1.3)))
  
  # Empirical cumulative distribution function (ECDF) plot
  ecdf_plot <- ggplot(myY, aes_string(x = trait)) +
    stat_ecdf(geom = "step", color = "blue", size = 1) +  # Plot cumulative density
    xlab(axis_label) +
    ylab("Cumulative Proportions") +
    ggtitle("ECDF Plot") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 10),
          axis.title = element_text(size = rel(1.3)))
  
  # Create the PDF with the four plots
  pdf("remix_gapit_pheno.pdf", width = 12, height = 8)
  grid.arrange(box_plot, violin_plot, hist_plot, ecdf_plot, ncol = 2)
  dev.off()
  
  message("The PDF file has been saved as: remix_gapit_pheno.pdf")
}

# Example usage
#create_plots_pdf(myY = myY, trait = "IP_Nplus", axis_label = "BLUP0")


################################################################################
#################### Chromosome Maps with Density ##############################
################################################################################

# Define the function
create_chromosome_map <- function(myG) {
  
  # Use the first row as header and remove it
  colnames(myG) <- as.character(myG[1, ])
  myG <- myG[-1, ]
  
  # Check if the necessary columns exist in myG
  if (!all(c("chrom", "pos") %in% colnames(myG))) {
    stop("The dataframe 'myG' does not contain the required columns 'chrom' and 'pos'.")
  }
  
  # Convert columns to appropriate types
  myG$chrom <- as.integer(as.factor(myG$chrom))
  myG$pos <- as.numeric(as.character(myG$pos)) / 1e6  # Convert positions to Mb
  
  # Prepare the data
  cat("Preparing the data...\n")
  data <- myG %>%
    dplyr::select(chrom, pos) %>%
    filter(!is.na(chrom) & !is.na(pos))
  
  # Calculate chromosome limits
  cat("Calculating chromosome limits...\n")
  chrom_lims <- data %>%
    group_by(chrom) %>%
    summarise(ymin = chrom - 0.5,
              ymax = chrom + 0.5,
              xmax = max(pos),  # Final position of each chromosome
              .groups = 'drop')
  
  # Create the chromosome map with marker density
  cat("Creating the chromosome map...\n")
  p <- ggplot(data, aes(x = pos, y = as.factor(chrom))) +
    geom_bin2d(bins = 100, aes(fill = after_stat(count))) +  # Add marker density
    scale_fill_viridis_c(option = "A", direction = -1, name = "Density") +  # Color scale
    labs(x = 'Position on the Chromosome (Mb)',
         y = 'Chromosome',
         title = 'Chromosome Map with Marker Density') +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 12))
  
  # Add segments for chromosomes
  p <- p +
    geom_segment(data = chrom_lims, 
                 aes(x = 0, xend = xmax, y = chrom, yend = chrom),
                 color = "black",
                 size = 1) +
    scale_x_continuous(labels = scales::comma)  # Add commas for positions in Mb
  
  # Save the map
  map_chromo <- "remix_gapit_chromosome_map.png"
  cat("Saving the map as", map_chromo, "...\n")
  tryCatch({
    ggsave(map_chromo, plot = p, width = 12, height = 8, dpi = 300, bg = "white")
  }, error = function(e) {
    stop("Error saving the file:", e$message)
  })
  
  # Return the name of the saved file
  return(map_chromo)
}

# Example usage
#map_chromo <- create_chromosome_map(myG)  # Ensure 'myG' is a dataframe available in your environment
# print(paste("Map saved as:", map_chromo))



################################################################################
################ Chromosome Map with Significant SNPs ##########################
################################################################################

# Function to create the chromosome map with significant SNPs
create_chromosome_map_snp <- function(model, trait_pheno, myG) {
  
  # Build the name of the GAPIT PVE file based on the model and the phenotypic trait
  gapit_file <- paste0("GAPIT.Association.PVE.", model, ".", trait_pheno, ".csv")
  
  # Load GAPIT data
  gapit_PVE <- read.csv(gapit_file)
  
  # Transform GAPIT data
  gapit_PVE <- gapit_PVE %>%
    rename(chrom = 2, pos = 3) %>%
    filter(!is.na(chrom) & !is.na(pos))
  
  # Prepare genome data from myG
  # Use the first row as the header
  colnames(myG) <- as.character(myG[1, ])
  myG <- myG[-1, ]
  
  # Convert columns to appropriate types
  myG <- myG %>%
    mutate(chrom = as.integer(chrom),  # Ensure chromosomes are integers
           pos = as.numeric(pos))  # Convert positions to numeric
  
  # Calculate minimum positions for each chromosome
  min_pos_per_chrom <- myG %>%
    group_by(chrom) %>%
    summarise(min_pos = min(pos), .groups = 'drop')
  
  # Adjust positions so that each chromosome starts at 0
  myG <- myG %>%
    left_join(min_pos_per_chrom, by = "chrom") %>%
    mutate(pos = (pos - min_pos) / 1e6) %>%  # Convert adjusted positions to Mb
    dplyr::select(-min_pos)  # Remove min_pos column after adjustment
  
  # Apply the same adjustment to the GAPIT data
  gapit_PVE <- gapit_PVE %>%
    left_join(min_pos_per_chrom, by = "chrom") %>%
    mutate(pos = (pos - min_pos) / 1e6) %>%  # Convert adjusted positions to Mb
    dplyr::select(-min_pos)  # Remove min_pos column after adjustment
  
  # Calculate maximum positions for each chromosome
  chrom_limits <- myG %>%
    group_by(chrom) %>%
    summarise(min_pos = min(pos), max_pos = max(pos), .groups = 'drop')
  
  # Create the chromosome map with rectangles and significant SNPs
  p <- ggplot() +
    # Draw rectangles for each chromosome (horizontally)
    geom_rect(data = chrom_limits, aes(ymin = min_pos, ymax = max_pos, xmin = chrom - 0.4, xmax = chrom + 0.4), 
              fill = "cyan4", color = "black", alpha = 0.5) +  # Chromosomes in blue with black border
    # Add significant SNPs with thicker red lines (horizontally)
    geom_segment(data = gapit_PVE, aes(y = pos, yend = pos, x = chrom - 0.5, xend = chrom + 0.5), 
                 color = "red", size = 1.5) +  # Increase segment size for QTLs
    labs(y = 'Position on Chromosome (Mb)',
         x = 'Chromosome',
         title = 'Chromosome Map with Significant SNPs') +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12),
          panel.grid.major = element_blank(),  # Remove major gridlines
          panel.grid.minor = element_blank(),  # Remove minor gridlines
          panel.spacing.x = unit(0.5, "lines")) +  # Add spacing between chromosomes
    scale_x_continuous(breaks = seq(min(myG$chrom), max(myG$chrom), by = 1), position = "top") +  # Chromosome numbers on top
    scale_y_reverse()  # Reverse the y-axis to have 0 at the top and increasing values downward
  
  # Display the plot on the screen
  print(p)
  
  # Define the output PDF filename based on the model and trait
  output_pdf <- paste0("remix_gapit_chromosome_map_", model, "_", trait_pheno, ".pdf")
  
  # Save the map to a PDF file
  ggsave(filename = output_pdf, plot = p, width = 8, height = 6, units = "in")
  
  # Message to indicate that the PDF has been saved
  message("Chromosome map saved to: ", output_pdf)
}

# Example usage of the function
# Assuming myG is already loaded in the environment
#create_chromosome_map_snp("MLMM", "IP_Nplus", myG)












################################################################################
############### Box Plot of Significant SNPs with Groups ######################
################################################################################

generate_box_plots_with_groups <- function(
    myG,
    myY,
    panel,
    model, 
    trait_pheno, 
    y_axis_title = "BLUP of AUDPC"  # Allow defining the y-axis title
) {
  # Construct the GAPIT PVE file name
  gapit_pve_file <- paste0("GAPIT.Association.PVE.", model, ".", trait_pheno, ".csv")
  
  # Load GAPIT data
  gwas_PVE <- read.csv(gapit_pve_file)
  
  # Rename the column "Phenotype_Variance_Explained(%)" to "variance"
  colnames(gwas_PVE)[ncol(gwas_PVE)] <- "variance"
  gwas_PVE$variance <- round(gwas_PVE$variance, 2)
  
  # Extract the header of the genotypes
  header <- myG[1, ]
  
  # Filter SNPs based on the SNPs listed in gwas_PVE
  myG_filtered <- myG %>%
    filter(!!sym(names(myG)[1]) %in% gwas_PVE$SNP)
  
  # Reintegrate the header
  myG_filtered <- rbind(header, myG_filtered)
  
  # Copy the coding of significant SNPs
  coding <- myG_filtered %>% dplyr::select(1:2)
  
  # Convert the first row into a header
  colnames(coding) <- as.character(coding[1, ])
  
  # Remove the first row now used as the header
  coding <- coding[-1, ]
  
  # Remove irrelevant columns (assuming columns 2 to 11 are irrelevant)
  myG_filtered <- myG_filtered[, -c(2:11)]
  
  # Transpose and reformat the data
  myG_filtered_transposed <- t(myG_filtered)
  
  # Remove leading/trailing whitespace from values
  myG_filtered_transposed <- apply(myG_filtered_transposed, 2, trimws)
  
  # Rename columns using the first row's names
  colnames(myG_filtered_transposed) <- as.character(unlist(myG_filtered_transposed[1, ]))
  myG_filtered_transposed <- myG_filtered_transposed[-1, ]
  
  # Convert to dataframe
  myG_filtered_transposed <- as.data.frame(myG_filtered_transposed, stringsAsFactors = FALSE)
  
  # Convert columns to factors
  myG_filtered_transposed[] <- lapply(myG_filtered_transposed, as.factor)
  
  # Rename the first column of myG_filtered_transposed to match myY
  names(myG_filtered_transposed)[1] <- names(myY)[1]
  
  # Merge genotype and phenotype data
  merged_data <- merge(myY, myG_filtered_transposed, by = names(myY)[1], all.y = TRUE)
  
  # Use the already loaded panel
  panel <- panel %>% dplyr::select(MUNQ, groupe_rapport_simple)
  
  # Merge panel data with merged_data using the MUNQ column
  merged_data <- merge(merged_data, panel, by = "MUNQ", all.x = TRUE)
  
  # Define colors and labels for the groups
  group_colors <- c("cidre" = "#FC4E07", "dessert_ancien" = "#009E73", "dessert_moderne" = "#E7B800", "sauvage" = "#CC79A7")
  group_labels <- c("cidre" = "cider", "dessert_ancien" = "old dessert", "dessert_moderne" = "modern dessert", "sauvage" = "wild")
  
   
  # Create box plots
  num_columns <- ncol(merged_data) - 1
  n_col_myY <- ncol(myY)
  k <- n_col_myY + 1
  
  for (i in k:num_columns) {
    column_name <- colnames(merged_data)[i]
    
    # Check if the column is present in merged_data
    if (!(column_name %in% colnames(merged_data))) {
      next
    }
    
    # Find the corresponding row in 'coding' for the current marker
    marker_row <- coding[coding$markers == column_name, ]
    
    # If the marker is found and it has allele coding
    if (nrow(marker_row) > 0 && !is.na(marker_row$alleles)) {
      alleles <- strsplit(marker_row$alleles, "/")[[1]]
      
      # Replace numeric codes with allele codes
      merged_data[[column_name]] <- ifelse(merged_data[[column_name]] == 0, paste0(alleles[1], alleles[1]),
                                           ifelse(merged_data[[column_name]] == 1, paste0(alleles[1], alleles[2]),
                                                  ifelse(merged_data[[column_name]] == 2, paste0(alleles[2], alleles[2]),
                                                         merged_data[[column_name]])))
    }
    
    # Create a box plot
    p <- ggplot(merged_data, aes(x = !!sym(column_name), y = !!sym(trait_pheno))) +
      geom_boxplot(fill = "slategray2", alpha = 0.7, outlier.shape = NA) +  # No coloring by group for boxplots
      geom_jitter(aes(fill = groupe_rapport_simple), width = 0.2, alpha = 0.7, 
                  shape = 21, size = 2, color = "black", stroke = 0.5) +  # Colored points with black outline
      scale_fill_manual(values = group_colors, labels = group_labels, name = "Groups") +  # Apply colors and labels for the legend
      labs(title = paste("Box Plot for SNP", column_name),
           x = column_name,
           y = y_axis_title) +  # Use the y-axis title defined in the function call
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal rotation of labels
        legend.position = "right",  # Show legend
        legend.title = element_blank()  # Remove the legend title
      )
    
    # Find information for the current marker
    marker_info <- gwas_PVE[gwas_PVE$SNP == column_name, ]
    chromosome <- ifelse(nrow(marker_info) > 0, marker_info$Chr, NA)
    position <- ifelse(nrow(marker_info) > 0, marker_info$Pos, NA)
    percentage <- ifelse(nrow(marker_info) > 0, paste0(round(marker_info$variance, 2), "%"), NA)
    
    # Add text annotations to the plot
    text_grob <- textGrob(
      label = paste("Chromosome:", chromosome, "\nPosition:", position, "\nPercentage of Phenotypic Variance:", percentage),
      gp = gpar(fontsize = 10),
      hjust = 0,
      x = unit(0.02, "npc")
    )
    
    # Define the PDF filename for each box plot
    pdf_file_name <- paste0("remix_gapit_box_plot_groups_", gsub(" ", "_", column_name), ".pdf")
    pdf(pdf_file_name, width = 8, height = 6)
    
    # Combine the plot and the text below and add to the PDF
    grid.arrange(p, text_grob, ncol = 1, heights = c(4, 1))
    
    # Close the PDF file
    dev.off()
    
    # Optional: Display a message for each generated file
    message(paste("PDF file created:", pdf_file_name))
  }
}

# Example function call
# generate_box_plots_with_groups(myG = myG, myY = myY, panel = panel, model = "MLMM", trait_pheno = "BLUP_AUDPC_reel", y_axis_title = "Your Custom Title Here")


################################################################################
##################### Structure with Groups on PCA #############################
################################################################################


# Define the function
plot_pca_scatter <- function(pca_file, panel, x_axis_title, y_axis_title) {
  # Load the PCA data from the CSV file
  gapit_pca <- read.csv(pca_file)
  
  # Rename the first column to "MUNQ"
  colnames(gapit_pca)[1] <- "MUNQ"
  
  # Ensure the panel data has the relevant columns
  if (!all(c("MUNQ", "groupe_rapport_simple") %in% colnames(panel))) {
    stop("The panel data must contain 'MUNQ' and 'groupe_rapport_simple' columns.")
  }
  
  # Select relevant columns from panel
  panel <- panel %>% dplyr::select(MUNQ, groupe_rapport_simple)
  
  # Merge the PCA data with the panel data based on "MUNQ"
  merged_data <- merge(gapit_pca, panel, by = "MUNQ")
  
  # Define the colors and labels for the groups
  group_colors <- c("cidre" = "#FC4E07", 
                    "dessert_ancien" = "#009E73", 
                    "dessert_moderne" = "#E7B800", 
                    "sauvage" = "#CC79A7")
  
  group_labels <- c("cidre" = "cider", 
                    "dessert_ancien" = "old dessert", 
                    "dessert_moderne" = "modern dessert", 
                    "sauvage" = "wild")
  
  # Create the scatter plot
  p <- ggplot(merged_data, aes(x = PC1, y = PC2, fill = groupe_rapport_simple)) +
    geom_point(size = 2, shape = 21, color = "black", stroke = 0.8) +  # Circles with black outline
    scale_fill_manual(values = group_colors, labels = group_labels, name = "Groups") +  # Custom colors and labels
    labs(title = "Representation of the first two components of the CPA",
         x = x_axis_title,
         y = y_axis_title) +
    theme_minimal() +  # Use a minimal theme for a clean look
    theme(legend.position = "right")  # Position the legend on the right
  
  # Define the output PDF file name
  pdf_filename <- "remix_gapit_PCA_Scatter_Plot.pdf"
  
  # Save the plot to a PDF
  ggsave(filename = pdf_filename, plot = p, device = "pdf", width = 8, height = 6)
  
  message("ACP with groups created")
}

# Example usage
# Assuming 'panel' is already loaded in the environment
# plot_pca_scatter("GAPIT.Genotype.PCA.csv", panel, "Principal Component 1 (PC1)", "Principal Component 2 (PC2)")

