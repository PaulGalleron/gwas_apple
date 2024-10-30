# fonction gwas


process_data <- function(pheno, geno, panel, filter_groups, pheno_percent, maf_threshold, pheno_column) {
  library(data.table)
  library(dplyr)

  ##############################################################################
  ########################### Preparation of myY ###############################
  ##############################################################################
  
  # Prepare panel_resume and pheno_groupe
  panel_resume <- panel[, c("MUNQ", "groupe_rapport_simple")]
  pheno_groupe <- merge(pheno, panel_resume, by = "MUNQ")

  # Filter groups
  filtered_pheno_groupe <- pheno_groupe %>%
    filter(groupe_rapport_simple %in% filter_groups)

  # Filter genotypes present in geno data
  MUNQ_320k <- data.frame(MUNQ = colnames(geno)[-c(1:5)])
  filtered_pheno_groupe_MUNQ_320k <- merge(filtered_pheno_groupe, MUNQ_320k, by = "MUNQ")

  # Function to select top percent of data based on a column
  select_top_percent <- function(df, column_name, percent) {
    if (percent <= 0 || percent > 100) stop("Invalid percentage")
    df_sorted <- df[order(df[[column_name]], decreasing = TRUE), ]
    num_to_select <- ceiling(nrow(df) * percent / 100)
    return(df_sorted[1:num_to_select, ])
  }

  # Filter based on phenotypic percentage
  filtered_pheno_groupe_MUNQ_320k_seuil <- select_top_percent(filtered_pheno_groupe_MUNQ_320k, pheno_column, pheno_percent)

  # Prepare myY by removing 'groupe_rapport_simple'
  myY <- filtered_pheno_groupe_MUNQ_320k_seuil[, !(names(filtered_pheno_groupe_MUNQ_320k_seuil) %in% "groupe_rapport_simple")]

  
  ##############################################################################
  ########################### Preparation of myG ###############################
  ##############################################################################
  
  geno_data <- as.data.table(geno)

  # Rename the first two columns
  setnames(geno_data, old = c(names(geno_data)[1], names(geno_data)[2]), new = c("markers", "chrom"))

  # Add the alleles column
  geno_data$alleles <- paste0(geno_data$A1, "/", geno_data$A2)

  # Reorder columns with alleles in the second position
  new_order <- c("markers", "alleles", "chrom", setdiff(names(geno_data), c("markers", "alleles", "chrom", "A1", "A2")))
  setcolorder(geno_data, new_order)

  # Remove columns A1 and A2
  geno_data <- geno_data[, !c("A1", "A2"), with = FALSE]

  # Create new columns filled with NA
  new_columns <- data.frame(
    strand = NA,
    assembly = NA,
    center = NA,
    protLSID = NA,
    assayLSID = NA,
    panel = NA,
    Qccode = NA
  )

  # Combine the new columns with geno_data
  geno_data <- cbind(geno_data[, 1:4], new_columns, geno_data[, 5:ncol(geno_data)])

  # Extract the first 11 columns
  first_cols <- geno_data[, 1:11]

  # Remove the first 11 columns from geno_data
  geno_data <- geno_data[, -(1:11)]

  # Select the MUNQ columns that match myY
  MUNQ_selection <- myY$MUNQ
  MUNQ_in_geno <- intersect(MUNQ_selection, colnames(geno_data))
  geno_selection <- geno_data[, ..MUNQ_in_geno, with = FALSE]


  ################### NA imputation

  # Fonction
  mode_function <- function(x) {
    unique_x <- unique(na.omit(x))
    return(unique_x[which.max(tabulate(match(x, unique_x)))])
  }

  geno_selection_imputed <- t(apply(geno_selection, 1, function(row) {
    row[is.na(row)] <- mode_function(row)
    return(row)
  }))

  geno_selection_imputed <- as.data.frame(geno_selection_imputed)


  # Combine first_cols with the selected geno columns
  data_geno_gapit <- cbind(first_cols, geno_selection_imputed)

  # Ensure columns 5 to 11 have NA values
  cols_to_na <- 5:11
  for (i in cols_to_na) {
    data_geno_gapit[, (i) := NA]
  }

  # Remove homomorphic markers
  row_variances <- apply(data_geno_gapit[, 12:ncol(data_geno_gapit)], 1, function(x) var(x, na.rm = TRUE))
  rows_to_keep <- which(row_variances != 0)
  data_geno_gapit_cleaned <- data_geno_gapit[rows_to_keep, ]

  # Calculate MAF and filter rows
  calculate_maf <- function(genotype_row) {
    allele_counts <- table(factor(genotype_row, levels = c(0, 1, 2)))
    total_alleles <- sum(allele_counts) * 2
    p_A <- (allele_counts["0"] * 2 + allele_counts["1"]) / total_alleles
    p_B <- (allele_counts["2"] * 2 + allele_counts["1"]) / total_alleles
    maf <- min(p_A, p_B)
    return(maf)
  }

  maf_values <- apply(data_geno_gapit_cleaned[, 12:ncol(data_geno_gapit_cleaned)], 1, calculate_maf)
  rows_to_keep <- which(maf_values >= maf_threshold)
  data_geno_gapit_cleaned_2 <- data_geno_gapit_cleaned[rows_to_keep, ]

  # Prepare myG by keeping column headers in the first row
  headers <- colnames(data_geno_gapit_cleaned_2)
  headers_df <- as.data.table(t(headers))

  # Renaming columns with V1, V2, ... for the rest of the data
  num_cols <- ncol(data_geno_gapit_cleaned_2)
  new_col_names <- paste0("V", 1:num_cols)
  setnames(data_geno_gapit_cleaned_2, new_col_names)

  # Bind the headers as the first row of the data
  myG_with_headers <- rbindlist(list(headers_df, as.data.table(data_geno_gapit_cleaned_2)), fill = TRUE)

  return(list(myY = myY, myG = myG_with_headers))
}



# Exemple d'appel de la fonction
# result <- process_data(
#  pheno = pheno,
#  geno = geno,
#  panel = panel,
#  filter_groups = c("cidre", "dessert_moderne", "dessert_ancien"),
#  pheno_percent = 100,
#  maf_threshold = 0.05,
#  pheno_column = "BLUP_AUDPC_reel"
# )
