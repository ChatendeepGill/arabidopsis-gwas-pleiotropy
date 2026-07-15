# Load packages
library(data.table)
library(dplyr)
library(purrr)
library(rhdf5)
library(stringr)

# Define inputs
# IMPORTANT: Replace these paths with the locations of the AraGWAS files on your computer. These files are downloaded separately and are not stored in this repository. 

# Directory containing the gene-specific SNP mapping CSV files downloaded from the AraGWAS Catalog. Each filename must begin with its TAIR10 gene identifier, followed by an underscore. Each file must contain: snp.chr, snp.position.
gene_dir <- "/path/to/downloaded/csv_files"

# Path to the compressed AraGWAS HDF5 database archive downloaded from the AraGWAS Catalog.
tar_file <- "/path/to/downloaded/aragwas_db.tar.gz"

# Define outputs
# Temporary directory used to extract one HDF5 file at a time
temp_dir <- "temp_extract"

# Main SNP-level output used in later analysis scripts
out_snp_file <- file.path("results", "aragwas_combined_results.csv")

# Processing log used to track progress and record errors
log_file <- file.path("results", "01_process_log.txt")

# Create output and temporary directories
dir.create("results", showWarnings = FALSE)
dir.create(temp_dir, showWarnings = FALSE)

# Define logging function
# Add a timestamped message to both the R console and the log file.
log_message <- function(msg, log_file) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  line <- paste(timestamp, msg)
  cat(line, "\n") # Print the message to the console
  write(line, file = log_file, append = TRUE) # Append the message to the log file
}

# Initialize log
write("", file = log_file) # Clear any log left from a previous run
log_message("=== GWAS processing started ===", log_file) # Record the beginning of processing

# Build SNP-to-gene lookup table
# Find all CSV files in the user-supplied SNP mapping directory
files <- list.files(gene_dir, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE)
log_message(paste("Found", length(files), "gene CSV files."), log_file) # Record the number of mapping files found

# Extract the TAIR10 gene identifier from the portion of each filename preceding the first underscore
extract_gene <- function(f) {
  str_split(basename(f), "_", simplify = TRUE) [1]
}

# Read each gene-specific SNP mapping file
mapping_list <- lapply(files, function(f) {
  gene <- extract_gene(f) # Extract the gene identifier from the filename
  dt <- fread(f, select = c("snp.chr", "snp.position")) # Read only the chromosome and SNP position columns
  dt[, gene := gene] # Add the gene identifier as a new column
  dt # Return the table for this gene
})

# Combine all gene-specific mapping tables
mapping_dt <- rbindlist(mapping_list, use.names = TRUE) # A SNP may appear in multiple rows if it maps to multiple genes

# Create one lookup row per chromosome and SNP position
# If a SNP maps to multiple genes, combine the unique gene identifiers into a semicolon-separated field
mapping_lookup <- mapping_dt[, .(gene = paste(unique(gene), collapse = ";")), by = .(snp.chr, snp.position)]

# OPTIONAL: Intermediate outputs
# These mapping tables are not used in later analysis scripts. Uncomment the following lines only if you want to inspect or retain them.

# write.csv(mapping_dt, file = file.path("results", "mapping_dt.csv"), row.names = FALSE)
# write.csv(mapping_lookup, file = file.path("results", "mapping_lookup.csv"), row.names = FALSE)
# log_message("Saved mapping tables (mapping_dt.csv and mapping_lookup.csv).", log_file)

# Standardize the lookup table before joining
lookup <- mapping_lookup %>%
  rename(chr = snp.chr, pos = snp.position) %>% # Rename columns to match the AraGWAS HDF5 association data
  mutate(chr = gsub("chr", "", chr), pos = as.numeric(pos)) # Remove the "chr" prefix and ensure positions are numeric

# Identify HDF5 files in the AraGWAS archive
# List all files inside the archive without extracting the entire archive
all_files <- untar(tar_file, list = TRUE)
# Retain only files ending in .hdf5
hdf5_files <- grep("\\.hdf5$", all_files, value = TRUE)
# Record the number of phenotype files found
log_message(paste("Found", length(hdf5_files), "HDF5 phenotype files."), log_file)

# Define function to process one HDF5 file
# Use tryCatch so that one problematic file does not stop the entire run
process_one <- function(hdf5_name, lookup, tar_file, temp_dir, file_index = NULL, total_files = NULL) {
  tryCatch({
    # Record which file is being processed
    if(!is.null(file_index) & !is.null(total_files)){
      log_message(paste("Processing file", file_index, "of", total_files, ":", hdf5_name), log_file)
    } else {
      log_message(paste("Processing:", hdf5_name), log_file)
    }
    
    # Extract only the current HDF5 file from the archive
    untar(tar_file, files = hdf5_name, exdir = temp_dir)
    # Construct the path to the extracted HDF5 file
    full_path <- file.path(temp_dir, hdf5_name)
    
    # Identify chromosome groups stored under /pvalues
    chr_groups <- h5ls(full_path) %>%
      filter(group == "/pvalues") %>%
      pull(name)
    
    # Extract association results for each chromosome
    snp_df <- map_dfr(chr_groups, function(chr) {
      chr_path <- paste0("/pvalues/", chr) # Construct the internal HDF5 path for the chromosome
      # Read SNP position, association score, and effect size
      data.frame(
        chr = gsub("chr", "", chr),
        pos = h5read(full_path, paste0(chr_path, "/positions")),
        score = h5read(full_path, paste0(chr_path, "/scores")),
        beta = h5read(full_path, paste0(chr_path, "/beta"))
      )
    }) %>% mutate(chr = as.character(chr), pos = as.numeric(pos)) # Standardize data types before joining
    
    # Match association results to the SNP-to-gene lookup using chromosome and genomic position
    joined <- inner_join(snp_df, lookup, by = c("chr", "pos")) %>%
      mutate(phenotype_file = basename(hdf5_name)) # Record the source phenotype HDF5 file
    
    # Record the number of matched rows
    log_message(paste("Matched", nrow(joined), "SNPs for file", hdf5_name), log_file)
    
    # Delete the extracted file after processing
    if(file.exists(full_path)) file.remove(full_path)
    
    # Return the matched association data
    return(joined)
  }, error = function(e) {
    # Record the filename and error message
    log_message(paste("ERROR processing file", hdf5_name, ":", e$message), log_file)
    return(data.frame()) # Return an empty data frame so loop can continue
  })
}

# Initialize main output file
# Create an empty CSV with expected column names
write.csv(data.frame(chr=character(), pos=numeric(), score=numeric(), beta=numeric(), gene=character(), phenotype_file=character()),
          out_snp_file, row.names = FALSE)

# Process HDF5 files incrementally
# Record the overall start time for calculating the estimated completion
start_time <- Sys.time()

# Process one phenotype file at a time to reduce memory use
for(i in seq_along(hdf5_files)) {
  # Record the start time for the current file
  iter_start <- Sys.time()
  
  # Process the current HDF5 file
  res <- process_one(hdf5_files[i], lookup, tar_file, temp_dir, i, length(hdf5_files))
  if(nrow(res) > 0) {
    # Append matched rows to the main output
    fwrite(res, out_snp_file, append = TRUE)
  }
  
  # Calculate the processing time for the current file
  iter_time <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))
  # Calculate the total elapsed time
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  # Calculate the average time per completed file
  avg_time <- elapsed / i
  # Estimate the remaining processing time
  remaining <- avg_time * (length(hdf5_files) - i)
  # Estimate the completion date and time
  eta <- Sys.time() + remaining
  
  # Record processing progress and estimated completion time  
  log_message(
    paste("Finished file", i, "of", length(hdf5_files),
          "- Time this file:", round(iter_time,1), "s;",
          "ETA:", format(eta, "%Y-%m-%d %H:%M:%S")),
    log_file
  )
}

# Record completion
log_message("Finished processing SNP-level results. Saved to aragwas_combined_results.csv.", log_file)
=