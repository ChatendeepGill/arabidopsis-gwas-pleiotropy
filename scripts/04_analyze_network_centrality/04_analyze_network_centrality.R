# Load packages
library(dplyr)
library(readr)

# Define inputs
# Co-expression network statistics
coexp_file <- file.path(
  "data",
  "external",
  "coexpression",
  "arabidopsis_coexpression_metrics.rds"
)

# Gene-level pleiotropy metrics
hill_file <- file.path(
  "results",
  "gene_pleiotropy_hill_metrics.csv"
)

# Define outputs
# Merged gene-level dataset used for the correlation analyses
out_merged_file <- file.path(
  "results",
  "merged_hill_coexp.csv"
)

# Create output directory
dir.create(
  "results",
  showWarnings = FALSE
)

# Load data 
coexp <- readRDS(
  coexp_file
)

hill <- read_csv(
  hill_file,
  show_col_types = FALSE
)

# Standardize gene identifiers
# Rename the TAIR10 gene identifier column so that all datasets use the same column name.
colnames(coexp)[colnames(coexp) == "TAIR_gene"] <- "gene"

# Retain genes with a valid, positive Hill q = 1 value
hill <- hill %>%
  filter(
    !is.na(hill_q1),
    hill_q1 > 0
  )

# Merge datasets
# Retain genes represented in both datasets.
merged <- hill %>%
  select(
    gene,
    trait_richness,
    total_effect_magnitude,
    H,
    hill_q1,
    magnitude,
    mean_of_means
  ) %>%
  inner_join(
    coexp %>%
      select(
        gene,
        node_degree,
        node_strength,
        node_closeness,
        node_betweenness
      ),
    by = "gene"
  ) %>%
  
  # Retain complete observations for all reported correlations.
  filter(
    !is.na(node_degree),
    !is.na(node_strength),
    !is.na(node_closeness),
    !is.na(node_betweenness)
  )

# Save merged dataset
write.csv(
  merged,
  out_merged_file,
  row.names = FALSE
)

# Calculate Spearman correlations
# Hill q = 1 versus co-expression network degree.
cor_hill_degree <- cor.test(
  merged$hill_q1,
  merged$node_degree,
  method = "spearman"
)

cor_hill_degree

# Hill q = 1 versus co-expression network strength.
cor_hill_strength <- cor.test(
  merged$hill_q1,
  merged$node_strength,
  method = "spearman"
)

cor_hill_strength

# Hill q = 1 versus co-expression network closeness.
cor_hill_closeness <- cor.test(
  merged$hill_q1,
  merged$node_closeness,
  method = "spearman"
)

cor_hill_closeness

# Hill q = 1 versus co-expression network betweenness.
cor_hill_betweenness <- cor.test(
  merged$hill_q1,
  merged$node_betweenness,
  method = "spearman"
)

cor_hill_betweenness
