# Load packages
library(dplyr)
library(readr)

# Define input
input_file <- file.path(
  "results",
  "gene_trait_collapsed.csv"
)

# Define output
output_dir <- "results"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_file <- file.path(
  output_dir,
  "gene_pleiotropy_hill_metrics.csv"
)

# Load gene-trait effect data 
# The input contains one row per unique gene-trait combination
df <- read_csv(
  input_file,
  show_col_types = FALSE
)

# Convert absolute effect sizes to proportions per gene
# Divide each trait effect by the total effect magnitude for that gene. The resulting proportions describe how the gene's effect magnitude is distributed across its associated traits.
df_props <- df %>%
  group_by(
    gene
  ) %>%
  mutate(
    proportion = abundance / sum(
      abundance,
      na.rm = TRUE
    )
  ) %>%
  ungroup()


# OPTIONAL CHECK: Confirm that trait proportions sum to approximately 1 for each gene.
# Uncomment this code to inspect the result.
# df_props %>%
#   group_by(gene) %>%
#   summarise(
#     sum_p = sum(proportion),
#     .groups = "drop"
#   )

# Calculate gene-level pleiotropy metrics 
# Shannon entropy:
df_hill <- df_props %>%
  group_by(
    gene
  ) %>%
  summarise(
    
    # Number of unique ontology traits associated with the gene.
    trait_richness = n_distinct(
      trait
    ),
    
    # Sum of absolute standardized effect sizes across traits.
    total_effect_magnitude = sum(
      abundance,
      na.rm = TRUE
    ),
    
    # Shannon entropy of the within-gene trait-effect distribution.
    H = -sum(
      proportion * log(proportion),
      na.rm = TRUE
    ),
    
    # Effective number of equally weighted traits.
    hill_q1 = exp(H),
    
    # Sum of squared effect sizes.
    magnitude = sum(
      abundance^2,
      na.rm = TRUE
    ),
    
    # Mean absolute standardized effect size across traits.
    mean_of_means = mean(
      abundance,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

# Save gene-level pleiotropy metrics
write.csv(
  df_hill,
  output_file,
  row.names = FALSE
)
