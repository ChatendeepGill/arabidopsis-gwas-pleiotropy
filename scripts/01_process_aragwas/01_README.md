**Purpose:** 
Extract SNP-level association results from the AraGWAS HDF5 archive and map SNPs to TAIR10 genes using gene-specific SNP mapping CSV files downloaded from the AraGWAS Catalog. 

The AraGWAS database and gene-specific CSV files are not included in this repository. Before running this script, download both inputs from the AraGWAS Catalog and update the input paths below. 

The gene-specific CSV files used in this analysis contain SNPs passing the AraGWAS Catalog significance threshold of p < 1 x 10^-4.

**Input:**
- `aragwas_db.tar.gz`: The AraGWAS HDF5 archive is not redistributed with this repository. It must be downloaded from the AraGWAS Catalog.
- `snp_gene_mappings/*.csv`: The gene-specific CSV files are not included in this repository. They must be downloaded from the AraGWAS Catalog. 

The TAIR10 gene identifier is extracted from the portion of the filename preceding the first underscore.

**Main Output:**
- `aragwas_combined_results.csv`: This file contains AraGWAS SNP-level association results that could be mapped to at least one TAIR10 gene. The output includes:

&emsp;&emsp; • `chr`: chromosome

&emsp;&emsp; • `pos`: genomic position

&emsp;&emsp; • `score`: association score reported in the AraGWAS HDF5 file

&emsp;&emsp; • `beta`: estimated SNP effect size

&emsp;&emsp; • `gene`: mapped TAIR10 gene identifier or identifiers

&emsp;&emsp; • `phenotype_file`: source AraGWAS phenotype HDF5 filename

When one SNP maps to multiple genes, gene identifiers are separated by semicolons. 

**Log Output:**
- `01_process_log.txt`: The log file records processing progress, matched SNP counts, estimated completion times, and any HDF5 files that could not be processed. 

**Optional Additional Outputs:**
- `mapping_dt.csv`: This file contains all SNP-gene mapping records before SNPs mapped to multiple genes are collapsed. 
- `mapping_lookup.csv`: This file contains one row per chromosome and genomic position. When a SNP maps to multiple genes, the gene identifiers are combined into a semicolon-separated field. 

The script contains commented code for saving these additional output files. These intermediate tables are not required by later analysis scripts. 

**Notes:**
- The script processes one HDF5 phenotype file at a time and appends its results to the output CSV. This reduces memory use when processing the complete AraGWAS archive. Temporary HDF5 files are removed after processing. 
