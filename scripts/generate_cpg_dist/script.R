#!/usr/bin/env Rscript

# Load necessary libraries first
suppressPackageStartupMessages({
  library(optparse) # Library for parsing command-line flags
  library(data.table)
})

################################################################################
## Argument Parsing
################################################################################

# Define the command-line options (flags)
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to the BED file containing RLM reads.", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Path for the output filtered BED file.", metavar="character")
)

# Create a parser object
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$input) || is.null(opt$output)){
  print_help(opt_parser)
  stop("All arguments (--input and --output) must be supplied.", call.=FALSE)
}

# Validate that input files exist
if (!file.exists(opt$input)) {
  stop("Error: Input RLM file does not exist: ", opt$input, call.=FALSE)
}

################################################################################
## Script Run
################################################################################

cat("Processing files:\n")
cat("  Input RLM file:", opt$input, "\n")
cat("  Output file:", opt$output, "\n")

# Read in the RLM bed file
rlm <- data.table::fread(
  opt$input,
  sep = "\t",
  data.table = TRUE
)

# Check if the RLM input file contains any data rows
if (nrow(rlm) == 0) {
  cat("Error: Input RLM file contains no data rows (only a header).")
  quit(status = 42)
}

# 1. Calculate base statistics per target in one pass
feature_stats <- rlm[, .(
    mean_methylation = mean(mean_methylation, na.rm = TRUE),
    mean_discordance = mean(discordance_score),
    mean_txscore = mean(transitions_score),
    cov = .N,
    chr = chr[1], # Assumes seqnames is constant per target
    start = min(start),
    end = max(end),
    n_CpGs = max(score)
), by = .(featureID)]

# 2. Count the occurrences of each methylation pattern
methyl_counts <- rlm[, .N, by = .(featureID, n_CpGs_methyl)]

# 3. Create a complete template for all possible methylation levels for each target
template <- feature_stats[, .(n_CpGs_methyl = 0:n_CpGs), by = .(featureID)]

# 4. Right-join the counts onto the template to ensure all levels are present
rlm_summary <- methyl_counts[template, on = .(featureID, n_CpGs_methyl)]
rlm_summary[is.na(N), N := 0] # Replace NA counts with 0 for patterns with no reads

# 5. Join the overall feature stats back in and calculate the proportion
rlm_summary <- feature_stats[rlm_summary, on = .(featureID)]
rlm_summary[, prop := N / cov]

# 6. Clean up columns and set the order for the final output
setorder(rlm_summary, featureID, n_CpGs_methyl)

# 7. Select, rename, and reorder columns for the final output
final_output <- rlm_summary[, .(
  chr = chr,
  start = start,
  end = end,
  featureID = featureID,
  n_CpGs = n_CpGs,
  n_CpGs_methyl = n_CpGs_methyl,
  N = N,
  cov = cov,
  prop = prop,
  mean_methylation = mean_methylation,
  mean_discordance = mean_discordance, 
  mean_txscore = mean_txscore
)]

# Check if the RLM input file contains any data rows
if (nrow(final_output) == 0) {
  cat("Error: Output file contains no data rows (only a header).")
  quit(status = 43)
}


# Write outfile
fwrite(
  final_output,
  opt$output, # Use the output path from the flag
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)

cat("Summarized BED file written to:", opt$output, "\n")