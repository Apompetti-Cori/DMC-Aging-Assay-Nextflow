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
  make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "A comma-separated list of paths to the BED files (e.g., 'file1.bed,file2.bed').",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "Path for the output filtered BED file.",
    metavar = "character"
  )
)

# Create a parser object
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$input) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Both --input and --output arguments must be supplied.", call. = FALSE)
}

# Split the comma-separated input string into a vector of file paths
input_files <- strsplit(opt$input, ",")[[1]]
names(input_files) <- input_files

# Validate that all input files exist
sapply(input_files, function(file) {
  if (!file.exists(file)) {
    stop(paste("Input file not found:", file), call. = FALSE)
  }
})

################################################################################
## Data Loading
################################################################################

# Read all specified BED files into a list of data.tables
# The `lapply` function applies the `fread` command to each file path in the `input_files` vector
message("Reading input files...")
list_of_data <- lapply(input_files, fread)
names(list_of_data) <- input_files

# Combine the list of data.tables into a single data.table
# `rbindlist` is a very fast way to row-bind data.tables
combined_data <- rbindlist(list_of_data, idcol = "file_path")

# Extract sample_id from file_path (e.g., "CB1234_cpg_dist.bed" -> "CB1234")
combined_data[, sample_id := sub("_cpg_dist\\.bed$", "", basename(file_path))]

message(paste(
  "Successfully loaded and combined",
  length(input_files),
  "files."
))

################################################################################
## Data Processing
################################################################################

# Step 1: Perform aggregation within each sample for each feature
sample_level_agg <- combined_data[,
  {
    x <- .SD

    # Perform the sub-aggregation to sum N by 'n_CpGs_methyl'
    sub_agg <- x[, .(N = sum(N)), by = .(n_CpGs_methyl)]

    # Find the peak CpG count robustly (handles ties)
    peak_val <- sub_agg$n_CpGs_methyl[which.max(sub_agg$N)]

    # Add additional columns describing the target-level stats for this sample
    sub_agg[, `:=`(
      chr = x$chr[1],
      start = min(x$start),
      end = max(x$end),
      peak = peak_val,
      n_CpGs = n_CpGs[1],
      mean_methylation = mean(x$mean_methylation),
      mean_discordance = mean(x$mean_discordance),
      mean_txscore = mean(x$mean_txscore)
    )]

    sub_agg
  },
  by = .(featureID, sample_id) # Group by both feature and sample
]

# Step 2: Average the results across samples for each feature
ref_cpg_dist <- sample_level_agg[,
  .(
    # Average metrics across samples
    N = as.integer(mean(N)), # Average count of reads
    mean_methylation = mean(mean_methylation),
    mean_discordance = mean(mean_discordance),
    mean_txscore = mean(mean_txscore),
    # Keep feature metadata (they are invariant for the feature)
    chr = chr[1],
    start = start[1],
    end = end[1],
    n_CpGs = n_CpGs[1],
    # Average the sample-specific peak value
    mean_peak = as.integer(mean(peak))
  ),
  by = .(featureID, n_CpGs_methyl) # Group by feature and methylation level
]

ref_cpg_dist[,
  # Final calculations for proportion and coverage
  `:=`(
    # Recalculate total coverage (sum N across n_CpGs_methyl for the feature)
    cov = sum(N),

    # Recalculate proportion based on averaged N
    prop = N / sum(N),

    # Calculate syn.prop using the *mean* peak value
    syn.prop = fcase(
      mean_peak < median(n_CpGs_methyl) & n_CpGs_methyl == 0                  , 1L ,
      mean_peak > median(n_CpGs_methyl) & n_CpGs_methyl == max(n_CpGs_methyl) , 1L ,
      default = 0L
    )
  ),
  by = .(featureID)
]

# Select, rename, and reorder columns for the final output
final_output <- ref_cpg_dist[, .(
  chr = chr,
  start = start,
  end = end,
  featureID = featureID,
  n_CpGs = n_CpGs,
  n_CpGs_methyl = n_CpGs_methyl,
  N = N,
  cov = cov,
  prop = prop,
  syn.prop,
  mean_methylation = mean_methylation,
  mean_discordance = mean_discordance,
  mean_txscore = mean_txscore
)]

# Write outfile
fwrite(
  final_output,
  opt$output, # Use the output path from the flag
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)

cat("Reference file written to:", opt$output, "\n")
