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
  make_option(c("-t", "--bed_file"), type="character", default=NULL,
              help="Path to the BED file containing target regions.\n\t(requires: chr, started, end, featureID, score)", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to the BED file containing RLM reads.", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Path for the output filtered BED file.", metavar="character"),
  make_option(c("-u", "--unannotated_output"), type="character", default=NULL,
              help="Optional path for the output file of unannotated reads.", metavar="character")
)

# Create a parser object
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$bed_file) || is.null(opt$input) || is.null(opt$output)){
  print_help(opt_parser)
  stop("All three arguments (--bed_file, --input, and --output) must be supplied.", call.=FALSE)
}

# Validate that input files exist
if (!file.exists(opt$bed_file)) {
  stop("Error: Target BED file does not exist: ", opt$bed_file, call.=FALSE)
}
if (!file.exists(opt$input)) {
  stop("Error: Input RLM file does not exist: ", opt$input, call.=FALSE)
}

################################################################################
## Script Run
################################################################################

cat("Processing files:\n")
cat("  Target BED file:", opt$bed_file, "\n")
cat("  Input RLM file:", opt$input, "\n")
cat("  Output file:", opt$output, "\n")

# Read in the regions of interest bed file
targets <- data.table::fread(
  opt$bed_file,
  sep = "\t",
  data.table = TRUE
)
colnames(targets) <- c("chr", "start", "end", "featureID", "score")

# Read in the RLM bed file
rlm <- data.table::fread(
  opt$input,
  sep = "\t",
  data.table = TRUE
)
colnames(rlm) <- gsub("#", "", colnames(rlm)) # remove "#" from the colnames

# Check if the RLM input file contains any data rows
if (nrow(rlm) == 0) {
  cat("Error: Input RLM file contains no data rows (only a header).")
  quit(status = 42)
}

# Check if the BED input file contains any data rows
if (nrow(targets) == 0) {
  cat("Error: Input BED file contains no data rows (only a header).")
  quit(status = 42)
}

# Set keys for the overlap operation
setkey(targets, chr, start, end)
setkey(rlm, chr, start, end)

# Overlap rlm and targets
rlm.sub <- foverlaps(rlm, targets, type = "any", nomatch = NA)

# If the optional flag is provided, write the unannotated reads
if (!is.null(opt$unannotated_output)) {
  cat("  Writing unannotated reads to:", opt$unannotated_output, "\n")
  # Isolate rows that did not overlap (where target columns are NA)
  unannotated <- rlm.sub[is.na(featureID)]

  # Select only the original columns from the input rlm file and write
  cols <- colnames(rlm)
  fwrite(
    unannotated,
    opt$unannotated_output,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
  )
}

cat("Total overlaps found:", nrow(rlm.sub[!is.na(featureID)]), "\n")

# Filter for reads that did overlap, then apply score filter
rlm.sub <- rlm.sub[!is.na(featureID)]
rlm.sub <- rlm.sub[n_CpGs == score]

# Check if the RLM input file contains any data rows
if (nrow(rlm.sub) == 0) {
  cat("Error: Output RLM file contains no data rows (only a header).")
  quit(status = 43)
}

# Write outfile
fwrite(
  rlm.sub,
  opt$output, # Use the output path from the flag
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)

cat("Filtered BED file written to:", opt$output, "\n")