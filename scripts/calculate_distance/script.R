#!/usr/bin/env Rscript

# Load necessary libraries first
suppressPackageStartupMessages({
  library(optparse) # Library for parsing command-line flags
  library(data.table)
  library(philentropy)
})

################################################################################
## Argument Parsing
################################################################################

# Define the command-line options (flags)
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="BED file containing the sample cpg distribution", metavar="character"),
  make_option(c("-r", "--reference"), type="character", default=NULL,
              help="BED file containing the reference cpg distribution", metavar="character"),    
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Path for the output distance file.", metavar="character")
)

# Create a parser object
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$input) || is.null(opt$reference) || is.null(opt$output)){
  print_help(opt_parser)
  stop("All arguments must be supplied.", call.=FALSE)
}

# Validate that input files exist
if (!file.exists(opt$input)) {
  stop("Error: Input BED file does not exist: ", opt$input, call.=FALSE)
}

if (!file.exists(opt$reference)) {
  stop("Error: Reference BED file does not exist: ", opt$reference, call.=FALSE)
}

################################################################################
## Script Run
################################################################################

cat("Processing files:\n")
cat("  Input RLM file:", opt$input, "\n")
cat("  Reference BED file:", opt$reference, "\n")
cat("  Output file:", opt$output, "\n")


# Read in the sample and reference bed files
sample <- fread(
  opt$input,
  sep = "\t"
)

ref <- fread(
  opt$reference,
  sep = "\t"
)

# Inner join on featureID and n_CpGs_methyl
# Then calculate JSD for each target region
# Also calculate coverage, mean methylation, and difference in mean methylation
merge_dt <- ref[sample, on = .(featureID, n_CpGs_methyl), allow.cartesian=TRUE, nomatch=NULL]

jsd_calculation <- merge_dt[,{
  j <- copy(.SD)
  x <- data.frame(
    ref = j$prop,
    sample = j$i.prop
  )
  x <- t(as.matrix(x))

  y <- data.frame(
    ref = j$syn.prop,
    sample = j$i.prop
  )
  y <- t(as.matrix(y))

  list(
    chr = unique(i.chr),
    start = unique(i.start),
    end = unique(i.end),
    cov = unique(j$i.cov),
    mean_methylation = unique(i.mean_methylation),
    diff_methylation = (unique(i.mean_methylation) - unique(mean_methylation)),
    jsd = suppressMessages({sqrt(philentropy::JSD(x))}),
    jsd.syn = suppressMessages({sqrt(philentropy::JSD(y))}),
    mean_txscore = unique(i.mean_txscore),
    mean_discordance = unique(i.mean_discordance)
  )
  
}, by = .(featureID, n_CpGs)]
data.table::setcolorder(jsd_calculation, c("chr", "start", "end"))

# Write the output to a file
fwrite(
  jsd_calculation, 
  opt$output,
  sep = "\t",
  col.names = T,
  row.names = F
)