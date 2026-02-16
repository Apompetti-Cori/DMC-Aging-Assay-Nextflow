library(magrittr)
library(tidyverse)
library(vroom)

args = commandArgs(trailingOnly = TRUE)
### Arguments script neeeds
# Rscript allele_frequency.R CpG_context_filename
# cpg_whitelist.bed filename_to_write_summarized_counts.tsv

### Read in files and whitelisted CpGs
vroom(
  args[1],
  delim = '\t',
  skip = 1,
  col_names = c('read_id', 'strand', 'chr', 'start', 'meth_status')
) %>%
  select(-strand) %>%
  unite(chr_base, c(chr, start), sep = '_', remove = FALSE) %>%
  distinct() -> cpg

vroom(
  args[2],
  delim = '\t',
  col_names = c('chr', 'start', 'end', 'target', 'target_cpg_count')
) %>%
  select(-end) %>%
  unite(chr_base, c(chr, start), sep = '_', remove = FALSE) %>%
  distinct() %>%
  filter(target != "C1766") -> cpg_whitelist

### Count the number of methylated CpGs per read
cpg %>%
  left_join(cpg_whitelist, by = 'chr_base') %>%
  na.omit() %>%
  distinct() %>%
  add_count(read_id, target) %>%
  filter(n == target_cpg_count) %>%
  select(-n) %>%
  dplyr::count(read_id, target, target_cpg_count, meth_status) %>%
  pivot_wider(names_from = meth_status, values_from = n) %>%
  replace_na(list(Z = 0, z = 0)) -> read_counts_wide

### Handle cases where the context file does not overlap the target regions
if (nrow(read_counts_wide) == 0) {
  write(
    "Error: CpG context file does not overlap the target regions.",
    stderr()
  )
  quit(save = "no", status = 2)
}

### Handle cases where no reads are methylated in target regions
if (!("Z" %in% names(read_counts_wide))) {
  read_counts_wide$Z = 0
}

read_counts_wide %>%
  select(read_id, target, target_cpg_count, count_meth_cpgs = Z) -> read_counts

### Summarize counts by target and number of methylated CpGs
tbl <- NULL
for (i in unique(cpg_whitelist$target)) {
  rbind(
    tibble(
      target = i,
      count_meth_cpgs = 0:unique(
        filter(cpg_whitelist, target == i)$target_cpg_count
      )
    ),
    tbl
  ) -> tbl
}

read_counts %>%
  add_count(target, name = 'target_reads') %>%
  dplyr::count(target, target_reads, count_meth_cpgs, name = 'read_count') %>%
  mutate(probability = read_count / target_reads) %>%
  full_join(tbl, by = c('target', 'count_meth_cpgs')) %>%
  arrange(target, count_meth_cpgs) %>%
  fill(target_reads) %>%
  replace_na(list(
    read_count = 0,
    probability = 0
  )) %>%
  group_by(target) %>%
  mutate(
    target_reads = sum(read_count, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(target_reads > 0) -> counts


write_tsv(counts, args[3])
