library(magrittr)

args = commandArgs(trailingOnly = TRUE)
### Arguments script neeeds
# Rscript bs_conversion_efficiency.R CHG_file.txt CHH_file.txt file_to_save

### read in files containing location of non-CpG methylation
# CHG context
rbind(vroom::vroom(
  args[1],
  delim = '\t',
  skip = 1,
  col_names = c('read_id', 'strand', 'chr', 'position', 'meth_status')
)) %>%
  dplyr::mutate(context = 'CHG') -> chg

# CHH context
rbind(
  vroom::vroom(
    args[2],
    delim = '\t',
    skip = 1,
    col_names = c('read_id', 'strand', 'chr', 'position', 'meth_status')
  )
) %>%
  dplyr::mutate(context = 'CHH') -> chh

### count percentage of "methylated" non-CpG cytosines
rbind(chg, chh) %>%
  # distinct() %>%
  dplyr::mutate(
    meth_status = ifelse(meth_status %in% c('x', 'h'), 'unmeth', 'meth')
  ) %>%
  dplyr::count(context, meth_status) %>%
  tidyr::pivot_wider(names_from = meth_status, values_from = n) -> context_wide

if (!("meth" %in% names(context_wide))) {
  context_wide$meth = 0
}

if (!("unmeth" %in% names(context_wide))) {
  context_wide$unmeth = 0
}

context_wide %>%
  tidyr::replace_na(list(meth = 0, unmeth = 0)) %>%
  dplyr::mutate(percent_meth = (meth / (meth + unmeth)) * 100) -> bs_efficiency

# save
readr::write_tsv(bs_efficiency, args[3])
