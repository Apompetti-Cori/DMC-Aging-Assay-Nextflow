library(magrittr)
library(tidyverse)
library(vroom)
library(entropy)

args <- commandArgs(trailingOnly = TRUE)
### Arguments script neeeds
# Rscript calc_summary/script.R reference_dist_filename allele_counts_filename
# bs_efficiency_filename file_to_save

### Read in files and whitelisted CpGs
# reference distribution for calculating JSD
vroom::vroom(args[1]) %>%
  dplyr::select(
    chr,
    start,
    end,
    target = featureID,
    count_meth_cpgs = n_CpGs_methyl,
    prop,
    syn.prop
  ) -> cb_ref

# allele count file
vroom::vroom(args[2]) %>%
  filter(!is.na(count_meth_cpgs)) %>%
  group_by(target) %>%
  # Add in the number of CpGs in the target for average percent meth calculation
  # Also need to recalculate total reads per target since that column is
  # inaccurate in multiple places, but the read counts for the allele frequencies
  # are correct and I can just sum them; can be removed in future
  mutate(
    target_cpg_count = max(count_meth_cpgs),
    target_reads = sum(read_count)
  ) %>%
  ungroup() %>%
  # drop targets missing reads
  na.omit() %>%
  # drop targets without reads
  filter(target_reads > 0) -> data

# bs efficiency file
vroom::vroom(args[3]) %>%
  dplyr::summarize(avg_nonCG_meth = mean(percent_meth)) -> bs_efficiency

### Calculate average percent methylation over the amplicon
data %>%
  # calculate average percent methylation
  group_by(target, target_reads) %>%
  dplyr::reframe(
    avg_perc_meth = (sum((count_meth_cpgs / target_cpg_count) * read_count) /
      target_reads)
  ) %>%
  ungroup() %>%
  distinct() -> perc_meth

### Calculate JSD for the amplicon
data %>%
  ### add in the cord blood reference information
  inner_join(cb_ref, by = c("target", "count_meth_cpgs")) %>%
  ### reformat for JSD function
  select(target, count_meth_cpgs, sample.prop = probability, prop, syn.prop) %>%
  pivot_longer(
    c(sample.prop, prop, syn.prop),
    names_to = 'set',
    values_to = 'prop'
  ) %>%
  pivot_wider(names_from = count_meth_cpgs, values_from = prop) %>%
  select(-set) %>%
  group_by(target) %>%
  nest() %>%
  ungroup() %>%
  ### reshape and use philentropy::JSD() to calculate the JSD
  mutate(
    data = map(data, ~ t(as.matrix(.))),
    data = map(data, ~ t(na.omit(.))),
    data.cord = map(data, ~ .[c(1, 2), ]),
    data.syn = map(data, ~ .[c(1, 3), ]),
    jsd = map(data.cord, ~ sqrt(philentropy::JSD(., unit = "log2"))),
    jsd.syn = map(data.syn, ~ sqrt(philentropy::JSD(., unit = "log2"))),
    entropy = map(
      data.cord,
      ~ {
        p <- na.omit(unlist(.[1, ]))
        philentropy::H(p)
      }
    ),
    entropy.relative = entropy <- map(
      data.cord,
      ~ {
        p <- na.omit(unlist(.[1, ]))
        philentropy::H(p) /
          log2(length(p))
      }
    ),
  ) %>%
  unnest(c(jsd, jsd.syn, entropy, entropy.relative)) %>%
  select(-c(data, data.cord, data.syn)) -> jsd

raw_data <- data %>%
  ### add in the cord blood reference information
  inner_join(cb_ref, by = c("target", "count_meth_cpgs")) %>%
  ### reformat for JSD function
  select(target, count_meth_cpgs, sample.prop = probability, prop, syn.prop) %>%
  pivot_longer(
    c(sample.prop, prop, syn.prop),
    names_to = 'set',
    values_to = 'prop'
  )


saveRDS(raw_data, file = gsub("_meth_jsd.tsv", ".rds", args[4]))

### Join together and save
left_join(perc_meth, jsd, by = "target") %>%
  dplyr::bind_cols(bs_efficiency) %>%
  write_tsv(
    x = .,
    file = args[4]
  )
