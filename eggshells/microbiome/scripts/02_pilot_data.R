# this script is for working with the pilot data
# Author: Anya Mueller
# Date: August 11, 2023

#### Set up ####
#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("rstatix")
#install.packages("janitor")
#install.packages("taxize")
#install.packages("eulerr")
#install.packages("readxl")
#install.packages("compositions")
#install.packages("emmeans")
#install.packages("BacDive", repos="http://R-Forge.R-project.org")
#install.packages("indicspecies")
#install.packages("gt")
#install.packages("devtools")
#library("devtools")
#install_github("microbiome/microbiome")
#install.packages("gtsummary")
#install.packages("ggdendro")
#install.packages("here")

## read in packages
library(tidyverse)
library(vegan)
library(rstatix)
library(janitor)
library(taxize)
library(eulerr)
library(readxl)
library(compositions)
library(emmeans)
library(BacDive)
library(indicspecies)
library(gt)
library(microbiome)
library(ggpubr)
library(gtsummary)
library(ggdendro)
library(here)
options(scipen = 999) #no scientific notation in numbers
set.seed(56739)

#read in data
pilot_asv_seq <- read_csv(here("eggshells", 
                               "microbiome",
                               "processed_data",
                               "pilot_asv_seq.csv"))

pilot_taxo <- read_csv(here("eggshells", 
                            "microbiome",
                            "processed_data",
                            "pilot_taxo.csv"))

pilot_sample_asv <- read_csv(here("eggshells", 
                                  "microbiome",
                                  "processed_data",
                                  "pilot_sample_asv.csv"))

pilot_swab_nests <- read_csv(here("eggshells", 
                                  "microbiome",
                                  "processed_data",
                                  "swab_nests.csv"),
                             guess_max = 2652) %>%
  filter(nest_id %in% c(pilot_sample_asv %>% 
                          pull(nest_id))) %>%
  select(-host_species) #fix dummy data for example purposes

pilot_metadata <- read_csv(here("eggshells", 
                                "microbiome",
                                "processed_data",
                                "pilot_metadata.csv"))

distance_mat <- read_csv(here("eggshells", 
                              "microbiome",
                              "processed_data",
                              "distance_between_nests_m.csv"))

#### data prep ####
#goal: prepare data for analysis

#identify the samples analysed, filter to nests of interest
pilot_swab_nests <- pilot_metadata %>%
  select(nest_id,
         date) %>%
  mutate(swab_analyzed = "yes") %>%
  full_join(pilot_swab_nests,
            .,
            by = c("nest_id",
                   "date"))

pilot_swab_nests <- pilot_metadata %>%
  select(nest_id,
         host_species) %>%
  full_join(pilot_swab_nests,
            .)

#they are all first swabs, filter to ones that are on first or second egg
pilot_swab_nests_reps <- pilot_swab_nests %>%
  filter(nest_id %in% c(pilot_swab_nests %>%
                          filter(swab_analyzed == "yes",
                                 num_host_egg < 3) %>%
                          pull(nest_id))) 
#make vectors of nests for data comparisons
yewa_nests <- pilot_swab_nests_reps %>%
  filter(host_species == "YEWA") %>%
  pull(nest_id) %>%
  unique()
dwf_nests <- pilot_swab_nests_reps %>%
  filter(site == "dwf") %>%
  pull(nest_id) %>%
  unique() 

#filter asv to nests of interest
pilot_sample_asv_reps <- pilot_sample_asv %>%
  filter(!is.na(nest_id)) %>%
  filter(nest_id %in% c(pilot_swab_nests_reps %>%
                          pull(nest_id) %>%
                          unique()))
#get rid of ASVS with no occurrences
pilot_sample_asv_reps <- pilot_sample_asv_reps %>%
  group_by(asv_name) %>%
  filter(sum(reads) > 0) %>%
  ungroup()

#filter taxonomy to ASVs from nests of interest
pilot_taxo_reps <- pilot_taxo %>% 
  filter(asv_name %in% c(pilot_sample_asv_reps %>%
                           pull(asv_name) %>%
                           unique()))

#how many reads do we have in each sample? (sampling depth)
pilot_sample_asv_reps %>% 
  full_join(.,
            pilot_swab_nests_reps %>%
              select(nest_id,
                     host_species,
                     site) %>%
              unique()) %>%
  unite(col = sample_group,
        site,
        host_species) %>%
  group_by(nest_id,
           sample_group) %>%
  summarize(read_depth = sum(reads)) %>%
  arrange(read_depth) %>%
  ggplot(data = .) +
  geom_bar(aes(y = read_depth,
               x = reorder(nest_id, 
                           -read_depth),
               fill = sample_group),
           stat = "identity") +
  labs(x = "Nest id",
       y = "Reads",
       title = "Sampling depth per nest",
       fill = "Sample group") +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 10))

ggsave(here("eggshells", 
              "microbiome",
              "outputs",
              "sampling_depth.jpeg"))

#make asv*sample matrix
asv_sample <- pilot_sample_asv_reps %>%
  select(asv_name,
         nest_id,
         reads) %>%
  pivot_wider(names_from = "asv_name",
              values_from = "reads") %>%
  column_to_rownames(var = "nest_id")

#transform raw read data in various ways:
#rarefy - https://www.youtube.com/watch?v=c7H8jLjTxSE
#are there any natural breaks in the data?
pilot_sample_asv_reps %>%
  group_by(nest_id) %>%
  summarize(sum_reads = sum(reads)) %>%
  ungroup() %>%
  ggplot(data = .) +
  geom_histogram(aes(x = sum_reads),
                 bins = 10)
#there is a break at 50000, so we could make a cut there
#but how well covered are all our samples?
#goods_coverage = 100(1-num otu w 1 seq/ total number seq)
pilot_sample_asv_reps %>% 
  left_join(.,
            pilot_swab_nests %>%
              select(nest_id,
                     site,
                     host_species) %>%
              unite(col = sample_group,
                    site,
                    host_species,
                    sep = "; ") %>%
              distinct()) %>%
  group_by(nest_id,
           sample_group) %>%
  summarize(sum_reads = sum(reads),
            num_sigs = sum(reads == 1),
            goods = 100*(1-num_sigs/sum_reads)) %>%
  ungroup() %>% 
  ggplot(data = .,
         aes(x = sum_reads,
             y = goods,
             colour = sample_group)) +
  geom_point(size = 3) +
  labs(x = "Number of reads",
       y = "Goods Coverage",
       title = "Goods Coverage of each sample",
       colour = "") +
  theme(text = element_text(size = 20))

ggsave(here("eggshells", 
            "microbiome",
            "outputs",
            "goods_coverage.jpeg"))
#all our samples have very high coverage!

#get min sampling depth
#ASV
min_reads <- pilot_sample_asv_reps %>% 
  group_by(nest_id) %>%
  summarize(sum_reads = sum(reads)) %>%
  pull(sum_reads) %>%
  min(.)

#check out how our samples are looking at our minimum read sampling size
rarefaction_curves <- rarecurve(x = asv_sample,
                                step = 100,
                                sample = min_reads,
                                ylab = "Amplicon sequence variants",
                                xlab = "Number of reads")
#make rarefaction curve for samples and sample groups
rarefaction_curves <- rarecurve(x = asv_sample,
                                step = 100,
                                tidy = TRUE,
                                ylab = "Amplicon sequence variants",
                                xlab = "Number of reads")
rarefaction_curves %>% 
  rename("nest_id" = "Site",
         "asv" = "Species",
         "reads" = "Sample") %>%
  left_join(.,
            pilot_swab_nests %>%
              select(nest_id,
                     site,
                     host_species) %>%
              unite(col = sample_group,
                    site,
                    host_species,
                    sep = "; ") %>%
              distinct()) %>%
  ggplot(data = .,
         aes(x = reads,
             y = asv,
             group = nest_id
             )) +
  geom_line() +
  geom_vline(xintercept=min_reads, linetype='dashed') + 
  facet_wrap(vars(sample_group)) +
  labs(y = "Amplicon sequence variants",
       x = "Number of reads",
       title = "Rarefaction curves of each sample group") +
  theme(text = element_text(size = 15),
        plot.margin = margin(15,15,15,15),
        panel.spacing = unit(20, "points"))

ggsave(here("eggshells", 
            "microbiome",
            "outputs",
            "rare_curves.jpeg"))
  
#the samples are saturated at 16209 reads, so we can rarefy to that level

#rarefy in loop so have a bunch of random subsets you can do work on

j <- c(1:100)
asv_sample_rared <- list()

for (i in j) {
  
  asv_sample_rared[[i]] <- rrarefy(x = asv_sample,
                                 sample = min_reads) %>% 
    as.data.frame()
  
}
asv_sample_rared_mean <- asv_sample_rared %>% 
  map(.x = .,
      .f = rownames_to_column,
      var = "nest_id") %>%
  map(.x = .,
      .f = pivot_longer,
      cols = starts_with("ASV"),
      names_to = "asv_name",
      values_to = "rare_reads") %>%
  reduce(.x = .,
         .f = full_join) %>%
  group_by(nest_id,
           asv_name) %>%
  summarize(mean_rare_reads = mean(rare_reads))

asv_sample_rared_mean <- asv_sample_rared_mean %>% 
  left_join(x = .,
            y = pilot_taxo_reps)

j <- c("asv_name",
       "kingdom",
       "phylum",
       "class",
       "order",
       "family",
       "genus",
       "species")
asv_sample_rared_mean_taxo <- list()
for (i in j) {
  asv_sample_rared_mean_taxo[[i]] <- asv_sample_rared_mean %>%
    group_by(nest_id,
             .data[[i]]) %>%
    summarize(mean_rare_reads = sum(mean_rare_reads),
              .groups = "keep") %>%
    ungroup()
}

#compute relative reads
rel_abund <- asv_sample %>%
  rownames_to_column(var = "nest_id") %>%
  pivot_longer(cols = starts_with("ASV"),
               names_to = "asv_name",
               values_to = "reads") %>%
  group_by(nest_id) %>%
  mutate(rel_abund = reads/sum(reads)) %>%
  ungroup()

rel_abund <- rel_abund %>% 
  left_join(x = .,
            y = pilot_taxo_reps)

j <- c("asv_name",
       "kingdom",
       "phylum",
       "class",
       "order",
       "family",
       "genus",
       "species")
rel_abund_taxo <- list()
for (i in j) {
  rel_abund_taxo[[i]] <- rel_abund %>%
    group_by(nest_id,
             .data[[i]]) %>%
    summarize(rel_abund = sum(rel_abund),
              .groups = "keep") %>%
    ungroup()
  #check that summing to 1
  rel_abund_taxo[[i]] %>%
    group_by(nest_id) %>%
    summarize(sumra = sum(rel_abund)) %>%
    print()
}

# filter distance matrix for examined nests and convert to matrix
distance_mat <- distance_mat %>%
  pivot_longer(cols = -nest_id,
               names_to = "nest_id_to",
               values_to = "geo_dist") %>%
  mutate(nest_id_to = str_remove(string = nest_id_to,
                                 pattern = "distance_to_ ")) %>% 
  filter(nest_id %in% c(pilot_swab_nests_reps %>%
                          pull(nest_id) %>%
                          unique()),
         nest_id_to %in% c(pilot_swab_nests_reps %>%
                             pull(nest_id) %>%
                             unique())) %>%
  pivot_wider(names_from = "nest_id_to",
              values_from = "geo_dist") %>%
  column_to_rownames(var = "nest_id") %>%
  as.matrix()
  
geo_distance_long <- distance_mat %>%
  replace_upper_triangle(x =.,
                         by = NA_character_,
                         diagonal = FALSE) %>% 
  as.data.frame() %>%
  rename("nest_id" = "rowname") %>%
  pivot_longer(cols = -nest_id,
               values_to = "geo_distance",
               names_to = "nest_id_compared") %>%
  drop_na()

#### alpha diversity ####

j <- c(1:100)
yewa_asv_sample_rared_alpha <- list()
dwf_asv_sample_rared_alpha <- list()
for (i in j) {
  yewa_df <- asv_sample_rared[[i]] %>% 
    rownames_to_column(var = "nest_id") %>% 
    filter(nest_id %in% c(yewa_nests)) %>%
    column_to_rownames(var = "nest_id")
  
  yewa_asv_sample_rared_alpha[[i]] <- list(estimateR(yewa_df) %>%
                                        t() %>%
                                        as.data.frame() %>%
                                        rownames_to_column(var = "nest_id") %>%
                                        rename("richness" = "S.obs",
                                               "chao1" = "S.chao1") %>%
                                        rename_with(.cols = - nest_id,
                                                    .fn = ~paste(.,
                                                                  i,
                                                                 sep = "_")),
                                      vegan::diversity(x = yewa_df,
                                                index = "shannon") %>%
                                        as.data.frame() %>%
                                        rename("shannon" = ".") %>%
                                        rownames_to_column(var = "nest_id") %>%
                                        rename_with(.cols = - nest_id,
                                                    .fn = ~paste(.,
                                                                 i,
                                                                 sep = "_")),
                                      vegan::diversity(x = yewa_df,
                                                index = "simpson") %>% 
                                        as.data.frame() %>%
                                        rename("simpson" = ".") %>%
                                        rownames_to_column(var = "nest_id") %>%
                                        rename_with(.cols = - nest_id,
                                                    .fn = ~paste(.,
                                                                 i,
                                                                 sep = "_")),
                                      tibble(nest_id = c(yewa_df %>% 
                                                           rownames()),
                                             evenness = vegan::diversity(x = yewa_df,
                                                                 index = "shannon") / 
                                               log(specnumber(x = yewa_df))) %>%
                                        rename_with(.cols = - nest_id,
                                                    .fn = ~paste(.,
                                                                 i,
                                                                 sep = "_"))
                                      ) %>%
    reduce(.x = .,
           .f = full_join)
  
  
  dwf_df <- asv_sample_rared[[i]] %>% 
    rownames_to_column(var = "nest_id") %>% 
    filter(nest_id %in% c(dwf_nests)) %>%
    column_to_rownames(var = "nest_id")
  
  dwf_asv_sample_rared_alpha[[i]] <- list(estimateR(dwf_df) %>%
                                             t() %>%
                                             as.data.frame() %>%
                                             rownames_to_column(var = "nest_id") %>%
                                             rename("richness" = "S.obs",
                                                    "chao1" = "S.chao1") %>%
                                             rename_with(.cols = - nest_id,
                                                         .fn = ~paste(.,
                                                                      i,
                                                                      sep = "_")),
                                           vegan::diversity(x = dwf_df,
                                                            index = "shannon") %>%
                                             as.data.frame() %>%
                                             rename("shannon" = ".") %>%
                                             rownames_to_column(var = "nest_id") %>%
                                             rename_with(.cols = - nest_id,
                                                         .fn = ~paste(.,
                                                                      i,
                                                                      sep = "_")),
                                           vegan::diversity(x = dwf_df,
                                                            index = "simpson") %>% 
                                             as.data.frame() %>%
                                             rename("simpson" = ".") %>%
                                             rownames_to_column(var = "nest_id") %>%
                                             rename_with(.cols = - nest_id,
                                                         .fn = ~paste(.,
                                                                      i,
                                                                      sep = "_")),
                                           tibble(nest_id = c(dwf_df %>% 
                                                                rownames()),
                                                  evenness = vegan::diversity(x = dwf_df,
                                                                              index = "shannon") / 
                                                    log(specnumber(x = dwf_df))) %>%
                                             rename_with(.cols = - nest_id,
                                                         .fn = ~paste(.,
                                                                      i,
                                                                      sep = "_"))
  ) %>%
    reduce(.x = .,
           .f = full_join)
}
#take mean of each measure
yewa_asv_sample_rared_alpha_sumd <- reduce(.x = yewa_asv_sample_rared_alpha,
                                      .f = full_join) %>%
  pivot_longer(cols = - nest_id,
               names_to = "a_div_type",
               values_to = "a_div") %>%
  mutate(a_div_type = str_remove(string = a_div_type,
                                 pattern = "_\\d+")) %>%
  group_by(nest_id,
           a_div_type) %>%
  summarize(mean_a_div = mean(a_div),
            median_a_div = median(a_div),
            sd_a_div = sd(a_div),
            min_a_div = min(a_div),
            max_a_div = max(a_div))

dwf_asv_sample_rared_alpha_sumd <- reduce(.x = dwf_asv_sample_rared_alpha,
                                           .f = full_join) %>%
  pivot_longer(cols = - nest_id,
               names_to = "a_div_type",
               values_to = "a_div") %>%
  mutate(a_div_type = str_remove(string = a_div_type,
                                 pattern = "_\\d+")) %>%
  group_by(nest_id,
           a_div_type) %>%
  summarize(mean_a_div = mean(a_div),
            median_a_div = median(a_div),
            sd_a_div = sd(a_div),
            min_a_div = min(a_div),
            max_a_div = max(a_div))

#### beta diversity ####
#take distance
#of all
asv_sample_rared_bray <- avgdist(x = asv_sample,
                                 sample = min_reads,
                                 meanfun = mean,
                                 iterations = 100,
                                 dmethod = "bray")

asv_sample_rared_jacc <- avgdist(x = asv_sample,
                                 sample = min_reads,
                                 meanfun = mean,
                                 iterations = 100,
                                 dmethod = "jaccard",
                                 binary = TRUE) 

#of data subsets

yewa_asv_sample <- asv_sample %>%
  rownames_to_column(var = "nest_id") %>%
  filter(nest_id %in% c(yewa_nests)) %>%
  column_to_rownames(var = "nest_id")

dwf_asv_sample <- asv_sample %>%
  rownames_to_column(var = "nest_id") %>%
  filter(nest_id %in% c(dwf_nests)) %>%
  column_to_rownames(var = "nest_id")

yewa_rared_bray <- avgdist(x = yewa_asv_sample,
                                 sample = min_reads,
                                 meanfun = mean,
                                 iterations = 100,
                                 dmethod = "bray")

dwf_rared_bray <- avgdist(x = dwf_asv_sample,
                          sample = min_reads,
                          meanfun = mean,
                          iterations = 100,
                          dmethod = "bray")
  
yewa_rared_jacc <- avgdist(x = yewa_asv_sample,
                                 sample = min_reads,
                                 meanfun = mean,
                                 iterations = 100,
                                 dmethod = "jaccard",
                                 binary = TRUE) #pa beta

dwf_rared_jacc <- avgdist(x = dwf_asv_sample,
                           sample = min_reads,
                           meanfun = mean,
                           iterations = 100,
                           dmethod = "jaccard",
                           binary = TRUE) #pa beta

####  Summary data ####
#how many reads in total?
pilot_sample_asv_reps %>% 
  summarise(n_reads = sum(reads))
#how many nest samples?
pilot_sample_asv_reps %>% 
  select(nest_id) %>%
  distinct() %>%
  summarise(n = n())
#how many nest samples in each group?
pilot_swab_nests_reps %>%
  group_by(site,
           host_species) %>%
  select(nest_id, 
         site,
         host_species) %>%
  distinct() %>%
  summarise(n = n())
#how many of each taxonomy?
j <- colnames(pilot_taxo_reps)
sumed <- list()
for (i in j) {
  sumed[[i]] <- pilot_taxo_reps %>% 
    filter(.data[[i]] != "Unknown Bacteria") %>%
    select(all_of(i)) %>%
    distinct() %>%
    summarise(n = n())
}
sumed

#### Differential abundance analysis ####
#https://cran.r-project.org/web/packages/indicspecies/vignettes/IndicatorSpeciesAnalysis.html
#do with rarefied data (reads as proxy for abundance) and relative abundance of raw reads

#need community data matrix
j <- c("asv_name",
       "kingdom",
       "phylum",
       "class",
       "order",
       "family",
       "genus",
       "species")
yewa_mean_rared <- list()
yewa_rel_abund <- list()
dwf_mean_rared <- list()
dwf_rel_abund <- list()
for (i in j) {
  yewa_mean_rared[[i]] <- asv_sample_rared_mean_taxo[[i]] %>%
    filter(nest_id %in% yewa_nests) %>%
    pivot_wider(names_from = all_of(i),
                values_from = "mean_rare_reads") %>%
    column_to_rownames(var = "nest_id")
    
  yewa_rel_abund[[i]] <- rel_abund_taxo[[i]] %>%
    filter(nest_id %in% yewa_nests) %>%
    pivot_wider(names_from = all_of(i),
                values_from = "rel_abund") %>%
    column_to_rownames(var = "nest_id")
  
  dwf_mean_rared[[i]] <- asv_sample_rared_mean_taxo[[i]] %>%
    filter(nest_id %in% dwf_nests) %>%
    pivot_wider(names_from = all_of(i),
                values_from = "mean_rare_reads") %>%
    column_to_rownames(var = "nest_id")
  
  dwf_rel_abund[[i]] <- rel_abund_taxo[[i]] %>%
    filter(nest_id %in% dwf_nests) %>%
    pivot_wider(names_from = all_of(i),
                values_from = "rel_abund") %>%
    column_to_rownames(var = "nest_id")
}

#need matrix splitting data into groups
j <- c("asv_name",
       "kingdom",
       "phylum",
       "class",
       "order",
       "family",
       "genus",
       "species")
yewa_mean_rared_groups <- list()
yewa_rel_abund_groups <- list()
dwf_mean_rared_groups <- list()
dwf_rel_abund_groups <- list()
for (i in j) {
  yewa_mean_rared_groups[[i]] <- yewa_mean_rared[[i]] %>%
    rownames_to_column(var = "nest_id") %>% 
    select(nest_id) %>%
    left_join(x = .,
              y = pilot_swab_nests_reps %>%
                select(nest_id,
                       site) %>%
                distinct()) %>%
    pull(site)
  
  yewa_rel_abund_groups[[i]] <- yewa_rel_abund[[i]] %>%
    rownames_to_column(var = "nest_id") %>% 
    select(nest_id) %>%
    left_join(x = .,
              y = pilot_swab_nests_reps %>%
                select(nest_id,
                       site) %>%
                distinct()) %>%
    pull(site)
  
  dwf_mean_rared_groups[[i]] <- dwf_mean_rared[[i]] %>%
    rownames_to_column(var = "nest_id") %>% 
    select(nest_id) %>%
    left_join(x = .,
              y = pilot_swab_nests_reps %>%
                select(nest_id,
                       host_species) %>%
                distinct()) %>%
    pull(host_species)
  
  dwf_rel_abund_groups[[i]] <- dwf_rel_abund[[i]] %>%
    rownames_to_column(var = "nest_id") %>% 
    select(nest_id) %>%
    left_join(x = .,
              y = pilot_swab_nests_reps %>%
                select(nest_id,
                       host_species) %>%
                distinct()) %>%
    pull(host_species)
}

#run test
j <- c("asv_name",
       "phylum",
       "class",
       "order",
       "family",
       "genus",
       "species")
yewa_mean_rared_indval <- list()
yewa_rel_abund_indval <- list()
dwf_mean_rared_indval <- list()
dwf_rel_abund_indval <- list()
for (i in j) {
  yewa_mean_rared_indval[[i]] <- multipatt(x = yewa_mean_rared[[i]],
                                           cluster = yewa_mean_rared_groups[[i]],
                                           control = how(nperm=999)) 
  print("yewa_mean_rared_indval")
  print(i)
  summary(yewa_mean_rared_indval[[i]],
          indvalcomp = TRUE)
  
  yewa_rel_abund_indval[[i]] <- multipatt(x = yewa_rel_abund[[i]],
                                          cluster = yewa_rel_abund_groups[[i]],
                                          control = how(nperm=999))
  print("yewa_rel_abund_indval")
  print(i)
  summary(yewa_rel_abund_indval[[i]],
          indvalcomp = TRUE)
  
  dwf_mean_rared_indval[[i]] <- multipatt(x = dwf_mean_rared[[i]],
                                          cluster = dwf_mean_rared_groups[[i]],
                                          control = how(nperm=999))
  print("dwf_mean_rared_indval")
  print(i)
  summary(dwf_mean_rared_indval[[i]],
          indvalcomp = TRUE)
  
  dwf_rel_abund_indval[[i]] <- multipatt(x = dwf_rel_abund[[i]],
                                         cluster = dwf_rel_abund_groups[[i]],
                                         control = how(nperm=999)) 
  print("dwf_rel_abund_indval")
  print(i)
  summary(dwf_rel_abund_indval[[i]],
          indvalcomp = TRUE)
}
#compA (specificity) is sample estimate of the probability that the surveyed site belongs to the target site group given the fact that the species has been found
# 1 = it occurs in sites belonging to this group only
#compB (fidelity) is sample estimate of the probability of finding the species in sites belonging to the site group
# 1 = it appears in all sites belonging to this group
dwf_rel_abund_indval[[1]]$sign %>%
  rownames_to_column(var = "asv_name") %>%
  filter(p.value <= 0.05) %>%
  clean_names() %>%
  mutate(indicator_for = case_when(s_rwbl == 0 & s_yewa == 1 ~ "YEWA",
                                   s_rwbl == 1 & s_yewa == 0 ~ "RWBL")) %>%
  select(asv_name,
         indicator_for,
         stat,
         p_value) %>%
  left_join(.,
            dwf_rel_abund[[1]] %>%
              rownames_to_column(var = "nest_id") %>%
              pivot_longer(cols = -nest_id,
                           names_to = "asv_name",
                           values_to = "rel_abund") %>%
              left_join(.,
                        pilot_swab_nests_reps %>%
                          select(nest_id,
                                 host_species) %>%
                          distinct()) %>%
              group_by(host_species,
                       asv_name) %>%
              summarize(mean_rel_abund = mean(rel_abund))
            )

#extract data and combine with summary data to get a sense of the importance of these taxa
j <- c("asv_name",
       "phylum",
       "class",
       "order",
       "family",
       "genus",
       "species")
yewa_mean_rared_indval_sum <- list()
yewa_rel_abund_indval_sum <- list()
dwf_mean_rared_indval_sum <- list()
dwf_rel_abund_indval_sum <- list()
for (i in j) {
  yewa_mean_rared_indval_sum[[i]] <- yewa_mean_rared_indval[[i]]$sign %>%
    rownames_to_column(var = i) %>%
    filter(p.value <= 0.05) %>%
    clean_names() %>%
    mutate(indicator_for = case_when(s_pcc == 0 & s_dwf == 1 ~ "dwf",
                                     s_pcc == 1 & s_dwf == 0 ~ "pcc")) %>%
    select(all_of(i),
           indicator_for,
           stat,
           p_value) %>%
    left_join(.,
              yewa_mean_rared[[i]] %>%
                rownames_to_column(var = "nest_id") %>%
                pivot_longer(cols = -nest_id,
                             names_to = i,
                             values_to = "rared") %>%
                left_join(.,
                          pilot_swab_nests_reps %>%
                            select(nest_id,
                                   site) %>%
                            distinct()) %>%
                group_by(site,
                         .data[[i]]) %>%
                summarize(mean_rared = mean(rared))
    ) %>%
    pivot_wider(names_from = site,
                values_from = mean_rared) %>%
    rename("taxon" = all_of(i)) %>%
    mutate(taxon_level = i,
           .before = taxon)
  
    yewa_rel_abund_indval_sum[[i]] <- yewa_rel_abund_indval[[i]]$sign %>%
    rownames_to_column(var = i) %>%
    filter(p.value <= 0.05) %>%
    clean_names() %>%
    mutate(indicator_for = case_when(s_pcc == 0 & s_dwf == 1 ~ "dwf",
                                     s_pcc == 1 & s_dwf == 0 ~ "pcc")) %>%
    select(all_of(i),
           indicator_for,
           stat,
           p_value) %>%
    left_join(.,
              yewa_rel_abund[[i]] %>%
                rownames_to_column(var = "nest_id") %>%
                pivot_longer(cols = -nest_id,
                             names_to = i,
                             values_to = "rel_abund") %>%
                left_join(.,
                          pilot_swab_nests_reps %>%
                            select(nest_id,
                                   site) %>%
                            distinct()) %>%
                group_by(site,
                         .data[[i]]) %>%
                summarize(mean_rel_abund = mean(rel_abund))
    ) %>%
      pivot_wider(names_from = site,
                  values_from = mean_rel_abund) %>%
      rename("taxon" = all_of(i)) %>%
      mutate(taxon_level = i,
             .before = taxon)
  
  dwf_mean_rared_indval_sum[[i]] <- dwf_mean_rared_indval[[i]]$sign %>%
    rownames_to_column(var = i) %>%
    filter(p.value <= 0.05) %>%
    clean_names() %>%
    mutate(indicator_for = case_when(s_rwbl == 0 & s_yewa == 1 ~ "YEWA",
                                     s_rwbl == 1 & s_yewa == 0 ~ "RWBL")) %>%
    select(all_of(i),
           indicator_for,
           stat,
           p_value) %>%
    left_join(.,
              dwf_mean_rared[[i]] %>%
                rownames_to_column(var = "nest_id") %>%
                pivot_longer(cols = -nest_id,
                             names_to = i,
                             values_to = "rared") %>%
                left_join(.,
                          pilot_swab_nests_reps %>%
                            select(nest_id,
                                   host_species) %>%
                            distinct()) %>%
                group_by(host_species,
                         .data[[i]]) %>%
                summarize(mean_rared = mean(rared))
    ) %>%
    pivot_wider(names_from = host_species,
                values_from = mean_rared) %>%
    rename("taxon" = all_of(i)) %>%
    mutate(taxon_level = i,
           .before = taxon)
  
  dwf_rel_abund_indval_sum[[i]] <- dwf_rel_abund_indval[[i]]$sign %>%
    rownames_to_column(var = i) %>%
    filter(p.value <= 0.05) %>%
    clean_names() %>%
    mutate(indicator_for = case_when(s_rwbl == 0 & s_yewa == 1 ~ "YEWA",
                                     s_rwbl == 1 & s_yewa == 0 ~ "RWBL")) %>%
    select(all_of(i),
           indicator_for,
           stat,
           p_value) %>%
    left_join(.,
              dwf_rel_abund[[i]] %>%
                rownames_to_column(var = "nest_id") %>%
                pivot_longer(cols = -nest_id,
                             names_to = i,
                             values_to = "rel_abund") %>%
                left_join(.,
                          pilot_swab_nests_reps %>%
                            select(nest_id,
                                   host_species) %>%
                            distinct()) %>%
                group_by(host_species,
                         .data[[i]]) %>%
                summarize(mean_rel_abund = mean(rel_abund))
    ) %>%
    pivot_wider(names_from = host_species,
                values_from = mean_rel_abund) %>%
    rename("taxon" = all_of(i)) %>%
    mutate(taxon_level = i,
           .before = taxon)
}

#make output tables for supplementary data
indicators_yewa_rel_abund <- yewa_rel_abund_indval_sum %>%
  reduce(.x = .,
         .f = full_join)
indicators_dwf_rel_abund <- dwf_rel_abund_indval_sum %>%
  reduce(.x = .,
         .f = full_join)

indicators_yewa_rel_abund_tbl <- indicators_yewa_rel_abund %>%
  mutate(taxon_level = str_to_title(taxon_level),
         taxon_level = str_replace(string = taxon_level,
                                   pattern = "Asv_name",
                                   replacement = "ASV"),
         stat = round(stat, digits = 3),
         p_value = round(p_value, digits = 3),
         dwf = round(c(dwf*100), digits = 3),
         pcc = round(c(pcc*100), digits = 3)) %>%
  gt(groupname_col = "taxon_level") %>%
  tab_style(
    style = list(cell_borders(
      sides = c("bottom"),
      style = "solid"),
      cell_text(weight = "bold",
                size = "x-small")),
    locations = cells_row_groups()) %>%
  tab_style(style = list(cell_borders(sides = "all",
                         style = "hidden"),
                         cell_text(size = "x-small",
                                   align = "left")),
            locations = cells_body()) %>%
  tab_style(style = list(cell_text(align = "left")),
            locations = cells_title()) %>%
  tab_spanner(label = "Relative abundance (%) ",
              columns = c(dwf,
                          pcc)) %>%
  tab_style(style = cell_text(align = "left",
                              weight = "bolder",
                              size = "small",
                              whitespace = "break-spaces"),
            locations = list(cells_column_labels(),
                             cells_column_spanners())) %>%
  tab_header(title = "Statistically significant indicator species for site comparison") %>%
  cols_label(taxon = md("Taxon"),
             indicator_for = md("Indicator of"),
             stat = md("Indval stat"),
             p_value = md("p value"))
indicators_yewa_rel_abund_tbl

gtsave(data = indicators_yewa_rel_abund_tbl,
       filename = here("eggshells",
                       "microbiome",
                       "outputs",
                       "indicators_yewa_rel_abund_tbl.png"))


indicators_dwf_rel_abund_tbl <- indicators_dwf_rel_abund %>% 
  mutate(taxon_level = str_to_title(taxon_level),
         taxon_level = str_replace(string = taxon_level,
                                   pattern = "Asv_name",
                                   replacement = "ASV"),
         stat = round(stat, digits = 3),
         p_value = round(p_value, digits = 3),
         YEWA = round(c(YEWA*100), digits = 3),
         RWBL = round(c(RWBL*100), digits = 3)) %>%
  gt(groupname_col = "taxon_level") %>%
  tab_style(
    style = list(cell_borders(
      sides = c("bottom"),
      style = "solid"),
      cell_text(weight = "bold",
                size = "x-small")),
    locations = cells_row_groups()) %>%
  tab_style(style = list(cell_borders(sides = "all",
                                      style = "hidden"),
                         cell_text(size = "x-small",
                                   align = "left")),
            locations = cells_body()) %>%
  tab_style(style = list(cell_text(align = "left")),
            locations = cells_title()) %>%
  tab_spanner(label = "Relative abundance (%) ",
              columns = c(YEWA,
                          RWBL)) %>%
  tab_style(style = cell_text(align = "left",
                              weight = "bolder",
                              size = "small",
                              whitespace = "break-spaces"),
            locations = list(cells_column_labels(),
                             cells_column_spanners())) %>%
  tab_header(title = "Statistically significant indicator species for species comparison") %>%
  cols_label(taxon = md("Taxon"),
             indicator_for = md("Indicator of"),
             stat = md("Indval stat"),
             p_value = md("p value"))
indicators_dwf_rel_abund_tbl

gtsave(data = indicators_dwf_rel_abund_tbl,
       filename = here("eggshells",
                       "microbiome",
                       "outputs",
                       "indicators_dwf_rel_abund_tbl.png"))

#euler plot
#genus
genus_euler <- list(pilot_sample_asv_reps %>%
                      filter(reads > 0) %>%
                      select(asv_name, 
                             nest_id) %>%
                      distinct(),
                    pilot_taxo_reps %>%
                      select(asv_name,
                             genus),
                    pilot_swab_nests_reps %>% 
                      select(nest_id,
                             site,
                             host_species) %>%
                      unite(col = sample_group,
                            site, 
                            host_species,
                            sep = "; ") %>%
                      distinct()) %>%
  reduce(.x = .,
         .f = full_join) %>%
  select(sample_group,
         genus) %>%
  distinct() %>%
  group_by(sample_group) %>%
  group_split()

names(genus_euler) <- c("dwf; RWBL",
                        "dwf; YEWA",
                        "pcc; RWBL",
                        "pcc; YEWA")

genus_euler <- genus_euler %>%
  map(.x = .,
      .f = pull,
      genus)

genus_euler_fit <- euler(genus_euler,
                     shape = "ellipse")

jpeg(here("eggshells",
          "microbiome",
          "outputs",
          "euler_genus.jpeg"))
plot(genus_euler_fit,
     quantities = list(type = 'counts', cex = .75, fontface = 2),
     main = "Genus between sample groups")
dev.off()

#species
species_euler <- list(pilot_sample_asv_reps %>%
                      filter(reads > 0) %>%
                      select(asv_name, 
                             nest_id) %>%
                      distinct(),
                    pilot_taxo_reps %>%
                      select(asv_name,
                             species),
                    pilot_swab_nests_reps %>% 
                      select(nest_id,
                             site,
                             host_species) %>%
                      unite(col = sample_group,
                            site, 
                            host_species,
                            sep = "; ") %>%
                      distinct()) %>%
  reduce(.x = .,
         .f = full_join) %>%
  select(sample_group,
         species) %>%
  distinct() %>%
  group_by(sample_group) %>%
  group_split()

names(species_euler) <- c("dwf; RWBL",
                          "dwf; YEWA",
                          "pcc; RWBL",
                          "pcc; YEWA")

species_euler <- species_euler %>%
  map(.x = .,
      .f = pull,
      species)

species_euler_fit <- euler(species_euler,
                         shape = "ellipse")

jpeg(here("eggshells",
          "microbiome",
          "outputs",
          "euler_species.jpeg"))
plot(species_euler_fit,
     quantities = list(type = 'counts', cex = .75, fontface = 2),
     main = "Species between sample groups")
dev.off()

#family
family_euler <- list(pilot_sample_asv_reps %>%
                        filter(reads > 0) %>%
                        select(asv_name, 
                               nest_id) %>%
                        distinct(),
                      pilot_taxo_reps %>%
                        select(asv_name,
                               family),
                      pilot_swab_nests_reps %>% 
                        select(nest_id,
                               site,
                               host_species) %>%
                        unite(col = sample_group,
                              site, 
                              host_species,
                              sep = "; ") %>%
                        distinct()) %>%
  reduce(.x = .,
         .f = full_join) %>%
  select(sample_group,
         family) %>%
  distinct() %>%
  group_by(sample_group) %>%
  group_split()

names(family_euler) <- c("dwf; RWBL",
                         "dwf; YEWA",
                         "pcc:RWBL",
                         "pcc; YEWA")

family_euler <- family_euler %>%
  map(.x = .,
      .f = pull,
      family)

family_euler_fit <- euler(family_euler,
                           shape = "ellipse")

jpeg(here("eggshells",
          "microbiome",
          "outputs",
          "euler_family.jpeg"))
plot(family_euler_fit,
     quantities = list(type = 'counts', cex = .75, fontface = 2),
     main = "Family between sample groups")
dev.off()

#order
order_euler <- list(pilot_sample_asv_reps %>%
                       filter(reads > 0) %>%
                       select(asv_name, 
                              nest_id) %>%
                       distinct(),
                     pilot_taxo_reps %>% 
                       select(asv_name,
                              order),
                     pilot_swab_nests_reps %>% 
                       select(nest_id,
                              site,
                              host_species) %>%
                       unite(col = sample_group,
                             site, 
                             host_species,
                             sep = "; ") %>%
                       distinct()) %>%
  reduce(.x = .,
         .f = full_join) %>%
  select(sample_group,
         order) %>%
  distinct() %>%
  group_by(sample_group) %>%
  group_split()

names(order_euler) <- c("dwf; RWBL",
                        "dwf; YEWA",
                        "pcc; RWBL",
                        "pcc; YEWA")

order_euler <- order_euler %>%
  map(.x = .,
      .f = pull,
      order)

order_euler_fit <- euler(order_euler,
                          shape = "ellipse")

jpeg(here("eggshells",
          "microbiome",
          "outputs",
          "euler_order.jpeg"))
plot(order_euler_fit,
     quantities = list(type = 'counts', cex = .75, fontface = 2),
     main = "Order between sample groups")
dev.off()

#class
class_euler <- list(pilot_sample_asv_reps %>%
                      filter(reads > 0) %>%
                      select(asv_name, 
                             nest_id) %>%
                      distinct(),
                    pilot_taxo_reps %>% 
                      select(asv_name,
                             class),
                    pilot_swab_nests_reps %>% 
                      select(nest_id,
                             site,
                             host_species) %>%
                      unite(col = sample_group,
                            site, 
                            host_species,
                            sep = "; ") %>%
                      distinct()) %>%
  reduce(.x = .,
         .f = full_join) %>%
  select(sample_group,
         class) %>%
  distinct() %>%
  group_by(sample_group) %>%
  group_split()

names(class_euler) <- c("dwf; RWBL",
                        "dwf; YEWA",
                        "pcc; RWBL",
                        "pcc; YEWA")

class_euler <- class_euler %>%
  map(.x = .,
      .f = pull,
      class)

class_euler_fit <- euler(class_euler,
                         shape = "ellipse")

jpeg(here("eggshells",
          "microbiome",
          "outputs",
          "euler_class.jpeg"))
plot(class_euler_fit,
     quantities = list(type = 'counts', cex = .75, fontface = 2),
     main = "Class between sample groups")
dev.off()

#phylum
phylum_euler <- list(pilot_sample_asv_reps %>%
                      filter(reads > 0) %>%
                      select(asv_name, 
                             nest_id) %>%
                      distinct(),
                    pilot_taxo_reps %>% 
                      select(asv_name,
                             phylum),
                    pilot_swab_nests_reps %>% 
                      select(nest_id,
                             site,
                             host_species) %>%
                      unite(col = sample_group,
                            site, 
                            host_species,
                            sep = "; ") %>%
                      distinct()) %>%
  reduce(.x = .,
         .f = full_join) %>%
  select(sample_group,
         phylum) %>%
  distinct() %>%
  group_by(sample_group) %>%
  group_split()

names(phylum_euler) <- c("dwf; RWBL",
                        "dwf; YEWA",
                        "pcc; RWBL",
                        "pcc; YEWA")

phylum_euler <- phylum_euler %>%
  map(.x = .,
      .f = pull,
      phylum)

phylum_euler_fit <- euler(phylum_euler,
                         shape = "ellipse")

jpeg(here("eggshells",
          "microbiome",
          "outputs",
          "euler_phylum.jpeg"))
plot(phylum_euler_fit,
     quantities = list(type = 'counts', cex = .75, fontface = 2),
     main = "Phylum between sample groups")
dev.off()

#asv_name
asv_name_euler <- list(pilot_sample_asv_reps %>%
                       filter(reads > 0) %>%
                       select(asv_name, 
                              nest_id) %>%
                       distinct(),
                     pilot_taxo_reps %>% 
                       select(asv_name,
                              asv_name),
                     pilot_swab_nests_reps %>% 
                       select(nest_id,
                              site,
                              host_species) %>%
                       unite(col = sample_group,
                             site, 
                             host_species,
                             sep = "; ") %>%
                       distinct()) %>%
  reduce(.x = .,
         .f = full_join) %>%
  select(sample_group,
         asv_name) %>%
  distinct() %>%
  group_by(sample_group) %>%
  group_split()

names(asv_name_euler) <- c("dwf; RWBL",
                         "dwf; YEWA",
                         "pcc;RWBL",
                         "pcc; YEWA")

asv_name_euler <- asv_name_euler %>%
  map(.x = .,
      .f = pull,
      asv_name)

asv_name_euler_fit <- euler(asv_name_euler,
                          shape = "ellipse")

jpeg(here("eggshells",
          "microbiome",
          "outputs",
          "euler_asv.jpeg"))
plot(asv_name_euler_fit,
     quantities = list(type = 'counts', cex = .75, fontface = 2),
     main = "ASV between sample groups")
dev.off()

#### how does geographic distance affect alpha diversity? ####
# to test the effect of geographic distance on alpha diversity we will do a 
#mantel test on the bray-curtis distance matrix of alpha diversity metrics

#take bray-curtis distance of alpha diversity metrics
yewa_alpha_diversity_mean_split <- yewa_asv_sample_rared_alpha_sumd %>%
  filter(a_div_type %in% c("chao1",
                           "evenness",
                           "richness",
                           "shannon",
                           "simpson")) %>%
  select(nest_id,
         a_div_type,
         mean_a_div) %>%
  group_by(a_div_type) %>%
  group_split() %>%
  map(.x = ., 
       .f = pivot_wider,
      names_from = "a_div_type",
      values_from = "mean_a_div") %>%
  map(.x = .,
      .f = ungroup) %>%
  setNames(c("chao1",
          "evenness",
          "richness",
          "shannon",
          "simpson"))

yewa_alpha_diversity_mean_split_bray <- yewa_alpha_diversity_mean_split %>%
  map(.x = ., 
      .f = column_to_rownames,
      var = "nest_id") %>%
  map(.x = .,
      .f = vegdist) %>%
  map(.x = .,
      .f = as.matrix)

dwf_alpha_diversity_mean_split <- dwf_asv_sample_rared_alpha_sumd %>%
  filter(a_div_type %in% c("chao1",
                           "evenness",
                           "richness",
                           "shannon",
                           "simpson")) %>%
  select(nest_id,
         a_div_type,
         mean_a_div) %>%
  group_by(a_div_type) %>%
  group_split() %>%
  map(.x = ., 
      .f = pivot_wider,
      names_from = "a_div_type",
      values_from = "mean_a_div") %>%
  map(.x = .,
      .f = ungroup) %>%
  setNames(c("chao1",
             "evenness",
             "richness",
             "shannon",
             "simpson"))

dwf_alpha_diversity_mean_split_bray <- dwf_alpha_diversity_mean_split %>%
  map(.x = ., 
      .f = column_to_rownames,
      var = "nest_id") %>%
  map(.x = .,
      .f = vegdist) %>%
  map(.x = .,
      .f = as.matrix)

#check data structure to help decide what correlation metric to use
yewa_alpha_dis_plot_data <- yewa_alpha_diversity_mean_split_bray %>% 
  map(.x = .,
      .f = replace_upper_triangle,
      by = NA_character_,
      diagonal = FALSE) %>%
  map(.x = .,
      .f = as.data.frame) %>%
  map(.x = .,
      .f = rename,
      "nest_id" = "rowname") %>%
  map(.x = .,
      .f = pivot_longer,
      cols = -nest_id,
      names_to = "nest_id_compared") %>%
  reduce(.x = .,
         .f = full_join,
         by = c("nest_id",
                "nest_id_compared")) %>% 
  drop_na() %>% 
  rename(c("chao1" = "value.x",
           "evenness" = "value.y",
           "richness" = "value.x.x",
           "shannon" = "value.y.y",
           "simpson" = "value")) %>% 
  pivot_longer(cols = -c(nest_id,
                         nest_id_compared),
               values_to = "distance",
               names_to = "distance_type") %>%
  full_join(.,
             geo_distance_long %>%
              filter(nest_id %in% c(yewa_nests),
                    nest_id_compared %in% c(yewa_nests))) %>% 
  left_join(.,
            pilot_swab_nests_reps %>%
              select(nest_id,
                     host_species,
                     site) %>%
              distinct()) %>%
  left_join(.,
            pilot_swab_nests_reps %>%
              select(nest_id,
                     host_species,
                     site) %>%
              rename_with(.data = .,
                          .fn = function(x){paste0(x,
                                         "_compared")}) %>%
              distinct()) %>%
  distinct() %>%
  mutate(distance = as.numeric(distance),
         geo_distance = as.numeric(geo_distance)) %>%
  drop_na()

yewa_alpha_dis_plot_data %>%
  full_join(.,
            geo_distance_long %>%
              filter(nest_id %in% c(yewa_nests),
                     nest_id_compared %in% c(yewa_nests)) %>%
              mutate(distance_type = "geographic distance (m)",
                     geo_distance = as.numeric(geo_distance)) %>%
              rename("distance" = "geo_distance")) %>%
  ggplot(data = .) +
  geom_histogram(aes(x = as.numeric(distance)),
                 bins = 5) + #adjust bins
  facet_wrap(vars(distance_type),
             scales = "free_x") +
  labs(title = "Distribution of pair-wise bray curtis or geographic distances",
       x = "Distance")
#data is normally distributed

#check linearity
yewa_alpha_dis_plot_data %>%
  ggplot(data = .,
         aes(x = as.numeric(geo_distance),
             y = as.numeric(distance))) +
  geom_smooth() +
  geom_point() +
  facet_wrap(vars(distance_type),
             scales = "free_x")
#not a linear relationship
#so can use Kendall's rank correlation https://www.phdata.io/blog/data-science-stats-review/

dwf_alpha_dis_plot_data <- dwf_alpha_diversity_mean_split_bray %>% 
  map(.x = .,
      .f = replace_upper_triangle,
      by = NA_character_,
      diagonal = FALSE) %>%
  map(.x = .,
      .f = as.data.frame) %>%
  map(.x = .,
      .f = rename,
      "nest_id" = "rowname") %>%
  map(.x = .,
      .f = pivot_longer,
      cols = -nest_id,
      names_to = "nest_id_compared") %>%
  reduce(.x = .,
         .f = full_join,
         by = c("nest_id",
                "nest_id_compared")) %>% 
  drop_na() %>% 
  rename(c("chao1" = "value.x",
           "evenness" = "value.y",
           "richness" = "value.x.x",
           "shannon" = "value.y.y",
           "simpson" = "value")) %>% 
  pivot_longer(cols = -c(nest_id,
                         nest_id_compared),
               values_to = "distance",
               names_to = "distance_type") %>% 
  full_join(.,
            geo_distance_long %>%
              filter(nest_id %in% c(dwf_nests),
                     nest_id_compared %in% c(dwf_nests))) %>% 
  left_join(.,
            pilot_swab_nests_reps %>%
              select(nest_id,
                     host_species,
                     site) %>%
              distinct()) %>%
  left_join(.,
            pilot_swab_nests_reps %>%
              select(nest_id,
                     host_species,
                     site) %>%
              rename_with(.data = .,
                          .fn = function(x){paste0(x,
                                                   "_compared")}) %>%
              distinct()) %>%
  distinct() %>%
  mutate(distance = as.numeric(distance),
         geo_distance = as.numeric(geo_distance)) %>%
  drop_na()

dwf_alpha_dis_plot_data %>%
  full_join(.,
            geo_distance_long %>%
              filter(nest_id %in% c(dwf_nests),
                     nest_id_compared %in% c(dwf_nests)) %>%
              mutate(distance_type = "geographic distance (m)",
                     geo_distance = as.numeric(geo_distance)) %>%
              rename("distance" = "geo_distance")) %>%
  ggplot(data = .) +
  geom_histogram(aes(x = as.numeric(distance)),
                 bins = 10) + #adjust bins
  facet_wrap(vars(distance_type),
             scales = "free_x") +
  labs(title = "Distribution of pair-wise bray curtis or geographic distances",
       x = "Distance")
#not normally distributed
#check linearity
dwf_alpha_dis_plot_data %>%
  ggplot(data = .,
         aes(x = as.numeric(geo_distance),
             y = as.numeric(distance))) +
  geom_smooth() +
  geom_point() +
  facet_wrap(vars(distance_type),
             scales = "free_x")
#not a linear relationship
#just yewa
yewa_distance_mat <- distance_mat[yewa_nests, yewa_nests]
j <- names(yewa_alpha_diversity_mean_split_bray)
yewa_alpha_mantel <- list()
for (i in j) {
  yewa_alpha_mantel[[i]] <- mantel(xdis = yewa_alpha_diversity_mean_split_bray[[i]],
                              ydis = yewa_distance_mat,
                              method = "kendall")
}

#make table with results
j <- names(yewa_alpha_mantel)
yewa_alpha_mantel_rez <- list()
for (i in j) {
  x <- tibble(distance = i,
              mantel_statistic.site = yewa_alpha_mantel[[i]]$statistic,
              significance.site = yewa_alpha_mantel[[i]]$signif)
  yewa_alpha_mantel_rez[[i]] <- x
}
yewa_alpha_mantel_rez <- yewa_alpha_mantel_rez %>%
  reduce(.x = .,
         .f = full_join)

#just dwf
dwf_distance_mat <- distance_mat[dwf_nests, dwf_nests]
j <- names(dwf_alpha_diversity_mean_split_bray)
dwf_alpha_mantel <- list()
for (i in j) {
  dwf_alpha_mantel[[i]] <- mantel(xdis = dwf_alpha_diversity_mean_split_bray[[i]],
                                   ydis = dwf_distance_mat,
                                   method = "kendall")
}

#make table with results
j <- names(dwf_alpha_mantel)
dwf_alpha_mantel_rez <- list()
for (i in j) {
  x <- tibble(distance = i,
              mantel_statistic.species = dwf_alpha_mantel[[i]]$statistic,
              significance.species = dwf_alpha_mantel[[i]]$signif)
  dwf_alpha_mantel_rez[[i]] <- x
}
dwf_alpha_mantel_rez <- dwf_alpha_mantel_rez %>%
  reduce(.x = .,
         .f = full_join)

alpha_mantel_rez <- full_join(yewa_alpha_mantel_rez,
                              dwf_alpha_mantel_rez)

alpha_mantel_rez_tbl <- alpha_mantel_rez %>% 
  mutate(across(.cols = where(is.numeric),
                .fns = ~ round(x = .x,
                               digits = 3)),
         distance = case_when(distance == "chao1" ~ "Chao1 Index",
                               distance == "evenness" ~ "Species evenness",
                               distance == "richness" ~ "Species richness",
                               distance == "shannon" ~ "Shannon's Diversity Index", 
                               distance == "simpson" ~ "Simpson's Diversity Index")) %>%
  gt() %>%
  tab_style(style = list(cell_borders(sides = "all",
                                      style = "hidden"),
                         cell_text(size = "x-small",
                                   align = "left")),
            locations = cells_body()) %>%
  tab_spanner(label = "Species",
              columns = c(ends_with(".species"))) %>% 
  tab_spanner(label = "Site",
              columns = c(ends_with(".site"))) %>%
  tab_style(style = list(cell_text(align = "left")),
            locations = cells_title()) %>%
  tab_style(style = cell_text(align = "left",
                              weight = "bolder",
                              size = "small",
                              whitespace = "break-spaces"),
            locations = list(cells_column_labels())) %>%
  tab_header(title = "Mantel test of alpha diversity and geographic distance") %>%
  cols_label(distance = "Alpha diversity metric",
             mantel_statistic.site = md("Mantel statistic"),
             significance.site = md("*p* value"),
             mantel_statistic.species = md("Mantel statistic"),
             significance.species = md("*p* value"))
alpha_mantel_rez_tbl
gtsave(data = alpha_mantel_rez_tbl,
       filename = here("eggshells",
                       "microbiome",
                       "outputs",
                       "alpha_mantel_rez_tbl.png"))

#plot results
full_join(dwf_alpha_dis_plot_data %>%
            mutate(comparison = "Species"),
          yewa_alpha_dis_plot_data %>%
            mutate(comparison = "Site")) %>%
  mutate(sample_group = case_when(comparison == "Species" & 
                                    host_species == "YEWA" & 
                                    host_species_compared == "YEWA" ~ "YEWA",
                                  comparison == "Species" & 
                                    host_species == "RWBL" & 
                                    host_species_compared == "RWBL" ~ "RWBL",
                                  comparison == "Species" & 
                                    host_species != host_species_compared ~ "different",
                                  comparison == "Site" & 
                                    site == "dwf" & 
                                    site_compared == "dwf" ~ "dwf",
                                  comparison == "Site" & 
                                    site == "pcc" & 
                                    site_compared == "pcc" ~ "pcc",
                                  comparison == "Site" & 
                                    site != site_compared ~ "different"),
         distance_type = case_when(distance_type == "chao1" ~ "Chao1 Index",
                              distance_type == "evenness" ~ "Species evenness",
                              distance_type == "richness" ~ "Species richness",
                              distance_type == "shannon" ~ "Shannon's\nDiversity Index", 
                              distance_type == "simpson" ~ "Simpson's\nDiversity Index")) %>% 
  ggplot(data = .,
         aes(x = geo_distance,
             y = distance,
             colour = sample_group,
             fill = sample_group
         )) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "glm",
              alpha = 0.2) +
  facet_grid(facets = comparison ~ distance_type,
             scales = "free_x") + #for some reason is not working
  labs(y = "Bray Curtis distance",
       x = "Geographic distance (m)",
       colour = "Sample group",
       fill = "Sample group") +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 5),
        strip.text = element_text(size = 7))

ggsave(here("eggshells", 
            "microbiome",
            "outputs",
            "alpha_mantel_plot.jpeg"))

#### how does geographic distance affect beta diversity? ####
# to test the effect of geographic distance on beta diversity we will do a 
#mantel test on the bray-curtis distance matrix of alpha diversity metrics

beta_dis_plot_data <- list(yewa_rared_bray %>%
       as.matrix() %>%
       replace_upper_triangle(x = .,
                              by = NA_character_,
                              diagonal = FALSE) %>%
       as.data.frame() %>%
       rename("nest_id" = "rowname") %>%
       pivot_longer(cols = - nest_id,
                    names_to = "nest_id_compared",
                    values_to = "bray") %>%
       drop_na() %>%
         mutate(comparison = "Site"),
       yewa_rared_jacc %>%
       as.matrix() %>%
       replace_upper_triangle(x = .,
                              by = NA_character_,
                              diagonal = FALSE) %>%
       as.data.frame() %>%
       rename("nest_id" = "rowname") %>%
       pivot_longer(cols = - nest_id,
                    names_to = "nest_id_compared",
                    values_to = "jacc") %>%
       drop_na() %>%
         mutate(comparison = "Site"),
       dwf_rared_bray %>%
         as.matrix() %>%
         replace_upper_triangle(x = .,
                                by = NA_character_,
                                diagonal = FALSE) %>%
         as.data.frame() %>%
         rename("nest_id" = "rowname") %>%
         pivot_longer(cols = - nest_id,
                      names_to = "nest_id_compared",
                      values_to = "bray") %>%
         drop_na() %>%
         mutate(comparison = "Species"),
       dwf_rared_jacc %>%
         as.matrix() %>%
         replace_upper_triangle(x = .,
                                by = NA_character_,
                                diagonal = FALSE) %>%
         as.data.frame() %>%
         rename("nest_id" = "rowname") %>%
         pivot_longer(cols = - nest_id,
                      names_to = "nest_id_compared",
                      values_to = "jacc") %>%
         drop_na() %>%
         mutate(comparison = "Species"),
     geo_distance_long) %>%
  reduce(.x = .,
         .f = full_join) %>%
  pivot_longer(cols = -c(nest_id,
                         nest_id_compared,
                         geo_distance,
                         comparison),
               names_to = "distance_type",
               values_to = "distance") %>%
  mutate(distance = as.numeric(distance),
         geo_distance = as.numeric(geo_distance)) %>% 
  left_join(.,
            pilot_swab_nests_reps %>%
              select(nest_id,
                     host_species,
                     site) %>%
              distinct()) %>%
  left_join(.,
            pilot_swab_nests_reps %>%
              select(nest_id,
                     host_species,
                     site) %>%
              rename_with(.data = .,
                          .fn = function(x){paste0(x,
                                                   "_compared")}) %>%
              distinct()) %>%
  distinct() %>%
  mutate(sample_group = case_when(comparison == "Species" & 
                                    host_species == "YEWA" & 
                                    host_species_compared == "YEWA" ~ "YEWA",
                                  comparison == "Species" & 
                                    host_species == "RWBL" & 
                                    host_species_compared == "RWBL" ~ "RWBL",
                                  comparison == "Species" & 
                                    host_species != host_species_compared ~ "different",
                                  comparison == "Site" & 
                                    site == "dwf" & 
                                    site_compared == "dwf" ~ "dwf",
                                  comparison == "Site" & 
                                    site == "pcc" & 
                                    site_compared == "pcc" ~ "pcc",
                                  comparison == "Site" & 
                                    site != site_compared ~ "different")) %>%
  drop_na()

beta_dis_plot_data %>%
  full_join(.,
            geo_distance_long %>%
              mutate(distance_type = "geographic distance (m)",
                     geo_distance = as.numeric(geo_distance)) %>%
              rename("distance" = "geo_distance")) %>%
  ggplot(data = .) +
  geom_histogram(aes(x = as.numeric(distance))) +
  facet_wrap(vars(distance_type,
                  comparison),
             scales = "free_x") +
  labs(title = "Distribution of pair-wise biological distance or geographic distances",
       x = "Distance")
#not normally distributed
#all
jacc_mantel <- mantel(xdis = asv_sample_rared_jacc,
                      ydis = distance_mat,
                      method = "kendall")

bray_mantel <- mantel(xdis = asv_sample_rared_bray,
                      ydis = distance_mat,
                      method = "kendall")
#yewa
yewa_jacc_mantel <- mantel(xdis = yewa_rared_jacc,
                      ydis = yewa_distance_mat,
                      method = "kendall")

yewa_bray_mantel <- mantel(xdis = yewa_rared_bray,
                      ydis = yewa_distance_mat,
                      method = "kendall")
#dwf
dwf_jacc_mantel <- mantel(xdis = dwf_rared_jacc,
                           ydis = dwf_distance_mat,
                           method = "kendall")

dwf_bray_mantel <- mantel(xdis = dwf_rared_bray,
                           ydis = dwf_distance_mat,
                           method = "kendall")

#make table with results

beta_mantel_rez <- tibble(distance = c("jaccards",
                                       "bray curtis"),
                          mantel_statistic.site = c(yewa_jacc_mantel$statistic,
                                                    yewa_bray_mantel$statistic),
                          significance.site = c(yewa_jacc_mantel$signif,
                                                yewa_bray_mantel$signif),
                          mantel_statistic.species = c(dwf_jacc_mantel$statistic,
                                                       dwf_bray_mantel$statistic),
                          significance.species = c(dwf_jacc_mantel$signif,
                                                   dwf_bray_mantel$signif))
beta_mantel_rez_tbl <- beta_mantel_rez %>% 
  mutate(across(.cols = where(is.numeric),
                .fns = ~ round(x = .x,
                               digits = 3)),
         distance = case_when(distance == "jaccards" ~ "Jaccard Index",
                              distance == "bray curtis" ~ "Bray Curtis dissimilarity")) %>%
  gt() %>%
  tab_style(style = list(cell_borders(sides = "all",
                                      style = "hidden"),
                         cell_text(size = "x-small",
                                   align = "left")),
            locations = cells_body()) %>%
  tab_spanner(label = "Species",
              columns = c(ends_with(".species"))) %>% 
  tab_spanner(label = "Site",
              columns = c(ends_with(".site"))) %>%
  tab_style(style = list(cell_text(align = "left")),
            locations = cells_title()) %>%
  tab_style(style = cell_text(align = "left",
                              weight = "bolder",
                              size = "small",
                              whitespace = "break-spaces"),
            locations = list(cells_column_labels())) %>%
  tab_header(title = "Mantel test of beta diversity and geographic distance") %>%
  cols_label(distance = "Beta diversity metric",
             mantel_statistic.site = md("Mantel statistic"),
             significance.site = md("*p* value"),
             mantel_statistic.species = md("Mantel statistic"),
             significance.species = md("*p* value"))
beta_mantel_rez_tbl
gtsave(data = beta_mantel_rez_tbl,
       filename = here("eggshells",
                       "microbiome",
                       "outputs",
                       "beta_mantel_rez_tbl.png"))
#plot results
beta_dis_plot_data %>%
  mutate(distance_type = case_when(distance_type == "jacc" ~ "Jaccard Index",
                              distance_type == "bray" ~ "Bray Curtis dissimilarity")) %>%
  ggplot(data = .,
         aes(x = geo_distance,
             y = distance,
             colour = sample_group,
             fill = sample_group
         )) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "glm",
              alpha = 0.2) +
  facet_grid(facets = comparison ~ distance_type,
             scales = "free") + # not working for some reason
  labs(y = "Beta diversity distance",
       x = "Geographic distance (m)",
       colour = "Sample group",
       fill = "Sample group") +
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 7))

ggsave(here("eggshells", 
            "microbiome",
            "outputs",
            "beta_mantel_plot.jpeg"))

#### how does site damage, nesting niche (species) and clutch initiation date effect eggshell microbiome alpha diversity? ####

alpha_plot <- yewa_alpha_diversity_mean_split %>%
  reduce(.x = .,
         .f = full_join) %>%
  mutate(comparison = "Site") %>%
  full_join(.,
          dwf_alpha_diversity_mean_split %>%
            reduce(.x = .,
                   .f = full_join) %>%
            mutate(comparison = "Species")) %>%
  pivot_longer(cols = c("chao1",
                        "evenness",
                        "richness",
                        "simpson",
                        "shannon"),
               names_to = "alpha_type",
               values_to = "alpha_val") %>%
  left_join(.,
            pilot_swab_nests_reps %>%
              select(nest_id,
                   host_species,
                   site) %>%
              distinct()) %>%
  mutate(x_axis = case_when(comparison == "Species" ~ host_species,
                            comparison == "Site" ~ site))
alpha_plot <- alpha_plot %>%
  mutate(alpha_type = case_when(alpha_type == "chao1" ~ "Chao1",
                                alpha_type == "evenness" ~ "Evenness",
                                alpha_type == "richness" ~ "Richness",
                                alpha_type == "shannon" ~ "Shannon", 
                                alpha_type == "simpson" ~ "Simpson"))

alpha_plot %>%
  unite(col = sample_group,
        host_species,
        site) %>%
  ggplot() +
  geom_boxplot(aes(y = sample_group,
                   x = alpha_val)) +
  geom_point(aes(y = sample_group,
                 x = alpha_val),
             alpha = 0.5) +
  facet_wrap(facets = "alpha_type",
             scales = "free_x",
             ncol = 5) +
  theme_classic() +
  labs(y = element_blank(),
       x = "Alpha diversity value") +
  theme(text = element_text(size = 15),
        strip.text = element_text(size = 15),
        axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.ticks.y = element_blank()) +
  scale_y_discrete(position = "right")

ggsave(here("eggshells", 
            "microbiome",
            "outputs",
            "alpha_div_plot.jpeg"))

##check normality##
#http://www.sthda.com/english/wiki/normality-test-in-r
alpha_diversity_mean <- alpha_plot %>%
  unite(col = alpha_type_host,
        alpha_type,
        host_species,
        remove = FALSE) %>%
  unite(col = alpha_type_site,
        alpha_type,
        site,
        remove = FALSE)
  
j <- alpha_diversity_mean %>%
  filter(site == "dwf") %>%
  pull(alpha_type_host) %>%
  unique()
alpha_normal_dwf <- list() 
for (i in j) {
  
  alpha_normal_dwf[[i]] <- alpha_diversity_mean %>%
    filter(site == "dwf",
           alpha_type_host == i) %>%
    pull(alpha_val) %>%
    shapiro_test()

  }
#if value of pvalue is >0.05 then the distribution is not significantly different than the normal distribution
#note that small sample sizes often pass
alpha_normal_dwf
#shannon_RWBL not normal, thius is because we are using the dummy data, we will proceed as though it is normally distributed

j <- alpha_diversity_mean %>%
  filter(host_species == "YEWA") %>%
  pull(alpha_type_site) %>%
  unique()
alpha_normal_yewa <- list() 
for (i in j) {
  alpha_normal_yewa[[i]] <- alpha_diversity_mean %>%
    filter(host_species == "YEWA",
           alpha_type_site == i) %>%
    pull(alpha_val) %>%
    shapiro_test()
}
alpha_normal_yewa
# all is normal

##check for same variance##
j <- alpha_diversity_mean %>%
  pull(alpha_type) %>%
  unique()
alpha_var_dwf <- list() 
for (i in j) {
  
  alpha_var_dwf[[i]] <- alpha_diversity_mean %>%
    filter(site == "dwf",
           alpha_type == i) %>%
    var.test(alpha_val ~ host_species,
             data = .)
  
}
#if value of pvalue is >0.05 then the variance of the two groups is not significantly different
alpha_var_dwf # all good

j <- alpha_diversity_mean %>%
  pull(alpha_type) %>%
  unique()
alpha_var_yewa <- list() 
for (i in j) {
  
  alpha_var_yewa[[i]] <- alpha_diversity_mean %>%
    filter(host_species == "YEWA",
           alpha_type == i) %>%
    var.test(alpha_val ~ site,
             data = .)

  }
alpha_var_yewa #all good

##compute t test on all
#http://www.sthda.com/english/wiki/unpaired-two-samples-t-test-in-r

j <- alpha_diversity_mean %>%
  pull(alpha_type) %>%
  unique()
alpha_ttest_yewa <- list() 
for (i in j) {
  
  alpha_ttest_yewa[[i]] <- alpha_diversity_mean %>%
    filter(host_species == "YEWA",
           alpha_type == i) %>%
    t.test(alpha_val ~ site,
           data = .,
           var.equal = TRUE)
  
}
alpha_ttest_yewa # no significant differences in alpha diversity between sites

j <- names(alpha_ttest_yewa)
ttest_rez <- list()
for (i in j) {
  ttest_rez[[paste(i, "site")]] <- tibble(t_statistic = alpha_ttest_yewa[[i]]$statistic %>% as.vector(),
                           degrees_of_freedom = alpha_ttest_yewa[[i]]$parameter %>% as.vector(),
                           p_value = alpha_ttest_yewa[[i]]$p.value,
                           alpha_div = i,
                           comparison = "site")
}

j <- alpha_diversity_mean %>%
  pull(alpha_type) %>%
  unique()
alpha_ttest_dwf <- list() 
for (i in j) {
  
  alpha_ttest_dwf[[i]] <- alpha_diversity_mean %>%
    filter(site == "dwf",
           alpha_type == i) %>%
    t.test(alpha_val ~ host_species,
           data = .,
           var.equal = TRUE)
  
}
alpha_ttest_dwf # no significant differences in shannon, simpson, eveness between species

j <- names(alpha_ttest_dwf)
for (i in j) {
  ttest_rez[[paste(i, "species")]] <- tibble(t_statistic = alpha_ttest_dwf[[i]]$statistic %>% as.vector(),
                           degrees_of_freedom = alpha_ttest_dwf[[i]]$parameter %>% as.vector(),
                           p_value = alpha_ttest_dwf[[i]]$p.value,
                           alpha_div = i,
                           comparison = "species")
}
ttest_rez <- ttest_rez %>%
  reduce(.x =.,
         .f = full_join)

ttest_rez_table <- ttest_rez %>% 
  mutate(t_statistic = round(t_statistic, digits = 3),
       p_value = round(p_value, digits = 3),
       alpha_div = case_when(alpha_div == "chao1" ~ "Chao1 Index",
                             alpha_div == "evenness" ~ "Species evenness",
                             alpha_div == "richness" ~ "Species richness",
                             alpha_div == "shannon" ~ "Shannon's Diversity Index", 
                             alpha_div == "simpson" ~ "Simpson's Diversity Index")) %>%
  gt(groupname_col = "alpha_div") %>%
  tab_style(
    style = list(cell_borders(
      sides = c("bottom"),
      style = "solid"),
      cell_text(weight = "bold",
                size = "x-small")),
    locations = cells_row_groups()) %>%
  tab_style(style = list(cell_borders(sides = "all",
                                      style = "hidden"),
                         cell_text(size = "x-small",
                                   align = "left")),
            locations = cells_body()) %>%
  tab_style(style = list(cell_text(align = "left")),
            locations = cells_title()) %>%
  tab_style(style = cell_text(align = "left",
                              weight = "bolder",
                              size = "small",
                              whitespace = "break-spaces"),
            locations = list(cells_column_labels())) %>%
  tab_header(title = "Two sample t-test for alpha diversity differences between species and sites") %>%
  cols_label(t_statistic = md("*t* statistic"),
             degrees_of_freedom = md("Degrees of freedom"),
             p_value = md("*p* value"),
             comparison = md("Comparison"))
ttest_rez_table

gtsave(data = ttest_rez_table,
       filename = here("eggshells",
                       "microbiome",
                       "outputs",
                       "ttest_rez_table.png"))


## look for covariation with day of year
#prep data 
meta_dat <- left_join(pilot_metadata,
                      pilot_swab_nests %>%
                        select(nest_id, day_of_year, date))
alpha_diversity_mean <- left_join(alpha_diversity_mean,
                                  meta_dat)
#asses assumptions
#linearity
alpha_diversity_mean %>%
  filter(site == "dwf") %>%
  ggplot(data = .,
         aes(x = day_of_year,
             y = alpha_val,
             colour = host_species)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(facets = "alpha_type",
             scales = "free_y") +
  stat_cor()
#all strong relationships for yewa

alpha_diversity_mean %>%
  filter(host_species == "YEWA") %>%
  ggplot(data = .,
         aes(x = day_of_year,
             y = alpha_val,
             colour = site)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(facets = "alpha_type",
             scales = "free_y") +
  stat_cor()
#a strong linear relationship for pcc for all measures, not so for dwf
#but note there is only 3 points for pcc and two of them are very close to each other

#homogeneity fo regression slopes
#no significant interaction between covariate and grouping variable
j <- alpha_diversity_mean %>%
  pull(alpha_type) %>%
  unique()
spe_interact <- list()
for (i in j) {
  spe_interact[[i]] <- alpha_diversity_mean %>%
    filter(site == "dwf",
           alpha_type == i) %>%
    anova_test(alpha_val ~ host_species*day_of_year)  
}
spe_interact
#significant interaction for simpson
j <- alpha_diversity_mean %>%
  pull(alpha_type) %>%
  unique()
site_interact <- list()
for (i in j) {
  site_interact[[i]] <- alpha_diversity_mean %>%
    filter(host_species == "YEWA",
           alpha_type == i) %>%
    anova_test(alpha_val ~ site*day_of_year)  
}
site_interact
#significant interaction for all

#normality of residuals
j <- alpha_diversity_mean %>%
  pull(alpha_type) %>%
  unique()
spe_lm <- list()
spe_rez <- list()
for (i in j) {
  spe_lm[[i]] <- alpha_diversity_mean %>%
    filter(site == "dwf",
           alpha_type == i) %>%
    lm(alpha_val ~ host_species+day_of_year, 
               data = .) %>%
    augment()
  
  spe_rez[[i]] <- spe_lm[[i]] %>%
    pull(.resid) %>%
    shapiro_test()
}
spe_rez
#eveness residuals are significantly different from the normal distribution, but for demonstration purposes will assume they are not

j <- alpha_diversity_mean %>%
  pull(alpha_type) %>%
  unique()
site_lm <- list()
site_rez <- list()
for (i in j) {
  site_lm[[i]] <- alpha_diversity_mean %>%
    filter(host_species == "YEWA",
           alpha_type == i) %>%
    lm(alpha_val ~ site + day_of_year,
       data = .) %>%
    augment()
  
  site_rez[[i]] <- site_lm[[i]] %>%
    pull(.resid) %>%
    shapiro_test() 
}
site_rez
#none of residuals are significantly different from the normal distribution

#homogeneity of variance
#assumption of equal variance of residuals for each group
j <- alpha_diversity_mean %>%
  pull(alpha_type) %>%
  unique()
spe_homg <- list()
for (i in j) {
  spe_homg[[i]] <- spe_lm[[i]] %>%
    levene_test(.resid ~ host_species)
}
spe_homg 
#non significant so can assume homogeneity of variance
j <- alpha_diversity_mean %>%
  pull(alpha_type) %>%
  unique()
site_homg <- list()
for (i in j) {
  site_homg[[i]] <- site_lm[[i]] %>%
    levene_test(.resid ~ site)
}
site_homg
#non significant so can assume homogeneity of variance

#outliers
j <- alpha_diversity_mean %>%
  pull(alpha_type) %>%
  unique()
spe_out <- list()
for (i in j) {
  spe_out[[i]] <- spe_lm[[i]] %>%
    filter(abs(.std.resid) > 3)
}
spe_out 
#no outliers
j <- alpha_diversity_mean %>%
  pull(alpha_type) %>%
  unique()
site_out <- list()
for (i in j) {
  site_out[[i]] <- site_lm[[i]] %>%
    filter(abs(.std.resid) > 3)
}
site_out
#no outliers

#ancova computation
j <- alpha_diversity_mean %>%
  pull(alpha_type) %>%
  unique()
ancova_species <- list() 
for (i in j) {
  
  ancova_species[[i]] <- alpha_diversity_mean %>%
    filter(site == "dwf",
           alpha_type == i) %>%
    anova_test(alpha_val ~ day_of_year * host_species, #covariate goes first and interaction required - if interaction is significant can't use ANCOVA https://r.qcbs.ca/workshop04/book-en/analysis-of-covariance-ancova.html
       data = .) 
}
ancova_species
#DFn Degrees of Freedom in the numerator (i.e. DF effect)
#DFd Degrees of Freedom in the denominator (i.e., DF error).
#F F-value.
#p p-value (probability of the data given the null hypothesis).
#ges Generalized Eta-Squared measure of effect size.

j <- names(ancova_species)
ancova_rez_spe <- list()
for (i in j) {
  ancova_rez_spe[[paste(i, "species")]]<- tibble(effect = ancova_species[[i]]$Effect,
                                          DFn.spe = ancova_species[[i]]$DFn,
                                          DFd.spe = ancova_species[[i]]$DFd,
                                          F_stat.spe = ancova_species[[i]]$F,
                                          p_val.spe = ancova_species[[i]]$p,
                                          effect_size.spe = ancova_species[[i]]$ges,
                                          alpha_div = i)
}

j <- alpha_diversity_mean %>%
  pull(alpha_type) %>%
  unique()
ancova_site <- list() 
for (i in j) {
  
  ancova_site[[i]] <- alpha_diversity_mean %>%
    filter(host_species == "YEWA",
           alpha_type == i) %>%
    anova_test(alpha_val ~ day_of_year * site,
       data = .)
}
ancova_site
#all significant interactions

j <- names(ancova_site)
ancova_rez_site <- list()
for (i in j) {
  ancova_rez_site[[paste(i, "site")]]<- tibble(effect = ancova_site[[1]]$Effect,
         DFn.site = ancova_site[[1]]$DFn,
         DFd.site = ancova_site[[1]]$DFd,
         F_stat.site = ancova_site[[1]]$F,
         p_val.site = ancova_site[[1]]$p,
         effect_size.site = ancova_site[[1]]$ges,
         alpha_div = i)
}
ancova_rez <- full_join(ancova_rez_spe %>%
                          reduce(.x = .,
                                 .f = full_join),
                        ancova_rez_site %>%
                          reduce(.x = .,
                                 .f = full_join)) 
#plot
alpha_plot %>%
  filter(comparison == "Site") %>%
  left_join(.,
            pilot_swab_nests_reps %>%
              select(nest_id,
                     host_clutch_initiation_date) %>%
              distinct()) %>%
  ggplot(data = .,
         aes(y = host_clutch_initiation_date,
             x = alpha_val,
             colour = x_axis)) +
  geom_line(data = alpha_plot %>%
              filter(comparison == "Site") %>%
              left_join(.,
                        pilot_swab_nests_reps %>%
                          select(nest_id,
                                 host_clutch_initiation_date) %>%
                          distinct()) %>%
              group_by(x_axis,
                       alpha_type,
                       host_clutch_initiation_date) %>%
              summarise(alpha_val = median(alpha_val)) %>%
              ungroup(),
            linewidth = 1,
            show.legend = FALSE,
            alpha = 0.7) +
  geom_jitter(size = 0.5,
              width = 0.1,
              height = 0) +
  facet_wrap("alpha_type",
             scales = "free_x",
             ncol = 5) +
  labs(y = "Clutch initaion date",
       x = "Alpha diversity value",
       colour = "Sample group") +
  theme_classic() +
  theme(text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  scale_colour_manual(values = c("#F8766D","#7CAE00")) 

ggsave(here("eggshells", 
            "microbiome",
            "outputs",
            "alpha_div_by_date_site_plot.jpeg"))

alpha_plot %>%
  filter(comparison == "Species") %>%
  left_join(.,
            pilot_swab_nests_reps %>%
              select(nest_id,
                     host_clutch_initiation_date) %>%
              distinct()) %>%
  ggplot(data = .,
         aes(y = host_clutch_initiation_date,
             x = alpha_val,
             colour = x_axis)) +
  geom_line(data = alpha_plot %>% 
              filter(comparison == "Species") %>%
              left_join(.,
                        pilot_swab_nests_reps %>%
                          select(nest_id,
                                 host_clutch_initiation_date) %>%
                          distinct()) %>% 
              group_by(x_axis,
                       alpha_type,
                       host_clutch_initiation_date) %>%
              summarise(alpha_val = median(alpha_val)) %>%
              ungroup(),
            linewidth = 1,
            show.legend = FALSE,
            alpha = 0.7) +
  geom_jitter(size = 0.5,
              width = 0.1,
              height = 0) +
  facet_wrap("alpha_type",
             scales = "free_x",
             ncol = 5) +
  labs(y = "Clutch initaion date",
       x = "Alpha diversity value",
       colour = "Sample group") +
  theme_classic() +
  theme(text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  scale_colour_manual(values = c("#00A9FF","#F8766D"))

ggsave(here("eggshells", 
            "microbiome",
            "outputs",
            "alpha_div_by_date_spe_plot.jpeg"))
  
#table
ancova_rez_tbl <- ancova_rez %>% 
  mutate(effect = str_replace(string = effect,
                              pattern = "day_of_year",
                              replacement = "clutch initiation date")) %>%
  gt(groupname_col = "alpha_div") %>%
  tab_style(
    style = list(cell_borders(
      sides = c("bottom"),
      style = "solid"),
      cell_text(weight = "bold",
                size = "x-small")),
    locations = cells_row_groups()) %>%
  tab_style(style = list(cell_borders(sides = "all",
                                      style = "hidden"),
                         cell_text(size = "x-small",
                                   align = "left")),
            locations = cells_body()) %>%
  tab_style(style = list(cell_text(align = "left")),
            locations = cells_title()) %>%
  tab_spanner(label = "Sepecies",
              columns = c(ends_with(".spe"))) %>% 
  tab_spanner(label = "Site",
              columns = c(ends_with(".site"))) %>%
  tab_style(style = cell_text(align = "left",
                              weight = "bolder",
                              size = "small",
                              whitespace = "break-spaces"),
            locations = list(cells_column_labels(),
                             cells_column_spanners())) %>%
  tab_header(title = "ANCOVA results for alpha diversity differences between species and sites") %>%
  cols_label(effect= md("Effect"),
             DFn.spe = md("Numerator Degrees of Freedom"),
             DFd.spe = md("Denomenator Degrees of Freedom"),
             F_stat.spe = md("F-value"),
             p_val.spe = md("*p*-value"),
             effect_size.spe = md("Effect size"),
             DFn.site = md("Numerator Degrees of Freedom"),
             DFd.site = md("Denomenator Degrees of Freedom"),
             F_stat.site = md("F-value"),
             p_val.site = md("*p*-value"),
             effect_size.site = md("Effect size"))
ancova_rez_tbl

gtsave(data = ancova_rez_tbl,
       filename = here("eggshells", 
                       "microbiome",
                       "outputs",
                       "ancova_rez_tbl.png"))


#### how does site damage, nesting niche (species), and clutch initiation data affect eggshell microbiome beta diversity ####
# use complete linkage clustering to look for discontinuities in the data
comclust_yewa_rared_jacc <- hclust(d = yewa_rared_jacc,
       method = "complete")
plot(comclust_yewa_rared_jacc)

comclust_yewa_rared_bray <- hclust(d = yewa_rared_bray,
                                   method = "complete")
plot(comclust_yewa_rared_bray)

comclust_dwf_rared_jacc <- hclust(d = dwf_rared_jacc,
                                   method = "complete")
plot(comclust_dwf_rared_jacc)

comclust_dwf_rared_bray <- hclust(d = dwf_rared_bray,
                                   method = "complete")
plot(comclust_dwf_rared_bray)
#can make dendrograms nicer http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning#ggdendro-package-ggplot2-and-dendrogram

# Build dendrogram object from hclust results
yewa_rared_jacc_dend <- as.dendrogram(comclust_yewa_rared_jacc)
yewa_rared_bray_dend <- as.dendrogram(comclust_yewa_rared_bray)
dwf_rared_jacc_dend <- as.dendrogram(comclust_dwf_rared_jacc)
dwf_rared_bray_dend <- as.dendrogram(comclust_dwf_rared_bray)
# Extract the data (for rectangular lines)
# Type can be "rectangle" or "triangle"
yewa_rared_jacc_dend_data <- dendro_data(yewa_rared_jacc_dend, type = "rectangle")
yewa_rared_bray_dend_data <- dendro_data(yewa_rared_bray_dend, type = "rectangle")
dwf_rared_jacc_dend_data <- dendro_data(dwf_rared_jacc_dend, type = "rectangle")
dwf_rared_bray_dend_data <- dendro_data(dwf_rared_bray_dend, type = "rectangle")
#prep data for plotting
clust_labels <- list(yewa_rared_jacc_dend_data[["labels"]] %>%
                       mutate(group = "Site Jaccards"),
                     yewa_rared_bray_dend_data[["labels"]] %>%
                       mutate(group = "Site Bray Curtis"),
                     dwf_rared_jacc_dend_data[["labels"]] %>%
                       mutate(group = "Species Jaccards"),
                     dwf_rared_bray_dend_data[["labels"]] %>%
                       mutate(group = "Species Bray Curtis")) %>%
  reduce(.x = .,
         .f = full_join) %>%
  full_join(x = .,
            y = pilot_swab_nests_reps %>%
              select(nest_id,
                     site,
                     host_species) %>%
              distinct(),
            by = c("label" = "nest_id")) %>%
  mutate(group_legend = case_when(str_detect(string = group,
                                             pattern = "Species") == TRUE ~ host_species,
                                  str_detect(string = group,
                                             pattern = "Site") == TRUE ~ site)) %>%
  mutate(comparison = case_when(str_detect(string = group,
                                           pattern = "Site Jaccards") == TRUE ~ "Site",
                                str_detect(string = group,
                                           pattern = "Site Bray Curtis") == TRUE ~ "Site" ,
                                str_detect(string = group,
                                           pattern = "Species Jaccards") == TRUE ~ "Species",
                                str_detect(string = group,
                                           pattern = "Species Bray Curtis") == TRUE ~ "Species"),
         dist_metric = case_when(str_detect(string = group,
                                            pattern = "Site Jaccards") == TRUE ~ "Jaccards",
                                 str_detect(string = group,
                                            pattern = "Site Bray Curtis") == TRUE ~ "Bray Curtis" ,
                                 str_detect(string = group,
                                            pattern = "Species Jaccards") == TRUE ~ "Jaccards",
                                 str_detect(string = group,
                                            pattern = "Species Bray Curtis") == TRUE ~ "Bray Curtis"))

clust_segments <- list(yewa_rared_jacc_dend_data[["segments"]] %>%
                         mutate(group = "Site Jaccards"),
                       yewa_rared_bray_dend_data[["segments"]] %>%
                         mutate(group = "Site Bray Curtis"),
                       dwf_rared_jacc_dend_data[["segments"]] %>%
                         mutate(group = "Species Jaccards"),
                       dwf_rared_bray_dend_data[["segments"]] %>%
                         mutate(group = "Species Bray Curtis")) %>%
  reduce(.x = .,
         .f = full_join) %>%
  mutate(comparison = case_when(str_detect(string = group,
                                           pattern = "Site Jaccards") == TRUE ~ "Site",
                                str_detect(string = group,
                                           pattern = "Site Bray Curtis") == TRUE ~ "Site" ,
                                str_detect(string = group,
                                           pattern = "Species Jaccards") == TRUE ~ "Species",
                                str_detect(string = group,
                                           pattern = "Species Bray Curtis") == TRUE ~ "Species"),
         dist_metric = case_when(str_detect(string = group,
                                            pattern = "Site Jaccards") == TRUE ~ "Jaccards",
                                 str_detect(string = group,
                                            pattern = "Site Bray Curtis") == TRUE ~ "Bray Curtis" ,
                                 str_detect(string = group,
                                            pattern = "Species Jaccards") == TRUE ~ "Jaccards",
                                 str_detect(string = group,
                                            pattern = "Species Bray Curtis") == TRUE ~ "Bray Curtis"))
# Plot line segments and add labels
dend_plot <- clust_segments %>%
  ggplot() + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = clust_labels %>%
              drop_na(), 
            aes(x = x, 
                y = y, 
                label = label,
                colour = group_legend),
            hjust = 1, angle = 90, size = 3) +
  ylim(-0.15, 1) +
  labs(x = element_blank(),
       y = "Height",
       title = "Complete-linkage clustering",
       colour = element_blank()) +
  theme(axis.line.x = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size = 14)) +
  facet_wrap(vars(group))
dend_plot
ggsave(here("eggshells", 
            "microbiome",
            "outputs",
            "dend_plot.jpeg"))

#do kmeans clustering to see if data will group by site or species
#kmeans needs to be done on euclidean distances -need to compute a rectangular data table with n rows by PCoA
#within-cluster sum of squares is the measure of variability of the observations within each cluster (smaller SS = more compact cluster)
kmeans_rez <- list()
pcoa_yewa_rared_jacc <- cmdscale(yewa_rared_jacc, #is metric so don't need to change
                                 k = (length(yewa_nests)-1),
                                 eig = TRUE)
kmeans_rez[["yewa_rared_jacc"]] <- kmeans(pcoa_yewa_rared_jacc$points,
       centers = 2,
       nstart = 100)

pcoa_yewa_rared_bray <- cmdscale(sqrt(yewa_rared_bray),  #square root transform bray curtis so it becomes metric
                                 k = (length(yewa_nests)-1),
                                 eig = TRUE)
kmeans_rez[["yewa_rared_bray"]] <- kmeans(pcoa_yewa_rared_bray$points,
                                 centers = 2,
                                 nstart = 100)

pcoa_dwf_rared_jacc <- cmdscale(dwf_rared_jacc,
                                 k = (length(dwf_nests)-1),
                                 eig = TRUE)
kmeans_rez[["dwf_rared_jacc"]] <- kmeans(pcoa_dwf_rared_jacc$points,
                                 centers = 2,
                                 nstart = 100)

pcoa_dwf_rared_bray <- cmdscale(sqrt(dwf_rared_bray),
                                 k = (length(dwf_nests)-1),
                                 eig = TRUE)
kmeans_rez[["dwf_rared_bray"]] <- kmeans(pcoa_dwf_rared_bray$points,
                                 centers = 2,
                                 nstart = 100)

#add k-means cluster groupings to dendrogram
j <- names(kmeans_rez)
k_clust <- list()
for (i in j) {
  k_clust[[i]] <- kmeans_rez[[i]]$cluster %>%
    as.data.frame() %>%
    rownames_to_column(var = "nest_id") %>%
    rename("kmeans_clust" = ".") %>%
    mutate(cluster = i)
}
clust_labels_w_k <- k_clust %>% 
  reduce(.x = .,
         .f = full_join) %>%
  mutate(group = case_when(str_detect(string = cluster,
                                      pattern = "yewa_rared_jacc") == TRUE ~ "Site Jaccards",
                           str_detect(string = cluster,
                                      pattern = "yewa_rared_bray") == TRUE ~ "Site Bray Curtis" ,
                           str_detect(string = cluster,
                                      pattern = "dwf_rared_jacc") == TRUE ~ "Species Jaccards",
                           str_detect(string = cluster,
                                      pattern = "dwf_rared_bray") == TRUE ~ "Species Bray Curtis")) %>%
  select(-cluster) %>%
  full_join(x =.,
            y = clust_labels,
            by = c("group",
                   "nest_id" = "label"))  %>%
  mutate(comparison = case_when(str_detect(string = group,
                                           pattern = "Site Jaccards") == TRUE ~ "Site",
                                str_detect(string = group,
                                           pattern = "Site Bray Curtis") == TRUE ~ "Site" ,
                                str_detect(string = group,
                                           pattern = "Species Jaccards") == TRUE ~ "Species",
                                str_detect(string = group,
                                           pattern = "Species Bray Curtis") == TRUE ~ "Species"),
         dist_metric = case_when(str_detect(string = group,
                                            pattern = "Site Jaccards") == TRUE ~ "Jaccards",
                                 str_detect(string = group,
                                            pattern = "Site Bray Curtis") == TRUE ~ "Bray Curtis" ,
                                 str_detect(string = group,
                                            pattern = "Species Jaccards") == TRUE ~ "Jaccards",
                                 str_detect(string = group,
                                            pattern = "Species Bray Curtis") == TRUE ~ "Bray Curtis"))
clust_segments %>%
  filter(comparison == "Site") %>%
  ggplot() + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = clust_labels %>%
              filter(comparison == "Site"), 
            aes(x = x, 
                y = y, 
                label = label,
                colour = group_legend),
            hjust = 1, angle = 90, size = 3) +
  geom_point(data = clust_labels_w_k %>%
               filter(comparison == "Site"),
             aes(x = x,
                 y = -0.3,
                 shape = as.character(kmeans_clust)),
             size = 3) +
  ylim(-0.4, 1) +
  labs(x = element_blank(),
       y = "Height",
       title = "Complete-linkage clustering site",
       colour = "Sample group",
       shape = "K-means\ncluster") +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size = 20)) +
  facet_wrap(vars(dist_metric),
             nrow = 2,
             strip.position = "right") +
  scale_colour_manual(values = c("#F8766D","#7CAE00"))

ggsave(here("eggshells", 
            "microbiome",
            "outputs",
            "dend_plot_kmean_site.jpeg"))

clust_segments %>%
  filter(comparison == "Species") %>%
  ggplot() + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = clust_labels %>%
              filter(comparison == "Species"), 
            aes(x = x, 
                y = y, 
                label = label,
                colour = group_legend),
            hjust = 1, angle = 90, size = 3) +
  geom_point(data = clust_labels_w_k %>%
               filter(comparison == "Species"),
             aes(x = x,
                 y = -0.3,
                 shape = as.character(kmeans_clust)),
             size = 3) +
  ylim(-0.4, 1) +
  labs(x = element_blank(),
       y = "Height",
       title = "Complete-linkage clustering species",
       colour = "Sample group",
       shape = "K-means\ncluster") + 
  theme_classic() +
  theme(axis.line.x = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size = 20)) +
  facet_wrap(vars(dist_metric),
             nrow = 2,
             strip.position = "right") +
  scale_colour_manual(values = c("#00A9FF","#F8766D"))

ggsave(here("eggshells", 
            "microbiome",
            "outputs",
            "dend_plot_kmean_species.jpeg"))

#make NMDS plots to visualize the data structure

beta <- list(yewa_rared_jacc,
             yewa_rared_bray,
             dwf_rared_jacc,
             dwf_rared_bray)
names(beta) <- c("yewa_rared_jacc",
                 "yewa_rared_bray",
                 "dwf_rared_jacc",
                 "dwf_rared_bray")
nmds_beta <- beta %>%
  map(.x = .,
      .f = metaMDS,
      k = 2) #note stress is really low, probs insuficient data

nmds_beta_scores <- nmds_beta %>%
  map(.x = .,
      .f = scores) %>%
  map(.x = .,
      .f = as.data.frame) %>%
  map(.x = .,
      .f = rownames_to_column,
      var = "nest_id") 

j <- c("yewa_rared_jacc",
       "yewa_rared_bray",
       "dwf_rared_jacc",
       "dwf_rared_bray")
hulls <- list()
centroids <- list()
for (i in j) {
  nmds_beta_scores[[i]] <- nmds_beta_scores[[i]] %>%
    mutate(item = i,
           stress = nmds_beta[[i]]$stress,
           item_type = "points") %>%
    full_join(.,
              pilot_metadata) %>%
    drop_na()
  
  #calculate hulls
  grp <- if_else(condition = str_detect(i,
                                        pattern = "dwf") == TRUE,
                 true = "host_species",
                 false = "site")
  
  hulls[[i]] <- nmds_beta_scores[[i]] %>%
    group_by(.data[[grp]]) %>%
    slice(chull(NMDS1, NMDS2)) %>%
    mutate(item_type = "hull",
           item = i)
  
  #calculate centroid
  centroids[[i]] <- nmds_beta_scores[[i]] %>%
    group_by(.data[[grp]]) %>%
    summarize(NMDS1 = mean(NMDS1),
              NMDS2 = mean(NMDS2)) %>%
    mutate(item = i)
  }

nmds_beta_scores_plotting <- nmds_beta_scores %>% 
  reduce(.x = .,
         .f = full_join) %>%
  mutate(var_of_interest = if_else(condition = str_detect(item,
                                                          pattern = "dwf") == TRUE,
                                   true = host_species,
                                   false = site),
         comparison = case_when(str_detect(string = item,
                                           pattern = "yewa_rared_jacc") == TRUE ~ "Site",
                                str_detect(string = item,
                                           pattern = "yewa_rared_bray") == TRUE ~ "Site" ,
                                str_detect(string = item,
                                           pattern = "dwf_rared_jacc") == TRUE ~ "Species",
                                str_detect(string = item,
                                           pattern = "dwf_rared_bray") == TRUE ~ "Species"),
         dist_metric = case_when(str_detect(string = item,
                                            pattern = "yewa_rared_jacc") == TRUE ~ "Jaccards",
                                 str_detect(string = item,
                                            pattern = "yewa_rared_bray") == TRUE ~ "Bray Curtis" ,
                                 str_detect(string = item,
                                            pattern = "dwf_rared_jacc") == TRUE ~ "Jaccards",
                                 str_detect(string = item,
                                            pattern = "dwf_rared_bray") == TRUE ~ "Bray Curtis"),
         item = case_when(str_detect(string = item,
                                     pattern = "yewa_rared_jacc") == TRUE ~ "Site Jaccards",
                          str_detect(string = item,
                                     pattern = "yewa_rared_bray") == TRUE ~ "Site Bray Curtis" ,
                          str_detect(string = item,
                                     pattern = "dwf_rared_jacc") == TRUE ~ "Species Jaccards",
                          str_detect(string = item,
                                     pattern = "dwf_rared_bray") == TRUE ~ "Species Bray Curtis"))

nmds_hulls <- hulls %>%
  reduce(.x =.,
         .f = full_join) %>%
  mutate(var_of_interest = if_else(condition = str_detect(item,
                                                          pattern = "dwf") == TRUE,
                                   true = host_species,
                                   false = site),
         comparison = case_when(str_detect(string = item,
                                           pattern = "yewa_rared_jacc") == TRUE ~ "Site",
                                str_detect(string = item,
                                           pattern = "yewa_rared_bray") == TRUE ~ "Site" ,
                                str_detect(string = item,
                                           pattern = "dwf_rared_jacc") == TRUE ~ "Species",
                                str_detect(string = item,
                                           pattern = "dwf_rared_bray") == TRUE ~ "Species"),
         dist_metric = case_when(str_detect(string = item,
                                            pattern = "yewa_rared_jacc") == TRUE ~ "Jaccards",
                                 str_detect(string = item,
                                            pattern = "yewa_rared_bray") == TRUE ~ "Bray Curtis" ,
                                 str_detect(string = item,
                                            pattern = "dwf_rared_jacc") == TRUE ~ "Jaccards",
                                 str_detect(string = item,
                                            pattern = "dwf_rared_bray") == TRUE ~ "Bray Curtis"),
         item = case_when(str_detect(string = item,
                                     pattern = "yewa_rared_jacc") == TRUE ~ "Site Jaccards",
                          str_detect(string = item,
                                     pattern = "yewa_rared_bray") == TRUE ~ "Site Bray Curtis" ,
                          str_detect(string = item,
                                     pattern = "dwf_rared_jacc") == TRUE ~ "Species Jaccards",
                          str_detect(string = item,
                                     pattern = "dwf_rared_bray") == TRUE ~ "Species Bray Curtis"))

nmds_centroids <- centroids %>%
  reduce(.x =.,
         .f = full_join) %>%
  mutate(var_of_interest = if_else(condition = str_detect(item,
                                                          pattern = "dwf") == TRUE,
                                   true = host_species,
                                   false = site),
         comparison = case_when(str_detect(string = item,
                                           pattern = "yewa_rared_jacc") == TRUE ~ "Site",
                                str_detect(string = item,
                                           pattern = "yewa_rared_bray") == TRUE ~ "Site" ,
                                str_detect(string = item,
                                           pattern = "dwf_rared_jacc") == TRUE ~ "Species",
                                str_detect(string = item,
                                           pattern = "dwf_rared_bray") == TRUE ~ "Species"),
         dist_metric = case_when(str_detect(string = item,
                                            pattern = "yewa_rared_jacc") == TRUE ~ "Jaccards",
                                 str_detect(string = item,
                                            pattern = "yewa_rared_bray") == TRUE ~ "Bray Curtis" ,
                                 str_detect(string = item,
                                            pattern = "dwf_rared_jacc") == TRUE ~ "Jaccards",
                                 str_detect(string = item,
                                            pattern = "dwf_rared_bray") == TRUE ~ "Bray Curtis"),
         item = case_when(str_detect(string = item,
                                     pattern = "yewa_rared_jacc") == TRUE ~ "Site Jaccards",
                          str_detect(string = item,
                                     pattern = "yewa_rared_bray") == TRUE ~ "Site Bray Curtis" ,
                          str_detect(string = item,
                                     pattern = "dwf_rared_jacc") == TRUE ~ "Species Jaccards",
                          str_detect(string = item,
                                     pattern = "dwf_rared_bray") == TRUE ~ "Species Bray Curtis"))

#spe_site
#base plot
plot_nmds_site <- nmds_beta_scores_plotting %>%
  filter(comparison == "Site") %>%
  ggplot(data = .,
         aes(x = NMDS1,
             y = NMDS2,
             colour = var_of_interest)) +
  geom_point() +
  #stat_ellipse() + # alternative to hull
  facet_wrap(vars(dist_metric),
             nrow = 2,
             strip.position = "right") +
    labs(colour = "Sample group",
         title = "Site") +
  theme_classic() +
  scale_colour_manual(values = c("#F8766D","#7CAE00"))

plot_nmds_spe <- nmds_beta_scores_plotting %>%
  filter(comparison == "Species") %>%
ggplot(data = .,
       aes(x = NMDS1,
           y = NMDS2,
           colour = var_of_interest)) +
  geom_point() +
  #stat_ellipse() + # alternative to hull
  facet_wrap(vars(dist_metric),
              nrow = 2,
              strip.position = "right") +
  labs(colour = "Sample group",
       title = "species") +
  theme_classic() +
  scale_colour_manual(values = c("#00A9FF","#F8766D")) 

#plot hulls
plot_nmds_site +
  aes(fill = var_of_interest) +
  geom_polygon(data = nmds_hulls %>%
                 filter(comparison == "Site"),
               alpha = 0.25,
               show.legend = FALSE) +
  geom_point(data = nmds_centroids %>%
               filter(comparison == "Site"),
             size = 3, 
             shape = 21,
             colour = "black",
             show.legend = FALSE) +
  guides(fill = "none") +
  theme(text = element_text(size = 20)) +
  scale_colour_manual(values = c("#F8766D","#7CAE00")) +
  scale_fill_manual(values = c("#F8766D","#7CAE00"))
ggsave(here("eggshells", 
            "microbiome",
            "outputs",
            "nmds_site.jpeg"))

plot_nmds_spe +
  aes(fill = var_of_interest) +
  geom_polygon(data = nmds_hulls %>%
                 filter(comparison == "Species"),
               alpha = 0.25,
               show.legend = FALSE) +
  geom_point(data = nmds_centroids %>%
               filter(comparison == "Species"),
             size = 3, 
             shape = 21,
             colour = "black",
             show.legend = FALSE) +
  guides(fill = "none") +
  theme(text = element_text(size = 20)) +
  scale_colour_manual(values = c("#00A9FF","#F8766D")) +
  scale_fill_manual(values = c("#00A9FF","#F8766D"))
ggsave(here("eggshells", 
            "microbiome",
            "outputs",
            "nmds_spe.jpeg"))

# look at multivariate dispersion and centroids of the species and host groups

#look at dispersion
j <- c("yewa_rared_jacc",
       "yewa_rared_bray",
       "dwf_rared_jacc",
       "dwf_rared_bray")
bdisp_anov <- list()
bdisp_tukey <- list()
for (i in j) {
  grp <- if_else(condition = str_detect(i,
                                        pattern = "dwf") == TRUE,
                 true = "host_species",
                 false = "site")
  #make sure groups are in right order
  x <- nmds_beta_scores[[i]] %>%
    right_join(.,
               beta[[i]] %>%
                 as.matrix() %>%
                 as.data.frame() %>%
                 rownames_to_column(var = "nest_id")) %>%
    pull(.data[[grp]]) %>%
    as.factor()
    
  y <- betadisper(d = beta[[i]],
                  group = x,
                  type = "centroid") #use group centroid
  
  bdisp_anov[[i]] <- anova(y)
  
  bdisp_tukey[[i]] <- TukeyHSD(y)

}
bdisp_anov # no significant difference in dispersion
bdisp_tukey # no significant difference in dispersion # note: tukey test is weaker than ANOVA
#NOTE permdisp and permanova are sensitive to unbalanced desing (might want to address this is future work - maybe subsample groups and run a bunch of tests with balanced groups)

#extract results
j <- c("yewa_rared_jacc",
       "yewa_rared_bray",
       "dwf_rared_jacc",
       "dwf_rared_bray")
disp_rez <- list()
for (i in j) {
  disp_rez[[i]] <- tibble(f_value_a = bdisp_anov[[i]]$`F value`[[1]],
                          p_value_a = bdisp_anov[[i]]$`Pr(>F)`[[1]],
                          diff_t = bdisp_tukey[[i]]$group[[1]] ,
                          p_value_t = bdisp_tukey[[i]]$group[[4]],
                          group = i)
}
disp_rez <- disp_rez %>%
  reduce(.x = .,
         .f = full_join)

#look at difference in centroids - permanova 
j <- c("yewa_rared_jacc",
       "yewa_rared_bray",
       "dwf_rared_jacc",
       "dwf_rared_bray")
cent_perm <- list()
for (i in j) {
  grp <- if_else(condition = str_detect(i,
                                        pattern = "dwf") == TRUE,
                 true = "host_species",
                 false = "site")
  #make sure groups are in right order
  x <- nmds_beta_scores[[i]] %>%
    right_join(.,
               beta[[i]] %>%
                 as.matrix() %>%
                 as.data.frame() %>%
                 rownames_to_column(var = "nest_id")) %>%
    pull(.data[[grp]]) %>%
    as.factor()
  
  cent_perm[[i]] <- adonis2(beta[[i]] ~ x)
}
cent_perm # no significant difference
#extract results
j <- c("yewa_rared_jacc",
       "yewa_rared_bray",
       "dwf_rared_jacc",
       "dwf_rared_bray")
perm_rez <- list()
for (i in j) {
  perm_rez[[i]] <- tibble(
    df_p = paste(cent_perm[[i]]$Df[[1]], "and", cent_perm[[i]]$Df[[2]]), #df
    r2_p = cent_perm[[i]]$R2[[1]], #R2
    pseudo_F_p = cent_perm[[i]]$F[[1]], #pseudo-F 
    p_value_p = cent_perm[[i]]$`Pr(>F)`[[1]], #p-value
    group = i)
}
perm_rez <- perm_rez %>%
  reduce(.x = .,
         .f = full_join)

perm_disp_rez <- full_join(disp_rez,
                           perm_rez) %>%
  mutate(dist = case_when(str_detect(string = group,
                                     pattern = "jacc") == TRUE ~ "Jaccards",
                          str_detect(string = group,
                                     pattern = "bray") == TRUE ~ "Bray Curtis"),
         comparison = case_when(str_detect(string = group,
                                           pattern = "yewa") == TRUE ~ "Site",
                                str_detect(string = group,
                                           pattern = "dwf") == TRUE ~ "Species"),
         across(.cols = where(is.numeric),
                .fns = ~ round(x = .x,
                               digits = 3)
                )) %>%
  select(-group)


perm_disp_rez_tbl <- perm_disp_rez %>% 
  gt(groupname_col = "dist") %>%
  tab_style(
    style = list(cell_borders(
      sides = c("bottom"),
      style = "solid"),
      cell_text(weight = "bold",
                size = "x-small")),
    locations = cells_row_groups()) %>%
  tab_style(style = list(cell_borders(sides = "all",
                                      style = "hidden"),
                         cell_text(size = "x-small",
                                   align = "left")),
            locations = cells_body()) %>%
  tab_style(style = list(cell_text(align = "left")),
            locations = cells_title()) %>%
  tab_spanner(label = "Dispersion ANOVA",
              columns = c(ends_with("_a"))) %>% 
  tab_spanner(label = "Dispersion Tukey",
              columns = c(ends_with("_t"))) %>%
  tab_spanner(label = "PERMANOVA",
              columns = c(ends_with("_p"))) %>%
  tab_style(style = cell_text(align = "left",
                              weight = "bolder",
                              size = "small",
                              whitespace = "break-spaces"),
            locations = list(cells_column_labels(),
                             cells_column_spanners())) %>%
  tab_header(title = "Difference in centroids and dispersion of species and sites") %>%
  cols_label(f_value_a = md("F-value"),
             p_value_a = md("*p*-value"),
             diff_t = md("Difference in means"),
             p_value_t = md("*p*-value"),
             df_p = md("Degrees of freedom"),
             r2_p = md("R^2^"),
             pseudo_F_p = md("Pseudo F-value"),
             p_value_p = md("*p*-value"),
             comparison = md("Comparison")) %>%
  cols_move_to_start(comparison)
perm_disp_rez_tbl
gtsave(data = perm_disp_rez_tbl,
       filename = here("eggshells",
                       "microbiome",
                       "outputs",
                       "perm_disp_rez_tbl.png"))

#try different visualization method - https://www.youtube.com/watch?v=jLVKJ_n6Qd0&t=119s
j <- c("yewa_rared_jacc",
       "yewa_rared_bray",
       "dwf_rared_jacc",
       "dwf_rared_bray")
beta_jitter <- list()
for (i in j) {
  beta_jitter[[i]] <- beta[[i]] %>%
    as.matrix() %>%
    replace_triangle(.,
                     #replace upper triangle
                     triangle = "upper",
                     #replace with NA
                     by = NA,
                     #exclude diagonal
                     diagonal = FALSE) %>%
    as.data.frame() %>%
    rename("nest_id" = "rowname") %>%
    pivot_longer(cols = - nest_id,
                 names_to = "nest_id_compared",
                 values_to = "dist") %>%
    drop_na() %>%
    mutate(item = i) %>%
    left_join(.,
              pilot_metadata) %>%
    left_join(.,
              pilot_metadata %>% 
                rename_with(~ paste0(., "_compared")))
}
beta_jitter <- beta_jitter %>%
  reduce(.x = .,
         .f = full_join) %>%
  mutate(item_sim = case_when(str_detect(string = item, 
                                         pattern = "yewa") == TRUE &
                                #site == site_compared ~ "same",
                                site == site_compared &
                                site == "dwf"~ "dwf",
                              str_detect(string = item, 
                                         pattern = "yewa") == TRUE &
                                site == site_compared &
                                site == "pcc"~ "pcc",
                              str_detect(string = item, 
                                         pattern = "yewa") == TRUE &
                                site != site_compared ~ "different",
                              str_detect(string = item, 
                                         pattern = "dwf") == TRUE &
                                #host_species == host_species_compared ~ "same",
                                host_species == host_species_compared &
                                host_species == "YEWA" ~ "YEWA",
                              str_detect(string = item, 
                                         pattern = "dwf") == TRUE &
                                host_species == host_species_compared &
                                host_species == "RWBL" ~ "RWBL",
                              str_detect(string = item, 
                                         pattern = "dwf") == TRUE &
                                host_species != host_species_compared ~ "different")) %>%
  mutate(comparison = case_when(str_detect(string = item,
                                           pattern = "yewa_rared_jacc") == TRUE ~ "Site",
                                str_detect(string = item,
                                           pattern = "yewa_rared_bray") == TRUE ~ "Site" ,
                                str_detect(string = item,
                                           pattern = "dwf_rared_jacc") == TRUE ~ "Species",
                                str_detect(string = item,
                                           pattern = "dwf_rared_bray") == TRUE ~ "Species"),
         dist_metric = case_when(str_detect(string = item,
                                            pattern = "yewa_rared_jacc") == TRUE ~ "Jaccards",
                                 str_detect(string = item,
                                            pattern = "yewa_rared_bray") == TRUE ~ "Bray Curtis" ,
                                 str_detect(string = item,
                                            pattern = "dwf_rared_jacc") == TRUE ~ "Jaccards",
                                 str_detect(string = item,
                                            pattern = "dwf_rared_bray") == TRUE ~ "Bray Curtis"),
         item = case_when(str_detect(string = item,
                                     pattern = "yewa_rared_jacc") == TRUE ~ "Site Jaccards",
                          str_detect(string = item,
                                     pattern = "yewa_rared_bray") == TRUE ~ "Site Bray Curtis" ,
                          str_detect(string = item,
                                     pattern = "dwf_rared_jacc") == TRUE ~ "Species Jaccards",
                          str_detect(string = item,
                                     pattern = "dwf_rared_bray") == TRUE ~ "Species Bray Curtis"))

beta_jitter %>%
  filter(comparison == "Site") %>%
  ggplot(data = .,
         aes(x = item_sim,
             y = dist)) +
  geom_jitter(width = 0.25,
              colour = "gray") +
  stat_summary(fun.data = median_hilow,
               colour = "red",
               size = 1,
               fun.args = list(conf.int = 0.5)) +
  theme_classic() +
  facet_wrap(vars(dist_metric),
             nrow = 2,
             strip.position = "right"
  ) +
  labs(y = "Distance",
       x = element_blank(),
       title = "site") +
  theme(text = element_text(size = 20))

ggsave(here("eggshells", 
            "microbiome",
            "outputs",
            "multivar_var_site.jpeg"))

beta_jitter %>%
  filter(comparison == "Species") %>%
  ggplot(data = .,
         aes(x = item_sim,
             y = dist)) +
  geom_jitter(width = 0.25,
              colour = "gray") +
  stat_summary(fun.data = median_hilow,
               colour = "red",
               size = 1,
               fun.args = list(conf.int = 0.5)) +
  theme_classic() +
  facet_wrap(vars(dist_metric),
             nrow = 2,
             strip.position = "right") +
  labs(y = "Distance",
       x = element_blank(),
       title = "species") +
  theme(text = element_text(size = 20))

ggsave(here("eggshells", 
            "microbiome",
            "outputs",
            "multivar_var_spe.jpeg"))

## do an RDA to see the effect of clutch initiation date and host species or site on beta diversity
#counts/pa = clutch initiation date * spe or site

#get rarefied data
j <- 1:100
asv_sample_rared_long <- list()
for (i in j) {
  asv_sample_rared_long[[i]] <- asv_sample_rared[[i]] %>%
    rownames_to_column(var = "nest_id") %>%
    pivot_longer(cols = -nest_id,
                 names_to = "asv",
                 values_to = "reads") %>%
    mutate(rared_it = i) 
}
asv_sample_rared_mean <- asv_sample_rared_long %>%
  reduce(.x = .,
         .f = full_join) %>%
  group_by(nest_id,
           asv) %>%
  summarise(reads = mean(reads)) %>%
  pivot_wider(names_from = "asv",
              values_from = "reads")

asv_sample_rared_mean_resp <- list(asv_sample_rared_mean %>%
                                filter(nest_id %in% dwf_nests) %>%
                                  column_to_rownames(var = "nest_id") %>%
                                  select_if(~sum(.) > 0),
                                asv_sample_rared_mean %>%
                                  filter(nest_id %in% dwf_nests) %>%
                                  column_to_rownames(var = "nest_id") %>%
                                  decostand(method = "pa") %>%
                                  select_if(~sum(.) > 0),
                                asv_sample_rared_mean %>%
                                  filter(nest_id %in% yewa_nests) %>%
                                  column_to_rownames(var = "nest_id") %>%
                                  select_if(~sum(.) > 0),
                                asv_sample_rared_mean %>%
                                  filter(nest_id %in% yewa_nests) %>%
                                  column_to_rownames(var = "nest_id") %>%
                                  decostand(method = "pa") %>%
                                  select_if(~sum(.) > 0)) %>%
  map(.x = .,
      .f = decostand,
      method = "hellinger") #data needs to be hellinger transformed before piping to rda

names(asv_sample_rared_mean_resp) <- c("dwf_count",
                                       "dwf_pa",
                                       "yewa_count",
                                       "yewa_pa")
j <- c("dwf_count",
       "dwf_pa",
       "yewa_count",
       "yewa_pa")
asv_sample_rared_mean_explan <- list()
for (i in j) {
  grp <- if_else(condition = str_detect(i,
                                        pattern = "dwf"),
                 true = "host_species",
                 false = "site")
  asv_sample_rared_mean_explan[[i]] <- asv_sample_rared_mean_resp[[i]] %>%
    rownames_to_column(var = "nest_id") %>%
    select("nest_id") %>%
    left_join(x = .,
              y = pilot_swab_nests_reps %>%
                mutate(host_species = as.factor(host_species),
                       site = as.factor(site)) %>%
                select(nest_id,
                       all_of(grp),
                       host_clutch_initiation_day_of_year) %>%
                distinct())
  }

j <- c("dwf_count",
       "dwf_pa")
tbrda_rez <- list()
for (i in j) {
  #remove samples missing data
  x <- asv_sample_rared_mean_explan[[i]] %>% 
    drop_na()
  nests <- x %>% pull(nest_id)
  y <- asv_sample_rared_mean_resp[[i]] %>%
    rownames_to_column(var = "nest_id") %>%
    filter(nest_id %in% nests) %>%
    column_to_rownames(var = "nest_id")
  
  tbrda_rez[[i]] <- rda(y ~ host_species*host_clutch_initiation_day_of_year,
                   data = x,
                   na.action = na.omit) 
  tbrda_rez[[paste0(i, "_coefs")]] <- tbrda_rez[[i]] %>% coef(.) 
  tbrda_rez[[paste0(i, "_r2")]] <- tbrda_rez[[i]] %>% RsquareAdj(.) #multiply adjusted r2 by proportion of accumulated constrained eigen values 
  tbrda_rez[[paste0(i, "_anovrda")]] <- tbrda_rez[[i]] %>% anova(.,
                                                                 permutations = how(nperm = 999))
  tbrda_rez[[paste0(i, "_anovaaxis")]] <- tbrda_rez[[i]] %>% anova(.,
                                                                 by = "axis", #variation in explanatory variables in fewer dimensions
                                                                 permutations = how(nperm = 999))
  tbrda_rez[[paste0(i, "_anovaterm")]] <- tbrda_rez[[i]] %>% anova(.,
                                                                   by = "terms",
                                                                   permutations = how(nperm = 999))
}
j <- c("yewa_count",
       "yewa_pa")
for (i in j) {
  #remove samples missing data
  x <- asv_sample_rared_mean_explan[[i]] %>% 
    drop_na()
  nests <- x %>% pull(nest_id)
  y <- asv_sample_rared_mean_resp[[i]] %>%
    rownames_to_column(var = "nest_id") %>%
    filter(nest_id %in% nests) %>%
    column_to_rownames(var = "nest_id")
  
  tbrda_rez[[i]] <- rda(y ~ site*host_clutch_initiation_day_of_year,
                        data = x,
                        na.action = na.omit) 
  tbrda_rez[[paste0(i, "_coefs")]] <- tbrda_rez[[i]] %>% coef(.) 
  tbrda_rez[[paste0(i, "_r2")]] <- tbrda_rez[[i]] %>% RsquareAdj(.) 
  tbrda_rez[[paste0(i, "_anovrda")]] <- tbrda_rez[[i]] %>% anova(.,
                                                              permutations = how(nperm = 999))
  tbrda_rez[[paste0(i, "_anovaaxis")]] <- tbrda_rez[[i]] %>% anova(.,
                                                                 by = "axis", #variation in explanatory variables in fewer dimensions
                                                                 permutations = how(nperm = 999))
  tbrda_rez[[paste0(i, "_anovaterm")]] <- tbrda_rez[[i]] %>% anova(.,
                                                                   by = "terms",
                                                                   permutations = how(nperm = 999))
}
tbrda_rez
#constrained variation is proportion variation explained by explanatory variables (baised R^2)
#how do constrained egein values vs unconstrained eigen values do w explaining the variance?
# dwf_count RDA is slightly larger, dwf_pa the PCA is larger, yewa_count PCA is larger, yewa_pa PCA is slightly larger
#significant for site

#make table for results
tbrda_rez_4tbl <- list(tibble(distance = "Bray-curtis",
                              variance_type = c("Total",
                                                "Constrained",
                                                "Unconstrained"),
                              variance_dwf = c(0.3414,
                                               0.1478 ,
                                               0.1937),
                              prop_var_explained_dwf = c(1.000,
                                                     0.4328,
                                                     0.5672),
                              summary_stat_type_dwf = c("adjusted R2", # measures the strength of the relationship between Y and X, but applies a correction of the R2 to take into account the number of explanatory variables
                                                    "F-statistic", #overall test of significance of an RDA by comparing the computed model to a null model
                                                    "p-value"),
                              summary_stat_dwf = c(0.1896458,
                                                   1.7801,
                                                   0.068),
                              term_dwf = c("species",
                                           "clutch initiation day of year",
                                           "species:clutch initiation day of year"),
                              term_f_stat_dwf = c(3.5135,
                                                  1.1253,
                                                  0.7015),
                              term_p_val_dwf = c(0.016,
                                                 0.299,
                                                 0.668)),
                       tibble(distance = "Jaccards",
                              variance_type = c("Total",
                                                "Constrained",
                                                "Unconstrained"),
                              variance_dwf = c(0.6032,
                                               0.1956,
                                               0.4076),
                              prop_var_explained_dwf = c(1.000,
                                                     0.3243,
                                                     0.6757),
                              summary_stat_type_dwf = c("adjusted R2", # measures the strength of the relationship between Y and X, but applies a correction of the R2 to take into account the number of explanatory variables
                                                    "F-statistic", #overall test of significance of an RDA by comparing the computed model to a null model
                                                    "p-value"),
                              summary_stat_dwf = c(0.03464573,
                                                   1.1196,
                                                   0.166),
                              term_dwf = c("species",
                                           "clutch initiation day of year",
                                           "species:clutch initiation day of year"),
                              term_f_stat_dwf = c(1.4011,
                                                  0.9229,
                                                  1.0349),
                              term_p_val_dwf = c(0.044,
                                                 0.667,
                                                 0.429)),
                       full_join(tibble(distance = "Bray-curtis",
                              variance_type = c("Total",
                                                "Constrained",
                                                "Unconstrained"),
                              variance_yewa = c(0.2387,
                                                0.1079,
                                                0.1308),
                              prop_var_explained_yewa = c(1.000,
                                                     0.4520,
                                                     0.5480),
                              summary_stat_type_yewa = c("adjusted R2", # measures the strength of the relationship between Y and X, but applies a correction of the R2 to take into account the number of explanatory variables
                                                    "F-statistic", #overall test of significance of an RDA by comparing the computed model to a null model
                                                    "p-value"),
                              summary_stat_yewa = c(0.2172117,
                                                    1.9249,
                                                    0.006),
                              term_yewa = c("site",
                                       "clutch initiation day of year",
                                       "site:clutch initiation day of year"),
                              term_f_stat_yewa = c(3.0089,
                                                   1.4737,
                                                   1.2923),
                              term_p_val_yewa = c(0.001,
                                                  0.140,
                                                  0.206)),
                       tibble(distance = "Jaccards",
                              variance_type = c("Total",
                                                "Constrained",
                                                "Unconstrained"),
                              variance_yewa = c(0.6097,
                                                0.1819,
                                                0.4278),
                              prop_var_explained_yewa = c(1.000,
                                                     0.2984,
                                                     0.7016),
                              summary_stat_type_yewa = c("adjusted R2", # measures the strength of the relationship between Y and X, but applies a correction of the R2 to take into account the number of explanatory variables
                                                    "F-statistic", #overall test of significance of an RDA by comparing the computed model to a null model
                                                    "p-value"),
                              summary_stat_yewa = c(-0.002321499,
                                                    0.9923,
                                                    0.522),
                              term_yewa = c("site",
                                       "clutch initiation day of year",
                                       "site:clutch initiation day of year"),
                              term_f_stat_yewa = c(1.1143,
                                                   0.8660,
                                                   0.9965),
                              term_p_val_yewa = c(0.16,
                                                  0.89,
                                                  0.52)))
       
       ) %>%
  reduce(.x = .,
         .f = full_join) %>%
  mutate(across(.cols = where(is.numeric),
                .fns = ~ round(x = .x,
                               digits = 3))) %>%
  unite(col = model_summary_yewa,
        summary_stat_type_yewa,
        summary_stat_yewa,
        sep = ": ") %>%
  unite(col = model_summary_dwf,
        summary_stat_type_dwf,
        summary_stat_dwf,
        sep = ": ")


tbrda_rez_tbl <- tbrda_rez_4tbl %>%
  gt(groupname_col = "distance") %>%
  tab_style(style = list(cell_borders(sides = "all",
                                      style = "hidden"),
                         cell_text(size = "x-small",
                                   align = "left")),
            locations = cells_body()) %>%
  tab_spanner(label = "Species",
              columns = c(ends_with("_dwf"))) %>% 
  tab_spanner(label = "Site",
              columns = c(ends_with("_yewa"))) %>%
  tab_style(style = list(cell_text(align = "left")),
            locations = cells_title()) %>%
  tab_style(style = cell_text(align = "left",
                              weight = "bolder",
                              size = "small",
                              whitespace = "break-spaces"),
            locations = list(cells_column_labels())) %>%
  tab_header(title = "tb-RDA for relationship between clutch initiation date and sample group") %>%
  cols_label(variance_type = md("Variance type"),
             variance_dwf = md("Variance"),
             variance_yewa = md("Variance"),
             prop_var_explained_dwf = md("Proportion variance explained"),
             prop_var_explained_yewa = md("Proportion variance explained"),
             model_summary_dwf = md("Model sumary statistics"),
             model_summary_yewa = md("Model sumary statistics"),
             term_dwf = md("Model term"),
             term_yewa = md("Model term"),
             term_f_stat_dwf = md("Model term F-statistic"),
             term_f_stat_yewa = md("Model term F-statistic"),
             term_p_val_dwf = md("Model term *p* value"),
             term_p_val_yewa = md("Model term *p* value"))
tbrda_rez_tbl

gtsave(data = tbrda_rez_tbl,
       filename = here("eggshells",
                       "microbiome",
                       "outputs",
                       "tbrda_rez_tbl.png"))

#plot results
beta_jitter <- beta_jitter %>%
  left_join(.,
            pilot_swab_nests_reps %>%
              select(nest_id, 
                   host_clutch_initiation_day_of_year) %>%
              distinct()) %>% 
  left_join(.,
            pilot_swab_nests_reps %>%
              select(nest_id, 
                   host_clutch_initiation_day_of_year) %>%
              distinct() %>%
            rename_with(~ paste0(., "_compared"))) %>%
  mutate(days_between_clutch_initiation = abs(host_clutch_initiation_day_of_year - 
           host_clutch_initiation_day_of_year_compared))


ggplot() +
  aes(x = days_between_clutch_initiation,
      y = dist,
      colour = item_sim) +
  geom_point(data = beta_jitter,
             size = .5,
             alpha = 0.5) +
  geom_line(data = beta_jitter %>% 
              group_by(item,
                       days_between_clutch_initiation,
                       item_sim) %>%
              summarise(median_dist = median(dist)) %>%
              ungroup() %>%
              rename("dist" = "median_dist"),
            linewidth = 1) +
  facet_wrap(vars(item),
            # scales = "free_y"
  ) +
  labs(x = "Days between clutch initiation",
       y = "Community distance",
       colour = "Paired\nSample\nGroup",
       title = "Change in beta diversity with time") +
  theme_classic() +
  theme(text = element_text(size = 20))

ggsave(here("eggshells", 
            "microbiome",
            "outputs",
            "multivar_var_time.jpeg"))