#### Purpose: This script is for cleaning and prepping the raw data
#### Author: Anya Mueller
#### Date: August 9, 2023

##### Set up ####
## install packages if needed
# install.packages("tidyverse")
# install.packages("janitor")
# install.packages("readxl")
# install.packages("geosphere")
# install.packages("here")

## read in packages
library(tidyverse)
library(janitor)
library(readxl)
library(geosphere)
library(here)

## read in data
all_field_data <- read_csv(here("eggshells",
                                "field",
                                "processed_data",
                                "master_and_nest_observations.csv"))

pilot_micro <- read_csv(here("eggshells", 
                             "microbiome", 
                             "dummy_data",
                             "pilot_micro_dummy.csv")) %>%
  clean_names()

#### Measure distances between nests ####

#make longitude negative
all_field_data <- all_field_data %>%
  mutate(longitude = -longitude)

#make matrix - long, lat
long <- all_field_data %>% 
  select(nest_id, longitude) %>%
  distinct() %>%
  pull(longitude)
lat <- all_field_data %>%
  select(nest_id, latitude) %>%
  distinct() %>%
  pull(latitude)

long_lat <- data.frame(long,lat) %>%
  as.matrix()

#calculate shortest distance between sites in meters
distance <- distm(long_lat,
                  fun = distGeo) #highly accurate method of shortest distance between two points on an ellipsoid (default is WGS84 ellipsoid)
#add in names
colnames(distance) <- paste("distance_to_", c(all_field_data %>% 
                                                select(nest_id) %>% 
                                                distinct() %>% 
                                                pull(nest_id)))
rownames(distance) <- all_field_data %>% 
  select(nest_id) %>% 
  distinct() %>% 
  pull(nest_id)

#make into dataframe
distance_mat_long <- distance %>% 
  as.data.frame(.) %>%
  rownames_to_column(var = "nest_id") 

#export
write_csv(distance_mat_long,
          here("eggshells", 
               "microbiome",
               "processed_data",
               "distance_between_nests_m.csv"))

#### Field data ####
#how many eggs on first day of monitoring?
first_day_monitoring <- all_field_data %>% 
  group_by(nest_id) %>%
  arrange(day_of_year) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(num_host_egg_first_day = num_host_egg)

all_field_data <- full_join(all_field_data,
                            first_day_monitoring %>%
                              select(nest_id, 
                                     num_host_egg_first_day))

not_swabbed <- all_field_data %>% 
  filter(!nest_id %in% c(all_field_data %>%
                           filter(host_swab == "yes") %>%
                           pull(nest_id) %>%
                           unique()))

swab_nests <- all_field_data %>%
  filter(nest_id %in% c(all_field_data %>%
                          filter(host_swab == "yes") %>%
                          pull(nest_id) %>%
                          unique()))
#number of eggs swabbed
first_swab_n_eggs <- swab_nests %>%
  filter(host_swab == "yes") %>%
  group_by(nest_id) %>%
  arrange(day_of_year) %>%
  filter(row_number() == 1) %>%
  mutate(num_eggs_first_host_swab = num_host_egg) %>%
  ungroup() %>%
  select(nest_id,
         num_eggs_first_host_swab)

second_swab_n_eggs <- swab_nests %>%
  filter(host_swab == "yes") %>%
  group_by(nest_id) %>%
  arrange(day_of_year) %>%
  filter(row_number() == 2) %>%
  mutate(num_eggs_second_host_swab = num_host_egg)%>%
  ungroup() %>%
  select(nest_id,
         num_eggs_second_host_swab)

third_swab_n_eggs <- swab_nests %>%
  filter(host_swab == "yes") %>%
  group_by(nest_id) %>%
  arrange(day_of_year) %>%
  filter(row_number() == 3) %>%
  mutate(num_eggs_third_host_swab = num_host_egg) %>%
  ungroup() %>%
  select(nest_id,
         num_eggs_third_host_swab)

swab_nests <- list(swab_nests,
        first_swab_n_eggs,
        second_swab_n_eggs,
        third_swab_n_eggs) %>%
  reduce(full_join, 
         by = "nest_id")

#1st, 2nd, or 3rd swab

first_swab_date <- swab_nests %>%
  filter(host_swab == "yes") %>%
  group_by(nest_id) %>%
  arrange(day_of_year) %>%
  filter(row_number() == 1) %>%
  mutate(host_swab_order = "first") %>%
  ungroup() %>%
  select(nest_id,
         host_swab_order,
         date)

second_swab_date <- swab_nests %>%
  filter(host_swab == "yes") %>%
  group_by(nest_id) %>%
  arrange(day_of_year) %>%
  filter(row_number() == 2) %>%
  mutate(host_swab_order = "second") %>%
  ungroup() %>%
  select(nest_id,
         host_swab_order,
         date)

third_swab_date <- swab_nests %>%
  filter(host_swab == "yes") %>%
  group_by(nest_id) %>%
  arrange(day_of_year) %>%
  filter(row_number() == 3) %>%
  mutate(host_swab_order = "third") %>%
  ungroup() %>%
  select(nest_id,
         host_swab_order,
         date)

swab_nests <- list(first_swab_date,
                   second_swab_date,
                   third_swab_date) %>%
  reduce(full_join, 
         by = c("nest_id", "date", "host_swab_order")) %>%
  full_join(swab_nests, 
            .,
            by = c("nest_id", "date"))

#### Microbiome data ####
#initial filtering
taxo <- pilot_micro %>% 
  select(asv_name:species) %>%
  colnames()
  
pilot_micro <- pilot_micro %>% 
  unite(all_of(taxo),
        col = taxo,
        sep = "-") %>%
  pivot_longer(cols = c(contains("mul1806"), contains("cou1805")),
               names_to = "sample",
               values_to = "reads") %>%
  filter(reads > 0) %>%
  group_by(taxo) %>%
  #calculate variables to filter against
  summarize(occurrence = n(),
            max_taxo_reads = max(reads),
            total_taxo_reads = sum(reads),
            sparsity = max_taxo_reads/total_taxo_reads) %>% 
    ungroup() %>%
  right_join(.,
            pilot_micro %>% 
              unite(all_of(taxo),
                    col = taxo,
                    sep = "-",
                    remove = FALSE)) %>%
    select(-taxo) %>%
  filter(
    #filter to occurrence threshold of 3 samples containing the ASV - make sure sequence is not random chance
    occurrence > 3,
    #filter to sparsity threshold of less than 90% - make sure reads are spread out
    sparsity <= 0.9
         )

pilot_taxo <- pilot_micro %>%
  select(asv_name,
         kingdom,
         phylum,
         class,
         order,
         family,
         genus,
         species) %>%
  distinct() %>%
  #filter for bacteria only
  filter(kingdom == "Bacteria") %>%
  #fix species names
  mutate(species = str_replace(string = species,
                              pattern = "_",
                              replacement = " ")) %>%
  mutate(across(.fns = str_replace_na)) %>%
  mutate(across(.fns = str_replace,
                pattern = "NA",
                replacement = "Unknown Bacteria"))
  

pilot_asv_seq <- pilot_micro %>%
  select(asv_name, asv_sequences) %>%
  #filter to bacterial ASVs
  filter(asv_name %in% c(pilot_taxo %>%
                           pull(asv_name) %>%
                           unique())) %>%
  distinct()

pilot_sample_asv <- pilot_micro %>% 
  pivot_longer(cols = c(contains("mul1806"), contains("cou1805")),
               names_to = "sample",
               values_to = "reads") %>%
  select(asv_name,
         sample,
         reads) %>%
  #filter to bacterial ASVs
  filter(asv_name %in% c(pilot_taxo %>%
                           pull(asv_name) %>%
                           unique())) %>%
  distinct() %>%
  mutate(nest_id = str_remove(string = sample,
                               pattern = "_mul.*") %>%
           str_replace_all(string = .,
                           pattern = "^dogsample",
                           replacement = NA_character_) %>%
           str_replace_all(string = .,
                           pattern = "^water",
                           replacement = NA_character_) %>%
           str_replace_all(string = .,
                           pattern = "^x00",
                           replacement = NA_character_),
         control_type1 = str_remove(string = sample,
                                     pattern = "._mul.*") %>%
           str_replace_all(string = .,
                          pattern = "_cou.*",
                          replacement = NA_character_) %>%
           str_replace_all(string = .,
                           pattern = "^x00.*",
                           replacement = "cloacal") %>%
           str_replace_all(string = .,
                           pattern = "^a",
                           replacement = NA_character_) %>%
           str_replace_all(string = .,
                           pattern = "^e",
                           replacement = NA_character_) %>%
           str_replace_all(string = .,
                           pattern = "^j",
                           replacement = NA_character_),
         control_type2 = str_remove(string = sample,
                                     pattern = "._cou.*")  %>%
           str_replace_all(string = .,
                           pattern = "_mul.*",
                           replacement = NA_character_) %>%
           str_replace_all(string = .,
                           pattern = "^x00",
                           replacement = NA_character_) %>%
           str_replace_all(string = .,
                           pattern = "^a",
                           replacement = NA_character_) %>%
           str_replace_all(string = .,
                           pattern = "^e",
                           replacement = NA_character_) %>%
           str_replace_all(string = .,
                           pattern = "^j",
                           replacement = NA_character_)) %>%
  unite(col = control_type,
        starts_with("control_type"),
        remove = TRUE) %>%
  mutate(control_type = str_replace(string = control_type,
                                    pattern = "NA_NA",
                                    replacement = NA_character_) %>%
           str_remove(string =.,
                      pattern = "^NA_") %>%
           str_remove(string =.,
                      pattern = "_NA$"))

#pilot metadata
pilot_metadata <- swab_nests %>%
  filter(nest_id %in% 
           c(pilot_sample_asv %>% 
               pull(nest_id))) %>%
  select(nest_id,
         date,
         day_of_year,
         time,
         num_host_egg,
         host_nest_outcome,
         host_swab,
         host_species,
         site) %>%
  filter(host_swab == "yes") %>%
  group_by(nest_id) %>%
  arrange(day_of_year) %>%
  filter(row_number() == 1) %>%
  select(nest_id,
         date,
         time,
         num_host_egg,
         host_nest_outcome,
         host_species,
         site) %>%
  mutate(host_species = case_when(host_species == "YEWA" ~ "YEWA",
                                  TRUE ~ "RWBL"),
         site = case_when(host_species == "YEWA" ~ site,
                          host_species == "RWBL" ~ "dwf")) #fix dummy data for example purposes

#### Export data ####
write_csv(swab_nests, 
          file = here("eggshells", 
                      "microbiome",
                      "processed_data",
                      "swab_nests.csv"))

write_csv(not_swabbed, 
          file = here("eggshells", 
                      "microbiome",
                      "processed_data",
                      "not_swabbed.csv"))

write_csv(pilot_asv_seq,
          file = here("eggshells", 
                      "microbiome",
                      "processed_data",
                      "pilot_asv_seq.csv"))

write_csv(pilot_taxo,
          file = here("eggshells", 
                      "microbiome",
                      "processed_data",
                      "pilot_taxo.csv"))

write_csv(pilot_sample_asv,
          file = here("eggshells", 
                      "microbiome",
                      "processed_data",
                      "pilot_sample_asv.csv"))

write_csv(pilot_metadata,
          file = here("eggshells", 
                      "microbiome",
                      "processed_data",
                      "pilot_metadata.csv"))
