#### Purpose: This script is for cleaning the raw data
#### Author: Anya Mueller
#### Date: August 12, 2022

# Set up -----------------------------------------------------------------
## install packages if needed
# install.packages("tidyverse")
# install.packages("janitor")
# install.packages("readxl")
# install.packages("lubridate")
# install.packages("here")

## read in packages
library(tidyverse)
library(janitor)
library(readxl)
library(lubridate)
library(here)


## Read in data ##
master <- read_csv(here("eggshells",
                        "field",
                        "dummy_data",
                        "master_dummy.csv")) %>% 
  clean_names()

nest_log <- read_csv(here("eggshells",
                          "field",
                          "dummy_data",
                          "nest_log_dummy.csv")) %>% 
  clean_names() 

#remove old nests
old_nests <- nest_log %>% 
  mutate(old_nest = case_when(str_detect(string = notes,
                                         pattern = "(?i)old nest") == TRUE ~ "yes",
                              str_detect(string = notes,
                                         pattern = "maybe old") == TRUE ~ "no",
                              str_detect(string = notes,
                                         pattern = "confirmed old") == TRUE ~ "yes",
                              str_detect(string = notes,
                                         pattern = "nest old") == TRUE ~ "yes",
                              TRUE ~ "no")) %>%
  filter(old_nest == "yes") %>%
  pull(nest_id) %>%
  unique()

master <- master %>%
  filter(!nest_id %in% old_nests)

nest_log <- nest_log %>%
  filter(!nest_id %in% old_nests)

#isolate nests that weren't monitored
non_monitored_nests <- nest_log %>% 
  mutate(nest_monitored = case_when(str_detect(string = notes, pattern = "(?i)outside search area") == TRUE ~ "no",
                             notes == "too high to swab if parasitized, no ladder access" ~ "yes",
                             notes == "need ladder to swab" ~ "yes",
                             str_detect(string = notes, pattern = "(?i)ladder") == TRUE ~ "no",
                             notes == "Not accessible at all" ~ "no",
                             notes == "Too high to monitor" ~ "no",
                             notes == "Too high to collect data (More than 5 meters)" ~ "no",
                             notes == "nest too high to acess" ~ "no",
                             TRUE ~ "yes")) %>%
  filter(nest_monitored == "no") %>%
  pull(nest_id) %>%
  unique()

# Master nest list --------------------------------------------------------

## Clean up data values ##

master <- master %>%
  mutate(tree_species = case_when(tree_species == "dwarf raspberry" ~ "raspberry",
                                  tree_species == "dwarf raspberrry" ~ "raspberry",
                                  tree_species == "choke cherry" ~ "chokecherry",
                                  tree_species == "Red-berried elder" ~ "red-berried elder",
                                  tree_species == "Peach willow" ~ "peachleaf willow",
                                  tree_species == "peach willow" ~ "peachleaf willow",
                                  tree_species == "Box elder" ~ "box elder",
                                  tree_species == "N/A" ~ NA_character_,
                                  tree_species == "na" ~ NA_character_, 
                                  tree_species == "NA" ~ NA_character_, 
                                  tree_species == "Cat tail" ~ "cat tail",
                                  tree_species == "dead tree" ~ "dead plant",
                                  tree_species == "dead bush" ~ "dead plant",
                                  TRUE ~ tree_species))

## Add new variables ##

master <- master %>%
  mutate(site = case_when(str_detect(nest_id, "^a") == TRUE ~ "dwf",
                          str_detect(nest_id, "^j") == TRUE ~ "dwf",
                          str_detect(nest_id, "^e") == TRUE ~ "pcc",
                          str_detect(nest_id, "^r02$") == TRUE ~ "pcc",
                          str_detect(nest_id, "^r01") == TRUE ~ "dwf"),
         target_species = case_when(str_detect(host_species, "AMRO$|RWBL$|GRCA$|YEWA$") == TRUE ~ "yes",
                                    TRUE ~ "no")) %>%
  #alter for specific nests
  mutate(site = case_when(str_detect(nest_id, "^e6$") == TRUE ~ "dwf",
                          str_detect(nest_id, "^a1$") == TRUE ~ "pcc",
                          str_detect(nest_id, "^a2$") == TRUE ~ "pcc",
                          str_detect(nest_id, "^j2$") == TRUE ~ "pcc",
                          str_detect(nest_id, "^j6$") == TRUE ~ "outside study area",
                          str_detect(nest_id, "^e5$") == TRUE ~ "outside study area",
                          str_detect(nest_id, "^a4$") == TRUE ~ "outside study area",
                          TRUE ~ site))

non_target_master <- master %>%
  select(-master_number) %>%
  filter(target_species == "no")

not_monitored_master <- master %>%
  filter(nest_id %in% non_monitored_nests)

master <- master %>%
  filter(target_species == "yes") %>%
  filter(site %in% c("pcc", "dwf")) %>%
  filter(!nest_id %in% non_monitored_nests) %>%
  mutate(master_number = 1:n())

## Export data ##

write_csv(master,
          file = here("eggshells",
                      "field",
                      "processed_data",
                      "master_nest_list.csv"))

# Nest observations -------------------------------------------------------

not_monitored_nest_log <- nest_log %>%
  filter(nest_id %in% non_monitored_nests)

## Clean up data values ##

nest_log <- nest_log %>%
  filter(!nest_id %in% non_monitored_nests) %>%
  rename("field_notes" = "notes") %>%
  mutate(time = str_remove(string = time, 
                           pattern = "1899-12-31 "),
         adult = case_when(adult == "Yes" ~ "yes",
                           adult == "yes h" ~ "yes",
                           adult == "NA" ~ NA_character_,
                           adult == "n/a" ~ NA_character_,
                           adult == "?" ~ NA_character_,
                           adult == "unknown" ~ NA_character_,
                           TRUE ~ adult),
         stage = case_when(stage == "Laying/Incubation" ~ "laying/incubating",
                           stage == "Ready to Lay" ~ "ready to lay",
                           stage == "Do Not Monitor (note)" ~ "do not monitor (note)",
                           stage == "Incubating" ~ "incubating",
                           stage == "Laying" ~ "laying",
                           stage == "Building" ~ "building",
                           stage == "Hatched" ~ "hatching",
                           stage == "Unknown (note)" ~ "unknown (note)",
                           stage == "DNM" ~ "do not monitor (note)",
                           stage == "ready to Lay" ~ "ready to lay",
                           stage == "NA" ~ NA_character_,
                           stage == "hatched" ~ "hatching",
                           stage == "laying/incubation" ~ "laying/incubating",
                           stage == "REady to Lay" ~ "ready to lay",
                           stage == "do Not Monitor (note)" ~ "do not monitor (note)",
                           stage == "unknown" ~ "unknown (note)",
                           stage == "Unknown" ~ "unknown (note)",
                           stage == "REady to lay" ~ "ready to lay",
                           stage == "do not Monitor (note)" ~ "do not monitor (note)",
                           stage == "incubation" ~ "incubating",
                           stage == "Do not Monitor (note)" ~ "do not monitor (note)",
                           stage == "do not monitor" ~ "do not monitor (note)",
                           stage == "see notes" ~ "unknown (note)",
                           stage == "? see notes" ~ "unknown (note)",
                           stage == "incubating/hatched" ~ "incubating/hatching",
                           stage == "incubating?" ~ "incubating",
                           TRUE ~ stage),
         new_parasitism = case_when(new_parasitism == "n/a" ~ NA_character_,
                                    new_parasitism == "na" ~ NA_character_,
                                    new_parasitism == "u" ~ NA_character_,
                                    new_parasitism == "bo" ~ "no",
                                    new_parasitism == "n" ~ "no",
                                    new_parasitism == "NA" ~ NA_character_,
                                    new_parasitism == "unknown" ~ NA_character_,
                               TRUE ~ new_parasitism),
         bhco_swab = case_when(bhco_swab == "n/a" ~ NA_character_,
                               bhco_swab == "na" ~ NA_character_,
                               bhco_swab == "u" ~ NA_character_,
                               bhco_swab == "n" ~ "no",
                               bhco_swab == "NA" ~ NA_character_,
                               TRUE ~ bhco_swab),
         host_swab = case_when(host_swab == "n/a" ~ NA_character_,
                               host_swab == "na" ~ NA_character_,
                               host_swab == "u" ~ NA_character_,
                               host_swab == "n" ~ "no",
                               host_swab == "NA" ~ NA_character_,
                               host_swab == "Yes" ~ "yes",
                               TRUE ~ host_swab),
         predated = case_when(
           predated == "n/a" ~ NA_character_,
           predated == "na" ~ NA_character_,
           predated == "u" ~ NA_character_,
           predated == "n" ~ "no",
           predated == "NA" ~ NA_character_,
           predated == "Unknown" ~ NA_character_,
           predated == "unknown" ~ NA_character_,
           predated == "yes (see notes)" ~ "yes",
           predated == "yes(see notes)" ~ "yes",
                              TRUE ~ predated)) %>%
  separate(col = time,
           into = c("hr", "min", "sec"),
           sep = ":") %>%
  unite(col = time,
        c(hr, min),
        sep = ":") %>%
  mutate(time = str_replace(string = time,
                            pattern = "NA:NA",
                            replacement = NA_character_)) %>%
  select(-sec) %>%
  separate(col = date,
           into = c("month", "day"),
           sep = 1) %>%
  mutate(year = "2022") %>%
  unite(col = date,
        c(year, month, day),
        sep = "-") %>%
  mutate(date = str_replace(string = date,
                            pattern = "NA/NA",
                            replacement = NA_character_) %>%
           as.Date(.)) %>%
  drop_na(date)

#fill in empty days
x <- nest_log %>% pull(nest_id) %>% unique()
y <- list()
for (i in x) {
  temp <- nest_log[nest_log$nest_id == i,] %>% mutate(original = "original")
  dates <- data.frame(date = seq(min(temp$date), max(temp$date), by="days"))
  y[[i]]<- merge(x = dates,
                 y = temp,
                 all.x = TRUE) %>%
    mutate(original = str_replace_na(original),
           nest_id = {{i}}) %>%
    mutate(data_cleaning_notes = if_else(condition = str_detect(string = original,
                                                                pattern = "NA"),
                                         true = "nest not monitored on this date",
                                         false = NA_character_)) %>%
    select(-original)
}
nest_log <- reduce(.x = y,
                   .f = full_join)


## Add new variables ##

nest_log <- nest_log %>%
  left_join(.,
            master %>%
              select(nest_id,
                     target_species))

non_target_nest_log <- nest_log %>%
  filter(target_species == "no")

nest_log <- nest_log %>%
  filter(target_species == "yes")

nest_log <- nest_log %>%
  mutate(num_host_swabs = nest_log %>%
           filter(host_swab == "yes") %>%
           select(nest_id,
                  date) %>% 
           distinct() %>% 
           group_by(nest_id) %>% 
           mutate(number_swabs = n()) %>%
           ungroup() %>%
           select(nest_id,
                  number_swabs) %>% 
           distinct() %>%
           left_join(nest_log,
                     .) %>% 
           pull(number_swabs),
         day_of_year = yday(date)) %>% 
  mutate(num_host_swabs = case_when(is.na(num_host_swabs) == TRUE ~ 0 %>% as.double(),
                                    TRUE ~ num_host_swabs %>% as.double())) %>%
  left_join(.,
            nest_log %>%
              drop_na(num_host_egg) %>%
              group_by(nest_id) %>%
              summarise(max_host_egg = max(num_host_egg))) %>%
  mutate(max_host_egg = case_when(nest_id == "a49" ~ 5,
                                       nest_id == "e30" ~ 5,
                                       nest_id == "a65" ~ 5,
                                       nest_id == "a25" ~ 2,
                                       nest_id == "e100" ~ 7,
                                       TRUE ~ max_host_egg),
         data_cleaning_notes1 = case_when(nest_id == "a49" ~ "5 eggs were laid with one removed during laying",
                                          nest_id == "e30" ~ "5 eggs were laid with three lost during laying",
                                          nest_id == "a65" ~ "5 eggs were laid with one removed during laying",
                                          nest_id == "a25" ~ "2 eggs were laid with one predated during laying",
                                          nest_id == "e100" ~ "7 eggs were laid with one predated during laying",
                                          TRUE ~ NA_character_)) %>%
  unite(col = data_cleaning_notes,
        starts_with("data_cleaning_notes"),
        sep = ", ",
        na.rm = TRUE) %>%
  left_join(.,
            nest_log %>%
              drop_na(num_bhco_egg) %>%
              group_by(nest_id) %>%
              summarise(max_bhco_egg = max(num_bhco_egg))) %>%
  left_join(.,
            nest_log %>%
              drop_na(num_host_hatch) %>%
              group_by(nest_id) %>%
              summarise(max_host_hatch = max(num_host_hatch))) %>%
  mutate(max_host_hatch = case_when(nest_id == "a92" ~ 5,
                                       TRUE ~ max_host_hatch),
         data_cleaning_notes1 = case_when(nest_id == "a92" ~ "5 eggs hatched with 3 hatchlings predated",
                                          TRUE ~ NA_character_)) %>%
  unite(col = data_cleaning_notes,
        starts_with("data_cleaning_notes"),
        sep = ", ",
        na.rm = TRUE) %>%
  left_join(.,
            nest_log %>%
              drop_na(num_bhco_hatch) %>%
              group_by(nest_id) %>%
              summarise(max_bhco_hatch = max(num_bhco_hatch))) %>% 
  group_by(nest_id) %>%
  mutate(host_proportion_hatch = max_host_hatch/max_host_egg,
         bhco_proportion_hatch = max_bhco_hatch/max_bhco_egg) %>%
  mutate(host_proportion_hatch = case_when(max_host_egg == 0 ~ NA_character_,
                                           nest_id == "a11" ~ NA_character_,
                                           TRUE ~ host_proportion_hatch %>% as.character()),
         bhco_proportion_hatch = case_when(max_bhco_egg == 0 ~ NA_character_,
                                           TRUE ~ bhco_proportion_hatch %>% as.character()),
         data_cleaning_notes1 = case_when(nest_id == "a11" ~ "host proportion hatched unknown because nest found with 1 egg 1 hatch",
                                          max_host_egg == 0 ~ "host proportion hatch unknown because nest not initiated",
                                          TRUE ~ NA_character_),
         data_cleaning_notes2 = case_when(max_bhco_egg == 0 ~ "bhco proportion hatch unknown because nest not initiated",
                                          TRUE ~ NA_character_)) %>%
  mutate(host_proportion_hatch = host_proportion_hatch %>% as.numeric(),
         bhco_proportion_hatch = bhco_proportion_hatch %>% as.numeric()) %>%
  unite(col = data_cleaning_notes,
        starts_with("data_cleaning_notes"),
        sep = ", ",
        na.rm = TRUE) %>%
  mutate(first_monitoring_day = if_else(condition = date == min(date),
                             true = "yes",
                             false = "no"),
         last_monitoring_day = if_else(condition = date == max(date),
                           true = "yes",
                           false = "no")) %>%
  ungroup() %>%
  mutate(weather_impact = case_when(field_notes == "3 hatched, 2 eggs. Swab day 10 missed on rain day, no access." ~ "no",
                                    str_detect(string = field_notes,
                                               pattern = "(?i)wind(?!y)") == TRUE ~ "yes",
                                    str_detect(string = field_notes,
                                               pattern = "(?i)rain") == TRUE ~ "yes",
                                    str_detect(string = field_notes,
                                               pattern = "(?i)storm") == TRUE ~ "yes",
                                    TRUE ~ "no"),
         nest_damage = case_when(str_detect(string = field_notes,
                                            pattern = "(?i)no nest destruction") == TRUE ~ "no",
                                 str_detect(string = field_notes,
                                            pattern = "(?i)no nest damage") == TRUE ~"no",
                                 field_notes == "I don't know what happen but I wrote that the nest was damaged without eggs the 0612 and found a nest in the same GPS point with 2 eggs today. New nest?" ~ "no", 
                                 str_detect(string = field_notes,
                                            pattern = "(?i)damaged") == TRUE ~ "yes",
                                 str_detect(string = field_notes,
                                            pattern = "(?i)destroyed") == TRUE  ~ "yes", 
                                 str_detect(string = field_notes,
                                            pattern = "(?i)turned on side") == TRUE  ~ "yes", 
                                 str_detect(string = field_notes,
                                            pattern = "(?i)messed up") == TRUE  ~ "yes", 
                                 str_detect(string = field_notes,
                                            pattern = "(?i)pulled apart") == TRUE  ~ "yes", 
                                 str_detect(string = field_notes,
                                            pattern = "(?i)disheveled") == TRUE  ~ "yes", 
                                 str_detect(string = field_notes,
                                            pattern = "(?i)detsroyed") == TRUE  ~ "yes", 
                                 str_detect(string = field_notes,
                                            pattern = "(?i)torn off") == TRUE  ~ "yes", 
                                 str_detect(string = field_notes,
                                            pattern = "(?i)pulled up") == TRUE  ~ "yes", 
                                 str_detect(string = field_notes,
                                            pattern = "(?i)nest tipped over") == TRUE ~ "yes",
                                 TRUE ~ "no"),
         day_before = date - 1,
         day_after = date + 1)

nest_log <- nest_log %>%
  left_join(.,
            nest_log %>%
              mutate(nest_deserted = case_when(str_detect(string = field_notes,
                                                     pattern = "likely deserted") == TRUE ~ "yes",
                                          str_detect(string = field_notes,
                                                     pattern = "nest deserted") == TRUE ~ "yes",
                                          TRUE ~ NA_character_)) %>%
              filter(nest_deserted == "yes") %>%
              select(nest_id,
                     nest_deserted)) %>%
  mutate(nest_deserted = case_when(nest_deserted == "yes" ~ "yes",
                                   is.na(nest_deserted) == TRUE ~ "no"))

nest_log <- nest_log %>%
  left_join(.,
            nest_log %>%
              mutate(nest_buried = case_when(str_detect(string = field_notes,
                                                     pattern = "nest buried") == TRUE ~ "yes",
                                          TRUE ~ NA_character_)) %>%
              filter(nest_buried == "yes") %>%
              select(nest_id,
                     nest_buried)) %>%
  mutate(nest_buried = case_when(nest_buried == "yes" ~ "yes",
                                 is.na(nest_buried) == TRUE ~ "no"))

#add variables for number eggs/ hatchlings the day before and day after
x <- nest_log %>% 
  group_by(nest_id) %>%
  summarise(number_observations = n()) %>%
  filter(number_observations > 1) %>%
  pull(nest_id) %>%
  unique()
l1 <- list()
for(i in x){
  temp <- nest_log %>% filter(nest_id == {{i}})
  y <- temp %>% filter(first_monitoring_day == "no") %>% 
    drop_na(date) %>% pull(date) %>% as.character()
  l2 <- list()
  for (j in y) {
    db <- temp %>% filter(date == {{j}}) %>% pull(day_before) %>% str_remove(string = ., pattern = " UTC")
    a <- temp %>% filter(date == db) %>% pull(num_host_egg)
    b <- temp %>% filter(date == db) %>% pull(num_bhco_egg)
    c <- temp %>% filter(date == db) %>% pull(num_host_hatch)
    d <- temp %>% filter(date == db) %>% pull(num_bhco_hatch)
    m <- temp %>% filter(date == db) %>% pull(data_cleaning_notes)
    l2[[j]] <- data.frame(nest_id = i,
                          date = j,
                          day_before = db,
                          num_host_egg_day_before = a,
                          num_bhco_egg_day_before = b,
                          num_host_hatch_day_before = c,
                          num_bhco_hatch_day_before = d,
                          monitored_day_before = m)
  }
  l1[[i]] <- reduce(.x = l2,
                    .f = full_join)
}
nest_log <- full_join(nest_log,
                  reduce(.x = l1,
                         .f = full_join) %>%
                    mutate(date = as.Date(date),
                           day_before = as.Date(day_before))) %>%
  mutate(monitored_day_before = case_when(first_monitoring_day == "yes" ~ "no",
                                          is.na(num_host_egg_day_before) == TRUE ~ "no",
                                          monitored_day_before == "nest not monitored on this date" ~ "no",
                                          is.na(monitored_day_before) == TRUE ~ "yes"))


x <- nest_log %>% 
  group_by(nest_id) %>%
  summarise(number_observations = n()) %>%
  filter(number_observations > 1) %>%
  pull(nest_id) %>%
  unique()
l1 <- list()
for(i in x){
  temp <- nest_log %>% filter(nest_id == {{i}})
  y <- temp %>% filter(last_monitoring_day == "no") %>% 
    drop_na(date) %>% pull(date) %>% as.character()
  l2 <- list()
  for (j in y) {
    da <- temp %>% filter(date == {{j}}) %>% pull(day_after) %>% str_remove(string = ., pattern = " UTC")
    e <- temp %>% filter(date == da) %>% pull(num_host_egg)
    f <- temp %>% filter(date == da) %>% pull(num_bhco_egg)
    g <- temp %>% filter(date == da) %>% pull(num_host_hatch)
    h <- temp %>% filter(date == da) %>% pull(num_bhco_hatch)
    l2[[j]] <- data.frame(nest_id = i,
                          date = j,
                          day_after = da,
                          num_host_egg_day_after = e,
                          num_bhco_egg_day_after = f,
                          num_host_hatch_day_after = g,
                          num_bhco_hatch_day_after = h)
  }
  l1[[i]] <- reduce(.x = l2,
                    .f = full_join)
}
nest_log <- full_join(nest_log,
                      reduce(.x = l1,
                             .f = full_join) %>%
                        mutate(date = as.Date(date),
                               day_after = as.Date(day_after)))

nest_log <- nest_log %>%
  mutate(day_of_year_host_egg_found = case_when(num_host_egg >= 1 ~ day_of_year %>%
                                             as.integer(),
                                           TRUE ~ NA_integer_),
         day_of_year_bhco_egg_found = case_when(num_bhco_egg >= 1 ~ day_of_year %>%
                                                  as.integer(),
                                                TRUE ~ NA_integer_)) %>%
  group_by(nest_id) %>%
  mutate(day_of_year_first_host_egg_found = case_when(max_host_egg == 0 ~ NA_integer_,
                                                 max_host_egg >= 1 ~ min(day_of_year_host_egg_found,
                                                                              na.rm = TRUE) %>%
                                                   as.integer()),
         day_of_year_first_bhco_egg_found = case_when(max_bhco_egg == 0 ~ NA_integer_,
                                                      max_bhco_egg >= 1 ~ min(day_of_year_bhco_egg_found,
                                                                                   na.rm = TRUE) %>%
                                                        as.integer())) %>%
  ungroup() %>%
  mutate(first_host_egg_found = case_when(day_of_year_host_egg_found == day_of_year_first_host_egg_found ~ "yes",
                                     TRUE ~ "no"),
         first_bhco_egg_found = case_when(day_of_year_bhco_egg_found == day_of_year_first_bhco_egg_found ~ "yes",
                                          TRUE ~ "no")) %>%
  select(-day_of_year_host_egg_found,
         -day_of_year_bhco_egg_found) 

nest_log <- nest_log %>%
  mutate(first_host_egg_witnessed = if_else(condition = str_detect(string = nest_id,
                                                              pattern = c(nest_log %>%
                                                                            filter(first_host_egg_found == "yes" & num_host_egg == 1 & num_host_egg_day_before == 0) %>%
                                                                            pull(nest_id) %>%
                                                                            unique() %>%
                                                                            paste(., collapse = "$|"))),
                                       true = "yes",
                                       false = "no"),
         first_bhco_egg_witnessed = if_else(condition = str_detect(string = nest_id,
                                                                   pattern = c(nest_log %>%
                                                                                 filter(first_bhco_egg_found == "yes" & num_bhco_egg == 1 & num_bhco_egg_day_before == 0) %>%
                                                                                 pull(nest_id) %>%
                                                                                 unique() %>%
                                                                                 paste(., collapse = "$|"))),
                                            true = "yes",
                                            false = "no")) %>%
  mutate(
    host_nest_outcome = case_when(max_host_egg == 0 ~ NA_character_,
                                  nest_id %in% c(
                                    "e128",  #nest that was not passed to roxanne
                                    #nests that we lost access to             
                                    "j42", 
                                    "j44",
                                    "j58") ~ NA_character_,
                                  #host hatch on final day == 0 - nest was predated and is total failure
                                  host_proportion_hatch == 1 ~ "complete sucess",
                                  host_proportion_hatch == 0 ~ "total failure",
                                  host_proportion_hatch < 1 ~ "partial sucess",
                                  TRUE ~ NA_character_),
    bhco_nest_outcome = case_when(max_bhco_egg == 0 ~ NA_character_,
                                  bhco_proportion_hatch == 1 ~ "complete sucess",
                                  bhco_proportion_hatch == 0 ~ "total failure",
                                  bhco_proportion_hatch < 1 ~ "partial sucess",
                                  TRUE ~ NA_character_)) %>%
  mutate(host_eggs_increasing = if_else(condition = num_host_egg < num_host_egg_day_after,
                                        true = "yes",
                                        false = "no"),
         host_eggs_decreasing = if_else(condition = num_host_egg > num_host_egg_day_after,
                                        true = "yes",
                                        false = "no"),
         bhco_eggs_increasing = if_else(condition = num_bhco_egg < num_bhco_egg_day_after,
                                        true = "yes",
                                        false = "no"),
         bhco_eggs_decreasing = if_else(condition = num_bhco_egg > num_bhco_egg_day_after,
                                        true = "yes",
                                        false = "no"))

nest_log <- nest_log %>%
  mutate(addled_host_egg = if_else(condition = str_detect(string = nest_id,
                                                          pattern = nest_log %>% 
                                                            mutate(addled_egg = case_when(str_detect(string = field_notes,
                                                                                                     pattern = "(?i)dud") == TRUE ~ "yes",
                                                                                          str_detect(string = field_notes,
                                                                                                     pattern = "(?i)addeled") == TRUE ~ "yes", 
                                                                                          str_detect(string = field_notes,
                                                                                                     pattern = "(?i)addled") == TRUE ~ "yes",
                                                                                          TRUE ~ "no")) %>% 
                                                            filter(addled_egg == "yes") %>% 
                                                            pull(nest_id) %>%
                                                            unique() %>%
                                                            paste(., collapse = "$|")),
                                   true = "yes",
                                   false = "no"))

nest_log <- nest_log %>%
  left_join(.,
            nest_log %>%
              mutate(num_addled_host_egg = case_when(addled_host_egg == "yes" & last_monitoring_day == "yes" ~ num_host_egg %>% as.integer(),
                                                     addled_host_egg == "no" ~ 0 %>% as.integer(),
                                                     TRUE ~ NA_integer_)) %>%
              select(nest_id,
                     num_addled_host_egg) %>%
              drop_na(num_addled_host_egg) %>%
              distinct())

nest_log <- nest_log %>% 
  mutate(day_of_year_host_hatch_found = case_when(num_host_hatch >= 1 ~ day_of_year %>%
                                                    as.integer(),
                                                  TRUE ~ NA_integer_),
         day_of_year_bhco_hatch_found = case_when(num_bhco_hatch >= 1 ~ day_of_year %>%
                                                    as.integer(),
                                                  TRUE ~ NA_integer_)) %>%
  group_by(nest_id) %>%
  mutate(host_first_hatch_day_of_year = case_when(max_host_hatch == 0 ~ NA_integer_,
                                                        max_host_hatch >= 1 ~ min(day_of_year_host_hatch_found,
                                                                                  na.rm = TRUE) %>%
                                                          as.integer()),
         bhco_first_hatch_day_of_year = case_when(max_bhco_hatch == 0 ~ NA_integer_,
                                                        max_bhco_hatch >= 1 ~ min(day_of_year_bhco_hatch_found,
                                                                                  na.rm = TRUE) %>%
                                                          as.integer())) %>%
  ungroup() %>%
  mutate(first_host_hatch_found = case_when(day_of_year_host_hatch_found == host_first_hatch_day_of_year ~ "yes",
                                            TRUE ~ "no"),
         first_bhco_hatch_found = case_when(day_of_year_bhco_hatch_found == bhco_first_hatch_day_of_year ~ "yes",
                                            TRUE ~ "no")) %>%
  select(-day_of_year_host_hatch_found,
         -day_of_year_bhco_hatch_found)

nest_log <- nest_log %>%
  mutate(first_host_hatch_witnessed = if_else(condition = str_detect(string = nest_id,
                                                                     pattern = c(nest_log %>%
                                                                                   filter(first_host_hatch_found == "yes" & num_host_hatch >= 1 & num_host_hatch_day_before == 0) %>%
                                                                                   pull(nest_id) %>%
                                                                                   unique() %>%
                                                                                   paste(., collapse = "$|"))),
                                              true = "yes",
                                              false = "no"),
         first_bhco_hatch_witnessed = if_else(condition = str_detect(string = nest_id,
                                                                     pattern = c(nest_log %>%
                                                                                   filter(first_bhco_hatch_found == "yes" & num_bhco_hatch >= 1 & num_bhco_hatch_day_before == 0) %>%
                                                                                   pull(nest_id) %>%
                                                                                   unique() %>%
                                                                                   paste(., collapse = "$|"))),
                                              true = "yes",
                                              false = "no")) %>%
  mutate(host_first_hatch_day_of_year = case_when(first_host_hatch_witnessed == "no" ~ NA_integer_,
                                                        first_host_hatch_witnessed == "yes" ~ host_first_hatch_day_of_year),
         bhco_first_hatch_day_of_year = case_when(first_bhco_hatch_witnessed == "no" ~ NA_integer_,
                                                        first_bhco_hatch_witnessed == "yes" ~ bhco_first_hatch_day_of_year)) %>%
  mutate(data_cleaning_notes1 = case_when(max_host_egg == 0 ~ NA_character_, #note from another part of code
                                          max_host_hatch == 0 ~ "host did not hatch",
                                          first_host_hatch_witnessed == "no" ~ "first host hatch was not witnessed"),
         data_cleaning_notes2 = case_when(max_bhco_egg == 0 ~ NA_character_, #note from another part of code,
                                          max_bhco_hatch == 0 ~ "bhco did not hatch",
                                          first_bhco_hatch_witnessed == "no" ~ "first bhco hatch was not witnessed")) %>%
  unite(col = "data_cleaning_notes",
        starts_with("data_cleaning_notes"),
        sep = ", ",
        na.rm = TRUE) %>%
  mutate(host_first_hatch_date = case_when(is.na(host_first_hatch_day_of_year) == FALSE ~ 
                                                   as.Date(host_first_hatch_day_of_year, origin = "2021-12-31")),
         bhco_first_hatch_date = case_when(is.na(bhco_first_hatch_day_of_year) == FALSE ~ 
                                             as.Date(bhco_first_hatch_day_of_year, origin = "2021-12-31"))) 
nest_log <- nest_log %>%
mutate(host_swab_first_egg = if_else(condition = str_detect(string = nest_id,
                                                     pattern = c(nest_log %>% 
                                                                   filter(host_swab == "yes" & num_host_egg == 1 & first_host_hatch_witnessed == "yes") %>%
                                                                   pull(nest_id) %>% 
                                                                   unique() %>% 
                                                                   paste(., collapse = "$|"))),
                              true = "yes",
                              false = "no"))

nest_log <- nest_log %>%
  left_join(.,
            nest_log %>%
              filter(first_host_hatch_witnessed == "yes" &
                       first_host_hatch_found == "yes") %>%
              select(nest_id, num_host_egg_day_before) %>%
              rename("num_host_egg_before_hatch" = "num_host_egg_day_before")) %>%
  left_join(.,
            nest_log %>%
              filter(first_bhco_hatch_witnessed == "yes" &
                       first_bhco_hatch_found == "yes") %>%
              select(nest_id, num_bhco_egg_day_before) %>%
              rename("num_bhco_egg_before_hatch" = "num_bhco_egg_day_before")) %>%
  mutate(host_hatching_sucess = max_host_hatch/num_host_egg_before_hatch,
         bhco_hatching_sucess = max_bhco_hatch/num_bhco_egg_before_hatch)

nest_log <- nest_log %>%
  left_join(.,
            master %>% 
              select(nest_id,
                     host_species) %>%
              distinct())

nest_log <- nest_log %>%
  mutate(average_host_incubation_period = case_when(host_species == "AMRO" ~ 12,
                                                    host_species == "GRCA" ~ 13,
                                                    host_species == "YEWA" ~ 11,
                                                    host_species == "RWBL" ~ 12 ),
         citation_average_host_incubation_period = case_when(host_species == "AMRO" ~ "Ospina, E. A., Merrill, L., &amp; Benson, T. J. (2018). Incubation temperature impacts nestling growth and survival in an open-cup nesting passerine. Ecology and Evolution, 8(6), 3270–3279. https://doi.org/10.1002/ece3.3911",
                                                        host_species == "GRCA" ~ "Johnson, E. J. and L. B. Best. (1980). Breeding biology of the Gray Catbird in Iowa. Iowa State J. Res. 55:171-183.",
                                                        host_species == "YEWA" ~ "McMaster, D. G., & Sealy, S. G. (1998). Short incubation periods of Brown-headed Cowbirds: How do cowbird eggs hatch before Yellow Warbler eggs?. The Condor, 100(1), 102-111.",
                                                        host_species == "RWBL" ~ "Birds of the world: Yasukawa, K. and W. A. Searcy (2020). Red-winged Blackbird (Agelaius phoeniceus), version 1.0. In Birds of the World (P. G. Rodewald, Editor). Cornell Lab of Ornithology, Ithaca, NY, USA. https://doi-org.proxy1.lib.uwo.ca/10.2173/bow.rewbla.01"),
         average_bhco_incubation_period = 12,
         citation_average_bhco_incubation_period = "McMaster, D. G., & Sealy, S. G. (1998). Short incubation periods of Brown-headed Cowbirds: How do cowbird eggs hatch before Yellow Warbler eggs?. The Condor, 100(1), 102-111.")

nest_log <- nest_log %>%
  left_join(.,
            nest_log %>% 
              mutate(across(.cols = c(day_of_year_first_host_egg_found, num_host_egg, average_host_incubation_period),
                            .fns = as.integer)) %>%
              mutate(host_clutch_initiation_day_of_year = case_when(first_host_egg_witnessed == "yes" ~ 
                                                                      day_of_year_first_host_egg_found %>% as.integer(),
                                                                    #nest found with eggs being laid
                                                                    first_host_egg_found == "yes" & 
                                                                      first_host_egg_witnessed == "no" & 
                                                                      host_eggs_increasing == "yes" ~ 
                                                                      day_of_year_first_host_egg_found - (num_host_egg - 1) %>% as.integer(),
                                                                    #nest found with eggs already laid
                                                                    first_host_egg_found == "yes" & 
                                                                      first_host_egg_witnessed == "no" & 
                                                                      host_eggs_increasing == "no" &
                                                                      first_monitoring_day == "yes"  &
                                                                      first_host_hatch_witnessed == "yes" ~ 
                                                                      host_first_hatch_day_of_year - average_host_incubation_period - (num_host_egg - 1) %>% as.integer(),
                                                                    first_host_egg_found == "yes" & 
                                                                      first_host_egg_witnessed == "no" & 
                                                                      is.na(host_eggs_increasing) == TRUE &
                                                                      first_host_hatch_witnessed == "yes" ~ 
                                                                      host_first_hatch_day_of_year - average_host_incubation_period - (num_host_egg - 1) %>% as.integer(),
                                                                    TRUE ~ NA_integer_)) %>%
              mutate(host_clutch_initiation_estimated = case_when(first_host_egg_witnessed == "yes" ~ "no",
                                                             first_host_egg_found == "yes" & 
                                                               first_host_egg_witnessed == "no" & 
                                                               host_eggs_increasing == "yes" ~ "yes",
                                                             first_host_egg_found == "yes" & 
                                                               first_host_egg_witnessed == "no" & 
                                                               host_eggs_increasing == "no" &
                                                               first_monitoring_day == "yes"  &
                                                               first_host_hatch_witnessed == "yes" ~ "yes",
                                                             first_host_egg_found == "yes" & 
                                                               first_host_egg_witnessed == "no" & 
                                                               is.na(host_eggs_increasing) == TRUE &
                                                               first_monitoring_day == "yes"  &
                                                               first_host_hatch_witnessed == "yes" ~ "yes"),
                     data_cleaning_notes_1 = case_when(first_host_egg_found == "yes" & 
                                                         first_host_egg_witnessed == "no" & 
                                                         host_eggs_increasing == "yes" ~ "clutch initiation estimated based on hosts laying every 24hrs",
                                                       first_host_egg_found == "yes" & 
                                                         first_host_egg_witnessed == "no" & 
                                                         host_eggs_increasing == "no" &
                                                         first_monitoring_day == "yes"  &
                                                         first_host_hatch_witnessed == "yes" ~ "clutch initiation estimated based on average incubation time and hosts laying every 24hrs",
                                                       first_host_egg_found == "yes" & 
                                                         first_host_egg_witnessed == "no" & 
                                                         is.na(host_eggs_increasing) == TRUE &
                                                         first_monitoring_day == "yes"  &
                                                         first_host_hatch_witnessed == "yes" ~ "clutch initiation estimated based on average incubation time and hosts laying every 24hrs")
              ) %>%
              select(nest_id,
                     host_clutch_initiation_day_of_year,
                     host_clutch_initiation_estimated,
                     data_cleaning_notes_1) %>%
              distinct() %>%
              drop_na(host_clutch_initiation_day_of_year),
            by = "nest_id") %>%
  left_join(.,
            nest_log %>% 
              mutate(across(.cols = c(day_of_year_first_bhco_egg_found, num_bhco_egg, average_bhco_incubation_period),
                            .fns = as.integer)) %>%
              mutate(bhco_clutch_initiation_day_of_year = case_when(first_bhco_egg_witnessed == "yes" ~ 
                                                                      day_of_year_first_bhco_egg_found %>% as.integer(),
                                                                    #nest found with eggs being laid
                                                                    first_bhco_egg_found == "yes" & 
                                                                      first_bhco_egg_witnessed == "no" & 
                                                                      bhco_eggs_increasing == "yes" ~ 
                                                                      day_of_year_first_bhco_egg_found - (num_bhco_egg - 1) %>% as.integer(),
                                                                    #nest found with eggs already laid
                                                                    first_bhco_egg_found == "yes" & 
                                                                      first_bhco_egg_witnessed == "no" & 
                                                                      bhco_eggs_increasing == "no" &
                                                                      first_monitoring_day == "yes"  &
                                                                      first_bhco_hatch_witnessed == "yes" ~ 
                                                                      bhco_first_hatch_day_of_year - average_bhco_incubation_period - (num_bhco_egg - 1) %>% as.integer(),
                                                                    first_bhco_egg_found == "yes" & 
                                                                      first_bhco_egg_witnessed == "no" & 
                                                                      is.na(bhco_eggs_increasing) == TRUE &
                                                                      first_bhco_hatch_witnessed == "yes" ~ 
                                                                      bhco_first_hatch_day_of_year - average_bhco_incubation_period - (num_bhco_egg - 1) %>% as.integer(),
                                                                    TRUE ~ NA_integer_)) %>%
              mutate(bhco_clutch_initiation_estimated = case_when(first_bhco_egg_witnessed == "yes" ~ "no",
                                                             first_bhco_egg_found == "yes" & 
                                                               first_bhco_egg_witnessed == "no" & 
                                                               bhco_eggs_increasing == "yes" ~ "yes",
                                                             first_bhco_egg_found == "yes" & 
                                                               first_bhco_egg_witnessed == "no" & 
                                                               bhco_eggs_increasing == "no" &
                                                               first_monitoring_day == "yes"  &
                                                               first_bhco_hatch_witnessed == "yes" ~ "yes",
                                                             first_bhco_egg_found == "yes" & 
                                                               first_bhco_egg_witnessed == "no" & 
                                                               is.na(bhco_eggs_increasing) == TRUE &
                                                               first_monitoring_day == "yes"  &
                                                               first_bhco_hatch_witnessed == "yes" ~ "yes"),
                     data_cleaning_notes_2 = case_when(first_bhco_egg_found == "yes" & 
                                                         first_bhco_egg_witnessed == "no" & 
                                                         bhco_eggs_increasing == "yes" ~ "clutch initiation estimated based on bhcos laying every 24hrs",
                                                       first_bhco_egg_found == "yes" & 
                                                         first_bhco_egg_witnessed == "no" & 
                                                         bhco_eggs_increasing == "no" &
                                                         first_monitoring_day == "yes"  &
                                                         first_bhco_hatch_witnessed == "yes" ~ "clutch initiation estimated based on average incubation time and bhcos laying every 24hrs",
                                                       first_bhco_egg_found == "yes" & 
                                                         first_bhco_egg_witnessed == "no" & 
                                                         is.na(bhco_eggs_increasing) == TRUE &
                                                         first_monitoring_day == "yes"  &
                                                         first_bhco_hatch_witnessed == "yes" ~ "clutch initiation estimated based on average incubation time and bhcos laying every 24hrs")
              ) %>%
              select(nest_id,
                     bhco_clutch_initiation_day_of_year,
                     bhco_clutch_initiation_estimated,
                     data_cleaning_notes_2) %>%
              distinct() %>%
              drop_na(bhco_clutch_initiation_day_of_year),
            by = "nest_id")

nest_log <- nest_log %>% 
  left_join(.,
            nest_log %>%
              mutate(data_cleaning_notes_3 = case_when(is.na(host_clutch_initiation_day_of_year) == TRUE & 
                                                         max_host_egg == 0 ~ "host clutch not initiated",
                                                       is.na(host_clutch_initiation_day_of_year) == TRUE & 
                                                         first_host_egg_found == "yes" & 
                                                         first_host_egg_witnessed == "no" & 
                                                         is.na(host_first_hatch_date) == TRUE &
                                                         first_host_hatch_witnessed == "no" ~ "unable to estimate clutch inititaion because first host hatch not witnessed",
                                                       is.na(host_clutch_initiation_day_of_year) == TRUE & 
                                                         first_host_egg_found == "yes" & 
                                                         first_host_egg_witnessed == "no" & 
                                                         max_host_hatch == 0 ~ "unable to estimate clutch inititaion because host did not hatch")) %>%
              select(nest_id,
                     data_cleaning_notes_2) %>%
              distinct() %>%
              drop_na(data_cleaning_notes_2)) %>%
  unite(col = "data_cleaning_notes",
        starts_with("data_cleaning_notes"),
        sep = ", ",
        na.rm = TRUE) %>%
  mutate(host_clutch_initiation_date = case_when(is.na(host_clutch_initiation_day_of_year) == FALSE ~ 
                                                   as.Date(host_clutch_initiation_day_of_year, origin = "2021-12-31"))) %>%
  left_join(.,
            nest_log %>%
              mutate(data_cleaning_notes_2 = case_when(is.na(bhco_clutch_initiation_day_of_year) == TRUE & 
                                                         max_bhco_egg == 0 ~ "bhco clutch not initiated",
                                                       is.na(bhco_clutch_initiation_day_of_year) == TRUE & 
                                                         first_bhco_egg_found == "yes" & 
                                                         first_bhco_egg_witnessed == "no" & 
                                                         is.na(bhco_first_hatch_date) == TRUE &
                                                         first_bhco_hatch_witnessed == "no" ~ "unable to estimate clutch inititaion because first bhco hatch not witnessed",
                                                       is.na(bhco_clutch_initiation_day_of_year) == TRUE & 
                                                         first_bhco_egg_found == "yes" & 
                                                         first_bhco_egg_witnessed == "no" & 
                                                         max_bhco_hatch == 0 ~ "unable to estimate clutch inititaion because bhco did not hatch")) %>%
              select(nest_id,
                     data_cleaning_notes_2) %>%
              distinct() %>%
              drop_na(data_cleaning_notes_2)) %>%
  unite(col = "data_cleaning_notes",
        starts_with("data_cleaning_notes"),
        sep = ", ",
        na.rm = TRUE) %>%
  mutate(bhco_clutch_initiation_date = case_when(is.na(bhco_clutch_initiation_day_of_year) == FALSE ~ 
                                                   as.Date(bhco_clutch_initiation_day_of_year, origin = "2021-12-31")))

nest_log <- nest_log %>%
  mutate(host_last_hatch_day_of_year_temp = case_when(num_host_hatch >= 1 & 
                                                       num_host_egg == 0 & 
                                                       num_host_egg_day_before > 0 ~ day_of_year %>% as.integer(),
                                                     num_host_hatch >= 1 & 
                                                       addled_host_egg == "yes" &
                                                       num_host_egg_day_before > num_addled_host_egg &
                                                       num_host_egg == num_addled_host_egg & 
                                                       num_host_egg_day_after == num_addled_host_egg ~ day_of_year %>% as.integer(),
                                                     TRUE ~ NA_integer_),
         bhco_last_hatch_day_of_year_temp = case_when(num_bhco_hatch >= 1 & 
                                                       num_bhco_egg == 0 & 
                                                       num_bhco_egg_day_before > 0 ~ day_of_year %>% as.integer(),
                                                     TRUE ~ NA_integer_)) %>%
  group_by(nest_id) %>%
  mutate(host_last_hatch_day_of_year = case_when(max_host_hatch == 0 ~ NA_integer_,
                                                max_host_hatch >= 1 ~ min(host_last_hatch_day_of_year_temp,
                                                                          na.rm = TRUE) %>%
                                                  as.integer()),
         bhco_last_hatch_day_of_year = case_when(max_bhco_hatch == 0 ~ NA_integer_,
                                                     max_bhco_hatch >= 1 ~ min(bhco_last_hatch_day_of_year_temp,
                                                                               na.rm = TRUE) %>%
                                                       as.integer())) %>%
  ungroup() %>%
  mutate(last_host_hatch_found = case_when(host_last_hatch_day_of_year == host_last_hatch_day_of_year_temp ~ "yes",
                                           TRUE ~ "no"),
         last_bhco_hatch_found = case_when(bhco_last_hatch_day_of_year == bhco_last_hatch_day_of_year_temp ~ "yes",
                                           TRUE ~ "no")) %>%
  select(-host_last_hatch_day_of_year_temp,
         -bhco_last_hatch_day_of_year_temp) 

nest_log <- nest_log %>%
  mutate(last_host_hatch_witnessed = if_else(condition = str_detect(string = nest_id,
                                                                    pattern = c(nest_log %>%
                                                                                  filter((last_host_hatch_found == "yes" & 
                                                                                            num_host_egg == 0 & 
                                                                                            num_host_egg_day_before > 0) |
                                                                                           (last_host_hatch_found == "yes" & 
                                                                                              addled_host_egg == "yes" &
                                                                                              num_host_egg_day_before > num_addled_host_egg &
                                                                                              num_host_egg == num_addled_host_egg & 
                                                                                              num_host_egg_day_after == num_addled_host_egg )) %>%
                                                                                  pull(nest_id) %>%
                                                                                  unique() %>%
                                                                                  paste(., collapse = "$|"))),
                                             true = "yes",
                                             false = "no"),
         last_bhco_hatch_witnessed = if_else(condition = str_detect(string = nest_id,
                                                                    pattern = c(nest_log %>%
                                                                                  filter((last_bhco_hatch_found == "yes" & 
                                                                                            num_bhco_egg == 0 & 
                                                                                            num_bhco_egg_day_before > 0)) %>%
                                                                                  pull(nest_id) %>%
                                                                                  unique() %>%
                                                                                  paste(., collapse = "$|"))),
                                             true = "yes",
                                             false = "no")) %>%
  mutate(host_last_hatch_day_of_year = case_when(last_host_hatch_witnessed == "no" ~ NA_integer_,
                                                last_host_hatch_witnessed == "yes" ~ host_last_hatch_day_of_year),
         bhco_last_hatch_day_of_year = case_when(last_bhco_hatch_witnessed == "no" ~ NA_integer_,
                                                last_bhco_hatch_witnessed == "yes" ~ bhco_last_hatch_day_of_year)) %>%
  mutate(data_cleaning_notes1 = case_when(max_host_hatch == 0 ~ NA_character_, #already note from above
                                          max_host_egg == 0 ~ NA_character_, #already note from above
                                          last_host_hatch_witnessed == "no" ~ "last host hatch was not witnessed"),
         data_cleaning_notes2 = case_when(max_bhco_hatch == 0 ~ NA_character_, #already note from above
                                          max_bhco_egg == 0 ~ NA_character_, #already note from above
                                          last_bhco_hatch_witnessed == "no" ~ "last bhco hatch was not witnessed")) %>%
  unite(col = "data_cleaning_notes",
        starts_with("data_cleaning_notes"),
        sep = ", ",
        na.rm = TRUE) %>%
  mutate(host_last_hatch_date = case_when(is.na(host_last_hatch_day_of_year) == FALSE ~ 
                                                as.Date(host_last_hatch_day_of_year, origin = "2021-12-31")),
         bhco_last_hatch_date = case_when(is.na(bhco_last_hatch_day_of_year) == FALSE ~ 
                                                as.Date(bhco_last_hatch_day_of_year, origin = "2021-12-31")))

nest_log <- nest_log %>%
  mutate(predation_weather_effect = if_else(condition = str_detect(string = nest_id,
                                                                   pattern = nest_log %>% 
                                                                     filter(predated == "yes" | weather_impact == "yes" ) %>%
                                                                     pull(nest_id)  %>%
                                                                     unique() %>%
                                                                     paste(., collapse = "$|")),
                                            true = "yes",
                                            false = "no"),
         desertion_burial_effect = if_else(condition = str_detect(string = nest_id,
                                                                  pattern = nest_log %>%
                                                                    filter(nest_deserted == "yes" | nest_buried == "yes" ) %>%
                                                                    pull(nest_id) %>%
                                                                    unique() %>%
                                                                    paste(., collapse = "$|")),
                                           true = "yes",
                                           false = "no"),
         natural_parasitism = if_else(condition = str_detect(string = nest_id,
                                                             pattern = nest_log %>% 
                                                               filter(num_bhco_egg > 0) %>%
                                                               pull(nest_id) %>%
                                                               unique() %>%
                                                               paste(., collapse = "$|")),
                                      true = "yes",
                                      false = "no"))

nest_log <- nest_log %>%
  mutate(host_time_to_hatch = host_last_hatch_day_of_year - host_clutch_initiation_day_of_year,
         host_hatching_period = host_last_hatch_day_of_year - host_first_hatch_day_of_year,
         bhco_time_to_hatch = bhco_last_hatch_day_of_year - bhco_clutch_initiation_day_of_year,
         bhco_hatching_period = bhco_last_hatch_day_of_year - bhco_first_hatch_day_of_year)

nest_log <- nest_log %>%
  select(nest_id,
         date,
         day_of_year,
         time,
         adult,
         num_host_egg,
         num_host_hatch,
         num_bhco_egg,
         num_bhco_hatch,
         new_parasitism,
         predated,
         nest_damage,
         weather_impact,
         natural_parasitism,
         nest_deserted,
         nest_buried,
         desertion_burial_effect,
         predation_weather_effect,
         host_nest_outcome,
         bhco_nest_outcome,
         max_host_egg,
         addled_host_egg,
         num_addled_host_egg,
         max_host_hatch,
         num_host_egg_before_hatch,
         host_proportion_hatch,
         host_hatching_sucess,
         host_clutch_initiation_date,
         host_clutch_initiation_day_of_year,
         host_clutch_initiation_estimated,
         host_first_hatch_date,
         host_first_hatch_day_of_year,
         host_last_hatch_date,
         host_last_hatch_day_of_year,
         host_time_to_hatch,
         host_hatching_period,
         max_bhco_egg,
         max_bhco_hatch,
         num_bhco_egg_before_hatch,
         bhco_proportion_hatch,
         bhco_hatching_sucess,
         bhco_clutch_initiation_date,
         bhco_clutch_initiation_day_of_year,
         bhco_clutch_initiation_estimated,
         bhco_first_hatch_date,
         bhco_first_hatch_day_of_year,
         bhco_last_hatch_date,
         bhco_last_hatch_day_of_year,
         bhco_time_to_hatch,
         bhco_hatching_period,
         bhco_swab,
         host_swab,
         host_swab_first_egg,
         num_host_swabs,
         average_host_incubation_period,
         citation_average_host_incubation_period,
         average_bhco_incubation_period,
         citation_average_bhco_incubation_period,
         field_notes,
         data_cleaning_notes)

## Export data ##

write_csv(nest_log,
          here("eggshells",
               "field",
               "processed_data",
          "nest_observations.csv"))

# Master and nest log -------------------------------------------------------------------------

master_and_nest <- full_join(nest_log,
                                master)
## Export data ##
write_csv(master_and_nest,
          here("eggshells",
               "field",
               "processed_data",
               "master_and_nest_observations.csv"))

# Non-monitored nests -------------------------------------------------------------------------

not_monitored_nests <- full_join(not_monitored_master,
                                not_monitored_nest_log)
## Export data ##
write_csv(not_monitored_nests,
          here("eggshells",
               "field",
               "processed_data",
               "not_monitored_nests.csv"))

# Summary table -------------------------------------------------------------------------

summary_table_numbers <- master_and_nest %>%
  group_by(host_species,
           site) %>%
  summarize(min_num_host_egg = min(num_host_egg, na.rm = TRUE),
            max_num_host_egg = max(num_host_egg, na.rm = TRUE),
            mean_num_host_egg = mean(num_host_egg, na.rm = TRUE),
            sd_num_host_egg = sd(num_host_egg, na.rm = TRUE),
            min_num_host_hatch = min(num_host_hatch, na.rm = TRUE),
            max_num_host_hatch = max(num_host_hatch, na.rm = TRUE),
            mean_num_host_hatch = mean(num_host_hatch, na.rm = TRUE),
            sd_num_host_hatch = sd(num_host_hatch, na.rm = TRUE),
            min_host_hatching_sucess = min(host_hatching_sucess, na.rm = TRUE),
            max_host_hatching_sucess = max(host_hatching_sucess, na.rm = TRUE),
            mean_host_hatching_sucess = mean(host_hatching_sucess, na.rm = TRUE),
            sd_host_hatching_sucess = sd(host_hatching_sucess, na.rm = TRUE),
            min_host_proportion_hatch = min(host_proportion_hatch, na.rm = TRUE),
            max_host_proportion_hatch = max(host_proportion_hatch, na.rm = TRUE),
            mean_host_proportion_hatch = mean(host_proportion_hatch, na.rm = TRUE),
            sd_host_proportion_hatch = sd(host_proportion_hatch, na.rm = TRUE),
            min_host_time_to_hatch = min(host_time_to_hatch, na.rm = TRUE),
            max_host_time_to_hatch = max(host_time_to_hatch, na.rm = TRUE),
            mean_host_time_to_hatch = mean(host_time_to_hatch, na.rm = TRUE),
            sd_host_time_to_hatch = sd(host_time_to_hatch, na.rm = TRUE),
            min_host_hatching_period = min(host_hatching_period, na.rm = TRUE),
            max_host_hatching_period = max(host_hatching_period, na.rm = TRUE),
            mean_host_hatching_period = mean(host_hatching_period, na.rm = TRUE),
            sd_host_hatching_period = sd(host_hatching_period, na.rm = TRUE),
            min_num_bhco_egg = min(num_bhco_egg, na.rm = TRUE),
            max_num_bhco_egg = max(num_bhco_egg, na.rm = TRUE),
            mean_num_bhco_egg = mean(num_bhco_egg, na.rm = TRUE),
            sd_num_bhco_egg = sd(num_bhco_egg, na.rm = TRUE),
            min_num_bhco_hatch = min(num_bhco_hatch, na.rm = TRUE),
            max_num_bhco_hatch = max(num_bhco_hatch, na.rm = TRUE),
            mean_num_bhco_hatch = mean(num_bhco_hatch, na.rm = TRUE),
            sd_num_bhco_hatch = sd(num_bhco_hatch, na.rm = TRUE),
            min_bhco_hatching_sucess = min(bhco_hatching_sucess, na.rm = TRUE),
            max_bhco_hatching_sucess = max(bhco_hatching_sucess, na.rm = TRUE),
            mean_bhco_hatching_sucess = mean(bhco_hatching_sucess, na.rm = TRUE),
            sd_bhco_hatching_sucess = sd(bhco_hatching_sucess, na.rm = TRUE),
            min_bhco_proportion_hatch = min(bhco_proportion_hatch, na.rm = TRUE),
            max_bhco_proportion_hatch = max(bhco_proportion_hatch, na.rm = TRUE),
            mean_bhco_proportion_hatch = mean(bhco_proportion_hatch, na.rm = TRUE),
            sd_bhco_proportion_hatch = sd(bhco_proportion_hatch, na.rm = TRUE),
            min_bhco_time_to_hatch = min(bhco_time_to_hatch, na.rm = TRUE),
            max_bhco_time_to_hatch = max(bhco_time_to_hatch, na.rm = TRUE),
            mean_bhco_time_to_hatch = mean(bhco_time_to_hatch, na.rm = TRUE),
            sd_bhco_time_to_hatch = sd(bhco_time_to_hatch, na.rm = TRUE),
            min_bhco_hatching_period = min(bhco_hatching_period, na.rm = TRUE),
            max_bhco_hatching_period = max(bhco_hatching_period, na.rm = TRUE),
            mean_bhco_hatching_period = mean(bhco_hatching_period, na.rm = TRUE),
            sd_bhco_hatching_period = sd(bhco_hatching_period, na.rm = TRUE)) %>%
  mutate(across(.cols = min_num_host_egg:sd_bhco_hatching_period,
                 .fns = as.character)) %>%
  mutate(across(.cols = min_num_host_egg:sd_bhco_hatching_period,
                .fns = ~ str_replace(string = .,
                                     pattern = "Inf",
                              replacement = "NA")
  )) %>%
  mutate(across(.cols = min_num_host_egg:sd_bhco_hatching_period,
                .fns = ~ str_replace(string = .,
                              pattern = "-NA",
                              replacement = "NA")
  )) %>%
  mutate(across(.cols = min_num_host_egg:sd_bhco_hatching_period,
                .fns =  ~ str_replace(string = .,
                              pattern = "NaN",
                              replacement = "NA")
  )) %>%
  mutate(across(.cols = min_num_host_egg:sd_bhco_hatching_period,
                .fns =  ~str_replace(string = .,
                              pattern = "NA",
                              replacement = NA_character_)
  ))

summary_table <- master_and_nest %>%
  filter(predation_weather_effect == "yes") %>%
  select(host_species,
         site, 
         nest_id,
         predation_weather_effect) %>%
  distinct() %>%
  group_by(host_species,
           site) %>%
  summarize(num_nest_affected_by_predation_weather = n()
            ) %>%
  left_join(.,
            summary_table_numbers)

summary_table <- master_and_nest %>%
  filter(desertion_burial_effect == "yes") %>%
  select(host_species,
         site, 
         nest_id,
         desertion_burial_effect) %>%
  distinct() %>%
  group_by(host_species,
           site) %>%
  summarize(num_nest_affected_by_desertion_burial = n()
  ) %>%
  full_join(.,
            summary_table)

summary_table <- master_and_nest %>%
  filter(natural_parasitism == "yes") %>%
  select(host_species,
         site, 
         nest_id,
         natural_parasitism) %>%
  distinct() %>%
  group_by(host_species,
           site) %>%
  summarize(num_nest_parasitized = n()
  ) %>%
  full_join(.,
            summary_table)

summary_table <- master_and_nest %>%
  filter(addled_host_egg == "yes") %>%
  select(host_species,
         site, 
         nest_id,
         addled_host_egg) %>%
  distinct() %>%
  group_by(host_species,
           site) %>%
  summarize(num_nest_with_addled_host_egg= n()
  ) %>%
  full_join(.,
            summary_table)

summary_table <- master_and_nest %>%
  filter(addled_host_egg == "yes") %>%
  select(host_species,
         site, 
         nest_id,
         addled_host_egg) %>%
  distinct() %>%
  group_by(host_species,
           site) %>%
  summarize(num_nest_with_addled_host_egg = n()
  ) %>%
  full_join(.,
            summary_table)

summary_table <- master_and_nest %>%
  filter(host_nest_outcome == "complete sucess") %>%
  select(host_species,
         site, 
         nest_id,
         host_nest_outcome) %>%
  distinct() %>%
  group_by(host_species,
           site) %>%
  summarize(num_complete_sucess_host_nest = n()
  ) %>%
  full_join(.,
            summary_table)

summary_table <- master_and_nest %>%
  filter(host_nest_outcome == "partial sucess") %>%
  select(host_species,
         site, 
         nest_id,
         host_nest_outcome) %>%
  distinct() %>%
  group_by(host_species,
           site) %>%
  summarize(num_partial_sucess_host_nest = n()
  ) %>%
  full_join(.,
            summary_table)

summary_table <- master_and_nest %>%
  filter(host_nest_outcome == "total failure") %>%
  select(host_species,
         site, 
         nest_id,
         host_nest_outcome) %>%
  distinct() %>%
  group_by(host_species,
           site) %>%
  summarize(num_total_failure_host_nest = n()
  ) %>%
  full_join(.,
            summary_table)

summary_table <- master_and_nest %>%
  filter(bhco_nest_outcome == "complete sucess") %>%
  select(host_species,
         site, 
         nest_id,
         bhco_nest_outcome) %>%
  distinct() %>%
  group_by(host_species,
           site) %>%
  summarize(num_complete_sucess_bhco_nest = n()
  ) %>%
  full_join(.,
            summary_table)

summary_table <- master_and_nest %>%
  filter(bhco_nest_outcome == "partial sucess") %>%
  select(host_species,
         site, 
         nest_id,
         bhco_nest_outcome) %>%
  distinct() %>%
  group_by(host_species,
           site) %>%
  summarize(num_partial_sucess_bhco_nest = n()
  ) %>%
  full_join(.,
            summary_table)

summary_table <- master_and_nest %>%
  filter(bhco_nest_outcome == "total failure") %>%
  select(host_species,
         site, 
         nest_id,
         bhco_nest_outcome) %>%
  distinct() %>%
  group_by(host_species,
           site) %>%
  summarize(num_total_failure_bhco_nest = n()
  ) %>%
  full_join(.,
            summary_table)

summary_table <- summary_table %>% 
  arrange(host_species)

## Export data ##
write_csv(summary_table,
          here("eggshells",
               "field",
               "processed_data",
               "summary_table.csv"))
