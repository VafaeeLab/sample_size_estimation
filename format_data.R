library(tidyverse)

relative_abundance <- read.table("data/species_abundance.txt", header = TRUE, sep = "\t")

relative_abundance <- relative_abundance %>%
  select(-c(2)) %>%
  column_to_rownames("clade_name")
colnames(relative_abundance) <- gsub(".", "_", 
                                     gsub("_profile", "", colnames(relative_abundance), fixed = TRUE),
                                     fixed = TRUE)
relative_abundance <- data.frame(t(relative_abundance))

write.csv(relative_abundance, "data/formatted_data.csv")



all_level_relative_abundance <- read.table("data/merged_abundance_table.txt", header = TRUE, sep = "\t")
all_level_relative_abundance <- all_level_relative_abundance %>%
  select(-c(2)) %>%
  column_to_rownames("clade_name")
colnames(all_level_relative_abundance) <- gsub(".", "_", 
                                               gsub("_profile", "", colnames(all_level_relative_abundance), fixed = TRUE),
                                               fixed = TRUE)
all_level_relative_abundance <- data.frame(t(all_level_relative_abundance))

write.csv(all_level_relative_abundance, "data/all_level_formatted_data.csv")