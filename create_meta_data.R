library(tidyverse)

file_list <- read.table("file_list.txt", header = FALSE)
colnames(file_list) <- "file_name"

downloaded_files <- read.table("fastq_downloaded.txt", header = FALSE)
colnames(downloaded_files) <- "file_name"

download_missing <- file_list %>%
  anti_join(downloaded_files)

samples <- downloaded_files %>%
  mutate(sample = gsub("_1.fastq.gz", "", file_name, fixed = TRUE)) %>%
  mutate(sample = gsub("_2.fastq.gz", "", sample, fixed = TRUE)) %>%
  select(sample) %>%
  unique() %>%
  separate(col = "sample", into = c("NCBI_accession", "sample_id"), 
           sep = "_", extra = "merge", remove = FALSE) %>%
  select(c(sample, sample_id, NCBI_accession)) %>%
  arrange(sample_id)

original_meta_data <- read.table("LeeKA_2022_metadata.tsv", header = TRUE, sep = "\t")

missing_data <- original_meta_data %>%
  anti_join(samples)

meta_data <- samples %>%
  inner_join(original_meta_data) %>%
  select(sample, sample_id, NCBI_accession, study_condition, disease, age, age_category, 
         subcohort, treatment, toxicity_above_zero, RECIST, ORR) %>%
  rename(c("ICIresponder" = "ORR"))

write.csv(meta_data, "data/meta_data.csv", row.names = FALSE)
