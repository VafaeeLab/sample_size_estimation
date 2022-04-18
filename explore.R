library(curatedMetagenomicData)
library(tidyverse)
library(ExperimentHub)


colnames(sampleMetadata)

head(sampleMetadata)

length(unique(sampleMetadata$study_name))
unique(sampleMetadata$study_name)
length(unique(sampleMetadata$NCBI_accession))

getMetada

study_df <- data.frame(study_name = unique(sampleMetadata$study_name)) %>%
  separate(col = "study_name", into = c("author", "year"), sep = "_", remove = FALSE) %>%
  mutate(year = strtoi(gsub("a|b", "", year))) %>%
  arrange(desc(year))
write.csv(study_df, "study_df.csv", row.names = FALSE)

#AsnicarF_2021
sample_metadata_AF2021 <- sampleMetadata %>%
  filter(study_name == "AsnicarF_2021")

unique(sampleMetadata$DNA_extraction_kit)
unique(sampleMetadata$disease)

summary(sample_metadata_AF2021$age)

length(unique(sample_metadata_AF2021$sample_id))
length(unique(sample_metadata_AF2021$subject_id))
levels(factor(sample_metadata_AF2021$disease))


sample_metadata_AF2017 <- sampleMetadata %>%
  filter(study_name == "AsnicarF_2017")
length(unique(sample_metadata_AF2017$sample_id))


sample_metadata <- sampleMetadata %>%
  filter(disease == "melanoma;metastases")


curatedMetagenomicData("AsnicarF_2017.relative_abundance", dryrun = FALSE, rownames = "short")

sample_metadata <- sampleMetadata %>%
  filter(study_name == "RubelMA_2020")

data <- sample_metadata_AF2021 %>%
  returnSamples("relative_abundance")

assayNames(data)
reads <- assays(data)$relative_abundance
meta_data <- data.frame(colData(data))


leeka_2022 <- read.table("LeeKA_2022_metadata.tsv", header = TRUE, sep = "\t")
levels(factor(leeka_2022$RECIST))
levels(factor(leeka_2022$ORR))
levels(factor(leeka_2022$PFS12))

sample_metadata_gk <- sampleMetadata %>%
  filter(study_name == "GopalakrishnanV_2018")



##LeeKA_2022
sample_metadata_LeeKA2022 <- sampleMetadata %>%
  filter(study_name == "LeeKA_2022")

data <- sample_metadata_LeeKA2022 %>%
  returnSamples("relative_abundance")

assayNames(data)
reads <- assays(data)$relative_abundance
meta_data <- data.frame(colData(data))



#AsnicarF_2021
sample_metadata <- sampleMetadata %>%
  filter(study_name == "AsnicarF_2021")

length(unique(sample_metadata$sample_id))
length(unique(sample_metadata$subject_id))
levels(factor(sample_metadata$disease))

data <- sample_metadata %>%
  returnSamples("relative_abundance")

assayNames(data)
reads <- assays(data)$relative_abundance
meta_data <- data.frame(colData(data))


data <- curatedMetagenomicData("AsnicarF_2021.relative_abundance", 
                               dryrun = FALSE, counts = FALSE, rownames = "long") %>%
  mergeData()

print(resourceTitles)

str_subset(resourceTitles, "AsnicarF_2021.relative_abundance")
eh_data <- ExperimentHub() %>%
  query("AsnicarF_2021.relative_abundance")

#WindTT_2020
#gives error
sample_metadata <- sampleMetadata %>%
  filter(study_name == "WindTT_2020")

length(unique(sample_metadata$sample_id))
length(unique(sample_metadata$subject_id))
levels(factor(sample_metadata$disease))

data <- sample_metadata %>%
  returnSamples("relative_abundance")
# snapshotDate(): 2021-05-18
# Error: dataType of list elements is different

assayNames(data)
reads <- assays(data)$relative_abundance
meta_data <- data.frame(colData(data))


#LeeKA_2022
#gives error
sample_metadata <- sampleMetadata %>%
  filter(study_name == "LeeKA_2022")

length(unique(sample_metadata$sample_id))
length(unique(sample_metadata$subject_id))
levels(factor(sample_metadata$disease))

data <- sample_metadata %>%
  returnSamples("relative_abundance")
# snapshotDate(): 2021-05-18
# Error: dataType of list elements is different

eh_data <- ExperimentHub() %>%
  query("LeeKA_2022.relative_abundance")
#nothing in this

assayNames(data)
reads <- assays(data)$relative_abundance
meta_data <- data.frame(colData(data))
