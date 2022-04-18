BiocManager::install("curatedMetagenomicData")

#https://askubuntu.com/questions/890027/installing-gsl-libraries-in-ubuntu-16-04-via-terminal
BiocManager::install("DirichletMultinomial")

BiocManager::install("waldronlab/curatedMetagenomicData", 
                     dependencies = TRUE, build_vignettes = TRUE,
                     force = TRUE)

# remove.packages("curatedMetagenomicData")

BiocManager::install("ExperimentHub")
