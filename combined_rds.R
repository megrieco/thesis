library(dplyr)
library(purrr)

project_dir <- "/home/mgrieco/Thesis/Test-231003/"

combined_results <- list.files(path=paste0(project_dir,"Results"),
                                  pattern = ".Rds$", full.names = TRUE, recursive=TRUE) %>%
  map_dfr(readRDS) 

#save combined results
saveRDS(object=combined_results,file=paste0(project_dir,"Results/combined_results.Rds"))