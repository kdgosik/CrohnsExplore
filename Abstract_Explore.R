library(pubmed.mineR)
library(magrittr)
library(igraph)
library(dplyr)
library(purrr)
library(networkD3)

# first 400 abstracts
abstracts <- readabs("pubmed_result_crohns_search1.txt")

pmids <- abstracts@PMID

pubtator_output <- lapply(pmids, pubtator_function)

pubtator_output[[1]]

pubtator_output_t <- transpose(pubtator_output)

disease_data <- Map(function(a,b) data.frame(cbind(id = a, disease = b)), pubtator_output_t$PMID, pubtator_output_t$Diseases) %>%
  Filter(function(x) NROW(x) != 0, .) %>%
  Filter(function(x) NCOL(x) == 2, .) %>%
  do.call(rbind, .)

disease_data$disease <- tolower(disease_data$disease)
disease_data$disease <- gsub("diseases", "disease", disease_data$disease)
disease_data <- unique(disease_data)


genetic_data <- Map(function(a,b) data.frame(cbind(id = a, gene = b)), pubtator_output_t$PMID, pubtator_output_t$Genes) %>%
  Filter(function(x) NROW(x) != 0, .) %>%
  Filter(function(x) NCOL(x) == 2, .) %>%
  do.call(rbind, .)



# net <- graph_from_data_frame(genetic_data)
# plot(net)
simpleNetwork(genetic_data, zoom = TRUE)


