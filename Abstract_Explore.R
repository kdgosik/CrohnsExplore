library(pubmed.mineR)
library(magrittr)
library(igraph)
library(dplyr)
library(purrr)
library(networkD3)
library(foreach)
library(doParallel)

## helper functions ####

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

# first 400 abstracts
abstracts <- readabs("pubmed_result_crohns_search1.txt")
pmids <- abstracts@PMID


cl <- makeCluster(4)
registerDoParallel(cl)

pubmed_list <- foreach( i = pmids, .export = "pubtator_function" ) %dopar% {
  pubtator_function(i)
}

stopCluster(cl)

disease_data <- pubmed_list %>%
  transpose %$%
  Map(function(a,b) data.frame(cbind(id = a, disease = b)), PMID, Diseases) %>%
  Filter(function(x) NROW(x) != 0, .) %>%
  Filter(function(x) NCOL(x) == 2, .) %>%
  do.call(rbind, .) %>%
  mutate(disease = tolower(disease), 
         disease = gsub("diseases", "disease", disease)) %>%
  unique

for( i in seq_along(disease_data$disease) ) {
  disease_data$disease <- gsub(paste0(disease_data$disease[i], "s"),
                               disease_data$disease[i], disease_data$disease)
}
disease_data$disease[nchar(disease_data$disease) < 5] <- "REMOVE"
disease_data <- disease_data[disease_data$disease != "REMOVE", ]
disease_data <- unique(disease_data)


genetic_data <- pubmed_list %>%
  transpose %$%
  Map(function(a,b) data.frame(cbind(id = a, gene = b)), PMID, Genes) %>%
  Filter(function(x) NROW(x) != 0, .) %>%
  Filter(function(x) NCOL(x) == 2, .) %>%
  do.call(rbind, .) %>%
  mutate(gene = as.character(gene)) %>%
  unique

data("HGNCdata")
HGNCdata %<>% 
  mutate(Approved.Symbol = as.character(Approved.Symbol),
         Approved.Name = as.character(Approved.Name),
         Previous.Symbols = as.character(Previous.Symbols),
         Previous.Names = as.character(Previous.Names),
         Synonyms = as.character(Synonyms))

for( i in 1 : nrow(HGNCdata) ) {
  
  name_list <- as.character(HGNCdata[i, c("Approved.Name", "Previous.Symbols", "Previous.Names", "Synonyms")])
  genetic_data$gene[tolower(genetic_data$gene) %in% tolower(name_list)] <- paste(HGNCdata[i, "Approved.Symbol"])

}

for( i in 1 : nrow(HGNCdata) ) {
  
  genetic_data$gene[tolower(genetic_data$gene) %in% tolower(genetic_data$gene)[i]] <- genetic_data$gene[i]
  
}

tnf_list <- c("TNF", "TNFa", "TNF-a", "Tumor necrosis factor (TNF)-alpha",
              "TNF-alpha", "TNF alpha", "TNFA", "tumor necrosis factor",
              "tumor necrosis factor-a", "tumor necrosis factor alpha",  
              "tumor necrosis factor-alpha", "tumor necrosis factor a",
              "Tumor Necrosis Factor Alpha", "tumour necrosis factor",
              "tumour necrosis factor-a")
genetic_data$gene[genetic_data$gene %in% tnf_list] <- "TNF"

anti_tumor <- c("anti-tumour necrosis factor a", "anti-tumour necrosis factor",
                "anti-tumor necrosis factor", "Anti-tumor necrosis factor-a",
                "anti-tumour necrosis factor")
genetic_data$gene[genetic_data$gene %in% anti_tumor] <- "anti-tumour necrosis factor"

genetic_data$gene[genetic_data$gene %in% "interleukin-1"] <- "IL-1"
genetic_data$gene[genetic_data$gene %in% "interleukin (IL)-6"] <- "IL-6"
genetic_data$gene[genetic_data$gene %in% c("interleukin-10", "Interleukin-10")] <- "IL-10"
genetic_data$gene[genetic_data$gene %in% "interleukin 17"] <- "IL-17"
genetic_data$gene[genetic_data$gene %in% "Interleukin-18"] <- "IL-18"
genetic_data$gene[genetic_data$gene %in% c("interleukin-23", "interleukin 23", "IL23", "Interleukin 23")] <- "IL-23"
genetic_data$gene[genetic_data$gene %in% c("interleukin 28A", "IL28A", "IL-28a")] <- "IL-28a"
genetic_data$gene[genetic_data$gene %in% "Interleukin (IL)-32"] <- "IL-32"
genetic_data$gene[genetic_data$gene %in% "TLR)2"] <- "TLR2"
genetic_data$gene[genetic_data$gene %in% "VDR"] <- "vitamin D receptor"
genetic_data$gene[genetic_data$gene %in% c("stat1", "Stat1IEC-KO")] <- "STAT1"
genetic_data$gene[genetic_data$gene %in% "signal transducer and activator of transcription 3"] <- "STAT3"
genetic_data$gene[genetic_data$gene %in% c("stat1", "Stat1IEC-KO")] <- "STAT1"
genetic_data$gene[genetic_data$gene %in% "P-gp"] <- "P-glycoprotein" 
genetic_data$gene[genetic_data$gene %in% "CRP  <  3"] <- "CRP"

genetic_data <- unique(genetic_data)


chemical_data <- pubmed_list %>%
  transpose %$%
  Map(function(a,b) data.frame(cbind(id = a, chemicals = b)), PMID, Chemicals) %>%
  Filter(function(x) NROW(x) != 0, .) %>%
  Filter(function(x) NCOL(x) == 2, .) %>%
  do.call(rbind, .) %>%
  mutate(chemicals = tolower(chemicals)) %>%
  unique

for( i in seq_along(chemical_data$chemicals) ) {
  chemical_data$chemicals <- gsub(paste0(chemical_data$chemicals[i], "s"),
                                     chemical_data$chemicals[i], chemical_data$chemicals)
}

chemical_data <- unique(chemical_data)


# net <- graph_from_data_frame(genetic_data)
# plot(net)
simpleNetwork(disease_data, zoom = TRUE)
simpleNetwork(genetic_data, zoom = TRUE)
simpleNetwork(chemical_data, zoom = TRUE)


disease_gene <- full_join(disease_data, genetic_data) %>%
  select(-id) %>%
  na.omit

disease_gene %>% 
  group_by(disease, gene) %>% 
  summarise(N = n()) %>% 
  arrange(-N)

disease_chemical <- full_join(disease_data, chemical_data) %>%
  select(-id) %>%
  na.omit

chemical_gene <- full_join(chemical_data, genetic_data) %>%
  select(-id) %>%
  na.omit


tmp <- disease_gene %>%
  select(from = disease, to = gene) %>%
  rbind(. , {disease_chemical %>% select(from = disease, to = chemicals)}) %>%
  rbind(. , {chemical_gene %>% select(from = chemicals, to = gene)})
  
  