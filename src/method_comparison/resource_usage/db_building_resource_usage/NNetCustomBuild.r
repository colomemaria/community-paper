library(OmnipathR)
library(nichenetr)
library(tidyverse)
library(mlrMBO)
library(parallelMap)
library(dplyr)
library(community)

interactionFormatTransf <- function(InputDf, InteractionType){
  
    OutputInt <- tibble(from = character(), to = character(), 
        source = character(), database = character())  
    
    n <- nrow(InputDf)
    sources <- dplyr::pull(InputDf, sources)
    sourceNodes <- dplyr::pull(InputDf, from)
    targetNodes <- dplyr::pull(InputDf, to)
    
    for (i in seq(n)){
        currentSources <- unlist(strsplit(sources[i],";"))
        for (j in seq(length(currentSources))){
            OutputInt <- add_row(OutputInt, 
                from = sourceNodes[i] , 
                to = targetNodes[i],  
                # source = paste(currentSources[j], InteractionType, sep="_"),
                source = currentSources[j],
                database = currentSources[j]) 
        }
    }
    
    return(OutputInt)
}

data(LR_database)

lr <- LR_database



lr <- lr %>%
    select(Ligand,Receptor,sources) %>%
    rename(from=Ligand, to=Receptor) %>% 
    filter(from != to) %>% 
    distinct()

lr_Network_Omnipath <- 
    lr %>%
    interactionFormatTransf(InteractionType="LigrecExtra") %>%
    dplyr::distinct() 

saveRDS(lr_Network_Omnipath, 
    "NNET_Custom/lig_rec_sources.rds")



## We next get protein-protein interactions from the different datasets availabe
## in Omnipath
AllInteractions <- 
    import_post_translational_interactions(exclude = "ligrecextra") %>% 
    dplyr::select(source_genesymbol, target_genesymbol, sources) %>% 
    dplyr::rename(from=source_genesymbol, to=target_genesymbol) %>% 
    dplyr::filter(from != to) %>% 
    dplyr::distinct() 

ligands <- unique(pull(lr, from))

# sig_Network_Omnipath <- sig_Network_Omnipath_raw

# ## Original Omnipath interactions
# sig_Network_Omnipath <- 
#     interactionFormatTransf(AllInteractions, InteractionType="Signalling") %>%
#     dplyr::distinct() 

# ## I have to remove self-interactions in the signaling network
# sig_Network_Omnipath <- sig_Network_Omnipath %>% 
#     dplyr::filter(from != to)

# # ## I also have to remove interactions going to ligands. See Methods Nichenet 
# # ## paper
# # sig_Network_Omnipath <- sig_Network_Omnipath %>% 
# #     dplyr::filter(!(to %in% ligands))

# ## There are in addition some records containing not input gene, we remove them
# ## since they are giving problems with running the model.
# sig_Network_Omnipath <- sig_Network_Omnipath %>% 
#     dplyr::filter(from != "") %>% 
#     dplyr::filter(to != "")


# ## We also remove signaling interactions that are already in the lig-receptor 
# ## network. 
# sig_Network_Omnipath <- dplyr::anti_join(
#   sig_Network_Omnipath, 
#   lr_Network_Omnipath, 
#   by = c("from" = "from", "to" = "to"))

# nrow(sig_Network_Omnipath)

## Original Omnipath interactions
sig_Network_Omnipath <- 
    interactionFormatTransf(AllInteractions, InteractionType="Signalling") %>%
    dplyr::distinct() 

## I have to remove self-interactions in the signaling network
sig_Network_Omnipath <- sig_Network_Omnipath %>% 
    dplyr::filter(from != to)

# ## I also have to remove interactions going to ligands. See Methods Nichenet 
# ## paper
sig_Network_Omnipath <- sig_Network_Omnipath %>% 
    dplyr::filter(!(to %in% ligands))

## There are in addition some records containing not input gene, we remove them
## since they are giving problems with running the model.
sig_Network_Omnipath <- sig_Network_Omnipath %>% 
    dplyr::filter(from != "") %>% 
    dplyr::filter(to != "")


## We also remove signaling interactions that are already in the lig-receptor 
## network. 
sig_Network_Omnipath <- dplyr::anti_join(
  sig_Network_Omnipath, 
  lr_Network_Omnipath, 
  by = c("from" = "from", "to" = "to"))

nrow(sig_Network_Omnipath)

202063 == 163533

saveRDS(sig_Network_Omnipath, 
    "NNET_Custom/sig_Network.rds")

gr_Interactions_Omnipath <- 
    dorothea(dorothea_levels = c("A","B","C")) %>%  
    select(source_genesymbol, target_genesymbol, sources) %>%
    rename(from=source_genesymbol, to=target_genesymbol) %>% 
    filter(from != to) %>%
    distinct()  

gr_Network_Omnipath <- 
    interactionFormatTransf(
        gr_Interactions_Omnipath, 
        InteractionType="Dorothea") %>%
    dplyr::distinct() 
nrow(gr_Network_Omnipath)
## [1] 113897

82767

saveRDS(gr_Network_Omnipath,
    "NNET_Custom/gr_Network.rds")

expression_settings_validation <- readRDS(url("https://zenodo.org/record/8010790/files/expression_settings"))

# index <- which(!unlist(lapply(expression_settings_validation, 
#     function(x) any(x$from != "IFNA1"))))

# expression_settings_validation <- expression_settings_validation[-index]

ls()

save(list=c("expression_settings_validation", "lr_Network_Omnipath", "sig_Network_Omnipath", "gr_Network_Omnipath"), 
     file="myvariables.RData")
