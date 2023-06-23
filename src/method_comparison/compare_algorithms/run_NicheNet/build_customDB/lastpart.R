library(OmnipathR)
library(nichenetr)
library(tidyverse)
library(mlrMBO)
library(parallelMap)
library(dplyr)

load("myvariables.RData")

ls()

All_sources <- unique(c(lr_Network_Omnipath$source,
    sig_Network_Omnipath$source, gr_Network_Omnipath$source))

my_source_weights_df <- 
     tibble(source = All_sources, weight = rep(1,length(All_sources)))

additional_arguments_topology_correction <- 
    list(source_names = my_source_weights_df$source %>% unique(), 
        algorithm = "PPR", 
        correct_topology = FALSE,
        lr_network = lr_Network_Omnipath, 
        sig_network = sig_Network_Omnipath, 
        gr_network = gr_Network_Omnipath, 
        settings = lapply(expression_settings_validation, 
            convert_expression_settings_evaluation), 
        secondary_targets = FALSE, 
        remove_direct_links = "no", 
        cutoff_method = "quantile")

nr_datasources <- additional_arguments_topology_correction$source_names %>% 
    length()

obj_fun_multi_topology_correction = makeMultiObjectiveFunction(name = "nichenet_optimization",
    description = "data source weight and hyperparameter optimization: expensive black-box function", 
    fn = model_evaluation_optimization, 
    par.set = makeParamSet(
        makeNumericVectorParam("source_weights", len = nr_datasources, 
            lower = 0, upper = 1, tunable = FALSE), 
        makeNumericVectorParam("lr_sig_hub", len = 1, lower = 0, upper = 1, 
            tunable = TRUE),  
        makeNumericVectorParam("gr_hub", len = 1, lower = 0, upper = 1, 
            tunable = TRUE),  
        makeNumericVectorParam("ltf_cutoff", len = 1, lower = 0.9, 
            upper = 0.999, tunable = TRUE),  
        makeNumericVectorParam("damping_factor", len = 1, lower = 0.01, 
            upper = 0.99, tunable =TRUE)), 
    has.simple.signature = FALSE,
    n.objectives = 4, 
    noisy = FALSE,
    minimize = c(FALSE,FALSE,FALSE,FALSE))

optimization_results = 
    lapply(1,mlrmbo_optimization, obj_fun = obj_fun_multi_topology_correction, 
           niter = 8, ncores = 8, nstart = 160, 
           additional_arguments = additional_arguments_topology_correction)

saveRDS(optimization_results, "../../../../../results/method_comparison/build_customDB/NicheNet/Optimization_results.rds")
