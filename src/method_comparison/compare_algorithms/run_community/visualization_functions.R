##### Mean weights vs number of interactions

# helper function for number of interactions vs mean interaction weight plot
# claclulate mean weight of (good) interactions for each pair of interactions cell types. 
# E.g. mean weight of interactions for "T to B", mean weight of interactions "DC to T", etc.

mean_weights <- function(my_interactions){
        
        # extract health status
        health_status <- unique(my_interactions$anno_samples$health_status)
        
        # extract cell types
        cell_types <- unique(my_interactions$anno_cells$cell_type)
        
        # extract good interactions
        idx_good <- my_interactions$anno_interactions$passed_QC_filter
        
        my_weights <- lapply(health_status
                             ,function(hs){
                                     
                                     idx_hs <- my_interactions$anno_samples$health_status == hs
                                     
                                     # number of interactions as heatmap
                                     df <- as.data.frame(matrix(,nrow = length(cell_types)
                                                                ,ncol = length(cell_types)
                                     )
                                     )
                                     rownames(df) <- cell_types
                                     colnames(df) <- cell_types
                                     
                                     # populate the data frame
                                     for(send in cell_types){
                                             idx_send <- my_interactions$anno_interactions$sending_cell_type == send
                                             
                                             for(rec in cell_types){
                                                     idx_rec <- my_interactions$anno_interactions$receiving_cell_type == rec
                                                     
                                                     # extract weights
                                                     weights <- my_interactions$weights[idx_good & idx_send & idx_rec,idx_hs]
                                                     
                                                     # remove interactions with all-zero weights
                                                     idx_all_zero <- rowSums(weights) == 0
                                                     weights <- weights[!idx_all_zero,]
                                                     
                                                     # calculate mean wieghts of interactions
                                                     df[send,rec] <- mean(unlist(weights))
                                                     
                                             }
                                     }
                                     df
                             }
        )
        
        names(my_weights) <- health_status
        
        return(my_weights)
}

# helper function for number of interactions vs mean interaction weight plot
# claclulate number of (good) interactions for each pair of interactions cell types. 
# E.g. number of interactions for "T to B", number of interactions "DC to T", etc.

number_interactions <- function(my_interactions){
        
        # extract health status
        health_status <- unique(my_interactions$anno_samples$health_status)
        
        # extract cell types
        cell_types <- unique(my_interactions$anno_cells$cell_type)
        
        # extract good interactions
        idx_good <- my_interactions$anno_interactions$passed_QC_filter
        
        my_weights <- lapply(health_status
                             ,function(hs){
                                     
                                     idx_hs <- my_interactions$anno_samples$health_status == hs
                                     
                                     # number of interactions as heatmap
                                     df <- as.data.frame(matrix(,nrow = length(cell_types)
                                                                ,ncol = length(cell_types)
                                     )
                                     )
                                     rownames(df) <- cell_types
                                     colnames(df) <- cell_types
                                     
                                     # populate the data frame
                                     for(send in cell_types){
                                             idx_send <- my_interactions$anno_interactions$sending_cell_type == send
                                             
                                             for(rec in cell_types){
                                                     idx_rec <- my_interactions$anno_interactions$receiving_cell_type == rec
                                                     
                                                     # extract weights
                                                     weights <- my_interactions$weights[idx_good & idx_send & idx_rec,idx_hs]
                                                     
                                                     # remove interactions with all-zero weights
                                                     idx_all_zero <- rowSums(weights) == 0
                                                     weights <- weights[!idx_all_zero,]
                                                     
                                                     # calculate number of interactions
                                                     df[send,rec] <- nrow(weights)
                                                     
                                             }
                                     }
                                     df
                             }
        )
        
        names(my_weights) <- health_status
        
        return(my_weights)
}

# plot mumber of interactions vs mean interaction weights per cell type to cell type interaction

plot_nrInt_vs_meanW_perCellType <- function(my_interactions
                                            ,interaction_type # named vector: values are interaction type, names are intereaction IDs
                                            ,colors = NA # named vector
                                            ,add_area = TRUE
                                            ,font_size = 26 # size of the axes text
                                            ,label_font_size = 6 # size of the text labels on the plots (i.e. "T to B", etc.)
                                            ,ylim = c(-4,0)
                                            ,xlim = c(-20,700)
){
        
        # calculate mean weight of interactions
        mean_weights_goodInteractions <- mean_weights(my_interactions)
        
        # calculate number of interactions
        number_goodInteractions <- number_interactions(my_interactions)
        
        # extract health status
        health_status <- unique(my_interactions$anno_samples$health_status)
        
        flatten_df <- function(df){
                df_flat <- do.call(c,df)
                names(df_flat) <- unlist(lapply(colnames(df)
                                                ,function(c){
                                                        lapply(rownames(df)
                                                               ,function(r){paste(r,c,sep = " to ")})
                                                })
                )
                df_flat
        }
        
        lapply(health_status
               ,function(hs){
                       
                       my_df <- data.frame(mean_interaction_weight = flatten_df(mean_weights_goodInteractions[[hs]])
                                           ,number_of_interactions = flatten_df(number_goodInteractions[[hs]])
                                           ,interaction_ID = names(flatten_df(number_goodInteractions[[hs]]))
                       )
                       my_df$interaction_type <- interaction_type[my_df$interaction_ID]
                       
                       print(str(my_df))
                       
                       p <- ggplot(data = my_df
                                   ,aes(x = number_of_interactions
                                        ,y = log10(mean_interaction_weight)
                                        ,color = interaction_type
                                   )
                       )+
                               geom_point(size = 2)+
                               ylim(ylim)+
                               xlim(xlim)+
                               xlab("number of interactions")+
                               ylab("log10 mean w")+
                               ggtitle(hs)+
                               theme_bw()+
                               theme(text = element_text(size = font_size) # axes
                                     ,legend.position = 'bottom'
                                     ,plot.title = element_text(hjust = 0.5)) 
                       
                       # add custom colors
                       if(typeof(colors) != "logical"){
                               p <- p+
                                       scale_color_manual("interacting cell types"
                                                          ,values = colors
                                       )+
                                       scale_fill_manual(values = colors)
                       }
                       
                       # add ellipse area
                       if(add_area){
                               p <- p+
                                       stat_ellipse(aes(fill = interaction_type)
                                                    ,geom = "polygon"
                                                    ,type = "t"
                                                    ,alpha = 0.2
                                                    ,show.legend = FALSE
                                                    ,lwd = 0
                                       )
                       }
                       
                       p+
                               geom_text_repel(aes(label = interaction_ID)
                                               ,color = "black"
                                               ,alpha = 0.75
                                               ,size = label_font_size
                               )
                       
               })
}

##### Volcano plot
plot_vulcano <- function(my_interactions
                         ,colors = c("red3" # upregulated
                                     ,"gray90" # unsignificant
                                     ,"lightslateblue" # downregulated
                         )
                         ,font_size = 18)
{
        # good quality interactions
        idx_good <- my_interactions$anno_interactions$passed_QC_filter
        
        # threshold log2FC
        threshold_log2FC <- my_interactions$thresholds$threshold_log2FC
        
        # threshold p.adj
        threshold_p.adj <- interactions$thresholds$threshold_fdr
        
        df <- data.frame(log2FC = my_interactions$anno_interactions$log2FC_weights[idx_good]
                         ,y = -log10(my_interactions$anno_interactions$p.adj[idx_good])
                         ,significant = my_interactions$anno_interactions$sign[idx_good]
        )
        
        df$direction <- "unchanged"
        df$direction[df$significant & (df$log2FC > threshold_log2FC)] <- "up"
        df$direction[df$significant & (df$log2FC < -threshold_log2FC)] <- "down"
        df$direction <- factor(df$direction
                               ,levels = c("up"
                                           ,"unchanged"
                                           ,"down")
                               ,ordered = TRUE)
        
        xlab <- "log2 fold change"
        ylab <- "-log10 p.adj"
        
        xlim <- c(-max(abs(df$log2FC))
                  ,max(abs(df$log2FC))
        )
        ylim <- c(0, max(df$y))
        
        p <- ggplot(data = df
                    ,aes(x = log2FC
                         ,y = y
                         ,color = direction
                         ,size = significant
                    ))+
                geom_point(alpha = 0.5
                           ,show.legend = FALSE)+
                scale_color_manual(values = colors)+
                scale_size_manual(values = c(0.5, 1.5))+
                guides(size = "none"
                       ,shape = "none")+
                xlab(xlab)+
                ylab(ylab)+
                xlim( xlim )+
                ylim(ylim)+
                theme_bw()+
                theme(text = element_text(size=font_size))+
                geom_vline(xintercept = threshold_log2FC
                           ,lty = 2
                           ,color = "gray")+
                geom_vline(xintercept = -threshold_log2FC
                           ,lty = 2
                           ,color = "gray")+
                geom_hline(yintercept = -log10(threshold_p.adj)
                           ,lty = 2
                           ,color = "gray")
        p
}

##### Heatmap 

# NEW:
# added option to add parameters to the Heatmap function

#' @title plot_heatmap
#' 
#' @description plots heatmap of chosen parameter
#' 
#' @param comm_result A dataframe of common interactions
#' @param which_interactions A logical vector with length of all interactions, indicating which interactions to include in the heatmap. NULL if top and by_param are set or "all"
#' @param by_param  A string indicating which column of the common interactions dataframe to use for filtering interactions, e.g. "passed_QC_filter", "passed_FDR_threshold", "passed_log2FC_threshold", "sign"
#' @param values_to_plot A string indicating which column of the common interactions dataframe to use as the values in the heatmap. e.g. "interactions_weights", "expr_l_s_active", "expr_r_r_active", "nr_l_s_active", "nr_r_r_active", "phi", "phi_l_s", "phi_r_r", "p", "p_l_s", "p_r_r"
#' @param row_font_size Font size for the row labels
#' @param color_case A string indicating the color for case samples in the heatmap
#' @param color_control A string indicating the color for control samples in the heatmap
#' @param ... Parameters passed to the ComplexHeatmap::Heatmap function
#' 
#' @export
#' @examples
#' # load example data
#' data("comm_result")
#' # calculate general statistics
#' comm_result <- general_stat(comm_result)
#' # plot heatmap
#' plot_heatmap(comm_result = comm_result
#'             ,which_interactions = NULL
#'            ,by_param = "passed_QC_filter"
#'           ,values_to_plot = "interactions_weights"
#'         ,row_font_size = 8
#'      ,color_case = "#7C001F" # "darkred"
#'  ,color_control = "#7AC5CD" # "CadetBlue3"
#' )
#'
plot_heatmap <- function(comm_result
                         ,which_interactions = NULL 
                         # NULL if top and by_param are set 
                         # or "all" 
                         # or boolean vector with length of all interactions
                         # e.g. "ADAM10" %in% Single_Cell_Result$anno_interactions$ligand_gene_name 
                         # or "T" %in% Single_Cell_Result$anno_interactions$sending_cell_type
                         #,top = NULL 
                         ,by_param 
                         # passed_QC_filter
                         # passed_FDR_threshold
                         # passed_log2FC_threshold
                         # sign
                         ,values_to_plot 
                         # weights, 
                         # expr_l_s_active, 
                         # expr_r_r_active, 
                         # nr_l_s_active, 
                         # nr_r_r_active, 
                         # phi, 
                         # phi_l_s, 
                         # phi_r_r, 
                         # p, 
                         # p_l_s, 
                         # p_r_r
                         ,centered = FALSE
                         ,row_font_size = 8
                         ,column_font_size = 8
                         ,legend_title_font_size = 8
                         ,labels_font_size = 8
                         ,color_case = "#7C001F" # "darkred"
                         ,color_control = "#7AC5CD" # "CadetBlue3"
                         ,color_values = NULL # or 
                         ,...
){
        
        df <- as.matrix(comm_result[[values_to_plot]])
        #print(str(df))
        if(is.null(which_interactions)){
                idx_interactions <- comm_result$anno_interactions[,by_param] & (!is.na(comm_result$anno_interactions[,by_param]))
        } else if(class(which_interactions) == 'logical'){
                idx_interactions <- which_interactions
        } else if(which_interactions == "all"){
                idx_interactions <- rep(TRUE
                                        ,nrow(df))
        } else {
                stop("ERROR: parameter which_interactions can be either NULL or 'all' or a logical vector.")
        }
        #print(str(idx_interactions) 
        
        col_samples <- sapply(comm_result$anno_samples$case_or_control ## add to object; create column with "case" and "control"
                              ,function(i){
                                      ifelse(grepl("case"
                                                   ,i)
                                             ,color_case
                                             ,color_control
                                      )
                              })
        #print(str(col_samples))
        
        
        if(centered){
                my_means <- matrix(nrow = nrow(df)
                                   ,ncol = ncol(df)
                                   ,rowMeans(df)
                )
                #print(my_means)
                df <- df - my_means
                
                my_max <- matrix(nrow = nrow(df)
                                 ,ncol = ncol(df)
                                 ,apply(df
                                        ,1
                                        ,function(x)max(abs(x)))
                )
                df <- df / my_max
        }
        
        if(is.null(color_values)){
                if(centered){
                        my_color <- circlize::colorRamp2(c(-1,0,1), c("lightslateblue","white", "red3"))
                } else my_color <- circlize::colorRamp2(seq(0, 1, length = 5), c("white", "red", "red4",  "darkred", "black"))
        } else my_color <-  color_values
        
        if(centered){
                title  <- paste("centered",values_to_plot)
        } else title <- values_to_plot
        
        h <- Heatmap(df[idx_interactions,]
                     ,name = title
                     ,col = my_color
                     ,column_names_side = "top"
                     ,column_names_gp = gpar(col = col_samples
                                             ,fontsize = column_font_size)
                     ,heatmap_legend_param = list(direction = "horizontal"
                                                  ,title_position = "lefttop"
                                                  ,title_gp = gpar(fontsize = legend_title_font_size)
                                                  ,labels_gp = gpar(fontsize = labels_font_size)
                     )
                     ,row_names_gp = gpar(fontsize = row_font_size)
                     ,row_names_max_width = max_text_width(rownames(df[idx_interactions,])
                                                           ,gp = gpar(fontsize = row_font_size))
                     ,...
        )
        plot(h
             ,heatmap_legend_side = "bottom")
}

##### Barplot down- and up- by interacting cell types

# stacked bar interaction classes: up- and downregulated
plot_stacked_bar <- function(my_interactions
                             ,idx_up
                             ,idx_down
                             ,idx_sign
                             ,interaction_type
                             ,colors
                             ,font_size = 18 # size of the axes fonts
){
        my_df <- my_interactions$anno_interactions[,c("interaction_ID","sign")]
        
        # assign directionality
        my_df$direction <- NA
        my_df$direction[idx_up & idx_sign]  <- "up" #"up-regulated"
        my_df$direction[idx_down & idx_sign] <- "down" #"down-regulated"
        
        # add interaction type
        my_df$interaction_type <- interaction_type[my_df$interaction_ID]
        
        # subset to only significant interactions
        my_df <- my_df[idx_sign,]
        my_df <- as.data.frame(table(my_df[,c("direction"
                                              ,"interaction_type")]))
        p <- ggplot(my_df
                    ,aes(fill=interaction_type
                         ,y=Freq
                         ,x=direction)
        ) + 
                geom_bar(position="stack", stat="identity")+
                scale_fill_manual("interacting cell types"
                                  ,values = colors)+
                ylab("number of interactions")+
                xlab("")+
                theme_bw()+
                theme(text = element_text(size=font_size))
        
        p
}

##### Cell type network

# cell type network plot of differential interactions
plor_celltype_network <- function(my_interactions
                                  ,amplify_edgeWidth = 50 # magnification for the edge width
                                  ,amplify_colorResolution = 20 # if color_palette is used
                                  ,nr_colors = 10 # if color_palette is used
                                  ,color_palette = c("red", "red4", "black")
                                  ,edge.color = NULL # custom edge color (named vector of colors, names are edge IDs)
                                  ,vertex.label.cex = 2
                                  ,vertex.shape="none"
                                  ,vertex.size = 10
                                  ,edge.arrow.size = 1
                                  ,title_cex = 1
                                  ,verbose = FALSE
                                  ,...
){
        
        # extract health status
        health_status <- unique(my_interactions$anno_samples$health_status)
        
        # prepare the data for plotting
        data <- lapply(health_status
                       ,function(hs){
                               
                               idx_hs <- my_interactions$anno_samples$health_status == hs
                               
                               # create a matrix with mean weights
                               mean_w_mat <- matrix(NA
                                                    ,nrow = length(cell_types)
                                                    ,ncol = length(cell_types))
                               rownames(mean_w_mat) <- cell_types
                               colnames(mean_w_mat) <- cell_types
                               
                               # populate the matrix with mean weights
                               for(i in cell_types){
                                       for(j in cell_types){
                                               idx_send <- my_interactions$anno_interactions$sending_cell_type == i
                                               idx_rec <- my_interactions$anno_interactions$receiving_cell_type == j
                                               
                                               # calculate the mean for each significant interaction between the cell types of interest in the cohort
                                               mean_weight <- rowMeans(my_interactions$weights[idx_sign & idx_send & idx_rec
                                                                                               ,idx_hs]
                                               )
                                               # calculate the mean of non-zero values
                                               my_mean <- mean(mean_weight[mean_weight != 0])
                                               
                                               # populate the matrix
                                               ifelse(is.na(my_mean)
                                                      ,mean_w_mat[i,j] <- 0.000000001 # if an edge is completely missing, the plotting function throws an error
                                                      ,mean_w_mat[i,j] <- my_mean)
                                               
                                       }
                               }
                               
                               # create a graph object
                               g <- igraph::graph(unlist(lapply(rownames(mean_w_mat)
                                                                , function(i){
                                                                        lapply(colnames(mean_w_mat)
                                                                               ,function(j){
                                                                                       c(i,j)
                                                                               })
                                                                        
                                                                })
                               ))
                               
                               # collapse the graph ogject to unique edges
                               g_simp <- igraph::simplify(g, remove.loops = F)
                               
                               
                               
                               
                               
                               
                               summary <- data.frame(edge_ID = {# extract edges from my_graph object
                                       edges <- as_edgelist(g, names = TRUE)
                                       edges <- as.data.frame(sapply(1:nrow(edges),function(i)paste(edges[i,])))
                                       
                                       ID <- sapply(1:ncol(edges)
                                                    ,function(j) paste(edges[1,j],edges[2,j],sep = " to "))
                                       ID
                               }
                               ,edge.width = {
                                       # define edge width as the desired magnification of the edge weight
                                       edge.width <- unlist(lapply(1:nrow(mean_w_mat),function(i) mean_w_mat[i,]))
                                       edge.width <- edge.width*amplify_edgeWidth
                                       edge.width
                               }
                               ,color_resolution = {
                                       # define resolution for the colors
                                       color_resolution <- round(edge.width*amplify_colorResolution
                                                                 ,digits = 0)
                               }
                               ) 
                               
                               # return the list
                               list(g_simp = g_simp
                                    ,mean_w_mat = mean_w_mat
                                    ,summary = summary)
                               
                       })
        names(data) <- health_status
        
        if(verbose){
                print(str(data))
                lapply(health_status
                       ,function(hs){
                               print(hs)
                               print(data[[hs]]$summary)
                       })
        }
        
        # plot
        for(hs in health_status){
                
                # define maximum number of colors if color_palette is used
                nr_colors <- max(unlist(sapply(health_status,function(h){data[[h]]$summary$color_resolution})))
                
                # add color scheme for the edges 
                if(is.null(edge.color)){
                        edge.col <- colorRampPalette(color_palette)
                        edge.color <- edge.col(nr_colors)[df$colorResolution+1]
                } else {
                        edge.color <- edge.color[data[[hs]]$summary$edge_ID] # make sure it is sorted correctly
                }
                
                # plot network
                plot(data[[hs]]$g_simp
                     ,edge.width = data[[hs]]$summary$edge.width 
                     ,edge.arrow.size=edge.arrow.size
                     ,edge.alpha=0.5
                     ,edge.curved=0.05
                     ,edge.color = edge.color
                     ,edge.attr.comb=c(weight="sum", type="ignore")
                     ,vertex.size=vertex.size
                     ,vertex.label.cex = vertex.label.cex
                     ,vertex.shape=vertex.shape
                     ,main = ""
                     ,layout = layout_in_circle(data[[hs]]$g_simp) 
                     ,...
                )
                title(hs, cex.main = title_cex)
        }
}

##### Forest plot functions

#' @title order_interactions_for_forests
#' 
#' @description This function takes in a vector of indices and subsets the interactions dataframe based on the indices. Then it sorts the interactions based on log2FC_weights, log2FC_p_s_l, log2FC_phi_s_l, log2FC_rho_s, log2FC_p_r_r, log2FC_phi_r_r, log2FC_rho_r and interaction_ID columns. It also has helper functions like cluster_interactions and pick_param which are used in the sorting process.
#' 
#' @param my_idx A vector of indices for subsetting the interactions dataframe
#' 
#' @return A sorted dataframe of interactions
#' 
#' @export
#' @examples
#' # load example data
#' data("comm_result")
#' # calculate general statistics
#' comm_result <- general_stat(comm_result)
#' # order interactions
#' order_interactions_for_forests(my_idx = 1:10)
#' 
order_interactions_for_forests <- function(my_anno_interactions
                                           ,threshold = 1
){
        order_by_direction <- function(values
                                       ,threshold = threshold
        ){ # reurns vector of numeric indices
                
                idx_greater <- values >= threshold
                idx_smaller <- values <= -threshold
                idx_within <- !(idx_greater | idx_smaller)
                
                num_idx_greater <- which(idx_greater)
                num_idx_smaller <- which(idx_smaller)
                num_idx_within <- which(idx_within)
                
                c(num_idx_greater
                  ,num_idx_smaller
                  ,num_idx_within)
        } 
        
        # sort by log2FC p_r_r
        my_anno_interactions <- my_anno_interactions[order_by_direction(my_anno_interactions$log2FC_p_r_r
                                                                        ,threshold = threshold
        )
        ,]
        
        # sort by log2FC p_s_l
        my_anno_interactions <- my_anno_interactions[order_by_direction(my_anno_interactions$log2FC_p_s_l
                                                                        ,threshold = threshold
        )
        ,]
        
        # sort by log2FC phi_r_r
        my_anno_interactions <- my_anno_interactions[order_by_direction(my_anno_interactions$log2FC_phi_r_r
                                                                        ,threshold = threshold
        )
        ,]
        
        # sort by log2FC phi_s_l
        my_anno_interactions <- my_anno_interactions[order_by_direction(my_anno_interactions$log2FC_phi_s_l
                                                                        ,threshold = threshold
        )
        ,]
        
        # sort by log2FC rho_r
        my_anno_interactions <- my_anno_interactions[order_by_direction(my_anno_interactions$log2FC_rho_r
                                                                        ,threshold = threshold
        )
        ,]
        
        # sort by log2FC rho_s
        my_anno_interactions <- my_anno_interactions[order_by_direction(my_anno_interactions$log2FC_rho_s
                                                                        ,threshold = threshold
        )
        ,]
        
        # sort by log2FC weights
        my_anno_interactions <- my_anno_interactions[order_by_direction(my_anno_interactions$log2FC_weights
                                                                        ,threshold = threshold
        )
        ,]
        
        
        # make as factor such that the plotting functino does not sort them alphabetically again
        my_anno_interactions$interaction_ID <- factor(my_anno_interactions$interaction_ID
                                                      ,levels = my_anno_interactions$interaction_ID
                                                      ,ordered = TRUE)
        
        # return result
        return(my_anno_interactions)
        
}

plot_cell_type_annotation <- function(my_df
                                      ,which_cell_type # "sending_cell_type" or "receiving cell type"
                                      ,title = NA
                                      ,color_cell_type = NA
){
        
        my_df$interaction_ID <- factor(my_df$interaction_ID
                                       ,levels=my_df$interaction_ID
                                       ,ordered = TRUE
        )
        
        my_df[[which_cell_type]] <- factor(my_df[[which_cell_type]]
                                           ,levels=unique(my_df[[which_cell_type]])
                                           ,ordered = TRUE
        )
        
        my_df$dummy_value <- -0.1
        my_df$cell_type_to_plot <- my_df[[which_cell_type]]
        
        p <- ggplot(my_df)+
                geom_bar(mapping = aes(x = interaction_ID
                                       , y = dummy_value
                                       , fill = cell_type_to_plot)
                         ,stat = "identity"
                         , width = 1
                )+
                ylim(c(-1,0.5))+
                theme_void()+
                coord_flip()+
                ggtitle(title)+
                theme(legend.position = 'bottom'
                      ,legend.title=element_blank()
                      ,plot.title = element_text(hjust = 0.5))
        
        if(class(color_cell_type) == 'logical'){
                p
        } else {
                p+
                        scale_fill_manual(values=color_cell_type)
                
        }
        
}

#' @title plot_forest
#'
#' @description creates a horizontal bar chart of log2FC expression values for a set of interactions. The function also allows to customize the color scale and the range of the log2FC values.
#'
#' @param my_df: a data frame containing the following columns: 'interaction_ID' and 'log2FC'. The interaction_ID column must be unique.
#' @param my_title: character: title for the plot
#' @param plot_legend: logical: should the legend be plotted (default is TRUE)
#' @param min: numeric: minimal log2FC value to be plotted (default is min of my_df$log2FC)
#' @param max: numeric: maximal log2FC value to be plotted (default is max of my_df$log2FC)
#'
#' @return ggplot object
#'
#' @export
#' @examples
#' # plot_forest
#' plot_forest(my_df
#'           ,my_title
#'          ,plot_legend = TRUE
#'         ,min
#'        ,max
#' )
#' 
plot_forest <- function(my_df
                        , my_title
                        , plot_legend = TRUE
                        , min
                        , max
                        ,legend_title_size
                        ,legend_text_size
) {
        
        my_values <- c(min, -1.1, -0.5, 0, 0.5, 1.1, max)
        my_colors <- c("lightslateblue", "lightslateblue", "aliceblue", "gray90", "lavenderblush", "red3",
                       "red3")
        names(my_colors) <- my_values
        
        idx_max <- max <= my_values
        
        if (sum(!idx_max) != length(my_values) - 1) {
                ifelse(min == max, {
                        my_colors <- my_colors[c(rep(TRUE, sum(!idx_max) + 2), rep(FALSE, sum(idx_max) -
                                                                                           2))]
                }, {
                        my_colors <- my_colors[c(rep(TRUE, sum(!idx_max) + 1), rep(FALSE, sum(idx_max) -
                                                                                           1))]
                })
                
        }
        my_length <- length(my_colors)
        
        idx_min <- min >= my_values[1:my_length]
        
        if (sum(!idx_min) != my_length - 1) {
                ifelse(min == max, {
                        my_colors <- my_colors[c(rep(FALSE, sum(idx_min) - 2), rep(TRUE, sum(!idx_min) +
                                                                                           2))]
                }, {
                        my_colors <- my_colors[c(rep(FALSE, sum(idx_min) - 1), rep(TRUE, sum(!idx_min) +
                                                                                           1))]
                })
                
        }
        
        resc_values <- rescale(as.numeric(names(my_colors)))
        
        my_p <- ggplot() + 
                geom_bar(data = my_df
                         , aes(x = interaction_ID
                               , y = log2FC
                               ,color = log2FC, fill = log2FC)
                         , stat = "identity"
                         , position = "identity") +
                scale_colour_gradientn(colours = my_colors
                                       , values = resc_values) + 
                scale_fill_gradientn(colours = my_colors
                                     ,values = resc_values) + 
                theme_void() + 
                geom_hline(yintercept = 0)+
                theme(axis.text.x = element_text(size = 0)
                      ,axis.text.y = element_blank()
                      ,axis.ticks.y = element_blank()
                      ,legend.position = "bottom"
                      ,plot.title = element_text(hjust = 0.5)
                      ,legend.title=element_text(size=legend_title_size)
                      ,legend.text=element_text(size=legend_text_size)
                      ,legend.box.spacing = unit(20, "pt")
                ) +
                xlab("") + 
                ylab("") + 
                ggtitle(my_title) + 
                coord_flip()
        
        ifelse(plot_legend, return(my_p), return(my_p + theme(legend.position = "none")))
}

#' @title plot_all_forests
#'
#' @description creates 7 forest plots for different parameters.
#'
#' @param my_idx Index of the interactions to be plotted
#' @param my_anno_interactions dataframe: Dataframe containing the interactions, having the rows as interactions and columns as the different parameters and the interaction_ID.
#' @param keep_order a logical value. If TRUE then keeps oder of my_anno_interactions ignoring the my_idx, default is FALSE
#' @param show_labels a logical value. If true interaction_IDs are plotted, default is FALSE
#' @param plot_legend a logical value. If true the legend is plotted, default is TRUE
#'
#' @return a list of 7 forest plots
#'
#' @export
#' @examples
#' # plot_all_forests
#' plot_all_forests(my_idx
#'                ,my_anno_interactions
#' )
#' 
plot_all_forests <- function(my_idx
                             ,my_anno_interactions
                             ,keep_order=FALSE
                             ,show_labels=FALSE
                             ,plot_legend=TRUE
                             ,color_cell_type = NA
                             ,annotate_interaction_types = FALSE
                             ,interaction_types = NA # named vector: names are interaction IDs, values are interaction classes
                             ,color_interaction_types = NA # named vector: names are interaction classes, values are colors
                             ,threshold = 1
                             ,legend_title_size = 10
                             ,legend_text_size = 9
){
        
        # define order
        if(!keep_order){
                my_anno_interactions <- order_interactions_for_forests(my_anno_interactions[my_idx,]
                                                                       ,threshold = threshold)  
        } else{
                my_anno_interactions$interaction_ID <- factor(my_anno_interactions$interaction_ID
                                                              ,levels = my_anno_interactions$interaction_ID
                                                              ,ordered = TRUE)
        }
        
        
        params <- c("log2FC_weights"
                    ,"log2FC_rho_s"
                    ,"log2FC_phi_s_l"
                    ,"log2FC_p_s_l"
                    ,"log2FC_rho_r"
                    ,"log2FC_phi_r_r"
                    ,"log2FC_p_r_r"
        )
        my_data <- lapply(params
                          ,function(i){
                                  
                                  test_df <- my_anno_interactions[,c(i,"interaction_ID")]
                                  colnames(test_df) <- c("log2FC","interaction_ID")
                                  
                                  test_df
                          }
        )
        
        names(my_data) <- params
        
        
        
        p_ct_s <- plot_cell_type_annotation(my_df = my_anno_interactions
                                            ,which_cell_type = "sending_cell_type"
                                            ,title = "cell_type_s"
                                            ,color_cell_type = color_cell_type)
        
        p_ct_r <- plot_cell_type_annotation(my_df = my_anno_interactions
                                            ,which_cell_type = "receiving_cell_type"
                                            ,title = "cell_type_r"
                                            ,color_cell_type = color_cell_type)
        
        
        p_w <- plot_forest(my_data$log2FC_weights
                           ,my_title = "w"
                           ,min = min(my_data$log2FC_weights$log2FC)
                           ,max = max(my_data$log2FC_weights$log2FC)
                           ,plot_legend = plot_legend
                           ,legend_title_size = legend_title_size
                           ,legend_text_size = legend_text_size
        )
        
        p_rho_s <- plot_forest(my_data$log2FC_rho_s
                               ,my_title = "rho_s"
                               ,min = min(my_data$log2FC_rho_s$log2FC)
                               ,max = max(my_data$log2FC_rho_s$log2FC)
                               ,plot_legend = plot_legend
                               ,legend_title_size = legend_title_size
                               ,legend_text_size = legend_text_size
        )
        
        p_phi_s_l <- plot_forest(my_data$log2FC_phi_s_l
                                 ,my_title = "phi_s_l"
                                 ,min = min(my_data$log2FC_phi_s_l$log2FC)
                                 ,max = max(my_data$log2FC_phi_s_l$log2FC)
                                 ,plot_legend = plot_legend
                                 ,legend_title_size = legend_title_size
                                 ,legend_text_size = legend_text_size
        )
        p_p_s_l <- plot_forest(my_data$log2FC_p_s_l
                               ,my_title = "p_s_l"
                               ,min = min(my_data$log2FC_p_s_l$log2FC)
                               ,max = max(my_data$log2FC_p_s_l$log2FC)
                               ,plot_legend = plot_legend
                               ,legend_title_size = legend_title_size
                               ,legend_text_size = legend_text_size
        )
        p_rho_r <- plot_forest(my_data$log2FC_rho_r
                               ,my_title = "rho_r"
                               ,min = min(my_data$log2FC_rho_r$log2FC)
                               ,max = max(my_data$log2FC_rho_r$log2FC)
                               ,plot_legend = plot_legend
                               ,legend_title_size = legend_title_size
                               ,legend_text_size = legend_text_size
        )
        
        
        p_phi_r_r <- plot_forest(my_data$log2FC_phi_r_r
                                 ,my_title = "phi_r_r"
                                 ,min = min(my_data$log2FC_phi_r_r$log2FC)
                                 ,max = max(my_data$log2FC_phi_r_r$log2FC)
                                 ,plot_legend = plot_legend
                                 ,legend_title_size = legend_title_size
                                 ,legend_text_size = legend_text_size
        )
        p_p_r_r <- plot_forest(my_data$log2FC_p_r_r
                               ,my_title = "p_r_r"
                               ,min = min(my_data$log2FC_p_r_r$log2FC)
                               ,max = max(my_data$log2FC_p_r_r$log2FC)
                               ,plot_legend = plot_legend
                               ,legend_title_size = legend_title_size
                               ,legend_text_size = legend_text_size
        )
        
        my_data$empty_values <- data.frame(log2FC = rep(0,(nrow(my_data[[1]])+1))
                                           ,interaction_ID = c("",as.character(my_data[[1]]$interaction_ID))
        )
        my_data$empty_values$interaction_ID <- factor(my_data$empty_values$interaction_ID
                                                      ,levels = c("",as.character(my_data[[1]]$interaction_ID))
                                                      ,ordered = TRUE)
        
        p_IDs <- ggplot(my_data$empty_values
                        ,aes(y = interaction_ID
                             ,x = log2FC
                        )
        )+ theme_classic()
        p_IDs <- p_IDs + theme(axis.line.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),
                               axis.title.x=element_blank(),
                               panel.grid.minor.x=element_blank(),
                               panel.grid.major.x=element_blank()
                               ,axis.title.y=element_blank()
                               ,axis.ticks.y=element_blank()
                               ,axis.line.y=element_blank()
        )
        
        
        
        margin = theme(plot.margin = unit(c(0
                                            ,-0.25
                                            ,0
                                            ,-0.25
        )
        , "cm")
        )
        
        if(annotate_interaction_types){
                df_classes <- cbind(my_anno_interactions
                                    ,interaction_type = interaction_types[as.character(my_anno_interactions$interaction_ID)])
                
                p_int_class <- plot_cell_type_annotation(my_df = df_classes
                                                         ,which_cell_type = "interaction_type"
                                                         ,title = "interaction_type"
                                                         ,color_cell_type = color_interaction_types)+guides(fill=guide_legend(nrow=2,byrow=TRUE))
        }
        
        p_empty <- ggplot()+ theme_void()+margin
        
        if(!show_labels){
                if(!annotate_interaction_types){
                        grid.arrange(p_w
                                     ,p_empty
                                     ,p_rho_s
                                     ,p_phi_s_l
                                     ,p_p_s_l
                                     ,p_empty
                                     ,p_rho_r
                                     ,p_phi_r_r
                                     ,p_p_r_r
                                     ,nrow = 1)
                }else{
                        grid.arrange(p_int_class
                                     ,p_w
                                     ,p_empty
                                     ,p_ct_s
                                     ,p_rho_s
                                     ,p_phi_s_l
                                     ,p_p_s_l
                                     ,p_empty
                                     ,p_ct_r
                                     ,p_rho_r
                                     ,p_phi_r_r
                                     ,p_p_r_r
                                     ,nrow = 1)
                }
        }else {
                if(!annotate_interaction_types){
                        
                        grid.arrange(p_int_class
                                     ,p_w
                                     ,p_empty
                                     ,p_ct_s
                                     ,p_rho_s
                                     ,p_phi_s_l
                                     ,p_p_s_l
                                     ,p_empty
                                     ,p_ct_r
                                     ,p_rho_r
                                     ,p_phi_r_r
                                     ,p_p_r_r
                                     ,p_IDs
                                     ,nrow = 1)
                }else {
                        grid.arrange(p_w
                                     ,p_empty
                                     ,p_ct_s
                                     ,p_rho_s
                                     ,p_phi_s_l
                                     ,p_p_s_l
                                     ,p_empty
                                     ,p_ct_r
                                     ,p_rho_r
                                     ,p_phi_r_r
                                     ,p_p_r_r
                                     ,p_IDs
                                     ,nrow = 1)
                }
        }
        
        
        
}

save(ls())