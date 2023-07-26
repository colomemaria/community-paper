#!/usr/bin/env python
# coding: utf-8

# # Part 4

# ## Correct batch effect

# Scgen batch correction.
# 
# - **INPUT:**
#     -  ```counts_norm.mtx``` 
#     -  ```anno_cells_norm.txt``` 
#     -  ```anno_samples_norm.txt``` 
#     -  ```anno_genes_norm.txt``` 
#     -  ```cell_relabelling.csv``` file containing unified cell type annotations. Stored in additional_input_files sub-directory.
#     
#     
# - **OUTPUT:**
# 
#     -  ```counts_corr.csv.gz``` 
#     -  ```anno_cells_corr.txt``` 
#     -  ```anno_samples_corr.txt``` 
#     -  ```anno_genes_corr.txt``` 
#     -  ```adata_afterCorrection.h5ad```

# ### load data

# In[1]:


import numpy as np
import anndata as ad
import scgen
import scanpy as sc
import pandas as pd
import numba as nb


# In[2]:


path_to_additional_files = "../../../data/data_preprocessing/vanGalen_Hourigan/additional_input_files/"
path_in = "../../../results/data_preprocessing/vanGalen_Hourigan/normalized/"
path_out = "../../../results/data_preprocessing/vanGalen_Hourigan/batch_corrected/"


# In[3]:


# read in count table and create an adata object
adata_beforeCorrection = sc.read_mtx(path_in + "counts_norm.mtx",'float64').T


# In[4]:


# read in cell annotation file
anno_cells = pd.read_csv(path_in + "anno_cells_norm.txt", sep = "\t") 


# In[5]:


# read in sample annotation file
anno_samples = pd.read_csv(path_in + "anno_samples_norm.txt", sep = "\t") 


# In[6]:


# read in gene annotation file
anno_genes = pd.read_csv(path_in + "anno_genes_norm.txt", sep = "\t") 


# In[7]:


# load cell subtype and color keys
cell_relabel=pd.read_csv(path_to_additional_files + "cell_relabelling.csv"
                        ,sep=';')


# ### process data

# In[8]:


adata_beforeCorrection.var_names=anno_genes['gene_symbol']
print("adata_beforeCorrection.var_names[1:10]")
print(adata_beforeCorrection.var_names[1:10])
adata_beforeCorrection.var = anno_genes
adata_beforeCorrection.var.index = adata_beforeCorrection.var.index.astype(str)


# In[9]:


adata_beforeCorrection.obs_names=anno_cells['cell_ID']
print("adata_beforeCorrection.obs_names[1:10]")
print(adata_beforeCorrection.obs_names[1:10])
adata_beforeCorrection.obs = anno_cells
adata_beforeCorrection.obs.index = adata_beforeCorrection.obs.index.astype(str)


# ### prepare for the visualization

# In[10]:


# define colors
colors_cell_type=dict(zip(cell_relabel["cell_type"],cell_relabel["cell_type_color_hex"]))
print(colors_cell_type)


# In[11]:


#there is an error at line 9, cell_subtype is only 2 values?
adata_beforeCorrection.obs["cell_subtype"].value_counts()


# In[12]:


# add cell_subtype column
cell_subtype = []
for t in adata_beforeCorrection.obs["cell_type_original"]:
    idx=cell_relabel["cell_type_original"]==t
    cell_subtype.append(cell_relabel["cell_subtype"][idx].values[0]) #
adata_beforeCorrection.obs["cell_subtype"] = pd.DataFrame(cell_subtype).values
# re-order such that the cell subtypes apper in the correct order in the legend
adata_beforeCorrection.obs["cell_subtype"]=adata_beforeCorrection.obs["cell_subtype"].astype("category")
adata_beforeCorrection.obs['cell_subtype'].cat.reorder_categories(['HSPC'
                                                  ,'early_Ery'
                                                  ,'late_Ery'
                                                  ,'pro_Mono'
                                                  ,'Mono'
                                                  ,'cDC'
                                                  ,'pDC'
                                                  ,'B'
                                                  ,'Plasma'
                                                  ,'CD4_T'
                                                  ,'CD8_T'], inplace=True)


# In[38]:


adata_beforeCorrection.obs["cell_subtype"]=adata_beforeCorrection.obs["cell_subtype"].astype("category")
adata_beforeCorrection.obs["cell_subtype"]=adata_beforeCorrection.obs['cell_subtype'].cat.reorder_categories(['HSPC'
                                                  ,'early_Ery'
                                                  ,'late_Ery'
                                                  ,'pro_Mono'
                                                  ,'Mono'
                                                  ,'cDC'
                                                  ,'pDC'
                                                  ,'B'
                                                  ,'Plasma'
                                                  ,'CD4_T'
                                                  ,'CD8_T'])


# In[39]:


# rename "malignangt" columns into "bares_mutatoin"
adata_beforeCorrection.obs["bares_mutation"] = adata_beforeCorrection.obs["malignant"].astype('str')


# In[40]:


# define color schemes
colors_cell_subtype=dict(zip(cell_relabel["cell_subtype"],cell_relabel["cell_subtype_color_hex"]))
print(colors_cell_subtype)

colors_cell_type=dict(zip(cell_relabel["cell_type"],cell_relabel["cell_type_color_hex"]))
print(colors_cell_type)


# In[41]:


adata_beforeCorrection.uns["health_status_colors"] = ["#7C001F" # bordeau for AML
                                                     , "#7ac5cd" # CadetBlue3 for healthy
                                                     ]

adata_beforeCorrection.uns["bares_mutation_colors"] = ["#A6ACAF" # grey for FALSE
                                                  , "#C0392B" # red for TRUE
                                                 ]

# ### visualize before batch correction

# In[43]:


# logtransform before HVG calculation
sc.pp.log1p(adata_beforeCorrection)

sc.pp.highly_variable_genes(adata_beforeCorrection)
sc.pl.highly_variable_genes(adata_beforeCorrection)


# In[44]:


# print how many HVGs we have:
print("Nr of HVGs:")
print(sum(adata_beforeCorrection.var.highly_variable))


# In[45]:


sc.tl.pca(adata_beforeCorrection, svd_solver='arpack')


# In[46]:


sc.pl.pca(adata_beforeCorrection
          , color=["cell_type"]
          ,palette=colors_cell_type
          ,save="_beforeCorrection_cell_type.pdf"
         )


# In[47]:


sc.pl.pca(adata_beforeCorrection
          , color=["dataset"]
          ,save="_beforeCorrection_cell_type_dataset.pdf"
         )


# In[48]:


sc.pp.neighbors(adata_beforeCorrection)


# In[49]:


sc.tl.umap(adata_beforeCorrection)



# ### batch correction

# In[57]:


# preprocess for batch correction
scgen.SCGEN.setup_anndata(adata_beforeCorrection, batch_key="dataset", labels_key="cell_type")


# In[ ]:


# create model
network = scgen.SCGEN(adata_beforeCorrection)


# In[ ]:


# train model
network.train(max_epochs=100,
            batch_size=32,
            early_stopping=True,
            early_stopping_patience=5,
)


# In[ ]:


# remove batch effect
adata_afterCorrection =  network.batch_removal()


# In[ ]:


# restore annotations
adata_afterCorrection.var = anno_genes
adata_afterCorrection.var.index = adata_afterCorrection.var.index.astype(str)
adata_afterCorrection.var_names=anno_genes['gene_symbol']
print("adata_afterCorrection.var_names[1:10]")
print(adata_afterCorrection.var_names[1:10])


# In[ ]:


adata_afterCorrection.obs = anno_cells
adata_afterCorrection.obs.index = adata_afterCorrection.obs.index.astype(str)
adata_afterCorrection.obs_names=anno_cells['cell_ID']
print("adata_afterCorrection.obs_names[1:10]")
print(adata_afterCorrection.obs_names[1:10])


# In[ ]:


print("adata_afterCorrection")
print(adata_afterCorrection)


# ### visualize after batch correction

# In[ ]:


sc.pp.highly_variable_genes(adata_afterCorrection)
sc.pl.highly_variable_genes(adata_afterCorrection)


# In[ ]:


print("Nr of HVGs:")
print(sum(adata_afterCorrection.var.highly_variable))


# In[ ]:


sc.tl.pca(adata_afterCorrection, svd_solver='arpack')


# In[ ]:


sc.pl.pca(adata_afterCorrection
          , color=["cell_type"]
          ,palette=colors_cell_type
          ,save="_afterCorrection_cell_type.pdf"
         )


# In[ ]:


sc.pl.pca(adata_afterCorrection
          , color=["dataset"]
          ,save="_afterCorrection_cell_type_dataset.pdf"
         )


# In[ ]:


sc.pp.neighbors(adata_afterCorrection)


# In[ ]:


sc.tl.umap(adata_afterCorrection)


# In[ ]:


adata_afterCorrection.uns["health_status_colors"] = ["#7C001F" # bordeau for AML
                                                     , "#7ac5cd" # CadetBlue3 for healthy
                                                     ]


# In[ ]:


adata_afterCorrection.uns["bares_mutation_colors"] = ["#A6ACAF" # grey for FALSE
                                                  , "#C0392B" # red for TRUE
                                                 ]


# In[ ]:


sc.pl.umap(adata_afterCorrection
           , color=["cell_type"]
          ,palette=colors_cell_type
          ,save="_afterCorrection_cell_type.pdf"
          )


# In[ ]:


sc.pl.umap(adata_afterCorrection
           , color=["cell_subtype"]
          ,palette=colors_cell_subtype
          ,save="_afterCorrection_cell_subtype.pdf"
          )


# In[ ]:


sc.pl.umap(adata_afterCorrection
           , color=["dataset"]
          ,save="_afterCorrection_dataset.pdf"
          )


# In[ ]:


sc.pl.umap(adata_afterCorrection
           , color=["patient_ID"]
          ,save="_afterCorrection_patient_ID.pdf"
          )


# In[ ]:


sc.pl.umap(adata_afterCorrection
           , color=["health_status"]
          ,save="_afterCorrection_health_status.pdf"
          )


# In[ ]:


sc.pl.umap(adata_afterCorrection
           , color=["bares_mutation"]
          ,save="_afterCorrection_malignant.pdf"
          )


# In[ ]:


sc.pl.umap(adata_afterCorrection
           , color=["cell_type_original"]
          ,save="_afterCorrection_cell_type_original.pdf"
          )


# # Export

# In[53]:


# export counts as csv.gz
print("save counts_corr.csv.gz")

counts_corr=pd.DataFrame(adata_afterCorrection.X
                        ,columns=adata_afterCorrection.var_names 
                        ,index=adata_afterCorrection.obs_names
                        ).transpose()
counts_corr.to_csv(path_out + "counts_corr.csv.gz"
                   ,index=True
                   ,compression="gzip"
                   )


# In[54]:


# export anno_cells_corr
print("save anno_cells_corr.txt")
adata_afterCorrection.obs.to_csv(path_out + "anno_cells_corr.txt"
                    ,sep = "\t"
                    ,index = True)


# In[55]:


# export anno_samples_corr
print("save anno_samples_corr.txt")
anno_samples.to_csv(path_out + "anno_samples_corr.txt"
                    ,sep = "\t"
                    ,index = True)


# In[56]:


# export anno_genes_corr
print("save anno_genes_corr.txt")
adata_afterCorrection.var.to_csv(path_out + "anno_genes_corr.txt"
                    ,sep = "\t"
                    ,index = True)


# In[57]:


# export adata object
print("save adata_afterCorrection.h5ad")
adata_afterCorrection.write(path_out + "adata_afterCorrection.h5ad")


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




