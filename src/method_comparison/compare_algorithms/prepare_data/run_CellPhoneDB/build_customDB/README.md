# Overview

CellPhoneDB stores ligand-receptor interactions as well as other properties of the interacting partners, including their subunit architecture and gene and protein identifiers. In order to create the content of the database, four main .csv data files are required: “gene_input.csv”, “protein_input.csv”, “ complex_input.csv” and “interaction_input.csv”.

In the 'build_customDB' notebook, we first generate the necessary files required by CellPhoneDB. Following this, the database is built using the provided command below. Detailed information can be found within the notebook. The custom database files needed for CellPhoneDB will be generated under the 'CPDB_Custom' directory as 'CPDB_Custom.db'

'''bash
cellphonedb database generate --user-interactions interactions.csv --user-interactions-only --user-protein prot_user.csv --user-gene gene_user.csv --result-path CPDB_Custom
'''
