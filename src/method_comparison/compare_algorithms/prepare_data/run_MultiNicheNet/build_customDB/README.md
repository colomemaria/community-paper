# Overview

This section is divided into two subsections:

1. **NNetCustomBuild**: To generate the custom database for NicheNet, we need to construct 3 networks, ligand-receptor, signaling and gene regulatory networks. We use 'community' database as ligand receptor network and then build signaling and gene regulatory networks upon it. This notebook requires internet access. In this step, 3 files are produced ('lig_rec_sources', 'sig_Network' and 'gr_Network') and used in the 'Model Construction'.

2. **ModelConstruction_Custom**: Building the weighted matrices. 
