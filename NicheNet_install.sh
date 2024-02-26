#!/bin/bash

# get the the conda path
CONDA_PATH=$(dirname $(dirname $(which conda)))


mamba create --name NicheNet python=3.12 -y
source $CONDA_PATH/bin/activate NicheNet
mamba env update --file NicheNet_env.yml
R -e 'options(repos = c(CRAN = "https://cloud.r-project.org/")); install.packages("BiocManager"); BiocManager::install("metaboliteIDmapping"); q()'
R -e 'devtools::install_github("satijalab/seurat-object@v4.1.3"); devtools::install_github("satijalab/seurat@v4.3.0"); q()';
R -e 'devtools::install_github("saezlab/OmnipathR"); devtools::install_github("saeyslab/nichenetr"); devtools::install_github("SoloveyMaria/community", upgrade = "always"); q()';
R -e 'IRkernel::installspec(); q()';


conda deactivate

mamba create --name cellphonedb python=3.8 r-base rpy2 -y
source $CONDA_PATH/bin/activate cellphonedb
pip install cellphonedb==3.1.0

conda deactivate

