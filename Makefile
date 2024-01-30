# Settings
CONDA_ENV=community_paper
SHELL=bash
MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

ifeq ($(OS),Windows_NT)
    CONDA := $(strip $(shell where.exe conda))
else
    CONDA := $(strip $(shell which conda))
endif

ifeq ($(CONDA),"")
    BASE := ${HOME}/miniconda3/bin/conda
else
    BASE := $(shell dirname $(shell dirname ${CONDA}))
endif

ACTIVATE=${BASE}/bin/activate


default: help

install-conda: ## install Miniconda
	curl -L $(MINICONDA_URL) -o miniconda.sh
	bash miniconda.sh -b
.PHONY: install-conda


create-env: ## create conda environment
	if ${CONDA} env list | grep -w ${CONDA_ENV}; then \
	   echo "Updating existing conda environment"; \
	   mamba env update -n ${CONDA_ENV} -f environment.yml; \
	   echo "Activating environment and installing community package from GitHub"; \
	   source ${ACTIVATE} ${CONDA_ENV} && R -e 'devtools::install_github("SoloveyMaria/community", upgrade = "always"); q()'; \
	else \
	    echo "Installing Mamba in base environment"; \
	    conda install -n base -c conda-forge mamba && \
	    echo "Activating base environment"; \
	    source ${ACTIVATE} base && \
	    echo "Creating new conda environment from environment.yml"; \
	    mamba env create -f environment.yml && \
	    echo "Activating new environment and installing R packages"; \
	    source ${ACTIVATE} ${CONDA_ENV} && R -e 'options(repos = c(CRAN = "https://cloud.r-project.org/")); install.packages("BiocManager"); BiocManager::install("metaboliteIDmapping"); q()' && \
	    echo "Installing OmnipathR from GitHub"; \
	    R -e 'install.packages("prettyunits"); devtools::install_github("saezlab/OmnipathR@v3.7.0"); q()' && \
	    echo "Installing community package from GitHub"; \
	    R -e 'devtools::install_github("SoloveyMaria/community", upgrade = "always"); q()' && \
	    echo "Installing Seurat object and Seurat from GitHub"; \
	    R -e 'devtools::install_github("satijalab/seurat-object@v4.1.3"); devtools::install_github("satijalab/seurat@v4.3.0"); q()'&& \
	    echo "Deactivating conda environment"; \
	    conda deactivate && \
	    echo "Running NicheNet install script"; \
	    ./NicheNet_install.sh; \
	fi
.PHONY: create-env


download-data: ## download preprocessed data
	./download_raw_data.sh --lasry && ./download_raw_data.sh --hourigan
.PHONY: download-data


run-jupyter: ## run jupyter notebooks
	 source ${ACTIVATE} ${CONDA_ENV} && \
		jupyter notebook

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help
