# Settings
CONDA_ENV=community_paper_latest
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
	if ${CONDA} env list | grep ${CONDA_ENV}; then \
	   mamba env update -n ${CONDA_ENV} -f environment.yml; \
	   source ${ACTIVATE} ${CONDA_ENV} && R -e 'devtools::install_github("SoloveyMaria/community", upgrade = "always"); q()'; \
	else \
	    conda install -n base -c conda-forge mamba && \
	    source ${ACTIVATE} base && \
	    mamba env create -f environment.yml && \
	    source ${ACTIVATE} ${CONDA_ENV} && R -e 'options(repos = c(CRAN = "https://cloud.r-project.org/")); install.packages("BiocManager"); BiocManager::install("metaboliteIDmapping"); q()' && R -e 'devtools::install_github("saezlab/OmnipathR"); devtools::install_github("SoloveyMaria/community", upgrade = "always"); devtools::install_github("satijalab/seurat-object@v4.1.3"); devtools::install_github("satijalab/seurat@v4.3.0"); q()'; \
	    ./NicheNet_install.sh \
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
