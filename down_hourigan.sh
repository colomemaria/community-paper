#!/bin/bash

# Create a directory to download
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116256/suppl/GSE116256_RAW.tar -P data/data_preprocessing/vanGalen_Hourigan/raw_data/

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120221/suppl/GSE120221_RAW.tar -P data/data_preprocessing/vanGalen_Hourigan/raw_data/

mkdir -p data/data_preprocessing/vanGalen_Hourigan/raw_data/GSE116256_RAW
mkdir -p data/data_preprocessing/vanGalen_Hourigan/raw_data/GSE120221_RAW

tar -xvf data/data_preprocessing/vanGalen_Hourigan/raw_data/GSE116256_RAW.tar --directory data/data_preprocessing/vanGalen_Hourigan/raw_data/GSE116256_RAW
tar -xvf data/data_preprocessing/vanGalen_Hourigan/raw_data/GSE120221_RAW.tar --directory data/data_preprocessing/vanGalen_Hourigan/raw_data/GSE120221_RAW

gunzip -d data/data_preprocessing/vanGalen_Hourigan/raw_data/GSE116256_RAW/*
gunzip -d data/data_preprocessing/vanGalen_Hourigan/raw_data/GSE120221_RAW/*
