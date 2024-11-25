#!/bin/bash

# usage ./download_raw_data.sh --lasry or ./download_raw_data.sh --hourigan


# Define a function to download and process data for a specific dataset
download_raw_data() {
    dataset="$1"
    mkdir -p "data/data_preprocessing/$dataset/raw_data"
    
    if [ "$dataset" == "Lasry" ]; then
        # Download page containing GSE185381 information
	wget https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185381 -qO- | 

	# Get the lines that contain "(scRNA-seq)" and grab the link
	grep -B2 "(CITE-seq)" | grep -Eoi '<a [^>]+>' | 
	grep -E 'href="/geo/query/acc.cgi' | 

	# Cut the string at the "=" sign and take the 3rd field
	cut -d\= -f3 | 

	# Remove the trailing characters after the first double quote
	sed 's/".*//' > barcodes.txt

	# Loop through the IDs in the barcodes.txt file
	while read ID; do 
	    # Download the directory containing the ID
	    wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5613nnn/$ID/suppl/ -qO- | 

	    # Get the lines that contain "processed" or "meta" and grab the link
	    grep -Eoi '<a [^>]+>' | grep -E 'processed|meta' | 

	    # Cut the string at the "=" sign and take the 2nd field
	    cut -d\= -f2 | 

	    # Remove the double quotes and ">" characters
	    tr -d '\"|>' | 

	    # Loop through the resulting file names and download them
	    while read -r line; do 
		wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5613nnn/$ID/suppl/$line -P data/Lasry/raw_data/
	    done
	done < barcodes.txt
	
	mkdir GSE185381_RAW
	
	for f in *.gz; do
	  STEM=$(basename "${f}" .gz)
	  gunzip -c "${f}" > ./GSE185381_RAW/"${STEM}"
	done	

	gunzip -d data/Lasry/raw_data/*.gz
	
	rm *.gz
	
	
	
    elif [ "$dataset" == "Hourigan" ]; then
        # Create a directory to download
	wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116256/suppl/GSE116256_RAW.tar -P data/vanGalen_Oetjen/raw_data/

	wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120221/suppl/GSE120221_RAW.tar -P data/vanGalen_Oetjen/raw_data/

	mkdir -p data/vanGalen_Oetjen/raw_data/GSE116256_RAW
	mkdir -p data/vanGalen_Oetjen/raw_data/GSE120221_RAW

	tar -xvf data/vanGalen_Oetjen/raw_data/GSE116256_RAW.tar --directory data/vanGalen_Oetjen/raw_data/GSE116256_RAW
	tar -xvf data/vanGalen_Oetjen/raw_data/GSE120221_RAW.tar --directory data/vanGalen_Oetjen/raw_data/GSE120221_RAW

	gunzip -d data/vanGalen_Oetjen/raw_data/GSE116256_RAW/*
	gunzip -d data/vanGalen_Oetjen/raw_data/GSE120221_RAW/*
	
	rm *.tar
	rm placeholder.yml
	
	
    else
        echo "Invalid dataset argument: $dataset"
        exit 1
    fi
}

# Check for the argument and call the function accordingly
if [ "$1" == "--lasry" ]; then
    download_raw_data "Lasry"
elif [ "$1" == "--vanGalen_Oetjen" ]; then
    download_raw_data "Hourigan"
else
    echo "Usage: $0 [--lasry | --vanGalen_Oetjen]"
    exit 1
fi


