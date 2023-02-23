#!/bin/bash

# Create a directory to download
mkdir downloaded_files

# Download page containing GSE185381 information
wget https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185381 -qO- | 

# Get the lines that contain "(scRNA-seq)" and grab the link
grep -B2 "(scRNA-seq)" | grep -Eoi '<a [^>]+>' | 
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
        wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5613nnn/$ID/suppl/$line -P downloaded_files/
    done
done < barcodes.txt

gunzip -d downloaded_files/*.gz