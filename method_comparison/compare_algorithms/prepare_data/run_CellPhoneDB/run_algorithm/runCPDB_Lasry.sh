#conda activate cpdb2
# Set the directory path to the directory containing the DEG samples
samples_dir=outs/samples_DEGs/

# Get a list of unique variable names (i.e., sample names) by extracting the first part of the file names in the directory and removing duplicates
my_vars=$(ls "$samples_dir" | cut -d_ -f1 | uniq)

# Set the path to the custom database file
custom_db=../build_customDB/CPDB_Custom/CPDB_Custom.db

# Loop over each sample variable name
for sample in $my_vars;
do
  # Create a subdirectory for the sample results
  mkdir ${samples_dir}${sample}_results;

  # Run CellPhoneDB's DEG analysis method on the sample using the custom database, with input files in the sample directory and output files in the sample results subdirectory
  cellphonedb method degs_analysis --threshold 0.1 ${samples_dir}${sample}_meta.tsv ${samples_dir}${sample}_counts.tsv ${samples_dir}${sample}_DEGs.tsv --database $custom_db --counts-data hgnc_symbol --output-path ${samples_dir}${sample}_results/;
done;

