samples_dir=../../../../../results/method_comparison/compare_algorithms/CPDB/samples_DEGs/
my_vars=$(ls "$samples_dir" | cut -d_ -f1 | uniq)
#output_dir=../../../../../results/method_comparison/compare_algorithms/CPDB/
custom_db=../../../../../results/method_comparison/build_customDB/CPDB/custom_cellphone.db

for sample in $my_vars;
do
mkdir ${samples_dir}${sample}_results; 
cellphonedb method degs_analysis ${samples_dir}${sample}_meta.tsv ${samples_dir}${sample}_counts.tsv ${samples_dir}${sample}_DEGs.tsv --database $custom_db --counts-data hgnc_symbol --output-path ${samples_dir}${sample}_results/;
done;
