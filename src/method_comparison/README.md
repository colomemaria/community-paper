# Database and Tool Comparison in Cell-Cell Communication Analysis

It is well recognized that communication tools and pipelines are constrained by the data available in their respective ligand-receptor databases. If there these databases contain bias towards a certain ligand-receptor spectrum, so too are results. Therefore, it's crucial to compare the default databases with unified custom databases, considering both ligand-receptor interactions confirmed in the literature and those predicted by computational tools. `compare_databases` section focuses on the comparative analysis of the databases.

## Database Structure and Customization
Comparing databases can be challenging due to the diverse structures used by different tools. For instance, `community` uses a single table,`CellPhoneDB` utilizes `SQLAlchemy` with multiple tables,  and `NicheNet` employs three different networks to link ligands to their target genes. The `community` tool is designed to simplify the creation of custom databases through user-friendly interfaces, allowing modifications of the default database as a csv/tsv file.

##### Dataset Used for Comparison
The comparison of databases is facilitated through the [`/compare_databases`](./method_comparison/compare_databases) section. This subsection details the datasets used and the process for comparing the default databases of `CellPhoneDB`, `NicheNet`, and `community`. It also provides instructions for generating a unified structure for each tool's database, aiming to standardize the comparison process.

##### Unified Structure for Comparing Results
To ensure a fair comparison between different communication tools, we standardize the output result from each tool. The pipeline under the [/compare_algorithms/](./method_comparison/compare_algorithms) directory is used for running each tool with the same input dataset and database. The databases are generated using the [/build_customDB](./method_comparison/compare_algorithms) section for each tool. The outputs are then transformed into a unified format, enabling a thorough comparison of the results. This approach ensures that differences in the output are attributable to the inherent capabilities of each tool, rather than to variations in output format or database.


# Subsampling and Measuring Resource Usage

As previously mentioned, our evaluation includes subsampling three datasets across a range from 6 to 32 samples. This method is designed to assess the adaptability of each tool to varying data sizes and complexities. In this section, we provide the subsampling scripts along with the resource measurements. Notably, all workflows were executed on a dedicated single node, equipped with 36 cores and 256GB of RAM.



























