In this section, we compare the default databases for CellPhoneDB, NicheNet, and our community database. 
We begin by comparing the L-R pairs in a Venn diagram. Since NicheNet has more than 10,000 predicted pairs, 
we also compare their database without predicted pairs. For our community database, we compare the curated pairs, 
while for CellPhoneDB, there are no predicted or curated pairs, so it is taken as is.

In the final section, we compare the gene space of these databases.

The comparison figures are saved under `ComparisonFigures/`, while the CSV files can be found in `output_csv/`.

CellPhoneDB employs SQLAlchemy and multiple tables to construct their database. To compare the databases, it is necessary to standardize the structure. To review the process we used to achieve this, please refer to the notebook in [restructure_CPDB](./restructure_CPDB/) directory, which includes comprehensive explanations.
