IMPC Data Collation, Cleaning, Parameter Grouping | RShiny Dashboard | SQL Database

This project provides an R-based workflow for collating and cleaning IMPC (International Mouse Phenotype Consortium) data, grouping parameters to reduce parameter space, and generating normalized SQL tables with automated SQL insert statements. An RShiny dashboard is included for exploration of gene–phenotype associations and clustering.

Features
1. Data Collation & Preprocessing

•	Merges multiple IMPC measurement data files (uses datafileNames.txt on the HPC) into a consolidated dataset. 

•	Standardizes text encoding, corrects case-sensitive naming inconsistencies and typos, and trims whitespaces. 

•	Deduplicates and exports cleaned CSV files.

Run: DataCleaning_ParameterGroupingCode.qmd
________________________________________
2. Parameter Grouping

•	Integrates analysis data, parameter descriptions, and procedure information.

•	Produces a unified dataset for parameter grouping. 

•	Grouping categories include: Brain, Weight, Images, Clinical Chemistry, Heart, Eye, Embryonic Development, Immunophenotyping, Hematology, Grip, Open Field, Organ Weight, Other

•	Output saved as: parametergroupings_working.csv
________________________________________
3. SQL Database Normalization & Automated Inserts

•	Generates normalized SQL tables aligned to the database schema

•	Assigns primary keys and establishes parent–child relationships using R joins and match().

•	Automatically creates SQL insert statements for MySQL execution

•	Outputs stored in: /SQLdb_Tables_CSVs_Normalized/

Run: SQL_DatabaseCode.qmd

Database Dump: impc_grp9_sqldb_dumpfile.sql

Sample Queries: /sql_queries/
________________________________________
4. RShiny Visualization Dashboard

•	Explore p-values by parameter for a selected gene

•	Explore p-values by gene knockout for a selected phenotype

•	View clusters of genes with similar phenotype scores

•	Visualize phenotype scores for the Genes of Interest across parameters

Run: /ProjectAppTrial/app.R


Technologies Used:
R, RShiny, MySQL, HPC Environment
