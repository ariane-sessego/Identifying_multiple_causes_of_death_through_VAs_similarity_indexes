This contains the script to compute the similarity indexes replicate the the results presented in  "Studying multiple causes of death through verbal autopsies: The contribution of a similarity index"

#### The scripts are in R version : 2024.04.04

#### Authors : anonymous for this review 

## This contains: 

- Computing__the_theoretical_index_of_similarity_InterVA4.R : script to compute the theoretical index from the probbase provided by InterVA and conducting robustness checks of the index using different norms or age-sex refinement. 
- Computing_the_individual_based_index_of_similarity.R : script to compute the empirical index from example detailed verbal autopsies provided.

Folders with input necessary for these computations:
- Input_probbase: The probability matrixes used by InterVA-4 to define causes of death, in letter and number format, with different HIV and malaria prevalence.
- Input_VA_example: Detailed_VA_examples.csv, simulated example of detailed verbal autopsies.
- probbaseVLVL_letters.csv : T

The rest of the documents are results/output created by the codes:
- Output_indexes: the computed indexes in different formats. 
 VLVL_Edgelist_indexes.csv: is the database with in line each association of causes and in column their corresponding indexes
 VLVL_Matrix_Eucl_index.csv, VLVL_Matrix_norm1_index.csv, VLVL_Matrix_normalise_index.csv: , VLVL_Matrix_scalar_product_index.csv: : is a matrix where each column and each line represent a cause. The value in each cell is the value of the index for the association of it’s line and column causes, each matrix representing a different index. 
Age- sex-specific Euclidean indexes and HIV and malaria specific indexes can also be found there.

- Output_figures: where output figures are stored, especially index heat maps.
