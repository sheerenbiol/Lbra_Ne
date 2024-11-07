# Lbra_Ne
Analysis on effective population size of L. braziliensis populations across South America.
Analyses were part of the following pre-print:
Heeren, S., Sanders, M., Shaw, J. J., Brandão-Filho, S. P., Côrtes Boité, M., Cantanhêde, L. M., Chourabi, K., Maes, I., Llanos-Cuentas, A., Arevalo, J., Marco, J. D., Lemey, P., Cotton, J. A., Dujardin, J.-C., Cupolillo, E., & Van den Broeck, F. (2024). Evolutionary genomics of a zoonotic parasite across the Neotropical Realm. bioRxiv. https://doi.org/10.1101/2024.06.06.597691

Analyses were performed with:
- (G-PhoCS)[https://github.com/gphocs-dev/G-PhoCS]
- (MSMC2)[]


## G-PhoCS
R scripts are available for the downstream analyses based on the output from G-PhoCS.

Two scripts are available:
- *gphocs_Rfunctions.R* --> pretty self-explanatory: it contains all functions defined for data analysis and visualization
- *gphocs_Analysis.R* --> the script with the actual data analysis and visualization.

Data for the analyses is provided in the *GPhoCS_DATA* folder.
Scripts should be able to run as is.

G-PhoCS analyses were done in R version 4.0.5. Information on the Session and the packages can be found in: *gphocs_Rsession.txt*

## MSMC2


## Points to elaborate on:

- Data to run G-phocs is added.
- script and functions for downstream analyses in R are provided.

- Data is an subset of the entire data that is analysed
	- only used on 1 subset of individuals and two models (No migration and bi-directional)
	- Perhaps only a few or ALL chromosomes are added.
