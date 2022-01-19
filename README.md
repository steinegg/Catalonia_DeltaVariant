# Catalonia_DeltaVariant
The simulations can be run by executing the file 'runStanShort.R'. The necessary data is available in the folder 'data'. 
In the file 'stanModel.stan' the epidemiological model and the inference process are implemented. Parameters and data sources are defined in 'initializeValues.R'.
The file 'sequencing.R' analyzes the raw sequencing data to provide the necessary input for 'runStanShort.R'. However, its output is already present in the folder 'data'.
The stan file 'binomial_test.stan' implements a binomial test to find the credible interval for the % of cases that stem from the Delta variant. 

All the plots from the main analysis are created in 'plots.R'. The plots for the sensitivity analysis are in 'plotsSens.R'. 
