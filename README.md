# CircaN

Nonlinear least squares model for accurate detection of circadian gene expression.

## Setup
If you don't already have, download and install the devtools package.

```
install.packages("devtools") 
```
Then, install CircaN from github directly to R.
```
library("devtools")
install_github("AndreaRP/CircaN")
```

## Running Full mode analysis

This will run the specified algorithms for circadian detection from JTK, MetaCycle meta2d and CircaN.

To help you get an idea of the type of data CircaN uses as input, we have included a toy dataset with 200 features, along with it's metadata file in the package. Please be mindful of the format of the data array, we have included the simplest way to convert a typical expression data into data CircaN can analyze by deleting the feature name column.

```
library("CircaN")
```
Load expression data and metadata and run circan function.

```
expression_example <- CircaN::expression_example
rownames(expression_example) <- expression_example$feature
expression_example <- expression_example[,-1]
metadata_example <-CircaN::metadata_example

results <- full_mode_analysis(data=expression_example, s2c=metadata_example)
```

This will run all three algorithms on your data with default parameters. Depending on your analysis you may want to change those to fit your needs. You can sepcify:

* data: Dataframe containing the expression data. Please observe that the name of the features must be specified as rownames, all columns in the dataframe must be numeric.
* s2c: Dataframe containing the metadata for the samples. Must have at least a 'sample' column with the sample name as it appears in the data matrix; a 'time' column with the time point the sample was collected; and an 'ind' column containing information for the individual the sample comes from. For an example see data(metadata_example).
* algorithms: Character vector with at least one of "circan", "jtk" or "metacycle", indicating the algorithms to use for the analysis. 
* circan_mode: Algorithm to use in the NLS regression. Must be one of 'default' for Gauss-Newton, 'plinear' for the Golub-Pereyra algorithm for partially linear least-squares models and 'port' for the ‘nl2sol’ algorithm from the Port library. Default is "port". See nls documentation for extended info.
* circan_init_value: Initial value for the period. Default is set to 24.
* max_per: Maximum period to regress. Default is set to 28.
* min_per: Minimum period to regress. Default is set to 20.
* mc_cycMethod: algorithms to test in the MetaCycle meta2d.

## Workflow

This will run each separate algorithm on your data, and then merge the results into a single array for convenience.
CircaN will fit each feature in your data to 7 different oscilating patterns and the it will determine which better fits each feature (selecting the minimum AIC). Then it will combine the period and amplitude p.values using Fisher's method and finally will get the corresponding BH q.value. The selection criteria we recommend is using and adjusted p.value < 0.05 combined with a relatively high R (>0.7).
For details on JTK (DOI: 10.1177/0748730410379711) and MetaCycle (DOI: 10.1093/bioinformatics/btw405), please refer to:
JTK: https://openwetware.org/wiki/HughesLab:JTK_Cycle    
MetaCycle: https://cran.r-project.org/web/packages/MetaCycle/vignettes/implementation.html    
