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

## Running CircaN

To help you get an idea of the type of data CircaN uses as input, we have included a toy dataset with 200 features, 
along with it's metadata file in the package.

```
library("CircaN")
```
Load expression data and metadata and run circan function.
```
expression_example <- CircaN::expression_example
metadata_example <-CircaN::metadata_example

p <- circan(data=expression_example, s2c=metadata_example)
```


This will run CircaN algorithm on your data with default parameters. Depending on your analysis you may want to change
those to fit your needs. You can sepcify:

* data: Dataframe containing the expression data. Samples must be in columns and genes in rows. For an example see data(expression_example).
* meta: Dataframe containing the metadata for the samples. Must have at least a 'sample' column with the sample name as it appears in the data matrix; a 'time' column with the time point the sample was collected; and an 'ind' column containing information for the individual the sample comes from. For an example see data(metadata_example).
* shiny: Is the package running in a shiny app? default to FALSE.
* mode: Algorithm to use in the NLS regression. Must be one of 'default' for Gauss-Newton, 'plinear' for the Golub-Pereyra algorithm for partially linear least-squares models and 'port' for the ‘nl2sol’ algorithm from the Port library. Default is default. See nls documentation for extended info.
* init_value: Initial value for the period. Default is set to 24.
* max_per: Maximum period to regress. Default is set to Inf.
* min_per: Minimum period to regress. Default is set to -Inf.

## Workflow

CircaN will fit each feature in your data to 7 different oscilating patterns and the it will determine which better fits each feature (selecting the minimum AIC). Then it will combine the period and amplitude p.values using Fisher's method and finally will get the corresponding BH q.value. The selection criteria we recommend is using and adjusted p.value < 0.05 combined with a relatively high R (>0.7).
