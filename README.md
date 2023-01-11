# roplspvs: R orthogonal projections of latent structures with permutation over variable selection
This script uses the Bioconductor package ropls  to produce OPLS models. It performs variable selection using p(corr) and optionally VIP. In addition to permutation post variable selection included in the ropls package also permutations pre variable selection with proceeding variable selection of every permutation resulting in p-values for R2 and Q2 including the variable selection procedure. It produces tables of all group comparisons including optional stratification by metadata.
Input is one matrix, one metafile with sampleID including the groups to be compared and two Configure files with desired settings.
The output is one html file per comparison of groups containing models using five different model startegies and one summary html file and .Rdata files of all comparisons. Also .txt files with tables with loadings showing p(corr) are created.

# Installation
To be able to contribute to the repository install git by running the following code on command line. All other commands are in R.
 ```
 apt-get install git
 git clone https://github.com/pulmonomics-lab/roplspvs.git
 
```
To run the package without contributing download it from https://github.com/pulmonomics-lab/roplspvs.git

to load the package to get help files and load the function oplspvs run the following code (this is not necessary to run the workflow since all help is also in the configuration files)

```
require(tools)
devtools::load_all(<path to roplspvs>)

```

This package is tested on R 3.6.2, 4.0.0, 4.0.3, 4.0.5 and 4.1.2
Install needed packages BiocManager, Ropls, ggplot2, ggrepel, kableExtr, gridEstra, ggpubr, matrixStats, stringr, tryCatchLog, tools, devtools and DescTools. These are installed if running the file dependencies.R using the following code
 
```
source("<path to roplspvs>/dependencies.R")

 ```

# Run example data
To run example-data download folders and files in roplspvs. Set working directory to your project folder where the results will be stored. To start the analysis run the following code

```
setwd("<path to your project folder>")
source("<path to roplspvs>/Roplspvs_Run.R")

```

# Run your own data
## Preparations:

1.	Prepare the datamatrix with subjectID in the first column and variable names in the first row. 
2.	All names has to be unique. Do not use the following symbols in subjectID or variable names; ?, $,%, ^, &, *, (, ),-,#, ?,,,<,>, /, |, \, [ ,] ,{, and };
3.	Missing values should be indicated by "", "NA" or "Inf"
4.	The datamatrix may contain numeric data, integers or categorical data as characters. Original categorical variable will be removed and replaced with dummy variables.
5. A sampleID file with sampleID in the first column and containing one column with the groups to be compared and one column with secondaryID (for example gender) if stratification is desired. Also this file should not contain the symbols above and the samleID´s should agree with the subjectID´s in the datamatrix.
6. Both datamatrix and sampleID files should be saved as Tab delimited ("*.txt") files. 
7. Create a project folder for the analysis manually. In this folder a subdirectory will be created called "outputR" where the results and the tables called model_table_to_analyse and the reordered_model_table_to_analyse describing which models are run will be stored. 
8. Create a folder and save your datamatrix and sampleID file here. Default name is "data" and is created by default in the project folder. You may also enter the path and foldername to your own input data folder.  
9. Save the roplspvs_Configure_Get_Started.R and if advanced parameters are to be changed roplspvs_Configure_Advanced.R in your project folder. The names of the files may be changed to your choice "Any_name_Advanced.R" and "Any_name_Get_Started.R".
10. Edit basic settings of parameters in file "roplspvs_Configure_Get_Started.R" which contains the parameters that has to be entered including file names, folder names and column names describing the groups to be compared. Default settings and advanced parameters may be altered in file "roplspvs_Configure_Advanced.R".

## Run analysis:
When running roplspvs_Run.R the data in Configure files are loaded and the function oplspvs is run which sends sent to file "roplspvs_Models_of_each_comparison.Rmd" which renders one HTML file for each comparison containing score plots, loading plots, permutation pre, post and over variable selection plots and model statistics. When all comparisons have been completed "roplspvs_Summary_of_models.Rmd" renders a summary HTML-file of all models containing tables of all selected models and all significant models.
Set working directory to the location of your project folder where the results will be stored and run your analysis using the following code

```
setwd("<path to your project folder>"")
source("<path to roplspvs>/roplspvs_Run.R")

```

## Features

### Preprocessing
Preprocessing includes optionally replacing 0 with NA or lower limit of detection (LLD) and replacing NA with LLD. A userset value for LLD may be chosen. Data filtering is performed, allowing for a userset missing value tolerance for each group in each comparison. THe data is optionally logtransformed. Character data is replaced with dummy variables. If NA remains in the dataset NIPAL (Nonlinear Iterative Partial Least Squares) is used to impute values for the missing data, mean centering and Unit Variance scaling is performed using the ropls package.

### Variable selection
The features that contribute the most to the model are selected using p(corr) and optionally also using VIP. P(corr) of a variable is the Pearson correlation between the raw data and the scores of the model, i.e., a measure of how well a variable correlates with the model. The cutoff for p(corr) is set by the user to either a value or to correspond to a p-value for the Pearson correlation. VIP is a relative measure of how much the variable contributes to the model.

### Five model strategies for selecting variables
The script generates five models using different strategies for selection variables: 
Model strategy 0 show models pre variable selection
Model strategy 1 uses p(corr) cutoff which is userset either to numeric value or corresponding to userset p[Pearson_pcorr_cutoff] and userset number of orthogonals with default 0. 
Model strategy 2 uses p(corr) cutoff corresponding to userset p[Pearson_pcorr_cutoff] resulting in best performing model after variable selection. 
Model strategy 3 uses both p(corr) cutoff resulting in best performing model after variable selection. Only pcorr cutoff higher than p[Pearson_pcorr_cutoff] is used.
Model strategy 4 sets p(corr) cutoff and number of orthogonals creating in best performing model after variable selection while keeping the amount of variables to a minimum. They are limited by only adding variables if a decrease in p(corr) cutoff less than userset delta-p(corr) cutoff results in an increase in Q2 of more than 1%. The minimum p(corr) cutoff for model strategy 3 and 4 is set by variable p[Pearson_pcorr_cutoff]. 
Model strategy 5 is an iteration model using increasing p(corr) cutoff stepwise as long as Q2 of the model post variable selection is increased more than 1%. 

### Best performing model during selection of amount of orthogonal variables and p(corr) cutoff
Best performing models are defined as models that after variable selection give high Q2, low difference between R2 and Q2 and also low p[Q2_perm_ post_vs] and low p[R2_perm_ post_vs]. The weight between low difference and low p[Q2_perm_ post_vs] is given by user-set preferred_pR2_and_pQ2_permutated_post_vs with lower preferred_pR2_and_pQ2_permutated_post_vs giving more weight to p[Q2_perm_ post_vs] compared to high Q2 and low difference between R2 and Q2.

### Detailed description of method for selecting amount of orthogonal variables 
1)	Maximum number of orthogonals are user-set by variable max_no_of_ortho with default setting 5 for Model2-4 and by setting no_of_ortho with default 0 for Model1.
2)	The number of orthogonals is set using ropls default method adding orthogonals as long as Q2 is increased by 1 %.
3)  It is checked if there is a model with higher Q2 using fewer orthogonals using the following proceedure

  3.1)	Selects models with p[R2_perm_ post_vs] and p[Q2_perm_ post_vs]< preferred_pR2_and_pQ2_permutated_post_vs and diff between R2 and Q2 < 0.2
  
  3.2)	Selects models with max Q2 as long as adding an orthogonal increases Q2 more than 1%.
  
  3.3)	If no model is found p[R2_perm_ post_vs] and p[Q2_perm_ post_vs] limit is increased by preferred_pR2_and_pQ2_permutated_post_vs and diff between R2 and Q2 is increased by 0.1 and selection is rerun.
  
  3.4)	If there is more than one model with the same number of orthogonals and the same Q2, the model with lowest number of orthogonals after variable selection is chosen.

### Detailed description of method for selecting p(corr) cutoff
1)	From the start, all variables with |p(corr)|< p[Pearson_pcorr_cutoff] are removed.
2)	Model stratgy 2: pcorr cutoff set to correspond to p[Pearson_pcorr_cutoff] set by user.
3)	Model strategy 3 and 4: Selects p(corr) cutoff resulting in models with p[R2_perm_ post_vs] and p[Q2_perm_ post_vs] < preferred_pR2_and_pQ2_permutated_post_vs, difference between R2 and Q2 < 0.2 and max Q2
4)	Model strategy 3: Decreases Q2 by maximum 1% if fewer features or fewer orthogonals may be used.
5)	Model strategy 4: Decreases Q2 as long as increasing p(corr) cutoff more than user-set pcorr_diff increases Q2 more than 1% for each pcorr_diff step resulting in fewer variables.
6)	If no model is found p[R2_perm_ post_vs] and p[Q2_perm_ post_vs] limit is increased by preferred_pR2_and_pQ2_permutated_post_vs and diff between R2 and Q2 is increased by 0.1 searching again for model with highest Q2.

### SUS plots of models post variable selection using the SUSplot function
To compare two models from roplspvs, correlation plots of p(corr) for each variable in each model are plotted in a so-called shared and unique structure (SUS) plot, Wiklund et al in Anal Chem 2008. These are easy to create for models pre variable selection where all variables are present in both models. Models post variable selection are trickier to compare, as the variables selected in the models differ. Therefore, to compare models post variable selection in a SUS plot, the variables from both models were used to fit a new model. For the visualization in the SUS plot, p(corr) is selected firstly from the original model and secondly from the newly fitted model. The variables are colored for being shared or unique to each model. The names in the plot can be substituted for new names. Alternatively, the length and part of the name can be userset.

## Model statistics
### R2 and Q2
R2 is derived from the ropls package and is a measure of the fitness of the model and shows how much of the variation is explained by the model with a maximum value of 1. Q2 is a measure of the predictability of the model determined by 7 fold cross-validation which is the default in the ropls package. The package gives R2 and Q2 pre and post variable selection.

### Permutations
In addition to the permutations pre and post variable selection which is calculated by default by the ropls package resulting in  p[R2_perm_ sans_vs], p[Q2_perm_ sans_vs], p[R2_perm_ post_vs] and p[Q2_perm_ post_vs]  this package also calculates p-values for permutations over variable selection p[R2_perm_ over_vs] and p[Q2_perm_ over_vs]. Permutations over variable selection is performed by  randomization of subject labels followed by variable selection and fitting of model post variable selection to obtain R2 and Q2 for permutated models. R2 and Q2 for the permutated models are compared to R2 and Q2 for the unpermutated models to generate p[R2_perm_ over_vs] and p[Q2_perm_ over_vs]. The number of orthogonals is set by ropls default for permutated models as well as in the unpermutated model under investigation.  The number of permutations performed is user-set by setting parameter no_permutations_sans_vs, no_permutations_post_vs, no_permutations_post_vs_selected_models and no_permutations_over_vs with default 20.

### Establishing the size of overestimate in models
As suggested by Lindgren et al. in Journal of Chemometrics 1996, the overfit may be established for a specific dataset by the difference between the average Q2 of the permuted models over variable selection and the average Q2 of the permutated models pre variable selection. This difference establishes the average increase in Q2 during variable selection in the specific dataset. This difference in Q2 is a measure of the overfit of the model but may also contain an actual enhancement of the model by removing unrelated variance.

### Adjusting Q2 to not include the overestimate of the model
By removing the calculated overfit from the Q2 post variable selection, we create an adjusted Q2 which should represent a prediction without overfit. If this adjusted Q2 is negative, the model should be rejected and considered insignificant. A limit for how low it can be and still represent a model remains to be evaluated.

## RMSE
RMSE is derived using ropls package and is the square root of the mean error between actual and predicted responses.

