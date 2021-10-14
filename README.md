# feamiR: Classification and feature selection for microRNA/mRNA interactions #

[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/feamiR)](https://github.com/r-hub/cranlogs.app)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/feamiR)](https://github.com/r-hub/cranlogs.app)

**microRNAs** play a key role in RNA interference, the sequence-driven targeting of mRNAs that regulates their translation to proteins, through translation inhibition or the degradation of the mRNA. Around ~30% of animal genes may be tuned by microRNAs. The prediction of miRNA/mRNA interactions is hindered by the short length of the interaction (seed) region (~7-8nt). We collate several large datasets overviewing validated interactions and propose *fea*mi**R**, a novel pipeline comprising optimised classification approaches (Decision Trees/Random Forests and an efficient feature selection based on embryonic Genetic Algorithms used in conjunction with Support Vector Machines) aimed at identifying discriminative nucleotide features, on the seed, compensatory and flanking regions, that increase the prediction accuracy for interactions.

Common and specific combinations of features illustrate differences between reference organisms, validation techniques or tissue/cell localisation. *fea*mi**R** revealed new key positions that drive the miRNA/mRNA interactions, leading to novel questions on the mode-of-action of miRNAs.

*fea*mi**R** was presented at EMBL Symposium: Non-Coding Genome in October 2021, the poster is uploaded [here](https://github.com/Core-Bioinformatics/feamiR/blob/main/docs/figures/feamiR%20Poster.pdf) and is described in this video:

https://user-images.githubusercontent.com/33287710/137338503-53e23622-1a4a-4c92-a540-0699bb8c86a9.mp4




The *fea*mi**R** package provides two categories of functions: [a] Dataset pre-processing functions and [b] analysis functions.

<img src="https://github.com/Core-Bioinformatics/feamiR/blob/main/docs/figures/workflow.png" width="400">

*Overview of the *fea*mi**R** workflow.*


## Dataset preparation functions ##

The pre-processing function is called *preparedataset*. For generating the input mRNA (3'UTR) dataset two options are implemented: (1) one that takes as input a reference genome and an annotation file and (2) one that works directly with the 3' UTR sequences.

1. The reference genome is expected in fasta format; the corresponding annotation file is expected in gtf format (parameters **fullchromosomes** and **annotations** respectively). With these inputs, the 3'UTR sequences are generated, based on the genomic coordinates and used for the microRNA alignment. For example, for *H. sapiens*, these files can be retrieved from Ensembl [1] (Homo_sapiens.GRCh38.dna.toplevel.fa and Homo_sapiens.GRCh38.100.chr.gtf). Before using this option for generating the 3' UTRs, check whether the naming of chromosomes is consistent between the two files and the IDs from the interactions file you intend to use (e.g. "chr1" is different from "1"); also the number of chromosomes should be specified using the chr parameter (e.g. 23 for *H sapiens*).
2. The 3'UTR input is expected in fasta format; it can be specified using the **mRNA_3pUTR** parameter. In a similar approach as for (1) the consistency of IDs across files should be checked. 

The input miRNA file (**miRNA_full**) is expected in fasta format; it contains mature miRNA sequences (e.g. from [miRBase](http://www.mirbase.org/) [2]). Check the consistency of miRNA IDs across datasets (e.g. the interaction dataset). Using the mature sequences as input, the seed and compensatory regions will be extracted and stored in separate fasta files.

The mRNA (3'UTR) and miRNA datasets will be used to determine **all possible interactions** (using the PatMaN [4] alignment tool). 

The resulting set of interactions (PatMaN output) is then split into a positive dataset (**validated interactions** e.g. using the miRTarBase [3] database) and negative dataset (**non-validated interactions**); the interaction dataset must contain a 'miRNA' column and 'Target Gene' column; additional information includes an Experiments column detailing the biological experiment used to validate the interaction. If this column is supplied, additional analyses can be performed on classes with >**minvalidationentries** (default 40) unique values e.g. the statistical analyses for determining discriminative features can be performed on the dataset split by experiment type. Also a 'Support Type' column may be included; if miRTarBase is used as interactions dataset, the expected values are 'Functional MTI', 'Functional MTI (Weak)', 'Non-Functional MTI' and 'Non-Functional MTI (Weak)'; interactions likely to yield robust, reproducible results are 'Functional MTI' entries.

Post alignment, two types of analyses may be performed: (1) statistical analyses based on feature frequencies and (2) ML based feature selection and classification. For the former, the default comprises only the seed features; this set can be extended using the nonseed_miRNA and flankingmRNA parameters which focuses on full miRNA features and mRNA flanking features, respetively. The chi-squared and Fisher exact p-values, calculated on the whole set of positive and negative interactions, are stored in csv files; heatmaps are also generated and saved as jpgs.

<img src="https://github.com/Core-Bioinformatics/feamiR/blob/main/docs/figures/statistical_heatmap.png" width="500">

*Heatmaps illustrating chi-squared and Fisher exact *p*-values when assessing differences in positional proportions of single nucleotides (left) and dinucleotides (right) for each position in the positive and negative datasets, using a subset of *D. melanogaster* miRNAs and mRNAs.* 

This preprocessing step can be bypassed by supplying a positive and negative set using positiveset and negativeset parameters; the formatting of these user-defined sets has to match the *fea*mi**R** setup (column names and formatting).

A pre-processing step for the ML and feature selection component is the subsampling of the negative set to a number of entries comparable to the positive set. The default is set to 100 representative subsamples (checked by chi-squared tests). 

A prefix for all output files can be supplied using the o parameter.

## Machine Learning functions ##

Based on the subsamples created using the preparedataset function, using *fea*mi**R**  several miRNA/mRNA classifiers can be tested. In addition, functions for selecting discriminative features are also available.

### Classifiers ###

The classifier functions are: *decisiontree*, *randomforest* and *svm*. For the Random Forest (RF), the optimal number of trees can be selected as a hyperparameter (*selectrfnumtrees*); for the Support Vector Machines (SVMs) the optimal kernel can be selected as a hyperparameter (*selectsvmkernel*). Both optimisation functions rely on robustness on cross validation and will produce plots from which an appropriate number of trees and kernel can be identified. These functions should be tested on multiple subsamples to evaluate the robustness of the parameter selection.

<img src="https://github.com/Core-Bioinformatics/feamiR/blob/main/docs/figures/kernel_hyperparameter.png" width="400">

*Performance of SVM model with a variety of kernels on seed features using 10-fold cross-validation on a set of *H. sapiens* miRNA/mRNA interactions validated by Immunoprecipitation experiments. The distribution of training and test accuracy across cross-validation runs for each kernel is shown.*

<img src="https://github.com/Core-Bioinformatics/feamiR/blob/main/docs/figures/numtrees_hyperparameter.png" width="800">

*Performance of RF model with a range of number of trees, from 1 to 100, using 10-fold cross-validation on a set of *H. sapiens* miRNA/mRNA interactions validated by Immunoprecipitation experiments. The distribution of test and training accuracy across cross-validation runs for each number of trees is shown.*

Next the *runallmodels* function creates the results for Decision Trees, Random Forests and SVMs on all 100 subsamples. The selected hyperparameters (determined using *selectsvmkernel* and *selectrfnumtrees*) should be input as parameters. The function will output a data frame containing test and training accuracy, sensitivity and specificity for each model on each subsample; summary box plots are also produced. Additional output comprises **dtreevote** containing the selection of features from the decision trees, for each subsample, and their level within the tree, **ongoingginis** containing the Gini index for each feature based on the Random Forest model. These methods are described in more detail in the next section.

<img src="https://github.com/Core-Bioinformatics/feamiR/blob/main/docs/figures/performance.png" width="400">

*Performance of all models across 100 different subsamples of the full H. sapiens dataset. The plot shows the distribution of test and training accuracy.*

### Feature selection ###

The *fea*mi**R** package implements 6 feature selection approaches:
* *dtreevoting*: Uses trained Decision Trees (DTs) for feature selection, by maintaining a record of frequent features and the level of their first occurrence. Selection focuses only on features which appear in the top **num_levels** levels (from the root, default 10), as these are expected to have higher discriminative power. The selection of features is performed over **num_runs** runs (user-defined parameter with a default value of 100) on different subsamples output by *preparedataset*. The first column of the output table contains the number of runs for which each feature was used. The output table can then be ordered by frequency and highest-performing features selected.

<img src="https://github.com/Core-Bioinformatics/feamiR/blob/main/docs/figures/dtree_voting.png" width="800">

*Boxplot showing distribution, on DTs, over 100 subsamples of the full H. sapiens dataset, of depth from the root for each feature (considering only features appearing in the top 10 levels). The number of trees where each feature was used is shown above.*

* *rfgini*: Using Random Forest (RF) models, we assess variable importance with entropy-based measures. The importance of predictor variables is measured using ‘mean decrease in node impurity’, based on the Gini index. The function calculates the cumulative mean decrease in the Gini index across **num_runs** (default: 100) samples output by *preparedataset*. The first column of the output table contains the cumulative Gini index for each feature across all runs.

<img src="https://github.com/Core-Bioinformatics/feamiR/blob/main/docs/figures/rf_gini.png" width="800">

*Boxplot showing distribution of mean decrease in Gini index of features across 100  runs  on  different  subsamples of the full *H. sapiens* dataset, ordered  by  cumulative  mean  decrease  in Gini.*

* *forwardfeatureselection*: Forward Feature Selection uses a greedy approach to select the most discriminative **k** features. Features are ordered by discriminative power by selecting the feature which increases test accuracy at each iteration. The accuracy with respect to the selection order of features is plotted and smoothed using LOESS. When the curve describing the distribution of accuracies plateaus i.e. adding extra features does not improve the accuracy significantly, the selection process stops. The accuracy is assessed using a specified **model**, e.g. linear SVM. The function outputs an ordered list of features, along with accuracy, sensitivity and specificity achieved using these features.

<img src="https://github.com/Core-Bioinformatics/feamiR/blob/main/docs/figures/ffs.png" width="400">

*Improvement in performance when adding features, ordered on the x-axis, selected by FFS on a subset of the full H. sapiens dataset. The line type indicates the training (dotted line) and test (continuous line) accuracy, sensitivity and specificity for the first 20 features.*

* *geneticalgorithm*: Implements a standard genetic algorithm using GA package [5] with a fitness function specialised for feature selection.
* *eGA*: Feature selection based on Embryonic Genetic Algorithms. It performs feature selection by maintaining an ongoing set of 'good' features which are improved run by run. This is achieved by randomly selecting new features to combine with the ‘good’ features and performing *forwardfeatureselection*. It outputs training and test accuracy, sensitivity and specificity and a list of <=k features.

<img src="https://github.com/Core-Bioinformatics/feamiR/blob/main/docs/figures/eGA.png" width="800">

*Summary of features selected by eGA on full H. sapiens positive and negative sets. Distribution of ranks (top) and accuracies (bottom) across 100 subsamples for the training data). On top of the boxplots we present the number of selections per feature.*

It is recommended to use combinations of these feature selection approaches across multiple subsamples, corroborated with statistical analysis to select discriminative features.

The feature selected using these various methods can be summarized in heatmaps.

<img src="https://github.com/Core-Bioinformatics/feamiR/blob/main/docs/figures/feature_heatmap.png" width="500">

*Heatmap showing the discriminative power of features across the 5 feature-selection methods on the full H. sapiens positive set: Fisher exact p-values, DT voting scheme, RF cumulative mean decrease in Gini, FFS and eGA. The top 10 features across methods were included; all scores were quantile normalised on all 144 features per method, for comparability. The stronger features are shown in darker red.*

## Installation prerequisites ##

### R packages ###

The following R packages are required:
* stringr
* randomForest
* rpart
* rpart.plot
* GA
* e1071
* ggplot2
* magrittr
* tibble
* dplyr
* reticulate

To use parallelisation within the *geneticalgorithm* function, the following packages are required:
* parallel
* doParallel

### Python packages ###

The following Python packages must be installed to use *preparedataset*:
* Bio
* gtfparse
* pandas
* numpy
* math
* scipy.stats
* matplotlib.pyplot
* seaborn
* statistics
* logging
* os

To use these packages, Python >=3.4 must be installed. The path to the Python executable can be set using pythonversion parameter. If the correct Python version is recognised throughout the OS as ‘python’, there is no need to specify the path. 

### PatMaN ###

PatMaN [4] is a DNA pattern matcher for short sequences. It can be downloaded and installed from https://bioinf.eva.mpg.de/patman/. The path to the PatMaN executable can be set using the patmanpath parameter. If patman is recognised throughout the OS, there is no need to specify the path.

### sreformat ###

sreformat is a tool used for removing hidden characters and converting between sequence formats (DNA/RNA alphabet). It can be downloaded as part of SQUID from http://eddylab.org/software.html. The path to the sreformat executable can be set using the sreformatpath parameter. If sreformat is recognised throughout the OS, there is no need to specify the path.

## References ##
[1] A.D. Yates et al. Ensembl 2020. Nucleic Acids Research, 48(D1):D682–D688, 11 2019. 

[2] A. Kozomara, M. Birgaoanu, and S. Griffiths-Jones. miRBase: from microRNA sequences to function. Nucleic Acids Research, 47(D1):D155–D162, 2018.

[3] H.Y. Huang et al. miRTarBase 2020: updates to the experimentally validated microRNA–target interaction database. Nucleic Acids Research, 48(D1):D148–D154, 2019.

[4] K. Prufer et al. PatMaN: Rapid alignment of short sequences to large databases. Bioinformatics, 24:1530–1, 2008.

[5] L. Scrucca. GA: A package for genetic algorithms in R. Journal of Statistical Software, 53(4):1–37, 2013.
