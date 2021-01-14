# feamiR: Classification and feature selection for microRNA/mRNA interactions #

**microRNAs** play a key role in RNA interference, the sequence-driven targeting of mRNAs that regulates their translation to proteins, through translation inhibition or the degradation of the mRNA. Around ~30% of animal genes may be tuned by microRNAs. The prediction of miRNA/mRNA interactions is hindered by the short length of the interaction (seed) region (~7-8nt). We collate several large datasets overviewing validated interactions and propose *fea*mi**R**, a novel pipeline comprising optimised classification approaches (Decision Trees/Random Forests and an efficient feature selection based on embryonic Genetic Algorithms used in conjunction with Support Vector Machines) aimed at identifying discriminative nucleotide features, on the seed, compensatory and flanking regions, that increase the prediction accuracy for interactions.

Common and specific combinations of features illustrate differences between reference organisms, validation techniques or tissue/cell localisation. feamiR revealed new key positions that drive the miRNA/mRNA interactions, leading to novel questions on the mode-of-action of miRNAs.


The *fea*mi**R** package provides two categories of functions: [a] Dataset pre-processing functions and [b] analysis functions.

## Dataset preparation functions ##

The pre-processing function is called *preparedataset*. For generating the input mRNA (3'UTR) dataset two options are implemented: [1] one that takes as input a reference genome and an annotation file and [2] one that works directly with the 3' UTR sequences.

1. The reference genome is expected in fasta format; the corresponding annotation file  is expected in gtf format. Using these inputs, the 3'UTR sequences are generated, using the genomic coordinates and used for the microRNA alignment e.g. for *H. sapiens* these files can be retreived from Ensembl  Homo_sapiens.GRCh38.dna.toplevel.fa and Homo_sapiens.GRCh38.100.chr.gtf. Before using this option for generating the 3' UTRs, check whether the naming of chromosomes is consistent between the two files and the IDs from the interactions file you intend to use (e.g. "chr1" is different from "1"); also the number of chromosomes should be specified using the chr parameter (e.g. 23 for *H sapiens*).
2. The 3'UTR input is expected in fasta format; it can be specified suing the mRNA_3pUTR parameter. In a similar approach as for [1] the consistency of IDs across files should be checked. 

The **input miRNA file** is expeted in fasta format; it contains mature miRNA sequences (e.g. from miRBase http://www.mirbase.org/). Check the consisteny of miRNA IDs across datasets (e.g. the interaction dataset). Using the mature sequences as input, the seed and compensatory regions, will be extracted and stored in separate fasta files.

The mRNA (3'UTR) and miRNA datasets will be used to determine **all possible interactions** (using the PaTMaN alignment tool). For reformating and PaTMaN alignment the sreformat and PaTMaN tools are instalation pre-requisites; the paths to the executables can be set using the sreformatpath and patmanpath hyperparameters. If sreformat and patman are recognised throughout the OS, there is no need to specify the path. 

The resulting set of interactions (PaTMaN output) is then split into a positive dataset (**validated interactions** e.g. using the miRTarBase database) and negative dataset (**non-validated interactions**); the interaction dataset must contain a 'miRNA' column and 'Target Gene' column; additioanl information includes an Experiments column detailing the biological experiment used to validate the interaction. If this column is supplied, additional analyses can be performed on classes with >10 unique values e.g. the statistical analyses for determining discriminative features can be performed on the dataset split by experiment type. Also a 'Support Type' column may be included; if miRTarBase is used as interactions dataset, the expected values are 'Functional MTI', 'Functional MTI (Weak)', 'Non-Functional MTI' and 'Non-Functional MTI (Weak)'; interactions likely to yield robust, reproducible results are 'Functional MTI' entries.

Post alignment, two types of analyses may be performed: [1] statistical analyses based on feature frequencies and [2] ML based feature selection and classification. For the former, the default comprises only the seed features; this set can be extended using the nonseed_miRNA and flankingmRNA parameters which focuses on full miRNA features and mRNA flanking features. The chi-squared and Fisher exact p-values, calculated on the whole set of positive and negative interactions, are stored in csv files; heatmaps are also generated and saved as jpgs. 

A pre-processing step for the ML and feature selection component is the subsampling of the negative set to a number of entries comparable to the positive set. The default is set to 100 representative subsamples (checked by chi-squared tests). This pre-processing step can be bypassed by supplying a positive and negative set using positiveset and negativeset parameters; the formating of these user-defined sets has to match the feamir setup (column names and formating).

A prefix for all output files can be supplied using the -o parameter.

PLEASE NOTE: an installation prerequisite is Python (>=3.4). The following Python libraries must also be installed: os, Bio, gtfparse, pandas, numpy, math, scipy.stats, matplotlib.pyplot, seaborn as sns, statistics, logging.

## ML and feature selection functions ##

Based on the subsamples created using the preparedataset function, using *fea*mi**R**  several miRNA/mRNA classifiers can be used. In addition, functions for selecting discriminative features are also available.

The classifier functions are: decisiontree, randomforest and svm. For the RF the optimal number of trees can be selected as a hyperparameter (selectrfnumtrees); for the SVM the optimal kernel can be selected as a hyperparameter (selectsvmkernel). Both optimisation functions are base on cross validation and will produce plots from which an appropriate number of trees and kernel can be identified. These function should be tested on multiple subsamples to evaluate the robustness of the parameter selection.

Next the runallmodels function creates the results for Decision Trees, Random Forests and SVMs on all 100 subsamples. The selected hyperparameters (determined using selectsvmkernel and selectrfnumtrees) should be input as parameters. The function will output a data frame containing test and training accuracy, sensitivity and specificity for each model on each subsample; summary boxplots are also be produced. Additional output comprises **dtreevote** containing the selection of features from the decision trees, for each subsample, and their level within the tree, **ongoingginis** containing the Gini index for each feature based on the Random Forest model. The first column of dtreevote contains the number of runs for which each feature was used. The first column of ongoingginis contains the cumulative Gini index for each feature across the 100 runs.

As well as using the Decision Tree voting scheme and Random Forest cumulative Gini index measure, *fea*mi**R** implements three further feature selection approaches: traditional forwardfeatureselection, a genetic algorithm approach and a novel approach based on embryonic Genetic Algorithms (accessed through the eGA function).
It is recommended to use combinations of these feature selection appraoches across multiple subsamples, corroborated with statistical analysis to select discriminative features; summaries across all methods are presented using summary heatmaps.
