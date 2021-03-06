% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyperparameters.R
\name{selectsvmkernel}
\alias{selectsvmkernel}
\title{Tuning SVM kernel.
Trains SVMs with a range of kernels (linear, polynomial degree 2, 3 and 4, radial and sigmoid) using cross validation so the optimal kernel can be chosen (using the resulting plots). If specified (by showplots=F) the plots are saved as jpegs.}
\usage{
selectsvmkernel(data, title, showplots = T, output_prefix = "")
}
\arguments{
\item{data}{Dataset: dataframe containing classification column and all other column features. Both the training and test datasets will be taken from this dataset.}

\item{title}{Title to be used for the resulting boxplot}

\item{showplots}{T if plots should be shown in standard output, F is plots should be saved as jpg files.}

\item{output_prefix}{Prefix used for saving plots. If showplots==F then plots are saved here. Otherwise, standard output.}
}
\value{
Dataframe containing test and training accuracy, sensitivity and specificity
}
\description{
Tuning SVM kernel.
Trains SVMs with a range of kernels (linear, polynomial degree 2, 3 and 4, radial and sigmoid) using cross validation so the optimal kernel can be chosen (using the resulting plots). If specified (by showplots=F) the plots are saved as jpegs.
}
\examples{
data = read.csv(paste(system.file('samples/subsamples', package = "feamiR"),'/sample0.csv',sep=''))
data = rbind(head(data,50),tail(data,50))
data$classification = as.factor(data$classification)
data = data[,2:ncol(data)]
selectsvmkernel(data,'SVM boxplots')
}
\keyword{SVM}
\keyword{hyperparameter}
\keyword{kernel}
\keyword{tuning}
