#' @importFrom magrittr %>%
NULL

classifier.performance <- function(classifier,data,label,includeplot=TRUE){
  predicted.label <- stats::predict(classifier,data, type = 'class')
  true.label = data[, 1]
  pred.table <- table(true.label, predicted.label)
  if (nrow(pred.table)==1){
    if (rownames(pred.table) %in% c(0,'0')){
      print(1)
      pred.table=data.frame('pred_0'=c(pred.table[1,1],0),'pred_1'=c(pred.table[1,2],0))
      rownames(pred.table)=c('true_0','true_1')
    }
    if (rownames(pred.table) %in% c(1,'1')){
      print(2)
      pred.table=data.frame('pred_0'=c(0,pred.table[1,1]),'pred_1'=c(0,pred.table[1,2]))
      rownames(pred.table)=c('true_0','true_1')
    }
  }
  if (ncol(pred.table)==1){
    if (colnames(pred.table) %in% c(0,'0')){
      pred.table=data.frame('pred_0'=c(pred.table[1,1],pred.table[2,1]),'pred_1'=c(0,0))
      rownames(pred.table)=c('true_0','true_1')
    }
    if (colnames(pred.table) %in% c(1,'1')){
      print(2)
      pred.table=data.frame('pred_0'=c(0,0),'pred_1'=c(pred.table[1,1],pred.table[2,1]))
      rownames(pred.table)=c('true_0','true_1')
    }
  }
  colnames(pred.table)=c('pred_0','pred_1')
  rownames(pred.table)=c('true_0','true_1')
  accuracy <- sum(diag(pred.table)) / sum(pred.table)
  sensitivity = pred.table['true_1','pred_1']/(pred.table['true_1','pred_0']+pred.table['true_1','pred_1'])
  specificity = pred.table['true_0','pred_0']/(pred.table['true_0','pred_0']+pred.table['true_0','pred_1'])
  return_list = list('accuracy'=accuracy,'sensitivity'=sensitivity,'specificity'=specificity)
  
  
  results = tibble::tibble(y_real  = true.label %>% factor(levels=c(0,1)),
                           y_pred  = predicted.label %>% factor(levels=c(0,1)),
                           Correct = ifelse(y_real == y_pred,"yes","no") %>% factor(levels=c('yes','no')))
  if (includeplot==TRUE){
    title = paste0(label)
    xlab  = 'True label'
    ylab  = 'Predicted label'
    graphics::plot(ggplot2::ggplot(results,ggplot2::aes(x = results$y_pred, y = results$y_real, colour = results$Correct)) +
                     ggplot2::geom_point() +
                     ggplot2::ggtitle(label = title, subtitle = paste0("Accuracy = ", 100*round(accuracy,3),"%")) +
                     ggplot2::xlab(xlab) +
                     ggplot2::ylab(ylab) +
                     ggplot2::scale_color_manual(labels = c('Yes', 'No'),
                                                 values = c('cornflowerblue','tomato')) +
                     ggplot2::geom_jitter() +
                     ggplot2::theme_bw()+
                     ggplot2::labs(colour = 'Correct'))
  }
  return(return_list)
}

#' Decision tree
#' Trains a decision on the given training dataset and uses it to predict classification for test dataset. The resulting accuracy, sensitivity and specificity are returned, as well as a tree summary.
#' @param data_train Training set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model is trained.
#' @param data_test Test set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model in tested.
#' @param includeplot Show performance scatter plot (default:FALSE)
#' @param showtree Show trained decision tree graphically (default:FALSE)
#' @return List containing performance summary, accessed using training (training accuracy), test (test accuracy), trainsensitivity, testsensitivity, trainspecificity, testspecificity. Also accessed using fit is the trained model produced. This can be used to find the features which appear at each level of the tree.
#' @keywords decision tree
#' @keywords dtree
#' @keywords classification
#' @export
#' @examples
#' data_train = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,0,0,1,1)),
#'       A=c(1,1,1,0,0,0,1,1,1,0),
#'       B=c(0,1,1,0,1,1,0,1,1,0),
#'       C=c(0,0,1,0,0,1,0,0,1,0))
#' data_test = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,1,0)),
#'       A=c(0,0,0,1,0,0,0,1),
#'       B=c(1,1,1,0,0,1,1,1),
#'       C=c(0,0,1,1,0,0,1,1))
#' decisiontree(data_train,data_test)

decisiontree<-function(data_train,data_test,includeplot=FALSE,showtree=FALSE){
  fit <- rpart::rpart(stats::as.formula("classification~."), data = data_train, method = 'class',minsplit=6)
  if (showtree==TRUE){print(rpart.plot::prp(fit))}
  training = suppressWarnings(classifier.performance(fit,data_train,label = 'Decision tree - training performance',includeplot = includeplot))
  test = suppressWarnings(classifier.performance(fit,data_test,label = 'Decision tree - test performance',includeplot = includeplot))
  return_list <- list("training"=training$accuracy,"test" = test$accuracy,"testsensitivity"=test$sensitivity,"testspecificity"=test$specificity,"trainsensitivity"=training$sensitivity,"trainspecificity"=training$specificity,"fit"=fit$frame)
  return(return_list)
}

#' Random Forest.
#' Trains a random forest on the training dataset and uses it to predict the classification of the test dataset. The resulting accuracy, sensitivity and specificity are returned, as well as a summary of the importance of features in the dataset.
#' @param data_train Training set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model is trained.
#' @param data_test Test set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model in tested.
#' @param numoftrees Number of trees used in the random forest (default:10)
#' @param includeplot Show performance scatter plot (default:FALSE)
#' @return List containing performance percentages, accessed using training (training accuracy), test (test accuracy), trainsensitivity, testsensitivity, trainspecificity, testspecificity. Also accessed using importance is the vector of Mean Decrease in Gini Index. This can be used to find the features which contribute most to classification.
#' @keywords rf
#' @keywords classification
#' @export
#' @examples
#' data_train = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,0,0,1,1)),
#'       A=c(1,1,1,0,0,0,1,1,1,0),
#'       B=c(0,1,1,0,1,1,0,1,1,0),
#'       C=c(0,0,1,0,0,1,0,0,1,0))
#' data_test = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,1,0)),
#'       A=c(0,0,0,1,0,0,0,1),
#'       B=c(1,1,1,0,0,1,1,1),
#'       C=c(0,0,1,1,0,0,1,1))
#' randomforest(data_train,data_test,numoftrees=5)

randomforest<-function(data_train,data_test,numoftrees=10,includeplot=FALSE){
  training.rf <- randomForest::randomForest(stats::as.formula("classification~."),data=data_train,ntree=numoftrees,proximity=TRUE,importance=TRUE,replace=TRUE)
  training = suppressWarnings(classifier.performance(training.rf,data_train,label = 'Random forest - training performance',includeplot = includeplot))
  test = suppressWarnings(classifier.performance(training.rf,data_test,label = 'Random forest - test performance',includeplot = includeplot))
  return_list <- list("training"=training$accuracy,"test" = test$accuracy,"testsensitivity"=test$sensitivity,"testspecificity"=test$specificity,"trainsensitivity"=training$sensitivity,"trainspecificity"=training$specificity,"importance"=randomForest::importance(training.rf))
  return(return_list)
}

#' SVM
#  Trains an SVM with the specified kernel (or default linear) on the supplied training dataset. The resulting accuracy, sensitivity and specificity are returned.
#' @param data_train Training set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model is trained.
#' @param data_test Test set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model in tested.
#' @param kernel Type of kernel to use for SVM model (default:linear)
#' @param degree Degree for kernel used (in polynomial or radial case)
#' @param poly Binary parameter stating whether the chosen kernel is polynomial of degree greater than 1 (default:0)
#' @param includeplot Show performance scatter plot (default:FALSE)
#' @return List containing performance percentages, accessed using training (training accuracy), test (test accuracy), trainsensitivity, testsensitivity, trainspecificity, testspecificity.
#' @keywords svm
#' @keywords classification
#' @export
#' @examples
#' data_train = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,0,0,1,1)),
#'       A=c(1,1,1,0,0,0,1,1,1,0),
#'       B=c(0,1,1,0,1,1,0,1,1,0),
#'       C=c(0,0,1,0,0,1,0,0,1,0))
#' data_test = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,1,0)),
#'       A=c(0,0,0,1,0,0,0,1),
#'       B=c(1,1,1,0,0,1,1,1),
#'       C=c(0,0,1,1,0,0,1,1))
#' svm(data_train,data_test,kernel='radial',degree=3)
#' svm(data_train,data_test,kernel='sigmoid')
#' svm(data_train,data_test,kernel='poly',degree=4,poly=1)

svm<-function(data_train,data_test,kernel='linear',degree=3,poly=0,includeplot=FALSE){
  if (poly==1){
    classifier = e1071::svm(formula = classification ~ .,
                            data = data_train,
                            type = 'C-classification',
                            kernel = kernel,degree=degree)
  }
  else if (kernel=='radial'){
    classifier = e1071::svm(formula = classification ~ .,
                            data = data_train,
                            type = 'C-classification',
                            kernel = 'radial',coef0=0, degree = 3)
  }
  else {
    classifier = e1071::svm(formula = classification ~ .,
                            data = data_train,
                            type = 'C-classification',
                            kernel = kernel)
  }
  training = suppressWarnings(classifier.performance(classifier,data_train,label = 'SVM - training performance',includeplot = includeplot))
  test = suppressWarnings(classifier.performance(classifier,data_test,label = 'SVM - test performance',includeplot = includeplot))
  
  return_list <- list("test" = test$accuracy,"testsensitivity"=test$sensitivity,"testspecificity"=test$specificity,"training" = training$accuracy, "trainsensitivity"=training$sensitivity,"trainspecificity"=training$specificity)
  return(return_list)
}

#' Linear SVM
#' Implements a linear SVM using the general svm function (for ease of use in feature selection)
#' @param data_train Training set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model is trained.
#' @param data_test Test set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model in tested.
#' @param includeplot Show performance scatter plot (default:FALSE)
#' @return List containing performance percentages, accessed using training (training accuracy), test (test accuracy), trainsensitivity, testsensitivity, trainspecificity, testspecificity.
#' @keywords svm
#' @keywords classification
#' @export
#' @examples
#' data_train = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,0,0,1,1)),
#'       A=c(1,1,1,0,0,0,1,1,1,0),
#'       B=c(0,1,1,0,1,1,0,1,1,0),
#'       C=c(0,0,1,0,0,1,0,0,1,0))
#' data_test = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,1,0)),
#'       A=c(0,0,0,1,0,0,0,1),
#'       B=c(1,1,1,0,0,1,1,1),
#'       C=c(0,0,1,1,0,0,1,1))
#' svmlinear(data_train,data_test)
svmlinear<-function(data_train,data_test,includeplot=FALSE){
  svm(data_train,data_test,kernel='linear',poly=0,includeplot=includeplot)
}

#' Polynomial degree 2 SVM
#' Implements a polynomial degree 2 SVM using the general svm function (for ease of use in feature selection)
#' @param data_train Training set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model is trained.
#' @param data_test Test set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model in tested.
#' @param includeplot Show performance scatter plot (default:FALSE)
#' @return List containing performance percentages, accessed using training (training accuracy), test (test accuracy), trainsensitivity, testsensitivity, trainspecificity, testspecificity.
#' @keywords svm
#' @keywords classification
#' @export
#' @examples
#' data_train = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,0,0,1,1)),
#'       A=c(1,1,1,0,0,0,1,1,1,0),
#'       B=c(0,1,1,0,1,1,0,1,1,0),
#'       C=c(0,0,1,0,0,1,0,0,1,0))
#' data_test = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,1,0)),
#'       A=c(0,0,0,1,0,0,0,1),
#'       B=c(1,1,1,0,0,1,1,1),
#'       C=c(0,0,1,1,0,0,1,1))
#' svmpolynomial2(data_train,data_test)
svmpolynomial2<-function(data_train,data_test,includeplot=FALSE){
  svm(data_train,data_test,kernel='polynomial',degree=2,poly=1,includeplot=includeplot)
}

#' Polynomial degree 3 SVM
#' Implements a polynomial degree 3 SVM using the general svm function (for ease of use in feature selection)
#' @param data_train Training set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model is trained.
#' @param data_test Test set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model in tested.
#' @param includeplot Show performance scatter plot (default:FALSE)
#' @return List containing performance percentages, accessed using training (training accuracy), test (test accuracy), trainsensitivity, testsensitivity, trainspecificity, testspecificity.
#' @keywords svm
#' @keywords classification
#' @export
#' @examples
#' data_train = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,0,0,1,1)),
#'       A=c(1,1,1,0,0,0,1,1,1,0),
#'       B=c(0,1,1,0,1,1,0,1,1,0),
#'       C=c(0,0,1,0,0,1,0,0,1,0))
#' data_test = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,1,0)),
#'       A=c(0,0,0,1,0,0,0,1),
#'       B=c(1,1,1,0,0,1,1,1),
#'       C=c(0,0,1,1,0,0,1,1))
#' svmpolynomial3(data_train,data_test)
svmpolynomial3<-function(data_train,data_test,includeplot=FALSE){
  svm(data_train,data_test,kernel='polynomial',degree=3,poly=1,includeplot=includeplot)
}

#' Polynomial degree 4 SVM
#' Implements a polynomial degree 4 SVM using the general svm function (for ease of use in feature selection)
#' @param data_train Training set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model is trained.
#' @param data_test Test set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model in tested.
#' @param includeplot Show performance scatter plot (default:FALSE)
#' @return List containing performance percentages, accessed using training (training accuracy), test (test accuracy), trainsensitivity, testsensitivity, trainspecificity, testspecificity.
#' @keywords svm
#' @keywords classification
#' @export
#' @examples
#' data_train = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,0,0,1,1)),
#'       A=c(1,1,1,0,0,0,1,1,1,0),
#'       B=c(0,1,1,0,1,1,0,1,1,0),
#'       C=c(0,0,1,0,0,1,0,0,1,0))
#' data_test = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,1,0)),
#'       A=c(0,0,0,1,0,0,0,1),
#'       B=c(1,1,1,0,0,1,1,1),
#'       C=c(0,0,1,1,0,0,1,1))
#' svmpolynomial4(data_train,data_test)
svmpolynomial4<-function(data_train,data_test,includeplot=FALSE){
  svm(data_train,data_test,kernel='polynomial',degree=4,poly=1,includeplot=includeplot)
}

#' Radial SVM
#' Implements a radial SVM using the general svm function (for ease of use in feature selection)
#' @param data_train Training set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model is trained.
#' @param data_test Test set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model in tested.
#' @param includeplot Show performance scatter plot (default:FALSE)
#' @return List containing performance percentages, accessed using training (training accuracy), test (test accuracy), trainsensitivity, testsensitivity, trainspecificity, testspecificity.
#' @keywords svm
#' @keywords classification
#' @export
#' @examples
#' data_train = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,0,0,1,1)),
#'       A=c(1,1,1,0,0,0,1,1,1,0),
#'       B=c(0,1,1,0,1,1,0,1,1,0),
#'       C=c(0,0,1,0,0,1,0,0,1,0))
#' data_test = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,1,0)),
#'       A=c(0,0,0,1,0,0,0,1),
#'       B=c(1,1,1,0,0,1,1,1),
#'       C=c(0,0,1,1,0,0,1,1))
#' svmradial(data_train,data_test)
svmradial<-function(data_train,data_test,includeplot=FALSE){
  svm(data_train,data_test,kernel='radial',degree=3,poly=0,includeplot=includeplot)
}

#' Sigmoid SVM
#' Implements a sigmoid SVM using general svm function (for ease of use in feature selection)
#' @param data_train Training set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model is trained.
#' @param data_test Test set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model in tested.
#' @param includeplot Show performance scatter plot (default:FALSE)
#' @return List containing performance percentages, accessed using training (training accuracy), test (test accuracy), trainsensitivity, testsensitivity, trainspecificity, testspecificity.
#' @keywords svm
#' @keywords classification
#' @export
#' @examples
#' data_train = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,0,0,1,1)),
#'       A=c(1,1,1,0,0,0,1,1,1,0),
#'       B=c(0,1,1,0,1,1,0,1,1,0),
#'       C=c(0,0,1,0,0,1,0,0,1,0))
#' data_test = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,1,0)),
#'       A=c(0,0,0,1,0,0,0,1),
#'       B=c(1,1,1,0,0,1,1,1),
#'       C=c(0,0,1,1,0,0,1,1))
#' svmsigmoid(data_train,data_test)
svmsigmoid<-function(data_train,data_test,includeplot=FALSE){
  svm(data_train,data_test,kernel='sigmoid',poly=0,includeplot=includeplot)
}
