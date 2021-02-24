#' @importFrom magrittr %>%
NULL

specfunc<-function(results){
  temp = results %>% dplyr::filter(results$y_real == 0)
  spectib <- as.data.frame(temp %>% dplyr::count(temp$Correct))
  if (nrow(spectib)==1) {
    if (spectib[1,'temp$Correct']=='yes'){
      specificity <- 1
    }
    else {
      specificity <- 0
    }}
  if (nrow(spectib)==2) {
    specificity <- spectib[1,2]/(spectib[1,2]+spectib[2,2])
  }
  else{specificity=0}
  specificity
}
sensfunc<-function(results){
  temp = results %>% dplyr::filter(results$y_real == 1)
  senstib <- as.data.frame(temp %>% dplyr::count(temp$Correct))
  if (nrow(senstib)==1) {
    if (senstib[1,'temp$Correct']=='yes'){
      sensitivity <- 1
    }
    else {
      sensitivity <- 0
    }}
  if (nrow(senstib)==2) {
    sensitivity<-senstib[1,2]/(senstib[2,2]+senstib[1,2])
  }
  else {sensitivity=0}
  sensitivity
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

  testdtree_y_pred <- stats::predict(fit,data_test, type = 'class')
  testdtree <- table(data_test[, 1], testdtree_y_pred)
  testdtree_accuracy_Test <- sum(diag(testdtree)) / sum(testdtree)

  traindtree_y_pred = stats::predict(fit,data_train, type = 'class')
  traindtree <- table(data_train[, 1], traindtree_y_pred)
  print(traindtree)
  traindtree_accuracy_Test <- sum(diag(traindtree)) / sum(traindtree)

  if (nlevels(testdtree_y_pred %>% factor(levels=c(0,1)))==1){
    if (levels(testdtree_y_pred %>% factor(levels=c(0,1)))=='1'){
      testsensitivity<-nrow(data_test %>% dplyr::filter(data_test$classification==1))/nrow(data_test)
      testspecificity<-0
    }
    if (levels(testdtree_y_pred %>% factor(levels=c(0,1)))=='0'){
      testsensitivity<-0
      testspecificity<-nrow(data_test %>% dplyr::filter(data_test$classification==0))/nrow(data_test)
    }
  }
  if (nlevels(testdtree_y_pred %>% factor(levels=c(0,1)))!=1){
    y_real  = data_test$classification %>% factor(levels=c(0,1))
    y_pred  = testdtree_y_pred %>% factor(levels=c(0,1))
    results = tibble::tibble(y_real  = data_test$classification %>% factor(levels=c(0,1)),
                     y_pred  = testdtree_y_pred %>% factor(levels=c(0,1)),
                     Correct = ifelse(y_real == y_pred,"yes","no") %>% factor(levels=c('yes','no')))
    testspecificity<-specfunc(results)
    testsensitivity<-sensfunc(results)
  }
  if (nlevels(traindtree_y_pred %>% factor(levels=c(0,1)))==1){
    if (levels(traindtree_y_pred %>% factor(levels=c(0,1)))=='1'){
      trainsensitivity<-nrow(data_train %>% dplyr::filter(data_train$classification==1))/nrow(data_train)
      trainspecificity<-0
    }
    if (levels(traindtree_y_pred %>% factor(levels=c(0,1)))=='0'){
      trainsensitivity<-0
      trainspecificity<-nrow(data_train %>% dplyr::filter(data_train$classification==0))/nrow(data_train)
    }
  }
  if (nlevels(traindtree_y_pred %>% factor(levels=c(0,1)))!=1){
    y_real  = data_train$classification %>% factor(levels=c(0,1))
    y_pred  = traindtree_y_pred %>% factor(levels=c(0,1))
    results = tibble::tibble(y_real  = data_train$classification %>% factor(levels=c(0,1)),
                     y_pred  = traindtree_y_pred %>% factor(levels=c(0,1)),
                     Correct = ifelse(y_real == y_pred,"yes","no") %>% factor(levels=c('yes','no')))
    trainspecificity<-specfunc(results)
    trainsensitivity<-sensfunc(results)
    if (includeplot==TRUE){
      y_real  = data_test$classification %>% factor(levels=c(0,1))
      y_pred  = testdtree_y_pred %>% factor(levels=c(0,1))
      results = tibble::tibble(y_real  = data_test$classification %>% factor(levels=c(0,1)),
                       y_pred  = testdtree_y_pred %>% factor(levels=c(0,1)),
                       Correct = ifelse(y_real == y_pred,"yes","no") %>% factor(levels=c('yes','no')))

      title = 'Performance on 10% unseen data - decision tree'
      xlab  = 'Measured (Real class)'
      ylab  = 'Predicted (Class assigned by decisiontree)'
      graphics::plot(ggplot2::ggplot(results,ggplot2::aes(x = testdtree_y_pred, y = data_test$classification, colour = results$Correct)) +
             ggplot2::geom_point() +
             ggplot2::ggtitle(label = title, subtitle = paste0("Accuracy = ", 100*round(testdtree_accuracy_Test,3),"%")) +
             xlab(xlab) +
             ylab(ylab) +
             ggplot2::scale_color_manual(labels = c('Yes', 'No'),
                                values = c('cornflowerblue','tomato')) +
             ggplot2::geom_jitter() +
             ggplot2::theme_bw())
    }
    return_list <- list("training"=traindtree_accuracy_Test,"test" = testdtree_accuracy_Test,"testsensitivity"=testsensitivity,"testspecificity"=testspecificity,"trainsensitivity"=trainsensitivity,"trainspecificity"=trainspecificity,"fit"=fit$frame)
    return(return_list)

  }}


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
  testrf_y_pred <- stats::predict(training.rf,newdata=data_test[-1])
  testrf <- table(data_test[, 1], testrf_y_pred)
  testrf_accuracy_Test <- sum(diag(testrf)) / sum(testrf)

  trainrf_y_pred = stats::predict(training.rf,newdata=data_train[-1])
  trainrf <- table(data_train[, 1], trainrf_y_pred)
  trainrf_accuracy_Test <- sum(diag(trainrf)) / sum(trainrf)

  if (nlevels(testrf_y_pred %>% factor(levels=c(0,1)))==1){
    if (levels(testrf_y_pred %>% factor(levels=c(0,1)))=='1'){
      testsensitivity<-nrow(data_test %>% dplyr::filter(data_test$classification==1))/nrow(data_test)
      testspecificity<-0
    }
    if (levels(testrf_y_pred %>% factor(levels=c(0,1)))=='0'){
      testsensitivity<-0
      testspecificity<-nrow(data_test %>% dplyr::filter(data_test$classification==0))/nrow(data_test)
    }
  }
  if (nlevels(testrf_y_pred %>% factor(levels=c(0,1)))!=1){
    y_real  = data_test$classification %>% factor(levels=c(0,1))
    y_pred  = testrf_y_pred %>% factor(levels=c(0,1))
    results = tibble::tibble(y_real  = data_test$classification %>% factor(levels=c(0,1)),
                     y_pred  = testrf_y_pred %>% factor(levels=c(0,1)),
                     Correct = ifelse(y_real == y_pred,"yes","no") %>% factor(levels=c('yes','no')))
    testspecificity<-specfunc(results)
    testsensitivity<-sensfunc(results)
  }
  if (nlevels(trainrf_y_pred %>% factor(levels=c(0,1)))==1){
    if (levels(trainrf_y_pred %>% factor(levels=c(0,1)))=='1'){
      trainsensitivity<-nrow(data_train %>% dplyr::filter(data_train$classification==1))/nrow(data_train)
      trainspecificity<-0
    }
    if (levels(trainrf_y_pred %>% factor(levels=c(0,1)))=='0'){
      trainsensitivity<-0
      trainspecificity<-nrow(data_train %>% dplyr::filter(data_train$classification==0))/nrow(data_train)
    }
  }
  if (nlevels(trainrf_y_pred %>% factor(levels=c(0,1)))!=1){
    y_real  = data_train$classification %>% factor(levels=c(0,1))
    y_pred  = trainrf_y_pred %>% factor(levels=c(0,1))
    results = tibble::tibble(y_real  = data_train$classification %>% factor(levels=c(0,1)),
                     y_pred  = trainrf_y_pred %>% factor(levels=c(0,1)),
                     Correct = ifelse(y_real == y_pred,"yes","no") %>% factor(levels=c('yes','no')))
    trainspecificity<-specfunc(results)
    trainsensitivity<-sensfunc(results)
  }
  if (includeplot==TRUE){
    y_real  = data_test$classification %>% factor(levels=c(0,1))
    y_pred  = testrf_y_pred %>% factor(levels=c(0,1))
    results = tibble::tibble(y_real  = data_test$classification %>% factor(levels=c(0,1)),
                     y_pred  = testrf_y_pred %>% factor(levels=c(0,1)),
                     Correct = ifelse(y_real == y_pred,"yes","no") %>% factor(levels=c('yes','no')))

    title = 'Performance on 10% unseen data - random forest'
    xlab  = 'Measured (Real class)'
    ylab  = 'Predicted (Class assigned by random forest)'
    graphics::plot(ggplot2::ggplot(results,ggplot2::aes(x = testrf_y_pred, y = data_test$classification, colour = results$Correct)) +
           ggplot2::geom_point() +
           ggplot2::ggtitle(label = title, subtitle = paste0("Accuracy = ", 100*round(testrf_accuracy_Test,3),"%")) +
           xlab(xlab) +
           ylab(ylab) +
           ggplot2::scale_color_manual(labels = c('Yes','No'),
                              values = c('cornflowerblue','tomato')) +
           ggplot2::geom_jitter() +
           ggplot2::theme_bw())
  }
  return_list <- list("training"=trainrf_accuracy_Test,"test" = testrf_accuracy_Test,"testsensitivity"=testsensitivity,"testspecificity"=testspecificity,"trainsensitivity"=trainsensitivity,"trainspecificity"=trainspecificity,"importance"=randomForest::importance(training.rf))
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
    classifier = e1071::svm(formula = stats::as.formula("classification~."),
                            data = data_train,
                            type = 'C-classification',
                            kernel = kernel,degree=degree,scale=FALSE)
  }
  else if (kernel=='radial'){
    classifier = e1071::svm(formula = stats::as.formula("classification~."),
                            data = data_train,
                            type = 'C-classification',
                            kernel = 'radial',coef0=0, degree = degree,scale=FALSE)
  }
  else {
    classifier = e1071::svm(formula = stats::as.formula("classification~."),
                            data = data_train,
                            type = 'C-classification',
                            kernel = kernel,scale=FALSE)
  }
  svmy_pred_train = stats::predict(classifier, newdata = subset(data_train,select=(colnames(data_train)[colnames(data_train)!='classification'])))
  svmcm_train = table(data_train[, 1], svmy_pred_train)
  svmaccuracy_Train <- sum(diag(svmcm_train)) / sum(svmcm_train)
  if (nlevels(svmy_pred_train %>% factor(levels=c(0,1)))==1){
    if (levels(svmy_pred_train %>% factor(levels=c(0,1)))=='1'){
      sensitivity_train<-nrow(data_train %>% dplyr::filter(data_train$classification==1))/nrow(data_train)
      specificity_train<-0
    }
    if (levels(svmy_pred_train %>% factor(levels=c(0,1)))=='0'){
      sensitivity_train<-0
      specificity_train<-nrow(data_train %>% dplyr::filter(data_train$classification==0))/nrow(data_train)
    }
  }
  if (nlevels(svmy_pred_train %>% factor(levels=c(0,1)))!=1){
    y_real  = data_train$classification %>% factor(levels=c(0,1))
    y_pred  = svmy_pred_train %>% factor(levels=c(0,1))
    results_train = tibble::tibble(y_real  = data_train$classification %>% factor(levels=c(0,1)),
                           y_pred  = svmy_pred_train %>% factor(levels=c(0,1)),
                           Correct = ifelse(y_real == y_pred,"yes","no") %>% factor(levels=c('yes','no')))
    sensitivity_train<-sensfunc(results_train)
    specificity_train<-specfunc(results_train)
  }


  svmy_pred = stats::predict(classifier, newdata = subset(data_test,select=(colnames(data_test)[colnames(data_test)!='classification'])))
  svmcm = table(data_test[, 1], svmy_pred)
  svmaccuracy_Test <- sum(diag(svmcm)) / sum(svmcm)
  if (nlevels(svmy_pred %>% factor(levels=c(0,1)))==1){
    if (levels(svmy_pred %>% factor(levels=c(0,1)))=='1'){
      sensitivity<-nrow(data_test %>% dplyr::filter(data_test$classification==1))/nrow(data_test)
      specificity<-0
    }
    if (levels(svmy_pred %>% factor(levels=c(0,1)))=='0'){
      sensitivity<-0
      specificity<-nrow(data_test %>% dplyr::filter(data_test$classification==0))/nrow(data_test)
    }
  }
  if (nlevels(svmy_pred %>% factor(levels=c(0,1)))!=1){
    y_real  = data_test$classification %>% factor(levels=c(0,1))
    y_pred  = svmy_pred %>% factor(levels=c(0,1))
    results = tibble::tibble(y_real  = data_test$classification %>% factor(levels=c(0,1)),
                     y_pred  = svmy_pred %>% factor(levels=c(0,1)),
                     Correct = ifelse(y_real == y_pred,"yes","no") %>% factor(levels=c('yes','no')))
    sensitivity<-sensfunc(results)
    specificity<-specfunc(results)
  }

  if (includeplot==TRUE){
    y_real  = data_test$classification %>% factor(levels=c(0,1))
    y_pred  = svmy_pred %>% factor(levels=c(0,1))
    results = tibble::tibble(y_real  = data_test$classification %>% factor(levels=c(0,1)),
                     y_pred  = svmy_pred %>% factor(levels=c(0,1)),
                     Correct = ifelse(y_real == y_pred,"yes","no") %>% factor(levels=c('yes','no')))

    title = 'Performance on 10% unseen data - SVM'
    xlab  = 'Measured (Real class)'
    ylab  = 'Predicted (Class assigned by SVM)'
    graphics::plot(ggplot2::ggplot(results,ggplot2::aes(x = svmy_pred, y = data_test$classification, colour = results$Correct)) +
           ggplot2::geom_point() +
           ggplot2::ggtitle(label = title, subtitle = paste0("Accuracy = ", 100*round(svmaccuracy_Test,3),"%")) +
           xlab(xlab) +
           ylab(ylab) +
           ggplot2::scale_color_manual(labels = c('Yes', 'No'),
                              values = c('cornflowerblue','tomato')) +
           ggplot2::geom_jitter() +
           ggplot2::theme_bw())
  }


  return_list <- list("test" = svmaccuracy_Test,"testsensitivity"=sensitivity,"testspecificity"=specificity,"training" = svmaccuracy_Train, "trainsensitivity"=sensitivity_train,"trainspecificity"=specificity_train)
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
