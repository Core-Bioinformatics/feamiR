#' Tuning number of trees hyperparameter.
#' Trains random forests with a range of number of trees so the optimal number can be identified (using the resulting plot) with cross validation
#' @param data Dataset: dataframe containing classification column and all other column features. Both the training and test datasets will be taken from this dataset.
#' @param maxnum Maximum number of trees to be considered. All numbers between 1 and maxnum will be considered. Default: 100.
#' @param title Title to be used for the resulting boxplot
#' @param showplots TRUE if plots should be shown in standard output, FALSE is plots should be saved as jpg files. Default: TRUE.
#' @param output_prefix Prefix used for saving plots. If showplots==FALSE then plots are saved here. Otherwise, standard output.
#' @return Dataframe containing test and training accuracy, sensitivity and specificity
#' @keywords hyperparameter
#' @keywords random forest
#' @keywords tuning
#' @export
#' @examples
#' data = read.csv(paste(system.file('samples/subsamples', package = "feamiR"),'/sample0.csv',sep=''))
#' data = rbind(head(data,50),tail(data,50))
#' data$classification = as.factor(data$classification)
#' data = data[,2:ncol(data)]
#' selectrfnumtrees(data,5,'RF boxplots')
selectrfnumtrees<-function(data,maxnum=100,title='',showplots=TRUE, output_prefix=''){
  ker<-data.frame(matrix(nrow = 120, ncol = 5))
  colnames(ker)<-list('num_trees','accuracy','sensitivity','specificity','type')
  if (nrow(data)<=100){
    numiter=5
  }
  else {numiter=10}
  partition <- sample(numiter,nrow(data),replace=TRUE)
  lenofker = 0
  for (i in c(1:numiter)){
    data_train <- data[partition!=i,]
    data_test <- data[partition==i,]
    for (num_trees in c(1:maxnum)){
      results <- randomforest(data_train,data_test,numoftrees=num_trees)
      if (num_trees<10){
        ker[lenofker+1,1]<-paste('0',toString(num_trees),sep='')}
      else {ker[lenofker+1,1] = toString(num_trees)}
      ker[lenofker+1,2]<-results$test
      ker[lenofker+1,3]<-results$testsensitivity
      ker[lenofker+1,4]<-results$testspecificity
      ker[lenofker+1,5]<-'test'
      if (num_trees<10){
        ker[lenofker+2,1]<-paste('0',toString(num_trees),sep='')}
      else {ker[lenofker+2,1] = toString(num_trees)}
      ker[lenofker+2,2]<-results$training
      ker[lenofker+2,3]<-results$trainsensitivity
      ker[lenofker+2,4]<-results$trainspecificity
      ker[lenofker+2,5]<-'train'
      lenofker<-lenofker+2
    }
  }
  ker=stats::na.omit(ker)
  accuracy<-suppressWarnings(ggplot2::ggplot(ker, ggplot2::aes(x=ker$num_trees, y=ker$accuracy,colour=ker$type,fill=ker$type))+ggplot2::ylim(0.5,1)+ggplot2::scale_fill_grey(start=0.8, end=0.2)+
    ggplot2::geom_boxplot()+ggplot2::ggtitle(paste(title))+ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+ggplot2::labs(colour = "type",fill='type')+ggplot2::ylab('accuracy')+ggplot2::xlab('number of trees'))

  sensitivity<-suppressWarnings(ggplot2::ggplot(ker, ggplot2::aes(x=ker$num_trees, y=ker$sensitivity,colour=ker$type,fill=ker$type))+ggplot2::ylim(0.5,1)+ggplot2::scale_fill_grey(start=0.8, end=0.2)+
    ggplot2::geom_boxplot()+ggplot2::ggtitle(paste(title))+ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+ggplot2::labs(colour = "type",fill='type')+ggplot2::ylab('sensitivity')+ggplot2::xlab('number of trees'))

  specificity<-suppressWarnings(ggplot2::ggplot(ker, ggplot2::aes(x=ker$num_trees, y=ker$specificity,colour=ker$type,fill=ker$type))+ggplot2::ylim(0.5,1)+ggplot2::scale_fill_grey(start=0.8, end=0.2)+
    ggplot2::geom_boxplot()+ggplot2::ggtitle(paste(title))+ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+ggplot2::labs(colour = "type",fill='type')+ggplot2::ylab('specificity')+ggplot2::xlab('number of trees'))

  if (showplots){
    suppressWarnings(graphics::plot(accuracy))
    suppressWarnings(graphics::plot(sensitivity))
    suppressWarnings(graphics::plot(specificity))}
  else {
    grDevices::jpeg(paste(output_prefix,'select_rf_numtrees_accuracy.jpg',sep=''))
    suppressWarnings(graphics::plot(accuracy))
    grDevices::dev.off()
    grDevices::jpeg(paste(output_prefix,'select_rf_numtrees_sensitivity.jpg',sep=''))
    suppressWarnings(graphics::plot(sensitivity))
    grDevices::dev.off()
    grDevices::jpeg(paste(output_prefix,'select_rf_numtrees_specificity.jpg',sep=''))
    suppressWarnings(graphics::plot(specificity))
    grDevices::dev.off()}
  return(ker)
}

#' Tuning SVM kernel.
#' Trains SVMs with a range of kernels (linear, polynomial degree 2, 3 and 4, radial and sigmoid) using cross validation so the optimal kernel can be chosen (using the resulting plots). If specified (by showplots=FALSE) the plots are saved as jpegs.
#' @param data Dataset: dataframe containing classification column and all other column features. Both the training and test datasets will be taken from this dataset.
#' @param title Title to be used for the resulting boxplot
#' @param showplots TRUE if plots should be shown in standard output, FALSE is plots should be saved as jpg files.
#' @param output_prefix Prefix used for saving plots. If showplots==FALSE then plots are saved here. Otherwise, standard output.
#' @return Dataframe containing test and training accuracy, sensitivity and specificity
#' @keywords hyperparameter
#' @keywords SVM
#' @keywords kernel
#' @keywords tuning
#' @export
#' @examples
#' data = read.csv(paste(system.file('samples/subsamples', package = "feamiR"),'/sample0.csv',sep=''))
#' data = rbind(head(data,50),tail(data,50))
#' data$classification = as.factor(data$classification)
#' data = data[,2:ncol(data)]
#' selectsvmkernel(data,'SVM boxplots')
selectsvmkernel<-function(data,title='',showplots=TRUE,output_prefix=''){
  ker<-data.frame(matrix(nrow = 120, ncol = 5))
  colnames(ker)<-list('kernel','accuracy','sensitivity','specificity','type')
  if (nrow(data)<=100){
    numiter=5
  }
  else {numiter=10}
  partition <- sample(numiter,nrow(data),replace=TRUE)
  kernellist<-c('linear','radial','sigmoid')
  lenofker<-0
  for (i in c(1:numiter)){
    data_train <- data[partition!=i,]
    data_train$classification = as.factor(data_train$classification)
    data_test <- data[partition==i,]
    data_test$classification = as.factor(data_test$classification)
    for (kernel in kernellist){
      results<-svm(data_train,data_test,kernel=kernel,degree=0,poly=0)
      ker[lenofker+1,1]<-kernel
      ker[lenofker+1,2]<-results$test
      ker[lenofker+1,3]<-results$testsensitivity
      ker[lenofker+1,4]<-results$testspecificity
      ker[lenofker+1,5]<-'test'
      ker[lenofker+2,1]<-kernel
      ker[lenofker+2,2]<-results$training
      ker[lenofker+2,3]<-results$trainsensitivity
      ker[lenofker+2,4]<-results$trainspecificity
      ker[lenofker+2,5]<-'train'
      lenofker<-lenofker+2
    }
    for (degree in c(2:4)){
      kernel<-paste('poly',toString(degree))
      results<-svm(data_train,data_test,kernel='poly',degree=degree,poly=1)
      ker[lenofker+1,1]<-kernel
      ker[lenofker+1,2]<-results$test
      ker[lenofker+1,3]<-results$testsensitivity
      ker[lenofker+1,4]<-results$testspecificity
      ker[lenofker+1,5]<-'test'
      ker[lenofker+2,1]<-kernel
      ker[lenofker+2,2]<-results$training
      ker[lenofker+2,3]<-results$trainsensitivity
      ker[lenofker+2,4]<-results$trainspecificity
      ker[lenofker+2,5]<-'train'
      lenofker<-lenofker+2
    }
  }
  ker = stats::na.omit(ker)
  accuracy<-suppressWarnings(ggplot2::ggplot(ker, ggplot2::aes(x=ker$kernel, y=ker$accuracy,colour=ker$kernel,fill=ker$type))+ggplot2::ylim(0,1)+ggplot2::scale_fill_grey(start=0.8, end=0.2)+
    ggplot2::geom_boxplot()+ggplot2::ggtitle(paste(title))+ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+ggplot2::ylab('accuracy')+ggplot2::xlab('kernel')+ggplot2::labs(colour = "kernel",fill='type'))

  sensitivity<-suppressWarnings(ggplot2::ggplot(ker, ggplot2::aes(x=ker$kernel, y=ker$sensitivity,colour=ker$kernel,fill=ker$type))+ggplot2::ylim(0,1)+ggplot2::scale_fill_grey(start=0.8, end=0.2)+
    ggplot2::geom_boxplot()+ggplot2::ggtitle(paste(title))+ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+ggplot2::ylab('sensitivity')+ggplot2::xlab('kernel')+ggplot2::labs(colour = "kernel",fill='type'))

  specificity<-suppressWarnings(ggplot2::ggplot(ker, ggplot2::aes(x=ker$kernel, y=ker$specificity,colour=ker$kernel,fill=ker$type))+ggplot2::ylim(0,1)+ggplot2::scale_fill_grey(start=0.8, end=0.2)+
    ggplot2::geom_boxplot()+ggplot2::ggtitle(paste(title))+ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+ggplot2::ylab('specificity')+ggplot2::xlab('kernel')+ggplot2::labs(colour = "kernel",fill='type'))
  if (showplots){
  suppressWarnings(graphics::plot(accuracy))
  suppressWarnings(graphics::plot(sensitivity))
  suppressWarnings(graphics::plot(specificity))}
  else {
    grDevices::jpeg(paste(output_prefix,'select_svm_kernel_accuracy.jpg',sep=''))
    suppressWarnings(graphics::plot(accuracy))
    grDevices::dev.off()
    grDevices::jpeg(paste(output_prefix,'select_svm_kernel_sensitivity.jpg',sep=''))
    suppressWarnings(graphics::plot(sensitivity))
    grDevices::dev.off()
    grDevices::jpeg(paste(output_prefix,'select_svm_kernel_specificity.jpg',sep=''))
    suppressWarnings(graphics::plot(specificity))
    grDevices::dev.off()}
  return(ker)
}
