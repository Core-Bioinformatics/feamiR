
#' Tuning number of trees hyperparameter
#' This function trains random forests with a range of number of trees so the optimal number can be identified (using the resulting plot) with cross validation
#' @param data Dataset: dataframe containing classification column and all other column features. Both the training and test datasets will be taken from this dataset.
#' @param maxnum Maximum number of trees to be considered. All numbers between 1 and maxnum will be considered. Default: 100.
#' @param title Title to be used for the resulting boxplot
#' @return Dataframe containing test and training accuracy, sensitivity and specificity
#' @keywords hyperparameter
#' @keywords random forest
#' @keywords tuning
#' @export
#' @examples
#' selectrfnumtrees(data,100,'Boxplots showing distribution of training and test performance for random forest using number of trees 1 to 100')
selectrfnumtrees<-function(data,maxnum=100,title='',output_prefix=''){
  ker<-data.frame(matrix(nrow = 120, ncol = 5))
  colnames(ker)<-list('num_trees','accuracy','sensitivity','specificity','type')
  partition <- sample(10,nrow(data),replace=TRUE)
  lenofker = 0
  for (i in c(1:10)){
    data_train <- data[partition!=i,]
    data_test <- data[partition==i,]
    for (numtrees in c(1:maxnum)){
      results <- randomforest(data_train,data_test,numoftrees=num_trees)
      if (numtrees<10){
        ker[lenofker+1,1]<-paste('0',toString(numtrees),sep='')}
      else {ker[lenofker+1,1] = toString(numtrees)}
      ker[lenofker+1,2]<-results$test
      ker[lenofker+1,3]<-results$testsensitivity
      ker[lenofker+1,4]<-results$testspecificity
      ker[lenofker+1,5]<-'test'
      if (numtrees<10){
        ker[lenofker+2,1]<-paste('0',toString(numtrees),sep='')}
      else {ker[lenofker+2,1] = toString(numtrees)}
      ker[lenofker+2,2]<-results$training
      ker[lenofker+2,3]<-results$trainsensitivity
      ker[lenofker+2,4]<-results$trainspecificity
      ker[lenofker+2,5]<-'train'
      lenofker<-lenofker+2
    }
  }
  accuracy<-ggplot(na.omit(ker), aes(x=num_trees, y=accuracy,colour=type,fill=type))+ylim(0.5,1)+scale_fill_grey(start=0.8, end=0.2)+
    geom_boxplot()+ggtitle(paste(title))+theme(plot.title = element_text(hjust = 0.5))
  plot(accuracy)
  sensitivity<-ggplot(na.omit(ker), aes(x=num_trees, y=sensitivity,colour=type,fill=type))+ylim(0.5,1)+scale_fill_grey(start=0.8, end=0.2)+
    geom_boxplot()+ggtitle(paste(title))+theme(plot.title = element_text(hjust = 0.5))
  plot(sensitivity)
  specificity<-ggplot(na.omit(ker), aes(x=num_trees, y=specificity,colour=type,fill=type))+ylim(0.5,1)+scale_fill_grey(start=0.8, end=0.2)+
    geom_boxplot()+ggtitle(paste(title))+theme(plot.title = element_text(hjust = 0.5))
  plot(specificity)
  return(ker)
}

#' Tuning SVM kernel
#' This function trains SVMs with a range of kernels (linear, polynomial degree 2, 3 and 4, radial and sigmoid) using cross validation so the optimal kernel can be chosen (using the resulting plots). If specified (by showplots=F) the plots are saved as jpegs.
#' @param data Dataset: dataframe containing classification column and all other column features. Both the training and test datasets will be taken from this dataset.
#' @param title Title to be used for the resulting boxplot
#' @return Dataframe containing test and training accuracy, sensitivity and specificity
#' @keywords hyperparameter
#' @keywords SVM
#' @keywords kernel
#' @keywords tuning
#' @export
#' @examples
#' selectsvmkernel(data,'Boxplots showing distribution of training and test performance for random forest using number of trees 1 to 100')
selectsvmkernel<-function(data,title,showplots=T,output_prefix=''){
  ker<-data.frame(matrix(nrow = 120, ncol = 5))
  colnames(ker)<-list('kernel','accuracy','sensitivity','specificity','type')
  partition <- sample(10,nrow(data),replace=TRUE)
  kernellist<-c('linear','radial','sigmoid')
  lenofker<-0
  for (i in c(1:10)){
    data_train <- data[partition!=i,]
    data_test <- data[partition==i,]
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
  accuracy<-ggplot(na.omit(ker), aes(x=kernel, y=accuracy,colour=kernel,fill=type))+ylim(0.5,1)+scale_fill_grey(start=0.8, end=0.2)+
    geom_boxplot()+ggtitle(paste(title))+theme(plot.title = element_text(hjust = 0.5))

  sensitivity<-ggplot(na.omit(ker), aes(x=kernel, y=sensitivity,colour=kernel,fill=type))+ylim(0.5,1)+scale_fill_grey(start=0.8, end=0.2)+
    geom_boxplot()+ggtitle(paste(title))+theme(plot.title = element_text(hjust = 0.5))
  
  specificity<-ggplot(na.omit(ker), aes(x=kernel, y=specificity,colour=kernel,fill=type))+ylim(0.5,1)+scale_fill_grey(start=0.8, end=0.2)+
    geom_boxplot()+ggtitle(paste(title))+theme(plot.title = element_text(hjust = 0.5))
  if (showplots){
  plot(accuracy)
  plot(sensitivity)
  plot(specificity)}
  else {
    jpeg(paste(output_prefix,'select_svm_kernel_accuracy.jpg',sep=''))
    plot(accuracy)
    dev.off()
    jpeg(paste(output_prefix,'select_svm_kernel_sensitivity.jpg',sep=''))
    plot(sensitivity)
    dev.off()
    jpeg(paste(output_prefix,'select_svm_kernel_specificity.jpg',sep=''))
    plot(specificity)
    dev.off()}
  return(ker)
}