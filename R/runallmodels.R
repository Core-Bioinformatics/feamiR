
#' Run all models
#' This function trains and tests Decision Tree, Random Forest and SVM models on 100 subsamples and provides a summary of the results, to select the best model. The number of trees and kernel chosen by select_svm_kernel and select_rf_numtrees should be used for SVM and Random Forest respectively. We can use this function to inform feature selection, using a Decision Tree voting scheme and a Random Forest measure based on the Gini index.
#' @param num_trees Number of trees for random forest (selected using select_rf_numtrees)
#' @param kernel Kernel for SVM (select using select_svm_kernel)
#' @param degree Degree for SVM kernel (not necessary for linear or sigmoid functions)
#' @param poly 1 if polynomial kernel is used, 0 if linear, radial or sigmoid.
#' @param file_path Where the 100 subsample files are found (e.g. if sample 10 is at 'subsamples/sample10.csv' then file_path should be 'subsamples/sample')
#' @return The function will output a data.frame of the achieved test and training accuracy, sensitivity and specificity for each model on each subsample. Summary boxplots showing accuracy, sensitivity and specificity for each model will be produced. The function will also output dtreevote containing the features used in the decision trees for each subsample and the level of the tree at which they appear. Finally, the function outputs ongoingginis which contains the Gini index for each feature in the Random Forest for each subsample. The first column of dtreevote contains the number of runs for which each feature was used which can be used for feature selection. The first column of ongoingginis contains the cumulative Gini index for each feature across the 100 runs which can be used for feature selection. 
#' @keywords decision tree
#' @keywords dtree
#' @keywords classification
#' @export
#' @examples
#' "runallmodels(num_trees=50,kernel='radial',degree=3,poly=1,file_path='/subsamples/sample')"

runallmodels<-function(num_trees=20,kernel='linear',degree=3,poly=0,file_path=file_path){
  l<-data.frame(matrix(nrow = 100, ncol = 5))
  colnames(l)<-list('model','accuracy','sensitivity','specificity','type')
  rf<-data.frame(matrix(nrow = 100, ncol = 5))
  colnames(rf)<-list('model','accuracy','sensitivity','specificity','type')
  svm <- data.frame(matrix(nrow = 100, ncol = 5))
  colnames(svm)<-list('model','accuracy','sensitivity','specificity','type')
  currentrun = 1
  
  
  #Initialise dtreevote which will be voting scheme for dtree
  coltraining <- utils::read.csv(paste(file_path,toString(1),".csv", sep=""), header = TRUE)
  coltraining = subset(coltraining, select = -c(X) )
  y<-data.frame(matrix(nrow = ncol(coltraining),ncol=101))
  rownames(y) <- colnames(coltraining)
  for (names in rownames(y)){
    y[names,1]<-0
  }
  
  ongoingginis<-data.frame(matrix(nrow=ncol(coltraining),ncol=101))
  rownames(ongoingginis)<-colnames(coltraining)
  for (names in rownames(ongoingginis)){
    ongoingginis[names,1]<-0
  }
  
  
  
  for (i in (0:99)){
    training <- utils::read.csv(paste(file_path,toString(i),".csv", sep=""), header = TRUE,colClasses=c('classification'='factor'))
    training = subset(training,select = -c(X))
    listoffeatures <- colnames(subset(training, select = -c(classification)))
    shuffledtrain <- training[sample(nrow(training)),]
    shuffledtrain$classification <- as.factor(shuffledtrain$classification)
    ind <- sample(2,nrow(shuffledtrain),replace=TRUE,prob=c(0.8,0.2))
    data_train <- shuffledtrain[ind==1,]
    data_test <- shuffledtrain[ind==2,]
    
    
    #Decision trees
    dtree <- decisiontree(data_train,data_test)
    names <- data.frame(index=row.names(dtree$fit),var=dtree$fit['var'])
    names <- subset(names,strtoi(index)<(2^10))
    for (p in (10:1)){
      subsetnames <- subset(names,strtoi(index)<2^p)
      subsetnames<-subset(subsetnames,strtoi(index)>=2^(p-1))
      for (n in subsetnames$var){
        if (n != "<leaf>"){
          y[n,i+1]<-p
        }
      }
    }
    for (n in colnames(training)) {
      if (n %in% (names$var)){
        y[n,1]<-y[n,1]+1
      }
    }
    l[currentrun,5]='train'
    l[currentrun,4]<-dtree$trainspecificity
    l[currentrun,3]<-dtree$trainsensitivity
    l[currentrun,2]<-dtree$training
    l[currentrun,1]<-'dtree'
    
    l[currentrun+1,5]='test'
    l[currentrun+1,4]<-dtree$testspecificity
    l[currentrun+1,3]<-dtree$testsensitivity
    l[currentrun+1,2]<-dtree$test
    l[currentrun+1,1]<-'dtree'
    
    rforest <- randomforest(data_train,data_test,num_trees)
    ginis<-rforest$importance[,4]

    for (feat in listoffeatures){
      ongoingginis[feat,1]<-ongoingginis[feat,1]+ginis[feat]
      ongoingginis[feat,i+1]<-ginis[feat]
    }
    rf[currentrun,5]='train'
    rf[currentrun,4]<-rforest$trainspecificity
    rf[currentrun,3]<-rforest$trainsensitivity
    rf[currentrun,2]<-rforest$training
    rf[currentrun,1]<-'rf'
    rf[currentrun+1,5]='test'
    rf[currentrun+1,4]<-rforest$testspecificity
    rf[currentrun+1,3]<-rforest$testsensitivity
    rf[currentrun+1,2]<-rforest$test
    rf[currentrun+1,1]<-'rf'
    
    svm_res <- svm(data_train,data_test,kernel=kernel,degree=degree,poly=poly)
    
    svm[currentrun,5]='train'
    svm[currentrun,4]<-svm_res$trainspecificity
    svm[currentrun,3]<-svm_res$trainsensitivity
    svm[currentrun,2]<-svm_res$train
    svm[currentrun,1]<-'svm'
    svm[currentrun+1,5]='test'
    svm[currentrun+1,4]<-svm_res$testspecificity
    svm[currentrun+1,3]<-svm_res$testsensitivity
    svm[currentrun+1,2]<-svm_res$test
    svm[currentrun+1,1]<-'svm'
    currentrun = currentrun +2
      }
  
  
  accuracies0<-rbind(l,rf)
  accuracies1<-rbind(accuracies0,svm)
  

  p<-ggplot2::ggplot(stats::na.omit(accuracies1), aes(x=model, y=accuracy,colour=model,fill=stats::reorder(type,-accuracy)))+ylim(0.5,1)+ylab('accuracy')+scale_fill_grey(start=0.2, end=0.8)+
    geom_boxplot()+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  p<-ggplot2::ggplot(stats::na.omit(accuracies1), aes(x=model, y=sensitivity,colour=model,fill=stats::reorder(type,-accuracy)))+ylim(0.5,1)+scale_fill_grey(start=0.2, end=0.8)+
    geom_boxplot()+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  p<-ggplot2::ggplot(stats::na.omit(accuracies1), aes(x=model, y=specificity,colour=model,fill=stats::reorder(type,-accuracy)))+ylim(0.5,1)+scale_fill_grey(start=0.2, end=0.8)+
    geom_boxplot()+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  return_list<-list("resultdf"=accuracies1,"dtreevote"=y,"ongoingginis"=ongoingginis)
  return(return_list)
}

#' Decision tree voting scheme
#' This function implements a feature selection approach based on Decision Trees, using a voting scheme across the top levels on trees trained on multiple subsamples.
#' @param num_runs Number of subsamples to use for voting scheme (default: 100)
#' @param num_levels Number of levels in each tree to consider. Only the features which appear in the top num_levels levels of the trees (from the root) will be counted
#' @param file_path Where the num_runs subsample files are found (e.g. if sample 10 is at 'subsamples/sample10.csv' then file_path should be 'subsamples/sample'). There must be enough samples to fulfill num_runs runs.
#' @return Outputs a dataframe containing (first column) total number of appearances of each feature (each row is a feature). The rest of the columns represent 1 run each and contain the level at which the feature appears. 
#' @keywords decision tree
#' @keywords dtree
#' @keywords voting scheme
#' @export
#' @examples
#' "dtreevoting(num_runs=100,num_levels=10,file_path='/subsamples/sample')"

dtreevoting<-function(num_runs=100,num_levels=10,file_path=file_path){
  coltraining <- utils::read.csv(paste(file_path,toString(1),".csv", sep=""), header = TRUE)
  coltraining = subset(coltraining, select = -c(X) )
  y<-data.frame(matrix(nrow = ncol(coltraining),ncol=num_runs+1))
  rownames(y) <- colnames(coltraining)
  for (names in rownames(y)){
    y[names,1]<-0
  }
  
  for (i in (0:num_runs-1)){
    training <- utils::read.csv(paste(file_path,toString(i),".csv", sep=""), header = TRUE,colClasses=c('classification'='factor'))
    training = subset(training,select = -c(X))
    listoffeatures <- colnames(subset(training, select = -c(classification)))
    shuffledtrain <- training[sample(nrow(training)),]
    shuffledtrain$classification <- as.factor(shuffledtrain$classification)
    ind <- sample(2,nrow(shuffledtrain),replace=TRUE,prob=c(0.8,0.2))
    data_train <- shuffledtrain[ind==1,]
    data_test <- shuffledtrain[ind==2,]

    dtree <- decisiontree(data_train,data_test)
    names <- data.frame(index=row.names(dtree$fit),var=dtree$fit['var'])
    names <- subset(names,strtoi(index)<(2^num_levels))
    for (p in (num_levels:1)){
      subsetnames <- subset(names,strtoi(index)<2^p)
      subsetnames<-subset(subsetnames,strtoi(index)>=2^(p-1))
      for (n in subsetnames$var){
        if (n != "<leaf>"){
          y[n,i+1]<-p
        }
      }
    }
    for (n in colnames(training)) {
      if (n %in% (names$var)){
        y[n,1]<-y[n,1]+1
      }
    }
  }
  return(y)
}



#' Random Forest cumulative MeanDecreaseGini feature selection
#' This function implements a feature selection approach based on cumulative MeanDecreaseGini using Random Forests trained on multiple subsamples.
#' @param num_runs Number of subsamples to use for voting scheme (default: 100)
#' @param num_trees Number of trees for random forest (selected using select_rf_numtrees)
#' @param file_path Where the num_runs subsample files are found (e.g. if sample 10 is at 'subsamples/sample10.csv' then file_path should be 'subsamples/sample'). There must be enough samples to fulfill num_runs runs.
#' @return The function will output a data.frame with cumulative mean decrease in Gini for each feature in the first columns (each row is a feature) and the rest of the column containing mean decrease in Gini for each of the num_runs runs.
#' @keywords random forest
#' @keywords rf
#' @keywords feature selection
#' @keywords Gini
#' @export
#' @examples
#' "rfgini(num_runs=100,num_trees=30,file_path='/subsamples/sample')"

rfgini<-function(num_runs=100,num_trees=30,file_path=file_path){
  coltraining <- utils::read.csv(paste(file_path,toString(1),".csv", sep=""), header = TRUE)
  coltraining = subset(coltraining, select = -c(X) )
  ongoingginis<-data.frame(matrix(nrow=ncol(coltraining),ncol=101))
  rownames(ongoingginis)<-colnames(coltraining)
  for (names in rownames(ongoingginis)){
    ongoingginis[names,1]<-0
  }

  for (i in (0:num_runs-1)){
    training <- utils::read.csv(paste(file_path,toString(i),".csv", sep=""), header = TRUE,colClasses=c('classification'='factor'))
    training = subset(training,select = -c(X))
    listoffeatures <- colnames(subset(training, select = -c(classification)))
    shuffledtrain <- training[sample(nrow(training)),]
    shuffledtrain$classification <- as.factor(shuffledtrain$classification)
    ind <- sample(2,nrow(shuffledtrain),replace=TRUE,prob=c(0.8,0.2))
    data_train <- shuffledtrain[ind==1,]
    data_test <- shuffledtrain[ind==2,]
    
    
    rforest <- randomforest(data_train,data_test,num_trees)
    ginis<-rforest$importance[,4]
    
    for (feat in listoffeatures){
      ongoingginis[feat,1]<-ongoingginis[feat,1]+ginis[feat]
      ongoingginis[feat,i+1]<-ginis[feat]
    }
  }

  return(ongoingginis)
}
