hybridfeatureselection <- function(k,training,test,features,ongoing=c(),lengthofongoing=list(),goodnessoffeature,runssofar=0,model=feamiR::svm_linear,mutprob=0.05,includePlot=FALSE,maxnumruns=50){

  #pick new features to make up k in total
  kfeatures <-sample(features,k,replace=TRUE)
  allfeatures<-union(features,ongoing)
  listoflengths<-lengthofongoing
  k_<-k+length(ongoing)
  combinedfeatures<-union(ongoing,kfeatures)

  #take the columns of the training and test sets which correspond to the chosen features
  dfonlyfeatures_training<-takefeaturecolumns(training,combinedfeatures)
  dfonlyfeatures_test<-takefeaturecolumns(test,combinedfeatures)

  #run forward feature selection to order the features and take those before the plateau. If includePlot is set to TRUE then a plot of accuracy for each feature up to the plateau is shown.
  ffs<-forwardfeatureselection(model,dfonlyfeatures_training,dfonlyfeatures_test,combinedfeatures)
  ongoingfeatures<-ffs$feature_list
  runssofar <- runssofar + 1

  #update goodness of feature measure
  for (feat in allfeatures){
    if (is.element(feat,ongoingfeatures)){
      goodnessoffeature[1,feat] <- goodnessoffeature[1,feat] + 1
    }
    else{
      goodnessoffeature[1,feat]<-max(goodnessoffeature[1,feat]-1,0)
    }
  }

  #add in features which were not chosen by forward feature selection but which have appeared in more than half of previous runs
  for (feat in combinedfeatures){
    if (!(is.element(feat,ongoingfeatures))&(goodnessoffeature[1,feat]>=runssofar/2)){
      ongoingfeatures<-append(ongoingfeatures,feat)
    }
  }

  #update accuracy, sensitivity and specificity values
  trainaccuracy <- ffs$accuracy
  testaccuracy<-ffs$testacc
  trainsens<-ffs$trainsens
  testsens<-ffs$testsens
  trainspec<-ffs$trainspec
  testspec<-ffs$testspec
  l<-length(ongoingfeatures)
  featureslessongoing<-allfeatures[!allfeatures %in% ongoingfeatures]

  #crossover
  pos<-sample(1:l,1,replace=TRUE)
  alt_set<-sample(featureslessongoing,l,replace=TRUE)
  ab<-union(utils::head(ongoingfeatures,pos),utils::tail(alt_set,l-pos))
  ba<-union(utils::head(alt_set,pos),utils::tail(ongoingfeatures,l-pos))
  training_bb<-takefeaturecolumns(training,alt_set)
  test_bb<-takefeaturecolumns(test,alt_set)
  training_ab<-takefeaturecolumns(training,ab)
  test_ab<-takefeaturecolumns(test,ab)
  training_ba<-takefeaturecolumns(test,ba)
  test_ba<-takefeaturecolumns(test,ba)
  model_bb<-model(training_bb,test_bb)
  model_ab<-model(training_ab,test_ab)
  model_ba<-model(training_ba,test_ba)
  acc_bb<-model_bb$training
  acc_ab<-model_ab$training
  acc_ba<-model_ba$training

  comparelist<-c(acc_bb,acc_ba,acc_ab,trainaccuracy)
  if (trainaccuracy!=max(comparelist)){
    if (max(comparelist)==acc_ba){
      trainaccuracy<-acc_ba
      ongoingfeatures<-ba
      lengthofongoing<-length(ba)
      testaccuracy<-model_ba$test
      trainspec<-model_ba$trainspecificity
      testspec<-model_ba$testspecificity
      trainsens<-model_ba$trainsensitivity
      testsens<-model_ba$testsensitivity
    }
    else if (max(comparelist)==acc_ab){
      trainaccuracy<-acc_ab
      ongoingfeatures<-ab
      lengthofongoing<-length(ab)
      testaccuracy<-model_ab$test
      trainspec<-model_ab$trainspecificity
      testspec<-model_ab$testspecificity
      trainsens<-model_ab$trainsensitivity
      testsens<-model_ab$testsensitivity
    }
    else {
      trainaccuracy<-acc_bb
      ongoingfeatures<-alt_set
      lengthofongoing<-length(alt_set)
      testaccuracy<-model_bb$test
      trainspec<-model_bb$trainspecificity
      testspec<-model_bb$testspecificity
      trainsens<-model_bb$trainsensitivity
      testsens<-model_bb$testsensitivity
    }
  }

  #mutation
  bernoulli<-sample(0:1,1,replace=TRUE,prob=c(1-mutprob,mutprob))
  if (bernoulli==1){
    #1. Pick a position j from 1 to l
    pos<-sample(1:l,1,replace=TRUE)
    #2. Pick a feature not in the current set and swap it for the other one
    newfeat<-sample(featureslessongoing,1,replace=TRUE)
    #3. Test the two accuracies and pick the better one
    altongoing<-ongoingfeatures
    altongoing[[pos]]<-newfeat
    alt_training<-takefeaturecolumns(training,altongoing)
    alt_test<-takefeaturecolumns(test,altongoing)
    alt_result<-model(alt_training,alt_test)
    alt_trainaccuracy<-alt_result$training
    if (alt_trainaccuracy>trainaccuracy){
      ongoingfeatures<-altongoing
      trainaccuracy<-alt_trainaccuracy
      testaccuracy<-alt_result$test
      trainsens<-alt_result$trainsensitivity
      testsens<-alt_result$testsensitivity
      trainspec<-alt_result$trainspecificity
      testspec<-alt_result$testspecificity
    }

  }
  featureslessongoing<-allfeatures[!allfeatures %in% ongoingfeatures]
  listoflengths<-append(listoflengths,l)

  #Termination conditions
  #1. If the feature set hasn't changed from the last run or we have reached the maximum number of runs then terminate
  if (identical(ongoing,ongoingfeatures)||runssofar>=maxnumruns){
    return_list <- list("feature_list" = ongoingfeatures, "training" = trainaccuracy,"test"=testaccuracy,"trainspecificity"=trainspec,"testspecificity"=testspec,"trainsensitivity"=trainsens,"testsensitivity"=testsens,"listofongoing"=listoflengths)
    return(return_list)
  }
  #2. If k good features have been found and 10 runs have passed then terminate
  if (l>=k_ & runssofar>=10){
    return_list <- list("feature_list" = ongoingfeatures, "training" = trainaccuracy,"test"=testaccuracy,"trainspecificity"=trainspec,"testspecificity"=testspec,"trainsensitivity"=trainsens,"testsensitivity"=testsens,"listofongoing"=listoflengths)
    return(return_list)
  }
  else{
    #3. Otherwise keep going
    hybridfeatureselection(k_-l,training,test,featureslessongoing,ongoingfeatures,listoflengths,goodnessoffeature,runssofar,model=model,mutprob=mutprob,includePlot=includePlot,maxnumruns=maxnumruns)}
}

#' Embryonic Genetic Algorithm.
#' Feature selection based on Embryonic Genetic Algorithms. It performs feature selection by maintaining an ongoing set of 'good' set of features which are improved run by run. It outputs training and test accuracy, sensitivity and specificity and a list of <=k features.
#' @param k Maximum number of features in the output feature set (default:30)
#' @param data_train Training set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model is trained.
#' @param data_test Test set: dataframe containing classification column and all other columns features. This is the dataset on which the decision tree model in tested.
#' @param mutprob Probability that mutation will be performed for each produced feature set from forward feature selection (default:0.05)
#' @param maxnumruns Maximum number of iterations after which the feature set will be output, if no other termination conditions have been met (default:50)
#' @param includePlot Show performance scatter plot (default:FALSE)
#' @return List containing (ordered list of) selected features, performance percentages, accessed using training (training accuracy), test (test accuracy), trainsensitivity, testsensitivity, trainspecificity, testspecificity. Also accessed using listofongoing is a list containing the length of the ongoing set at each stage.
#' @keywords eGA
#' @keywords embryonic
#' @keywords feature selection
#' @keywords genetic
#' @export
#' @examples
#' data_train = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,0,0,1,1)),
#'       A=c(1,1,1,0,0,0,1,1,1,0),
#'       B=c(0,1,1,0,1,1,0,1,1,0),
#'       C=c(0,0,1,0,0,1,0,0,1,0),
#'       D=c(0,1,1,0,0,0,1,0,0,0),
#'       E=c(1,0,1,0,0,1,0,1,1,0))
#' data_test = data.frame(
#'       classification=as.factor(c(1,1,0,0,1,1,1,0)),
#'       A=c(0,0,0,1,0,0,0,1),
#'       B=c(1,1,1,0,0,1,1,1),
#'       C=c(0,0,1,1,0,0,1,1),
#'       D=c(0,0,1,1,0,1,0,1),
#'       E=c(0,0,1,0,1,0,1,1))
#' data = read.csv(paste(system.file('samples/subsamples', package = "feamiR"),'/sample0.csv',sep=''))
#' data = rbind(head(data,50),tail(data,50))
#' data$classification = as.factor(data$classification)
#' ind <- sample(2,nrow(data),replace=TRUE,prob=c(0.8,0.2))
#' data_train <- data[ind==1,]
#' data_test <- data[ind==2,]
#' eGA(k=7,data_train,data_test,maxnumruns=3)
eGA <- function(k=30,data_train,data_test,mutprob=0.05,includePlot=FALSE,maxnumruns=50){
  listoffeatures <- colnames(data_train)[colnames(data_train)!='classification']
  goodnessoffeature<-data.frame(matrix(nrow = 1, ncol = length(listoffeatures)))
  for (i in (1:length(listoffeatures))){
    goodnessoffeature[1,i]<-0
  }
  colnames(goodnessoffeature)<-listoffeatures
  return(hybridfeatureselection(k=k,training=data_train,test=data_test,listoffeatures,c(),list(),goodnessoffeature,0,model=feamiR::svmlinear,mutprob=mutprob,includePlot=includePlot,maxnumruns=maxnumruns))
}
