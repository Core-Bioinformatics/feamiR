takefeaturecolumns<-function(sam,feat){
  featuretrain = subset(sam, select=c('classification'))
  feat<-feat[feat!='classification']
  for (i in 1:(length(feat))){
    row=feat[i]
    featuretrain[row]<-sam[row]
  }
  featuretrain
}

#Selects feature which maximises accuracy combined with ongoing features
nextbestfeature<-function(model,trainingsam,testsam,featuresleft,ongoingfeatures){
  training<-takefeaturecolumns(trainingsam,union(featuresleft,ongoingfeatures))
  test<-takefeaturecolumns(testsam,union(featuresleft,ongoingfeatures))
  maxaccuracy<-0
  maxfeat<-utils::head(featuresleft,1)
  for (feat in featuresleft){
    feat_trainingsample<-takefeaturecolumns(trainingsam,append(ongoingfeatures,feat))
    feat_testsample<-takefeaturecolumns(testsam,append(ongoingfeatures,feat))
    featresults<-model(feat_trainingsample,feat_testsample)
    feataccuracy<-featresults$training
    if (maxaccuracy<=feataccuracy){
      maxaccuracy<-feataccuracy
      max_train<-featresults$training
      max_test<- featresults$test
      max_testsens<-featresults$testsensitivity
      max_testspec<-featresults$testspecificity
      max_trainsens<-featresults$trainsensitivity
      max_trainspec<-featresults$trainspecificity
      maxfeat<-feat
    }
  }
  return_list <- list("feature" = maxfeat, "training_accuracy" = max_train,"test_accuracy"= max_test,"trainsens"=max_trainsens,"testsens"=max_testsens,"trainspec"=max_trainspec,"testspec"=max_testspec)
  return(return_list)
}

#Tests whether a position position is at a plateau
beforeandafter<-function(vect,centre,places,threshold){
  flag<-FALSE
  for (i in 1:places){
    if (centre+i<length(vect) & (abs(vect[centre]-vect[centre-i])>threshold || abs(vect[centre]-vect[centre+i])>threshold)){
      flag<-TRUE
    }
    else if (centre+i>=length(vect) & abs(vect[centre]-vect[centre-i])>threshold){
      flag<-TRUE
    }
  }
  flag
}


#' Forward Feature Selection.
#' Performs forward feature selection on the given list of features, placing them in order of discriminative power using a given model on the given dataset up to the accuracy plateau.
#' @param model The ML models used to classify the data, typically SVM with a given kernel
#' @param training Training dataset as a data.frame with classification column and column for each feature.
#' @param test Test dataset with matching columns to training.
#' @param featurelist List of features to order
#' @param includePlot Show number of features vs accuracy line plot (default:FALSE)
#' @return Ordered list of most discriminative features when classifying the dataset along with training and test accuracy, sensitivity and specificity
#' @keywords forward
#' @keywords feature selection
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
#' listoffeatures = colnames(data_train)[colnames(data_train)!='classification']
#' forwardfeatureselection(feamiR::svmlinear,data_train,data_test,listoffeatures)
forwardfeatureselection <-function(model=feamiR::svmlinear,training,test,featurelist,includePlot=FALSE){
  ongoing<-c()
  featureaccuracy<-data.frame(matrix(nrow=length(featurelist), ncol = 0))
  i<-1
  while (i<=length(featurelist)){
    b<-nextbestfeature(model,training,test,featurelist[!featurelist %in% ongoing],ongoing)
    ongoing<-c(ongoing,b$feature)
    featureaccuracy[i,'training_accuracy']<-b$training_accuracy
    featureaccuracy[i,'test_accuracy']<-b$test_accuracy
    featureaccuracy[i,'accuracy']<-b$training_accuracy
    featureaccuracy[i,'number of features']<-i
    featureaccuracy[i,'test_specificity']<-b$testspec
    featureaccuracy[i,'train_specificity']<-b$trainspec
    featureaccuracy[i,'test_sensitivity']<-b$testsens
    featureaccuracy[i,'train_sensitivity']<-b$trainsens
    i<-i+1
  }
  x<-featureaccuracy$'number of features'
  y<-featureaccuracy$accuracy
  if (length(x)<=5){
    loess30 = suppressWarnings(stats::loess(y ~ x, degree = 1, span = 0.5))
  }
  else{
  loess30<-suppressWarnings(stats::loess(y ~ x,degree=1,span=0.2))}
  y1 <- stats::predict(loess30,newdata=x,se=FALSE)
  flag<-TRUE
  i<-3
  numfeat<-length(x)
  while (i<=length(x) & flag==TRUE){
    if (length(x)<20 & i>3){
      flag<-beforeandafter(y1,i,3,0.01)
    }
    else {
      if (i>5){
        flag<-beforeandafter(y1,i,5,0.01)
      }
    }
    if (!flag){numfeat<-x[i]}
    i<-i+1
  }

  y2_<-featureaccuracy$test_accuracy
  loess30_2<-suppressWarnings(stats::loess(y2_ ~ x,degree=1,span=0.2))
  y2 <- stats::predict(loess30_2,newdata=x,se=FALSE)

  y4_<-featureaccuracy$training_accuracy
  loess30_4<-suppressWarnings(stats::loess(y4_ ~ x,degree=1,span=0.2))
  y4 <- stats::predict(loess30_4,newdata=x,se=FALSE)

  if (includePlot == TRUE){
    graphics::plot(utils::head(x,numfeat),utils::head(y,numfeat),col='red',main='Number of Features Against Accuracy up to Plateau',xlab='Number of features',ylab='Accuracy',ylim=c(0:1),type='b') + graphics::text(utils::head(x,numfeat), utils::head(y1,numfeat), utils::head(ongoing,numfeat), cex=0.6, pos=4,adj=0.5, col="red",srt=90)
    graphics::points(utils::head(x,numfeat),utils::head(y2_,numfeat),col='green')
    graphics::lines(utils::head(x,numfeat),utils::head(y1,numfeat),col='red')
    graphics::lines(utils::head(x,numfeat),utils::head(y2,numfeat),col='green')}
  return_list <- list("feature_list" = utils::head(ongoing,numfeat), "accuracy" = y1[numfeat],"testacc"=featureaccuracy$test_accuracy[numfeat],"trainacc"=featureaccuracy$training_accuracy[numfeat],"trainsens"=featureaccuracy$train_sensitivity[numfeat],"trainspec"=featureaccuracy$train_specificity[numfeat],"testsens"=featureaccuracy$test_sensitivity[numfeat],"testspec"=featureaccuracy$test_specificity[numfeat])
  return(return_list)
}
