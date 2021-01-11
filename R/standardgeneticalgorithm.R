
#' Standard Genetic Algorithm.
#' Implements a standard genetic algorithm using GA package (\link[GA]{ga}) with a fitness function specialised for feature selection.
#' @param model The ML models used to classify the data, typically SVM with a given kernel
#' @param k Maximum number of features to be output.
#' @param training Training dataset as a data.frame with classification column and column for each feature.
#' @param test Test dataset with matching columns to training.
#' @param parallel Specifies whether GA should be run sequentially or in parallel (default: T)
#' @param mutprob The probability that an individual undergoes mutation in a particular iteration (default: 0.1)
#' @param crossprob The probability of crossover between pairs of individuals (default: 0.8)
#' @param popsize The size of the solution population (default:20)
#' @param maxiter The maximum number of iterations to run before termination (default: 1000)
#' @param maxiter_withoutimprovement The maximum number of consecutive iterations without improvement to fitness before termination (default: 300)
#' @param numberpassedon The number of best fitness individuals to be passed on to the next generation in each iteration (default: 3)
#' @param plot Specifies whether GA plot should be shown (default: F)
#' @return Set (unordered) of <=k features and training and test accuracy, sensitivity and specificity using these features.
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
#' geneticalgorithm(
#'   feamiR::svmlinear,
#'   k=2,
#'   data_train,
#'   data_test,
#'   parallel=FALSE,
#'   maxiter=5,
#'   maxiter_withoutimprovement=5,
#'   popsize=10)

geneticalgorithm <- function(model=feamiR::svmlinear,k=30,training,test,parallel=T,mutprob=0.1,crossprob=0.8,popsize=20,maxiter=1000,maxiter_withoutimprovement=300,numberpassedon=3,plot=F){
  param_nBits=ncol(training)-1
  col_names=utils::tail(colnames(training),param_nBits)
  custom_fitness <- function(vars){
    currentvarnames=col_names[stats::setNames(vars,col_names)==1]
    if (length(currentvarnames)==0){
      return(0)}
    varfilt<-currentvarnames[!(currentvarnames%in%col_names)]
    if (length(varfilt)>0){
      return(0)
    }
    tr<-takefeaturecolumns(training,currentvarnames)
    te<-takefeaturecolumns(test,currentvarnames)

    if (length(currentvarnames)<=k) {
      acc<-model(tr,te)
      return(acc$training)
    }
    else{
      acc<-model(tr,te)
      return(acc$training/length(currentvarnames))
    }
  }
  ga_GA_1 = GA::ga(fitness = function(vars) custom_fitness(vars = vars),
               type = "binary", # optimization data type  # cross-over method
               elitism = numberpassedon, # best N indiv. to pass to next iteration
               pmutation = mutprob, # mutation rate prob
               pcrossover = crossprob, #crossover rate prob
               popSize = popsize, # the number of indivduals/solutions
               nBits = param_nBits, # total number of variables
               names=col_names, # variable name
               run=maxiter_withoutimprovement, # max iter without improvement (stopping criteria)
               maxiter = maxiter, # total runs or generations
               monitor=plot, # plot the result at each iteration
               keepBest = TRUE, # keep the best solution at the end
               parallel = parallel, # allow parallel procesing
               seed=84211 # for reproducibility purposes
  )
  sol=(ga_GA_1@solution[1,])
  acc=(ga_GA_1@fitnessValue)
  solfeatures=col_names[sol==1]
  data_train<-takefeaturecolumns(training,solfeatures)
  data_test<-takefeaturecolumns(test,solfeatures)
  run<-model(data_train,data_test)

  return_list <- list("feature_list" = solfeatures, "fitness" = acc,"testaccuracy"=run$test,"trainspec"=run$trainspecificity,"testspec"=run$testspecificity,"trainsens"=run$trainsensitivity,"testsens"=run$testsensitivity)
  return(return_list)
}
