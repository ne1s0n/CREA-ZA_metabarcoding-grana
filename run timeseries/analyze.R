library(e1071)
library(plyr)
library(madbito)
setwd('~/research/CREA-ZA_metabarcoding-grana/run timeseries/')
df = read.table('filtrato_finale_caseifici.csv', stringsAsFactors = FALSE, sep=';', header = TRUE)


# DATA SETUP --------------------------------------------------------------
#separating meta data from species abundances
meta = df[, c("Samples", "Dairy", "Province", "Month")]
abu = as.matrix(df[, setdiff(colnames(df), colnames(meta))])

#relative abundances
abu_rel = abu
for(i in 1:nrow(abu)){
  abu_rel[i,] = abu[i,] / sum(abu[i,]) 
}

#dominant species, i.e. all species that represent at least 1%
#of the abundances in a single sample
dom = colSums(abu_rel > 0.01) > 0
dom_names = colnames(abu)[dom]

# SUPPORT FUNCTIONS -------------------------------------------------------
KNN_classifier = function(X, Y, Ks, folds, reps){
  acc = NULL
  
  for (rep in 1:reps){
    folds_list = stratified_crossvalidation_folds(classes = Y, folds.num = folds)
    for (K in Ks){
      tmp = NULL
      for (f in 1:folds){
        #separating abundances
        X_train = X[folds_list != f,]
        X_test  = X[folds_list == f,]
        
        #separating labels
        Y_train = factor(Y[folds_list != f])
        Y_test  = factor(Y[folds_list == f])
        
        #training a classifier
        m = gknn(x = X_train, y = Y_train, k = K, type = 'class')
        Y_hat = predict(m, X_test)
        
        #accuracy
        tmp = c(tmp, sum(Y_hat == Y_test) / length(Y_test))
      }
      acc = rbind(acc, data.frame(
        K = K, 
        accuracy = mean(tmp)
      ))
    }}
  
  #averaging over repetitions
  acc = ddply(acc, .(K), function(x){
    return(data.frame(accuracy = mean$accuracy))
  })
  
  #and we are done
  return(acc)
}
SVM_classifier = function(X, Y, kernel, gammas, Cs, folds, reps){
  acc = NULL
  
  for (rep in 1:reps){
    folds_list = stratified_crossvalidation_folds(classes = Y, folds.num = folds)
    for (gamma in gammas){
    for (C in Cs){
      for (f in 1:folds){
        #separating abundances
        X_train = X[folds_list != f,]
        X_test  = X[folds_list == f,]
        
        #separating labels
        Y_train = factor(Y[folds_list != f])
        Y_test  = factor(Y[folds_list == f])
        
        #training a classifier
        if (kernel == 'linear'){
          m = svm(x = X_train, y = Y_train, kernel = kernel, cost = C, scale = FALSE)  
        }else{
          m = svm(x = X_train, y = Y_train, kernel = kernel, gamma = gamma, cost = C, scale = FALSE)
        }
        Y_hat = predict(m, X_test)
        
        #accuracy
        acc = rbind(acc, data.frame(
          kernel = kernel,
          gamma = gamma,
          cost = C,
          accuracy = sum(Y_hat == Y_test) / length(Y_test))
        )
      }
    }#looping over Cs
    }#looping over gammas
    }#looping over reps
  
  #averaging over repetitions and folds
  acc = ddply(acc, .(kernel, gamma, cost), function(x){
    return(data.frame(accuracy = mean(x$accuracy)))
  })
  
  #and we are done
  return(acc)
}

# SVM CLASSIFIER ----------------------------------------------------------
acc_prov = SVM_classifier(X=abu_rel, Y = meta$Province, kernel = 'linear', gammas = NA, Cs = c(0.01, 0.1, 1, 10), reps = 50, folds = 5)
write.csv(acc_prov, file = 'classifiers/SVM-lin_prov.csv', row.names = FALSE)
acc_prov = SVM_classifier(X=abu_rel, Y = meta$Province, kernel = 'radial', gammas = 10^(-3:+3), Cs = 10^(-3:+3), reps = 50, folds = 5)
write.csv(acc_prov, file = 'classifiers/SVM-rbf_prov.csv', row.names = FALSE)
stop('here')

acc_dairy = SVM_classifier(X=abu_rel, Y = meta$Dairy, kernel = 'linear', gammas = NA, Cs = c(0.01, 0.1, 1, 10), reps = 50, folds = 5)
write.csv(acc_dairy, file = 'classifiers/SVM-lin_dairy.csv', row.names = FALSE)
acc_dairy = SVM_classifier(X=abu_rel, Y = meta$Dairy, kernel = 'radial', gammas = 10^(-3:+3), Cs = 10^(-3:+3), reps = 50, folds = 5)
write.csv(acc_dairy, file = 'classifiers/SVM-rbf_dairy.csv', row.names = FALSE)
stop('here')

# KNN CLASSIFIER --------------------------------------------------------------
if (FALSE){
  acc_dairy = KNN_classifier(X=abu_rel, Y = meta$Dairy,    Ks = 1:10, reps = 100, folds = 5)
  write.csv(acc_dairy, file = 'classifiers/KNN_dairy.csv', row.names = FALSE)
  acc_prov  = KNN_classifier(X=abu_rel, Y = meta$Province, Ks = 1:10, reps = 100, folds = 5)
  write.csv(acc_dairy, file = 'classifiers/KNN_prov.csv', row.names = FALSE)
}




