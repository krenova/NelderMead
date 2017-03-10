source('lib/lib_NelderMead.R')
data(iris)

#Creating the Various data table types
library("Matrix")
dat = iris
dat$Sepal.LenWidRatio = dat[,"Sepal.Length"] / dat[,"Sepal.Width"]
dat$Petal.LenWidRatio = dat[,"Petal.Length"] / dat[,"Petal.Width"]
dat$SepalPetal.Length = dat[,"Sepal.Width"] + dat[,"Petal.Length"]
dat$SepalPetal.Width = dat[,"Petal.Length"]/dat[,"Petal.Width"]
dat$SepalPetal.LenWidRatio = dat[,"SepalPetal.Length"]/dat[,"SepalPetal.Width"]
dat$Len.SepPetRatio = dat[,"Sepal.Length"]/dat[,"Petal.Length"]
dat$Width.SepPetRatio = dat[,"Sepal.Width"]/dat[,"Petal.Width"]
dat$SepLenPetWidRatio = dat[,"Sepal.Length"]/dat[,"Petal.Width"]
dat$SepWidPetLenRatio = dat[,"Sepal.Width"]/dat[,"Petal.Length"]

dat_mod = model.matrix(Species~.-1, data = dat, sparse = TRUE) #
dat_sparse = sparse.model.matrix(Species~.-1, data = dat)
dat_y = as.numeric(dat[,"Species"])-1


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#----------------------Examplary Application of Nelder Mead (xgboost)--------------------------------------
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

xgb_wrap_obj_ex <- function(param_fixed, vtx, dat_x, dat_y, nfolds = 3, pred = FALSE, verb = FALSE ){ 
  return(   xgb.cv( param=param_fixed, nround = vtx["nround"],
                    max_depth = vtx["max_depth"], eta = vtx["eta"],
                    data= dat_x, label = dat_y, nfold = nfolds,
                    stratified = TRUE, prediction=pred, verbose=verb, missing = NA)$evaluation_log$test_mlogloss_mean[vtx["nround"]] )
                    #*****************************************************************************************************#
                    # -if using the version of xgboost in March 2016                                                      #
                    # stratified = TRUE, prediction=pred, verbose=verb, missing = NA)$test.mlogloss.mean[vtx["nround"]] ) #
                    #*****************************************************************************************************#
}

xgb_bdry_ex <- function(vtx = vector(0), ind_int = character(0), ind_num = character(0), min_int = integer(0)) {
  vtx[ind_int[1]] <- exp(vtx[ind_int[1]])
  vtx[ind_int[2]] <- exp(vtx[ind_int[2]])
  #to ensure that integer outputs are in their appropriate integer format
  vtx[ind_num] <- exp(vtx[ind_num])/(1+exp(vtx[ind_num]))
  vtx[ind_int][ vtx[ind_int] < min_int ] <- min_int
  vtx[ind_int] <- round(vtx[ind_int])
  return(vtx)
}

vtcs0 <- NM_opt( vtcs_init = cbind(nround = log(c(100,200,300,50)), 
                                  max_depth = log(c(7,3,2,20)), 
                                  eta = -log(1/c(0.2, 0.05, 0.09, 0.005) - 1) ),
                obj_fun = xgb_wrap_obj_ex,
                fxd_obj_param  = list(param = list("nthread" = 20,   # number of threads to be used 
                                                   "objective" = "multi:softmax",    # multi-class classification 
                                                   "eval_metric" ="mlogloss",    # evaluation metric
												                           "num_class" = 3
                                                   ),
                                      dat_x = dat_sparse,
                                      dat_y = dat_y ),
                bdry_fun = xgb_bdry_ex,
                fxd_bdry_param = list( ind_int = c("nround","max_depth"),
                                       ind_num = "eta", 
                                       min_int = 1 ), 
                a_e = 2.5,
                max_0prgrss = 50)


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo







#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#----------------------Examplary Application of Nelder Mead (randomForest)---------------------------------
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

library(caret)
K_Folds_K = 3
K_Folds = createFolds(dat_y, k = K_Folds_K, list = TRUE, returnTrain = FALSE)

library(doParallel)
cl <- makeCluster(K_Folds_K-1)
registerDoParallel(cl)



rf_wrap_ll_par <- function(vtx, ipt_target, ipt_data, idx_folds , ipt_ntrees = 100, ipt_replace = TRUE, ipt_do.trace = FALSE )	{
  k_folds <- length(idx_folds)
  ll_ttl <- foreach(i_cv = 1:k_folds, .combine=sum, .packages=c("randomForest")) %dopar% {
    N_i <- length(ipt_target) - length(idx_folds[[i_cv]])
    rf <- randomForest( ipt_target[-idx_folds[[1]]]~., data = ipt_data[-idx_folds[[1]],], importance=FALSE,
                        ntrees = vtx["ntrees"], replace = ipt_replace, mtry = vtx["mtry"], do.trace = ipt_do.trace )
    prob1 <- matrix(predict(rf, ipt_data[-idx_folds[[1]],],type = "prob", predict.all=TRUE)$aggregate, ncol=3)
    return( -sum(log(prob1[matrix(as.logical(model.matrix( ~ipt_target[-idx_folds[[1]]]-1, sparse = TRUE)),ncol=3)]), na.rm = TRUE) / N_i)
  }
  return(ll_ttl/k_folds)
}

vtcs0 <- NM_opt( vtcs_init = cbind(ntrees = log(c(50,200,400)),
                                   mtry = log(c(7,20,5))),
                 obj_fun = rf_wrap_ll_par,
                 fxd_obj_param  = list( ipt_target = factor(dat_y),
                                        ipt_data = dat_mod,
                                        idx_folds = K_Folds),
                 bdry_fun = rf_bdry,
                 fxd_bdry_param = list( ind_int = c("ntrees", "mtry"), min_int = 1),
                 a_e = 2.5,
                 max_0prgrss = 30)

#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo








#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#----------------------Examplary Calling of Nelder Mead Function for optimization on 8 paramaters --------
# [1] Does not work well for this example due to the small data size which causes problems optimizing  
#     xgboost on 'colsample_bytree' and 'colsample_bylevel'.
# [2] However, the following wrapper function can be easily modified to take in your own data
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo.oooooo
# 
# 
# xgb_wrap_obj3 <- function(param_fixed, vtx, dat_x, dat_y, nfolds = 3, pred = FALSE, verb = FALSE ){ 
#   return(   xgb.cv( param=param_fixed, 
#                     nround = vtx["nround"],
#                     max_depth = vtx["max_depth"], 
#                     eta = vtx["eta"],
#                     gamma = vtx["gamma"],
#                     subsample = vtx["subsample"],
#                     colsample_bytree = vtx["colsample_bytree"], 
#                     colsample_bylevel = vtx["colsample_bylevel"], 
#                     min_child_weight = vtx["min_child_weight"],
#                     data= dat_x, label = dat_y, nfold = nfolds, 
#                     stratified = TRUE, prediction=pred, verbose=verb, missing = NA)$evaluation_log$test_mlogloss_mean[vtx["nround"]] )
#                     #*****************************************************************************************************#
#                     # -if using the version of xgboost in March 2016                                                      #
#                     # stratified = TRUE, prediction=pred, verbose=verb, missing = NA)$test.mlogloss.mean[vtx["nround"]] ) #
#                     #*****************************************************************************************************#
# }
# 
# 
# xgb_bdry3 = function(vtx = vector(0), ind_int = character(0), ind_num = character(0), min_int = integer(0)) {
#   #mapping integer outputs
#   vtx[ind_int[1]] <- exp(vtx[ind_int[1]])             #nround [1,inf)
#   vtx[ind_int[2]] <- exp(vtx[ind_int[2]])             #max_depth [1,inf)
#   vtx[ind_int][ vtx[ind_int] < min_int ] <- min_int   
#   vtx[ind_int] <- round(vtx[ind_int])
#   #mapping numeric outputs
#   vtx[ind_num[1]] <- exp(vtx[ind_num[1]])/(1+exp(vtx[ind_num[1]]))  #eta [0,1]
#   vtx[ind_num[2]] <- exp(vtx[ind_num[2]])                           #gamma [0,inf)
#   vtx[ind_num[3]] <- exp(vtx[ind_num[3]])/(1+exp(vtx[ind_num[3]]))  #subsample (0,1]
#   vtx[ind_num[4]] <- exp(vtx[ind_num[4]])/(1+exp(vtx[ind_num[4]]))  #colsample_bytree (0,1]
#   vtx[ind_num[5]] <- exp(vtx[ind_num[5]])/(1+exp(vtx[ind_num[5]]))  #colsample_bylevel (0,1]
#   vtx[ind_num[6]] <- exp(vtx[ind_num[6]])                           #min_child_weigth [0,inf)
#   return(vtx)
# }
# 
# 
# # Data set too small to run a sufficiently big enough variety of parameters
# vtcs0 <- NM_opt( vtcs_init = cbind( nround = log(c(200,300,400,250,90,410,190,440,420)), 
#                                     max_depth = log(c(14,14,6,10,4,9,5,6,8)),
#                                     eta = -log(1/c(0.2, 0.05, 0.09, 0.005,0.3,0.02,0.5,0.3,0.7) - 1 +1e-05),
#                                     gamma = log(c(0.01,10,2,5,5,0.5,2,1,7)),
#                                     subsample = -log(1/runif(9) - 1 +1e-05),
#                                     colsample_bytree = -log(1/c(0.7, 0.1, 0.3, 0.15,0.8,0.14,0.7,0.3,0.5) - 1 +1e-05),
#                                     colsample_bylevel = -log(1/c(0.5, 0.7, 0.9, 0.4,0.9,0.5,0.7,0.5,0.5) - 1 +1e-05),
#                                     min_child_weight = log(c(4,3,4,9,2,10,3,1,7)) ),
#                                     obj_fun = xgb_wrap_obj3,
#                 fxd_obj_param  = list(param = list("nthread" = 20,   # number of threads to be used
#                                                    "objective" = "multi:softmax",    # multi-class classification 
#                                                    "eval_metric" ="mlogloss",    # evaluation metric
#                                                    "num_class" = 3
#                                                     ),
#                                     dat_x = dat_sparse,
#                                     dat_y = dat_y ),
#                 bdry_fun = xgb_bdry3,
#                 fxd_bdry_param = list(ind_int = c("nround","max_depth"),
#                                       ind_num = c("eta","gamma", "subsample", "colsample_bytree",
#                                                   "colsample_bylevel", "min_child_weight"),
#                                       min_int = 1 ),
#                 a_e = 2.7,
#                 max_0prgrss = 20
# )
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo