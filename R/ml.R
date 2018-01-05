#' @export
cvgroup <- function(data,cvnum,seed=1234){
  set.seed(seed)
  cv <- list()
  data <-  sample(data,length(data),replace = F)
  for(i in 1:cvnum){
    datastart <- (i-1)*round(length(data)/cvnum,0)+1
    dataend <- (i)*round(length(data)/cvnum,0)
    if(i==cvnum){dataend <- length(data)}
    a <- data[datastart:dataend]
    cv[[i]] <- list(train=setdiff(data,a),test=a)

  }
  cv
}

.ml_m5c <- function(cvnum=5,repeatTimes=5,posnames,negnames,pos_sample,neg_sample,
                   idependent_pos,idependent_neg,times=5,modeltype,over_sampling=F, perc.over = 500,perc.under=120,...){


  cv_pos <- cvgroup(data = posnames,cvnum = cvnum)
  cv_neg <- cvgroup(data = negnames,cvnum = cvnum)

  if(cvnum==1){
    cv_pos[[1]]$train <- posnames
    cv_neg[[1]]$train <- negnames
    cv_pos[[1]]$test <- idependent_pos
    cv_neg[[1]]$test <- idependent_neg
  }
  res_cv <- list()
  for (j in 1:cvnum){
  cat("cv part",j,"\t")
    #require(parallel)
    #cl <- makeCluster(cores, type="FORK")
    #fun <- function(j){
    require(randomForest)
    res_cv[[j]] <- list()

    for (i in 1:repeatTimes){
      pos <- cv_pos[[j]]$train
      set.seed(i)
      neg <- sample(cv_neg[[j]]$train,length(pos)*times)

      fmat <- rbind(pos_sample[c(pos),],neg_sample[c(neg),])
      label <- as.factor(c(rep(1,length(pos)),rep(0,length(neg))))
      pp <- data.frame(fmat,label)

      pos_test <- cv_pos[[j]]$test
      set.seed(1234)
      neg_test <- sample(cv_neg[[j]]$test,length(pos_test)*times)
      ################ over_sampling ##########
      if(over_sampling==T){
        newData <- SMOTE(label~., pp, perc.over = perc.over,perc.under=perc.under,...)
      }else{
        newData <- pp #SMOTE(label~., pp, perc.over = 500,perc.under=120)
      }
      ################ machine learning #########
      if(modeltype=="SVC"){
        model <-  svm(label~.,newData, probability = TRUE)
        pre_pos_test <- attr( predict(model,pos_sample[c(pos_test),],probability = T),"probabilities")[,1]
        pre_neg_test <- attr( predict(model,neg_sample[c(neg_test),],probability = T),"probabilities")[,1]

        pre_pos_idenpendnt <- attr( predict(model,pos_sample[c(idependent_pos),],probability = T),"probabilities")[,1]
        pre_neg_idenpendnt <- attr( predict(model,neg_sample[c(idependent_neg),],probability = T),"probabilities")[,1]
      }
      if(modeltype=="RFC"){
        #require(RWeka)
        #require(RWekajars)
        require(randomForest)

        # RT <- make_Weka_classifier("weka/classifiers/trees/RandomForest")
        newData <- as.data.frame(newData)
        # model <- RT(label ~ .,data=newData)
        model <-  randomForest(label~., newData,...)
        pre_pos_test <-predict(model,as.data.frame(pos_sample[c(pos_test),]),type=c("vote"))[,2]
        pre_neg_test <-predict(model,as.data.frame(neg_sample[c(neg_test),]),type=c("vote"))[,2]
        pre_pos_idenpendnt <-  predict(model,as.data.frame(pos_sample[c(idependent_pos),]),type=c("vote"))[,2]
        pre_neg_idenpendnt <- predict(model,as.data.frame(neg_sample[c(idependent_neg),]),type=c("vote"))[,2]
      }
      #print("modeling")
      # if(modeltype=="RFCE"){
      #   print(dim(newData))
      #
      #   model1 <-  randomForest(label~., newData[,c(1:164,271)])
      #   model2 <-  randomForest(label~., newData[,c(165: 248,271)])
      #   model3 <-  randomForest(label~., newData[,c(248:270,271)])
      #   model4 <-  randomForest(label~., newData[,c(70:100,271)])
      #   model5 <-  randomForest(label~., newData[,c(1:270,271)])
      #
      #   testdata <- rbind(pos_sample[c(pos_test),],neg_sample[c(neg_test),])
      #   print("ensembling")
      #   type1_1 <- predict(model1,testdata[,c(1:164)],type = "vote")[,2]
      #   type1_2 <- predict(model2,testdata[,c(165: 248)],type = "vote")[,2]
      #   type1_3 <- predict(model3,testdata[,c(248:270)],type = "vote")[,2]
      #   type1_4 <- predict(model4,testdata[,c(70:100)],type = "vote")[,2]
      #   type1_5 <- predict(model5,testdata[,c(1:270)],type = "vote")[,2]
      #
      #   type1 <- data.frame(type1_1,type1_2,type1_3,type1_4,type1_5,
      #                       label=as.factor(c(rep(1,length(pos_test)),rep(0,length(neg_test)))))
      #   model_type1 <-  randomForest(label~., type1)
      #
      #   pre_pos_test <- predict(model_type1,type1[1:length(pos_test),-6],type = "vote")[,2]
      #   pre_neg_test <-predict(model_type1,type1[-c(1:length(pos_test)),-6],type = "vote")[,2]
      #   ###########################################
      #   testdata <- rbind(pos_sample[c(idependent_pos),],neg_sample[c(idependent_neg),])
      #
      #   type1_1 <- predict(model1,testdata[,c(1:164)],type = "vote")[,2]
      #   type1_2 <- predict(model2,testdata[,c(165: 248)],type = "vote")[,2]
      #   type1_3 <- predict(model3,testdata[,c(248:270)],type = "vote")[,2]
      #   type1_4 <- predict(model4,testdata[,c(70:100)],type = "vote")[,2]
      #   type1_5 <- predict(model5,testdata[,c(1:270)],type = "vote")[,2]
      #
      #
      #
      #   pre_pos_idenpendnt <-  predict(model_type1,data.frame(type1_1,type1_2,type1_3,type1_4,type1_5)[1:length(idependent_pos),],
      #                                  type = "vote")[,2]
      #   pre_neg_idenpendnt <- predict(model_type1,data.frame(type1_1,type1_2,type1_3,type1_4,type1_5)[-c(1:length(idependent_pos)),],
      #                                 type = "vote")[,2]
      # }
      # if(modeltype=="RFCE2"){
      #   print(dim(newData))
      #
      #   model1 <-  randomForest(label~., newData[,c(1:164,271)])
      #   model2 <-  randomForest(label~., newData[,c(165: 248,271)])
      #   model3 <-  randomForest(label~., newData[,c(248:270,271)])
      #   model4 <-  randomForest(label~., newData[,c(70:100,271)])
      #   model5 <-  randomForest(label~., newData[,c(1:270,271)])
      #
      #   print("ensembling")
      #   type1_1 <- predict(model1,newData[,c(1:164)],type = "vote")[,2]
      #   type1_2 <- predict(model2,newData[,c(165: 248)],type = "vote")[,2]
      #   type1_3 <- predict(model3,newData[,c(248:270)],type = "vote")[,2]
      #   type1_4 <- predict(model4,newData[,c(70:100)],type = "vote")[,2]
      #   type1_5 <- predict(model5,newData[,c(1:270)],type = "vote")[,2]
      #
      #   type1 <- data.frame(type1_1,type1_2,type1_3,type1_4,type1_5,
      #                       label=newData[,c(271)])
      #   model_type1 <-  randomForest(label~., type1)
      #   #################################### test
      #   testdata <- rbind(pos_sample[c(pos_test),],neg_sample[c(neg_test),])
      #
      #   type1_1 <- predict(model1,testdata[,c(1:164)],type = "vote")[,2]
      #   type1_2 <- predict(model2,testdata[,c(165: 248)],type = "vote")[,2]
      #   type1_3 <- predict(model3,testdata[,c(248:270)],type = "vote")[,2]
      #   type1_4 <- predict(model4,testdata[,c(70:100)],type = "vote")[,2]
      #   type1_5 <- predict(model5,testdata[,c(1:270)],type = "vote")[,2]
      #
      #   type1 <- data.frame(type1_1,type1_2,type1_3,type1_4,type1_5)
      #
      #
      #   pre_pos_test <- predict(model_type1,type1[1:length(pos_test),],type = "vote")[,2]
      #   pre_neg_test <-predict(model_type1,type1[-c(1:length(pos_test)),],type = "vote")[,2]
      #   ########################################### idependent
      #   testdata <- rbind(pos_sample[c(idependent_pos),],neg_sample[c(idependent_neg),])
      #
      #   type1_1 <- predict(model1,testdata[,c(1:164)],type = "vote")[,2]
      #   type1_2 <- predict(model2,testdata[,c(165: 248)],type = "vote")[,2]
      #   type1_3 <- predict(model3,testdata[,c(248:270)],type = "vote")[,2]
      #   type1_4 <- predict(model4,testdata[,c(70:100)],type = "vote")[,2]
      #   type1_5 <- predict(model5,testdata[,c(1:270)],type = "vote")[,2]
      #
      #
      #
      #   pre_pos_idenpendnt <-  predict(model_type1,data.frame(type1_1,type1_2,type1_3,type1_4,type1_5)[1:length(idependent_pos),],
      #                                  type = "vote")[,2]
      #   pre_neg_idenpendnt <- predict(model_type1,data.frame(type1_1,type1_2,type1_3,type1_4,type1_5)[-c(1:length(idependent_pos)),],
      #                                 type = "vote")[,2]
      # }
      require(PRROC)
      pr <- roc.curve(pre_pos_test,pre_neg_test,curve = T)

      #hist(a[c(1:22)])
      #hist(a[-c(1:22)])
      res_cv[[j]][[i]] <- list(pos_train=pos,neg_train=neg,model=model,
                               positives.test=pos_test,negatives.test=neg_test,
                               positives.test.score=pre_pos_test,negatives.test.score=pre_neg_test,auc= pr$auc,
                               pre_pos_idenpendnt=pre_pos_idenpendnt,
                               pre_neg_idenpendnt=pre_neg_idenpendnt)
    }
    cvAUCmax <- which.max(sapply(res_cv[[j]][1:repeatTimes],function(x){a<-x[["auc"]];a}))
    posscore1 <- apply(sapply(res_cv[[j]],function(x) x$pre_pos_idenpendnt), 1, function(x) mean(x))#x[cvAUCmax])
    negscore1 <- apply(sapply(res_cv[[j]],function(x) x$pre_neg_idenpendnt), 1, function(x) mean(x))#x[cvAUCmax])

    posscore2 <- apply(sapply(res_cv[[j]],function(x) x$positives.test.score), 1, function(x) mean(x))#x[cvAUCmax])
    negscore2 <- apply(sapply(res_cv[[j]],function(x) x$negatives.test.score), 1, function(x) mean(x))#x[cvAUCmax])

    res_cv[[j]][["positives.test.score.id"]] <- posscore1
    res_cv[[j]][["negatives.test.score.id"]] <- negscore1

    res_cv[[j]][["positives.test.score"]] <- posscore2
    res_cv[[j]][["negatives.test.score"]] <- negscore2

    res_cv[[j]][["positives.test"]] <- res_cv[[j]][[cvAUCmax]]$positives.test
    res_cv[[j]][["negatives.test"]] <- res_cv[[j]][[cvAUCmax]]$negatives.test

    res_cv[[j]][["auc_test"]] <- max(sapply(res_cv[[j]][1:repeatTimes],function(x){a<-x[["auc"]];a}))
    res_cv[[j]][["auc_test_id"]] <- roc.curve(posscore1,negscore1,curve = T)$auc
    res_cv[[j]]
  }
  #res_cv <- parLapply(cl, 1:cvnum,  fun)
  #stopCluster(cl)

  return(res_cv)
}




############ extra ig res ####
#' @export
PEA_ml <- function(pos_sample,neg_sample,independent_num=100,ig="ALL",
                   times = 1,modeltype = "RFC",cvnum = 5,repeatTimes = 1, ntree=200,over_sampling = F){

  Root_cdna_pos <- rownames(pos_sample)
  #set.seed(1)
  #neg_sample <- neg_sample[sample(rownames(neg_sample),nrow(pos_sample)),]

  set.seed(1234)
  idependent_Root <- list()
  idependent_Root[["pos"]] <- sample(Root_cdna_pos,independent_num)
  set.seed(1234)
  idependent_Root[["neg"]] <- sample(rownames(neg_sample),independent_num)
  ################train data
  train_data <- list()
  train_data[["pos"]] <- setdiff(Root_cdna_pos,idependent_Root[["pos"]])
  train_data[["neg"]] <- setdiff(rownames(neg_sample),idependent_Root[["neg"]])


  ################
  require(randomForest)
  #############  feature object
  print("")
  R <- seq(0.2,1,by = 0.2)

  res_cv <- list()
  alldata <- rbind(pos_sample[train_data[["pos"]],],neg_sample[train_data[["neg"]],])
  alldata <- cbind(alldata,c(rep(1,nrow(pos_sample[train_data[["pos"]],])),rep(0,nrow(neg_sample[train_data[["neg"]],]))))
  colnames(alldata)[ncol(alldata)] <- "label"
  alldata <- as.data.frame(alldata)
  alldata$label <- as.factor(alldata$label)
  ############# information.gain correlation
  cat("information.gain correlation \n")
  require(FSelector)
  weights <- information.gain(label~., alldata)
  #feature_object <- cutoff.k.percent(weights,0.05)
  feature_num <- ncol(pos_sample)
  AUC <- c()
  length(AUC) <- feature_num
  ###############
  if (ig == "ALL"){
    first_num <- ncol(alldata)-1
  }else{
    if (ig=="gradient") {
      first_num <- 3:(ncol(alldata)-1)
    }else{
      first_num <- ig
    }
  }
  cat("modeling and prediction with different features based information gain rank \n")
  for ( i in first_num){
    feature_object <- cutoff.k(weights,i)
    res_cv[[i]] <- .ml_m5c(cvnum = cvnum,cores=5,repeatTimes = repeatTimes,posnames = train_data[["pos"]] ,negnames =train_data[["neg"]]  ,
                           neg_sample = neg_sample[,feature_object],pos_sample = pos_sample[,feature_object],idependent_pos = idependent_Root[["pos"]],
                           idependent_neg = idependent_Root[["neg"]], ntree=ntree,times = times,modeltype = modeltype,over_sampling = over_sampling)
    names(res_cv[[i]]) <- paste0("cv_",1:cvnum)
    AUC[i] <- mean(sapply(res_cv[[i]],function(x) x$auc_test))
    cat("\n This part AUC mean value is ",AUC[i],"(cv test).",'\n')
    cat("res of information gain with",i,"feature over",'\n')
  }

  AUC <- (cbind(first_num,AUC))
  #
  # save(weights,res_cv,train_data,idependent_Root,neg_sample,AUC,first_num,second_num,third_num,
  #      file = paste0("../res/",deal,"information"))
  res_cv <- res_cv[first_num]
  names(res_cv) <- paste0("ig_",first_num)
  res_cv[["weights"]] <- weights
  res_cv[["AUC"]] <- AUC

  res_cv
}

# load("/home/malab5/R/PEAm5c/root_36_46_sample.RData")
# aaa <- PEA_ml(pos_sample = pos_sample,neg_sample = neg_sample, ig=50 )

############ extra model ####
#' @export
extra_model <- function(res,ignum=150) {
  models <- res[[paste0("ig_",ignum)]]
  weights <- res$weights
  length(models)
  res <- list()
  for(i in 1:length(models)){
    res[[i]] <- models[[i]][[1]]$model
  }
  require(FSelector)
  res[["weights"]] <- cutoff.k(attrs = weights,k = ignum)
  return(res)
}
#
# load("/home/malab5/R/PEAm5c/root_36_46_sample.RData")
# aaa <- PEA_ml(pos_sample = pos_sample,neg_sample = neg_sample)
# ddd <- extra_model(res = aaa,ignum=50)
# ddd

####### self-defined prediction ###
#' @export
predict_self_model <- function(models,sequence_dir,end = 5,up = 5){

  aaa <- extra_motif_seq(input_seq_dir = sequence_dir,end = end,up = up)

  aaa <- lapply(aaa, c2s)
  bbb <- FeatureExtract(aaa)

  res <- 0
  for(i in 1:(length(models)-1)){
    alldata <- bbb
    alldata <- as.data.frame(alldata)
    model <-  models[[i]]
    res <- predict(model,alldata[,models$weights],type=c("vote"))[,2]/(length(models)-1)+res
    cat("predicting part ",i," ...\t")

  }
  res_table <- cbind(str_split(names(res),"_",simplify = T),res,"non-CMR")
  res_table[which(res>0.5),4] <- "CMR"
  colnames(res_table) <- c("transcript","position","score","res")
  cat("over")
  return(res_table)
}
# load("/home/malab5/R/PEAm5c/root_36_46_sample.RData")
# aaa <- PEA_ml(pos_sample = pos_sample,neg_sample = neg_sample)
# ddd <- extra_model(res = aaa)
# ddd
#
# eee <- predict_self_model(models = ddd,sequence_dir = "~/work/ArabidopsisMethylation_m5c/newdata/all/cdna1_6.fa")
# table(eee[,4])
