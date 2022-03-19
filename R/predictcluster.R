


#' Title predictcluster
#'
#' @param blcadata blcadata is a gene expression matrix form rna-seq or array. rowname is gene name and colname is sample's name. The value need log2(data+1) process before predicted.
#'
#' @return
#' @export
#'
#' @examples cluster <- predictcluster(test)
predictcluster <- function(blcadata){
   library(openxlsx)
   library(Matrix)
   library(impute)
   library(limma)
   library(sva)
   library(rpart)
   library(partykit)
   library(randomForest)
   library(xgboost)




  rt <- train
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)

  train <- data
  train=as.matrix(train)



  rt <- blcadata
  rt=as.matrix(rt)
  exp=rt[,1:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)

  test <- data
  test=as.matrix(test)



  sameSite=intersect(row.names(train),row.names(test))
  trainout=train[sameSite,]
  testout=test[sameSite,]

  batchType=c(rep(1,ncol(train)), rep(2,ncol(test)))
  all=cbind(trainout,testout)
  outTab=ComBat(all, batchType, par.prior=TRUE)

  trainout=outTab[,(1:ncol(trainout))]

  testout=outTab[,((ncol(trainout)+1):ncol(all))]



  ######准备训练集文


  sameSample=intersect(row.names(trainout),row.names(marker))
  trainout=trainout[sameSample,]




  #########xgboost


  #读取训练集输入文

  data=t(trainout)


  sameSample=intersect(row.names(data),row.names(cli))
  data=data[sameSample,]
  cli=cli[sameSample,]


  group <- cli
  group<-as.factor(group)

  letter=c(0,1,2,3,3,5,6)
  uniqClu=levels(group)
  group=letter[match(group, uniqClu)]

  data <- as.matrix(data)
  datatt <- data
  data<- Matrix(data,sparse=T)
  dtrain <- xgb.DMatrix(data = data, label = group)

  #输出重要基因的表达量
  sigExp1=t(datatt)
  sigExp <- as.matrix(sigExp1)


  #读取test集输入文
  data2<- testout

  sameSample3=intersect(row.names(sigExp1),row.names(data2))
  data2=data2[sameSample3,]

  data2=t(data2)



  data2 <- as.matrix(data2)
  data2 <- Matrix(data2,sparse=T)
  dtest <- xgb.DMatrix(data = data2)
  #预测xgbpre的建
  sigExp <- t(sigExp)
  sigExp <- Matrix(sigExp,sparse=T)
  dtrain2 <- xgb.DMatrix(data = sigExp, label = group)

  xgbpre <- xgboost( data = dtrain2,max_depth=6  #树的深度自己
                     , eta=0.4 #调节每棵树的权重，越小则约避免overfit
                     ,booster="gbtree" #gblinear#二分类，线性回
                     , objective='multi:softmax', num_class=3 #二分类，线性回 需要修
                     ,nround=100)




  pre_test<-predict(xgbpre,newdata = dtest)

  ##输出


  letter2=c("C1","C2","C3","C4","C5","C6","C7")
  uniqClu2=levels(factor(pre_test))
  pre_test=letter2[match(pre_test, uniqClu2)]


  out <- cbind(ID=rownames(data2),Cluster=pre_test)


  print(out)}
