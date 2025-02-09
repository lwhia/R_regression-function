#Build a multiple linear regression function, use the matrix method.

Proregression <- function(dataset, outcome, predictors){
  dataset <- data.frame(dataset)
  
  #define the original value of predictors' matrix, including the 1 as "offset"
  xmatrix <- data.frame(intercept =rep(1,nrow(dataset)))
  y1 <- dataset[[outcome]]
  
  #go through each variable name of predictors, if it exist in dataset, 
  #then add it's value to x matrix, 
  for (variablename in predictors){
    
    if (variablename %in% colnames(dataset)){
      xmatrix[[variablename]] <- dataset[[variablename]]} 
    else {print (paste("Predictor", variablename, "not found in dataset"))}  }

  #use the predictors and outcome to generate a new dataset2 (remove not provided predictors)
  dataset2 <- cbind(xmatrix,y1)

  #remove missing values if needed, and add warning message
  if (any(is.na(dataset2))){
    rowsnum1 <- nrow(dataset2)
    dataset3 <- na.omit(dataset2)
    rowsnum2 <- nrow(dataset3)
    diff <- rowsnum1-rowsnum2
    dataset2 <- dataset3
    print (paste(diff,"Rows with Missing values in predictors or outcome have been removed during modeling. If you want to include them, please impute data and provide new dataset"))  }
  
  y <- dataset2$y1
  x <- as.matrix(dataset2[, -ncol(dataset2)])

  #below is to calculate the regression result, eg. beta estimate, SE, RSS, t statistic, etc.
  x <- as.matrix(x)
  y <- as.matrix(y)
  xtransx <- t(x) %*% x
  betamatrix <- solve(xtransx) %*% t(x) %*% y
  betamatrix
  yhat <- x %*% betamatrix
  SYY <- sum((y - mean(y))^2)
  RSS <- t(y) %*% y - t(betamatrix) %*% t(x) %*% y 
  RSSeach <- y-yhat
  variance <- RSS/(nrow(dataset2)-length(predictors)-1)
  varbeta <- sqrt(diag((as.numeric(variance) * solve(xtransx))))
  tvalue <- betamatrix/varbeta
  degree <-  nrow(dataset2) - length(predictors) - 1
  pvalue <- 2 * (1 - pt(abs(tvalue),degree))
  
  #below is to calculate R square
  RSQUARE <- round(1-RSS/SYY,5)
  adjustedRsquare <- round(((nrow(x)-1)*RSQUARE-length(predictors))/(nrow(x)-length(predictors)-1),5)
  
  #below is to calculate the model's F statistic and p value, to see the whole model's fitness
  SSR <- SYY-RSS
  MSR <- SSR/length(predictors)
  MSE <- RSS/degree
  fstatistic <- round(MSR/MSE,5)
  fpvalue <- pf(fstatistic,df1=length(predictors),df2=degree,lower.tail = FALSE)
  
  #to put all the regresson results together in a table
  output <- data.frame(Estimate = round(as.numeric(betamatrix),5),StdError = round(varbeta,5),T_Statistic= round(tvalue,5),P_value= format(round(pvalue,5),scientific=TRUE))
  
  #below is the display output: residuals, R square, F statistic,  coefficient estimate, 
  print("Residuals Summary")
  print (summary(RSSeach))
  print (paste("R square is", RSQUARE, ",    Adjusted R square is", adjustedRsquare))
  
  print(paste("F Statistic:",fstatistic,  "on ",length(predictors), "and ", degree, "Degrees of freedom", "p-value is ",format(fpvalue,scientific=TRUE)))
  print("Coefficients:")
  return(output)
}

#Test the regression result
testfile <- read.csv("C:/Users/wangx/iCloudDrive/BS803 Statisitical programming (R)/final project/gas_vapor_with_missing.csv")
#testfile <- read.csv("gas_vapor_with_missing.csv")
attach(testfile)
#below is the function's output
Proregression(testfile,"y",c("x1","x3","x4"))

#below is the R default regression output
model1 <- lm(y~x1+x3+x4,data=testfile)
summary(model1)
testfile2 <- read.csv("C:/Users/wangx/iCloudDrive/BS852mthods in epi/project/framdat4.csv")
Proregression(testfile2,"FVC4",c("CHOL4","AGE4"))
model2 <- lm(FVC4~CHOL4+AGE4,data=testfile2)
summary(model2)