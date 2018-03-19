#read dataset
dataset = read.csv('dataset22-26(fix).csv')

#selecting related attributes in building model
selected = dataset[, 5:8]

#selecting c02 value from dataframe x
co2 = selected[, 3]
lonlat = selected[, 1:2]

#finding 9 nearest point from will-be-predicted point
library(FNN)
nearest = knn.index(lonlat, k=9, algorithm = "kd_tree")

#making blank dataframe to store concentration value from nearest point(by index)
dataset2 = as.data.frame(matrix(nrow = 25789, ncol = 9))

#storing nearest concentration value to dataset2
for (i in 1:9){
  for (j in 1: 25789){
    dataset2[j, i] = co2[nearest[j, i]]
  }
}

x =  dataset2
Y = co2

id = 1:25789
datafull = cbind(cbind(id, x),  Y)

library(caTools)
set.seed(123)
split = sample.split(datafull$Y, SplitRatio = 0.8)
training_set = subset(datafull, split == TRUE)
test_set = subset(datafull, split == FALSE)

# using SVR
library(e1071)
regressor = svm(x = training_set[, 2:10] ,
                y = training_set[, 11],
                type = 'eps-regression',
                kernel = 'radial')

predictedY = predict(regressor, test[, 2:10])

library(ggplot2)
ggplot() +
  geom_point(aes(x = test_set$id , y = test_set[, 11]),
             colour = 'red') +
  geom_point(aes(x = test_set$id , y = predict(regressor, newdata = test_set[, 2:10])),
             colour = 'blue') 
  
  

library(Metrics)
predictionRMSE = rmse(test_set$Y, predict(regressor, newdata = test_set[, 2:10]))

pacf(co2, length(co2) - 1, pl = TRUE)

library(forecast)
library(tseries)
adf.test(co2, alternative = 'stationary', lag.max = 20)

#making blank dataframe to store concentration value from nearest point(by index)
dataset2 = as.data.frame(matrix(nrow = 25789, ncol = 30))
  


