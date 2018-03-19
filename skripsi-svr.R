library(ggplot2)
library(e1071)
library(caTools)
library(FNN)
library(Metrics)
#read dataset
dataset = read.csv('dataset22-26(fix).csv')

#####     CO2       #####
#finding 9 nearest point from will-be-predicted point
nearest = knn.index(data = dataset[, 3:6],  k=9, algorithm = "kd_tree")

#making blank dataframe to store concentration value from nearest point(by index)
dataset2 = as.data.frame(matrix(nrow = 25789, ncol = 9))

#storing nearest concentration value to dataset2
for (i in 1:9){
  for (j in 1: 25789){
    dataset2[j, i] = dataset[, 7][nearest[j, i]]
  }
}

#storing features and target to var x and Y
x =  dataset2
Y = dataset[, 7]

#binding x and Y with its id
id = 1:25789
datafull = cbind(cbind(id, x),  Y)

#data sampling (0.8 train, 0.2 test)
set.seed(123)
split = sample.split(datafull$Y, SplitRatio = 0.8)
training_set = subset(datafull, split == TRUE)
test_set = subset(datafull, split == FALSE)

# using SVR
regressor = svm(x = training_set[, 2:10] ,
                y = training_set[, 11],
                type = 'eps-regression',
                kernel = 'radial')

predictedY = predict(regressor, test_set[, 2:10])

#plotting the actual and predcited data
ggplot() +
  geom_point(aes(x = test_set$id , y = test_set[, 11]),
             colour = 'red') +
  geom_point(aes(x = test_set$id , y = predict(regressor, newdata = test_set[, 2:10])),
             colour = 'blue') 
  
#Root Mean Squared Error of CO2 Spatial Prediction
predictionRMSE = rmse(test_set$Y, predict(regressor, newdata = test_set[, 2:10]))

##### CO #####


  


