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
y = co2

# using SVR
library(e1071)
regressor = svm(x = x, 
                y = y,
                type = 'eps-regression',
                kernel = 'sigmoid')


rmse = function(error)
{
  sqrt(mean(error^2))
}

error = regressor$residuals

predictionRMSE = rmse(error)



