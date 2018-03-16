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

## 75% of the sample size
smp_size <- floor(0.8 * nrow(x))

## set the seed to make your partition reproductible
set.seed(123)
train_ind <- sample(seq_len(nrow(x)), size = smp_size)

train <- x[train_ind, ]
test <- x[-train_ind, ]

# using SVR
library(e1071)
regressor = svm(x = x, 
                y = Y,
                type = 'eps-regression',
                kernel = 'radial')


rmse = function(error)
{
  sqrt(mean(error^2))
}

error = regressor$residuals

predictionRMSE = rmse(error)

predictedY = predict(regressor, x)

library(ggplot2)
ggplot() +
  geom_point(aes(x = datafull$id , y = datafull$Y),
             colour = 'red') +
  geom_point(aes(x = datafull$id , y = predict(regressor, newdata = x)),
            colour = 'blue')


id = 1:25789
datafull = cbind(cbind(id, x),  Y)

  


