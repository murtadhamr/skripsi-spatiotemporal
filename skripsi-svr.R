library(ggplot2)
library(e1071)
library(caTools)
library(FNN)
library(Metrics)
library(caret)
library(RSNNS)
library(EnvStats)
library(car)
library(RANN)
#read dataset
dataset = read.csv('dataset22-26(fix).csv')

#checking missing value
sum(is.na(dataset))

#normalizing dataset into range 0-1
dataset_normalized = dataset

dataset_normalized$CO201500 = normalizeData(dataset$CO201500, type='0_1')
dataset_normalized$CO01500 = normalizeData(dataset$CO01500, type='0_1')

#####     Spatial Model     #####
#####     CO2       #####
#finding 8 nearest point from will-be-predicted point
## nearest = knn.index(data = dataset[, 3:6],  k=8, algorithm = "kd_tree")
nearest = nn2(data = dataset[, 3:6], k = min(9, nrow(dataset)), treetype = "kd", searchtype = "radius", radius = 0.5)
nearest_index = nearest[["nn.idx"]]

##check which one has value 0 or no neighbours
less_neighbour = c()
blank = 0


for(j in 1:nrow(nearest_index)){
  for(i in 1:ncol(nearest_index)){
    if(nearest_index[j, i] == 0)
      blank = blank + 1
  }
  if(blank > 3)
    less_neighbour = append(less_neighbour, nearest_index[j, 1])
  blank = 0
}

no_neighbour = as.data.frame(matrix(nrow = length(less_neighbour), ncol = 9))

for(j in 1:length(less_neighbour)){
  for(i in 1:9){
    no_neighbour[j, i] = nearest_index[less_neighbour[j], i]
  }
  nearest_index[less_neighbour[j], 1] = 0
}

nearest_index_fix = subset(nearest_index, nearest_index[, 1] > 0)

mean_co2 = mean(dataset_normalized$CO201500)
mean_co = mean(dataset_normalized$CO01500)

#making blank dataframe to store concentration value from nearest point(by index)
datafull_co2 = as.data.frame(matrix(nrow = nrow(nearest_index_fix), ncol = 9))

#storing nearest concentration value to dataset2
#interpolation with mean if index == 0
for (i in 1:9){
  for (j in 1: nrow(nearest_index_fix)){
    if(nearest_index_fix[j, i] == 0)
      datafull_co2[j, i] = mean_co2
    else
      datafull_co2[j, i] = dataset_normalized[, 7][nearest_index_fix[j, i]]
  }
}

id = nearest_index_fix[, 1]

datafull_co2 = cbind(id, datafull_co2)

#data sampling (0.8 train, 0.2 test)
set.seed(123)
split_co2 = sample.split(datafull_co2$V1, SplitRatio = 0.8)
training_set_co2 = subset(datafull_co2, split_co2 == TRUE)
test_set_co2 = subset(datafull_co2, split_co2 == FALSE)

# using SVR
regressor = svm(x = training_set_co2[, 3:10] ,
                y = training_set_co2[, 2],
                type = 'eps-regression',
                kernel = 'radial')

#storing predicted value in variable predictedY
predictedY1 = predict(regressor, test_set_co2[, 3:10])

test = test_set_co2$V1
write.csv(test, file = "test.csv")


#denormalize data
test_normalized = denormalizeData(test_set_co2$V1, getNormParameters(dataset_normalized$CO201500))
predictedY1_normalized = denormalizeData(predictedY1, getNormParameters(dataset_normalized$CO201500))

#Root Mean Squared Error of CO2 Spatial Prediction
predictionRMSE1 = rmse(test_normalized, predictedY1_normalized)

#correlation
correlation = cor(x = test_set_co2[, 2],
                  y = predict(regressor, newdata = test_set_co2[, 3:10]), 
                  method = c("pearson", "kendall", "spearman"))

#plotting the actual and predcited data
ggplot() +
  geom_line(aes(x = test_set_co2$id , y = test_normalized),
             colour = 'red') +
  geom_line(aes(x = test_set_co2$id , y = predictedY1_normalized),
             colour = 'blue') 

# perform a grid search
tuneResult2 = tune(svm, train.x = training_set_co2[, 3:10] ,
                  train.y = training_set_co2[, 2],
                  ranges = list(epsilon = seq(0.01,0.0001), cost = 1:10)
)
print(tuneResult2)
# Draw the tuning graph
plot(tuneResult2)

####
#https://www.kdnuggets.com/2017/03/building-regression-models-support-vector-regression.html

##### CO #####














