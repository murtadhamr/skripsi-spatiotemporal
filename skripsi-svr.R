#read dataset
dataset = read.csv('dataset22-26(fix).csv')

#selecting related attributes in building model
x = dataset[, 5:8]

#selecting c02 value from dataframe x
co2 = x[, 3]
lonlat = x[, 1:2]

#finding 9 nearest point from will-be-predicted point
nearest = knn.index(lonlat, k=9, algorithm = "kd_tree")

#making blank dataframe to store concentration value from nearest point(by index)
dataset2 = as.data.frame(matrix(nrow = 25789, ncol = 9))

#storing nearest concentration value to dataset2
for (i in 1:9){
  for (j in 1: 25789){
    dataset2[j, i] = co2[nearest[j, i]]
  }
}

#merging attributes with target
dataco2 = merge(dataset2, co2)

#using SVR
library(e1071)



