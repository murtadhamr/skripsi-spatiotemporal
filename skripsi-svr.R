library(ggplot2)
library(e1071)
library(caTools)
library(Metrics)
library(RSNNS)
library(EnvStats)
library(RANN)
#read dataset
dataset = read.csv('dataset22-26(fix).csv')

#checking missing value
sum(is.na(dataset))

#normalizing dataset into range 0-1
dataset_normalized = dataset

dataset_normalized$CO201500 = normalizeData(dataset$CO201500, type='0_1')
dataset_normalized$CO01500 = normalizeData(dataset$CO01500, type='0_1')

#finding 8 nearest point from will-be-predicted point
## nearest = knn.index(data = dataset[, 3:6],  k=8, algorithm = "kd_tree")
nearest_point = nn2(data = dataset[, 3:6], k = min(9, nrow(dataset)), treetype = "kd", searchtype = "radius", radius = 0.45045)
nearest_index = nearest_point[["nn.idx"]]

##check which one has value 0 and not reach 70% minimum neighbors
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

#to store the qualified points that has less neighbours
less_neighbour_point = as.data.frame(matrix(nrow = length(less_neighbour), ncol = 9))


for(j in 1:length(less_neighbour)){
  for(i in 1:9){
    less_neighbour_point[j, i] = nearest_index[less_neighbour[j], i]
  }
  nearest_index[less_neighbour[j], 1] = 0
}

longlat_no_neighbour = as.data.frame(matrix(nrow = length(less_neighbour), ncol = 3))

#checking longlat of less neighbour points
for(j in 1:length(less_neighbour)){
  for(i in 1:3){
    if(i == 1)
      longlat_no_neighbour[j, i] = less_neighbour[j]
    else if(i == 2)
      longlat_no_neighbour[j, i] = dataset[j,5]
    else if(i == 3)
      longlat_no_neighbour[j, i] = dataset[j,6]
  }
}

nearest_index_fix = subset(nearest_index, nearest_index[, 1] > 0)

#row which should be interpolated
index_interpolated = as.data.frame(matrix(nrow = sum(nearest_index_fix[, 9] == 0), ncol = 9))

j = 1

for(i in 1:nrow(nearest_index_fix)){
  if(nearest_index_fix[i, 9] == 0){
    index_interpolated[j, ] = nearest_index_fix[i, ]
    j = j+1
  }
}

#long lat for determining random point position
long_lat_interpolation = {}

random_point_longlat = function(i, long_lat_interpolation, index_interpolated, dataset){
  long_lat_interpolation = as.data.frame(matrix(nrow = sum(index_interpolated[i, ]!=0), ncol = 2))
  #lat 1 long 2
  for(j in 1:nrow(long_lat_interpolation)){
    for(k in 1:2){
      #lat
      if(k == 1)
        long_lat_interpolation[j, k] = dataset[index_interpolated[i, j], 5]
      #long
      else if (k == 2)
        long_lat_interpolation[j, k] = dataset[index_interpolated[i, j], 6]
    }
  }
  return(list(min(long_lat_interpolation[, 1]), max(long_lat_interpolation[, 1]), min(long_lat_interpolation[, 2]), max(long_lat_interpolation[, 2])))
}

min_max_longlat = as.data.frame(matrix(nrow = nrow(less_neighbour_point), ncol = 4))

for(i in 1:nrow(index_interpolated)){
  min_max_longlat[i, ] = matrix(unlist(random_point_longlat(i, long_lat_interpolation, index_interpolated, dataset)), ncol = 4, byrow = TRUE)
}


trunc <- function(x, ..., prec = 0) base::trunc(x * 10^prec, ...) / 10^prec


random_point = as.data.frame(matrix(nrow = nrow(less_neighbour_point), ncol = 6))

for(i in 1:nrow(index_interpolated)){
  for(j in 1:ncol(random_point)){
    if(j == (sum(index_interpolated[i, ]==0)*2) || j < (sum(index_interpolated[i, ]==0)*2)){
      #lat
      if(j%%2 != 0)
        random_point[i, j] = trunc(runif(1, min_max_longlat[i, 1], min_max_longlat[i, 2]), prec = 3)
      #lon
      else if (j%%2 == 0)
        random_point[i, j] = trunc(runif(1, min_max_longlat[i, 3], min_max_longlat[i, 4]), prec = 3)
    }
    else
      break
  }
}

euclideanDistance = function(x, y, a, b){
  dist = sqrt(((x-a)^2)+((y-b)^2)) 
  return(dist)
}

## Spatial Interpolation ##

## id polutant 7 = co2, 8 = co
inverseDitanceWeight_3 = function(random_point, dataset_normalized, index_interpolated, row_interpolated, id_polutant){
  
    d1 = euclideanDistance(random_point[row_interpolated, 5], random_point[row_interpolated, 6], dataset_normalized[index_interpolated[row_interpolated, 2], 5], dataset_normalized[index_interpolated[row_interpolated, 2], 6])
    d2 = euclideanDistance(random_point[row_interpolated, 5], random_point[row_interpolated, 6], dataset_normalized[index_interpolated[row_interpolated, 3], 5], dataset_normalized[index_interpolated[row_interpolated, 3], 6])
    d3 = euclideanDistance(random_point[row_interpolated, 5], random_point[row_interpolated, 6], dataset_normalized[index_interpolated[row_interpolated, 4], 5], dataset_normalized[index_interpolated[row_interpolated, 4], 6])
    d4 = euclideanDistance(random_point[row_interpolated, 5], random_point[row_interpolated, 6], dataset_normalized[index_interpolated[row_interpolated, 5], 5], dataset_normalized[index_interpolated[row_interpolated, 5], 6])
    d5 = euclideanDistance(random_point[row_interpolated, 5], random_point[row_interpolated, 6], dataset_normalized[index_interpolated[row_interpolated, 6], 5], dataset_normalized[index_interpolated[row_interpolated, 6], 6])
  
  a1 = dataset_normalized[index_interpolated[row_interpolated, 2], id_polutant]
  a2 = dataset_normalized[index_interpolated[row_interpolated, 3], id_polutant]
  a3 = dataset_normalized[index_interpolated[row_interpolated, 4], id_polutant]
  a4 = dataset_normalized[index_interpolated[row_interpolated, 5], id_polutant]
  a5 = dataset_normalized[index_interpolated[row_interpolated, 6], id_polutant]
  sigma_wz = as.vector((a1/d1^2)+(a2/d2^2)+(a3/d3^2)+(a4/d4^2)+(a5/d5^2))
  sigma_w = as.vector((1/d1^2)+(1/d2^2)+(1/d3^2)+(1/d4^2)+(1/d5^2))
  return(sigma_wz/sigma_w)
}

inverseDitanceWeight_2_random = function(interpolated_point_normalized, random_point, dataset_normalized, index_interpolated, row_interpolated, id_polutant){
  
    d1 = euclideanDistance(random_point[row_interpolated, 3], random_point[row_interpolated, 4], dataset_normalized[index_interpolated[row_interpolated, 2], 5], dataset_normalized[index_interpolated[row_interpolated, 2], 6])
    d2 = euclideanDistance(random_point[row_interpolated, 3], random_point[row_interpolated, 4], dataset_normalized[index_interpolated[row_interpolated, 3], 5], dataset_normalized[index_interpolated[row_interpolated, 3], 6])
    d3 = euclideanDistance(random_point[row_interpolated, 3], random_point[row_interpolated, 4], dataset_normalized[index_interpolated[row_interpolated, 4], 5], dataset_normalized[index_interpolated[row_interpolated, 4], 6])
    d4 = euclideanDistance(random_point[row_interpolated, 3], random_point[row_interpolated, 4], dataset_normalized[index_interpolated[row_interpolated, 5], 5], dataset_normalized[index_interpolated[row_interpolated, 5], 6])
    d5 = euclideanDistance(random_point[row_interpolated, 3], random_point[row_interpolated, 4], dataset_normalized[index_interpolated[row_interpolated, 6], 5], dataset_normalized[index_interpolated[row_interpolated, 6], 6])
    d6 = euclideanDistance(random_point[row_interpolated, 3], random_point[row_interpolated, 4], random_point[row_interpolated, 5], random_point[row_interpolated, 6])
    #d6 = euclideanDistance(random_point[row_interpolated, 3], random_point[row_interpolated, 4], dataset_normalized[index_interpolated[row_interpolated, 7], 5], dataset_normalized[index_interpolated[row_interpolated, 7], 6])
    
    a1 = dataset_normalized[index_interpolated[row_interpolated, 2], id_polutant]
    a2 = dataset_normalized[index_interpolated[row_interpolated, 3], id_polutant]
    a3 = dataset_normalized[index_interpolated[row_interpolated, 4], id_polutant]
    a4 = dataset_normalized[index_interpolated[row_interpolated, 5], id_polutant]
    a5 = dataset_normalized[index_interpolated[row_interpolated, 6], id_polutant]
    a6 = interpolated_point_normalized[row_interpolated, 7]
    sigma_wz = as.vector((a1/d1^2)+(a2/d2^2)+(a3/d3^2)+(a4/d4^2)+(a5/d5^2)+(a6/d6^2))
    sigma_w = as.vector((1/d1^2)+(1/d2^2)+(1/d3^2)+(1/d4^2)+(1/d5^2)+(1/d6^2))
    return(sigma_wz/sigma_w)
}

inverseDitanceWeight_2_fix = function(interpolated_point_normalized, random_point, dataset_normalized, index_interpolated, row_interpolated, id_polutant){
  
  d1 = euclideanDistance(random_point[row_interpolated, 3], random_point[row_interpolated, 4], dataset_normalized[index_interpolated[row_interpolated, 2], 5], dataset_normalized[index_interpolated[row_interpolated, 2], 6])
  d2 = euclideanDistance(random_point[row_interpolated, 3], random_point[row_interpolated, 4], dataset_normalized[index_interpolated[row_interpolated, 3], 5], dataset_normalized[index_interpolated[row_interpolated, 3], 6])
  d3 = euclideanDistance(random_point[row_interpolated, 3], random_point[row_interpolated, 4], dataset_normalized[index_interpolated[row_interpolated, 4], 5], dataset_normalized[index_interpolated[row_interpolated, 4], 6])
  d4 = euclideanDistance(random_point[row_interpolated, 3], random_point[row_interpolated, 4], dataset_normalized[index_interpolated[row_interpolated, 5], 5], dataset_normalized[index_interpolated[row_interpolated, 5], 6])
  d5 = euclideanDistance(random_point[row_interpolated, 3], random_point[row_interpolated, 4], dataset_normalized[index_interpolated[row_interpolated, 6], 5], dataset_normalized[index_interpolated[row_interpolated, 6], 6])
  d6 = euclideanDistance(random_point[row_interpolated, 3], random_point[row_interpolated, 4], dataset_normalized[index_interpolated[row_interpolated, 7], 5], dataset_normalized[index_interpolated[row_interpolated, 7], 6])
  
  a1 = dataset_normalized[index_interpolated[row_interpolated, 2], id_polutant]
  a2 = dataset_normalized[index_interpolated[row_interpolated, 3], id_polutant]
  a3 = dataset_normalized[index_interpolated[row_interpolated, 4], id_polutant]
  a4 = dataset_normalized[index_interpolated[row_interpolated, 5], id_polutant]
  a5 = dataset_normalized[index_interpolated[row_interpolated, 6], id_polutant]
  a6 = interpolated_point_normalized[row_interpolated, 7]
  sigma_wz = as.vector((a1/d1^2)+(a2/d2^2)+(a3/d3^2)+(a4/d4^2)+(a5/d5^2)+(a6/d6^2))
  sigma_w = as.vector((1/d1^2)+(1/d2^2)+(1/d3^2)+(1/d4^2)+(1/d5^2)+(1/d6^2))
  return(sigma_wz/sigma_w)
}

inverseDitanceWeight_1_fix_fix = function(interpolated_point_normalized,random_point, dataset_normalized, index_interpolated, row_interpolated, id_polutant){
    d1 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 2], 5], dataset_normalized[index_interpolated[row_interpolated, 2], 6])
    d2 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 3], 5], dataset_normalized[index_interpolated[row_interpolated, 3], 6])
    d3 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 4], 5], dataset_normalized[index_interpolated[row_interpolated, 4], 6])
    d4 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 5], 5], dataset_normalized[index_interpolated[row_interpolated, 5], 6])
    d5 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 6], 5], dataset_normalized[index_interpolated[row_interpolated, 6], 6])
    #fix-fix
    d6 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 7], 5], dataset_normalized[index_interpolated[row_interpolated, 7], 6])
    d7 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 8], 5], dataset_normalized[index_interpolated[row_interpolated, 8], 6])
    
    a1 = dataset_normalized[index_interpolated[row_interpolated, 2], id_polutant]
    a2 = dataset_normalized[index_interpolated[row_interpolated, 3], id_polutant]
    a3 = dataset_normalized[index_interpolated[row_interpolated, 4], id_polutant]
    a4 = dataset_normalized[index_interpolated[row_interpolated, 5], id_polutant]
    a5 = dataset_normalized[index_interpolated[row_interpolated, 6], id_polutant]
    a6 = interpolated_point_normalized[row_interpolated, 7]
    a7 = interpolated_point_normalized[row_interpolated, 8]
    sigma_wz = as.vector((a1/d1^2)+(a2/d2^2)+(a3/d3^2)+(a4/d4^2)+(a5/d5^2)+(a6/d6^2)+(a7/d7^2))
    sigma_w = as.vector((1/d1^2)+(1/d2^2)+(1/d3^2)+(1/d4^2)+(1/d5^2)+(1/d6^2)+(1/d7^2))
    return(sigma_wz/sigma_w)
}

inverseDitanceWeight_1_fix_random = function(interpolated_point_normalized,random_point, dataset_normalized, index_interpolated, row_interpolated, id_polutant){
  d1 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 2], 5], dataset_normalized[index_interpolated[row_interpolated, 2], 6])
  d2 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 3], 5], dataset_normalized[index_interpolated[row_interpolated, 3], 6])
  d3 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 4], 5], dataset_normalized[index_interpolated[row_interpolated, 4], 6])
  d4 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 5], 5], dataset_normalized[index_interpolated[row_interpolated, 5], 6])
  d5 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 6], 5], dataset_normalized[index_interpolated[row_interpolated, 6], 6])
  #fix-random
  d6 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 7], 5], dataset_normalized[index_interpolated[row_interpolated, 7], 6])
  d7 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], random_point[row_interpolated, 3], random_point[row_interpolated, 4])

  a1 = dataset_normalized[index_interpolated[row_interpolated, 2], id_polutant]
  a2 = dataset_normalized[index_interpolated[row_interpolated, 3], id_polutant]
  a3 = dataset_normalized[index_interpolated[row_interpolated, 4], id_polutant]
  a4 = dataset_normalized[index_interpolated[row_interpolated, 5], id_polutant]
  a5 = dataset_normalized[index_interpolated[row_interpolated, 6], id_polutant]
  a6 = interpolated_point_normalized[row_interpolated, 7]
  a7 = interpolated_point_normalized[row_interpolated, 8]
  sigma_wz = as.vector((a1/d1^2)+(a2/d2^2)+(a3/d3^2)+(a4/d4^2)+(a5/d5^2)+(a6/d6^2)+(a7/d7^2))
  sigma_w = as.vector((1/d1^2)+(1/d2^2)+(1/d3^2)+(1/d4^2)+(1/d5^2)+(1/d6^2)+(1/d7^2))
  return(sigma_wz/sigma_w)
}

inverseDitanceWeight_1_random_random = function(interpolated_point_normalized,random_point, dataset_normalized, index_interpolated, row_interpolated, id_polutant){
  d1 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 2], 5], dataset_normalized[index_interpolated[row_interpolated, 2], 6])
  d2 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 3], 5], dataset_normalized[index_interpolated[row_interpolated, 3], 6])
  d3 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 4], 5], dataset_normalized[index_interpolated[row_interpolated, 4], 6])
  d4 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 5], 5], dataset_normalized[index_interpolated[row_interpolated, 5], 6])
  d5 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], dataset_normalized[index_interpolated[row_interpolated, 6], 5], dataset_normalized[index_interpolated[row_interpolated, 6], 6])
  #random-random
  d6 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], random_point[row_interpolated, 3], random_point[row_interpolated, 4])
  d7 = euclideanDistance(random_point[row_interpolated, 1], random_point[row_interpolated, 2], random_point[row_interpolated, 5], random_point[row_interpolated, 6])
 
  a1 = dataset_normalized[index_interpolated[row_interpolated, 2], id_polutant]
  a2 = dataset_normalized[index_interpolated[row_interpolated, 3], id_polutant]
  a3 = dataset_normalized[index_interpolated[row_interpolated, 4], id_polutant]
  a4 = dataset_normalized[index_interpolated[row_interpolated, 5], id_polutant]
  a5 = dataset_normalized[index_interpolated[row_interpolated, 6], id_polutant]
  a6 = interpolated_point_normalized[row_interpolated, 7]
  a7 = interpolated_point_normalized[row_interpolated, 8]
  sigma_wz = as.vector((a1/d1^2)+(a2/d2^2)+(a3/d3^2)+(a4/d4^2)+(a5/d5^2)+(a6/d6^2)+(a7/d7^2))
  sigma_w = as.vector((1/d1^2)+(1/d2^2)+(1/d3^2)+(1/d4^2)+(1/d5^2)+(1/d6^2)+(1/d7^2))
  return(sigma_wz/sigma_w)
}

interpolated_point_normalized_co2 = as.data.frame(matrix(nrow = nrow(index_interpolated), ncol = ncol(index_interpolated)))

for(i in 1:nrow(interpolated_point_normalized_co2)){
  for(j in 1:ncol(interpolated_point_normalized_co2)){
    if(index_interpolated[i, j] != 0)
      interpolated_point_normalized_co2[i, j] = dataset_normalized[index_interpolated[i, j], 7]
    else if (index_interpolated[i, j] == 0 & j == 7) 
      interpolated_point_normalized_co2[i, j] = inverseDitanceWeight_3(random_point, dataset_normalized, index_interpolated, i, 7)
    #j=8
    else if(index_interpolated[i, j] == 0 & j == 8 & index_interpolated[i, 7] != 0)
      interpolated_point_normalized_co2[i, j] = inverseDitanceWeight_2_fix(interpolated_point_normalized_co2,random_point, dataset_normalized, index_interpolated, i, 7)
    else if(index_interpolated[i, j] == 0 & j == 8 & index_interpolated[i, 7] == 0)
      interpolated_point_normalized_co2[i, j] = inverseDitanceWeight_2_random(interpolated_point_normalized_co2,random_point, dataset_normalized, index_interpolated, i, 7)
    #j=9
    else if(index_interpolated[i, j] == 0 & j == 9 & index_interpolated[i, 7] != 0 & index_interpolated[i, 8] != 0 )
        interpolated_point_normalized_co2[i, j] = inverseDitanceWeight_1_fix_fix(interpolated_point_normalized_co2, random_point, dataset_normalized, index_interpolated, i, 7)
    else if(index_interpolated[i, j] == 0 & j == 9 & index_interpolated[i, 7] != 0 & index_interpolated[i, 8] == 0)
        interpolated_point_normalized_co2[i, j] = inverseDitanceWeight_1_fix_random(interpolated_point_normalized_co2, random_point, dataset_normalized, index_interpolated, i, 7)
    else if(index_interpolated[i, j] == 0 & j == 9 & index_interpolated[i, 7] == 0 & index_interpolated[i, 8] == 0)
        interpolated_point_normalized_co2[i, j] = inverseDitanceWeight_1_random_random(interpolated_point_normalized_co2, random_point, dataset_normalized, index_interpolated, i, 7)
  }
}

interpolated_point_normalized_co = as.data.frame(matrix(nrow = nrow(index_interpolated), ncol = ncol(index_interpolated)))

for(i in 1:nrow(interpolated_point_normalized_co)){
  for(j in 1:ncol(interpolated_point_normalized_co)){
    if(index_interpolated[i, j] != 0)
      interpolated_point_normalized_co[i, j] = dataset_normalized[index_interpolated[i, j], 8]
    else if (index_interpolated[i, j] == 0 & j == 7) 
      interpolated_point_normalized_co[i, j] = inverseDitanceWeight_3(random_point, dataset_normalized, index_interpolated, i, 8)
    #j=8
    else if(index_interpolated[i, j] == 0 & j == 8 & index_interpolated[i, 7] != 0)
      interpolated_point_normalized_co[i, j] = inverseDitanceWeight_2_fix(interpolated_point_normalized_co,random_point, dataset_normalized, index_interpolated, i, 8)
    else if(index_interpolated[i, j] == 0 & j == 8 & index_interpolated[i, 7] == 0)
      interpolated_point_normalized_co[i, j] = inverseDitanceWeight_2_random(interpolated_point_normalized_co,random_point, dataset_normalized, index_interpolated, i, 8)
    #j=9
    else if(index_interpolated[i, j] == 0 & j == 9 & index_interpolated[i, 7] != 0 & index_interpolated[i, 8] != 0 )
      interpolated_point_normalized_co[i, j] = inverseDitanceWeight_1_fix_fix(interpolated_point_normalized_co, random_point, dataset_normalized, index_interpolated, i, 8)
    else if(index_interpolated[i, j] == 0 & j == 9 & index_interpolated[i, 7] != 0 & index_interpolated[i, 8] == 0)
      interpolated_point_normalized_co[i, j] = inverseDitanceWeight_1_fix_random(interpolated_point_normalized_co, random_point, dataset_normalized, index_interpolated, i, 8)
    else if(index_interpolated[i, j] == 0 & j == 9 & index_interpolated[i, 7] == 0 & index_interpolated[i, 8] == 0)
      interpolated_point_normalized_co[i, j] = inverseDitanceWeight_1_random_random(interpolated_point_normalized_co, random_point, dataset_normalized, index_interpolated, i, 8)
  }
}

#####     Spatial Model     #####
#####     CO2       #####
#making blank dataframe to store concentration value from nearest point(by index)
datafull_co2 = as.data.frame(matrix(nrow = nrow(nearest_index_fix), ncol = 9))

#storing nearest concentration value to datasetfull_co2
k = 1
for (i in 1:9){
  for (j in 1: nrow(nearest_index_fix)){
    if(nearest_index_fix[j, i] == 0){
      datafull_co2[j, i] = interpolated_point_normalized_co2[k, i]
      k = k + 1
    }
    else
      datafull_co2[j, i] = dataset_normalized[, 7][nearest_index_fix[j, i]]
  }
  k = 1
}

id = nearest_index_fix[, 1]

datafull_co2 = cbind(id, datafull_co2)

#data sampling (0.8 train, 0.2 test)
set.seed(123)
split_co2 = sample.split(datafull_co2$V1, SplitRatio = 0.8)
training_set_co2 = subset(datafull_co2, split_co2 == TRUE)
test_set_co2 = subset(datafull_co2, split_co2 == FALSE)

# using SVR
regressor_co2 = svm(x = training_set_co2[, 3:10] ,
                y = training_set_co2[, 2],
                type = 'eps-regression',
                kernel = 'radial',
                gamma = 0.125)

#storing predicted value in variable predictedY
predicted_co2 = predict(regressor_co2, test_set_co2[, 3:10])

#denormalize data
test_denormalized_co2 = denormalizeData(test_set_co2$V1, getNormParameters(dataset_normalized$CO201500))
predicted_denormalized_co2 = denormalizeData(predicted_co2, getNormParameters(dataset_normalized$CO201500))

# write.csv(test_normalized_co2, file = "test.csv")
# write.csv(predictedY1_normalized_co2, file = "predict.csv")

#Root Mean Squared Error of CO2 Spatial Prediction
predictionRMSE_co2 = rmse(test_denormalized_co2, predicted_denormalized_co2)

#correlation
correlation_co2 = cor(x = test_set_co2[, 2],
                  y = predict(regressor_co2, newdata = test_set_co2[, 3:10]), 
                  method = c("pearson", "kendall", "spearman"))

#plotting the actual and predcited data
ggplot() +
  geom_line(aes(x = test_set_co2$id , y = test_denormalized_co2),
             colour = 'red') +
  geom_line(aes(x = test_set_co2$id , y = predicted_denormalized_co2),
             colour = 'blue') 

# # perform a grid search
# tuneResult2 = tune(svm, train.x = training_set_co2[, 3:10] ,
#                   train.y = training_set_co2[, 2],
#                   ranges = list(epsilon = seq(0.01,0.0001), cost = 1:10)
# )
# print(tuneResult2)
# # Draw the tuning graph
# plot(tuneResult2)

####
#https://www.kdnuggets.com/2017/03/building-regression-models-support-vector-regression.html

##### CO #####
#making blank dataframe to store concentration value from nearest point(by index)
datafull_co = as.data.frame(matrix(nrow = nrow(nearest_index_fix), ncol = 9))

#storing nearest concentration value to dataset2
#interpolation with mean if index == 0
k = 1
for (i in 1:9){
  for (j in 1: nrow(nearest_index_fix)){
    if(nearest_index_fix[j, i] == 0){
      datafull_co[j, i] = interpolated_point_normalized_co[k, i]
      k = k + 1
    }
    else
      datafull_co[j, i] = dataset_normalized[, 8][nearest_index_fix[j, i]]
  }
  k = 1
}
id = nearest_index_fix[, 1]

datafull_co = cbind(id, datafull_co)

#data sampling (0.8 train, 0.2 test)
set.seed(123)
split_co = sample.split(datafull_co$V1, SplitRatio = 0.8)
training_set_co = subset(datafull_co, split_co == TRUE)
test_set_co = subset(datafull_co, split_co == FALSE)

# using SVR
regressor_co = svm(x = training_set_co[ ,3:10],
                y = training_set_co[ ,2],
                type = 'eps-regression',
                kernel = 'radial',
                gamma = 0.125 )

#storing predicted value in variable predictedY
predicted_co = predict(regressor_co, test_set_co[, 3:10])


#denormalize data
test_co_denormalized = denormalizeData(test_set_co$V1, getNormParameters(dataset_normalized$CO01500))
predicted_co_denormalized = denormalizeData(predicted_co, getNormParameters(dataset_normalized$CO01500))

#Root Mean Squared Error of CO2 Spatial Prediction
predictionRMSE_co = rmse(test_co_denormalized, predicted_co_denormalized)

#correlation
correlation_co = cor(x = test_set_co[, 2],
                  y = predict(regressor_co, newdata = test_set_co[, 3:10]), 
                  method = c("pearson", "kendall", "spearman"))

#plotting the actual and predcited data
ggplot() +
  geom_line(aes(x = test_set_co$id , y = test_co_denormalized),
            colour = 'red') +
  geom_line(aes(x = test_set_co$id , y = predicted_co_denormalized),
            colour = 'blue') 

# # perform a grid search
# tuneResult2 = tune(svm, train.x = training_set_co2[, 3:10] ,
#                    train.y = training_set_co2[, 2],
#                    ranges = list(epsilon = seq(0.01,0.0001), cost = 1:10)
# )
# print(tuneResult2)
# # Draw the tuning graph
# plot(tuneResult2)













