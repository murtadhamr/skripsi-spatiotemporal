library(tseries)

data_hotspot2 <- read.csv(file.choose(),header=TRUE)
## plot data titik panas bulanan ##
#data_hotspot2=data 2001-2013
#data_hotspot(didalamfolder data asli)=data 2001-2015
plot(data_hotspot2$hotspot, type = "l")

#normalisasi
input<- normalizeData(data_hotspot2[,3], type='0_1')

## Uji Augmented Dickey Fuller ##
adf.test(input)

## Uji Bartlett and Levene ##
##bartlett ini pake data yang dinormalisasi, atributnya bulan dan titik panas
bartlett.test(hotspot_norm)

data_hotspot2$hotspot<-data_hotspot2$hotspot+0.000001
## Mencari nilai koefisien ?? transformasi Boxcox ##
a<- boxcox(hotspot ~ month, data=data_hotspot2, plotit=T, 
           lambda = seq(0.08, 0.2, length = 10))

##mencari nilai terkecil pada plot boxcox untuk lambda
which.max(a$y)
a$x[60]

## Melakukan Transformasi Boxcox ##
hotspot_transform<-input^0.1515152

## Melakukan plot data hasil Transformasi Boxcox ##
plot(hotspot_transform)

## Melakukan plot ACF ##
acf(hotspot_transform)

## Melakukan plot PACF ##
pacf(hotspot_transform)

#untuk denormalisasi
a<-denormalizeData(values,getNormParameters(values))
