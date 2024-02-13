library("readxl")
colon_collagens<-read_excel("Path to colon_collagen_PCA.xlsx file")
library(missMDA)
library(ggplot2)
library(factoextra)

#Visualization of PCA plot for entire collagen subtypes without any 'NAN' in CPTAC colon cancer database
colon_collagens_PCA = prcomp(colon_collagens[,c(-1,-2)])
fviz_pca_ind(colon_collagens_PCA, habillage=colon_collagens$type, addEllipses=TRUE)
summary(colon_collagens_PCA)


