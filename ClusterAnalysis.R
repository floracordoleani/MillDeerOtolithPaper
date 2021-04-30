
#---------------------#
# Clustering analysis #
#---------------------#

### Code written by Anthony Malkassian and Flora Cordoleani 
### Use of clustering analysis to quantify Mill/Deer Creek juvenile life history strategies
# Load libraries needed to run this script
library(fda)
library(mclust)
library(ggplot2)
library(wesanderson)
library(ggpubr)


# Strontium profile smoothing with splines ----------------------------------------------

### Creating a spline function that will work within dplyr environments

spline.fun = function(x, xvar, yvar, ...) {
  smooth.spline(x = x[,xvar], y = x[,yvar], ...)
}

df <- 10 ##This value controls the "wigglyness" of the spline. Higher values of df (degrees of freedom) will make a wigglier curve that more closely matches the data. Lower values are a more smoothed fit. 

oto_data <- read.csv('Data/MillDeerOtoliths.csv',stringsAsFactors = FALSE)
oto_data <- oto_data[complete.cases(oto_data),]
oto_data_list <- oto_data %>% split(f = oto_data$sample) 
oto_data_list <- unique(oto_data_list)

### Fitting a spline for each Sr profile and storing the predicted Sr values 
x_vec <- seq(0,1000)
predicted_Sr <- data.frame(matrix(NA, nrow = length(unique(oto_data$sample)), ncol = length(x_vec)))
vwshed=vyear=vname=tab=tab2=vector()

for (i in 1:length(unique(oto_data$sample))){
  X <- data.frame(oto_data_list[i])
  oto_spl <- spline.fun(X, xvar = "distance", yvar = "otoSr", df = df)
  
  ## create a plateau for the ocean part of the profile
  # find the last Sr value of the profile
  maxSr <- X$otoSr[dim(X)[1]]
  
  # Fit a spline with a fixed plateau for the ocean phase. This helps for the clustering
  # if the total profile distance is < 1001, fix the Sr value for x_vec = last data point of the profile - 1001 with the max Sr value obtained from the spline
  if (X$distance[which(X$otoSr==maxSr)] < 1001){
    #predicted_Sr[i,1:(round(x_dist[length(x_dist)])+1)] <- predict(oto_spl,x_temp)$y
    predicted_Sr[i,] <- predict(oto_spl,x_vec)$y
    x_dist <- X$distance[which(X$otoSr==maxSr)]
    x_temp <- seq(0,round(x_dist[length(x_dist)])) # if there is more than one distance that has the max Sr value take the last distance with this value
    predicted_Sr[i,(round(x_dist[length(x_dist)])+2):1001] <- max(predict(oto_spl,x_temp)$y)
    
    # else if the total distance is > 1001 fit the spline to the entire profile (i.e x_vec = 0-1000)
  } else {
    predicted_Sr[i,] <- predict(oto_spl,x_vec)$y
  }
  
  vyear=c(vyear,X$year[1])
  vwshed=c(vwshed,X$watershed[1])
  vname=c(vname,X$sample[1])
  tab=cbind(tab,X$otoSr)
  tab2=cbind(tab2,X$distance)
  
}

# Create functionnal object using bspline basis system ---------------------------------------------
## (cf doc https://cran.r-project.org/web/packages/fda/fda.pdf p70)

mabasecut=create.bspline.basis(c(0,400),nbasis=50)
mabase=create.bspline.basis(c(0,1000),nbasis=50)

tfd=Data2fd(x_vec,t(predicted_Sr),basisobj=mabase)
tfdcut=Data2fd(x_vec[200:400],t(predicted_Sr)[200:400,],basisobj=mabasecut)

# Clustering analysis; step 1: focus on the 200 - 400 Sr profile zone ----------------------------------

## Apply functionnal pca to tfdcut object 
## (cf pca.fd() doc https://cran.r-project.org/web/packages/fda/fda.pdf p181)

tfpca=pca.fd(tfdcut,10)
# set of eigenvalues
tfpca$values
#proportion of variance explained by each eigenfunction
tfpca$varprop

#we keep only 3 harmonics regarding the proportion of variance explained by the first two eigenfunction
fpca=pca.fd(tfdcut,3)

#plot functionnal PCA
plot(fpca)

#plot functionnal pca harmonics
plot(fpca$harmonics)

## Then apply model based clustering with differents models (VVV, VEV, etc...)
##cf doc https://cran.r-project.org/web/packages/mclust/mclust.pdf p93)
clus1 <- Mclust(fpca$scores)

#summary function shows optimal clustering with VVE model with 4 components
summary(clus1,parameters = TRUE)

## best model given : Mclust EEE (ellipsoidal, equal volume, shape and orientation) model with 2 components:
#BIC or Bayesian Information Criterion for the specified mixture models and numbers of clusters
plot(clus1, what = "BIC")  # optimal clustering found for model VVV with 3 clusters

plot(clus1, what = "classification")

plot(clus1, what = "uncertainty")


# Clustering analysis; step 2: Separation based on the full profiles ------------------------------------
tfd2=tfd[clus1$classification==2,]

## Apply functionnal pca to tfd object 
##(cf doc https://cran.r-project.org/web/packages/fda/fda.pdf p181)

tfpca=pca.fd(tfd2,10)

# set of eigenvalues
tfpca$values

#proportion of variance explained by each eigenfunction
tfpca$varprop

#we keep only 1 harmonic regarding the proportion of variance explained by the first two eigenfunction
fpca=pca.fd(tfd2,1)

#plot functionnal PCA
plot(fpca)

#plot functionnal pca harmonics
plot(fpca$harmonics)

## Apply model based clustering on the remaining two population subset
clus <- Mclust(fpca$scores)

#summary function shows optimal clustering with EEV (ellipsoidal, equal volume and shape) model with 3 components:
summary(clus,parameters = TRUE)

#BIC or Bayesian Information Criterion for the specified mixture models and numbers of clusters
plot(clus, what = "BIC")  

plot(clus, what = "classification")

plot(clus, what = "uncertainty")

# Cluster Plots -------------------------------------------------------------------
### Assign each fish to associated cluster
FishID <- unique(oto_data$sample)
EarlyOutmigrant_ID <- FishID[which(clus1$classification==1)]
Other_ID <- FishID[which(clus1$classification==2)]
IntermediateOutmigrant_ID <- Other_ID[which(clus$classification==2)]
LateOutmigrant_ID <- Other_ID[which(clus$classification==1)]

cluster1 <- data.frame(sample = EarlyOutmigrant_ID) %>% 
            mutate(reartype ='EarlyOutmigrant')

cluster2 <- data.frame(sample = IntermediateOutmigrant_ID) %>% 
  mutate(reartype ='IntermediateOutmigrant')

cluster3 <- data.frame(sample =LateOutmigrant_ID) %>% 
  mutate(reartype ='LateOutmigrant')

subcluster <- data.frame(sample = Other_ID) %>% 
  mutate(reartype ='OtherOutmigrant')

Allclusters <- data.frame(rbind(cluster1,cluster2,cluster3))
clusters23 <- data.frame(rbind(cluster2,cluster3))
clusters1sub <- data.frame(rbind(cluster1,subcluster))

oto_clust1sub_merged <- merge(oto_data,clusters1sub,by="sample")
oto_clust23_merged <- merge(oto_data,clusters23,by="sample")
oto_Allclust_merged <- merge(oto_data,Allclusters,by="sample")

# Save clustering analysis results for use in MasterCode.R
write.csv(Allclusters, file = 'Data/MillDeerClusters.csv',row.names = FALSE)

C1 <- ggplot()+ 
  geom_line(data=oto_clust1sub_merged,aes(distance,otoSr,group=sample,colour=reartype))+ 
  labs(y=expression(paste({}^"87","Sr/",{}^"86","Sr")))+ xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20),legend.title = element_blank(),
        legend.position = "none")+
  scale_color_manual(values = c( "red","dodgerblue2"))

C1

C2 <- ggplot()+ 
  geom_line(data=oto_clust23_merged,aes(distance,otoSr,group=sample,colour=reartype))+ 
  labs(y=expression(paste({}^"87","Sr/",{}^"86","Sr")))+ xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20),legend.title = element_blank(),
        legend.position = "none")+
  scale_color_manual(values = c( "darkcyan","darkgoldenrod2"))


C2

C3 <- ggplot()+ 
  geom_line(data=oto_Allclust_merged,aes(distance,otoSr,group=sample,colour=reartype))+ 
  labs(x=(expression(paste('Otolith radius (',mu,'m)',sep = ''))))+
  labs(y=expression(paste({}^"87","Sr/",{}^"86","Sr")))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20),legend.title = element_blank(),
        legend.position = "none")+
  scale_color_manual(values = wes_palette("Darjeeling1", n = 3))

C3

### Supplemental Material Figure S1

S1 <- ggarrange(C1,C2,C3,
                labels = c("a","b","c"),
                font.label=list(size = 25,color="black"),
                nrow=3)

png("Figures/FigureS1.png", 
    family = "sans serif", width = 5, height=9, units = "in", res =300)

S1

dev.off()
