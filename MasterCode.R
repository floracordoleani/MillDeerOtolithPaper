############ Code written by Flora Cordoleani to analyze Mill/Deer Creek spring-run otolith data 

#------------------------#
# Packages and functions #
#------------------------#

# Load libraries needed to run this script
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggridges)
library(wesanderson)


# Load function
move.size <- function(x,Threshold_Sr8786){ # Function used to look at natal and freshwater exit
  if(sum(is.na(x$otoSr) > 0)) x<-x[!is.na(x$otoSr),]
  pre.x <- x[x$otoSr < Threshold_Sr8786, ] 
  #  if no otoSr measurements are below threshold value return NA
  if (nrow(pre.x) == 0) return(NA)
  #  Maximum otolith distance of < thresh data subset  
  distPre <- max(pre.x$distance, na.rm = TRUE)
  #  Last row before thresh exceeded  
  pre.entry <- pre.x[pre.x$distance == distPre, ]  
  #  All rows greater than the max distance of the < thresh data subset
  post.x <- x[x$distance > distPre, ]
  #  if first distance value post threshold is NA, return NA
  if (is.na(post.x$distance[1])) return(NA)
  #  Minimum otolith distance of post threshold data subset
  distPost <- min(post.x$distance, na.rm = TRUE)
  #  First row after threshold is exceeded for the last time
  post.entry <- post.x[post.x$distance == distPost, ]
  
  SrPre <- pre.entry$otoSr
  SrPost <- post.entry$otoSr
  
  distMix <- distPre + (Threshold_Sr8786 - SrPre) * ((distPost - distPre)/(SrPost - SrPre))
  return(distMix)
}


#---------------------------------#
# Strontium (Sr) Profile Analysis #
#---------------------------------#

# Sr profile figure ----------------------------------------------------------------------

### Read data files
Otodata <- read.csv('Data/MillDeerOtoliths.csv',stringsAsFactors = FALSE)
Clusterdata <- read.csv('Data/MillDeerClusters.csv',stringsAsFactors = FALSE)

Otodata_complete <- merge(Otodata,Clusterdata,by="sample")

### Prepare profiles for figure
SubsetEarly <- data.frame(Otodata_complete[Otodata_complete$reartype=="EarlyOutmigrant",])
SubsetIntermediate <- data.frame(Otodata_complete[Otodata_complete$reartype=="IntermediateOutmigrant",] )
SubsetLate <- data.frame(Otodata_complete[Otodata_complete$reartype=="LateOutmigrant",] )

earlyexample <- data.frame(SubsetEarly[SubsetEarly$sample=='FC10-9-13-3',])
intermediateexample <- data.frame(SubsetIntermediate[SubsetIntermediate$sample=='FC10-4-13-001',]) 
lateexample <- SubsetLate[SubsetLate$sample== 'DC18_68517',] 


x = c(0,1000)
position_trib <- data.frame(x=x,y=c(0.7035, 0.7049647)) #range of Mill/Deer Cr Sr values
position_sac <- data.frame(x=x,y=c(0.7049647,0.7061)) #range of Sacremento River Sr values
position_delta <- data.frame(x=x,y=c(0.7061, 0.7078)) #range of Delta Sr values
position_bay<- data.frame(x=x,y=c(0.7078, 0.710)) #range of Bay Sr values
position_ocean <- data.frame(x=x,y=c(0.70918, 0.715)) #range of Ocean Sr values
position_incub <- data.frame(x=c(0,200),y=c(0.7035,0.710)) # show incubation period

profile1 <- ggplot() +
  labs(y=expression(paste({}^"87","Sr/",{}^"86","Sr")))+ xlab("")+
  ggtitle('Early migrant') + 
  scale_x_continuous(limits=c(0,1200), breaks=seq(0,1200,200))+
  scale_y_continuous(limits=c(0.7035,0.710), breaks=seq(0.704,0.710,0.001))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=25),legend.title = element_blank(),
        plot.title = element_text(face = "bold"))+
  geom_rect(data= position_trib, inherit.aes = FALSE,
            aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]), 
            fill="#001B87", alpha=0.8)+ 
  geom_rect(data= position_sac, inherit.aes = FALSE,
            aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]), 
            fill= "#0064D6",alpha=0.5)+ 
  geom_rect(data= position_delta, inherit.aes = FALSE,
            aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]), 
            fill="#00A6D7", alpha=0.2)+ 
  geom_rect(data= position_bay, inherit.aes = FALSE,
            aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]), 
            fill="#C1E6E0", alpha=1)+
  geom_rect(data= position_incub, inherit.aes = FALSE, 
            aes(xmin=-Inf, xmax=200,ymin=0.7035, ymax=0.710), 
            fill="lightgrey", alpha=0.8)+
  annotate(geom="text",x=1100, y=0.7097, label=" Bay & Ocean", col="black", size = 7,fontface="bold")+
  annotate(geom="text",x=1100, y=0.707,label="Delta", col="black", size = 7,fontface="bold")+
  annotate(geom="text",x=1100, y=0.7055,label="Sacramento 
  River", col="black", size = 7,fontface="bold")+
  annotate(geom="text",x=1100, y=0.7042,label="Mill & Deer 
  Creeks", col="grey60", size = 7,fontface="bold")+
  annotate(geom="text",x=70, y=0.7097,label="Incubation", 
           col="black", size = 7,fontface="bold")+
  geom_line(data=SubsetEarly ,aes(distance,otoSr,group=sample),col="grey50" ,lwd=1,alpha=0.6) + 
  geom_line(data=earlyexample,aes(distance,otoSr,group=sample),col="black")+
  geom_point(data=earlyexample,aes(distance,otoSr,group=sample),col="black")+
  geom_segment(data=earlyexample,aes(x = distance, y = otoSr - 1.96*SE1, 
                                        xend = distance, yend =  otoSr + 1.96*SE1),col="black")+
  geom_segment(aes(x=200,xend=+Inf,y=0.7049647,yend=0.7049647), 
               col="black", lty="dashed",lwd = 0.8)+
  geom_segment(aes(x=200,xend=+Inf,y=0.7061,yend=0.7061), 
               col="black", lty="dashed",lwd = 0.8)+
  geom_segment(aes(x=200,xend=+Inf,y= 0.7078,yend= 0.7078), 
               col="black", lty="dashed",lwd = 0.8)+
  geom_vline(xintercept = 200,col="black", lty="dashed",lwd = 0.8)

profile1

profile2 <- ggplot() +
  labs(y=expression(paste({}^"87","Sr/",{}^"86","Sr"))) + xlab("")+
  ggtitle('Intermediate migrant') + 
  scale_x_continuous(limits=c(0,1200), breaks=seq(0,1200,200))+ #
  scale_y_continuous(limits=c(0.7035,0.710), breaks=seq(0.704,0.710,0.001))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=25),legend.title = element_blank(),
        plot.title = element_text(face = "bold"))+
  geom_rect(data= position_trib, inherit.aes = FALSE,
            aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]), 
            fill="#001B87", alpha=0.8)+ # "#F2AD00" "#FCD16B"
  geom_rect(data= position_sac, inherit.aes = FALSE,
            aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]), 
            fill= "#0064D6",alpha=0.5)+ # "#D1362F" ",#c03728"
  geom_rect(data= position_delta, inherit.aes = FALSE,
            aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]), 
            fill="#00A6D7", alpha=0.2)+ # "#0066CC"   "#2166AC"
  geom_rect(data= position_bay, inherit.aes = FALSE,
            aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]), 
            fill="#C1E6E0", alpha=1)+# "#009966" "#3c5e27"
  geom_rect(data= position_incub, inherit.aes = FALSE, 
            aes(xmin=-Inf, xmax=200,ymin=0.7035, ymax=0.710), 
            fill="lightgrey", alpha=0.8)+
  annotate(geom="text",x=1100, y=0.7097, label=" Bay & Ocean", col="black", size = 7,fontface="bold")+
  annotate(geom="text",x=1100, y=0.707,label="Delta", col="black", size = 7,fontface="bold")+
  annotate(geom="text",x=1100, y=0.7055,label="Sacramento 
  River", col="black", size = 7,fontface="bold")+
  annotate(geom="text",x=1100, y=0.7042,label="Mill & Deer 
  Creeks", col="grey60", size = 7,fontface="bold")+
  annotate(geom="text",x=70, y=0.7097,label="Incubation", 
           col="black", size = 7,fontface="bold")+
  geom_line(data=SubsetIntermediate,aes(distance,otoSr,group=sample),col="grey50" ,lwd=1,alpha=0.6) + #"grey30"
  geom_line(data=intermediateexample ,aes(distance,otoSr,group=sample),col="black")+
  geom_point(data=intermediateexample,aes(distance,otoSr,group=sample),col="black")+
  geom_segment(data=intermediateexample,aes(x = distance, y = otoSr - 1.96*SE1, 
                                                xend = distance, yend =  otoSr + 1.96*SE1),col="black")+
  geom_segment(aes(x=200,xend=+Inf,y=0.7049647,yend=0.7049647), 
               col="black", lty="dashed",lwd = 0.8)+
  geom_segment(aes(x=200,xend=+Inf,y=0.7061,yend=0.7061), 
               col="black", lty="dashed",lwd = 0.8)+
  geom_segment(aes(x=200,xend=+Inf,y= 0.7078,yend= 0.7078), 
               col="black", lty="dashed",lwd = 0.8)+
  geom_vline(xintercept = 200,col="black", lty="dashed",lwd = 0.8)

profile2

profile3 <- ggplot() +
  labs(x=(expression(paste('Otolith radius (',mu,'m)',sep = ''))),
       y=(expression(paste({}^"87","Sr/",{}^"86","Sr"))))+
  ggtitle('Late migrant') + 
  scale_x_continuous(limits=c(0,1200),breaks=seq(0,1200,200))+
  scale_y_continuous(limits=c(0.7035,0.710), breaks=seq(0.704,0.710,0.001))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=25),legend.title = element_blank(),
        plot.title = element_text(face = "bold"))+
  geom_rect(data= position_trib, inherit.aes = FALSE,
            aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]), 
            fill="#001B87", alpha=0.8)+ 
  geom_rect(data= position_sac, inherit.aes = FALSE,
            aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]), 
            fill= "#0064D6",alpha=0.5)+ 
  geom_rect(data= position_delta, inherit.aes = FALSE,
            aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]), 
            fill="#00A6D7", alpha=0.2)+
  geom_rect(data= position_bay, inherit.aes = FALSE,
            aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]), 
            fill="#C1E6E0", alpha=1)+
  geom_rect(data= position_incub, inherit.aes = FALSE, 
            aes(xmin=-Inf, xmax=200,ymin=0.7035, ymax=0.710), 
            fill="lightgrey", alpha=0.8)+
  annotate(geom="text",x=1100, y=0.7097, label=" Bay & Ocean", col="black", size = 7,fontface="bold")+
  annotate(geom="text",x=1100, y=0.707,label="Delta", col="black", size = 7,fontface="bold")+
  annotate(geom="text",x=1100, y=0.7055,label="Sacramento 
  River", col="black", size = 7,fontface="bold")+
  annotate(geom="text",x=1100, y=0.7042,label="Mill & Deer 
  Creeks", col="grey60", size = 7,fontface="bold")+
  annotate(geom="text",x=70, y=0.7097,label="Incubation", 
           col="black", size = 7,fontface="bold")+
  geom_line(data=SubsetLate,aes(distance,otoSr,group=sample),col="grey50" ,lwd=1,alpha=0.6) + 
  geom_line(data=lateexample,aes(distance,otoSr,group=sample),col="black")+
  geom_point(data=lateexample,aes(distance,otoSr,group=sample),col="black")+
  geom_segment(data=lateexample,aes(x = distance, y = otoSr - 1.96*SE1, 
                                       xend = distance, yend =  otoSr + 1.96*SE1),col="black")+
  geom_segment(aes(x=200,xend=+Inf,y=0.7049647,yend=0.7049647), 
               col="black", lty="dashed",lwd = 0.8)+
  geom_segment(aes(x=200,xend=+Inf,y=0.7061,yend=0.7061), 
               col="black", lty="dashed",lwd = 0.8)+
  geom_segment(aes(x=200,xend=+Inf,y= 0.7078,yend= 0.7078), 
               col="black", lty="dashed",lwd = 0.8)+
  geom_vline(xintercept = 200,col="black", lty="dashed",lwd = 0.8)

profile3

### Paper figure 1

p1 <- ggarrange(nrow=3,profile1, profile2,profile3,
          labels = c("a","b","c"),
          font.label=list(size = 30,color="black"))

png("Figures/Figure1.png", 
    family = "sans serif", width =9, height= 18, units = "in", res =200)

p1

dev.off()

#----------------------------------------------------------------------#
# Otolith radius (OR) distance at Natal and Freshwater Exit Estimation #
#----------------------------------------------------------------------#

# OR distance at Natal Exit ---------------------------------------------------------------------

### Read data files
NatalExitdata <- read.csv('Data/MillDeerNatalExit.csv',stringsAsFactors = FALSE) # Natal exit OR distance determined visually 
Incdata <- read.csv('Data/MillDeerIncrements.csv',stringsAsFactors = FALSE)

###  Estimate OR distance at natal exit
SacRiverThreshold <- 0.7046338 # Strontium threshold value used for the Sacramento River

for (i in 1:dim(NatalExitdata)[1]){
  subset <- subset(Otodata,sample == NatalExitdata$sample[i]) 
  NatalExitdata$NatalExit_OR_estimated[i] <-  move.size(subset,SacRiverThreshold) 
  
  if (is.na(NatalExitdata$NatalExit_OR_estimated[i])){ #Sacramento River Sr threshold was never crossed
    NatalExitdata$NatalExit_OR_final[i] <- NatalExitdata$NatalExit_OR[i] # OR distance identified visually is used as final distance at natal exit
  } else {
    NatalExitdata$NatalExit_OR_final[i] <- NatalExitdata$NatalExit_OR_estimated[i] # OR distance estimated is used as final distance at natal exit
  }
}


### Add number of days spent in natal tributary 
for (i in 1:dim(NatalExitdata)[1]){
  dataSubset<-subset(Incdata,sample == NatalExitdata$sample[i]) 
  if (dim(dataSubset)[1]!=0){
    NatalExit_OR_index <- which.min(abs(dataSubset$inc_distance - NatalExitdata$NatalExit_OR_final[i]))
    NatalExitdata$NatalExit_IncNum[i] <- dataSubset$inc_num[NatalExit_OR_index]  
  } else NatalExitdata$NatalExit_IncNum[i] <- 'NA' # Microchemistry not performed
}

### Add rearing strategy cluster information
NatalExitdata_complete <- merge(NatalExitdata,Clusterdata,by="sample")

# Change object category for future analysis
NatalExitdata_complete$reartype <- as.factor(NatalExitdata_complete$reartype)
NatalExitdata_complete$NatalExit_IncNum <- as.numeric(NatalExitdata_complete$NatalExit_IncNum)

### Summary statistics for manuscript
# Evaluating number of fish otolith analyzed each year in each trib
oto.summary <- NatalExitdata_complete %>% group_by(watershed,year) %>% 
  dplyr::summarise(count=n())

# Evaluating number of fish with each rearing strategy each year
strategy.summary <- NatalExitdata_complete %>% group_by(watershed,year,reartype) %>% 
  dplyr::summarise(count=n())

# Evaluating the proportion of rearing strategy expression
propstrategy.summary <- strategy.summary %>% group_by(reartype) %>% 
  dplyr::summarise(percent=round((sum(count)/123)*100))

# Evaluating OR mean and sd at natal exit for each rearing strategy
NatalOR.summary <- NatalExitdata_complete %>% group_by(reartype) %>% 
                   dplyr::summarise(Mean =round(mean(NatalExit_OR_final, na.rm=TRUE)),
                                    Sd = round(sd(NatalExit_OR_final, na.rm=TRUE)))

# Evaluating mean and sd increment number (=number of day) at natal exit for each rearing strategy
NatalInc.summary <- NatalExitdata_complete %>% group_by(reartype) %>% 
                    dplyr::summarise(Mean = round(mean(NatalExit_IncNum, na.rm=TRUE)),
                                     Sd = round(sd(NatalExit_IncNum, na.rm=TRUE)))

# Anova test
res.aov.NatalOR <- aov(NatalExit_OR_final ~ reartype , data = NatalExitdata_complete)
summary(res.aov.NatalOR)

NatalIncdata <-  NatalExitdata_complete %>% filter(NatalExit_IncNum!='NA')
res.aov.Natalinc <- aov(NatalExit_IncNum~ reartype , data = NatalIncdata)
summary(res.aov.Natalinc)

# Tukey test
tukey.res.NatalOR <- TukeyHSD(res.aov.NatalOR)
plot(tukey.res.NatalOR)

tukey.res.Natalinc <- TukeyHSD(res.aov.Natalinc)
plot(tukey.res.Natalinc)


# OR distance at Freshwater Exit -------------------------------------------------------------
ChippsThreshold <-0.7078 # Sr threshold value for Chipps Island, which corresponds to location of freshwater exit

###  Estimate OR distance at Freshwater (FW) exit
FWdata <- Otodata[,c("sample","distance","otoSr")]

# Remove columns with  Sample_ID= NA
FWdata <- FWdata[!is.na(FWdata$sample),]
FWdata$sample <- as.character(FWdata$sample)

# Generate unique fish IDs
FishID <-unique(FWdata$sample)

# Create empty data frame of FW exit to fill
FWExitdata <- data.frame(matrix(ncol=2,nrow=length(FishID), 
                                dimnames=list(NULL, c("sample", "FWExit_OR_estimated"))))

for (i in 1:length(FishID)){
  dataSubset<-subset(FWdata,FWdata$sample==FishID[i])
  dataSubset<-dataSubset[order(dataSubset$distance),]
  FWExitdata[i,"sample"] <- dataSubset$sample[1]
  FWExitdata[i,"FWExit_OR_estimated"] <-move.size(dataSubset,ChippsThreshold)
}


### Add number of days spent in freshwater for each fish
for (i in 1:length(FishID)){
  dataSubset<-subset(Incdata,sample == FWExitdata$sample[i]) 
  if (dim(dataSubset)[1]!=0){
    FWExit_OR_index <- which.min(abs(dataSubset$inc_distance -(FWExitdata$FWExit_OR_estimated[i])))
    FWExitdata$FWExit_IncNum[i] <- dataSubset$inc_num[FWExit_OR_index]  
  } else FWExitdata$FWExit_IncNum[i] <- 'NA' # Microchemistry not performed
}

FWExitdata_complete <- merge(x=FWExitdata,y=Clusterdata,by="sample")

# Change object category for future analysis
FWExitdata_complete$reartype <- as.factor(FWExitdata_complete$reartype)
FWExitdata_complete$FWExit_IncNum <- as.numeric(FWExitdata_complete$FWExit_IncNum)

### FW Exit statistics for manuscript
# Evaluating OR mean and sd at FW exit for each rearing strategy
FWOR.summary <- FWExitdata_complete %>% group_by(reartype) %>% 
                dplyr::summarise(Mean = round(mean(FWExit_OR_estimated, na.rm=TRUE)),
                                 Sd = round(sd(FWExit_OR_estimated, na.rm=TRUE)))

FWIncNum.summary <- FWExitdata_complete %>% group_by(reartype) %>% 
                    dplyr::summarise(Mean = round(mean(FWExit_IncNum, na.rm=TRUE)),
                                     Sd = round(sd(FWExit_IncNum, na.rm=TRUE)))

# Anova test
res.aov.FWOR <- aov(FWExit_OR_estimated ~ reartype, data = FWExitdata_complete)
summary(res.aov.FWOR)

FWIncdata <-  FWExitdata_complete %>% filter(FWExit_IncNum!='NA')
res.aov.FWinc <- aov(FWExit_IncNum ~ reartype, data = FWIncdata)
summary(res.aov.FWinc)

# Tukey test
tukey.res.FWOR <- TukeyHSD(res.aov.FWOR)
plot(tukey.res.FWOR)

tukey.res.FWinc <- TukeyHSD(res.aov.FWinc)
plot(tukey.res.FWinc)

# Natal & FW Exit density figures --------------------------------------------------------------

### Paper Figure 2
p2a <- ggplot(NatalExitdata_complete)+
  geom_density(aes(x=NatalExit_OR_final,y=..count../sum(..count..),
                   fill = reartype,color=reartype),alpha=0.3,size=1) +
  scale_y_continuous(limits=c(0,0.008), breaks=seq(0,0.008,0.002))+
  scale_x_continuous(limits=c(100,1000), breaks=seq(100,1000,200))+
  labs(x=(expression(paste('Otolith radius (',mu,'m)',sep = ''))))+
  ggtitle('Natal Exit') + 
  ylab("Density") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20),legend.title = element_blank(),
        legend.position = "none", plot.title = element_text(face = "bold"))+
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 3))+  
  scale_color_manual(values = wes_palette("Darjeeling1", n = 3))

p2a

p2b <- ggplot(FWExitdata_complete)+
  geom_density(aes(x=FWExit_OR_estimated, y=..count../sum(..count..),color=reartype, fill = reartype),
               alpha=0.3,adjust=1,size=1) +
  scale_y_continuous(limits=c(0,0.005), breaks=seq(0,0.005,0.001))+
  scale_x_continuous(limits=c(400,1000), breaks=seq(400,1000,100))+
  ylab("")+ 
  labs(x=(expression(paste('Otolith radius (',mu,'m)',sep = ''))))+
  ggtitle('Freshwater Exit') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20),legend.title = element_blank(),
        legend.position = "none", plot.title = element_text(face = "bold"))+
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 3))+  
  scale_color_manual(values = wes_palette("Darjeeling1", n = 3))

p2b

p2 <- ggarrange(p2a, p2b, 
          labels = c("a", "b"),
           font.label=list(size = 25,color="black"))

png("Figures/Figure2.png", 
    family = "sans serif", width = 9, height= 4, units = "in", res =200)

p2

dev.off()

### Paper figure 3
# Combine natal and freshwater exit dataframes
NatalFW_data <- merge(NatalExitdata_complete,FWExitdata_complete,by=c("sample","reartype"))

NatalFWdata_ridge <- melt(NatalFW_data, id.vars = "year",
                 measure.vars = c("NatalExit_OR_final", "FWExit_OR_estimated"))

p3a <- ggplot(NatalFWdata_ridge,aes(x = value,y=forcats::fct_rev(as.factor(year)),color = variable, 
                                 point_color = variable, fill = variable))+ 
  geom_density_ridges(alpha=0.3,adjust=0.8,size=1,scale = .95)+ 
  scale_x_continuous(limits=c(50,1000), breaks=seq(100,1000,200),expand = c(0,0))+
  scale_y_discrete(expand = c(0,0)) +
  labs(x=(expression(paste('Otolith radius (',mu,'m)',sep = ''))))+ ylab("")+ 
  scale_fill_manual(values = c("#001B87", "#00A6D7"), labels = c("Natal", "Freshwater"),name="") +
  scale_color_manual(values = c("#001B87", "#00A6D7"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#001B87", "#00A6D7"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#001B87", "#00A6D7"),
         color = NA, point_color = NA))) +
  theme_ridges(center = TRUE)+
  theme(text = element_text(size=20),
        axis.text.x = element_text(size=16),
        axis.text.y =element_text(size=16))

p3a 


NatalFWdata_cast = dcast(NatalFW_data, year ~ reartype, fun.aggregate = length)
NatalFWdata_melt = melt(NatalFWdata_cast, id.vars = "year",
                measure.vars = c("EarlyOutmigrant", "IntermediateOutmigrant","LateOutmigrant"))

p3b <- ggplot(NatalFWdata_melt, aes(x = as.factor(year),y = value, fill = variable)) + 
              geom_bar(stat = "identity" ,alpha=0.6,width = 12,position="stack")+
       labs(y = "Number of recovered otoliths", x="") +
       facet_grid(year~.,switch="y") +
       coord_flip()+
       theme(strip.background =element_rect(fill="white"),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank(),
             text = element_text(size=20),
             legend.title = element_blank())+
       scale_fill_manual(values = wes_palette("Darjeeling1",n=3),
                         labels = c("Early", "Intermediate", "Late"))

p3b


p3 <- cowplot::plot_grid(p3a+theme(legend.position="top"),
                p3b+ theme(legend.position="top"), 
                labels = c('a', 'b'), label_size = 20,
                nrow=1,
                rel_widths = c(0.6,1))

png("Figures/Figure3.png", 
    family = "sans serif", width = 10, height= 6, units = "in", res =200)

p3

dev.off()


#---------------------------#
# Increment Growth Analysis #
#---------------------------#

# Growth analysis on the first 15 and 30 rearing days after emergence ----------------------------------------------------------

### Assign rearing locations to increment data 
GrowthData <- merge(Incdata ,NatalExitdata_complete,by=c('sample','year','watershed'))
GrowthData$inc_num <- as.numeric(GrowthData$inc_num)
GrowthData$NatalExit_IncNum <- as.numeric(GrowthData$NatalExit_IncNum)

GrowthData_final <- GrowthData %>% 
  group_by(sample) %>%
  arrange(-desc(inc_num), .by_group = TRUE) %>%
  mutate(Habitat = ifelse(inc_num <= NatalExit_IncNum, 'Trib', 'Non-Trib')) %>%
  ungroup()

### Estimate growth during the first 15 or 30 days after emergence
growth15 <- GrowthData_final %>% filter(inc_num < 16)
growth15$year <- as.factor(growth15$year)

growth15_avg <- growth15 %>% group_by(sample) %>%
  summarize(mean_growth = mean(inc_width),
            cumgrowth = max(cuminc_width),
            dist_exit = mean(NatalExit_OR_final),
            inc_exit=mean(NatalExit_IncNum),
            Strat=unique(reartype))

growth30 <- GrowthData_final %>% filter(inc_num < 31)
growth30$year <- as.factor(growth30$year)

growth30_avg <- growth30 %>% group_by(sample) %>%
  summarize(mean_growth = mean(inc_width),
            cumgrowth = max(cuminc_width),
            dist_exit = mean(NatalExit_OR_final),
            inc_exit=mean(NatalExit_IncNum),
            Strat=unique(reartype))

growthday15 <- GrowthData_final  %>%  filter(inc_num == 15) %>% 
               filter(NatalExit_IncNum >= 15)
 
growthday30 <- GrowthData_final  %>%  filter(inc_num == 30) %>% 
  filter(NatalExit_IncNum >= 30)

### Growth statistics for manuscript 
growthday15.summary <- growthday15 %>% group_by(reartype) %>% 
                      dplyr::summarise(Mean =round(mean(cuminc_width, na.rm=TRUE)),
                                       Sd= round(sd(cuminc_width, na.rm=TRUE)))

growthday30.summary <- growthday30 %>% group_by(reartype) %>% 
  dplyr::summarise(Mean =round(mean(cuminc_width, na.rm=TRUE)),
                   Sd= round(sd(cuminc_width, na.rm=TRUE)))


res.aov_15 <- aov(mean_growth ~ Strat, data=growth15_avg[growth15_avg$inc_exit>=15,])
summary(res.aov_15)

tuk15 <- TukeyHSD(res.aov_15)
tuk15

plot(tuk15)


res.aov_30 <- aov(mean_growth ~ Strat, data=growth30_avg[growth30_avg$inc_exit>=30,])
summary(res.aov_30)

tuk30 <- TukeyHSD(res.aov_30)
tuk30
plot(tuk30)


### Paper figure 4
p4a <- ggplot(data=growth15_avg[growth15_avg$inc_exit>=15,],aes(x=inc_exit,y=mean_growth)) + 
  geom_point(aes(x=inc_exit,y=mean_growth,colour=Strat),size=4)+
  geom_smooth(color = 'black',method="lm", color="black") +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = 3,size=7,col="black")+
  geom_point(data=growth15_avg[growth15_avg$inc_exit<15,],aes(x=inc_exit,y=mean_growth),
             col="red",shape=1,size=4,stroke=2)+
  scale_x_continuous(limits=c(0,260),breaks = seq(0,250,by =50))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=45, hjust = 1),
        text = element_text(size=30),legend.title = element_blank())+
  ylab('Average increment width (µm)')+ xlab("Increment number at natal exit") +
  scale_fill_manual(values = wes_palette("Darjeeling1",n=3))+  
  scale_color_manual(values = wes_palette("Darjeeling1",n=3),
                     labels = c("Early", "Intermediate", "Late"))

p4a


value_med <- growthday15 %>% group_by(reartype) %>% summarize(med_value = median(cuminc_width))
res.aov_cumgrowth15 <- aov(cuminc_width ~ reartype, data=growthday15)
tuk_cumgrowth15 <- agricolae::HSD.test(res.aov_cumgrowth15, trt='reartype' , group = T)
sig.letters <- tuk_cumgrowth15$groups[order(row.names(tuk_cumgrowth15$groups)), ]

p4b <- ggplot(growthday15,aes(x = reartype, y = cuminc_width,fill=reartype)) +
  geom_boxplot(outlier.shape = NA,alpha=0.5,col="black")+
  geom_jitter(width=0.15,size=3,col="black")+
  geom_text(data = value_med, aes(x=reartype, y = 23 + med_value, 
                                  label = sig.letters$groups,colour=reartype),
                                  size=8, vjust=0,show.legend=FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        text = element_text(size=30),legend.title = element_blank(),
        legend.position = "none")+
  ylab('Cumulative width of first 15 increments')+ xlab("") + 
  scale_x_discrete(labels=c('Early','Intermediate','Late'))+
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 3),
                    labels = c("Early", "Intermediate", "Late"))+  
  scale_color_manual(values = wes_palette("Darjeeling1", n = 3))

p4b


p4 <- ggarrange(p4a,p4b, labels = c("a", "b"),
          nrow=2,
          font.label=list(size = 35,color="black"),
          common.legend =FALSE, legend = "right")


png("Figures/Figure4.png", 
    family = "sans serif", width = 13, height= 18, units = "in", res =300)

p4

dev.off()

### Supp Material figure S3

pS3a <- ggplot(data=growth30_avg[-which(growth30_avg$Strat=='EarlyOutmigrant'),],
             aes(x=inc_exit,y=mean_growth)) + 
  geom_point(aes(x=inc_exit,y=mean_growth,colour=Strat),size=4)+ 
  geom_smooth(color = 'black',method="lm", color="black") +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,label.x = 50,size=7,col="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=45, hjust = 1),
        text = element_text(size=30),legend.title = element_blank())+
  scale_x_continuous(limits=c(50,260),breaks = seq(50,250,by =50))+
  ylab('Average increment width (µm)')+ xlab("Increment number at natal exit") +
  scale_fill_manual(values = wes_palette("Darjeeling1")[2:3])+  
  scale_color_manual(values = wes_palette("Darjeeling1")[2:3],
                     labels = c("Intermediate", "Late"))

pS3a

pS3b <- ggplot(growthday30[-which(growthday30$reartype=='EarlyOutmigrant'),],
              aes(x = reartype, y = cuminc_width,fill=reartype)) +
  geom_boxplot(outlier.shape = NA,alpha=0.5,col="black")+
  geom_jitter(width=0.15,size=3,col="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        text = element_text(size=30),legend.title = element_blank())+
  ylab('Cumulative width of first 30 increments')+ xlab("") +
  scale_fill_manual(values = wes_palette("Darjeeling1")[2:3],
                    labels = c("Intermediate", "Late"))+  
  scale_color_manual(values = wes_palette("Darjeeling1")[2:3])

pS3b


pS3 <- ggarrange(pS3a, pS3b, labels = c("a", "b"),
          font.label=list(size = 35,color="black"),
          common.legend = FALSE, legend = "right",
          nrow=2)

png("Figures/FigureS3.png", 
    family = "sans serif", width = 13, height= 18, units = "in", res =300)

pS3

dev.off()


#-----------------------------------------#
# Rotary Screw Trap (RST) data  from CDFW #
#-----------------------------------------#

# Identify rearing strategies in RST data ---------------------------------------------------

RSTdata <- read.csv('Data/RSTChinMillDeer.csv',header=T)

RSTdata$MonthDay <-format(as.Date(RSTdata$Date), format="%m-%d")
RSTdata  <- RSTdata %>%
            mutate(Type = case_when(Month > 9 & Length > 50 |
                          Month <= 2 & Length > 60 |
                          Month == 3 & Length > 76 |
                          Month == 4 & Day >= 1 & Day <15 & Length > 85 |
                          Month == 4 & Day >= 15 & Day <30 & Length > 95 |
                          Month >= 4 & Month < 7 & Length > 100 ~ 'Yearling',
                          Month >= 1 & Month < 3 & Length > 45 & Length <= 60 |
                          Month >= 3 & Month < 7 & Length > 45 & Length <= 100  ~ 'Smolt',
                          Month > 10 & Length <= 45 |
                          Month <= 6 & Length <= 45  ~ 'Fry'))

# Plot RST data grouped by strategy ----------------------------------------------------------------------
# RST data summarized by month and day
date_ord <-  seq(as.Date("2001-10-01"), as.Date("2002-07-30"), by="days")
date_ord <- format(date_ord, format="%m-%d")


pRST <- ggplot(RSTdata,aes(x=factor(MonthDay,levels = date_ord),
                         y=Length, colour = Type))+
  geom_point()+xlab("")+ylab("Length (mm)")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_discrete(breaks=c("10-03","11-01" ,"12-01","01-01","02-01","03-01","04-01","05-01","06-01"),
                   labels=c("October","November","December","January","February","March", "April","May","June"))+
  scale_y_continuous(breaks=seq(20,160,20),limits = c(20, 150))+
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        text = element_text(size=25),
        legend.position = "none")+
  scale_color_manual(values = wes_palette("Darjeeling1", n = 3))

pRST


#----------------------------------------#
# Fish Fork Length (FL) back-calculation #
#----------------------------------------#

# Calculation of FL at natal and FW exit based on OR ------------------------------------------------

### Broken stick FL calibration model
calib = read.csv("Data/OR_FL_FINALforR.csv")
calib = subset(calib, select=c("Sample_ID","OR","FL"))
FL <- calib$FL; OR <- calib$OR
forced.intercept <- 30
res.lm <- lm(I(FL-forced.intercept) ~ 0 + OR)
res.bs <- segmented::segmented(res.lm, seg.Z = ~ 0 + OR)
FL.bs <- predict(res.bs) + forced.intercept

### Application on Mill/Deer Creek population
ORNataldata <- data.frame(OR=as.double(NatalExitdata_complete$NatalExit_OR_final))
NatalExitdata_complete$NatalExitFL <- segmented::predict.segmented(res.bs,newdata = ORNataldata) + forced.intercept

ORFWdata <- data.frame(OR=as.double(FWExitdata_complete$FWExit_OR_estimated))
FWExitdata_complete$FWExitFL <- segmented::predict.segmented(res.bs,newdata = ORFWdata) + forced.intercept

### FL summary statistics for Supp Mat
NatalExitFL.summary <- NatalExitdata_complete %>% group_by(reartype) %>%
                       dplyr::summarise(Mean = round(mean(NatalExitFL, na.rm=TRUE)),
                                        Max = round(max(NatalExitFL, na.rm=TRUE)),
                                        Min = round(min(NatalExitFL, na.rm=TRUE)),
                                        Sd = round(sd(NatalExitFL, na.rm=TRUE)))

FWExitFL.summary <- FWExitdata_complete %>% group_by(reartype) %>%
                    dplyr::summarise(Mean = round(mean(FWExitFL, na.rm=TRUE)),
                                     Max = round(max(FWExitFL, na.rm=TRUE)),
                                     Min = round(min(FWExitFL, na.rm=TRUE)),
                                     Sd = round(sd(FWExitFL, na.rm=TRUE)))

# Figure S2 for Supp Mat -------------------------------------------------------------

pS2b <- ggplot() + geom_point(data = calib,aes(OR,FL),col="black")+
  geom_point(data=NatalExitdata_complete,aes(NatalExit_OR_final,NatalExitFL),col='red',size=2)+
  labs(y=(expression(paste('Fork Length (mm)'))), fill="")+
  labs(x= "Otolith radius (µm)")+
  scale_y_continuous(limits=c(30,150),
                     breaks=c(30,60,90,120,150))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=25),legend.title = element_blank(),
        legend.position = "none", plot.title = element_text(face = "bold"))

pS2b

pS2c  <- ggplot() + geom_point(data = calib,aes(OR,FL),col="black")+
  geom_point(data=FWExitdata_complete,aes(FWExit_OR_estimated,FWExitFL),col='red',size=2)+
  labs(x= "Otolith radius (µm)",y="")+
  scale_y_continuous(limits=c(30,150),breaks=c(30,60,90,120,150))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=25),legend.title = element_blank(),
        legend.position = "none", plot.title = element_text(face = "bold"))

pS2c


pS2d <- ggplot(NatalExitdata_complete)+
  geom_density(aes(x= as.double(NatalExitFL),y=..count../sum(..count..),
                   fill = reartype,color=reartype),alpha=0.3,adjust=1.2, size=1)+  
  scale_x_continuous(limits=c(20,160), breaks=seq(30,150,30))+
  scale_y_continuous(limits=c(0,0.007), breaks=seq(0,0.006,0.002))+
  labs(x=(expression(paste('Fork Length (mm)'))), fill="")+
  labs(y="Density")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=25),legend.title = element_blank(),
        legend.position = "none", plot.title = element_text(face = "bold"))+
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 3))+  
  scale_color_manual(values = wes_palette("Darjeeling1", n = 3))

pS2d 

pS2e <- ggplot(FWExitdata_complete)+
  geom_density(aes(x=as.double(FWExitFL),y=..count../sum(..count..),fill = reartype,color=reartype),
               alpha=0.3,adjust=1.2, size=1) +
  labs(x=(expression(paste('Fork Length (mm)'))), fill="",y="")+
  scale_x_continuous(limits=c(50,160), breaks=seq(50,150,20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=25),legend.title = element_blank(),
        legend.position = "none", plot.title = element_text(face = "bold"))+
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 3))+  
  scale_color_manual(values = wes_palette("Darjeeling1", n = 3))

pS2e


pS2f <- ggplot(NatalExitdata_complete)+
  geom_density(aes(x=NatalExit_IncNum,y=..count../sum(..count..), 
                   fill = reartype,color=reartype),alpha=0.3,adjust=1.4,size=1) +
  scale_x_continuous(limits=c(0,300), breaks=seq(0,300,50))+
  scale_y_continuous(limits=c(0,0.005), breaks=seq(0,0.005,0.001))+
  labs(x='Increment Number')+
  ylab("Density") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=25),legend.title = element_blank(),
        legend.position = "none", plot.title = element_text(face = "bold"))+
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 3))+  
  scale_color_manual(values = wes_palette("Darjeeling1", n = 3))

pS2f

pS2g <- ggplot(FWExitdata_complete)+
  geom_density(aes(x=FWExit_IncNum,y=..count../sum(..count..), 
                   fill = reartype,color=reartype),alpha=0.3,adjust=1.4, size=1) +
  scale_x_continuous(limits=c(20,350), breaks=seq(50,350,50))+
  labs(x='Increment Number')+
  ylab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=25),legend.title = element_blank(),
        legend.position = "none", plot.title = element_text(face = "bold"))+
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 3))+  
  scale_color_manual(values = wes_palette("Darjeeling1", n = 3))

pS2g


pS2 <- ggarrange(ggarrange(ggExtra::ggMarginal(pRST,type = 'density',adjust = 4,margins = 'y',
                                               groupColour = TRUE, groupFill = TRUE),
                           font.label=list(size =25,color="black"),labels ='a'),
                 ggarrange(pS2b, pS2c, labels = c("b", "c"),
                           font.label=list(size =25,color="black")),
                 ggarrange(pS2d,pS2e, labels=c("d","e"),
                           font.label=list(size =25,color="black")),
                 ggarrange(pS2f,pS2g, labels=c("f","g"),
                           font.label=list(size =25,color="black")),
                 nrow=4)

png("Figures/FigureS2.png", 
    family = "sans serif", width = 12, height= 20, units = "in", res =200)

pS2

dev.off()



