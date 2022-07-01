library(vegan)
library(ggplot2)

#read and setup OTU data
data<-read.csv("2015_data.csv")

names<-data[,1]
substrate<-data[,2]
branch<-data[,3]
point<-data[,4]
BP<-data[,5]
trunk.dist<-data[,6]
rel.dist<-data[,7]
meta<-data[,1:7]

#order substrates from "top" to "bottom"
substrate<-factor(substrate,c("LB","DB","OB","IB"))

otus<-data[,8:5784]

#separate by substrate
#OB - bark, DB - dead bryophytes, LB - live bryophytes, IB - inner bark/wood (not used in analysis)
OB<-otus[1:121,]
DB<-otus[122:238,]
LB<-otus[239:349,]
IB<-otus[350:380,]

#read in spatial data
spaces<-read.csv("all_dists.csv",header=FALSE)
#and separate by substrate
bk.dist<-spaces[1:121,1:121]
dm.dist<-spaces[122:238,122:238]
lm.dist<-spaces[239:349,239:349]

#convert OTUs to relative abundance and create dissimilarity matrices using Bray-Curtis (BC)
BRKmat<-vegdist(decostand(OB,"total"),"bray")
LMOmat<-vegdist(decostand(LB,"total"),"bray")
DMOmat<-vegdist(decostand(DB,"total"),"bray")


#make data frames of BC dissimilarities and distances for each substrate
BRKall<-cbind(as.dist(bk.dist,upper=FALSE,diag=FALSE),BRKmat)
DMOall<-cbind(as.dist(dm.dist,upper=FALSE,diag=FALSE),DMOmat)
LMOall<-cbind(as.dist(lm.dist,upper=FALSE,diag=FALSE),LMOmat)

#combine them all
All<-as.data.frame(rbind(BRKall,DMOall,LMOall))
All$Substrate<-c(rep("bark",length(BRKall[,1])),
                 rep("dead bryophytes",length(DMOall[,1])),
                 rep("live bryophytes",length(LMOall[,1])))

#optionally turn into similarity instead of dissimilarity
#All[,2]<-(1-All[,2])

#re-subset by substrate for plotting separately
bark<-subset(All,Substrate=="bark")
deadbry<-subset(All,Substrate=="dead bryophytes")
livebry<-subset(All,Substrate=="live bryophytes")

#full distance decay using all distances, bark
bkplot<-ggplot(bark,aes(x=V1,y=BRKmat))
bkplot + geom_point(color="#00B8E7",alpha=0.7) + 
  geom_smooth(method="loess",se=FALSE) + theme_classic() +
  labs(title="Bark",x="Distance in cm",y="Bray-Curtis Dissimilarity")

#dead bryophytes
dbplot<-ggplot(deadbry,aes(x=V1,y=BRKmat))
dbplot + geom_point(color="#E68613",alpha=0.7) + 
  geom_smooth(method="loess",se=FALSE) + theme_classic() +
  labs(title="Dead Bryophytes",x="Distance in cm",y="Bray-Curtis Dissimilarity")

#live bryophytes
lbplot<-ggplot(livebry,aes(x=V1,y=BRKmat))
lbplot + geom_point(color="#0CB702",alpha=0.7) + 
  geom_smooth(method="loess",se=FALSE) + theme_classic() +
  labs(title="Live Bryophytes",x="Distance in cm",y="Bray-Curtis Dissimilarity")

#or plot all together, messy
dd<-ggplot(All,aes(x=V1,y=BRKmat,col=Substrate))
dd + geom_point(alpha=0.7) + 
  geom_smooth(method="loess") + facet_grid(vars(Substrate)) +
  scale_color_manual(values=c("royalblue4","darkmagenta","olivedrab4")) +
  labs(x="Distance in cm",y="Bray-Curtis Dissimilarity")


#extract the smallest distances <=10 cm and plot all substrates together
All10<-subset(All,All[,1] <=10)
dplot<-ggplot(All10,aes(x=V1,y=BRKmat,col=Substrate,shape=Substrate))
dplot + geom_point(alpha=0.8) + geom_smooth(method="lm",se=FALSE) + 
  labs(x="Distance in cm",y="Bray-Curtis Dissiilarity") +
  scale_color_manual(values=c("#00B8E7","#E68613","#0CB702")) + 
  theme_classic() + theme(legend.position = "top")
  
#get formulas for lines
modA<-lm(All10[,1]~All10[,2]*All10[,3])
summary(modA)
#can use these formulas to calculate dissimilarity at a distance of 0 cm (y-intercept)
#and the spatial distance at which inter-annual dissimilarity equals spatial dissimilarity

#what are the average dissimilarity values at large distances?
All100<-subset(All,All[,1]>=100)
summary(subset(All100,Substrate=="bark")[,2])
summary(subset(All100,Substrate=="dead bryophytes")[,2])
summary(subset(All100,Substrate=="live bryophytes")[,2])






#create new distance decay graphs with presence-absence data and Jaccard dissimilarity
BRKj<-vegdist(decostand(OB,"pa"),"jaccard")
LMOj<-vegdist(decostand(LB,"pa"),"jaccard")
DMOj<-vegdist(decostand(DB,"pa"),"jaccard")

BRKpa<-cbind(as.dist(bk.dist,upper=FALSE,diag=FALSE),BRKj)
DMOpa<-cbind(as.dist(dm.dist,upper=FALSE,diag=FALSE),DMOj)
LMOpa<-cbind(as.dist(lm.dist,upper=FALSE,diag=FALSE),LMOj)

AllPA<-as.data.frame(rbind(BRKpa,DMOpa,LMOpa))
AllPA$Substrate<-c(rep("bark",length(BRKpa[,1])),
                 rep("dead bryophytes",length(DMOpa[,1])),
                 rep("live bryophytes",length(LMOpa[,1])))

PA50<-subset(AllPA,V1<=50)

pplot<-ggplot(PA50,aes(x=V1,y=BRKj,col=Substrate,shape=Substrate))
pplot + geom_point(alpha=0.8) + geom_smooth(method="lm",se=FALSE) + 
  labs(x="Distance in cm",y="Jaccard Dissiilarity") +
  scale_color_manual(values=c("#00B8E7","#E68613","#0CB702")) + 
  theme_classic() + theme(legend.position = "top")
