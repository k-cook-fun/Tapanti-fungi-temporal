library(vegan)
library(ggplot2)
library(emmeans)
library(lme4)

#load data set, wood already removed
data<-read.csv("nw_data.csv")

otus<-data[,12:11421]
subs<-data$substrate
year<-data$year
BP<-data$BP
BPS<-data$BPS

#rename and reorder factors so they are in better order
library(plyr)
year<-revalue(year, c("Y15"="2015", "Y16"="2016", "Y17"="2017"))
subs<-revalue(subs, c("dead_bryos"="dead bryophytes","live_bryos"=
                        "live bryophytes","outer_bark"="bark"))
subs<-factor(subs,c("live bryophytes","dead bryophytes","bark"))


#get expected richness at 1000 reads
alpha<-rarefy(otus,sample=1000)

#for shannon, do repeated rarefactions  via for loop
#warning: SLOW
#first make empty data frames for the data
x<-data.frame(row.names=(1:408))
y<-data.frame(row.names=(1:408))

set.seed(62323222)
for (i in 1:1000) {
  rared<-rrarefy(otus,sample=1000)
  shan<-diversity(rared,"shannon")
  y<-cbind(y,shan)
}

#get average
Shannon<-rowMeans(y)

comb<-as.data.frame(cbind(Shannon,alpha))
comb<-cbind(comb,subs,year,BP,BPS)
colnames(comb)[colnames(comb)=="year"]<-"Year"


#boxplots
ggplot(comb,aes(x=subs,y=alpha,fill=Year)) + geom_boxplot() +
  xlab("Substrate") + ylab("OTU Richness")

ggplot(comb,aes(x=subs,y=Shannon,fill=Year)) + geom_boxplot() +
  xlab("Substrate") + ylab("Shannon")


#mixed effects models
#OTU richness
almod<-lmer(alpha~subs*year + (1|BP/BPS),data=comb)
#this random effect is the same as + (1|BP) + (1|BPS:BP)
plot(almod)
#high heterogeneity of variance
#use log transformation to improve this
logmod<-lmer(log(alpha) ~ subs*year + (1|BP/BPS), data = comb)
plot(logmod)
qqnorm(residuals(logmod))
hist(residuals(logmod))
car::Anova(logmod)
emmeans(logmod, pairwise ~ year | subs)
emmeans(logmod, pairwise ~ subs | year)


#Shannon
shamod <- lmer(Shannon ~ subs*Year + (1|BP/BPS), data = comb)
plot(shamod)
qqnorm(residuals(shamod))
hist(residuals(shamod))
car::Anova(shamod)
emmeans(shamod, pairwise ~Year | subs)
emmeans(shamod, pairwise ~ subs | Year)
