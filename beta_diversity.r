library(vegan)
library(ggplot2)
library(ggpubr)
library(emmeans)
library(lme4)

#read in dissimilarity data
#this contains dissimilarities between pairs of years within the same sampling location
#ex. dissimilarity between fungi in sample point 1 in 2015 and fungi in the same location in 2016
#these values were generated using vegdist() in vegan and beta.temp() in betapart and manually
#compiled into a table
data<-read.csv("coupled_dissims.csv")
first<-subset(data,years=="Y1516")
second<-subset(data,years=="Y1617")
third<-subset(data,years=="Y1517")

boxplot(first$BC.dist~first$substrate)
boxplot(second$BC.dist~second$substrate)
boxplot(third$BC.dist~third$substrate)

#remove un-needed inner bark/wood substrate
nowood<-subset(data,substrate!="inner_bark")


#omit the 2015-2017 comparison
ones<-subset(nowood,years!="Y1517")
ones<-droplevels(ones)

library(plyr)
ones$years<-revalue(ones$years,c("Y1516"="2015-2016","Y1617"="2016-2017"))
ones$substrate<-revalue(ones$substrate, c("dead_bryos"="dead bryophytes","live_bryos"=
                        "live bryophytes","outer_bark"="bark"))
ones$substrate<-factor(ones$substrate,c("live bryophytes","dead bryophytes","bark"))
colnames(ones)[colnames(ones)=="substrate"]<-"Substrate"

#relative abundance data with Bray-Curtis dissimilarity
ggplot(data=ones,aes(x=Substrate,y=BC.dist)) + theme_bw() +
  geom_boxplot() +
  facet_wrap(~years) + ylab("Paired Bray-Curtis dissimilarity") +
  theme(axis.title = element_text(size=11)) +
  theme(axis.text = element_text(size=10)) +
  scale_x_discrete(labels=c("live\nbryophytes", "dead\nbryophytes", "bark"))

#presence-absence data with Jaccard dissimilarity
ggplot(data=ones,aes(x=Substrate,y=jac.dist)) +
  geom_boxplot() +
  facet_wrap(~years)


#beta diversity partitioning with Jaccard
ggplot(data=ones,aes(x=substrate,y=beta.jtu)) +
  geom_boxplot() +
  facet_wrap(~years)
ggplot(data=ones,aes(x=substrate,y=beta.jne)) +
  geom_boxplot() +
  facet_wrap(~years)

Substrate <-ones$Substrate
years<-ones$years
Jac.turnover<-ones$beta.jtu
Jac.nestedness<-ones$beta.jne
Jac.total<-ones$beta.jac
jacs<-data.frame(Substrate,years,Jac.turnover,Jac.nestedness,Jac.total)

library(tidyverse)
long <- jacs %>% 
  gather(key=index,value=beta,-years,-Substrate)
long$index<-as.factor(long$index)

long$index<-revalue(long$index,c("Jac.turnover"="turnover","Jac.nestedness"=
                                   "nestedness","Jac.total"="total Jaccard"))
long$index<-factor(long$index,c("nestedness","turnover","total Jaccard"))
colnames(long)[colnames(long)=="index"]<-"Index"

ggplot(data=long,aes(x=Substrate,y=beta,color=Index)) + theme_bw() +
  geom_boxplot() + scale_color_viridis_d(option="viridis") +
  facet_wrap(~years) + ylab("Beta diversity") +
  scale_x_discrete(labels=c("live\nbryophytes", "dead\nbryophytes", "bark")) +
  theme(axis.text = element_text(size=10)) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size=10))



#mixed models
BCmod<-lmer(BC.dist ~ substrate*years + (1|BP/BPS), data = ones)
#this random effect is the same as + (1|BP) + (1|BPS:BP)
#significance test
car::Anova(BCmod, type=3)

Jacmod<-lmer(jac.dist ~ substrate*years + (1|BP/BPS), data = ones)
car::Anova(Jacmod, type=3)

emmeans(BCmod, pairwise ~ years)
emmeans(BCmod, pairwise ~ substrate)
emmeans(BCmod, pairwise ~ years | substrate)
emmeans(BCmod, pairwise ~ substrate | years)

emmeans(Jacmod, pairwise ~ years)
emmeans(Jacmod, pairwise ~ substrate)
emmeans(Jacmod, pairwise ~ years | substrate)
emmeans(Jacmod, pairwise ~ substrate | years)
