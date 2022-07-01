library(vegan)
library(adespatial)

#load otu data
data<-read.csv("nw_data.csv")
meta<-data[,1:11]
subs<-data[,2]
year<-data[,4]
otus<-data[,12:11421]

#spatial and temporal distances
spaces<-read.csv("nw_dists.csv",header = F)

#convert otus to relative abundance
rel<-decostand(otus,method="total")
#and get presence-absence (PA)
PA <- decostand(otus,method="pa")

#get PCNM vectors for spatial distances
sp.PCNM<-vegan::pcnm(spaces)
dmin<-sp.PCNM$threshold
sp.PCNM.data<-as.data.frame(sp.PCNM$vectors)
nb.ev<-length(which(sp.PCNM$values>0))


#RDA on relative abundance OTUs with spatial vectors + significance test
sp.rda<-rda(rel,sp.PCNM$vectors)
anova.cca(sp.rda)
R1<-RsquareAdj(sp.rda)$adj.r.squared


#RDA for PA OTUs with spatal vectors + significance test
pa.rda<-rda(PA,sp.PCNM$vectors)
anova.cca(pa.rda)
Ra<-RsquareAdj(pa.rda)$adj.r.squared


#forward select pcnm vectors
#WARNING - this step can be very slow
#nd.pcnm.fwd<-forward.sel(rel,as.matrix(sp.PCNM$vectors),adjR2thresh=R1)
#pcnm.sign<-sort(nd.pcnm.fwd[,2])
#these are the vectors I got, use the line below to skip the slow forward selection
pcnm.red<-sp.PCNM$vectors[,c(1,2,3,4,5,6,7,8,9,10,13,14,19,20,22,23,25,27,28,29,30,31,32,94)]

#repeat for PA data
#pa.fwd<-forward.sel(PA,as.matrix(sp.PCNM$vectors),adjR2thresh=Ra)
pcnm.pa<-sp.PCNM$vectors[,c(2,1,3,4,79,6,84,83,13,99,5,23,121,93,118,29,48)]

#rda with new subset of spatial vectors + significance tests
nd.rda2<-rda(rel ~ ., data=as.data.frame(pcnm.red))
anova.cca(nd.rda2)
pa.rda2<-rda(PA ~ ., data=as.data.frame(pcnm.pa))
anova.cca(pa.rda2)



#variance partitioning
#recode substrate variable as binary
substrate<-model.matrix(~subs)[,-1]

#test significance of substrate and year
sub.rda<-rda(rel,substrate)
anova.cca(sub.rda)
year.rda<-rda(rel,year)
anova.cca(year.rda)

#variace partitioning with forward selected PCNMs, substrates, and year
varpart<-varpart(rel,substrate,year,pcnm.red)
plot(varpart,digits=2)

#test significance of each partition
#WARNING - can be slow to run
sub.test<-anova.cca(rda(rel,substrate,cbind(year,pcnm.red)))
space.test<-anova.cca(rda(rel,pcnm.red,cbind(year,substrate)))
year.test<-anova.cca(rda(rel,year,cbind(substrate,pcnm.red)))


#repeat for PA
pa.sub<-rda(PA,substrate)
anova.cca(pa.sub)
pa.year<-rda(PA,year)
anova.cca(pa.year)
pa.varpart<-varpart(PA,substrate,year,pcnm.pa)
plot(pa.varpart,digits=2)

#test significance of each partition
sub.test.pa<-anova.cca(rda(PA,substrate,cbind(year,pcnm.pa)))
space.test.pa<-anova.cca(rda(PA,pcnm.pa,cbind(year,substrate)))
year.test.pa<-anova.cca(rda(PA,year,cbind(substrate,pcnm.pa)))