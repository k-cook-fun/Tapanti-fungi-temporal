library(vegan)

#load community data
data<-read.csv("nw_data.csv")

otus<-data[,12:11421]
row.names(otus)<-data[,1]

#do two trasformations, presence-absence and relative abundance
PA<-decostand(otus,method="pa")
relative<-decostand(otus,method="total")

#what are the maximum relative abundances of OTUs per sample?
maxed<-apply(relative,FUN=max,MARGIN=1)
hist(maxed,xlab="Relative Abundance",main="Most abundant OTU per Sample")

#calculate the number of samples each OTU occurs in
freqs<-colSums(PA)
hist(freqs)
ecdf(freqs)(3)
quantile(freqs,0.5)


#abundance
#transform all absences (0s in the OTU table) to NA
noz<-relative
noz[noz==0]<-NA
#calculate mean non-zero relative abundance for each OTU
abds<-colMeans(noz,na.rm=TRUE)
hist(abds)

#make a list of all non-zero abundances
abds2<-as.vector(as.matrix(noz))
abds2<-abds2[!is.na(abds2)]
hist(abds2)
hist(log(abds2),xlab="log Relative Abundance",main="All OTU Occurrences")

#find percentile of a given relative abundance level
ecdf(abds2)(0.05)
