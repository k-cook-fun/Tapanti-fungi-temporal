library(vegan)
library(tidyverse)

#load data
data<-read.csv("nw_data.csv")
#separate by year
Y15<-subset(data,year=="Y15")
Y16<-subset(data,year=="Y16")
Y17<-subset(data,year=="Y17")


#calculate rates at which fungi persist from one year to the next

#2015 to 2016
#not all samples/substrates exist in both years
#subset out only those present in both years
ax<-intersect(Y15$BPS,Y16$BPS)
Y15a <- Y15 %>%
  filter(BPS %in% c(ax)) 
Y16a <- Y16 %>%
  filter(BPS %in% c(ax)) 

#get otu relative abundance for the 1st year
Y15a<-decostand(Y15a[,12:11421],"total")
#and presence-absence for the 2nd year
Y16a<-decostand(Y16a[,12:11421],"pa")

#samples are in same order in both years
#substrates are:
#dead bryos - 1-35
#live bryos - 36-61
#bark - 62-99

#if somethings is present in 2015, how often is it also in 2016?
#set all otus present in 2016 to 0.5
PA16a<-Y16a/2
#set outs present in 2015 to 1
PA15a<-decostand(Y15a,"pa")

#add them together
#0 = otu missing in both years
#1 = otu present in 2015 only
#0.5 = otu present in 2016 only
#1.5 = present in both years
total<-PA15a+PA16a
sum(total[,]==1.5)/(sum(total[,]==1)+sum(total[,]==1.5))

#by substrate type
sum(total[1:35,]==1.5)/(sum(total[1:35,]==1)+sum(total[1:35,]==1.5))
sum(total[36:61,]==1.5)/(sum(total[36:61,]==1)+sum(total[36:61,]==1.5))
sum(total[62:99,]==1.5)/(sum(total[62:99,]==1)+sum(total[62:99,]==1.5))

#how does this vary by initial abundance of the OTU?
#very low abundance - <=0.01%
#low abundance - >0.01%, <=0.1%
#medium - >0.1%,<=1%
#high - >1%,<=10%  
#very high - >10%

#very low
vlow15<-Y15a
vlow15[vlow15 > 0.0001] <- 0
vlow15[vlow15 > 0] <- 1
vlow<-vlow15+PA16a
sum(vlow[,]==1.5)/(sum(vlow[,]==1)+sum(vlow[,]==1.5))
#by substrate type
sum(vlow[1:35,]==1.5)/(sum(vlow[1:35,]==1)+sum(vlow[1:35,]==1.5))
sum(vlow[36:61,]==1.5)/(sum(vlow[36:61,]==1)+sum(vlow[36:61,]==1.5))
sum(vlow[62:99,]==1.5)/(sum(vlow[62:99,]==1)+sum(vlow[62:99,]==1.5))

#low
low15<-Y15a
low15[low15 <= 0.0001] <- 0 
low15[low15 > 0.001] <- 0 
low15[low15 > 0] <- 1
low<-low15+PA16a
sum(low[,]==1.5)/(sum(low[,]==1)+sum(low[,]==1.5))
#by substrate type
sum(low[1:35,]==1.5)/(sum(low[1:35,]==1)+sum(low[1:35,]==1.5))
sum(low[36:61,]==1.5)/(sum(low[36:61,]==1)+sum(low[36:61,]==1.5))
sum(low[62:99,]==1.5)/(sum(low[62:99,]==1)+sum(low[62:99,]==1.5))

#medium
med15<-Y15a
med15[med15 <= 0.001] <- 0 
med15[med15 > 0.01] <- 0 
med15[med15 > 0] <- 1
med<-med15+PA16a
sum(med[,]==1.5)/(sum(med[,]==1)+sum(med[,]==1.5))
#by substrate type
sum(med[1:35,]==1.5)/(sum(med[1:35,]==1)+sum(med[1:35,]==1.5))
sum(med[36:61,]==1.5)/(sum(med[36:61,]==1)+sum(med[36:61,]==1.5))
sum(med[62:99,]==1.5)/(sum(med[62:99,]==1)+sum(med[62:99,]==1.5))

#high
high15<-Y15a
high15[high15 <= 0.01] <- 0 
high15[high15 > 0.1] <- 0 
high15[high15 > 0] <- 1
high<-high15+PA16a
sum(high[,]==1.5)/(sum(high[,]==1)+sum(high[,]==1.5))
#by substrate type
sum(high[1:35,]==1.5)/(sum(high[1:35,]==1)+sum(high[1:35,]==1.5))
sum(high[36:61,]==1.5)/(sum(high[36:61,]==1)+sum(high[36:61,]==1.5))
sum(high[62:99,]==1.5)/(sum(high[62:99,]==1)+sum(high[62:99,]==1.5))

#very high
vhigh15<-Y15a
vhigh15<-ifelse(Y15a>0.1, 1, 0)
vhigh<-vhigh15+PA16a
sum(vhigh[,]==1.5)/(sum(vhigh[,]==1)+sum(vhigh[,]==1.5))
#by substrate type
sum(vhigh[1:35,]==1.5)/(sum(vhigh[1:35,]==1)+sum(vhigh[1:35,]==1.5))
sum(vhigh[36:61,]==1.5)/(sum(vhigh[36:61,]==1)+sum(vhigh[36:61,]==1.5))
sum(vhigh[62:99,]==1.5)/(sum(vhigh[62:99,]==1)+sum(vhigh[62:99,]==1.5))







#2016-2017
ay<-intersect(Y16$BPS,Y17$BPS)
Y16b <- Y16 %>%
  filter(BPS %in% c(ay)) 
Y17b <- Y17 %>%
  filter(BPS %in% c(ay)) 

#get otu relative abundance for the 1st year
Y16b<-decostand(Y16b[,12:11421],"total")
#and presence-absence for the 2nd year
Y17b<-decostand(Y17b[,12:11421],"pa")

#set all otus present in 2017 to 0.5
PA17b<-Y17b/2
#set outs present in 2016 to 1
PA16b<-decostand(Y16b,"pa")

#samples are in same order in both years
#substrates are:
#dead bryos - 1-29
#live bryos - 30-49
#bark - 50-81

#add them together
#0 = otu missing in both years
#1 = otu present in 2016 only
#0.5 = otu present in 2017 only
#1.5 = present in both years
total<-PA17b+PA16b
sum(total[,]==1.5)/(sum(total[,]==1)+sum(total[,]==1.5))
#by substrate type
sum(total[1:29,]==1.5)/(sum(total[1:29,]==1)+sum(total[1:29,]==1.5))
sum(total[30:49,]==1.5)/(sum(total[30:49,]==1)+sum(total[30:49,]==1.5))
sum(total[50:81,]==1.5)/(sum(total[50:81,]==1)+sum(total[50:81,]==1.5))

#very low
vlow16<-Y16b
vlow16[vlow16 > 0.0001] <- 0
vlow16[vlow16 > 0] <- 1
vlow<-vlow16+PA17b
sum(vlow[,]==1.5)/(sum(vlow[,]==1)+sum(vlow[,]==1.5))
#by substrate type
sum(vlow[1:29,]==1.5)/(sum(vlow[1:29,]==1)+sum(vlow[1:29,]==1.5))
sum(vlow[30:49,]==1.5)/(sum(vlow[30:49,]==1)+sum(vlow[30:49,]==1.5))
sum(vlow[50:81,]==1.5)/(sum(vlow[50:81,]==1)+sum(vlow[50:81,]==1.5))

#low
low16<-Y16b
low16[low16 <= 0.0001] <- 0 
low16[low16 > 0.001] <- 0 
low16[low16 > 0] <- 1
low<-low16+PA17b
sum(low[,]==1.5)/(sum(low[,]==1)+sum(low[,]==1.5))
#by substrate type
sum(low[1:29,]==1.5)/(sum(low[1:29,]==1)+sum(low[1:29,]==1.5))
sum(low[30:49,]==1.5)/(sum(low[30:49,]==1)+sum(low[30:49,]==1.5))
sum(low[50:81,]==1.5)/(sum(low[50:81,]==1)+sum(low[50:81,]==1.5))

#medium
med16<-Y16b
med16[med16 <= 0.001] <- 0 
med16[med16 > 0.01] <- 0 
med16[med16 > 0] <- 1
med<-med16+PA17b
sum(med[,]==1.5)/(sum(med[,]==1)+sum(med[,]==1.5))
#by substrate type
sum(med[1:29,]==1.5)/(sum(med[1:29,]==1)+sum(med[1:29,]==1.5))
sum(med[30:49,]==1.5)/(sum(med[30:49,]==1)+sum(med[30:49,]==1.5))
sum(med[50:81,]==1.5)/(sum(med[50:81,]==1)+sum(med[50:81,]==1.5))

#high
high16<-Y16b
high16[high16 <= 0.01] <- 0 
high16[high16 > 0.1] <- 0 
high16[high16 > 0] <- 1
high<-high16+PA17b
sum(high[,]==1.5)/(sum(high[,]==1)+sum(high[,]==1.5))
#by substrate type
sum(high[1:29,]==1.5)/(sum(high[1:29,]==1)+sum(high[1:29,]==1.5))
sum(high[30:49,]==1.5)/(sum(high[30:49,]==1)+sum(high[30:49,]==1.5))
sum(high[50:81,]==1.5)/(sum(high[50:81,]==1)+sum(high[50:81,]==1.5))

#very high
vhigh16<-Y16b
vhigh16<-ifelse(Y16b>0.1, 1, 0)
vhigh<-vhigh16+PA17b
sum(vhigh[,]==1.5)/(sum(vhigh[,]==1)+sum(vhigh[,]==1.5))
#by substrate type
sum(vhigh[1:29,]==1.5)/(sum(vhigh[1:29,]==1)+sum(vhigh[1:29,]==1.5))
sum(vhigh[30:49,]==1.5)/(sum(vhigh[30:49,]==1)+sum(vhigh[30:49,]==1.5))
sum(vhigh[50:81,]==1.5)/(sum(vhigh[50:81,]==1)+sum(vhigh[50:81,]==1.5))





#2015-2017
az<-intersect(Y15$BPS,Y17$BPS)
Y15c <- Y15 %>%
  filter(BPS %in% c(az)) 
Y17c <- Y17 %>%
  filter(BPS %in% c(az)) 

#get otu relative abundance for the 1st year
Y15c<-decostand(Y15c[,12:11421],"total")
#and presence-absence for the 2nd year
Y17c<-decostand(Y17c[,12:11421],"pa")

#set all otus present in 2017 to 0.5
PA17c<-Y17c/2
#set outs present in 2015 to 1
PA15c<-decostand(Y15c,"pa")

#add them together
#0 = otu missing in both years
#1 = otu present in 2015 only
#0.5 = otu present in 2017 only
#1.5 = present in both years
total<-PA17c+PA15c
sum(total[,]==1.5)/(sum(total[,]==1)+sum(total[,]==1.5))
#by substrate type
sum(total[1:39,]==1.5)/(sum(total[1:39,]==1)+sum(total[1:39,]==1.5))
sum(total[40:75,]==1.5)/(sum(total[40:75,]==1)+sum(total[40:75,]==1.5))
sum(total[76:123,]==1.5)/(sum(total[76:123,]==1)+sum(total[76:123,]==1.5))

#samples are in same order in both years
#substrates are:
#dead bryos - 1-39
#live bryos - 40-75
#bark - 76-123

#very low
vlow15<-Y15c
vlow15[vlow15 > 0.0001] <- 0
vlow15[vlow15 > 0] <- 1
vlow<-vlow15+PA17c
sum(vlow[,]==1.5)/(sum(vlow[,]==1)+sum(vlow[,]==1.5))
#by substrate type
sum(vlow[1:39,]==1.5)/(sum(vlow[1:39,]==1)+sum(vlow[1:39,]==1.5))
sum(vlow[40:75,]==1.5)/(sum(vlow[40:75,]==1)+sum(vlow[40:75,]==1.5))
sum(vlow[76:123,]==1.5)/(sum(vlow[76:123,]==1)+sum(vlow[76:123,]==1.5))

#low
low15<-Y15c
low15[low15 <= 0.0001] <- 0 
low15[low15 > 0.001] <- 0 
low15[low15 > 0] <- 1
low<-low15+PA17c
sum(low[,]==1.5)/(sum(low[,]==1)+sum(low[,]==1.5))
#by substrate type
sum(low[1:39,]==1.5)/(sum(low[1:39,]==1)+sum(low[1:39,]==1.5))
sum(low[40:75,]==1.5)/(sum(low[40:75,]==1)+sum(low[40:75,]==1.5))
sum(low[76:123,]==1.5)/(sum(low[76:123,]==1)+sum(low[76:123,]==1.5))

#medium
med15<-Y15c
med15[med15 <= 0.001] <- 0 
med15[med15 > 0.01] <- 0 
med15[med15 > 0] <- 1
med<-med15+PA17c
sum(med[,]==1.5)/(sum(med[,]==1)+sum(med[,]==1.5))
#by substrate type
sum(med[1:39,]==1.5)/(sum(med[1:39,]==1)+sum(med[1:39,]==1.5))
sum(med[40:75,]==1.5)/(sum(med[40:75,]==1)+sum(med[40:75,]==1.5))
sum(med[76:123,]==1.5)/(sum(med[76:123,]==1)+sum(med[76:123,]==1.5))

#high
high15<-Y15c
high15[high15 <= 0.01] <- 0 
high15[high15 > 0.1] <- 0 
high15[high15 > 0] <- 1
high<-high15+PA17c
sum(high[,]==1.5)/(sum(high[,]==1)+sum(high[,]==1.5))
#by substrate type
sum(high[1:39,]==1.5)/(sum(high[1:39,]==1)+sum(high[1:39,]==1.5))
sum(high[40:75,]==1.5)/(sum(high[40:75,]==1)+sum(high[40:75,]==1.5))
sum(high[76:123,]==1.5)/(sum(high[76:123,]==1)+sum(high[76:123,]==1.5))

#very high
vhigh15<-Y15c
vhigh15<-ifelse(Y15c>0.1, 1, 0)
vhigh<-vhigh15+PA17c
sum(vhigh[,]==1.5)/(sum(vhigh[,]==1)+sum(vhigh[,]==1.5))
#by substrate type
sum(vhigh[1:39,]==1.5)/(sum(vhigh[1:39,]==1)+sum(vhigh[1:39,]==1.5))
sum(vhigh[40:75,]==1.5)/(sum(vhigh[40:75,]==1)+sum(vhigh[40:75,]==1.5))
sum(vhigh[76:123,]==1.5)/(sum(vhigh[76:123,]==1)+sum(vhigh[76:123,]==1.5))
