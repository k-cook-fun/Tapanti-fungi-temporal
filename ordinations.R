library(vegan)
library(ggplot2)
library(phyloseq)
library(microbiome)
library(grid)

#load OTU table formatted for phyloseq
otus<-read.csv("phylo_otu_tab.csv",header=TRUE,sep=",",row.names=1)
motus<-as.matrix(otus)

#read taxonomy
taxonomy<-read.csv("phylo_taxonomy.csv",header=TRUE,sep=",",row.names=1,na.strings="")
mtax<-as.matrix(taxonomy)

#read sample data
samples<-read.csv("phylo_metadata.csv",header=TRUE,row.names=1)

#merge these into a phyloseq object
otutab<-otu_table(motus,taxa_are_rows=TRUE)
taxtab<-tax_table(mtax)
samtab<-sample_data(samples)
physeq<-phyloseq(otutab,taxtab,samtab)

#adjust names for substrate and year variables
library(plyr)
sample_data(physeq)$substrate<-revalue(sample_data(physeq)$substrate,
                              c("dead_bryos"="dead bryophytes",
                              "live_bryos"="live bryophytes","outer_bark"="bark"))
sample_data(physeq)$year<-revalue(sample_data(physeq)$year,
                              c("Y15"="2015","Y16"="2016","Y17"="2017"))

#make absolute and relative abundance tables
otu.absolute <- abundances(physeq)
otu.relative <- phyloseq::transform_sample_counts(physeq, function(x) x/sum(x))

#drop wood/inner bark samples that are not used in this study
nw <- subset_samples(otu.relative,substrate!="inner_bark")
colnames(sample_data(nw))[colnames(sample_data(nw))=="substrate"]<-"Substrate"
colnames(sample_data(nw))[colnames(sample_data(nw))=="year"]<-"Year"


#ordination of samples; color by substrate type, shape by sampling year
set.seed(54648644)
ord <- ordinate(nw, "NMDS", "bray",try=100)
plot.ord<-plot_ordination(nw, ord, type="samples", color="Substrate", shape="Year") +
  theme_bw()

plot.ord + scale_color_viridis_d(option="viridis") +
  scale_shape_manual(values=c(16,3,17)) +
  geom_point(size=2) + theme(axis.title = element_text(size=13)) +
  theme(axis.text = element_text(size=13)) +
  theme(legend.text = element_text(size=13)) +
  theme(legend.title = element_text(size=13)) +
  theme(legend.key.size = unit(1,"cm"))



#ordinations of individual substrates to better see patterns by year
BRK<-subset_samples(nw,Substrate=="bark")
DMO<-subset_samples(nw,Substrate=="dead bryophytes")
LMO<-subset_samples(nw,Substrate=="live bryophytes")

set.seed(20222760)
BKord<-ordinate(BRK,"NMDS","bray",try=100)
DMord<-ordinate(DMO,"NMDS","bray",try=100)
LMord<-ordinate(LMO,"NMDS","bray",try=100)

BKplot<-plot_ordination(BRK, BKord, type="samples", color="Year", shape="Year") + theme_bw() + ggtitle("Bark") + stat_ellipse()
DMplot<-plot_ordination(DMO, DMord, type="samples", color="Year", shape="Year") + theme_bw() + ggtitle("Dead Bryophytes") + stat_ellipse()
LMplot<-plot_ordination(LMO, LMord, type="samples", color="Year", shape="Year") + theme_bw() + ggtitle("Live Bryophytes") + stat_ellipse()

BKplot + scale_color_viridis_d(option="viridis") +
  scale_shape_manual(values=c(16,3,17)) +
  geom_point(size=2) + theme(axis.title = element_text(size=13)) +
  theme(axis.text = element_text(size=13)) +
  theme(legend.text = element_text(size=13)) +
  theme(legend.title = element_text(size=13)) +
  theme(legend.key.size = unit(1,"cm")) +
  annotate("text",x=-0.685,y=0.55,label="a",size=7) +
  coord_cartesian(xlim=c(-0.45,0.85),ylim=c(-0.4,0.45),clip="off")

DMplot + scale_color_viridis_d(option="viridis") +
  scale_shape_manual(values=c(16,3,17)) +
  geom_point(size=2) + theme(axis.title = element_text(size=13)) +
  theme(axis.text = element_text(size=13)) +
  theme(legend.text = element_text(size=13)) +
  theme(legend.title = element_text(size=13)) +
  theme(legend.key.size = unit(1,"cm")) +
  annotate("text",x=-0.74,y=0.75,label="b",size=7) +
  coord_cartesian(xlim=c(-0.54,0.6),ylim=c(-0.61,0.61),clip="off")

LMplot + scale_color_viridis_d(option="viridis") +
  scale_shape_manual(values=c(16,3,17)) +
  geom_point(size=2) + theme(axis.title = element_text(size=13)) +
  theme(axis.text = element_text(size=13)) +
  theme(legend.text = element_text(size=13)) +
  theme(legend.title = element_text(size=13)) +
  theme(legend.key.size = unit(1,"cm")) +
  annotate("text",x=-0.236,y=0.205,label="c",size=7) +
  coord_cartesian(xlim=c(-0.18,0.15),ylim=c(-0.21,0.16),clip="off")
