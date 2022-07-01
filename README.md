# Tapanti-fungi-temporal

This repository contains scripts and accompanying data tables used in the study "Inter-annual persistence of canopy fungi driven by abundance despite high spatial turnover".

The purpose of this research project was to better understand how fungal communities change over time, how these changes differ among growth substrates, how temporal changes compare to spatial turnover within a single time point, and if a fungus' abundance at one time can be used to predict its presence in the furture.

This study took place in a low montane tropical rainforest in Parque Nacional Tapanti in central Costa Rica and focused on fungi in a tree canopy environment.  Samples were collected from five adjaced _Sauraia_ tree branches with a cork borer and manually dissected into up to three substrates: living bryphytes, dead bryophytes, and host tree bark.  Sampling first took place in July 2015, and the exact locations were revisited in 2016 and 2017 and resampled.  These repeated samples were taken immediately adjacent to the original ones.  High-throughput sequencing of the ITS region was used to determine fungi present in each sample.  Sequence reads can be access in NCBI GenBank SRA PRJNA762332.  The spatial study performed in the same system, to which data from this temporal study were compared, can be accessed here:  https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.16358




Summary of scripts:

uparse_seq_processing.txt - 
Scripts used for sequence processing and creation of OTUs (operational taxonomic units, roughly analogous to species) and OTU tables.
  
alpha_diversity.r - 
Calculating species richness and the Shannon index, models for contrasting them among years and substrates
  
beta_diversity.r - 
Assessing beta diversity using relative abundance/Bray-Curtis dissimilarity and presence-absence/Jaccard dissimilarity, comparisons among years and substrates
  
distance_decay.R - 
Creation of distance decay curves using spatial data, comparing temporal to spatial turnover
  
frequency_abundance.r - 
Assessing the frequency (commonness, number of samples a taxon is observed in) and abundance (relative abundance, % of a given sample that the taxon makes up) of OTUs
  
ordination.R - 
Constructing NMDS ordinations of samples to visualize differences in fungal community composition among years and substrates
  
reoccurrence_rates.R - 
Calculate how often OTUs persist from one year to the next and how relative abundance of OTUs impacts this
  
variance_partitioning.r - 
Variance partition analysis to determine effects of sampling year, substrate, and spatial relationships on fungal community composition
  
  
  
The remaining files are data tables used in the above R scripts.  These include spatial distance matrices (all_dists.csv and nw_dists.csv), in which values are spatial distances between pairs of samples measured in centimeters, a table of dissimilarity values (coupled_dissims.csv) between pairs of samples collected at the same location in different years, and OTU tables (species by site matrices) with accompanying sample metadata and OTU taxonomic information (nw_data.csv, phylo_otu_tab.csv, phylo_metadata.csv, phylo_taxonomy.csv).



Sample metadata includes the following:

sample - 
Unique sample name, incorporates year, branch, point, and substrate

substrate - 
Substrate type: live bryophytes, dead bryophytes, or bark (there is an additional wood/inner bark substrate type used in the spatial study but not in the temporal one)

year - 
Collection year:  Y15 for 2015, Y16 for 2016, Y17 for 2017

y.num - 
Collection year in a purely numeric format

branch - 
Tree branch that the sample was collected from:  B01-B05

point - 
Collection point along the branch, with each branch having up to 27 points, going from P01 closest to the trunk and P27 closest to the branch tip

BP - 
Branch + point, the full location that a sample was collected from; so all substrates dissected from a single core will have the same value, as will samples from the same location but collected in different years)

BPS - 
Branch + point + substrate

YBP - 
Year + branch + point; all substrates collected from the same location WITHIN the same year will have the same value

YS - 
Year + substrates; all substrates of the same kind collected in the same year will have the same value, regardless of the location they were collected from

old_name - 
Older names for samples, retained for use with older versions of the data no longer in use; this can be ignored

