library(ape)
library(tidyverse)


########### load data files #######################################################

deduplicated_masked_tree_Pat<-read.tree("deduplicated_masked_tree_Pat.tre")
sttab<-as_tibble(read.table(file="supplementary_data.csv", sep=",",  header=TRUE))
alignment<- read.dna(file="deduplicated_masked_aligned_Pat.fst", format="fasta")
outbreak_realigned<-read.dna(file="outbreakSeqAllMafft.fst", format="fasta")

##############################################################################################
# delete samples in non-outbreak clade
##############################################################################################

nonOutbreakClade<-extract.clade(deduplicated_masked_tree_Pat, 451)
tipsNonOutbreak<-nonOutbreakClade$tip.label


outbreak_tips<- labels(alignment) [!(labels(alignment) %in% tipsNonOutbreak)]
outbreak<-alignment[outbreak_tips,]
outbreak_short<-outbreak[,1:59712]
outbreak<-outbreak_short[,16:59712]
write.dna(outbreak, "outbreakSeq.fst", format="fasta")

####################################################################################
# these samples have been realigned (with reference) using MAFFT (standard parameters)
# mafft outbreakSeq.fst > outbreakSeqAllMafft.fst
# the alignment has then been masked for any mutation that happen within 50 bp from a gap
# using the script mask_mutations_mainclade_OXA48.R
# now I need to re-add the sequences I have excluded above and reallign using mafft


#####################################################################################
#isolate samples in non-outbreak clade
#####################################################################################

nonOutbreakClade<-extract.clade(deduplicated_masked_tree_Pat, 451)
tipsNonOutbreak<-nonOutbreakClade$tip.label

non_outbreak_tips<- labels(alignment) [(labels(alignment) %in% tipsNonOutbreak)]
non_outbreak<-alignment[non_outbreak_tips,]
non_outbreak_short<-non_outbreak[,1:59712]
non_outbreak<-non_outbreak_short[,16:59712]
write.dna(non_outbreak, "NON_outbreakSeq.fst", format="fasta")
