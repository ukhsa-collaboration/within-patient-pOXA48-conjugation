library(ape)
library(tidyverse)


####################### adding node labels to the tree to submit to augur ancestral ###################################
masked_tree_Pat<-read.tree("Final_alignment_masked_PAT_50bases_0.1samples.fst.treefile")
masked_tree_Pat$node.label<-paste("node",1:masked_tree_Pat$Nnode, sep="-")
write.tree(masked_tree_Pat, "Final_NODES_alignment_masked_PAT_50bases_0.1samples.fst.treefile" )
#####################################################################################################



deduplicated_masked_tree_Pat<-read.tree("Final_NODES_alignment_masked_PAT_50bases_0.1samples.fst.treefile")

sttab<-as_tibble(read.table(file="data/supplementary_data.csv", sep=",",  header=TRUE))  

#####################################################################################
# FUNCTION NEEDED
######################################################################################
# 3 regions have been cut out from plasmid allignemnts because sequences did not allign well there
# to plot back on the plasmid I need to add those regios back to the mutations I find
########################################################################################
fastafile<- read.dna(file="", format="fasta", as.character=TRUE)

referenceAll<-paste0(fastafile["JN626286.1",],collapse="")
fastafile<-0
gc()

tobecut<-read.csv("data/table_masked.csv")
  
backtopos<-function(loc){
  oldstring<-substring(referenceAll, 1, loc)
  noGapString<-str_remove_all(oldstring, "-")
  oldloc<-length(unlist(str_split(noGapString,"")))
  
  if(loc<tobecut[1,3]){newloc<-oldloc}
  else{
    if (loc<tobecut[2,3]){newloc<-oldloc+1500}
    else{
      if(loc<tobecut[3,3]){newloc<-oldloc+3500}
      else{
        newloc<-oldloc+4330
      }
    }
  }
  return(newloc)
}

######################################################################################
# process the augur mutations so that I can plot them
####################################################################################
################# Clades ##############################################################
#######################################################################################
#CladeA
CladeAsamples<- c("PAT45-S_100", "PAT45-S_101", "PAT45-S_102", "PAT240-S_519" ,"PAT240-S_520")
#CladeB
CladeBsamples<-c("PAT63-S_139", "PAT63-S_140","PAT229-S_495","PAT74-S_164")
#CladeC
CladeCsamples<-c("PAT157-S_341", "PAT157-S_340" ,"PAT157-S_342","PAT64-S_141", "PAT64-S_142",
                 "PAT48-S_107" ,"PAT48-S_108","PAT22-S_052", "PAT22-S_053","PAT235-S_509")
#cladeD
CladeDsamples<-c("PAT211-S_457","PAT111-S_243" ,"PAT111-S_244","PAT175-S_378", "PAT175-S_379",
                 "PAT224-S_485", "PAT224-S_486", "PAT70-S_155", "PAT70-S_156","PAT33-S_075", "PAT33-S_074",
                 "PAT148-S_322" ,"PAT148-S_323")
#CladeE
CladeEsamples<-c("PAT141-S_308", "PAT141-S_309","PAT219-S_474", "PAT219-S_475", "PAT93-S_205", "PAT93-S_206",
                 "PAT76-S_167" ,"PAT76-S_168")
#CladeF
CladeFsamples<-c("PAT09-S_019","PAT09-S_021","PAT18-S_043" , "PAT18-S_045", "PAT44-S_098", "PAT44-S_099")

# clade A
cladeA<-getMRCA(deduplicated_masked_tree_Pat,CladeAsamples)
pippoclade<-extract.clade(deduplicated_masked_tree_Pat, cladeA)
plot.phylo(pippoclade,show.node.label = TRUE)
AllMutCladeA<- allmutations["node-36"]
OnlyMutCladeA<-onlymutations["node-36"]
MutCladeA<-as.numeric(
  unlist(
    str_split(
      str_replace_all(
        str_replace_all(OnlyMutCladeA,pattern = "[A-Z]",replacement=""),pattern = "-",replacement=""),",")))



#clade B
cladeB<-getMRCA(deduplicated_masked_tree_Pat, CladeBsamples)
pippoclade<-extract.clade(deduplicated_masked_tree_Pat, cladeB)
plot.phylo(pippoclade,show.node.label = TRUE)
AllMutCladeB<- allmutations["node-147"]
OnlyMutCladeB<-onlymutations["node-147"]
MutCladeB<-as.numeric(
  unlist(
    str_split(
      str_replace_all(
        str_replace_all(OnlyMutCladeB,pattern = "[A-Z]",replacement=""),pattern = "-",replacement=""),",")))


#clade C
cladeC<-getMRCA(deduplicated_masked_tree_Pat, CladeCsamples)
pippoclade<-extract.clade(deduplicated_masked_tree_Pat, cladeC)
plot.phylo(pippoclade,show.node.label = TRUE)
AllMutCladeC<- allmutations["node-131"]
OnlyMutCladeC<-onlymutations["node-131"]
MutCladeC<-as.numeric(
  unlist(
    str_split(
      str_replace_all(
        str_replace_all(OnlyMutCladeC,pattern = "[A-Z]",replacement=""),pattern = "-",replacement=""),",")))


#clade D
cladeD<-getMRCA(deduplicated_masked_tree_Pat,  CladeDsamples)
pippoclade<-extract.clade(deduplicated_masked_tree_Pat, cladeD)
plot.phylo(pippoclade,show.node.label = TRUE)
AllMutCladeD<- allmutations["node-153"]
OnlyMutCladeD<-onlymutations["node-153"]
MutCladeD<-as.numeric(
  unlist(
    str_split(
      str_replace_all(
        str_replace_all(OnlyMutCladeD,pattern = "[A-Z]",replacement=""),pattern = "-",replacement=""),",")))

#clade E
cladeE<-getMRCA(deduplicated_masked_tree_Pat,  CladeEsamples)
pippoclade<-extract.clade(deduplicated_masked_tree_Pat, cladeE)
plot.phylo(pippoclade,show.node.label = TRUE)
AllMutCladeE<- allmutations["node-205"]
OnlyMutCladeE<-onlymutations["node-205"]
MutCladeE<-as.numeric(
  unlist(
    str_split(
      str_replace_all(
        str_replace_all(OnlyMutCladeE,pattern = "[A-Z]",replacement=""),pattern = "-",replacement=""),",")))

#clade F
cladeF<-getMRCA(deduplicated_masked_tree_Pat,  CladeFsamples)
pippoclade<-extract.clade(deduplicated_masked_tree_Pat, cladeF)
plot.phylo(pippoclade,show.node.label = TRUE)
AllMutCladeF<- allmutations["node-200"]
OnlyMutCladeF<-onlymutations["node-200"]
MutCladeF<-as.numeric(
  unlist(
    str_split(
      str_replace_all(
        str_replace_all(AllMutCladeF,pattern = "[A-Z]",replacement=""),pattern = "-",replacement=""),",")))


################ FIGURE 3b ####################################


############ the proper figure ###################

rectEnds<-read.table("data/AnonympOXA48GenesLocSmall.csv", sep=";",header=TRUE)
#### this is needed for the shade in the red points in the bottom part of the figure
library(scales)
colorshadekhaki<-alpha("khaki",alpha=0.7)
colorshadered<-alpha("firebrick",alpha=0.8)
axiscolore<-c("#582C83", "#00AB8E", "#E40046", "#00A5DF", "#FFB81C", "#003B5C")
############# This is the main part of the plot
par(fig=c(0,1,0,0.9))


plot(as.numeric(MutCladeA),rep(10.5,length(MutCladeA)),xlab="plasmid length",ylab="",pch="|", 
     ylim=c(0.2,11),xlim=c(0,62000),yaxt="n", cex=1.5)

for (i in 1:length(rectEnds[,2])){
  rect(rectEnds[i,2],0.2,rectEnds[i,3],11,col="grey60",border=NA,)
}

for (j in 1:length(tobecut[,2])){
  rect(tobecut[j,3],0.2,tobecut[j,4],11, col=colorshadekhaki, border=NA,)
}

rect(rectEnds[7,2],0.2,rectEnds[7,3],11,col=colorshadered,border=NA,)

points(as.numeric(MutCladeA),rep(10.5,length(MutCladeA)),pch="|", cex=1.5, col="#582C83")
points(as.numeric(MutCladeB),rep(8.5,length(MutCladeB)),pch="|", cex=1.5, col="#00AB8E")
points(as.numeric(MutCladeC),rep(6.5,length(MutCladeC)),pch="|", cex=1.5, col="#E40046")
points(as.numeric(MutCladeD),rep(4.5,length(MutCladeD)),pch="|", cex=1.5, col="#0036df")
points(as.numeric(MutCladeE),rep(2.5,length(MutCladeE)),pch="|", cex=1.5, col="#FFB81C")
points(as.numeric(MutCladeF),rep(0.7,length(MutCladeF)),pch="|", cex=1.5, col= "#003B5C")


axis(side=2,las=2, at=c(0.5,2.5,4.5,6.5,8.5,10.5),
     label=c("clade F","clade E","clade D", "clade C", "clade B", "clade A"),cex.axis=0.8,
     col=c("#582C83", "#00AB8E", "#E40046", "#00A5DF", "#FFB81C", "#003B5C"), font=2)

abline(h=1.5,col="grey40")
abline(h=3.5,col="grey40")
abline(h=5.5,col="grey40")
abline(h=7.5,col="grey40")
abline(h=9.5,col="grey40")

########### figure 3a #############
#########preparing for the plot ##################
env<-read.table("Deduplicate/AnonympOXA48GenesLocSmall.csv", sep=";",header=TRUE)
namesP <- env[,1]
startsP <- env[,2]
endsP <-env[,3]
strandsP <-env[,4]
dfP <- data.frame(name=namesP, start=startsP, end=endsP,strand=strandsP)

library(genoPlotR)

dna_segP <- dna_seg(dfP)
dna_segP[,"fill"]<-"grey80"
dna_segP[7,"fill"]<-"red"
dna_segP[4,"fill"]<-"grey28"
dna_segP[c(23:35),"fill"]<-"royalblue"

dna_segP[26,"fill"]<-"grey48"
dna_segP[c(36:38),"fill"]<-"orchid"
dna_segP[c(1,2,3,39,40),"fill"]<-"gold"
dna_segP[c(21,22),"fill"]<-"chartreuse"

dna_segP[,"col"]<-"black"

ENVP<-list(dna_segP)

#mid_pos <-ENV[[1]]$start+((ENV[[1]]$start-middle(ENV[[1]]))/2)
mid_pos<-middle(ENVP[[1]])

custom_pos<-mid_pos


annotUNO <- annotation(x1 = mid_pos, x2 = rep(NA, length(mid_pos)), text =ENVP[[1]]$name, rot = rep(80, length(mid_pos)), col = rep("black", length(mid_pos)))
annotDUE <- annotation(x1 = custom_pos, x2 = rep(NA, length(mid_pos)), text =ENVP[[1]]$name, rot = rep(80, length(mid_pos)), col = rep("black", length(mid_pos)))



par(fig=c(0.2,1,0.8,1))

plot_gene_map(ENVP,annotations=annotUNO,plot_new=TRUE, cex=0.8)
plot_gene_map(ENVP,annotations=annotDUE,plot_new=TRUE, cex=0.8)

custom_pos<-mid_pos
custom_pos<-c( 61881.0,   200.0 , 1502.5  ,2900.0 , 3800.0  ,4685.0 , 5843.5  ,7000.0,
               7900.0,  8815.0,  9633.0 ,17200.0 ,18500.0, 19400.0, 20391.5 ,23716.0,
               25891.5, 27655.0, 29587.0 ,30800.0, 33353.5 ,34490.0, 35500.0, 36378.5,
               37360.5, 39862.5, 41200.0, 42500.0, 43946.5, 45208.5, 46247.0, 47200.0,
               48200.0, 49621.0, 51747.0, 55000.0, 56000.0, 56921.5, 59000.0, 60500.0)


annotDUE <- annotation(x1 = custom_pos, x2 = rep(NA, length(mid_pos)), text =ENVP[[1]]$name, rot = rep(80, length(mid_pos)), col = rep("black", length(mid_pos)))
plot_gene_map(ENVP,annotations=annotDUE,plot_new=TRUE, cex=0.8)

dev.copy2pdf(file="plasmidAnnotNonOverl.pdf")

###############################################################################################
#################################################################################################



par(fig=c(0,1,0,0.9))


plot(as.numeric(MutCladeA),rep(10.5,length(MutCladeA)),xlab="plasmid length",ylab="",pch="|", 
     ylim=c(0.2,11),xlim=c(0,62000),yaxt="n", cex=1.5)

for (i in 1:length(rectEnds[,2])){
  rect(rectEnds[i,2],0.2,rectEnds[i,3],11,col="grey60",border=NA,)
}

for (j in 1:length(tobecut[,2])){
  rect(tobecut[j,3],0.2,tobecut[j,4],11, col=colorshadekhaki, border=NA,)
}

rect(rectEnds[7,2],0.2,rectEnds[7,3],11,col=colorshadered,border=NA,)

points(as.numeric(MutCladeA),rep(10.5,length(MutCladeA)),pch="|", cex=1.5, col="#582C83")
points(as.numeric(MutCladeB),rep(8.5,length(MutCladeB)),pch="|", cex=1.5, col="#00AB8E")
points(as.numeric(MutCladeC),rep(6.5,length(MutCladeC)),pch="|", cex=1.5, col="#E40046")
points(as.numeric(MutCladeD),rep(4.5,length(MutCladeD)),pch="|", cex=1.5, col="#0036df")
points(as.numeric(MutCladeE),rep(2.5,length(MutCladeE)),pch="|", cex=1.5, col="#FFB81C")
points(as.numeric(MutCladeF),rep(0.7,length(MutCladeF)),pch="|", cex=1.5, col= "#003B5C")


axis(side=2,las=2, at=c(0.5,2.5,4.5,6.5,8.5,10.5),
     label=c("clade F","clade E","clade D", "clade C", "clade B", "clade A"),cex.axis=0.8,
     col=c("#582C83", "#00AB8E", "#E40046", "#00A5DF", "#FFB81C", "#003B5C"), font=2)

abline(h=1.5,col="grey40")
abline(h=3.5,col="grey40")
abline(h=5.5,col="grey40")
abline(h=7.5,col="grey40")
abline(h=9.5,col="grey40")


par(fig=c(0,1,0.8,1))

plot_gene_map(ENVP,annotations=annotUNO, cex=0.8)





###############################################################################################
# delete samplesn in non-outbreak clade
##############################################################################################

nonOutbreakClade<-extract.clade(deduplicated_masked_tree_Pat, 451)
tipsNonOutbreak<-nonOutbreakClade$tip.label

alignment<- read.dna(file="deduplicated_masked_aligned_Pat.fst", format="fasta")
outbreak_tips<- labels(alignment) [!(labels(alignment) %in% tipsNonOutbreak)]
outbreak<-alignment[outbreak_tips,]
outbreak_short<-outbreak[,1:59712]
outbreak<-outbreak_short[,16:59712]
write.dna(outbreak, "outbreakSeq.fst", format="fasta")

############################################################################################
#isolate samples in non-outbreak clade
#####################################################################################

nonOutbreakClade<-extract.clade(deduplicated_masked_tree_Pat, 451)
tipsNonOutbreak<-nonOutbreakClade$tip.label
