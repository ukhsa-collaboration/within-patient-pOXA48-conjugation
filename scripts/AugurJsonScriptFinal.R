library(rjson)
library(tidyjson)
library(tidyverse)


filename<-"augur.output_dedu.data.json"
  
  file1<-fromJSON(file=filename)
  pippo<-file1$nodes
  nodenames<-attributes(pippo)$names
  tipsnames<-nodenames[1:262]
  
  allmutationstips<-c()

  onlymutationstips<-c()
  numbermutationstips<-c()
  for(i in 1:262){
    thistip<-tipsnames[i]
    itsmuts<-pippo[[i]]$muts
    allmutationstips[thistip]<-paste0(unlist(itsmuts), collapse=",")
    temp<-unlist(itsmuts)[grep("-", unlist(itsmuts), invert=TRUE ) ] 
    temp2<-temp[grep("N",temp, invert=TRUE ) ]
    onlymutationstips[thistip]<-paste0(temp2, collapse=",")
    numbermutationstips[thistip]<-length(temp2)
  }
  
  
  tibblemutations<-as_tibble(cbind(names(numbermutations),numbermutations)) %>% rename ("plasmid"="V1")
  
  allmutations<-c()
  
  onlymutations<-c()
  numbermutations<-c()
  for(i in 1:length(nodenames)){
    thistip<-nodenames[i]
    itsmuts<-pippo[[i]]$muts
    allmutations[thistip]<-paste0(unlist(itsmuts), collapse=",")
    temp<-unlist(itsmuts)[grep("-", unlist(itsmuts), invert=TRUE ) ] 
    temp2<-temp[grep("N",temp, invert=TRUE ) ]
    onlymutations[thistip]<-paste0(temp2, collapse=",")
    numbermutations[thistip]<-length(temp2)
  }
  
  #################################################################################
  # masked_tree_Pat is in TreePlot_Plasmids_dedu_Final.R
  ##################################################################################
  getMRCA(deduplicated_masked_tree_Pat, c("PAT22-S_052","PAT157-S_342" ),)
  minnie<-extract.clade(masked_tree_Pat, 424)
  plot.phylo(minnie)
  
  
  getMRCA(masked_tree_Pat, c("PAT175-S_379","PAT33-S_074" ))
  getMRCA(masked_tree_Pat, c("PAT240-S_519","PAT45-S_100" ))
  getMRCA(masked_tree_Pat, c("PAT76-S_167","PAT76-S_168" ))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  