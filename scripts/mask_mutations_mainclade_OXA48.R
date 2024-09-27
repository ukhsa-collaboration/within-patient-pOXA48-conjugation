library(ape)

source("scripts/DistDate.R")

###########################################################################################
# using ALL the samples in the outbreak clade
# # # # # # # # # # # # # # # # # # # # # # # 
# this is what we want to use as the mutations in the most distant samples look like bona fide mutation 
# (although probably from homologous recombination)
##################################################################################################

seq<-read.dna("outbreakSeqAllMafft.fst",format="fasta")

mdiff <- dist.dna(seq,model="N",as.matrix=T)

seqmainclade <- seq
maskedbases <- 50 # remove 50 bases from any gap in all directions
seqmaincladec <- as.character(seqmainclade)
seqmaincladec_masked <- seqmaincladec
for(i in 1:(dim(seqmaincladec)[1])){
  singleseqc<-seqmaincladec[i,]
  if(length(which(singleseqc=='-'))>0){ 
    seqmaincladec_masked[i,unique(pmin(dim(seqmaincladec)[2],pmax(0,as.vector(outer(which(singleseqc=='-'),-maskedbases:maskedbases,'+')))))]<-"n" 
  }
}
# main alignment masked (the reference IS included)
seqmainclade_masked <- as.DNAbin(seqmaincladec_masked)
#write.dna(seqmainclade_masked,"OutbreakCladeNoRefClean.fst", format="fasta")

# plot tree
#seq_final<-rbind(seq[c("JN626286.1","PAT18-S_045",  "PAT09-S_021", "PAT229-S_495",  
  #                     "PAT63-S_140",  "PAT63-S_139",  "PAT74-S_164", "PAT234-S_507"),],seqmainclade_masked)
seq_final<-seqmainclade_masked
#write.dna(seq_final,"OutbreakALLCladeWithREFClean.fst", format="fasta")

tree_final<-(bionjs(dist.dna(seq_final,as.matrix=T)))
tree_final<-root(tree_final,"JN626286.1")
plot(tree_final,cex=0.2)
# get frequency of the minor allele at each position
freqs<-apply(seqmaincladec_masked,2,function(x){y<-sort(table(x[x %in% c('a','c','g','t')]),decreasing=TRUE); if(length(y)>1){return(y[2])} else {return(0)}})
which(freqs>0)
which(freqs>2) # some of the sites in just 1-2 sequences are still a bit weird, though possibly correct

###########################################################################################
#removing the most distant samples on the outbreak clade
# not interesting
########################################################################################

seq<-read.dna("outbreakSeqAllMafft.fst",format="fasta")

mdiff <- dist.dna(seq,model="N",as.matrix=T)
mainclade <- which(mdiff[1,]<10) # by hand! this removes everything but the main clade
#these are the samples that have been removed
#c(PAT18-S_045,  PAT09-S_021, PAT229-S_495,  PAT63-S_140,  PAT63-S_139,  PAT74-S_164, PAT234-S_507, JN626286.1 )
seqmainclade <- seq[mainclade,]
#seqmainclade <- seq
maskedbases <- 50 # remove 50 bases from any gap in all directions
seqmaincladec <- as.character(seqmainclade)
seqmaincladec_masked <- seqmaincladec
for(i in 1:(dim(seqmaincladec)[1])){
  singleseqc<-seqmaincladec[i,]
  if(length(which(singleseqc=='-'))>0){ 
    seqmaincladec_masked[i,unique(pmin(dim(seqmaincladec)[2],pmax(0,as.vector(outer(which(singleseqc=='-'),-maskedbases:maskedbases,'+')))))]<-"n" 
  }
}
# main alignment masked
seqmainclade_masked <- as.DNAbin(seqmaincladec_masked)
#write.dna(seqmainclade_masked,"OutbreakCladeNoRefClean.fst", format="fasta")

# plot tree
seq_final<-rbind(seq[c("JN626286.1","PAT18-S_045",  "PAT09-S_021", "PAT229-S_495",  
                       "PAT63-S_140",  "PAT63-S_139",  "PAT74-S_164", "PAT234-S_507"),],seqmainclade_masked)
#seq_final<-seqmainclade_masked
#write.dna(seq_final,"OutbreakCladeWithREFClean.fst", format="fasta")

tree_final<-(bionjs(dist.dna(seq_final,as.matrix=T)))
tree_final<-root(tree_final,"JN626286.1")
plot(tree_final,cex=0.2)
# get frequency of the minor allele at each position
freqs<-apply(seqmaincladec_masked,2,function(x){y<-sort(table(x[x %in% c('a','c','g','t')]),decreasing=TRUE); if(length(y)>1){return(y[2])} else {return(0)}})
which(freqs>0)
which(freqs>2) # some of the sites in just 1-2 sequences are still a bit weird, though possibly correct



#############################################################################################################
# redoing the same for the alignment obtained using mafft --add to add the non outbreak samples 
# on the aligned and cleaned (with this same script) outbreak clade.
#############################################################################################################

seq<-read.dna("Re-aligned_allClade_allSEq_WithRef_Mafft.fst",format="fasta")

mdiff <- dist.dna(seq,model="N",as.matrix=T)
mainclade <- which(mdiff[1,]<10) # by hand! this removes everything but the main clade
seqmainclade <- seq[mainclade,]
maskedbases <- 50 # remove 50 bases from any gap in all directions
seqmaincladec <- as.character(seqmainclade)
seqmaincladec_masked <- seqmaincladec
for(i in 1:(dim(seqmaincladec)[1])){
  singleseqc<-seqmaincladec[i,]
  if(length(which(singleseqc=='-'))>0){ 
    seqmaincladec_masked[i,unique(pmin(dim(seqmaincladec)[2],pmax(0,as.vector(outer(which(singleseqc=='-'),-maskedbases:maskedbases,'+')))))]<-"n" 
  }
}
# main alignment masked
seqmainclade_masked <- as.DNAbin(seqmaincladec_masked)
#write.dna(seqmainclade_masked,"Re-aligned_all_seq_NOREF_Clean.fst", format="fasta")

# plot tree
seq_final<-rbind(seq["JN626286.1",],seqmainclade_masked)
#write.dna(seq_final, "Re-aligned_all_seq_NOREF_Clean.fst", format="fasta")

tree_final<-(bionjs(dist.dna(seq_final,as.matrix=T)))
tree_final<-root(tree_final,"JN626286.1")
plot(tree_final,cex=0.2)
# get frequency of the minor allele at each position
freqs<-apply(seqmaincladec_masked,2,function(x){y<-sort(table(x[x %in% c('a','c','g','t')]),decreasing=TRUE); if(length(y)>1){return(y[2])} else {return(0)}})
which(freqs>0)
which(freqs>2) # some of the sites in just 1-2 sequences are still a bit weird, though possibly correct

