library(tidyverse)
library(ggplot2)
library(ape)
library(igraph)
library(lubridate)
library(BiocManager)
library(Biostrings)
library(tidyverse)
library(rcartocolor)
library(RColorBrewer)
library(latex2exp)
library("dplyr")


source("scripts/DistDate.R")







# UKHSA colour scheme:
UKHSA_cols <- list(
  UKHSA.teal = "#007C91",
  midnight   = "#003B5C",
  plum       = "#582C83",
  moonlight  = "#1D57A5",
  wine       = "#8A1B61",
  cherry     = "#E40046",
  DHSC.green = "#00AB8E",
  ocean      = "#00A5DF",
  grass      = "#84BD00",
  tangerine  = "#FF7F32",
  sunny      = "#FFB81C",
  sand       = "#D5CB9F"
)

UKHSA_cols_ordered = c(UKHSA_cols$plum, UKHSA_cols$grass, UKHSA_cols$tangerine, UKHSA_cols$cherry, UKHSA_cols$moonlight, UKHSA_cols$sunny, UKHSA_cols$wine, UKHSA_cols$DHSC.green, UKHSA_cols$sand, UKHSA_cols$midnight)

custom_ukhsa_theme <- list(legend.position     ="none",
                           panel.grid.major.y  = element_line(color="grey"),
                           # AXIS
                           axis.text.x         = element_text(size=12, colour="black", margin=margin(t=-10)),
                           axis.text.y         = element_text(size=12, colour="black", hjust=0.5, margin=margin(r=10)),
                           axis.ticks.x.bottom = element_line(colour="grey"),
                           axis.ticks.y.left   = element_line(colour="grey"),
                           # axis.ticks.length.y = unit(c(0.5), "cm"),
                           axis.title.y        = element_text(face="bold", angle=90, size=12, colour="black"),
                           axis.title.x        = element_text(face="bold", size=12, colour="black"),
                           # PLOT
                           # plot.background = element_rect(colour = "black", fill=NA, size=1), # - border colour
                           plot.margin         = unit(c(1,1,1,1), "cm"),
                           plot.title          = element_text(face="bold", size=20, colour="black"),
                           plot.subtitle       = element_text(size=12, colour="black", margin=margin(b=10))
)

default_figure_width <- 26
default_figure_height <- 16










#####################################################################################
# READ DATA: EPI DATA
#####################################################################################
# Read in epi file (contains data on region, species, time, ...)


epi_table<-as_tibble(read.table(file="./data/finalEpiDataRef.csv", sep=",",  header=TRUE))

# Removing 3 Patients (see methods.rmd)
patients_to_remove <- c("PAT182", "PAT202", "PAT230")
epi_table <- epi_table %>% filter(! ORDPATNAME %in% patients_to_remove)


# computing the time distances between samples (in days)
epi_table<-epi_table %>%
  mutate(timeDist=unlist(pmap(list(day,month,year),distanceDate))) %>% 
  mutate(dateS=unlist(pmap(list(day,month,year),function(x,y,z){as.Date(str_c(c(z,y,x),collapse="-"))}))) %>%
  mutate(timeDec=unlist(pmap(list(day,month,year),function(x,y,z){decimal_date(as.Date(str_c(c(z,y,x),collapse="-")))}))) %>%
  mutate(PatSample=str_c(ORDPATNAME, SAMPLE_ID,sep="-")) %>% mutate(PatRes=str_c(ORDPATNAME, CARB,sep="-")) %>%
  mutate(PatSp=str_c(ORDPATNAME, Species,sep="-"))
# Removing timedist outlier: (<0 days - this was because of an error - date of birth was keyed in as a sample date)
epi_table<-epi_table %>% mutate(timeDist = ifelse(timeDist < 0, NA, timeDist))

#############################################################
OXApatients<- epi_table %>%  select (ORDPATNAME) %>% unique ()

timedOXA<- epi_table %>% filter (ORDPATNAME %in% as.vector(OXApatients$ORDPATNAME) )









#####################################################################################
# READ DATA: ALIGNED SEQUENCE DATA + CALCULATE GENETIC DISTANCES
#####################################################################################
cleaned_and_aligned_sequence_data<-read.dna("./data/deduplicated_id_seq_clean_aligned.fasta",format="fasta")
# Calculating SNPs
genetic_distance_matrix<-dist.dna(cleaned_and_aligned_sequence_data, model="N", as.matrix=TRUE)
# REMOVING LOWER TRIANGULAR (as genetic distance of a-b is the same as b-a; so we would have twice as many data points, needlessly) + DIAGONAL (genetic distance of a-a is 0, this will skew the distributions)
# genetic_distance_table <- genetic_distance_matrix[upper.tri(genetic_distance_matrix, diag = FALSE)]
gdt_ind <- which( upper.tri(genetic_distance_matrix,diag=FALSE) , arr.ind = TRUE )
genetic_distance_table <- data.frame( pl1 = dimnames(genetic_distance_matrix)[[2]][gdt_ind[,2]] ,
                                      pl2 = dimnames(genetic_distance_matrix)[[1]][gdt_ind[,1]] ,
                                      SNPs = genetic_distance_matrix[ gdt_ind ] )

# Converting data type to numeric
genetic_distance_table<-genetic_distance_table %>% 
  mutate(SNPs=as.numeric(SNPs))
# Joining patient id information:
genetic_distance_table <- left_join(genetic_distance_table, epi_table %>% select(MOLIS, ORDPATNAME) %>% mutate(ORDPATNAME=ifelse(is.na(ORDPATNAME), "Reference", ORDPATNAME)) %>% dplyr::rename(pl1=MOLIS, pl1_pat=ORDPATNAME), by = "pl1")
genetic_distance_table <- left_join(genetic_distance_table, epi_table %>% select(MOLIS, ORDPATNAME) %>% mutate(ORDPATNAME=ifelse(is.na(ORDPATNAME), "Reference", ORDPATNAME)) %>% dplyr::rename(pl2=MOLIS, pl2_pat=ORDPATNAME), by = "pl2")
# Adding between/within patient data:
genetic_distance_table$within_patient <- (genetic_distance_table$pl1_pat == genetic_distance_table$pl2_pat)
# genetic_distance_table[genetic_distance_table$pl1_pat == "Reference",]$within_patient <- NA
genetic_distance_table[genetic_distance_table$pl2_pat == "Reference",]$within_patient <- NA

write.csv(genetic_distance_table %>% select(c("pl1_pat", "pl2_pat", "SNPs", "within_patient")),"data/genetic_distance_table_removed_molis.csv", row.names = FALSE)
write.csv(genetic_distance_table,"data/genetic_distance_table.csv", row.names = FALSE)

















### Pairwise Single Nucleotide Polymorphisms Distribution GRAPHS ###

##########################################################################################
# PLOT: SNP Histogram (all sample-sample pairs) - CORRECTED
##########################################################################################

genetic_distance_table <- read.csv("data/genetic_distance_table.csv")
genetic_distance_table$within_patient <- as.logical(genetic_distance_table$within_patient)
genetic_distance_table <- genetic_distance_table %>%
  filter(pl1 != pl2) %>% # remove self-comparison
  filter(pl1 != "JN626286") %>% # remove any reference plasmid
  filter(pl2 != "JN626286") %>% # remove any reference plasmid
  filter(pl1 != "Reference") %>% # remove any reference plasmid
  filter(pl2 != "Reference") # remove any reference plasmid

ggplot(genetic_distance_table , aes(x=SNPs))+
  geom_histogram(col="white", bins=101)+
  # geom_histogram(col="white", width=1,breaks=c(-0.5, seq(1,10,1) - 0.5, seq(20,100,10) - 5, seq(200, 1000, 100) - 50, seq(2000, 10000, 1000) - 500, seq(20000,100000, 10000) - 5000))+
  # geom_histogram(col="white", breaks=c(seq(0,10,1), seq(15,95,5), seq(100,700, 50)))+
  # scale_x_log10() +
  # scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0,1,10,100,1000), minor_breaks=c(seq(2,9,1), seq(20,90,10), seq(200,900,100)), limits=c(NA, 1000)) +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0,1,10,100,700), minor_breaks=c(seq(2,9,1), seq(20,90,10), seq(200,700,100))) +
  ylim(0,12000) +
  # scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits=c(NA,30000), breaks=c(0,1,10,100,1000, 10000, 30000), minor_breaks=c(seq(2,9,1), seq(20,90,10), seq(200,900,100), seq(2000,9000,1000), seq(20000,30000,10000))) +
  theme_void() +
  do.call(theme, custom_ukhsa_theme) + # theme with custom_ukhsa_theme passed as arguments
  ggtitle("Distribution of SNPs", subtitle="\nAll Pairs of Samples (excluding Reference)") +
  xlab("SNPs") +
  ylab("Count\n") +
  geom_hline(yintercept=0, color="grey")+
  theme(
    panel.grid.minor.y  = element_line(color="gray75", linetype = "dotted"),
  )

ggsave("figures/snp_distribution_all_pairs.pdf", width=default_figure_width, height=default_figure_height, units = "cm")


##########################################################################################
# PLOT: SNP Histogram (all between-patient pairs)
##########################################################################################


ggplot(genetic_distance_table %>% filter(within_patient == FALSE) , aes(x=SNPs))+
  geom_histogram(col="white", bins=31)+
  # geom_histogram(col="white", width=1,breaks=c(-0.5, seq(1,10,1) - 0.5, seq(20,100,10) - 5, seq(200, 1000, 100) - 50, seq(2000, 10000, 1000) - 500, seq(20000,100000, 10000) - 5000))+
  # geom_histogram(col="white", breaks=c(seq(0,10,1), seq(15,95,5), seq(100,700, 50)))+
  # scale_x_log10() +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0,1,10,100,1000), minor_breaks=c(seq(2,9,1), seq(20,90,10), seq(200,900,100)), limits=c(NA, 1000)) +
  # scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0,1,10,100,700), minor_breaks=c(seq(2,9,1), seq(20,90,10), seq(200,700,100)), limits=c(NA, 700)) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0,1,10,100,500), limits=c(NA,500), minor_breaks=c(seq(2,9,1), seq(20,90,10), seq(200,400,100))) +
  theme_void() +
  do.call(theme, custom_ukhsa_theme) + # theme with custom_ukhsa_theme passed as arguments
  ggtitle("Distribution of SNPs", subtitle="\nPairs of Samples Between Patients (excluding Reference)") +
  xlab("SNPs") +
  ylab("Count\n") +
  geom_hline(yintercept=0, color="grey")+
  theme(panel.grid.minor.y  = element_line(color="gray75", linetype = "dotted"))

ggsave("figures/snp_distribution_within_patient_pairs.pdf", width=default_figure_width, height=default_figure_height, units = "cm")



##########################################################################################
# PLOT: SNP Histogram (all within-patient pairs)
##########################################################################################


ggplot(genetic_distance_table %>% filter(within_patient == TRUE) , aes(x=SNPs))+
  geom_histogram(col="white", bins=31)+
  # geom_histogram(col="white", width=1,breaks=c(-0.5, seq(1,10,1) - 0.5, seq(20,100,10) - 5, seq(200, 1000, 100) - 50, seq(2000, 10000, 1000) - 500, seq(20000,100000, 10000) - 5000))+
  # geom_histogram(col="white", breaks=c(seq(0,10,1), seq(15,95,5), seq(100,700, 50)))+
  # scale_x_log10() +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0,1,10,100,1000), minor_breaks=c(seq(2,9,1), seq(20,90,10), seq(200,900,100)), limits=c(NA, 1000)) +
  # scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0,1,10,100,700), minor_breaks=c(seq(2,9,1), seq(20,90,10), seq(200,700,100)), limits=c(NA, 700)) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0,1,10,100,500), limits=c(NA,500), minor_breaks=c(seq(2,9,1), seq(20,90,10), seq(200,400,100))) +
  theme_void() +
  do.call(theme, custom_ukhsa_theme) + # theme with custom_ukhsa_theme passed as arguments
  ggtitle("Distribution of SNPs", subtitle="\nPairs of Samples Within Patients (excluding Reference)") +
  xlab("SNPs") +
  ylab("Count\n") +
  geom_hline(yintercept=0, color="grey")+
  theme(panel.grid.minor.y  = element_line(color="gray75", linetype = "dotted"))

ggsave("figures/snp_distribution_within_patient_pairs.pdf", width=default_figure_width, height=default_figure_height, units = "cm")

