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
#library(ggtext)

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


deduplicated_epi_table<-as_tibble(read.table(file="data/deduplicatedEpiDataRef.csv", sep=",",  header=TRUE))
epi_table<-as_tibble(read.table(file="data/finalEpiDataRef.csv", sep=",",  header=TRUE))

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
# Plot: Number of Samples per Patient Count Plot
#####################################################################################
epi_table_grouped_by_patient <- epi_table %>% 
  filter(epi_table$MOLIS != "JN626286") %>%
  group_by(ORDPATNAME) %>% 
  nest() %>% 
  mutate(n = map_dbl(.x=data,.f=~ length (.x$SAMPLE_ID)),
         meanTime=map_dbl(.x=data,.f=~ mean( .x$timeDec)),
         timediff=map_dbl(.x = data,.f=~ diff(c(min(.x$timeDist), max(.x$timeDist)))   ) )


plot_1 <- ggplot(epi_table_grouped_by_patient, aes(x=n)) +
  # geom_bar(fill="gray50") +
  geom_histogram(col="white", size=2, binwidth=1, fill="gray50") +
  theme_void() +
  do.call(theme, custom_ukhsa_theme) +
  ggtitle("1A. Number of Samples per Patient") +
  xlab("Number of Samples per Patient") +
  ylab("Count") +
  ylim(0,100) +
  geom_hline(yintercept=0, color="grey")

ggsave("figures/number_of_samples_per_patient_histogram.pdf", width=default_figure_width, height=default_figure_height, units = "cm")








#####################################################################################
# Plot: Number of Samples per Region Countplot
###################################################################################

epi_table_grouped_by_region <- epi_table  %>%
  filter(epi_table$MOLIS != "JN626286") %>%
  group_by(ORDPATNAME,Region) %>%
  nest()
epi_table_grouped_by_region$Region.Category = factor(epi_table_grouped_by_region$Region, levels=c("LONDON", "YORK&HUM", "N EAST", "N WEST", "E MIDS", "W MIDS", "EAST", "S EAST", "S WEST", "Turkey" ))
epi_table_grouped_by_region <- epi_table_grouped_by_region[order(epi_table_grouped_by_region$Region.Category),]

plot_2 <- ggplot(epi_table_grouped_by_region, aes(x=Region.Category)) + 
  geom_bar(stat="count", fill="gray50") +

  theme_void() +
  do.call(theme, custom_ukhsa_theme) + # theme with custom_ukhsa_theme passed as arguments
  ggtitle("1B. Region of Patient") +
  xlab("\nRegion") +
  ylab("Count\n") +
  ylim(0,80)

ggsave("figures/sample_region_countplot.pdf", width=default_figure_width, height=default_figure_height, units = "cm")




#####################################################################################
# Plot: Time in between samples Histogram
#####################################################################################
plot_3 <- ggplot(epi_table_grouped_by_patient, aes(x=timediff)) +
  geom_histogram(col="white", size=1, fill="gray50") +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0,1,10,100,500), limits=c(NA,500)) +
  theme_void() + 
  do.call(theme, custom_ukhsa_theme) +
  ggtitle("1C. Time Between Patient's First and Last Sample") +
  xlab("Time (days)") +
  ylab("Count\n") +
  geom_hline(yintercept=0, color="grey")

ggsave("figures/sample_time_difference_histogram.pdf", width=default_figure_width, height=default_figure_height, units = "cm")



##########################################################################################
# PLOT: Sample Date Histogram
##########################################################################################


timedOXAallMY <- epi_table %>% filter(epi_table$MOLIS != "JN626286")
timedOXAallMY <- timedOXAallMY %>% mutate(sample.date=as.Date(with(timedOXAallMY, paste(year,month,day,sep="-")),"%Y-%m-%d"))

plot_4 <- ggplot(timedOXAallMY %>% drop_na(sample.date) %>% filter(sample.date != "1947-11-23 "), aes(x=sample.date)) +
  geom_histogram(bins=36, col="white", size=1, fill="gray50") +
  theme_void() + 
  do.call(theme, custom_ukhsa_theme) +
  ggtitle("1D. Date Sample Taken") +
  xlab("Date") +
  ylab("Count\n") +
  geom_hline(yintercept=0, color="grey") +
  xlim(as.Date("2014-1-1"), as.Date("2016-12-31"))

ggsave("figures/sample_date_histogram.pdf", width=default_figure_width, height=default_figure_height, units = "cm")







##########################################################################################
# PLOT: Species Countplot
##########################################################################################

# Adding second column of information (bacterial species) - grouping the species together...
species_count_table <- epi_table %>% 
  filter(epi_table$MOLIS != "JN626286") %>% 
  group_by(Species) %>% 
  count() %>% 
  arrange(desc(n))
species_count_table$SpeciesGroup <- species_count_table$Species
species_count_table <- species_count_table %>%
  mutate(SpeciesGroup=map_chr(SpeciesGroup, function(x){ ifelse(x %in% c("Escherichia_coli","Escherichia sp."), "Escherichia coli",x)})) %>%
  mutate(SpeciesGroup=map_chr(SpeciesGroup, function(x){ ifelse(x %in% c("Citrobacter amalonaticus","Citrobacter sp.","Citrobacter_freundii"),
                                                      "Citrobacter",x)})) %>%
  mutate(SpeciesGroup=map_chr(SpeciesGroup, function(x){ ifelse(x %in% c("Enterobacter_aerogenes","Enterobacter_cloacae" ,"Enterobacter_hormaechei" ),
                                                      "Enterobacter",x)})) %>%
  mutate(SpeciesGroup=map_chr(SpeciesGroup, function(x){ ifelse(x %in% c("Klebsiella_oxytoca","Klebsiella_variicola"),
                                                      "Other Klebsiella",x)})) %>%  
  mutate(SpeciesGroup=map_chr(SpeciesGroup, function(x){ ifelse(x %in% c("Morganella_morganii ","Raoultella planticola","Serratia marcescens" ),
                                                      "Other Enterobacterales",x)})) %>%
  mutate(SpeciesGroup=map_chr(SpeciesGroup, function(x){ ifelse(is.na(x),
                                                      "Other Enterobacterales",x)})) %>%
  mutate(SpeciesGroup=map_chr(SpeciesGroup, function(x){ ifelse(x %in% c("Proteus_mirabilis", "Providencia stuartii"), # <- FAN ADDITION, DOUBLE CHECK
                                                                "Other Enterobacterales",x)}))
species_count_table$Species <- species_count_table$Species %>% str_replace_all(c('_' = ' '))




plot_5 <- ggplot(species_count_table, aes(x=n, y=reorder(Species,n,decreasing=TRUE), fill=SpeciesGroup)) +
  geom_bar(stat="identity",color="white") +
  scale_fill_manual(
    name = "Bacterial Species Group",
    values = c("darkred", "red3", "lightcoral", "darkorange", "gold", "wheat"),
    breaks = c("Klebsiella_pneumoniae", "Other Klebsiella", "Enterobacter", "Other Enterobacterales", "Escherichia coli", "Citrobacter"),
    labels = c("Klebsiella pneumoniae", "Klebsiella (Other)", "Enterobacter", "Enterobacterales (Other)", "Escherichia coli", "Citrobacter"),
    na.value = "white"
  ) +
  scale_x_discrete(breaks=c(0,20,40,60,80,100)) +
  theme_void() +
  do.call(theme, custom_ukhsa_theme) +
  # ggtitle("Distribution of Species") +
  ggtitle("1E. Sample Species") +
  xlab("Count") +
  ylab("Species") +
  geom_text(aes(x=n+1, label=Species), size=4, hjust=0) +
  theme(axis.text.y = NULL) +
  theme(axis.title.x = element_text(face="bold", size=12, margin=margin(t=10))) +
  theme(axis.text.x = element_text(size=12, margin=margin(t=10))) +
  theme(panel.grid.major.y = NULL, panel.grid.major.x = element_line(color="grey"), axis.ticks.length.y = NULL ) +
  theme(plot.title = element_text(face="bold", size=20, colour="black", margin=margin(b=10))) +
  xlim(0, 100) +
  theme(# panel.background = element_rect(fill='transparent'), # LEGEND
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.position = "right",
        legend.position = c(0.835, 0.835),
        legend.margin=margin(c(5,-100,5,5)),
        # legend.position = c(1.09,.7895),
        legend.title = element_text(size = 10,face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(c(.7,.7), "cm"),
        #legend.box = "horizontal",
        # legend.margin = margin(),
        # legend.margin
        legend.background = element_rect(fill='white')
  )


ggsave("figures/species_countplot.pdf", width=default_figure_width, height=default_figure_height, units = "cm")





##########################################################################################
# PLOT: Strain Countplot
##########################################################################################

# Adding second column of information (bacterial species) - grouping the species together...
strain_count_table <- epi_table %>% group_by(Strain, Species) %>% count() %>% arrange(desc(n))
strain_count_table$SpeciesGroup <- strain_count_table$Species
strain_count_table <- strain_count_table %>%
  mutate(SpeciesGroup=map_chr(SpeciesGroup, function(x){ ifelse(x %in% c("Escherichia_coli","Escherichia sp."), "Escherichia coli",x)})) %>%
  mutate(SpeciesGroup=map_chr(SpeciesGroup, function(x){ ifelse(x %in% c("Citrobacter amalonaticus","Citrobacter sp.","Citrobacter_freundii"),
                                                                "Citrobacter",x)})) %>%
  mutate(SpeciesGroup=map_chr(SpeciesGroup, function(x){ ifelse(x %in% c("Enterobacter_aerogenes","Enterobacter_cloacae" ,"Enterobacter_hormaechei" ),
                                                                "Enterobacter",x)})) %>%
  mutate(SpeciesGroup=map_chr(SpeciesGroup, function(x){ ifelse(x %in% c("Klebsiella_oxytoca","Klebsiella_variicola"),
                                                                "Other Klebsiella",x)})) %>%  
  mutate(SpeciesGroup=map_chr(SpeciesGroup, function(x){ ifelse(x %in% c("Morganella_morganii ","Raoultella planticola","Serratia marcescens" ),
                                                                "Other Enterobacterales",x)})) %>%
  mutate(SpeciesGroup=map_chr(SpeciesGroup, function(x){ ifelse(is.na(x),
                                                                "Other Enterobacterales",x)}))

strain_count_table <- strain_count_table %>% filter (Strain!="-") %>% mutate(SpeciesStrain=paste0(Species, "-",Strain)) %>% arrange(desc(n)) %>% ungroup() 
strain_count_table$SpeciesStrain <- strain_count_table$SpeciesStrain %>% str_replace_all(c('_' = ' '))


plot_6 <- ggplot(strain_count_table %>% filter(n>2), aes(x=n, y=reorder(SpeciesStrain,n,decreasing=TRUE), fill=SpeciesGroup))+
  geom_bar(stat="identity",color="white") +
  scale_fill_manual(
    name = "Bacterial Species Group",
    values = c("darkred", "red3", "lightcoral", "darkorange", "gold", "wheat"),
    breaks = c("Klebsiella_pneumoniae", "Other Klebsiella", "Enterobacter", "Other Enterobacterales", "Escherichia coli", "Citrobacter"),
    labels = c("Klebsiella pneumoniae", "Klebsiella (Other)", "Enterobacter", "Enterobacterales (Other)", "Escherichia coli", "Citrobacter"),
    na.value = "white"
  ) +
  theme_void()+
  do.call(theme, custom_ukhsa_theme) +
  ggtitle(label="1F. Sample Strain", subtitle="Only strains with at least 3 samples shown") +
  # ggtitle(label="Distribution of Strains", subtitle="Strains with at least 3 samples (not all strain data is present in the dataset)") +
  xlab("Count") +
  ylab("Strain") +
  geom_text(aes(x=n+0.2,label=SpeciesStrain), hjust=0, size=4)+
  theme(axis.text.y = NULL) +
  theme(panel.grid.major.y = NULL, panel.grid.major.x = element_line(color="grey")) +
  theme(plot.title = element_text(face="bold", size=20, colour="black", margin=margin(b=10))) +
  theme(axis.text.x = element_text(size=12, margin=margin(t=10))) +
  theme(axis.title.x = element_text(face="bold", size=12, margin=margin(t=10))) +
  xlim(0,20) +
  theme(
    legend.position = c(0.835, 0.835),
    legend.margin=margin(c(5,-100,5,5)),
    legend.title = element_text(size = 10,face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(c(.7,.7), "cm"),
    legend.background = element_rect(fill='white')
  )


ggsave("figures/strain_countplot.pdf", width=default_figure_width, height=default_figure_height, units = "cm")








######################################################################################
# table 1 on the number of non-OXA-48 resistances
#########################################################################################

View(epi_table %>% filter (CARB!="OXA-48") %>% group_by(CARB) %>% count() %>% arrange(n,desc=TRUE))





