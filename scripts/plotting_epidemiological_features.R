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

custom_theme <- list(legend.position     ="none",
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


epi_table<-as_tibble(read.table(file="data/supplementary_data.csv", sep=",",  header=TRUE))

# Removing 3 Patients (see methods.rmd)
# patients_to_remove <- c("PAT182", "PAT202", "PAT230")
# epi_table <- epi_table %>% filter(! PATIENT_ID %in% patients_to_remove)

epi_table$SAMPLE_DATE <- as.Date(epi_table$SAMPLE_DATE)
epi_table$day <- as.numeric(format(epi_table$SAMPLE_DATE, format = "%d"))
epi_table$month <- as.numeric(format(epi_table$SAMPLE_DATE, format = "%m"))
epi_table$year <- as.numeric(format(epi_table$SAMPLE_DATE, format = "%Y"))


# computing the time distances between samples (in days)
epi_table<-epi_table %>%
  mutate(timeDist=unlist(pmap(list(day,month,year),distanceDate))) %>% 
  mutate(dateS=unlist(pmap(list(day,month,year),function(x,y,z){as.Date(str_c(c(z,y,x),collapse="-"))}))) %>%
  mutate(timeDec=unlist(pmap(list(day,month,year),function(x,y,z){decimal_date(as.Date(str_c(c(z,y,x),collapse="-")))}))) %>%
  mutate(PatSample=str_c(PATIENT_ID, SAMPLE_ID,sep="-")) %>% mutate(PatRes=str_c(PATIENT_ID, CARB,sep="-")) %>%
  mutate(PatSp=str_c(PATIENT_ID, SPECIES,sep="-"))
# Removing timedist outlier: (<0 days - this was because of an error - date of birth was keyed in as a sample date)
epi_table<-epi_table %>% mutate(timeDist = ifelse(timeDist < 0, NA, timeDist))

#############################################################
OXApatients<- epi_table %>%  select (PATIENT_ID) %>% unique ()

timedOXA<- epi_table %>% filter (PATIENT_ID %in% as.vector(OXApatients$PATIENT_ID) )



#####################################################################################
# Plot: Number of Samples per Patient Count Plot
#####################################################################################
epi_table_grouped_by_patient <- epi_table %>% 
  filter(epi_table$SAMPLE_ID != "JN626286") %>%
  group_by(PATIENT_ID) %>% 
  nest() %>% 
  mutate(n = map_dbl(.x=data,.f=~ length (.x$SAMPLE_ID)),
         meanTime=map_dbl(.x=data,.f=~ mean( .x$timeDec)),
         timediff=map_dbl(.x = data,.f=~ diff(c(min(.x$timeDist), max(.x$timeDist)))   ) )


plot_1 <- ggplot(epi_table_grouped_by_patient, aes(x=n)) +
  # geom_bar(fill="gray50") +
  geom_histogram(col="white", size=2, binwidth=1, fill="gray50") +
  theme_void() +
  do.call(theme, custom_theme) +
  ggtitle("1A. Number of Samples per Patient") +
  xlab("\nNumber of Samples per Patient") +
  ylab("Count") +
  ylim(0,100) +
  geom_hline(yintercept=0, color="grey")

ggsave("figures/number_of_samples_per_patient_histogram.pdf", width=default_figure_width, height=default_figure_height, units = "cm")








#####################################################################################
# Plot: Number of Samples per REGION Countplot
###################################################################################

epi_table_grouped_by_region <- epi_table  %>%
  filter(epi_table$SAMPLE_ID != "JN626286") %>%
  group_by(PATIENT_ID,REGION) %>%
  nest()
epi_table_grouped_by_region$REGION.Category = factor(epi_table_grouped_by_region$REGION, levels=c("LONDON", "YORK&HUM", "N EAST", "N WEST", "E MIDS", "W MIDS", "EAST", "S EAST", "S WEST", "Turkey" ))
epi_table_grouped_by_region <- epi_table_grouped_by_region[order(epi_table_grouped_by_region$REGION.Category),]

plot_2 <- ggplot(epi_table_grouped_by_region, aes(x=REGION.Category)) + 
  geom_bar(stat="count", fill="gray50") +

  theme_void() +
  do.call(theme, custom_theme) + # theme with custom_theme passed as arguments
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
  do.call(theme, custom_theme) +
  ggtitle("1C. Time Between Patient's First and Last Sample") +
  xlab("Time (days)") +
  ylab("Count\n") +
  geom_hline(yintercept=0, color="grey")

ggsave("figures/sample_time_difference_histogram.pdf", width=default_figure_width, height=default_figure_height, units = "cm")



##########################################################################################
# PLOT: Sample Date Histogram
##########################################################################################


timedOXAallMY <- epi_table %>% filter(epi_table$SAMPLE_ID != "JN626286")
timedOXAallMY <- timedOXAallMY %>% mutate(sample.date=as.Date(with(timedOXAallMY, paste(year,month,day,sep="-")),"%Y-%m-%d"))

plot_4 <- ggplot(timedOXAallMY %>% drop_na(sample.date) %>% filter(sample.date != "1947-11-23 "), aes(x=sample.date)) +
  geom_histogram(bins=36, col="white", size=1, fill="gray50") +
  theme_void() + 
  do.call(theme, custom_theme) +
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
  filter(epi_table$SAMPLE_ID != "JN626286") %>% 
  group_by(SPECIES) %>% 
  count() %>% 
  arrange(desc(n))
species_count_table$SpeciesGroup <- species_count_table$SPECIES
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
species_count_table$SPECIES <- species_count_table$SPECIES %>% str_replace_all(c('_' = ' '))




plot_5 <- ggplot(species_count_table, aes(x=n, y=reorder(SPECIES,n,decreasing=TRUE), fill=SpeciesGroup)) +
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
  do.call(theme, custom_theme) +
  # ggtitle("Distribution of Species") +
  ggtitle("1E. Sample Species") +
  xlab("Count") +
  ylab("Species") +
  geom_text(aes(x=n+1, label=SPECIES), size=4, hjust=0) +
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
strain_count_table <- epi_table %>% group_by(SEQUENCE_TYPE, SPECIES) %>% count() %>% arrange(desc(n))
strain_count_table$SpeciesGroup <- strain_count_table$SPECIES
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

strain_count_table <- strain_count_table %>% filter (SEQUENCE_TYPE!="-") %>% mutate(SpeciesStrain=paste0(SPECIES, "-",SEQUENCE_TYPE)) %>% arrange(desc(n)) %>% ungroup() 
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
  do.call(theme, custom_theme) +
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











##########################################################################################
# MULTI-PLOTS: Fig 1A - F
##########################################################################################
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
cm_to_inches <- 0.393701
pdf("figures/fig_1_a_to_f.pdf", width=default_figure_width*2*cm_to_inches, height=default_figure_height*3.75*cm_to_inches)
multiplot(plot_1, plot_3, plot_5, plot_2, plot_4, plot_6, cols=2)
dev.off()


png("figures/fig_1_a_to_f.png", width=default_figure_width*100, height=default_figure_height*100)
multiplot(plot_1, plot_3, plot_5, plot_2, plot_4, plot_6, cols=2)
dev.off()










######################################################################################
# table 1 on the number of non-OXA-48 resistances
#########################################################################################

View(epi_table %>% filter (CARB!="OXA-48") %>% group_by(CARB) %>% count() %>% arrange(n,desc=TRUE))





