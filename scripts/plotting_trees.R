library(tidyverse)
library(ape)
library(ggtree)
library(treeio)
library(ggnewscale)



#####################################################################################
# LOAD: UKHSA Colours and Templates
#####################################################################################
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
# READ: TREE FILE
#####################################################################################
# Read in Tree File:
patient_tree<-read.tree("data/Final_alignment_masked_PAT_50bases_0.1samples.fst.treefile")

# Assigning names to the node labels (just to see where clades begin etc...) - tip labels are already defined in the file
patient_tree$node.label<-paste("node",1:patient_tree$Nnode, sep="-")

#####################################################################################
# READ: EPI DATA
#####################################################################################
# Read in epi file (Contains all the info we need to colour the columns - region, species, time, ...)
epi_table<-as_tibble(read.table(file="data/supplementary_data.csv", sep=",",  header=TRUE))

# Adding second column of information (bacterial species) - grouping the species together...
epi_table <- epi_table %>%
  mutate(SPECIES=map_chr(SPECIES, function(x){ ifelse(x %in% c("Escherichia_coli","Escherichia sp."), "Escherichia coli",x)})) %>%
  mutate(SPECIES=map_chr(SPECIES, function(x){ ifelse(x %in% c("Citrobacter amalonaticus","Citrobacter sp.","Citrobacter_freundii"),
                                                      "Citrobacter",x)})) %>%
  mutate(SPECIES=map_chr(SPECIES, function(x){ ifelse(x %in% c("Enterobacter_aerogenes","Enterobacter_cloacae" ,"Enterobacter_hormaechei" ),
                                                      "Enterobacter",x)})) %>%
  mutate(SPECIES=map_chr(SPECIES, function(x){ ifelse(x %in% c("Klebsiella_oxytoca","Klebsiella_variicola"),
                                                      "Other Klebsiella",x)})) %>%  
  mutate(SPECIES=map_chr(SPECIES, function(x){ ifelse(x %in% c("Morganella_morganii ","Raoultella planticola","Serratia marcescens" ),
                                                      "Other Enterobacterales",x)})) %>%
  mutate(SPECIES=map_chr(SPECIES, function(x){ ifelse(is.na(x),
                                                      "Other Enterobacterales",x)}))

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


#####################################################################################
# Tree Plot: NODE-LABELLED TREE
#####################################################################################
# We will require this tree to manually find where clades 'start from'

node_labelled_tree_plot <- ggtree(patient_tree) %<+% epi_table +
  geom_tiplab(
    size = 1,
    linesize=.2,
    align=TRUE) +
  geom_text2(aes(subset=!isTip, label=node), # labels all the nodes in the tree
             size = 1,
             color = "darkred", 
             hjust = 0, 
             vjust = 0) +
  geom_text2(aes(subset=isTip, label=node), # labels all the nodes in the tree
             size = 1,
             color = "darkblue", 
             hjust = 1, 
             vjust = 0) 

node_labelled_tree_plot

ggsave("figures/node_labelled_tree.pdf", width=18, height=28, units = "cm")

# New tree labels when the non-epidemic "variable" clade is collapsed

collapsed_node_labelled_tree_plot <- collapse(ggtree(patient_tree), node=449) %<+% epi_table +
  geom_tiplab(
    size = 1,
    linesize=.2,
    align=TRUE) +
  geom_text2(aes(subset=!isTip, label=node), # labels all the nodes in the tree
             size = 1,
             color = "darkred", 
             hjust = 0, 
             vjust = 0) +
  geom_text2(aes(subset=isTip, label=node), # labels all the nodes in the tree
             size = 1,
             color = "darkblue", 
             hjust = 1, 
             vjust = 0) 

collapsed_node_labelled_tree_plot

ggsave("figures/collapsed_node_labelled_tree.pdf", width=18, height=28, units = "cm")





#####################################################################################
# UPDATING TREE WITH GROUPS AND METADATA
#####################################################################################

# To view as Tibble dataframe:
# as_tibble(patient_tree) %>% print(n = 466)
# View(as_tibble(patient_tree))

# Adding metadata to table: Grouping:
total_nodes <- 464
small_clade <- c(217:233, 449:464) # c(218:234,451:466) old
epidemic_clade_node <- 235 # 236 (old - this one had JN.... .1)
smaller_clade_node <- 449 # 451 (old)
reference_clade_node <- 1

reference_clade <- c(reference_clade_node)
group_labels <- replicate(total_nodes, "epidemic")
group_labels[small_clade] = "other"
group_labels[reference_clade] = "reference"
# Adding metadata to table: Colours: (similar to group, but we want root of these clades not to have a colour...)
group_colours <- group_labels
group_colours[c(234, 235, 449, 1)] = "reference" # old: group_colours[c(235, 236, 451, 1)] = "reference"
# Adding metadata to table: Type of node: (root, tip, branch)
# node_type <- c(replicate(234,"tip"), c("root"), replicate(231, "branch"))
tree_metadata <- data.frame(
  node = 1:total_nodes,
  # node_type = node_type,
  group = group_labels,
  colour = group_colours
)
# Combining tree and metadata:
tree_table <- full_join(as_tibble(patient_tree), tree_metadata, by = c("node"="node"))
full_tree <- as.treedata(tree_table)











#####################################################################################
# Tree Plot: BASE TREE
#####################################################################################
# We will add epi data & focusing etc... on top of the base tree plot

base_tree_plot <- ggtree(full_tree) %<+% epi_table
grouped_tree_plot <- ggtree(full_tree, aes(color=colour)) %<+% epi_table













#####################################################################################
# Tree Plot 1: 2 Main  (Epidemic and 'Variable' clades) + Reference Clade
#####################################################################################
# Adding focuss of the two main clades

# Consult Node-Labelled Tree for node numbers
# Plot the tree with node labels
tree_plot_1 <- grouped_tree_plot +
  scale_color_manual(
    values=c(UKHSA_cols$moonlight, UKHSA_cols$tangerine, "black"),
    guide="none",
    expand = expansion(mult = c(0, 0.1), add = c(1, 0))
  ) +
  geom_tiplab(
    size = 1,
    linesize=.2,
    colour="black",
    align=TRUE
  ) +
  geom_hilight( # focuss the epidemic clade
    node = epidemic_clade_node,
    fill = UKHSA_cols$moonlight,
    alpha=0.1,
    extend=0.0178 # + 0.0003
  ) +
  # geom_cladelabel(
  #   node=epidemic_clade_node,
  #   label="Epidemic Clade",
  #   angle=90,
  #   offset=-.0192,
  #   offset.text=-0.0008,
  #   align=T,
  #   color=UKHSA_cols$moonlight,
  #   barsize=0.5) +
  geom_hilight( # focuss the smaller clade
    node = smaller_clade_node,
    fill = UKHSA_cols$cherry,
    alpha=0.1,
    extend=0.0011 # + 0.0003
  ) +
  geom_hilight( # focuss JN... "Clade"
    node = reference_clade_node,
    fill = UKHSA_cols$grass,
    alpha=0.1,
    extend=0.0192 # + 0.0003
  ) +
  theme(panel.background = element_rect(fill = NA)) +
  theme(plot.margin = unit(c(.2,2.9,.2,-.7), "cm")) # t,r,b,l


tree_plot_1


#####################################################################################
# Tree Plot 2: Adding Region Information
#####################################################################################
# Setup heatmap table:
region_table<-data.frame("Region"=epi_table$REGION)
rownames(region_table) <- epi_table$SAMPLE_ID
# Tree plot
tree_plot_2<-
  gheatmap(
    tree_plot_1, region_table, 
    offset = 0.0009,
    width = 0.05, 
    color="black", 
    colnames = TRUE,
    font.size=2, 
    colnames_angle=0,
    hjust=0.5,
    colnames_position="top",
    colnames_offset_y=0.8
  ) +
  scale_fill_manual(
    name = "Region",                       # define the coloring scheme and legend for gender
    # values = c("skyblue1", "royalblue", "maroon", "maroon4", "maroon3", "skyblue3","antiquewhite3", "antiquewhite4", "darkblue", "gray100"),
    # "YORK&HUM", "N WEST", "N EAST", "W MIDS", "E MIDS", "EAST", "LONDON", "S WEST", "S EAST"
    values = c("mediumpurple1", "slateblue4", "dodgerblue4", "dodgerblue", "lightskyblue", "aquamarine2", "honeydew2", "yellowgreen", "seagreen4", "white"),
    breaks = c("LONDON", "S WEST", "S EAST", "YORK&HUM", "N WEST", "N EAST", "W MIDS", "E MIDS", "EAST", "Turkey" ),
    labels = c("LONDON", "S WEST", "S EAST", "YORK&HUM", "N WEST", "N EAST", "W MIDS", "E MIDS", "EAST", "Turkey (Reference)" )
  )+
  theme(
    legend.position = "bottom",
    legend.title = element_text(face="bold", size=12),
    legend.text = element_text(size=8),
    axis.title = element_blank(),
    axis.text = element_text(size=8, face="bold")
  )
# plot with first column of information
tree_plot_2








#####################################################################################
# Tree Plot 3: Adding Species Information
#####################################################################################
# Setup heatmap table:
species_table<-data.frame("Species"=epi_table$SPECIES)
rownames(species_table) <- epi_table$SAMPLE_ID
# Plot:
# actually adding it to the plot
tree_plot_2a <- tree_plot_2 + new_scale_fill()

tree_plot_3 <-
  gheatmap(
    tree_plot_2a, species_table, 
    offset = 0.0022,
    width = 0.05, 
    color="black", 
    colnames = TRUE,
    font.size=2,
    colnames_angle=0,
    hjust=0.5,
    colnames_position="top",
    colnames_offset_y=0.8
  ) +
  scale_fill_manual(
    name = "Bacterial Species",
    values = c("darkred", "red3", "lightcoral", "darkorange", "gold", "wheat"),
    breaks = c("Klebsiella_pneumoniae", "Other Klebsiella", "Enterobacter", "Other Enterobacterales", "Escherichia coli", "Citrobacter"),
    labels = c("Klebsiella pneumoniae", "Other Klebsiella", "Enterobacter", "Other Enterobacterales", "Escherichia coli", "Citrobacter"),
    na.value = "white"
  )+
  theme(#legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    #legend.box = "vertical", legend.margin = margin()
  )
#+
#  guides(fill = guide_legend(nrow = 2,byrow = TRUE))

tree_plot_3








#####################################################################################
# Tree Plot 4: Adding Time/Date Information
#####################################################################################

date_table<-data.frame("Time"=epi_table$timeDec)
rownames(date_table) <- epi_table$SAMPLE_ID

tree_plot_3a <-tree_plot_3 + new_scale_fill()

tree_plot_4 <- gheatmap(
  tree_plot_3a,
  date_table, 
  offset = 0.00351,
  width = 0.05, 
  color="black", 
  colnames = TRUE,
  font.size=2,
  colnames_angle=0,
  hjust=0.5,
  colnames_position="top",
  colnames_offset_y=0.8
) +
  scale_fill_continuous(name = "Date",
                        low = "gray90", high ="black",
                        breaks = c(2014,2014.5,2015,2015.5,2016,2016.5),
                        labels=c("2014","","2015","","2016",""),
                        na.value = "white")+
  guides(fill = guide_colourbar(barwidth = .8, barheight = 5))+
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        #legend.position = "bottom",
        legend.position = c(1.1,.765),
        legend.title = element_text(size = 6,face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(c(.4,.4), "cm"),
        #legend.box = "horizontal", legend.margin = margin(),
        legend.background = element_rect(fill='transparent')
  ) +
  ggtree::vexpand(.001, 1) # adds more space to top (so colnames aren't cut off)
#+
# guides(color = guide_legend(override.aes = list(size = 0.5),))
# final plot
# beware of ylim, it can cut out samples (but then it gives an error)
tree_plot_4



# ggsave("figures/main_phylo_tree.pdf", width=25, height=25, units = "cm")
ggsave("figures/main_phylo_tree.png", width=25, height=25, units = "cm")





#####################################################################################
# Tree Plot 5: Outbreak Tree - Clades Labelled...
#####################################################################################
# Now we collapse the tree and focus on the Epidemic clade
# We focus a few clades of interest

subclade_tree<- collapse(ggtree(full_tree), node=449) %<+% epi_table+
  geom_tiplab(size = 1,
              offset = 0.00004,
              align = TRUE,
              linesize=.2)+
  geom_hilight(  # focuss node 39 in blue, "extend =" allows us to define the length of the color block
    node = 269,
    fill = UKHSA_cols$plum,alpha=0.25,extend=0.00134) +
  geom_cladelab(
    node=269, label="A", offset = 0.00008, offset.text=0, fontsize=5, barcolour=rgb(red = 0, green = 0, blue = 0, alpha = 0),
  ) +
  geom_hilight(  # focuss node 39 in blue, "extend =" allows us to define the length of the color block
    node = 380,
    fill = UKHSA_cols$DHSC.green ,alpha=0.25,extend=0.00014) +
  geom_cladelab( # Using 381 instead of 380 as it will position better
    node=381, label="B", offset = 0.00008, offset.text=0, fontsize=5, barcolour=rgb(red = 0, green = 0, blue = 0, alpha = 0),
  ) +
  geom_hilight(  # focuss node 39 in blue, "extend =" allows us to define the length of the color block
    node = 364,
    fill = UKHSA_cols$cherry,alpha=0.25,extend=0.000949) +
  geom_cladelab(
    node=364, label="C", offset = 0.00008, offset.text=0, fontsize=5, barcolour=rgb(red = 0, green = 0, blue = 0, alpha = 0),
  ) +
  geom_hilight(  # focuss node 39 in blue, "extend =" allows us to define the length of the color block
    node = 386,
    fill = UKHSA_cols$ocean,alpha=0.25,extend=0.001273) +
  geom_cladelab(
    node=386, label="D", offset = 0.00008, offset.text=0, fontsize=5, barcolour=rgb(red = 0, green = 0, blue = 0, alpha = 0),
  ) +
  geom_hilight(  # focuss node 39 in blue, "extend =" allows us to define the length of the color block
    node = 438,
    fill = UKHSA_cols$sunny,alpha=0.25,extend=0.00136) +
  geom_cladelab(
    node=438, label="E", offset = 0.00008, offset.text=0, fontsize=5, barcolour=rgb(red = 0, green = 0, blue = 0, alpha = 0),
  ) +
  geom_hilight(  # focuss node 39 in blue, "extend =" allows us to define the length of the color block
    node = 433,
    fill = UKHSA_cols$midnight,alpha=0.25,extend=0.00088) +
  geom_cladelab(
    node=433, label="F", offset = 0.00008, offset.text=0, fontsize=5, barcolour=rgb(red = 0, green = 0, blue = 0, alpha = 0),
  ) +
  geom_point2(aes(subset=(node==449)), shape=18, size=3, fill=UKHSA_cols$tangerine) +
  theme(plot.margin= unit(c(.2,2.5,.2,-.7), "cm"))

subclade_tree

ggsave("figures/subclades_in_epidemic_phylo_tree.pdf", width=25, height=25, units = "cm")

#####################################################################################
# SUBCLADES + REGION
#####################################################################################

subclade_tree_2<-
  gheatmap(
    subclade_tree, region_table, 
    offset = 0.00012,
    width = 0.05, 
    color="black", 
    colnames = TRUE,
    font.size=2, 
    colnames_angle=0,
    hjust=0.5,
    colnames_position="bottom",
    colnames_offset_y=220
  ) +
  scale_fill_manual(
    name = "Region",
    values = c("mediumpurple1", "slateblue4", "dodgerblue4", "dodgerblue", "lightskyblue", "aquamarine2", "honeydew2", "yellowgreen", "seagreen4", "white"),
    breaks = c("LONDON", "S WEST", "S EAST", "YORK&HUM", "N WEST", "N EAST", "W MIDS", "E MIDS", "EAST", "Turkey" ),
    labels = c("LONDON", "S WEST", "S EAST", "YORK&HUM", "N WEST", "N EAST", "W MIDS", "E MIDS", "EAST", "Turkey (Reference)" )
  )

subclade_tree_2


#####################################################################################
# SUBCLADES + REGION + SPECIES
#####################################################################################

subclade_tree_2a <- subclade_tree_2 + new_scale_fill()

subclade_tree_3 <-
  gheatmap(
    subclade_tree_2a, species_table, 
    offset = 0.00023,
    width = 0.05, 
    color="black", 
    colnames = TRUE,
    font.size=2,
    colnames_angle=0,
    hjust=0.5,
    colnames_position="bottom",
    colnames_offset_y=220
  ) +
  scale_fill_manual(
    name = "Bacterial Species",
    values = c("darkred", "red3", "lightcoral", "darkorange", "gold", "wheat"),
    breaks = c("Klebsiella_pneumoniae", "Other Klebsiella", "Enterobacter", "Other Enterobacterales", "Escherichia coli", "Citrobacter"),
    labels = c("Klebsiella pneumoniae", "Other Klebsiella", "Enterobacter", "Other Enterobacterales", "Escherichia coli", "Citrobacter"),
    na.value = "white"
  )

subclade_tree_3

#####################################################################################
# SUBCLADES + REGION + SPECIES + TIME
#####################################################################################

subclade_tree_3a <-subclade_tree_3 + new_scale_fill()

subclade_tree_4 <- gheatmap(
  subclade_tree_3a,
  date_table, 
  offset = 0.000339,
  width = 0.05, 
  color="black", 
  colnames = TRUE,
  font.size=2,
  colnames_angle=0,
  hjust=0.5,
  colnames_position="bottom",
  colnames_offset_y=220
) +
  scale_fill_continuous(name = "Date",
                        low = "gray90", high ="black",
                        breaks = c(2014,2014.5,2015,2015.5,2016,2016.5),
                        labels=c("2014","","2015","","2016",""),
                        na.value = "white")+
  guides(fill = guide_colourbar(barwidth = .8, barheight = 5))+
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        #legend.position = "bottom",
        legend.position = c(1.09,.7895),
        legend.title = element_text(size = 5.5,face = "bold"),
        legend.text = element_text(size = 5.5),
        legend.key.size = unit(c(.3,.3), "cm"),
        #legend.box = "horizontal", legend.margin = margin(),
        legend.background = element_rect(fill='transparent')
  ) +
  ggtree::vexpand(.01, 1) # adds more space to top (so colnames aren't cut off)

subclade_tree_4

ggsave("figures/subclades_in_epidemic_phylo_tree_with_metadata.pdf", width=25, height=25, units = "cm")









##############################################################################################
#
# FOCUS PATIENTS....
#
##############################################################################################



focus_patients <- c("PAT09", "PAT106", "PAT142", "PAT148", "PAT157", "PAT160", "PAT175", "PAT18", "PAT185", "PAT192", "PAT219", "PAT221", "PAT229", "PAT235", "PAT241", "PAT28", "PAT52", "PAT90")

focus_patients_tree <- patient_tree
focus_patients_tree$tip.label <- vapply(strsplit(focus_patients_tree$tip.label,"-"), `[`, 1, FUN.VALUE=character(1))

focus_patients_1_shape_palette <- c(
  PAT106   = 15,
  PAT219   = 15,
  PAT229   = 15,
  PAT241   = 15,
  PAT52    = 15,
  PAT90    = 15,
  JN626286 = 15,
  other    = 16
)

focus_patients_1_fill_palette <- c(
  PAT106   = "#a6cee3",
  PAT219   = "#33a02c",
  PAT229   = "#e31a1c",
  PAT241   = "#ff7f00",
  PAT52    = "#6a3d9a",
  PAT90    = "#1f78b4",
  JN626286 = "gray30",
  other    = "#FFFFFF00"
)




#####################################################################################
# TREE PART 1: FOCUSING ON SPECIFIC PATIENTS WITH AT LEAST TWO GENETICALLY DIFFERENT SAMPLES
# Patients in variable clade
#####################################################################################
# Patients who appear in variable clade with at least 2 genetically distinct samples
focus_patients_1_without_reference <- c(
  "PAT219",
  "PAT52",
  "PAT241",
  "PAT90",
  "PAT106",
  "PAT229"
)
focus_patients_1 <- c(
  "PAT219",
  "PAT52",
  "PAT241",
  "PAT90",
  "PAT106",
  "PAT229",
  "JN626286"
)
# Making tree
focus_patients_tree_1 <- full_join(
  as_tibble(focus_patients_tree),
  data.frame(
    label=focus_patients_1,
    group=focus_patients_1),
  by = c("label"="label")
)

focus_patients_tree_1$group[is.na(focus_patients_tree_1$group)] <- "Patients"
focus_patients_tree_1$label <- "\u2014" 
focus_patients_tree_1 <- as.treedata(focus_patients_tree_1)


# Plotting tree
focus_patients_tree_1_plot <- ggtree(focus_patients_tree_1) %<+% epi_table
focus_patients_tree_1_plot +
  geom_tiplab(
    aes(subset = group != "Patients", color=group, fill=group),
    size=9,
    fontface=2,
    show.legend=FALSE,
    vjust=0.35
  ) +
  geom_tippoint(aes(color=group, shape=group), size=0.5, stroke = 0.1) +
  scale_shape_manual(values=focus_patients_1_shape_palette) +
  scale_fill_manual(values=focus_patients_1_fill_palette ) +
  scale_color_manual(values=focus_patients_1_fill_palette ) +
  guides(color = guide_legend(override.aes = list(size=5))) +
  theme(legend.position = c(0.8, 0.8)) +
  theme(legend.title=element_blank()) +
  ggtree::vexpand(.005, -1) # adds more space to top (so colnames aren't cut off)

ggsave("figures/focus_patients_tree_part_1", width=25, height=25, units = "cm")













#####################################################################################
# TREE PART 2: FOCUSING ON SPECIFIC PATIENTS WITH AT LEAST TWO GENETICALLY DIFFERENT SAMPLES
# Patients in epidemic clade
#####################################################################################
# Patients who appear in variable clade with at least 2 genetically distinct samples
# focus_patients_2 <- setdiff(focus_patients, focus_patients_1_without_reference)

focus_patients_2 <- c(
  "PAT09",
  "PAT142",
  "PAT148",
  "PAT157",
  "PAT160",
  "PAT175",
  "PAT18",
  "PAT185",
  "PAT192",
  "PAT221",
  "PAT235",
  "PAT28",
  "JN626286"
)


focus_patients_2_shape_palette <- c(
  PAT09  = 15,
  PAT142 = 15,
  PAT148 = 15,
  PAT157 = 15,
  PAT160 = 15,
  PAT175 = 15,
  PAT18  = 15,
  PAT185 = 15,
  PAT192 = 15,
  PAT221 = 15,
  PAT235 = 15,
  PAT28  = 15,
  JN626286 = 15,
  other    = 16
)

focus_patients_2_fill_palette <- c(
  PAT09  = "#a6cee3",
  PAT142 = "#1f78b4",
  PAT148 = "#b2df8a",
  PAT157 = "#33a02c",
  PAT160 = "#fb9a99",
  PAT175 = "#e31a1c",
  PAT18  = "#fdbf6f",
  PAT185 = "#ff7f00",
  PAT192 = "#cab2d6",
  PAT221 = "#6a3d9a",
  PAT235 = "#ffff99",
  PAT28  = "#b15928",
  JN626286 = "gray30",
  other    = "#FFFFFF00"
)

focus_patients_2_fill_list <- c(
  "#a6cee3",
  "#1f78b4",
  "#b2df8a",
  "#33a02c",
  "#fb9a99",
  "#e31a1c",
  "#fdbf6f",
  "#ff7f00",
  "#cab2d6",
  "#6a3d9a",
  "#ffff99",
  "#b15928",
  "black"
)

focus_patients_2_label_list <- c(
  "PAT09",
  "PAT142",
  "PAT148",
  "PAT157",
  "PAT160",
  "PAT175",
  "PAT18",
  "PAT185",
  "PAT192",
  "PAT221",
  "PAT235",
  "PAT28",
  "Reference"
)

focus_patients_2_colour_palette <- c(
  PAT09  = "black",
  PAT142 = "black",
  PAT148 = "black",
  PAT157 = "black",
  PAT160 = "black",
  PAT175 = "black",
  PAT18  = "black",
  PAT185 = "black",
  PAT192 = "black",
  PAT221 = "black",
  PAT235 = "black",
  PAT28  = "black",
  JN626286 = "black",
  other    = "#FFFFFF00"
)


# Making tree
focus_patients_tree_2 <- full_join(
  as_tibble(focus_patients_tree),
  data.frame(
    label=focus_patients_2,
    group=focus_patients_2),
  by = c("label"="label")
)
focus_patients_tree_2$group[is.na(focus_patients_tree_2$group)] <- "Patients"

focus_patients_tree_2$label <- "\u2014" # "-" # "--------" # "◄--" # "⟵"

focus_patients_tree_2 <- as.treedata(focus_patients_tree_2)


# Plotting tree
focus_patients_tree_2_plot <- collapse(ggtree(focus_patients_tree_2), node=449) %<+% epi_table # ggtree(focus_patients_tree_2) %<+% epi_table
focus_patients_tree_2_plot +
  geom_tiplab(
    aes(subset = group != "Patients", color=group, fill=group),
    size=9,
    fontface=2,
    show.legend=FALSE,
    vjust=0.35
  ) +
  # geom_label2(
  #   aes(subset = group != "Patients", color=group),
  #   label="-",
  #   label.padding = unit(0, "lines"),
  #   label.r = unit(0, "lines"),
  #   label.size = 0,
  #   ) +
  # geom_tiplab(aes(color=group),
  #             size = 1.5,
  #             linesize=.4,
  #             fontface = 2,
#             # align=TRUE
# ) +

# geom_tippoint(aes(color=group, fill=group, shape=group), size=1.4, stroke = .2) +

geom_tippoint(aes(color=group, shape=group), size=0.5, stroke = 0.1) +
  
  geom_cladelabel(node=449, label="Collapsed Diverse Clade", 
                  color="black", align=FALSE, offset = .000025) + 
  
  scale_shape_manual(values=focus_patients_2_shape_palette) +
  scale_fill_manual(values=focus_patients_2_fill_palette ) +
  scale_color_manual(values=focus_patients_2_fill_palette ) +
  
  # scale_color_manual(
  #   name = "legend",
  #   values = focus_patients_2_fill_list,
  #   limits = focus_patients_2_label_list
  # ) +
  
  geom_point2(aes(subset=(node==449)), shape=18, size=7, fill=UKHSA_cols$tangerine) +
  
  guides(color = guide_legend(override.aes = list(size=5))) +
  
  # theme(legend.position = "none")
  theme(legend.position = c(0.8, 0.8)) +
  # theme(legend.key.size = unit(1, "cm")) +
  theme(legend.title=element_blank()) +
  # scale_shape_manual(name = "", values=18, label="Collapsed Variable Clade") +
  ggtree::vexpand(.005, -1) # adds more space to top (so colnames aren't cut off)

ggsave("figures/focus_patients_tree_part_2.pdf", width=25, height=25, units = "cm")

