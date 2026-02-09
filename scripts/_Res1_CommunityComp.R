
#######################################################
 ## COMMUNITY STRUCTURE
#######################################################

## BACTERIA -- Figure 2a
plot1a = amp_heatmap(amp_subset_samples(
  amp.bac,sample_title!="NK" & time %in% c("0h","17h","49h")),   
  tax_aggregate = "Order",
  tax_show = c(
    "Micrococcales","Microtrichales",
    #"Puniceispirillales","PeM15","Synechococcales",
    "Flavobacteriales","Rhodobacterales",
    "Pseudomonadales","SAR11 clade",
    "Chitinophagales","Verrucomicrobiales",
    "Cytophagales","Enterobacterales"),
  order_y_by = c(
    "Verrucomicrobiales","SAR11 clade",
    "Pseudomonadales","Rhodobacterales",
    "Flavobacteriales","Micrococcales",
    "Microtrichales","Enterobacterales",
    "Cytophagales","Chitinophagales"),
  normalise = T,
  plot_values = F,
  facet_by = "condition",
  round = 1,
  max_abundance = 55,
  min_abundance = 0.5,
  plot_colorscale = "log10",
  plot_legendbreaks = c(1,5,10,25,40),
  group_by = "time", 
  color_vector = c( 
    "gray97","aliceblue",#"lightsteelblue2",
    "lightsteelblue3","#ab6c9e",#"#8f117d",#c230ad
    "#610053","#36002e")) +
  geom_text(aes(
      label = round(Abundance),
      color = ifelse(Abundance < 4.5,"black","white"))) +
  scale_color_identity() +
  theme_classic() +
  theme(
    axis.text = element_text(size=9.5),
    #axis.text.x = element_blank(),
    legend.position = "right",
    axis.ticks = element_blank())

##########################

## EUKARYOTES -- Figure 2b
plot1b = amp_heatmap(amp_subset_samples(
  amp.euk, sample_title!="NK" & time %in% c("0h","17h","49h")),   
  tax_aggregate = "Order",
  tax_show = c(
    "Chlamydomonadales","Mamiellales",
    "Bacillariophyta uc","Chloropicales",
    "Cryomonadida","Prymnesiales","Phaeocystales", 
    #"Pterocystida","Sarcinochrysidales","Telonemia uc",
    #"Dino-Group-I","Acanthoecida","Gymnodiniales",
    "Peridiniales","Suessiales","Pelagomonadales"),
  normalise = T,
  plot_values = F,
  facet_by = "condition",
  order_y_by = c(
    "Bacillariophyta uc",
    "Chloropicales","Mamiellales",
    "Prymnesiales", "Phaeocystales",
    #"Pterocystida","Sarcinochrysidales","Telonemia uc",
    #"Dino-Group-I","Acanthoecida","Gymnodiniales",
    "Chlamydomonadales","Cryomonadida",
    "Peridiniales","Suessiales","Pelagomonadales"),
  round = 0,
  max_abundance = 55,
  min_abundance = 0.5,
  plot_colorscale = "log10",
  plot_legendbreaks = c(1,5,10,25,40),
  group_by = "time", 
  color_vector = c( 
    "gray97","aliceblue",#"lightsteelblue2",
    "lightsteelblue3","#ab6c9e",#"#8f117d",#c230ad
    "#610053","#36002e")) +
  geom_text(aes(
    label = round(Abundance),
    color = ifelse(Abundance < 4.5,"black","white"))) +
  scale_color_identity() +
  theme_classic() +
  theme(
    axis.text.y = element_text(size=9.5),
    #axis.text.x = element_blank(),
    legend.position = "right",
    axis.ticks = element_blank())

# Export size 6x6
plot_grid(
  plot1a, plot1b,
  ncol = 1,
  align = "hv",
  axis = "tblr")

#######################################################

# ALL REPLICATES -- Figure S1
plot2a = amp_heatmap(amp_subset_samples(
  amp.bac,sample_title!="NK" & time %in% c("0h","17h","49h")),   
  tax_aggregate = "Order",
  tax_show = c(
    "Micrococcales","Microtrichales",
    "Cytophagales","Flavobacteriales",
    #"Synechococcales","Puniceispirillales","PeM15",
    "Rhodobacterales","Pseudomonadales",
    "SAR11 clade","Chitinophagales",
    "Verrucomicrobiales","Enterobacterales"),
  order_y_by = c(
    "Verrucomicrobiales","SAR11 clade",
    "Pseudomonadales","Rhodobacterales",
    "Flavobacteriales","Micrococcales",
    "Microtrichales","Enterobacterales",    
    "Cytophagales","Chitinophagales"),
  normalise = T,
  plot_values = F,
  facet_by = c("time","condition"),
  round = 1,
  max_abundance = 55,
  min_abundance = 0.5,
  plot_colorscale = "log10",
  plot_legendbreaks = c(1,5,10,25,40),
  group_by = "replicate", 
  color_vector = c( 
    "gray97","aliceblue",#"lightsteelblue2",
    "lightsteelblue3","#ab6c9e",#"#8f117d",#c230ad
    "#610053","#36002e")) +
  geom_text(aes(
    label = round(Abundance),
    color = ifelse(Abundance < 4.5,"black","white"))) +
  scale_color_identity() +
  labs(x="Replicate") +
  theme_classic() +
  theme(
    axis.text = element_text(size=9.5),
    #axis.text.x = element_blank(),
    legend.position = "right",
    axis.ticks = element_blank())

plot2b = amp_heatmap(amp_subset_samples(
  amp.euk, sample_title!="NK" & time %in% c("0h","17h","49h")),   
  tax_aggregate = "Order",
  tax_show = c(
    "Chlamydomonadales","Mamiellales",
    "Bacillariophyta uc","Chloropicales",
    "Prymnesiales","Phaeocystales", 
    #"Pterocystida","Sarcinochrysidales","Telonemia uc",
    #"Dino-Group-I","Acanthoecida","Gymnodiniales",
    "Peridiniales","Suessiales",   
    "Cryomonadida","Pelagomonadales"),
  normalise = T,
  plot_values = F,
  facet_by = c("time","condition"),
  order_y_by = c("Bacillariophyta uc",
    "Chloropicales","Mamiellales",
    "Prymnesiales", "Phaeocystales",
    "Chlamydomonadales","Cryomonadida",
    "Peridiniales","Suessiales","Pelagomonadales"),
  round = 0,
  max_abundance = 55,
  min_abundance = 0.5,
  plot_colorscale = "log10",
  plot_legendbreaks = c(1,5,10,25,40),
  group_by = "replicate", 
  color_vector = c( 
    "gray97","aliceblue",#"lightsteelblue2",
    "lightsteelblue3","#ab6c9e",#"#8f117d",#c230ad
    "#610053","#36002e")) +
  labs(x="replicate") +
  geom_text(aes(
    label = round(Abundance),
    color = ifelse(Abundance < 4.5,"black","white"))) +
  scale_color_identity() +
  theme_classic() +
  theme(
    axis.text.y = element_text(size=9.5),
    #axis.text.x = element_blank(),
    legend.position = "right",
    axis.ticks = element_blank())

# Export size 9x7
plot_grid(
  plot2a, plot2b,
  ncol=1,
  align = "hv",
  axis = "tblr")
