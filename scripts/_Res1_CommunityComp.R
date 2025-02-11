
#######################################################
 ## COMMUNITY STRUCTURE
#######################################################

## BACTERIA

amp_heatmap(amp_subset_samples(
  amp.bac,sample_title!="NK" & time %in% c("0h","17h","49h")),   
  tax_aggregate = "Order",
  tax_show = c(
    "Micrococcales","Microtrichales",
    "Cytophagales","Synechococcales",
    "Flavobacteriales","Rhodobacterales",
    "Pseudomonadales","SAR11 clade","PeM15",
    "Chitinophagales","Verrucomicrobiales",
    "Enterobacterales","Puniceispirillales"),
  order_y_by = c(
    "Puniceispirillales",
    "Synechococcales","Verrucomicrobiales",
    "SAR11 clade","Pseudomonadales",
    "Rhodobacterales",    "Flavobacteriales",
    "PeM15","Micrococcales","Microtrichales",
    "Enterobacterales",    
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

#######################################################

## EUKARYOTES
# 3.4 x 5.5
amp_heatmap(amp_subset_samples(
  amp.euk, sample_title!="NK" & time %in% c("0h","17h","49h")),   
  tax_aggregate = "Order",
  tax_show = c(
   "Chalomydomonadales",
    "Telonemia uc","Mamiellales",
    "Bacillariophyta uc","Chloropicales",
    "Prymnesiales","Phaeocystales", 
    "Pterocystida",
   #"Dino-Group-I", "Gymnodiniales",
    "Peridiniales","Acanthoecida",
    "Suessiales",   "Cryomonadida",
    "Sarcinochrysidales","Pelagomonadales"),
  normalise = T,
  plot_values = F,
  facet_by = "condition",
  order_y_by = c(
    "Telonemia uc","Bacillariophyta uc",
    "Chloropicales","Mamiellales",
    "Prymnesiales", "Phaeocystales",
    "Dino-Group-I","Gymnodiniales",
   "Chalomydomonadales","Pterocystida",
    "Peridiniales", "Suessiales","Cryomonadida",
   "Acanthoecida","Sarcinochrysidales","Pelagomonadales"),
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
    axis.text.y = element_text(size=9.5),
    #axis.text.x = element_blank(),
    legend.position = "right",
    axis.ticks = element_blank())
