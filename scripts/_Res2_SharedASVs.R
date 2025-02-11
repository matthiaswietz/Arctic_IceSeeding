#############################################################
 ## SHARED AND UNIQUE ASVs -- BAC
#############################################################

# org. seawater + control
core1 <- ASV.bac.rel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(clip_id=variable) %>%
  left_join(ENV.bac) %>%
  filter(condition %in% c("OrgSW","Ctr") & time %in% c("0h","17h","49h")) %>%
  group_by(asv) %>% 
  summarize(Abundance=mean(Abundance)) %>% 
  filter(Abundance >= 0.01) %>%
  distinct(asv, .keep_all=T)  %>%
  column_to_rownames("asv")

# org ice
core2 <- ASV.bac.rel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(clip_id=variable) %>%
  left_join(ENV.bac) %>%
  filter(condition=="OrgIce") %>%
  group_by(asv) %>% 
  summarize(Abundance=mean(Abundance)) %>% 
  filter(Abundance >= 0.01) %>%
  distinct(asv, .keep_all=T)  %>%
  column_to_rownames("asv")

# ice melt
core3 <- ASV.bac.rel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(clip_id=variable) %>%
  left_join(ENV.bac) %>%
  filter(condition=="IceMelt" & time %in% c("17h","49h")) %>%
  group_by(asv) %>% 
  summarize(Abundance=mean(Abundance)) %>% 
  filter(Abundance >= 0.01) %>%
  distinct(asv, .keep_all=T)  %>%
  column_to_rownames("asv")

# Combine in list
core.bac <- list()
core.bac[["SW"]] <- as.character(row.names(core1))
core.bac[["OrgIce"]] <- as.character(row.names(core2))
core.bac[["IceMelt"]] <- as.character(row.names(core3))

overlaps.bac <- presAbs(core.bac) %>%
  rownames_to_column("asv")

upset(
  overlaps.bac,
  number.angles = 0, 
  nsets=3,
  main.bar.color = "gray12",
  sets.bar.color = "gray12",
  matrix.color = "gray12",
  point.size = 2.44, line.size = 0.8, text.scale = 1.2,
  mainbar.y.label = "Bacterial ASVs",
  order.by = "freq") 

# Categorize by occurrence; add taxinfo
nicheDist.bac <- overlaps.bac %>%
  mutate(niche = case_when(
    OrgIce == 1 & IceMelt == 1 & SW == 0 ~"ice-exported",
    OrgIce == 1 & IceMelt == 0 & SW == 0 ~"ice-unique",
    OrgIce == 0 & IceMelt == 1 & SW == 0 ~"unclear",
    SW == 1 & OrgIce == 0 & IceMelt == 1  ~"seawater-unique",
    SW == 1 & OrgIce == 0 & IceMelt == 0  ~"seawater-unique",
    SW == 1 & OrgIce == 1 & IceMelt == 1  ~"seawater+ice",
    SW == 1 & OrgIce == 1 & IceMelt == 0  ~"seawater+ice")) %>%
  left_join(TAX.bac %>% mutate(
    asv = rownames(TAX.bac)), overlaps, by='asv') 

# Add abundances & metadata
nicheAbd.bac <- ASV.bac.rel %>%
  rownames_to_column(var = "asv") %>%
  pivot_longer(cols = -asv, names_to="sample", values_to="abundance") %>%
  right_join(nicheDist.bac) %>%
  left_join(ENV.bac, by =c("sample"="clip_id")) %>%
  drop_na(abundance) %>%
  mutate(type="Bacteria")

# Most common niche per Genus; ASVs per niche
nicheGenus.bac <- nicheDist.bac %>%
  left_join(topBac) %>%
  drop_na(mean) 


#############################################################
## SHARED AND UNIQUE ASVs -- EUK
#############################################################

# org. seawater & control
core4 <- ASV.euk.rel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(clip_id=variable) %>%
  left_join(ENV.euk) %>%
  filter(condition %in% c("OrgSW","Ctr") & time %in% c("0h","17h","49h")) %>%
  group_by(asv) %>% 
  summarize(Abundance=mean(Abundance)) %>% 
  filter(Abundance >= 0.01) %>%
  distinct(asv, .keep_all=T)  %>%
  column_to_rownames("asv")

# org ice
core5 <- ASV.euk.rel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(clip_id=variable) %>%
  left_join(ENV.euk) %>%
  filter(condition=="OrgIce") %>%
  group_by(asv) %>% 
  summarize(Abundance=mean(Abundance)) %>% 
  filter(Abundance >= 0.01) %>%
  distinct(asv, .keep_all=T)  %>%
  column_to_rownames("asv")

# icemelt
core6 <- ASV.euk.rel %>%
  rownames_to_column("asv") %>%
  reshape2::melt(value.name="Abundance") %>%
  dplyr::rename(clip_id=variable) %>%
  left_join(ENV.euk) %>%
  filter(condition=="IceMelt" & time %in% c("17h","49h")) %>%
  group_by(asv) %>% 
  summarize(Abundance=mean(Abundance)) %>% 
  filter(Abundance >= 0.01) %>%
  distinct(asv, .keep_all=T)  %>%
  column_to_rownames("asv")

# Combine in list
core.euk <- list()
core.euk[["SW"]] <- as.character(row.names(core4))
core.euk[["OrgIce"]] <- as.character(row.names(core5))
core.euk[["IceMelt"]] <- as.character(row.names(core6))

# Combine
overlaps.euk <- presAbs(core.euk) %>%
  rownames_to_column("asv")

upset(
  overlaps.euk,
  number.angles = 0, 
  nsets=4,
  main.bar.color = "gray12",
  sets.bar.color = "gray12",
  matrix.color = "gray12",
  point.size = 2.44, line.size = 0.8, text.scale = 1.2,
  mainbar.y.label = "Eukaryotic ASVs",
  order.by = "freq") 

#############################################

# Categorize by occurrence; add taxinfo
nicheDist.euk <- overlaps.euk %>%
  mutate(niche = case_when(
    OrgIce == 1 & IceMelt == 1 & SW == 0 ~"ice-exported",
    OrgIce == 1 & IceMelt == 0 & SW == 0 ~"ice-unique",
    OrgIce == 0 & IceMelt == 1 & SW == 0 ~"unclear",
    SW == 1 & OrgIce == 0 & IceMelt == 1  ~"seawater-unique",
    SW == 1 & OrgIce == 0 & IceMelt == 0  ~"seawater-unique",
    SW == 1 & OrgIce == 1 & IceMelt == 1  ~"seawater+ice",
    SW == 1 & OrgIce == 1 & IceMelt == 0  ~"seawater+ice")) %>%
  left_join(TAX.euk %>% mutate(
    asv = rownames(TAX.euk)), overlaps, by='asv') 

# Add abundances & metadata
nicheAbd.euk <- ASV.euk.rel %>%
  rownames_to_column(var = "asv") %>%
  pivot_longer(cols = -asv, names_to="sample", values_to="abundance") %>%
  right_join(nicheDist.euk) %>%
  left_join(ENV.euk, by =c("sample"="clip_id")) %>%
  drop_na(abundance) %>%
  mutate(type="Eukaryotes")

# Preferred niche per Genus; based on most abundant ASV
nicheGenus.euk <- nicheDist.euk %>%
  left_join(topEuk) %>%
  drop_na(mean)


#############################################################
 ## COMBINE # PLOT
#############################################################

# Keep genera with mean >0.5 in min. 1 condition
# Remove some uc taxa (especially if related, known genera are present)
# Focus on taxa with consistent detection
# export size 5x6
rbind(nicheAbd.bac, nicheAbd.euk) %>%
  filter(niche %in% c("ice-exported") & time %in% c("0h","17h","49h")) %>%
  group_by(asv, condition, time, Genus, type, niche) %>%
  summarize(sum = sum(abundance, na.rm=T), .groups="drop") %>%
  group_by(condition, time, Genus, type, niche) %>%
  summarize(mean = mean(sum, na.rm=T), .groups="drop") %>%
  group_by(Genus) %>%
  filter(any(mean > 0.2)) %>%  
  ungroup() %>%
  filter(!Genus %in% c(
    "JGI 0000069-P22 uc","ML602J-51","Flavobacteriaceae uc",
    "Cryomorphaceae uc","Saprospiraceae uc","Rhodobacteraceae uc",
    "Cellvibrionaceae uc","MAST-1C uc","Gymnodiniaceae uc",
    "Raphid-pennate uc","Naviculaceae uc","Botuliformida uc",
    "Dinophyceae uc","Filosa-Thecofilosea uc","Cercozoa uc",
    "Protaspa-lineage uc","MOCH-2 uc","Botuliformidae uc",
    "Cryothecomonas-lineage uc","Apocalathium","Amphidiniopsis",
    "Suessiaceae uc","Uronychia","Peritromus","Pedinellales uc",
    "Marimonadida uc","Mataza-lineage uc","Telonemia-Group-1 uc",
    "Pavlovaceae uc","Dino-Group-II uc","TAGIRI1-lineage uc")) %>% 
  ggplot(aes(x=time, y=reorder(Genus, desc(-mean)), fill=mean)) +  # reorder Genus by mean
  geom_tile(color = "white") +
  facet_wrap(type ~ ., scales = "free", ncol=2) +
  scale_fill_gradientn(colors = c(
    "whitesmoke","lemonchiffon1","lemonchiffon2",
    "lightsteelblue1","lightsteelblue3","deeppink4"),
    trans = "log10",
    na.value="white",
    breaks = c(0, 0.01, 0.05, 1, 5, 10)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title = element_blank(),
    axis.ticks = element_blank())

## Export ice-specific ASVs for SI table
rbind(nicheAbd.bac, nicheAbd.euk) %>%
  filter(niche %in% c("ice-exported") & time %in% c("0h","17h","49h")) %>%
  write.table()

######################################

# Contribution of each fraction
nicheSummary <- rbind(nicheAbd.bac, nicheAbd.euk) %>%
  mutate(Kingdom=case_when(
    grepl("Euk",Kingdom)~"Eukaryotes",
    TRUE~"Bacteria")) %>%
  group_by(niche, Kingdom, sample, condition) %>%
  summarize(sum=sum(abundance)) %>%
  group_by(niche, Kingdom, condition) %>%
  summarize(mean=mean(sum)) %>%
  ungroup()

nicheSummary %>%
  filter(niche!="unclear") %>%
  mutate(niche=factor(niche, levels=c(
    "seawater-unique","seawater+ice","ice-unique","ice-exported"))) %>%
  ggplot(aes(x = condition, y = mean, fill = niche)) +
  geom_bar(stat = "identity", position="stack") +
  labs(y = "Relative abundance") +
  scale_fill_manual(values = c(
    "seawater-unique" = "lightskyblue4", 
    "seawater+ice"="cadetblue2",
    "ice-exported" = "maroon4",
    "ice-unique"="palegoldenrod")) +  
  theme_classic() +
  facet_grid(Kingdom~.) +
  theme(
    axis.text.x = element_text(
      angle = 0, vjust=1, hjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank())

###############################

# Remove temporary data
rm(core1,core2,core3,core4,core5,core6,core7,core8)

