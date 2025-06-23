#############################################################
 ## SHARED AND UNIQUE ASVs -- BAC
#############################################################

# OrgSW + Ctr
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

# OrgIce
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

# IceMelt
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
core.bac[["OrgSW+Ctr"]] <- as.character(row.names(core1))
core.bac[["OrgIce"]] <- as.character(row.names(core2))
core.bac[["IceMelt"]] <- as.character(row.names(core3))

# Reformat
overlaps.bac <- presAbs(core.bac) %>%
  rownames_to_column("asv")

## Figure S2
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
nicheCat.bac <- overlaps.bac %>%
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

# Count 
nicheCat.bac %>%
  summarise(
    n_OrgIce = sum(OrgIce == 1 & SW == 0),
    shared_OrgIce_IceMelt = sum(OrgIce == 1 & IceMelt == 1 & SW == 0),
    perc_OrgIce_IceMelt = ifelse(
    n_OrgIce > 0,shared_OrgIce_IceMelt / n_OrgIce * 100,NA_real_))

# Add abundances & metadata
nicheAbd.bac <- ASV.bac.rel %>%
  rownames_to_column(var = "asv") %>%
  pivot_longer(cols = -asv, names_to="sample", values_to="abundance") %>%
  right_join(nicheCat.bac) %>%
  left_join(ENV.bac, by =c("sample"="clip_id")) %>%
  drop_na(abundance) %>%
  mutate(type="Bacteria")

# Most common niche per Genus; ASVs per niche
nichePref.bac <- nicheCat.bac %>%
  left_join(topBac) %>%
  drop_na(mean) 


#############################################################
## SHARED AND UNIQUE ASVs -- EUK
#############################################################

# OrgSW + Ctr
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

# OrgIce
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

# IceMelt
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
core.euk[["OrgSW+Ctr"]] <- as.character(row.names(core4))
core.euk[["OrgIce"]] <- as.character(row.names(core5))
core.euk[["IceMelt"]] <- as.character(row.names(core6))

# Reformat
overlaps.euk <- presAbs(core.euk) %>%
  rownames_to_column("asv")

## Figure S2
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
nicheCat.euk <- overlaps.euk %>%
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

# Count 
nicheCat.euk %>%
  summarise(
    n_OrgIce = sum(OrgIce == 1 & SW == 0),
    shared_OrgIce_IceMelt = sum(OrgIce == 1 & IceMelt == 1 & SW == 0),
    perc_OrgIce_IceMelt = ifelse(
    n_OrgIce > 0, shared_OrgIce_IceMelt / n_OrgIce * 100,NA_real_))

# Add abundances & metadata
nicheAbd.euk <- ASV.euk.rel %>%
  rownames_to_column(var = "asv") %>%
  pivot_longer(cols = -asv, names_to="sample", values_to="abundance") %>%
  right_join(nicheCat.euk) %>%
  left_join(ENV.euk, by =c("sample"="clip_id")) %>%
  drop_na(abundance) %>%
  mutate(type="Eukaryotes")

# Preferred niche per Genus; based on most abundant ASV
nichePref.euk <- nicheCat.euk %>%
  left_join(topEuk) %>%
  drop_na(mean)


#############################################################
 ## COMBINE + PLOT
#############################################################

## SEEEDED FRACTIONS -- FIGURE 3

## Contribution of each fraction
nicheSummary <- rbind(nicheAbd.bac, nicheAbd.euk) %>%
  mutate(Kingdom=case_when(
    grepl("Euk",Kingdom)~"Eukaryotes",
    TRUE~"Bacteria")) %>%
  group_by(niche, Kingdom, sample, condition) %>%
  summarize(sum=sum(abundance)) %>%
  group_by(niche, Kingdom, condition) %>%
  summarize(mean=mean(sum)) %>%
  ungroup()

## Ice-exported over time
nicheTime <- rbind(nicheAbd.bac, nicheAbd.euk) %>%
  mutate(Kingdom=case_when(
    grepl("Euk",Kingdom)~"Eukaryotes", TRUE~"Bacteria")) %>%
  filter(niche=="ice-exported", time %in% c("0h", "17h", "49h")) %>%
  group_by(niche, asv, Kingdom, sample, condition, time) %>%
  summarize(sum=sum(abundance)) %>%
  group_by(niche,asv, Kingdom, condition, time) %>%
  summarize(mean=mean(sum)) %>%
  ungroup() %>%
  filter(mean > 0.01) %>%
  group_by(Kingdom, time) %>%
  summarise(n_ASVs = n_distinct(asv)) %>%
  arrange(Kingdom, time)

## FIGURE 3b
# Export size 3x5
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

# Show all replicates
nicheAbd.bac %>%
  filter(
    condition %in% c("IceMelt", "Ctr"),
    niche != "unclear",
    time %in% c("17h", "49h")) %>%
  group_by(condition, time, replicate) %>%
  mutate(rel_abd = abundance / sum(abundance) * 100) %>%
  ungroup() %>%
  group_by(niche, condition, time, replicate) %>%
  summarize(rel_abd = sum(rel_abd), .groups = "drop") %>%
  mutate(label = paste0("Rep", replicate, "_", time)) %>%
  mutate(label = factor(label, levels = unique(label))) %>%
  ggplot(aes(x = label, y = rel_abd, fill = niche)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Relative abundance (%)") +
  scale_fill_manual(values = c(
    "seawater-unique"="lightskyblue4", 
    "seawater+ice"="cadetblue2",
    "ice-exported"= "maroon4",
    "ice-unique"="palegoldenrod")) +
  facet_grid(condition ~ .) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank())

nicheAbd.euk %>%
  filter(
    condition %in% c("IceMelt", "Ctr"),
    niche != "unclear",
    time %in% c("17h", "49h")) %>%
  group_by(condition, time, replicate) %>%
  mutate(rel_abd = abundance / sum(abundance) * 100) %>%
  ungroup() %>%
  group_by(niche, condition, time, replicate) %>%
  summarize(rel_abd = sum(rel_abd), .groups = "drop") %>%
  mutate(label = paste0("Rep", replicate, "_", time)) %>%
  mutate(label = factor(label, levels = unique(label))) %>%
  ggplot(aes(x = label, y = rel_abd, fill = niche)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Relative abundance (%)") +
  scale_fill_manual(values = c(
    "seawater-unique"="lightskyblue4", 
    "seawater+ice"="cadetblue2",
    "ice-exported"= "maroon4",
    "ice-unique"="palegoldenrod")) +
  facet_grid(condition ~ .) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank())

######################################

## TAXA-SPECIFIC DYNAMICS
## FIGURE 4

# Keep ASVs in >=2 replicates at least once in 17h/49h
asv.keep <- rbind(nicheAbd.bac, nicheAbd.euk) %>%
  filter(niche == "ice-exported", time %in% c("0h","17h","49h")) %>%
  group_by(asv, condition, time, replicate, Genus, type, niche) %>%
  summarize(sum_abd = sum(abundance, na.rm=T), .groups="drop") %>%
  group_by(asv, time) %>%
  summarize(n_reps = sum(sum_abd > 0.1), .groups="drop") %>%
  mutate(pass = case_when(
    time %in% c("17h", "49h") & n_reps >= 2 ~"true",
    time == "0h" & n_reps >= 1 ~"T0")) %>%
  filter(!is.na(pass)) %>%
  group_by(asv) %>%
  summarize(n = n_distinct(time), .groups="drop") %>%
  filter(n >= 2)

## FIGURE 4
# Keep genera with mean >0.5 in min. 1 condition
# further filtering after initial plotting
# e.g. remove some UC taxa (especially if related, known genera present)
# Focus on taxa with consistent detection
# Ensure that IceMelt abundance is greater than Ctr
# export size 5x6
rbind(nicheAbd.bac, nicheAbd.euk) %>%
  filter(niche == "ice-exported", time %in% c("0h","17h","49h")) %>%
  semi_join(asv.keep, by="asv") %>%
  group_by(condition, time, replicate, Genus, type) %>%
  summarize(genus_abd = sum(abundance, na.rm=T), .groups="drop") %>%
  group_by(condition, time, Genus, type) %>%
  summarize(mean = mean(genus_abd), .groups="drop") %>%
  ungroup() %>%
  filter(
    (condition=="IceMelt" & mean > max(mean[condition=="Ctr"])) |  
    condition=="OrgIce") %>%
  filter(condition %in% c("OrgIce","IceMelt")) %>%
  filter(!Genus %in% c(
    "Cellvibrionaceae uc","Protaspa-lineage uc",
    "Cercozoa uc","Chaetoceros","Pedinellales uc",
    "Pyramimonadales uc","Cryomonadida uc",
    "Cryothecomonas-lineage uc","Suessiaceae uc",
    "Rhodobacteraceae uc","Gymnodiniaceae uc",
    "Filosa-Thecofilosea uc","Mataza-lineage uc")) %>%
  ggplot(aes(x = time, y = reorder(Genus, desc(-mean)), fill=mean)) +  # reorder Genus by mean
  geom_tile(color="white") +
  facet_wrap(type ~., scales="free", ncol=2) +
  scale_fill_gradientn(
    colors = c("whitesmoke", "lemonchiffon1", "lemonchiffon2", 
               "lightsteelblue1", "lightsteelblue3", "deeppink4"),
    trans = "log10",
    na.value = "white",
    breaks = c(0, 0.01, 0.05, 1, 5, 12)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title = element_blank(),
    axis.ticks = element_blank())

# SUPPLEMENT -- Figure S3
# Plot all replicates
# Export size 5x8
rbind(nicheAbd.bac, nicheAbd.euk) %>%
  filter(niche=="ice-exported", time %in% c("0h","17h","49h")) %>%
  semi_join(asv.keep, by="asv") %>%
  group_by(asv, condition, time, replicate, Genus, type) %>%
  summarize(sum_abd = sum(abundance, na.rm=T), .groups="drop") %>%
  group_by(condition, time, replicate, Genus, type) %>%
  summarize(abundance = sum(sum_abd, na.rm=T), .groups="drop") %>%
  group_by(Genus) %>%
  mutate(mean = mean(abundance, na.rm=T)) %>%  # Calculate mean per Genus
  ungroup() %>%
  mutate(
    Genus = forcats::fct_reorder(Genus, mean, .desc=F),
    label = paste0("Rep", replicate, "_", time),
    label = factor(label, levels = unique(label))) %>%
  filter(
    (condition=="IceMelt" & abundance > max(abundance[condition == "Ctr"])) |  # IceMelt strictly greater than Ctr
    condition=="OrgIce") %>%
  filter(!Genus %in% c(
    #"JGI 0000069-P22 uc","Saprospiraceae uc",
    #"MAST-1C uc","Raphid-pennate uc","Botuliformida uc",
    # "Apocalathium", "Amphidiniopsis", #"TAGIRI1-lineage uc",
    "Cellvibrionaceae uc","Protaspa-lineage uc",
    "Cercozoa uc","Pedinellales uc","Chaetoceros",
    "Pyramimonadales uc","Cryomonadida uc",
    "Cryothecomonas-lineage uc","Suessiaceae uc",
    "Rhodobacteraceae uc","Gymnodiniaceae uc",
    "Filosa-Thecofilosea uc","Mataza-lineage uc")) %>%
  ggplot(aes(x = label, y = Genus, fill = abundance)) +
  geom_tile(color = "white") +
  facet_wrap(type ~., scales = "free", ncol=2) +
  scale_fill_gradientn(
    colors=c("whitesmoke","lemonchiffon1","lemonchiffon2",
             "lightsteelblue1","lightsteelblue3","deeppink4"),
    trans = "log10",
    na.value="white",
    breaks = c(0, 0.01, 0.05, 1, 5, 12)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

###############################

# Remove temporary data
rm(core1,core2,core3,core4,core5,core6,core7,core8)

