
#############################################################
  ## MICROBIAL SOURCE TRACKING 
#############################################################

## BACTERIA
meta <- ENV.bac %>%
  filter(sample_title!="NK" & time %in% c("0h","17h","49h")) %>%
  dplyr::select(c("clip_id","condition")) %>%
  #dplyr::rename("Env" = "condition") %>%
  mutate(SourceSink = case_when(
    condition == "IceMelt" ~ "Sink",
    condition == "Ctr" ~ "Sink",
    TRUE ~ "Source")) %>%
  mutate(Env = case_when(
    grepl("Ctr", condition)~"Ctr",
    grepl("OrgSW", condition)~"OrgSW",
    grepl("OrgIce", condition)~"OrgIce",
    TRUE ~ "IceMelt")) %>%
  group_by(SourceSink) %>%
  mutate(id = as.character(row_number())) %>%
  ungroup() %>%
  column_to_rownames("clip_id") 

# Reformat & match ASVs
asv <- ASV.bac %>% 
  dplyr::select(any_of(rownames(meta))) %>%
  #filter(rowMeans(across(everything(), ~ . >= 0.001)) >= 0.001) %>%
  as.matrix() %>% t()

# Reformat & match metadata
meta <- meta[match(rownames(asv), rownames(meta)), ]

# Calculate
sourceBac <- FEAST(
  asv, meta, EM_iterations = 1000,
  different_sources_flag = 0, 
  #dir_path="C:/Users/mwietz/ownCloud - mwietz@owncloud.mpi-bremen.de/AWI_MPI/collaborations/Beluga/Rstats",
  dir_path="~/",
  outfile="sourceSink_bac")

# Extract results
sourceBacRes <- do.call(rbind, lapply(names(sourceBac), function(name) {
  values <- sourceBac[[name]]
  if (length(values) == 22) {data.frame(
    source = name, t(as.data.frame(values)))} else {NULL}}))

# Rename sinks by original names
id <- meta %>%
  rownames_to_column("sink") %>%
  filter(Env %in% c("IceMelt","Ctr"))
names(sourceBacRes)[2:23] <- id$sink

# Re-add metadata
#sourceBacRes2 <- sourceBacRes %>%
  #filter(grepl("_", source)) %>%
 # mutate(clip_id = gsub("_orgSeawater|_orgCore*","", source)) 

# Average by sample type
sourceBacAvg <- sourceBacRes %>%
  reshape2::melt() %>%
  left_join(ENV.bac, by =c("variable"="clip_id")) %>%
  #group_by(source, condition) %>%
  #summarize(sum = sum(value, na.rm=T), .groups="drop") %>%
  #mutate(source=gsub(".*org", "org", source)) %>%
  group_by(source, condition) %>%
  summarize(mean = mean(value, na.rm=T), .groups="drop") %>%
  ungroup() %>%
  group_by(condition) %>%
  mutate(mean = ifelse(grepl("OrgIce", source), mean(mean[grepl("OrgIce", source)]), mean)) %>%
  ungroup() %>%
  distinct(mean, .keep_all=T) %>%
  mutate(type="Bacteria")

###################################################

## EUKARYOTES
meta <- ENV.euk %>%
  filter(sample_title!="NK" & time %in% c("0h","17h","49h")) %>%
  dplyr::select(c("clip_id","condition")) %>%
  #dplyr::rename("Env" = "condition") %>%
  mutate(SourceSink = case_when(
    condition == "IceMelt" ~ "Sink",
    condition == "Ctr" ~ "Sink",
    TRUE ~ "Source")) %>%
  mutate(Env = case_when(
    grepl("Ctr", condition)~"Ctr",
    grepl("OrgSW", condition)~"OrgSW",
    grepl("OrgIce", condition)~"OrgIce",
    TRUE ~ "IceMelt")) %>%
  group_by(SourceSink) %>%
  mutate(id = as.character(row_number())) %>%
  ungroup() %>%
  column_to_rownames("clip_id") 

# Reformat & match ASVs
asv <- ASV.euk %>% 
  dplyr::select(any_of(rownames(meta))) %>%
  #filter(rowMeans(across(everything(), ~ . >= 0.001)) >= 0.001) %>%
  as.matrix() %>% t()

# Reformat & match metadata
meta <- meta[match(rownames(asv), rownames(meta)), ]

# Calculate
sourceEuk <- FEAST(
  asv, meta, EM_iterations = 1000,
  different_sources_flag = 0, 
  #dir_path="C:/Users/mwietz/ownCloud - mwietz@owncloud.mpi-bremen.de/AWI_MPI/collaborations/Beluga/Rstats",
  dir_path="~/",
  outfile="sourceSink_euk")

# Extract results
sourceEukRes <- do.call(rbind, lapply(names(sourceEuk), function(name) {
  values <- sourceEuk[[name]]
  if (length(values) == 22) {data.frame(
    source = name, t(as.data.frame(values)))} else {NULL}}))

# Rename sinks by original names
id <- meta %>%
  rownames_to_column("sink") %>%
  filter(Env %in% c("IceMelt","Ctr"))
names(sourceEukRes)[2:23] <- id$sink

# Re-add metadata
#sourceBacRes2 <- sourceBacRes %>%
#filter(grepl("_", source)) %>%
# mutate(clip_id = gsub("_orgSeawater|_orgCore*","", source)) 

# Average by sample type
sourceEukAvg <- sourceEukRes %>%
  reshape2::melt() %>%
  left_join(ENV.euk, by =c("variable"="clip_id")) %>%
  #group_by(source, condition) %>%
  #summarize(sum = sum(value, na.rm=T), .groups="drop") %>%
  #mutate(source=gsub(".*org", "org", source)) %>%
  group_by(source, condition) %>%
  summarize(mean = mean(value, na.rm=T), .groups="drop") %>%
  ungroup() %>%
  group_by(condition) %>%
  mutate(mean = ifelse(grepl("OrgIce", source), mean(mean[grepl("OrgIce", source)]), mean)) %>%
  ungroup() %>%
  distinct(mean, .keep_all=T) %>%
  mutate(type="Eukaryotes")

###################################################

## MERGE AND PLOT
sourceSink <- rbind(sourceBacAvg, sourceEukAvg)  %>%
mutate(source=gsub(".*Org", "Org", source)) 
#%>%
  #mutate(type=factor(type, levels=c("SML","Seawater","Snow")))

# Heatmap
sourceSink %>%
  filter(source!="Unknown") %>%
ggplot() + 
  aes(x = source, y = condition, fill = mean) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colors = c("seashell","mistyrose4","lightsteelblue4"),
    values = scales::rescale(c(0, 0.005, 0.01, 0.05, 0.1))) +
  geom_text(aes(
    label = ifelse(is.na(mean), "NA", sprintf("%.1f", mean * 100)),
    color = ifelse(mean < 0.12 | is.na(mean),"black","white")),
    size = 4) +
  scale_color_identity() +
  theme_minimal() +
  facet_grid(type~., scales="free_x") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text = element_text(size=12),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12)) 

###################################################

# remove temp-data
rm(asv, meta, id)
