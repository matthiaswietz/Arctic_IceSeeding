#############################################################
 ## ENRICHMENT -- PROK
#############################################################

# Prepare metadata
# remove OrgSW: not fractionated 
meta <- ENV.bac %>%
  filter(!type %in% c("NK","OrgSW")) %>%
  remove_rownames %>%
  column_to_rownames("clip_id")

# Prepare ASVs
asv <- ASV.bac %>% 
  filter(rowMeans(across(everything(), ~ . >= 0.1)) >= 0.1) %>% 
  dplyr::select(which(names(.) %in% row.names(meta)))
tax <- TAX.bac[row.names(asv),] 

# Format tables
asv = otu_table(asv, taxa_are_rows=T)
tax = tax_table(tax)

# Rename 
rownames(tax) <- rownames(asv)
colnames(tax) <- c(
  "Kingdom","Phylum","Class","Order",
  "Family","Genus","Species")

# Create phyloseq object
pseq <- phyloseq(
  otu_table(asv, taxa_are_rows=F), 
  sample_data(meta), 
  tax_table(tax)) 

# Calculate enrichment
FracBac <- ancombc2(
  pseq,
  assay_name = "counts",
  tax_level = "Genus",
  fix_formula = "fraction",
  p_adj_method = "holm",
  pseudo = 0,
  pseudo_sens = T,
  prv_cut = 0.1,
  lib_cut = 0,
  s0_perc = 0.05,
  group = NULL,
  struc_zero=F,
  neg_lb = F ,
  alpha = 0.05,
  n_cl = 1,
  iter_control = list(tol = 0.01, max_iter=20, verbose=F),
  em_control = list(tol = 1e-05, max_iter=100),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(fwer_ctrl_method="holm", B=100),
  trend_control = list(contrast=NULL, node=NULL, solver="ECOS", B=100))

# Extract results
# Join the niche of the most abundant ASV per Genus
FracBacRes = FracBac$res %>%
  mutate(fraction = case_when(
    lfc_fractiongreater3um > 1 ~ "greater 3um",
    lfc_fractiongreater3um < -1 ~ "0.2-3um",
    diff_fractiongreater3um=="FALSE"~"both")) %>%
  mutate(enrichFrac = case_when(
    #fraction %in% c("FL","PA") ~ (abs(lfc_seasonSpring) + abs(lfc_fractionPA)) / 2,
    #fraction == "FL+PA" ~ abs(lfc_seasonSpring),
    fraction %in% c("0.2-3um","greater 3um","both") ~ abs(lfc_fractiongreater3um),
    TRUE ~ NA_real_)) %>%
  right_join(nicheGenus.bac, by=c("taxon"="Genus")) %>%
  mutate(niche = case_when(
    is.na(niche) ~"unclear", T~niche),
    tax="Bacteria") %>%
  drop_na(fraction) 


#############################################################
 ## ENRICHMENT -- EUK
#############################################################

# Prepare metadata
# remove OrgSW: not fractionated 
meta <- ENV.euk %>%
  filter(!type %in% c("NK","OrgSW")) %>%
  remove_rownames %>%
  column_to_rownames("clip_id")

# Prepare ASVs
asv <- ASV.euk %>% 
  filter(rowMeans(across(everything(), ~ . >= 0.1)) >= 0.1) %>% 
  dplyr::select(which(names(.) %in% row.names(meta)))
tax <- TAX.euk[row.names(asv),] 

# Format tables
asv = otu_table(asv, taxa_are_rows=T)
tax = tax_table(tax)

# Rename 
rownames(tax) <- rownames(asv)
colnames(tax) <- c(
  "Kingdom","Phylum","Class","Order",
  "Family","Genus","Species")

# Create phyloseq object
pseq <- phyloseq(
  otu_table(asv, taxa_are_rows=F), 
  sample_data(meta), 
  tax_table(tax)) 

# Calculate enrichment
FracEuk <- ancombc2(
  pseq,
  assay_name = "counts",
  tax_level = "Genus",
  fix_formula = "fraction",
  p_adj_method = "holm",
  pseudo = 0,
  pseudo_sens = T,
  prv_cut = 0.1,
  lib_cut = 0,
  s0_perc = 0.05,
  group = NULL,
  struc_zero=F,
  neg_lb = F ,
  alpha = 0.05,
  n_cl = 1,
  iter_control = list(tol = 0.01, max_iter=20, verbose=F),
  em_control = list(tol = 1e-05, max_iter=100),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(fwer_ctrl_method="holm", B=100),
  trend_control = list(contrast=NULL, node=NULL, solver="ECOS", B=100))

# Extract results
# Join the niche of the most abundant ASV per Genus
FracEukRes = FracEuk$res %>%
  mutate(fraction = case_when(
    lfc_fractiongreater3um > 1 ~ "greater 3um",
    lfc_fractiongreater3um < -1 ~ "0.2-3um",
    diff_fractiongreater3um=="FALSE"~"both")) %>%
  mutate(enrichFrac = case_when(
    fraction %in% c("0.2-3um","both","greater 3um") ~ abs(lfc_fractiongreater3um),
    TRUE ~ NA_real_)) %>%
  drop_na() %>%
  right_join(nicheGenus.euk, by=c("taxon"="Genus")) %>%
  mutate(niche = case_when(
    is.na(niche) ~"unclear", T~niche),
    tax="Eukaryotes") %>%
  drop_na(fraction) 


#########################################################
 ## COMBINE + PLOT
#########################################################

# Relative contribution of fractions per niche
rbind(FracBacRes, FracEukRes) %>%
  group_by(fraction, niche, tax) %>%
  summarise(nicheCount = n(), .groups = "drop") %>%
  group_by(fraction, tax) %>%
  mutate(relContri = nicheCount / sum(nicheCount)) %>% 
  filter(niche!="unclear") %>%
  mutate(niche=factor(niche, levels=c(
    "seawater-unique","seawater+ice","ice-unique",
    "ice-exported", "unclear"))) %>%
  mutate(fraction=factor(fraction, levels=c(
    "0.2-3um","greater 3um","both"))) %>%
  ggplot(aes(x = fraction, y = relContri, fill = niche)) +
  geom_bar(stat = "identity", position = "stack") + 
  facet_grid(tax~.) +
  scale_fill_manual(values = c(
    "seawater-unique" = "lightskyblue4", 
    "seawater+ice"="cadetblue2",
    "ice-exported" = "maroon4",
    "ice-unique"="palegoldenrod")) +  
  theme_classic() +
  ylab("Relative fraction") +
  theme(
    axis.text.x = element_text(
      angle = 34, hjust = 1, size=10),
    axis.title.x = element_blank(),
    axis.ticks.x=element_blank()) 

######################################

# Major taxa per fraction & niche
# omit some unclassifieds
# omit "double taxa" (eg keep Pelagomonas; remove uncult. clades)
# Export size 6x8
rbind(FracEukRes, FracBacRes) %>%
  filter(!taxon %in% c(
    "Hemidiscaceae uc", "Gymnodiniaceae uc",
    "Cand Puniceispirillum","MB11C04 marine group",
    "TAGIRI1-lineage uc","Botuliformidae uc",
    "Pyramimonadales uc","Prymnesiophyceae_Clade_B3 uc",
    "Stramenopiles ucX","Prorocentrum","Ochrophyta uc",
    "Dolichomastigaceae-B","Pyramimonadales uc",
    "Pelagomonadales_clade_B1","Pelagomonadales_clade_B2",
    "Pelagomonadales_clade_A uc","Raphid-pennate uc",
    "Ventricleftida uc","NPK2-lineage uc")) %>%
  filter(!is.na(enrichFrac) & fraction != "both" & niche != "unclear") %>%
  filter(p_fractiongreater3um < 0.05 & enrichFrac > 1.5) %>%
  mutate(taxon = gsub("uc","", gsub("Marinimicrobia \\(SAR406 clade\\)","SAR406", taxon))) %>%
  #mutate(taxon = gsub("uc", "", taxon)) %>%
  mutate(
    fraction = factor(fraction, levels = c("0.2-3um","greater 3um")),  # Ordering fraction
    taxon = factor(taxon, levels = taxon[order(niche, enrichFrac)])) %>%
  ggplot(aes(x = enrichFrac, y = taxon)) +
  # Lollipop plot with thin lines
  geom_segment(aes(xend=0, yend=taxon), size=0.2, color="gray89") +  # Lollipop lines
  geom_point(aes(fill=niche, shape=Kingdom), color="black", size=3) + 
  scale_shape_manual(values = c(21L, 22L)) +
  # Place taxon names to the left of the dots with better alignment
  geom_text(aes(label=taxon), size=3, hjust=1, vjust=0.4, color="black", nudge_x = -0.1) +  # Labels inside
  scale_x_continuous(
    name = "Enrichment Factor",
    expand = c(0,0),
    limits=c(0,4)) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),  
    axis.ticks = element_blank(),
    strip.text.y = element_blank(),  # Remove y-axis facet label
    strip.text.x = element_text(hjust = 0.5),  # Center fraction label
    strip.placement = "outside",  # Place fraction label outside plot area
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "right") +
  scale_fill_manual(values = c(
    "seawater-unique" = "lightskyblue4", 
    "seawater+ice"="cadetblue2",
    "ice-exported" = "maroon4",
    "ice-unique"="lightgoldenrod2")) +  
  facet_wrap(vars(fraction), ncol = 2, scales = "free_y") +  
  theme(strip.text.y = element_text(angle = 90))

######################################

# remove temp-data
rm(asv, tax, meta, pseq)
