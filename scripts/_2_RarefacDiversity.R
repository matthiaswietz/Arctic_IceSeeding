
#######################################################
 ### ALPHA-DIVERSITY ### 
#######################################################

# Calculated on full set (i.e. without subsetting)
# Run on AWI-VM for speed
# Extract from iNEXT objects -- continued here  
load("iNEXT_16S.Rdata")
load("iNEXT_18S.Rdata")

###################################

## RAREFACTION ##
rarefac.bac <- fortify(iNEXT.bac, type=1) %>%
  left_join(ENV.bac, by=c("Assemblage"="clip_id")) %>%
  drop_na(condition)
rarefac.line.bac <- rarefac.bac %>%
  filter(Method != "Observed") %>%
  mutate(Method = factor(Method, levels=c(
    "Rarefaction","Extrapolation")))

rarefac.euk <- fortify(iNEXT.euk, type=1) %>%
  left_join(ENV.euk, by=c("Assemblage"="clip_id")) %>%
  drop_na(condition)
rarefac.line.euk <- rarefac.euk %>%
  filter(Method != "Observed") %>%
  mutate(Method = factor(Method, levels=c(
    "Rarefaction","Extrapolation")))

###################################

## COVERAGE ##

cover.bac <- fortify(iNEXT.bac, type=2) %>%
  left_join(ENV.bac, by=c("Assemblage"="clip_id")) %>%
  drop_na(condition)
cover.line.bac <- cover.bac %>%
  filter(Method != "Observed") %>%
  mutate(Method = factor(Method, levels=c(
    "Rarefaction","Extrapolation")))

cover.euk <- fortify(iNEXT.euk, type=2) %>%
  left_join(ENV.euk, by=c("Assemblage"="clip_id")) %>%
  drop_na(condition)
cover.line.euk <- cover.euk %>%
  filter(Method != "Observed") %>%
  mutate(Method = factor(Method, levels=c(
    "Rarefaction","Extrapolation")))

###################################

## COMBINE + PLOT

coverage <- rbind(cover.bac,cover.euk)
cover.line <- rbind(cover.line.bac,cover.line.euk)
rarefaction <- rbind(rarefac.bac,rarefac.euk)
rarefac.line <- rbind(rarefac.line.bac,rarefac.line.euk)

ggplot(coverage, aes(
  x=x, y=y, colour=Assemblage))+ 
geom_line(
  aes(linetype=Method), 
  lwd = 0.5, data=cover.line) +
scale_colour_discrete(guide="none") +
scale_x_continuous(
  limits = c(0,1e+5)) +
scale_y_continuous(
  breaks = seq(0.9,1,0.05), 
  limits = c(0.9,1)) +
facet_grid(locus_tag~condition) +
labs(x="Sample size", y="Sample coverage") +
theme_bw(base_size = 12) + 
theme(legend.position="bottom")

ggplot(rarefaction, aes(
  x=x, y=y, colour=Assemblage)) +
geom_line(aes(
  linetype=Method), 
  lwd = 0.5, data=rarefac.line) +
scale_colour_discrete(guide="none") +
scale_x_continuous(limits = c(0,1e+5)) +
facet_grid(locus_tag~condition) +
labs(x="Sample size", y="Species richness") +
theme_bw() + 
theme(legend.position="bottom")


#######################################################
  ### REFORMATTING
#######################################################

# Extract metrices
richness <- iNEXT.bac$AsyEst[
  iNEXT.bac$AsyEst$Diversity=="Species richness",] %>%
  arrange(Assemblage) 
simpson <- iNEXT.bac$AsyEst[
  iNEXT.bac$AsyEst$Diversity=="Simpson diversity",] %>%
  arrange(Assemblage) 
shannon <- iNEXT.bac$AsyEst[
  iNEXT.bac$AsyEst$Diversity=="Shannon diversity",] %>%
  arrange(Assemblage) 

# compile; calculate evenness
div.bac <- data.frame(
  clip_id = richness$Assemblage,
  richness = richness$Observed,
  simpson = simpson$Observed,
  shannon = shannon$Observed) %>%
  mutate(evenness = shannon/log(richness))

###################################

richness <- iNEXT.euk$AsyEst[
  iNEXT.euk$AsyEst$Diversity=="Species richness",] %>%
  arrange(Assemblage) 
simpson <- iNEXT.euk$AsyEst[
  iNEXT.euk$AsyEst$Diversity=="Simpson diversity",] %>%
  arrange(Assemblage) 
shannon <- iNEXT.euk$AsyEst[
  iNEXT.euk$AsyEst$Diversity=="Shannon diversity",] %>%
  arrange(Assemblage) 

div.euk <- data.frame(
  clip_id = richness$Assemblage,
  richness = richness$Observed,
  simpson = simpson$Observed,
  shannon = shannon$Observed) %>%
  mutate(evenness = shannon/log(richness))

###################################

## MERGE ALL
AlphaDiv <- rbind(
  div.bac, div.euk) %>%
  left_join(ENV) %>%
  filter(!grepl("NK", sample_title))

# Remove temporary objects
rm(richness, simpson, shannon,
   div.bac, div.euk, cover.bac, cover.euk,
   cover.line.bac, cover.line.euk,
   rarefac.line.bac, rarefac.line.euk,
   rarefac.bac, rarefac.euk)
