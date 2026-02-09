#############################################################
 ## DOC
#############################################################

# Load measurements
DOC <- read.csv(
  "DOC.txt", h=T, sep="\t", check.names=F,
  stringsAsFactors=F, skipNul=T) 

## FIGURE 5
# Plot each replicate over time
# Use averaged T0 value
# Export size
DOC %>%
  mutate(time = as.character(time)) %>%
  mutate(coreid=paste0("Replicate ","",coreid))  %>%
  group_by(coreid, condition, time) %>%
  summarize(
    avg_NPOC = mean(`NPOC_umol_L-1`, na.rm=T),
    sd_NPOC  = sd(`NPOC_umol_L-1`, na.rm=T),
    .groups = "drop") %>%
  {
  # Compute shared 0h Ctr average
    t0_shared <- filter(., time=="0h", condition=="Ctr") %>%
      summarize(
        time = "0h", condition = "Ctr",
        avg_NPOC = mean(avg_NPOC, na.rm=T),
        sd_NPOC  = mean(sd_NPOC, na.rm=T))
  # Remove original 0h-Ctr rows and replace with shared across all coreids
    bind_rows(
      filter(., !(time=="0h" & condition=="Ctr")),
      crossing(coreid = unique(.$coreid), t0_shared))
  } %>%
  ggplot(aes(
    x = time, y = avg_NPOC, 
    group = condition, color = condition)) +
  geom_line() +
  geom_errorbar(aes(
    ymin = avg_NPOC - sd_NPOC, 
    ymax = avg_NPOC + sd_NPOC), width = 0.2) +
  geom_point(size = 2) +
  facet_grid(coreid ~ .) +
  labs(
    y = "DOC (µmol L-1)", color = "Condition") +
  scale_y_continuous(n.breaks = 3.4) +
  scale_color_manual(values = c(
    "Ctr" = "darkmagenta", 
    "IceMelt" = "orchid3")) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(linewidth = 0.3),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank())

#############################################################
 ## DOM / FT-ICR-MS
#############################################################

## FIRST: download full FT data from PANGAEA
# Open in browser: https://doi.pangaea.de/10.1594/PANGAEA.978551?format=textfile
# Rename the downloaded file to dataset.tab

# Prepare import
lines <- readLines("dataset.tab")
end_comment <- grep("\\*/", lines)[1]

# Load data
D<- read.delim(
  text = lines[(end_comment + 1):length(lines)],
  header = TRUE,
  check.names = FALSE)

# Load key table
Dkey<-read.csv("DOC.txt",sep="\t")

# filter out all isotopologues (otherwise some formulae remain multiple times)
D2<-D1[
  D1$`13C atoms [#]`==0 & 
  D1$`15N atoms [#]`==0 &
  D1$`34S atoms [#]`==0 & 
  D1$`18O atoms [#]`==0,]

# filter more, rename columns, give IDs
ID <- data.frame(
  sample=Dkey$SampleID,
  originalname=Dkey$original_sample_name,
  type=Dkey$condition,
  time=Dkey$time,
  core=Dkey$coreid)

# Define samples
seawater=ID$sample[ID$type=="Ctr" & ID$time=="0h"]
ice=ID$sample[ID$type=="IceMelt" & ID$time=="0h"]
D2$type<-D2$`Exp condition`
D2$type[D2$type=="OrgIce"]<-"IceMelt"
D2$type[D2$type=="OrgSW"]<-"Ctr"

# New DF with specific samples at time 0h
IceVsSeaReference <- D2[
  D2$type %in% c("Ctr","IceMelt") & D2$`Time point`=="0h",]

# sulfur compounds in ice and seawater references
ggplot(
  IceVsSeaReference,
  aes(x=`O atoms [#]` / `C atoms [#]`, y=`H atoms [#]`/`C atoms [#]`, color=as.factor(`S atoms [#]`))) +
  geom_point() +
  facet_grid(`S atoms [#]`~type) + 
  ggtitle("h=0")+scale_color_discrete("S")

# get formulae unique to ice
allice<-IceVsSeaReference%>%group_by(type,Formula)%>%filter(n()==4 &type=="IceMelt")
atleastoneseawater<-IceVsSeaReference%>%group_by(type,Formula)%>%filter(type=="Ctr")
uniqueallice<-IceVsSeaReference[IceVsSeaReference$Formula%in%setdiff(allice$Formula,atleastoneseawater$Formula),]
uniqueallice<-uniqueallice[!is.na(uniqueallice$Int),]

# van Krevelen plot -- OrgIce
ggplot(uniqueallice,aes(y=`H atoms [#]`/`C atoms [#]`,x=`O atoms [#]` / `C atoms [#]`,color=as.factor(`S atoms [#]`)))+geom_point()+scale_color_viridis_d("S")+theme_classic()+xlab("O/C")+ylab("H/C")

# Filter OrgIce but never in Ctr 
stillinotherseawater<-D2$Formula[D2$type=="Ctr"]
uniqueallice2<-uniqueallice[!uniqueallice$Formula%in%unique(stillinotherseawater),]

length(unique(uniqueallice2$Formula))

D7<-D2%>%group_by(`Time point`,ID,type,`Sample ID`)%>%dplyr::summarise(uniqueN=sum(Formula %in% uniqueallice2$Formula))
D7$timec<-0
D7$timec[D7$`Time point`=="15h"]<-1
D7$timec[D7$`Time point`=="17h"]<-2
D7$timec[D7$`Time point`=="27h"]<-3
D7$timec[D7$`Time point`=="33h"]<-4
D7$timec[D7$`Time point`=="49h"]<-5

# Plot sum of unique ice over time, total sum declines rapidly 
ggplot(D7[D7$type=="IceMelt",],aes(x=timec,y=uniqueN,color=as.factor(ID)))+geom_point()+ylab("sum of unique ice compounds present")+geom_smooth(method="gam",formula=(y)~s(x,bs="tp",k=4),method.args=list(family="poisson"))

# Van krevelen plot of unique ice formulae
ggplot(uniqueallice2,aes(y=`H atoms [#]`/`C atoms [#]`,x=`O atoms [#]` / `C atoms [#]`,color=as.factor(`S atoms [#]`)))+geom_point()+scale_color_viridis_d("S")+theme_classic()+xlab("O/C")+ylab("H/C")

uniqueallice2$type="IceMelt unique"
complete<-rbind(IceVsSeaReference,uniqueallice2)
complete$type[complete$type=="IceMelt"]<-"IceMelt reference"
complete$type[complete$type=="Ctr"]<-"Ctr reference"
# more van krevelen plots
ggplot(complete,aes(y=`H atoms [#]`/`C atoms [#]`,x=`O atoms [#]` / `C atoms [#]`,color=as.factor(`S atoms [#]`)))+geom_point()+scale_color_viridis_d("S")+theme_classic()+xlab("O/C")+ylab("H/C")+facet_wrap(~type)

# Plot for manuscript -- FIGURE 7 
PVK<-ggplot(IceVsSeaReference[IceVsSeaReference$type=="IceMelt",],aes(y=`H atoms [#]`/`C atoms [#]`,x=`O atoms [#]` / `C atoms [#]`))+geom_point(color="lightgray")+geom_point(data=uniqueallice2,aes(y=`H atoms [#]`/`C atoms [#]`,x=`O atoms [#]`/`C atoms [#]`,color=as.factor(`S atoms [#]`)),inherit.aes = FALSE)+scale_color_viridis_d("S")+theme_classic()+xlab("O/C")+ylab("H/C")

Pmass<-ggplot(complete, aes(y=type,x=`m/z`)) + geom_boxplot()+theme_classic()+ylab("")+xlab("m/z")+xlim(100,1000)
PAI<-ggplot(complete, aes(y=type,x=`AI mod`)) + geom_boxplot()+theme_classic()+ylab("")+xlab(bquote(AI[mod]))+xlim(0,1)

plot1<-ggplot(complete, aes(fill=as.factor(`S atoms [#]`),x=type)) + geom_bar(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..]), position = position_dodge2(width = 0.9, preserve = "single"))+theme_classic()+xlab("")+ylab("")+scale_fill_manual("S",values=c("#440154","#fde725","orange"))+scale_y_continuous(labels=scales::percent_format())
plot2<-ggplot(complete, aes(fill=as.factor(`N atoms [#]`),x=type)) + geom_bar(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..]), position = position_dodge2(width = 0.9, preserve = "single"))+theme_classic()+xlab("")+ylab("% molecular formulae")+scale_fill_manual("N",values=c("#440154","#fde725","orange","red","violet"))+scale_y_continuous(labels=scales::percent_format())
plot3<-ggplot(complete, aes(fill=as.factor(`P atoms [#]`),x=type)) + geom_bar(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..]), position = position_dodge2(width = 0.9, preserve = "single"))+theme_classic()+xlab("")+ylab("")+scale_fill_manual("P",values=c("#440154","#fde725","orange"))+scale_y_continuous(labels=scales::percent_format())

PZ1<-plot_grid(plot1,plot2,plot3,nrow=3)
PZ2<-plot_grid(Pmass,PAI,nrow=2)

# compose parts
plot_grid(PZ1,plot_grid(PZ2,PVK,nrow=2,labels=c("b","c")),nrow=1,scale=c(1,1),rel_heights = c(1,1),labels=c("a",""))

## van krevelen plots for unique ice formulae colored by S,N,P and AI mod
icex<-IceVsSeaReference[IceVsSeaReference$type=="IceMelt",]
seawaterx<-IceVsSeaReference[IceVsSeaReference$type=="Ctr",]
percentageiceinsea<-length(unique(uniqueallice2$formula))/length(unique(icex$formula))*100
nice<-length(unique(icex$formula))
nsea<-length(unique(seawaterx$formula))

a1<-ggplot(uniqueallice2%>%arrange(`P atoms [#]`),aes(y=`H atoms [#]`/`C atoms [#]`,x=`O atoms [#]`/`C atoms [#]`,color=as.factor(`P atoms [#]`)))+geom_point()+facet_wrap(~type)+scale_color_viridis_d("P")
a2<-ggplot(uniqueallice2%>%arrange(`S atoms [#]`),aes(y=`H atoms [#]`/`C atoms [#]`,x=`O atoms [#]`/`C atoms [#]`,color=as.factor(`S atoms [#]`)))+geom_point()+facet_wrap(~type)+scale_color_viridis_d("S")
a3<-ggplot(uniqueallice2%>%arrange(`N atoms [#]`),aes(y=`H atoms [#]`/`C atoms [#]`,x=`O atoms [#]`/`C atoms [#]`,color=as.factor(`N atoms [#]`)))+geom_point()+facet_wrap(~type)+scale_color_viridis_d("N")
a4<-ggplot(uniqueallice2%>%arrange(`N atoms [#]`),aes(y=`H atoms [#]`/`C atoms [#]`,x=`O atoms [#]`/`C atoms [#]`,color=`AI mod`))+geom_point()+facet_wrap(~type)+scale_color_viridis_c("AImod")

plot_grid(a1,a2,a3,a4,nrow=2)

########################################

# Comparing Ice and Ctr
modAIv<-IceVsSeaReference$`AI mod`
modmz<-IceVsSeaReference$`m/z`
typeT<-IceVsSeaReference$type

# compare by mass and AI mod
coin::wilcox_test(modAIv~as.factor(typeT))
coin::wilcox_test(modmz~as.factor(typeT))

# compare by S,N,P
hasSice<-sum(IceVsSeaReference[!is.na(IceVsSeaReference$Int)&IceVsSeaReference$type=="IceMelt",]$`S atoms [#]` >0)
hasSseawater<-sum(IceVsSeaReference[!is.na(IceVsSeaReference$Int)&IceVsSeaReference$type=="Ctr",]$`S atoms [#]`>0)
hasSiceonly<-sum(uniqueallice2$S>0)
hasnoSice<-sum(IceVsSeaReference[!is.na(IceVsSeaReference$Int)&IceVsSeaReference$type=="IceMelt",]$`S atoms [#]`==0)
hasnoSseawater<-sum(IceVsSeaReference[!is.na(IceVsSeaReference$Int)&IceVsSeaReference$type=="Ctr",]$`S atoms [#]`==0)
hasnoSiceonly<-sum(uniqueallice2$`S atoms [#]`==0)
dfS<-data.frame(S=c(hasSice,hasSseawater),noS=c(hasnoSice,hasnoSseawater),row.names=c("IceMelt","Ctr"))
dfS1<-data.frame(S=c(hasSiceonly,hasSseawater),noS=c(hasnoSiceonly,hasnoSseawater),row.names=c("IceMelt unique","Ctr"))
dfS2<-data.frame(S=c(hasSice,hasSiceonly),noS=c(hasnoSice,hasnoSiceonly),row.names=c("IceMelt","IceMelt unique"))

mosaicplot(dfS, color = TRUE) 
mosaicplot(dfS1, color = TRUE) 
mosaicplot(dfS2, color = TRUE) 

fisher.test(dfS)$p.value
#fisher.test(dfS1)$p.value
#fisher.test(dfS2)$p.value

hasNice<-sum(IceVsSeaReference[!is.na(IceVsSeaReference$Int)&IceVsSeaReference$type=="IceMelt",]$`N atoms [#]`>0)
hasNseawater<-sum(IceVsSeaReference[!is.na(IceVsSeaReference$Int)&IceVsSeaReference$type=="Ctr",]$`N atoms [#]`>0)
hasNiceonly<-sum(uniqueallice2$N>0)
hasnoNice<-sum(IceVsSeaReference[!is.na(IceVsSeaReference$Int)&IceVsSeaReference$type=="IceMelt",]$`N atoms [#]`==0)
hasnoNseawater<-sum(IceVsSeaReference[!is.na(IceVsSeaReference$Int)&IceVsSeaReference$type=="Ctr",]$`N atoms [#]`==0)
hasnoNiceonly<-sum(uniqueallice2$`N atoms [#]`==0)
dfN<-data.frame(N=c(hasNice,hasNseawater),noN=c(hasnoNice,hasnoNseawater),row.names=c("IceMelt","Ctr"))
dfN1<-data.frame(N=c(hasNiceonly,hasNseawater),noN=c(hasnoNiceonly,hasnoNseawater),row.names=c("IceMelt unique","Ctr"))
dfN2<-data.frame(N=c(hasNice,hasNiceonly),noN=c(hasnoNice,hasnoNiceonly),row.names=c("IceMelt","IceMelt unique"))

mosaicplot(dfN, color = TRUE) 
mosaicplot(dfN1, color = TRUE) 
mosaicplot(dfN2, color = TRUE) 

fisher.test(dfN)$p.value
#fisher.test(dfN1)$p.value
#fisher.test(dfN2)$p.value

hasPice<-sum(IceVsSeaReference[!is.na(IceVsSeaReference$Int)&IceVsSeaReference$type=="IceMelt",]$`P atoms [#]`>0)
hasPseawater<-sum(IceVsSeaReference[!is.na(IceVsSeaReference$Int)&IceVsSeaReference$type=="Ctr",]$`P atoms [#]`>0)
hasPiceonly<-sum(uniqueallice2$`P atoms [#]`>0)
hasnoPice<-sum(IceVsSeaReference[!is.na(IceVsSeaReference$Int)&IceVsSeaReference$type=="IceMelt",]$`P atoms [#]`==0)
hasnoPseawater<-sum(IceVsSeaReference[!is.na(IceVsSeaReference$Int)&IceVsSeaReference$type=="Ctr",]$`P atoms [#]`==0)
hasnoPiceonly<-sum(uniqueallice2$`P atoms [#]`==0)
dfP<-data.frame(P=c(hasPice,hasPseawater),noP=c(hasnoPice,hasnoPseawater),row.names=c("IceMelt","Ctr"))
dfP1<-data.frame(P=c(hasPiceonly,hasPseawater),noP=c(hasnoPiceonly,hasnoPseawater),row.names=c("IceMelt unique","Ctr"))
dfP2<-data.frame(P=c(hasPice,hasPiceonly),noP=c(hasnoPice,hasnoPiceonly),row.names=c("IceMelt","IceMelt unique"))

mosaicplot(dfP, color = TRUE) 
mosaicplot(dfP1, color = TRUE) 
mosaicplot(dfP2, color = TRUE) 

fisher.test(dfP)$p.value
#fisher.test(dfP1)$p.value
#fisher.test(dfP2)$p.value

pvaluespairwisecorrected=data.frame(
  case=c("ice,sea,S","iceonly,sea,S","ice,iceonly,S","ice,sea,N","iceonly,sea,N","ice,iceonly,N","ice,sea,P","iceonly,sea,P","ice,iceonly,P"),
  pvalue=p.adjust(c(fisher.test(dfS)$p.value,fisher.test(dfS1)$p.value,fisher.test(dfS2)$p.value,fisher.test(dfN)$p.value,fisher.test(dfN1)$p.value,fisher.test(dfN2)$p.value,fisher.test(dfP)$p.value,fisher.test(dfP1)$p.value,fisher.test(dfP2)$p.value)))
pvaluespairwisecorrected$signi05<-pvaluespairwisecorrected$pvalue<0.05

pvaluespairwisecorrected2=data.frame(
  case=c("ice,sea,S","ice,sea,N","ice,sea,P"),
  pvalue=p.adjust(c(fisher.test(dfS)$p.value,fisher.test(dfN)$p.value,fisher.test(dfP)$p.value)))
pvaluespairwisecorrected2$signi05<-pvaluespairwisecorrected2$pvalue<0.05
#result with p values
pvaluespairwisecorrected2

# relative abundance of S,N,P in unique ice formulae %
sum(uniqueallice2$`S atoms [#]`>0)/nrow(uniqueallice2)*100
sum(uniqueallice2$`P atoms [#]`>0)/nrow(uniqueallice2)*100
sum(uniqueallice2$`N atoms [#]`>0)/nrow(uniqueallice2)*100

########################################

## Comparison with other Arctic FT datasets
# Handle empty strings + whitespaces
# Harmonize sulfur formatting (S / S1)
# Ensure no duplicates per dataset
DOMcomp <- read.csv(
  "DOM-ArcticComp.txt", h=T, sep="\t", check.names=F,
  stringsAsFactors=F, skipNul=T) %>% 
  mutate_all(~ na_if(trimws(.), "")) %>%   # Remove leading/trailing whitespaces
  mutate_all(~ na_if(., "NA")) %>%         # Replace "NA" string with NA values
  mutate(across(c("HG-Sediment","UIW-DOM"), ~gsub("S$", "S1", .))) %>%
  mutate(across(everything(), ~{
    unique_vals <- unique(.x)  # Get unique values in the column
    # If the column has fewer unique values, pad with NA to match original length
    length(unique_vals) <- length(.x)  
    unique_vals})) %>%
  # Ensure no column has duplicate rows with the same formula (ignoring NaN/NA)
  distinct() %>%
  filter(if_any(everything(), ~ !is.na(.)))

# Load databases
db_npAtlas <- read.csv(
  "npAtlas.txt", h=T, sep="\t", check.names=F,
  stringsAsFactors=F, skipNul=T)  %>%
  mutate(formula = gsub("S$", "S1", formula))
db_Chebi <- read.csv(
  "Chebi.txt", h=T, sep="\t", check.names=F,
  stringsAsFactors=F, skipNul=T)  %>%
  mutate(formula = gsub("S$", "S1", formula))

# DOM comparison: count shared formula
DOMcount <- DOMcomp %>%
  pivot_longer(everything(), names_to="dataset", values_to="formula") %>%
  group_by(formula) %>%
  filter(n() > 1) %>%
  summarise(datasets = paste(sort(unique(dataset)), collapse = " / "), .groups="drop") %>%
  drop_na() %>%
  group_by(datasets) %>%
  summarise(shared = n(), .groups="drop")

# Summarize
formulaeTotal <- DOMcomp %>%
  summarise(across(everything(), ~n_distinct(.))) %>%
  pivot_longer(everything(), names_to="dataset", values_to="sum")  

# Set function
getSum <- function(combination, formulaeTotal) {
  sum(formulaeTotal$sum[formulaeTotal$dataset %in% combination])}

# Generate combinations of all datasets 
formulaeComb <- lapply(2:nrow(formulaeTotal), function(i) {
  combn(formulaeTotal$dataset, i, simplify = F)}) %>%
  unlist(., recursive = F)

# Create summary dataframe
formulaeComb <- data.frame(
  datasets = sapply(formulaeComb, function(x) paste(sort(x), collapse = " / ")),
  sum = sapply(formulaeComb, getSum, formulaeTotal = formulaeTotal))

# Calculate shared fractions
DOMcount <- left_join(DOMcount, formulaeComb) %>%
  mutate(percentage = (shared / sum) * 100)

######################################

# Specific formula patterns? 
DOMshared <- DOMcomp %>%
  pivot_longer(everything(), names_to="dataset", values_to="formula") %>%
  group_by(formula) %>%
  #filter(n_distinct(dataset) == 3) %>%  
  filter(n_distinct(dataset) == 3 & all(c("Unique-Ice") %in% dataset)) %>%  
  summarise(datasets = paste(sort(unique(dataset)), collapse=" / "), .groups="drop") %>%
  drop_na()

DOMnonshared <- DOMcomp %>%
  pivot_longer(everything(), names_to="dataset", values_to="formula") %>%
  group_by(formula) %>%
  filter(n_distinct(dataset) == 1) %>%  
  summarise(datasets = paste(sort(unique(dataset)), collapse = " / "), .groups="drop") %>%
  drop_na()

# Hits in compounds databases?
hits.chebi <- inner_join(DOMshared, db_Chebi)
hits.npAtlas <- inner_join(DOMshared, db_npAtlas)
