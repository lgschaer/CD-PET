#-----MAKING A PHYLOSEQ OBJECT-----#

library(ggpubr)
library(tidyverse)
library(phyloseq)
library(csv)

#load data into R for phyloseq analysis



#load sample data
sdata <- as.csv("/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/samples2/DCPET_metadata.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(sdata)

sdata2 <- sdata %>% 
  rownames_to_column(var = "Sample_ID") %>%
  mutate(SampleName = Sample_Name) %>%
  select(c("SampleName", "Sample_Name", "Sample_ID", "Enrichment", "Carbon", "Media", "Replicate", "Media_Carbon"))
head(sdata2)                                                       #view data to make sure everything is OK

sdata3 <- sdata2 %>%                                               #make data to use in phyloseq object with Sample_ID as rownames
  column_to_rownames(var = "Sample_Name")
head(sdata3)                                                       #view data to make sure everything is OK

#load sequence table
sequence_table <- readRDS("/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/samples2/dada2_output/seqtab.rds")
colnames(sequence_table) <- NULL                                   #remove column names

#make nonzero subset to remove all columns with zero taxa counts
sequence_table <- as.matrix(sequence_table)                                #change to matrix format
m <- (colSums(sequence_table, na.rm=TRUE) != 0)                        #T if colSum is not 0, F otherwise
nonzero <- sequence_table[, m]                                         #all the non-zero columns

#load taxa table
taxa_table <- readRDS("/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/samples2/dada2_output/taxa.rds")
taxa_table <- as.matrix(taxa_table)                                #change to matrix format
taxa_table[1:5,1:5]                                                #view to make sure everything looks good

#make phyloseq object
samdata = sample_data(sdata3)                                      #define sample data
colnames(nonzero) <- NULL                                          #remove column names from "nonzero"
seqtab = otu_table(nonzero, taxa_are_rows = FALSE)                 #define sequence table
taxtab = tax_table(taxa_table)                                     #define taxa table
rownames(taxtab) <- NULL                                           #remove rownames from taxa table


sample_names(samdata)
sample_names(seqtab)

phyloseq_object_all = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))
phyloseq_object_all

sample_sums(phyloseq_object_all)

#sample counts before rarefying
sample_counts <- sample_data(phyloseq_object_all) %>%
  group_by(Enrichment, Carbon, Media) %>%
  mutate(Count = 1) %>%
  summarise(SumCount = sum(Count))
sample_counts

#normalize data
#Delete samples with a mean of less than 1000
samplesover1000_all <- subset_samples(phyloseq_object_all, sample_sums(phyloseq_object_all) > 1000)

#Check if there are OTUs with no counts, if so how many?
any(taxa_sums(samplesover1000_all) == 0)
sum(taxa_sums(samplesover1000_all) == 0)

#Prune OTUs with no counts 
prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)
any(taxa_sums(prune_samplesover1000_all) == 0)

#make sure seed is set the same each time, set to 81 here
rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))

#number rarefied at:
min(sample_sums(prune_samplesover1000_all))

#filter out eukaryotes and mitochondria
justbacteria <- rarefy_samplesover1000_all %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  )
justbacteria

#saveRDS file
saveRDS(justbacteria, file = "/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/samples2/phyloseq_output/dcpet_rarefied_nochloroplasts.rds")

#-----ALPHA DIVERISTY-----#

#Violin plot of alpha diversity Observed OTUs and Shannon Diversity (with color)
colors <- c("PositiveControl" = "lightblue", "Emma2" = "purple", "Laura1" = "orange")

justbacteria %>%                                                     #phyloseq object
  plot_richness(
    x = "Enrichment",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_violin(aes(fill = Enrichment), show.legend = FALSE)+          #make violin plot, set fill aes to sampletype
  geom_boxplot(width=0.1) +                                          #add boxplot, set width
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = colors)+                         #set fill colors
  ggtitle("Alpha Diversity") +                                       #add title
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position

colors2 <- c("purple", "orange", "green", "lightblue", "violet")

justbacteria2 <- subset_samples(justbacteria, Enrichment!="PositiveControl")
justbacteria2

justbacteria2 %>%                                                     #phyloseq object
  plot_richness(
    x = "Media_Carbon",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  #geom_violin(aes(fill = Media_Carbon), show.legend = FALSE)+          #make violin plot, set fill aes to sampletype
  geom_boxplot(width=4) +                                          #add boxplot, set width
  geom_jitter(aes(fill = Enrichment),color = "black", size = 5, shape = 21)+
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = colors2)+                         #set fill colors
  ggtitle("Alpha Diversity") +                                       #add title
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position

#-----ALPHA DIVERISTY STATISTICS-----#

#power analysis to determime whether we have sufficent power use an ANOVA with non-normally distributed data.
#count samples
sample_counts <- sample_data(justbacteria2) %>%
  group_by(Media_Carbon) %>%
  mutate(Count = 1)%>%
  summarise(SumCount = sum(Count)) %>%
  mutate(countpercategory = mean(SumCount))
sample_counts
head(sample_data(justbacteria))

library(pwr)

pwr.anova.test(k=4, f=0.23, sig.level=.05, power=.8)
#this data set is too small to meet the assumptions of an ANOVA

#add alpha diversity data to a data frame
richness <- justbacteria2 %>%
  estimate_richness(measures = c("Observed", "Shannon")) %>%           #specify which measures
  rownames_to_column(var = "Sample_Name") %>%                          #add column name to SampleID column
  as_tibble() 
head(richness)

alphadiv <- richness %>%
  left_join(sdata2, by = "Sample_Name")
head(alphadiv)

#Kruskal-Wallis Test
set.seed(81)

##BY ENRICHMENT
#Observed
kruskal.test(Observed ~ Enrichment, data = alphadiv) 

#Shannon
kruskal.test(Shannon ~ Enrichment, data = alphadiv) 

##BY MEDIA/CARBON & ENRICHMENT
head(alphadiv)

L1_alphadiv <- filter(alphadiv, Enrichment == "Laura1")
E2_alphadiv <- filter(alphadiv, Enrichment == "Emma2")

#Observed
kruskal.test(Observed ~ Media_Carbon, data = L1_alphadiv) 
kruskal.test(Observed ~ Media_Carbon, data = E2_alphadiv) 

#Shannon
kruskal.test(Shannon ~ Media_Carbon, data = L1_alphadiv) 
kruskal.test(Shannon ~ Media_Carbon, data = E2_alphadiv) 

#Dunn test (post hoc)

##Observed
#dunnO <- dunnTest(Observed ~ sampletype,
#                  data=alphadiv,
#                  method="bh")
#dunnO

#dunnO <- dunnO$res
#View(dunnO)

##Shannon
dunnS <- dunnTest(Shannon ~ Media_Carbon,
                  data=E2_alphadiv,
                  method="bh")
dunnS

#dunnS <- dunnS$res
#View(dunnS)

#-----BETA DIVERSITY-----#

#ordination
distance <- ordinate(
  physeq = justbacteria2, 
  method = "PCoA", 
  distance = "bray"
)
summary(distance)
distance

#t-SNE plot
head(sdata3)

tsne <- tsne_phyloseq(justbacteria2, distance = "bray", perplexity = 15, dimensions = 2,
                      precomputed_distance = NULL, pseudocounts = 1, verbose = 1,
                      rng_seed = 81, philr_options = list(), control = list())
summary(tsne)

justbacteria2

#tSNE Plot
plot_tsne_phyloseq(justbacteria2, tsne, color = "Enrichment", shape = "Carbon") +
  geom_point(aes(color = Enrichment, fill = Enrichment), color = "black", size = 5, show.legend = TRUE) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  scale_fill_manual(values = colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(color = FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

#PCOA Plot

#ordination
all_pcoa <- ordinate(
  physeq = justbacteria2, 
  method = "PCoA", 
  distance = "bray"
)

colors3 <- c("green", "blue", "lightblue")

#plot
PlotB <- plot_ordination(
  physeq = justbacteria2,                                                          #phyloseq object
  ordination = all_pcoa)+                                                #ordination
  geom_point(aes(fill = Carbon, shape = Enrichment), size = 6) +                         #sets fill color to sampletype
  scale_fill_manual(values = colors3) +
  scale_shape_manual(values = c(21, 22, 23))+
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank())+                                      #removes legend title
    #legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.position = "bottom",
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
PlotB

#-----BETA DIVERSTIY STATISTICS-----#

#PERMANOVA
set.seed(81)
head(sample_data(justbacteria2))

#subset phyloseq object, all samples by datatype
one <- subset_samples(justbacteria2, Media_Carbon %in% c("DCPET_ARW", "DCPET_AMW"))
two <- subset_samples(justbacteria2, Media_Carbon %in% c("DCPET_ARW", "DCPET_BH"))
three <- subset_samples(justbacteria2, Media_Carbon %in% c("DCPET_ARW", "TPA_BH"))
four <- subset_samples(justbacteria2, Media_Carbon %in% c("DCPET_AMW", "DCPET_BH"))
five <- subset_samples(justbacteria2, Media_Carbon %in% c("DCPET_AMW", "TPA_BH"))
six <- subset_samples(justbacteria2, Media_Carbon %in% c("DCPET_BH", "TPA_BH"))

enrich_comp <- subset_samples(justbacteria2, Enrichment %in% c("Laura1", "Emma2"))

E1<-subset_samples(justbacteria2, Enrichment == "Emma2")
L1<-subset_samples(justbacteria2, Enrichment == "Laura1")

one <- subset_samples(L1, Media_Carbon %in% c("DCPET_ARW", "DCPET_AMW"))
two <- subset_samples(L1, Media_Carbon %in% c("DCPET_ARW", "DCPET_BH"))
three <- subset_samples(L1, Media_Carbon %in% c("DCPET_ARW", "TPA_BH"))
four <- subset_samples(L1, Media_Carbon %in% c("DCPET_AMW", "DCPET_BH"))
five <- subset_samples(L1, Media_Carbon %in% c("DCPET_AMW", "TPA_BH"))
six <- subset_samples(L1, Media_Carbon %in% c("DCPET_BH", "TPA_BH"))


# Calculate bray curtis distance matrix, all samples
bray1 <- phyloseq::distance(one, method = "bray")
bray2 <- phyloseq::distance(two, method = "bray")
bray3 <- phyloseq::distance(three, method = "bray")
bray4 <- phyloseq::distance(four, method = "bray")
bray5 <- phyloseq::distance(five, method = "bray")
bray6 <- phyloseq::distance(six, method = "bray")

bray_enrich_comp <- phyloseq::distance(enrich_comp, method = "bray")

# make a data frame from the sample_data, all samples
sam1 <- data.frame(sample_data(one))
sam2 <- data.frame(sample_data(two))
sam3 <- data.frame(sample_data(three))
sam4 <- data.frame(sample_data(four))
sam5 <- data.frame(sample_data(five))
sam6 <- data.frame(sample_data(six))

sam_enrich_comp <- data.frame(sample_data(enrich_comp))

# Adonis test, all samples
adonis(bray1 ~ Media_Carbon, data = sam1)
adonis(bray2 ~ Media_Carbon, data = sam2)
adonis(bray3 ~ Media_Carbon, data = sam3)
adonis(bray4 ~ Media_Carbon, data = sam4)
adonis(bray5 ~ Media_Carbon, data = sam5)
adonis(bray6 ~ Media_Carbon, data = sam6)

adonis(bray_enrich_comp ~ Enrichment, data = sam_enrich_comp)

#EXPLORING TAXA

phyloseq_object_all
justbacteria

#Summarize abundance of each class
genusabundance <- rarefy_samplesover1000_all %>%
  tax_glom(taxrank = "Genus") %>%                      # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus) 
head(genusabundance)

#just the positive control
pos <- genusabundance %>%
  select(Phylum, Class, Family, Genus, Sample, Abundance, Enrichment, Media_Carbon, Replicate) %>%
  filter(Abundance != 0) %>%
  filter(Enrichment == "PositiveControl") %>%
  mutate(
    Phylum = as.character(Phylum),
    Class = as.character(Class),
    Family = as.character(Family),
    Genus = as.character(Genus)) %>%
  group_by(Enrichment, Media_Carbon, Genus)%>%
  summarise(Abundance = sum(Abundance)) %>%
  mutate(Genus.2p = ifelse(Abundance < 0.02, "<2%", Genus))%>%
  arrange(desc(Abundance, Sample))                           #Arrange with descending abundances
head(pos)

View(pos)

#save color palatte
colors10 <- c(
  "black",   "darkcyan",     "orchid1",   "green",       "blue",   
  "grey47",  "cyan",    "coral1",     "yellow",    "darkgreen",   "palegoldenrod",    
  "grey77",  "darkblue",     "orange",    "red",         "mediumpurple1", "tan4",   "purple4",
  "dodgerblue",    "white", "firebrick", "yellowgreen", "magenta", "blue", "green", "red", "orchid", "lightblue"
) 

ggplot(pos)+
  geom_col(mapping = aes(x = Media_Carbon, y = Abundance, fill = Genus), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(rows = vars(Enrichment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors10) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        title = element_text(size = 25))


#Select and summarize necessary variables
all <- genusabundance %>%
  select(Phylum, Class, Family, Genus, Sample, Abundance, Enrichment, Media_Carbon, Replicate) %>%
  filter(Abundance != 0) %>%
  filter(Enrichment != "PositiveControl") %>%
  mutate(
    Phylum = as.character(Phylum),
    Class = as.character(Class),
    Family = as.character(Family),
    Genus = as.character(Genus)
    )
head(all)

phylum <- all %>%
  dplyr::group_by(Enrichment, Media_Carbon, Phylum)%>%
  summarise(Abundance = sum(Abundance)) %>%
  mutate(Phylum.1p = ifelse(Abundance < 0.01, "<1%", Phylum))%>%
  arrange(desc(Abundance, Sample))                           #Arrange with descending abundances
head(phylum)

class <- all %>%
  group_by(Enrichment, Media_Carbon, Class)%>%
  summarise(Abundance = sum(Abundance)) %>%
  mutate(Class = ifelse(Abundance < 0.01, "<1%", Class))%>%
  arrange(desc(Abundance, Sample))                           #Arrange with descending abundances
head(class)

family <- all %>%
  group_by(Enrichment, Media_Carbon, Family)%>%
  summarise(Abundance = sum(Abundance)) %>%
  mutate(Family = ifelse(Abundance < 0.01, "<1%", Family))%>%
  arrange(desc(Abundance, Sample))                           #Arrange with descending abundances
head(family)

genus <- all %>%
  group_by(Enrichment, Media_Carbon, Genus)%>%
  summarise(Abundance = sum(Abundance)) %>%
  mutate(Genus = ifelse(Abundance < 0.02, "<2%", Genus))%>%
  arrange(desc(Abundance, Sample))                           #Arrange with descending abundances
head(genus)

#MAKING A TAXA PLOT 

#save color palatte
colors10 <- c(
  "black",   "darkcyan",     "orchid1",   "green",       "blue",   
  "grey47",  "cyan",    "coral1",     "yellow",    "darkgreen",   "palegoldenrod",    
  "grey77",  "darkblue",     "orange",    "red",         "mediumpurple1", "tan4",   "purple4",
   "dodgerblue",    "white", "firebrick", "yellowgreen", "magenta", "blue", "green", "red", "orchid", "lightblue"
)  

length(colors10)

#plot
phy <- ggplot(phylum)+
  geom_col(mapping = aes(x = Media_Carbon, y = Abundance, fill = Phylum), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(rows = vars(Enrichment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors10) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 10),
    #    legend.title = element_blank(),
        title = element_text(size = 18))

cla <- ggplot(class)+
  geom_col(mapping = aes(x = Media_Carbon, y = Abundance, fill = Class), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(rows = vars(Enrichment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors10) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 10),
     #   legend.title = element_blank(),
        title = element_text(size = 18))

fam <- ggplot(family)+
  geom_col(mapping = aes(x = Media_Carbon, y = Abundance, fill = Family), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(rows = vars(Enrichment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors10) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 10),
       # legend.title = element_blank(),
        title = element_text(size = 18))

gen <- ggplot(genus)+
  geom_col(mapping = aes(x = Media_Carbon, y = Abundance, fill = Genus), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(rows = vars(Enrichment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors10) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 10),
      #  legend.title = element_blank(),
        title = element_text(size = 18))

ggarrange(phy, cla, fam, gen, 
          align = "hv", 
          legend = "right", 
          ncol = 2, nrow = 2)
         # hjust = -6, vjust = 0.7,
         # labels = c("Phylum", "Family", "Class", "Genus"))

#genus by replicate taxa plot
genus2 <- all %>%
  group_by(Enrichment, Media_Carbon, Genus, Replicate)%>%
  summarise(Abundance = sum(Abundance)) %>%
  mutate(Genus.2p = ifelse(Abundance < 0.02, "<2%", Genus))%>%
  arrange(desc(Abundance, Sample))                           #Arrange with descending abundances
head(genus2)

ggplot(genus2)+
  geom_col(mapping = aes(x = Replicate, y = Abundance, fill = Genus.2p), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(rows = vars(Enrichment), cols = vars(Media_Carbon))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors10) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        title = element_text(size = 25))
