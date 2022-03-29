####################################################################################################
# Abundance and species richness for Rhabdastrella_globostellata_microbiome
####################################################################################################
# By Pierre-Louis Stenger
# mail : Pierrelouis.stenger@gmail.com
# For all phylum, use the "all_phylum.txt" file
# For chloroflexi : on the cholroflexi.txt file  which is a subset on chloroflexi from the "all_data_Rhabastrella.xlsx", here we are looking for geus.

library(dplyr)
library(tidyr)
library(ggplot2)
#devtools::install_github("hms-dbmi/UpSetR")
library(UpSetR)
library(VennDiagram)

####################################################################################################
####################################################################################################
####################################################################################################
# ALL_DATA
####################################################################################################
####################################################################################################
####################################################################################################

####################################################################################################
# Import the data
####################################################################################################

setwd("~/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_01_Rhabdastrella_globostellata_and_Mycale_microbiome/Rhabdastrella_globostellata_and_Mycale_microbiome/09_Venn_Diagramm")

taxo <- read.table("taxonomy.tsv", sep="\t", header=T)
head(taxo)
str(taxo)

level <- read.csv("level-7.csv",  header=T)
head(level)

#level_02 <- subset(level, select = - c(Dataset, FullName, Acronyme, time, Site))   
level_02 <- subset(level, select = - c(Dataset, Acronyme, time, Site))   


#head(t(level_02))

str(level_02)


pb_long <- level_02 %>%
  gather(key = FullName, value = value, -index)
head(pb_long)
tail(pb_long)


pb_long_02 <- pb_long[!grepl("FullName", pb_long$FullName),]
tail(pb_long_02)
head(pb_long_02)
str(pb_long_02)

pb_long_02$value <- as.numeric(pb_long_02$value)


#################################################
# Kingdom 
word.list <- strsplit(pb_long_02$FullName, "p__")
n.obs <- sapply(word.list, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(word.list, "[", i = seq.max))
head(mat)

kingdom <- mat[,1]
phylum <- mat[,2]

#################################################
# phylum
word.list <- strsplit(phylum, "c__")
n.obs <- sapply(word.list, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(word.list, "[", i = seq.max))
head(mat)

phylum <- mat[,1]
class <- mat[,2]


#################################################
# class
word.list <- strsplit(class, "o__")
n.obs <- sapply(word.list, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(word.list, "[", i = seq.max))
head(mat)

class <- mat[,1]
order <- mat[,2]

#################################################
# order
word.list <- strsplit(order, "f__")
n.obs <- sapply(word.list, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(word.list, "[", i = seq.max))
head(mat)

order <- mat[,1]
family <- mat[,2]

#################################################
# family
word.list <- strsplit(family, "g__")
n.obs <- sapply(word.list, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(word.list, "[", i = seq.max))
head(mat)

family <- mat[,1]
genus <- mat[,2]

#################################################
# genus
word.list <- strsplit(genus, "s__")
n.obs <- sapply(word.list, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(word.list, "[", i = seq.max))
head(mat)

genus <- mat[,1]
species <- mat[,2]


taxo_split <- data.frame(kingdom, phylum, class, order, family, genus, species)
taxo_split[] <- lapply(taxo_split, gsub, pattern='d__', replacement='')
taxo_split[] <- lapply(taxo_split, gsub, pattern='__', replacement='')
taxo_split[] <- lapply(taxo_split, gsub, pattern='[.]', replacement='')

head(taxo_split)
asv_link <- paste(taxo_split$kingdom, "-", taxo_split$phylum, "-", taxo_split$class,"-",  taxo_split$order, "-", taxo_split$family, "-", taxo_split$genus, "-", taxo_split$pecies)

length(taxo_split$kingdom)
length(pb_long_02$FullName)

dat <- cbind(pb_long_02, taxo_split, asv_link)

head(dat)
dat$index<-gsub("_V1V3","",as.character(dat$index))
dat$index<-gsub("Bourake_","",as.character(dat$index))
dat$index<-gsub("_2019","",as.character(dat$index))
head(dat)

########################################################
metadata <- read.table("taxonomy.tsv", header=T, sep="\t")
head(metadata)


#################################################
# Kingdom 
word.list <- strsplit(metadata$Taxon, "; p__")
n.obs <- sapply(word.list, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(word.list, "[", i = seq.max))
head(mat)

kingdom <- mat[,1]
phylum <- mat[,2]

#################################################
# phylum
word.list <- strsplit(phylum, "; c__")
n.obs <- sapply(word.list, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(word.list, "[", i = seq.max))
head(mat)

phylum <- mat[,1]
class <- mat[,2]


#################################################
# class
word.list <- strsplit(class, "; o__")
n.obs <- sapply(word.list, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(word.list, "[", i = seq.max))
head(mat)

class <- mat[,1]
order <- mat[,2]

#################################################
# order
word.list <- strsplit(order, "; f__")
n.obs <- sapply(word.list, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(word.list, "[", i = seq.max))
head(mat)

order <- mat[,1]
family <- mat[,2]

#################################################
# family
word.list <- strsplit(family, "; g__")
n.obs <- sapply(word.list, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(word.list, "[", i = seq.max))
head(mat)

family <- mat[,1]
genus <- mat[,2]

#################################################
# genus
word.list <- strsplit(genus, "; s__")
n.obs <- sapply(word.list, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(word.list, "[", i = seq.max))
head(mat)

genus <- mat[,1]
species <- mat[,2]


taxo_split <- data.frame(kingdom, phylum, class, order, family, genus, species)
taxo_split[] <- lapply(taxo_split, gsub, pattern='d__', replacement='')
taxo_split[] <- lapply(taxo_split, gsub, pattern='__', replacement='')
taxo_split[] <- lapply(taxo_split, gsub, pattern='[.]', replacement='')

head(taxo_split)

asv_link <- paste(taxo_split$kingdom, "-", taxo_split$phylum, "-", taxo_split$class,"-",  taxo_split$order, "-", taxo_split$family, "-", taxo_split$genus, "-", taxo_split$pecies)
head(asv_link)
head(metadata)

asv_link_02 <- data.frame(cbind(asv_link, metadata$Feature.ID))
head(asv_link_02)

head(dat)

asv_link_03 <- merge(asv_link_02, dat, by.x='asv_link', by.y='asv_link')
head(asv_link_03)

dat_ASV <- data.frame(asv_link_03$V2, asv_link_03$index, asv_link_03$asv_link, asv_link_03$kingdom, asv_link_03$phylum, asv_link_03$class, asv_link_03$order, asv_link_03$family, asv_link_03$genus, asv_link_03$species, asv_link_03$value)
colnames(dat_ASV) <- c("ASV", "sites", "taxonomy", "kingdom", "phylum", "class", "order", "family", "genus", "species", "value")
head(dat_ASV)
str(dat_ASV)

unique(dat_ASV$sites)

dat_ASV_02 <- data.frame(dat_ASV$ASV, dat_ASV$sites, dat_ASV$value)
colnames(dat_ASV_02) <- c("ASV", "sites", "value")
head(dat_ASV_02)
str(dat_ASV_02)

pb_long_dat_ASV <- spread(dat_ASV,
                        key=sites,
                        value=value)
head(pb_long_dat_ASV)

length(pb_long_dat_ASV$ASV)

dat_ASV_03 <- data.frame(pb_long_dat_ASV$ASV, 
                         pb_long_dat_ASV$Mycale_1, 
                         pb_long_dat_ASV$Mycale_2, 
                         pb_long_dat_ASV$Mycale_3, 
                         pb_long_dat_ASV$Rhabdastrella_globostellata_1, 
                         pb_long_dat_ASV$Rhabdastrella_globostellata_2, 
                         pb_long_dat_ASV$Rhabdastrella_globostellata_3, 
                         pb_long_dat_ASV$SeaWater_downstream_0_2micro, 
                         pb_long_dat_ASV$SeaWater_upstream_0_2micro)


colnames(dat_ASV_03) <- c("ASV", 
                          "Mycale_1", 
                          "Mycale_2", 
                          "Mycale_3", 
                          "Rhabdastrella_globostellata_1", 
                          "Rhabdastrella_globostellata_2", 
                          "Rhabdastrella_globostellata_3", 
                          "SeaWater_downstream_0_2micro", 
                          "SeaWater_upstream_0_2micro"
)  

head(dat_ASV_03)


dat_ASV_04 <- dat_ASV_03 %>% mutate(Mycale_1  = replace(Mycale_1 , Mycale_1 != 0, 1))
dat_ASV_04 <- dat_ASV_04 %>% mutate(Mycale_2  = replace(Mycale_2 , Mycale_2 != 0, 1))
dat_ASV_04 <- dat_ASV_04 %>% mutate(Mycale_3  = replace(Mycale_3 , Mycale_3 != 0, 1))
dat_ASV_04 <- dat_ASV_04 %>% mutate(Rhabdastrella_globostellata_1  = replace(Rhabdastrella_globostellata_1 , Rhabdastrella_globostellata_1 != 0, 1))
dat_ASV_04 <- dat_ASV_04 %>% mutate(Rhabdastrella_globostellata_2  = replace(Rhabdastrella_globostellata_2 , Rhabdastrella_globostellata_2 != 0, 1))
dat_ASV_04 <- dat_ASV_04 %>% mutate(Rhabdastrella_globostellata_3  = replace(Rhabdastrella_globostellata_3 , Rhabdastrella_globostellata_3 != 0, 1))
dat_ASV_04 <- dat_ASV_04 %>% mutate(SeaWater_downstream_0_2micro  = replace(SeaWater_downstream_0_2micro , SeaWater_downstream_0_2micro != 0, 1))
dat_ASV_04 <- dat_ASV_04 %>% mutate(SeaWater_upstream_0_2micro  = replace(SeaWater_upstream_0_2micro , SeaWater_upstream_0_2micro != 0, 1))

head(dat_ASV_04)  
tail(dat_ASV_04)  

####################################################################################################
# UpSet plot
####################################################################################################

## WARNING, for all data, plot made in 2 hours !

#upset(dat_ASV_04, sets = c("Bourake_Mycale_2019_1_V1V3",
#                           "Bourake_Mycale_2019_2_V1V3",
#                           "Bourake_Mycale_2019_3_V1V3",
#                           "Bourake_Rhabdastrella_globostellata_2019_1_Chimeric_V1V3V4",
#                           "Bourake_Rhabdastrella_globostellata_2019_1_V1V3",
#                           "Bourake_Rhabdastrella_globostellata_2019_1_V4",
#                           "Bourake_Rhabdastrella_globostellata_2019_2_Chimeric_V1V3V4",
#                           "Bourake_Rhabdastrella_globostellata_2019_2_V1V3",
#                           "Bourake_Rhabdastrella_globostellata_2019_2_V4",
#                           "Bourake_Rhabdastrella_globostellata_2019_3_Chimeric_V1V3V4",
#                           "Bourake_Rhabdastrella_globostellata_2019_3_V1V3",
#                           "Bourake_Rhabdastrella_globostellata_2019_3_V4",
#                           "Bourake_SeaWater_downstream_0_2micro_V1V3",
#                           "Bourake_SeaWater_upstream_0_2micro_V1V3",
#                           "Cafard_island_autochtone_Rhabdastrella_globostellata_2019_V1V3",
#                           "Cafard_island_transplant_Rhabdastrella_globostellata_2019_12weeks_V1V3",
#                           "Yate_Rhabdastrella_globostellata_2019_1_V1V3",
#                           "Yate_Rhabdastrella_globostellata_2019_2_V1V3",
#                           "Yate_Rhabdastrella_globostellata_2019_3_V1V3"), sets.bar.color = "#56B4E9",
#order.by = "freq", empty.intersections = "on")

upset(dat_ASV_04, sets = c("Mycale_1",
                           "Mycale_2",
                           "Mycale_3",
                           "Rhabdastrella_globostellata_1",
                           "Rhabdastrella_globostellata_2",
                           "Rhabdastrella_globostellata_3"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")


upset(dat_ASV_04, sets = c("Mycale_1",
                           "Mycale_2",
                           "Mycale_3"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")




library(VennDiagram)

sum(dat_ASV_04$Mycale_1)
sum(dat_ASV_04$Mycale_2)
sum(dat_ASV_04$Mycale_3)

draw.triple.venn(area1 = sum(dat_ASV_04$Mycale_1), 
                 area2 = sum(dat_ASV_04$Mycale_2), 
                 area3 = sum(dat_ASV_04$Mycale_3), 
                 
                 n12 = 28+550, n23 = 8+550, n13 = 19+550, 
                 n123 = 550, category = c("Mycale 01", "Mycale 02", "Mycale 03"), lty = "blank", 
                 fill = c("skyblue", "cadetblue4", "dodgerblue3"))


####################################################################################################
# ONLY MYCALE & RHABA V1V3 & Sea Water Downstream and Upstream
####################################################################################################

####################################################################################################
# by taxonomic drawer
####################################################################################################

head(dat_ASV)

link <- data.frame(dat_ASV$ASV, dat_ASV$kingdom, dat_ASV$phylum, dat_ASV$class, dat_ASV$order, dat_ASV$family, dat_ASV$genus, dat_ASV$species)
colnames(link) <- c("ASV", "kingdom", "phylum", "class", "order", "family", "genus", "species")
head(link)

link_02 <- link[!duplicated(link), ]
head(link_02)
tail(link_02)
length(link_02$ASV)

head(dat_ASV_04)
tail(dat_ASV_04)
length(unique(dat_ASV_04$ASV))

dat_ASV_marker <- dat_ASV_04

dat_ASV_marker_03 <- merge(dat_ASV_marker, link_02, by="ASV",  all.x=T)
tail(dat_ASV_marker_03)
length(dat_ASV_marker_03)
length(unique(dat_ASV_marker_03$ASV))
length((dat_ASV_marker_03$ASV))

dat_ASV_marker_03_test <- dat_ASV_marker_03[!duplicated(dat_ASV_marker_03$ASV), ]
length(dat_ASV_marker_03_test$ASV)

dat_ASV_marker_03 <- dat_ASV_marker_03_test

head(dat_ASV_marker_03)
tail(dat_ASV_marker_03)
unique(dat_ASV_marker_03$phylum)

length(unique(dat_ASV_marker_03$ASV))

upset(dat_ASV_marker_03, sets = c("Mycale sp. individual 1", "Mycale sp. individual 2", "Mycale sp. individual 3","Rhabdastrella globostellata 1", "Rhabdastrella globostellata 2", "Rhabdastrella globostellata 3","SeaWater downstream 0.2 micro",
                                  "SeaWater upstream 0.2 micro"
), sets.bar.color = c( "coral3", "chocolate3", "chocolate1", "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod", "cadetblue3", "cadetblue3"),
order.by = "freq", empty.intersections = "on")


####################################################################################################
# ONLY MYCALE & RHABA V1V3
####################################################################################################

####################################################################################################
# by taxonomic drawer
####################################################################################################

head(dat_ASV)

link <- data.frame(dat_ASV$ASV, dat_ASV$kingdom, dat_ASV$phylum, dat_ASV$class, dat_ASV$order, dat_ASV$family, dat_ASV$genus, dat_ASV$species)
colnames(link) <- c("ASV", "kingdom", "phylum", "class", "order", "family", "genus", "species")
head(link)

link_02 <- link[!duplicated(link), ]
head(link_02)
tail(link_02)
length(link_02$ASV)

head(dat_ASV_marker)
tail(dat_ASV_marker)
length(unique(dat_ASV_marker$ASV))

dat_ASV_marker_03 <- merge(dat_ASV_marker, link_02, by="ASV",  all.x=T)
tail(dat_ASV_marker_03)
length(dat_ASV_marker_03)
length(unique(dat_ASV_marker_03$ASV))
length((dat_ASV_marker_03$ASV))

dat_ASV_marker_03_test <- dat_ASV_marker_03[!duplicated(dat_ASV_marker_03$ASV), ]
length(dat_ASV_marker_03_test$ASV)

dat_ASV_marker_03 <- dat_ASV_marker_03_test

head(dat_ASV_marker_03)
tail(dat_ASV_marker_03)
unique(dat_ASV_marker_03$phylum)

length(unique(dat_ASV_marker_03$ASV))

upset(dat_ASV_marker_03, sets = c("Mycale sp. individual 1", "Mycale sp. individual 2", "Mycale sp. individual 3","Rhabdastrella globostellata 1", "Rhabdastrella globostellata 2", "Rhabdastrella globostellata 3"
), sets.bar.color = c( "coral3", "chocolate3", "chocolate1", "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
order.by = "freq", empty.intersections = "on")



### CHLOROFLEXI 

dat_ASV_marker_03_chloroflexi <- dat_ASV_marker_03[ which(dat_ASV_marker_03$phylum=='Chloroflexi'), ]
head(dat_ASV_marker_03_chloroflexi)
tail(dat_ASV_marker_03_chloroflexi)
dat_ASV_marker_03_chloroflexi_02 <- dat_ASV_marker_03_chloroflexi[!duplicated(dat_ASV_marker_03_chloroflexi), ]
head(dat_ASV_marker_03_chloroflexi_02)
tail(dat_ASV_marker_03_chloroflexi_02)

upset(dat_ASV_marker_03_chloroflexi_02, sets = c("Mycale sp. individual 1", "Mycale sp. individual 2", "Mycale sp. individual 3","Rhabdastrella globostellata 1", "Rhabdastrella globostellata 2", "Rhabdastrella globostellata 3"
), sets.bar.color = c( "coral3", "chocolate3", "chocolate1", "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
order.by = "freq", empty.intersections = "on")


### SAR202

head(dat_ASV_marker_03)

dat_ASV_marker_03_SAR202_clade <- dat_ASV_marker_03[ which(dat_ASV_marker_03$order=='SAR202_clade'), ]
head(dat_ASV_marker_03_SAR202_clade)
tail(dat_ASV_marker_03_SAR202_clade)
dat_ASV_marker_03_SAR202_clade_02 <- dat_ASV_marker_03_SAR202_clade[!duplicated(dat_ASV_marker_03_SAR202_clade), ]
head(dat_ASV_marker_03_SAR202_clade_02)
tail(dat_ASV_marker_03_SAR202_clade_02)

upset(dat_ASV_marker_03_SAR202_clade_02,  sets = c("Mycale sp. individual 1", "Mycale sp. individual 2", "Mycale sp. individual 3","Rhabdastrella globostellata 1", "Rhabdastrella globostellata 2", "Rhabdastrella globostellata 3"
), sets.bar.color = c( "coral3", "chocolate3", "chocolate1", "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
order.by = "freq", empty.intersections = "on")





####################################################################################################
####################################################################################################
# ONLY MYCALE
####################################################################################################
####################################################################################################

####################################################################################################
# by taxonomic drawer
####################################################################################################

head(dat_ASV)

link <- data.frame(dat_ASV$ASV, dat_ASV$kingdom, dat_ASV$phylum, dat_ASV$class, dat_ASV$order, dat_ASV$family, dat_ASV$genus, dat_ASV$species)
colnames(link) <- c("ASV", "kingdom", "phylum", "class", "order", "family", "genus", "species")
head(link)

link_02 <- link[!duplicated(link), ]
head(link_02)
tail(link_02)
length(link_02$ASV)

head(dat_ASV_marker)
tail(dat_ASV_marker)
length(unique(dat_ASV_marker$ASV))

dat_ASV_marker_03 <- merge(dat_ASV_marker, link_02, by="ASV",  all.x=T)
tail(dat_ASV_marker_03)
length(dat_ASV_marker_03)
length(unique(dat_ASV_marker_03$ASV))
length((dat_ASV_marker_03$ASV))

dat_ASV_marker_03_test <- dat_ASV_marker_03[!duplicated(dat_ASV_marker_03$ASV), ]
length(dat_ASV_marker_03_test$ASV)

dat_ASV_marker_03 <- dat_ASV_marker_03_test

names(dat_ASV_marker_03) <- gsub(x = names(dat_ASV_marker_03), pattern = "Mycale sp. individual ", replacement = "Mycale_0")  

head(dat_ASV_marker_03)
tail(dat_ASV_marker_03)
unique(dat_ASV_marker_03$phylum)

length(unique(dat_ASV_marker_03$ASV))

upset(dat_ASV_marker_03, sets = c("Mycale_1", "Mycale_2", "Mycale_3"
), sets.bar.color = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03$Mycale_1)
area2 <- sum(dat_ASV_marker_03$Mycale_2)
area3 <- sum(dat_ASV_marker_03$Mycale_3)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 13+147, n23 = 2+147, n13 = 0+147, 
                 n123 = 147, category = c("Mycale 01", "Mycale 02", "Mycale 03"), lty = "blank", 
                 fill = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"))

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 29+128, n23 = 0+128, n13 = 10+128, 
                 n123 = 128, category = c("", "", ""), lty = "blank", 
                 fill = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"))


### CHLOROFLEXI 

dat_ASV_marker_03_chloroflexi <- dat_ASV_marker_03[ which(dat_ASV_marker_03$phylum=='Chloroflexi'), ]
head(dat_ASV_marker_03_chloroflexi)
tail(dat_ASV_marker_03_chloroflexi)
dat_ASV_marker_03_chloroflexi_02 <- dat_ASV_marker_03_chloroflexi[!duplicated(dat_ASV_marker_03_chloroflexi), ]
head(dat_ASV_marker_03_chloroflexi_02)
tail(dat_ASV_marker_03_chloroflexi_02)

upset(dat_ASV_marker_03_chloroflexi_02, sets = c("Mycale_1", "Mycale_2", "Mycale_3"
), sets.bar.color = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_chloroflexi_02$Mycale_1)
area2 <- sum(dat_ASV_marker_03_chloroflexi_02$Mycale_2)
area3 <- sum(dat_ASV_marker_03_chloroflexi_02$Mycale_3)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 1+73, n23 = 0+73, n13 = 0+73, 
                 n123 = 73, category = c("Mycale 01", "Mycale 02", "Mycale 03"), lty = "blank", 
                 fill = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"))

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 7+73, n23 = 0+73, n13 = 0+73,
                 n123 = 73, category = c("", "", ""), lty = "blank", 
                 fill = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"), cex=0)

length(unique(dat_ASV_marker_03_chloroflexi_02$ASV))

### SAR202

head(dat_ASV_marker_03)

dat_ASV_marker_03_SAR202_clade <- dat_ASV_marker_03[ which(dat_ASV_marker_03$order=='SAR202_clade'), ]
head(dat_ASV_marker_03_SAR202_clade)
tail(dat_ASV_marker_03_SAR202_clade)
dat_ASV_marker_03_SAR202_clade_02 <- dat_ASV_marker_03_SAR202_clade[!duplicated(dat_ASV_marker_03_SAR202_clade), ]
head(dat_ASV_marker_03_SAR202_clade_02)
tail(dat_ASV_marker_03_SAR202_clade_02)

upset(dat_ASV_marker_03_SAR202_clade_02, sets = c("Mycale_1", "Mycale_2", "Mycale_3"
), sets.bar.color = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_SAR202_clade_02$Mycale_1)
area2 <- sum(dat_ASV_marker_03_SAR202_clade_02$Mycale_2)
area3 <- sum(dat_ASV_marker_03_SAR202_clade_02$Mycale_3)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 0+58, n23 = 0+58, n13 = 0+58, 
                 n123 = 58, category = c("Mycale 01", "Mycale 02", "Mycale 03"), lty = "blank", 
                 fill = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"))

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 0+59, n23 = 0+59, n13 = 0+59, 
                 n123 = 59, category = c("", "", ""), lty = "blank", 
                 fill = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"))

length(unique(dat_ASV_marker_03_SAR202_clade_02$ASV))

################################################
## By phylum

head(dat_ASV_marker_03, n=100)
head(dat_ASV_marker_03)
     
dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$phylum), ]
head(dat_ASV_marker_03_phylum, n=100)
head(dat_ASV_marker_03_phylum)
length(unique(dat_ASV_marker_03_phylum$ASV))

upset(dat_ASV_marker_03_phylum, sets = c("Mycale_1", "Mycale_2", "Mycale_3"
), sets.bar.color = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$Mycale_01)
area2 <- sum(dat_ASV_marker_03_phylum$Mycale_02)
area3 <- sum(dat_ASV_marker_03_phylum$Mycale_03)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 19+358, n23 = 14+358, n13 = 10+358, 
                 n123 = 358, category = c("Mycale 01", "Mycale 02", "Mycale 03"), lty = "blank", 
                 fill = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"))

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 19+358, n23 = 14+358, n13 = 10+358, 
                 n123 = 358, category = c("", "", ""), lty = "blank", 
                 fill = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"))


################################################
## By class


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$class), ]


upset(dat_ASV_marker_03_phylum, sets = c("Mycale_01", "Mycale_02", "Mycale_03"
), sets.bar.color = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$Mycale_01)
area2 <- sum(dat_ASV_marker_03_phylum$Mycale_02)
area3 <- sum(dat_ASV_marker_03_phylum$Mycale_03)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 19+337, n23 = 14+337, n13 = 6+337, 
                 n123 = 337, category = c("Mycale 01", "Mycale 02", "Mycale 03"), lty = "blank", 
                 fill = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"))


################################################
## By order


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$order), ]


upset(dat_ASV_marker_03_phylum, sets = c("Mycale_01", "Mycale_02", "Mycale_03"
), sets.bar.color = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$Mycale_01)
area2 <- sum(dat_ASV_marker_03_phylum$Mycale_02)
area3 <- sum(dat_ASV_marker_03_phylum$Mycale_03)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 19+188, n23 = 14+188, n13 = 6+188, 
                 n123 = 188, category = c("Mycale_01", "Mycale_02", "Mycale_03"), lty = "blank", 
                 fill = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"))

################################################
## By family


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$family), ]


upset(dat_ASV_marker_03_phylum, sets = c("Mycale_01", "Mycale_02", "Mycale_03"
), sets.bar.color = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$Mycale_01)
area2 <- sum(dat_ASV_marker_03_phylum$Mycale_02)
area3 <- sum(dat_ASV_marker_03_phylum$Mycale_03)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 9+173, n23 = 14+173, n13 = 6+173, 
                 n123 = 173, category = c("Mycale_01", "Mycale_02", "Mycale_03"), lty = "blank", 
                 fill = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"))

################################################
## By genus


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$genus), ]


upset(dat_ASV_marker_03_phylum, sets = c("Mycale_01", "Mycale_02", "Mycale_03"
), sets.bar.color = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$Mycale_01)
area2 <- sum(dat_ASV_marker_03_phylum$Mycale_02)
area3 <- sum(dat_ASV_marker_03_phylum$Mycale_03)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 9+171, n23 = 14+171, n13 = 6+171, 
                 n123 = 171, category = c("Mycale_01", "Mycale_02", "Mycale_03"), lty = "blank", 
                 fill = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"))

################################################
## By species


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$species), ]


upset(dat_ASV_marker_03_phylum, sets = c("Mycale_01", "Mycale_02", "Mycale_03"
), sets.bar.color = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$Mycale_01)
area2 <- sum(dat_ASV_marker_03_phylum$Mycale_02)
area3 <- sum(dat_ASV_marker_03_phylum$Mycale_03)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 9+171, n23 = 14+171, n13 = 6+171, 
                 n123 = 171, category = c("Mycale_01", "Mycale_02", "Mycale_03"), lty = "blank", 
                 fill = c("darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"))


















####################################################################################################
# ONLY RHABA
####################################################################################################

####################################################################################################
# by taxonomic drawer
####################################################################################################

head(dat_ASV)
#head(dat_ASV_marker_02_dupli)

link <- data.frame(dat_ASV$ASV, dat_ASV$kingdom, dat_ASV$phylum, dat_ASV$class, dat_ASV$order, dat_ASV$family, dat_ASV$genus, dat_ASV$species)
colnames(link) <- c("ASV", "kingdom", "phylum", "class", "order", "family", "genus", "species")
head(link)

link_02 <- link[!duplicated(link), ]
head(link_02)
tail(link_02)
length(link_02$ASV)

head(dat_ASV_marker)
tail(dat_ASV_marker)
length(unique(dat_ASV_marker$ASV))

dat_ASV_marker_03 <- merge(dat_ASV_marker, link_02, by="ASV",  all.x=T)
tail(dat_ASV_marker_03)
length(dat_ASV_marker_03)
length(unique(dat_ASV_marker_03$ASV))
length((dat_ASV_marker_03$ASV))

dat_ASV_marker_03_test <- dat_ASV_marker_03[!duplicated(dat_ASV_marker_03$ASV), ]
length(dat_ASV_marker_03_test$ASV)

dat_ASV_marker_03 <- dat_ASV_marker_03_test

head(dat_ASV_marker_03)
tail(dat_ASV_marker_03)
unique(dat_ASV_marker_03$phylum)

length(unique(dat_ASV_marker_03$ASV))

upset(dat_ASV_marker_03, sets = c("Rhabdastrella_globostellata_1", "Rhabdastrella_globostellata_2", "Rhabdastrella_globostellata_3"
), sets.bar.color = c("coral3", "chocolate3", "chocolate1"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03$Rhabdastrella_globostellata_1)
area2 <- sum(dat_ASV_marker_03$Rhabdastrella_globostellata_2)
area3 <- sum(dat_ASV_marker_03$Rhabdastrella_globostellata_3)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 18+162, n23 = 8+162, n13 = 13+162, 
                 n123 = 162, category = c("R_globostellata 01", "R_globostellata 02", "R_globostellata 03"), lty = "blank", 
                 fill = c("coral3", "chocolate3", "chocolate1"))

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 2+154, n23 = 14+154, n13 = 6+154, 
                 n123 = 154, category = c("", "", ""), lty = "blank", 
                 fill = c("coral3", "chocolate3", "chocolate1"))


### CHLOROFLEXI 

dat_ASV_marker_03_chloroflexi <- dat_ASV_marker_03[ which(dat_ASV_marker_03$phylum=='Chloroflexi'), ]
head(dat_ASV_marker_03_chloroflexi)
tail(dat_ASV_marker_03_chloroflexi)
dat_ASV_marker_03_chloroflexi_02 <- dat_ASV_marker_03_chloroflexi[!duplicated(dat_ASV_marker_03_chloroflexi), ]
head(dat_ASV_marker_03_chloroflexi_02)
tail(dat_ASV_marker_03_chloroflexi_02)

upset(dat_ASV_marker_03_chloroflexi_02, sets = c("Rhabdastrella_globostellata_1", "Rhabdastrella_globostellata_2", "Rhabdastrella_globostellata_3"
), sets.bar.color = c("coral3", "chocolate3", "chocolate1"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_chloroflexi_02$Rhabdastrella_globostellata_1)
area2 <- sum(dat_ASV_marker_03_chloroflexi_02$Rhabdastrella_globostellata_2)
area3 <- sum(dat_ASV_marker_03_chloroflexi_02$Rhabdastrella_globostellata_3)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 0+84, n23 = 0+84, n13 = 1+84, 
                 n123 = 84, category = c("R_globostellata 01", "R_globostellata 02", "R_globostellata 03"), lty = "blank", 
                 fill = c("coral3", "chocolate3", "chocolate1"))

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 0+81, n23 = 0+81, n13 = 1+81, 
                 n123 = 81, category = c("", "", ""), lty = "blank", 
                 fill = c("coral3", "chocolate3", "chocolate1"))

length(unique(dat_ASV_marker_03_chloroflexi_02$ASV))

### SAR202

head(dat_ASV_marker_03)

dat_ASV_marker_03_SAR202_clade <- dat_ASV_marker_03[ which(dat_ASV_marker_03$order=='SAR202_clade'), ]
head(dat_ASV_marker_03_SAR202_clade)
tail(dat_ASV_marker_03_SAR202_clade)
dat_ASV_marker_03_SAR202_clade_02 <- dat_ASV_marker_03_SAR202_clade[!duplicated(dat_ASV_marker_03_SAR202_clade), ]
head(dat_ASV_marker_03_SAR202_clade_02)
tail(dat_ASV_marker_03_SAR202_clade_02)

upset(dat_ASV_marker_03_SAR202_clade_02, sets = c("Rhabdastrella_globostellata_1", "Rhabdastrella_globostellata_2", "Rhabdastrella_globostellata_3"
), sets.bar.color = c("coral3", "chocolate3", "chocolate1"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_SAR202_clade_02$Rhabdastrella_globostellata_1)
area2 <- sum(dat_ASV_marker_03_SAR202_clade_02$Rhabdastrella_globostellata_2)
area3 <- sum(dat_ASV_marker_03_SAR202_clade_02$Rhabdastrella_globostellata_3)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 0+58, n23 = 0+58, n13 = 0+58, 
                 n123 = 58, category = c("R_globostellata 01", "R_globostellata 02", "R_globostellata 03"), lty = "blank", 
                 fill = c("coral3", "chocolate3", "chocolate1"))

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 0+59, n23 = 0+59, n13 = 0+59, 
                 n123 = 59, category = c("", "", ""), lty = "blank", 
                 fill = c("coral3", "chocolate3", "chocolate1"))


draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 0+111, n23 = 0+111, n13 = 0+111, 
                 n123 = 111, category = c("", "", ""), lty = "blank", 
                 fill = c("chartreuse", "forestgreen", "forestgreen"))

length(unique(dat_ASV_marker_03_SAR202_clade_02$ASV))

################################################
## By phylum

head(dat_ASV_marker_03, n=100)
head(dat_ASV_marker_03)

dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$phylum), ]
head(dat_ASV_marker_03_phylum, n=100)
head(dat_ASV_marker_03_phylum)
length(unique(dat_ASV_marker_03_phylum$ASV))

upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_01", "R_globostellata_02", "R_globostellata_03"
), sets.bar.color = c("coral3", "chocolate3", "chocolate1"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$R_globostellata_01)
area2 <- sum(dat_ASV_marker_03_phylum$R_globostellata_02)
area3 <- sum(dat_ASV_marker_03_phylum$R_globostellata_03)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 0+452, n23 = 8+452, n13 = 10+452, 
                 n123 = 452, category = c("R_globostellata 01", "R_globostellata 02", "R_globostellata 03"), lty = "blank", 
                 fill = c("coral3", "chocolate3", "chocolate1"))


################################################
## By class


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$class), ]


upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_01", "R_globostellata_02", "R_globostellata_03"
), sets.bar.color = c("coral3", "chocolate3", "chocolate1"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$R_globostellata_01)
area2 <- sum(dat_ASV_marker_03_phylum$R_globostellata_02)
area3 <- sum(dat_ASV_marker_03_phylum$R_globostellata_03)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 0+414, n23 = 8+414, n13 = 6+414, 
                 n123 = 414, category = c("R_globostellata 01", "R_globostellata 02", "R_globostellata 03"), lty = "blank", 
                 fill = c("coral3", "chocolate3", "chocolate1"))


################################################
## By order


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$order), ]


upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_01", "R_globostellata_02", "R_globostellata_03"
), sets.bar.color = c("coral3", "chocolate3", "chocolate1"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$R_globostellata_01)
area2 <- sum(dat_ASV_marker_03_phylum$R_globostellata_02)
area3 <- sum(dat_ASV_marker_03_phylum$R_globostellata_03)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 0+249, n23 = 8+249, n13 = 6+249, 
                 n123 = 249, category = c("R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), lty = "blank", 
                 fill = c("coral3", "chocolate3", "chocolate1"))

################################################
## By family


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$family), ]


upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_01", "R_globostellata_02", "R_globostellata_03"
), sets.bar.color = c("coral3", "chocolate3", "chocolate1"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$R_globostellata_01)
area2 <- sum(dat_ASV_marker_03_phylum$R_globostellata_02)
area3 <- sum(dat_ASV_marker_03_phylum$R_globostellata_03)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 0+224, n23 = 5+224, n13 = 6+224, 
                 n123 = 224, category = c("R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), lty = "blank", 
                 fill = c("coral3", "chocolate3", "chocolate1"))

################################################
## By genus


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$genus), ]


upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_01", "R_globostellata_02", "R_globostellata_03"
), sets.bar.color = c("coral3", "chocolate3", "chocolate1"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$R_globostellata_01)
area2 <- sum(dat_ASV_marker_03_phylum$R_globostellata_02)
area3 <- sum(dat_ASV_marker_03_phylum$R_globostellata_03)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 0+224, n23 = 3+224, n13 = 6+224, 
                 n123 = 224, category = c("R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), lty = "blank", 
                 fill = c("coral3", "chocolate3", "chocolate1"))

################################################
## By species


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$species), ]


upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_01", "R_globostellata_02", "R_globostellata_03"
), sets.bar.color = c("coral3", "chocolate3", "chocolate1"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$R_globostellata_01)
area2 <- sum(dat_ASV_marker_03_phylum$R_globostellata_02)
area3 <- sum(dat_ASV_marker_03_phylum$R_globostellata_03)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 0+203, n23 = 3+203, n13 = 6+203, 
                 n123 = 203, category = c("R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), lty = "blank", 
                 fill = c("coral3", "chocolate3", "chocolate1"))













####################################################################################################
# ONLY Transplant and autochtone
####################################################################################################

head(dat_ASV_04) 


dat_ASV_marker <- data.frame(dat_ASV_04$ASV,
                             dat_ASV_04$Cafard_island_autochtone_Rhabdastrella_globostellata_2019_V1V3,
                             dat_ASV_04$Cafard_island_transplant_Rhabdastrella_globostellata_2019_12weeks_V1V3)


colnames(dat_ASV_marker) <- c("ASV", "R_globostellata_autochtone", "R_globostellata_transplant")
head(dat_ASV_marker)
tail(dat_ASV_marker)

length(unique(dat_ASV_marker$ASV))

upset(dat_ASV_marker, sets = c("R_globostellata_autochtone", "R_globostellata_transplant"
), sets.bar.color = c("chartreuse", "forestgreen"),
order.by = "freq", empty.intersections = "on")

sum(dat_ASV_marker$R_globostellata_autochtone)
sum(dat_ASV_marker$R_globostellata_transplant)

draw.pairwise.venn(area1 = sum(dat_ASV_marker$R_globostellata_autochtone)-65, area2 = sum(dat_ASV_marker$R_globostellata_transplant)-65, cross.area = 65, category = c("R_globostellata autochtone", "R_globostellata transplant"))

draw.triple.venn(area1 = sum(dat_ASV_marker$R_globostellata_autochtone), area2 = sum(dat_ASV_marker$R_globostellata_transplant), 
                 
                 n12 = 60+822, n23 = 14+822, n13 = 10+822, 
                 n123 = 822, category = c("R_globostellata autochtone", "R_globostellata transplant"), lty = "blank", 
                 fill = c("chartreuse", "forestgreen"))


####################################################################################################
# by taxonomic drawer
####################################################################################################

head(dat_ASV)
#head(dat_ASV_marker_02_dupli)

link <- data.frame(dat_ASV$ASV, dat_ASV$kingdom, dat_ASV$phylum, dat_ASV$class, dat_ASV$order, dat_ASV$family, dat_ASV$genus, dat_ASV$species)
colnames(link) <- c("ASV", "kingdom", "phylum", "class", "order", "family", "genus", "species")
head(link)

link_02 <- link[!duplicated(link), ]
head(link_02)
tail(link_02)
length(link_02$ASV)

head(dat_ASV_marker)
tail(dat_ASV_marker)
length(unique(dat_ASV_marker$ASV))

dat_ASV_marker_03 <- merge(dat_ASV_marker, link_02, by="ASV",  all.x=T)
tail(dat_ASV_marker_03)
length(dat_ASV_marker_03)
length(unique(dat_ASV_marker_03$ASV))
length((dat_ASV_marker_03$ASV))

dat_ASV_marker_03_test <- dat_ASV_marker_03[!duplicated(dat_ASV_marker_03$ASV), ]
length(dat_ASV_marker_03_test$ASV)

dat_ASV_marker_03 <- dat_ASV_marker_03_test

head(dat_ASV_marker_03)
tail(dat_ASV_marker_03)
unique(dat_ASV_marker_03$phylum)

length(unique(dat_ASV_marker_03$ASV))

upset(dat_ASV_marker_03, sets = c("R_globostellata_autochtone", "R_globostellata_transplant"
), sets.bar.color = c("chartreuse", "forestgreen"),
order.by = "freq", empty.intersections = "on")

draw.pairwise.venn(area1 = sum(dat_ASV_marker_03$R_globostellata_autochtone), area2 = sum(dat_ASV_marker_03$R_globostellata_transplant), 
                   cross.area = 602, category = c("R_globostellata autochtone", "R_globostellata transplant"),
                   fill = c("chartreuse", "forestgreen"))

draw.pairwise.venn(area1 = sum(dat_ASV_marker_03$R_globostellata_autochtone), area2 = sum(dat_ASV_marker_03$R_globostellata_transplant), 
                   cross.area = 602, category = c("", ""),
                   fill = c("chartreuse", "forestgreen"))



### CHLOROFLEXI 

dat_ASV_marker_03_chloroflexi <- dat_ASV_marker_03[ which(dat_ASV_marker_03$phylum=='Chloroflexi'), ]
head(dat_ASV_marker_03_chloroflexi)
tail(dat_ASV_marker_03_chloroflexi)
dat_ASV_marker_03_chloroflexi_02 <- dat_ASV_marker_03_chloroflexi[!duplicated(dat_ASV_marker_03_chloroflexi), ]
head(dat_ASV_marker_03_chloroflexi_02)
tail(dat_ASV_marker_03_chloroflexi_02)

upset(dat_ASV_marker_03_chloroflexi_02, sets = c("R_globostellata_autochtone", "R_globostellata_transplant"
), sets.bar.color = c("chartreuse", "forestgreen"),
order.by = "freq", empty.intersections = "on")




draw.pairwise.venn(area1 = sum(dat_ASV_marker_03_chloroflexi_02$R_globostellata_autochtone), area2 = sum(dat_ASV_marker_03_chloroflexi_02$R_globostellata_transplant), 
                   cross.area = 171, category = c("R_globostellata autochtone", "R_globostellata transplant"),
                   fill = c("chartreuse", "forestgreen"))

draw.pairwise.venn(area1 = sum(dat_ASV_marker_03_chloroflexi_02$R_globostellata_autochtone), area2 = sum(dat_ASV_marker_03_chloroflexi_02$R_globostellata_transplant), 
                   cross.area = 171, category = c("", ""),
                   fill = c("chartreuse", "forestgreen"))




length(unique(dat_ASV_marker_03_chloroflexi_02$ASV))

### SAR202

head(dat_ASV_marker_03)

dat_ASV_marker_03_SAR202_clade <- dat_ASV_marker_03[ which(dat_ASV_marker_03$order=='SAR202_clade'), ]
head(dat_ASV_marker_03_SAR202_clade)
tail(dat_ASV_marker_03_SAR202_clade)
dat_ASV_marker_03_SAR202_clade_02 <- dat_ASV_marker_03_SAR202_clade[!duplicated(dat_ASV_marker_03_SAR202_clade), ]
head(dat_ASV_marker_03_SAR202_clade_02)
tail(dat_ASV_marker_03_SAR202_clade_02)

upset(dat_ASV_marker_03_SAR202_clade_02, sets = c("R_globostellata_autochtone", "R_globostellata_transplant"
), sets.bar.color = c("chartreuse", "forestgreen"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_SAR202_clade_02$R_globostellata_autochtone)
area2 <- sum(dat_ASV_marker_03_SAR202_clade_02$R_globostellata_transplant)
area3 <- sum(dat_ASV_marker_03_SAR202_clade_02$R_globostellata_03)

draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 0+111, n23 = 0+111, n13 = 0+111, 
                 n123 = 111, category = c("R_globostellata autochtone", "R_globostellata transplant"), lty = "blank", 
                 fill = c("chartreuse", "forestgreen"))


draw.pairwise.venn(area1 = sum(dat_ASV_marker_03_SAR202_clade_02$R_globostellata_autochtone), area2 = sum(dat_ASV_marker_03_SAR202_clade_02$R_globostellata_transplant), 
                   cross.area = 111, category = c("", ""),
                   fill = c("chartreuse", "forestgreen"))


length(unique(dat_ASV_marker_03_SAR202_clade_02$ASV))















################################################
## By phylum

head(dat_ASV_marker_03, n=100)
head(dat_ASV_marker_03)

dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$phylum), ]
head(dat_ASV_marker_03_phylum, n=100)
head(dat_ASV_marker_03_phylum)
length(unique(dat_ASV_marker_03_phylum$ASV))

upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant"
), sets.bar.color = c("chartreuse", "forestgreen"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$R_globostellata_autochtone)
area2 <- sum(dat_ASV_marker_03_phylum$R_globostellata_transplant)


draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 19+358, n23 = 14+358, n13 = 10+358, 
                 n123 = 358, category = c("R_globostellata autochtone", "R_globostellata transplant"), lty = "blank", 
                 fill = c("chartreuse", "forestgreen"))


################################################
## By class


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$class), ]


upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant"
), sets.bar.color = c("chartreuse", "forestgreen"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$R_globostellata_autochtone)
area2 <- sum(dat_ASV_marker_03_phylum$R_globostellata_transplant)


draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 19+337, n23 = 14+337, n13 = 6+337, 
                 n123 = 337, category = c("R_globostellata autochtone", "R_globostellata transplant"), lty = "blank", 
                 fill = c("chartreuse", "forestgreen"))


################################################
## By order


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$order), ]


upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant"
), sets.bar.color = c("chartreuse", "forestgreen"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$R_globostellata_autochtone)
area2 <- sum(dat_ASV_marker_03_phylum$R_globostellata_transplant)


draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 19+188, n23 = 14+188, n13 = 6+188, 
                 n123 = 188, category = c("R_globostellata_autochtone", "R_globostellata_transplant"), lty = "blank", 
                 fill = c("chartreuse", "forestgreen"))

################################################
## By family


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$family), ]


upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant"
), sets.bar.color = c("chartreuse", "forestgreen"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$R_globostellata_autochtone)
area2 <- sum(dat_ASV_marker_03_phylum$R_globostellata_transplant)


draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 9+173, n23 = 14+173, n13 = 6+173, 
                 n123 = 173, category = c("R_globostellata_autochtone", "R_globostellata_transplant"), lty = "blank", 
                 fill = c("chartreuse", "forestgreen"))

################################################
## By genus


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$genus), ]


upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant"
), sets.bar.color = c("chartreuse", "forestgreen"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$R_globostellata_autochtone)
area2 <- sum(dat_ASV_marker_03_phylum$R_globostellata_transplant)


draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 9+171, n23 = 14+171, n13 = 6+171, 
                 n123 = 171, category = c("R_globostellata_autochtone", "R_globostellata_transplant"), lty = "blank", 
                 fill = c("chartreuse", "forestgreen"))

################################################
## By species


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$species), ]


upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant"
), sets.bar.color = c("chartreuse", "forestgreen"),
order.by = "freq", empty.intersections = "on")

area1 <- sum(dat_ASV_marker_03_phylum$R_globostellata_autochtone)
area2 <- sum(dat_ASV_marker_03_phylum$R_globostellata_transplant)


draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, 
                 n12 = 9+171, n23 = 14+171, n13 = 6+171, 
                 n123 = 171, category = c("R_globostellata_autochtone", "R_globostellata_transplant"), lty = "blank", 
                 fill = c("chartreuse", "forestgreen"))











####################################################################################################
# ONLY Transplant and autochtone + 3 RHABDA
####################################################################################################

head(dat_ASV_04) 


dat_ASV_marker <- data.frame(dat_ASV_04$ASV,
                             dat_ASV_04$Cafard_island_autochtone_Rhabdastrella_globostellata_2019_V1V3,
                             dat_ASV_04$Cafard_island_transplant_Rhabdastrella_globostellata_2019_12weeks_V1V3,
                             dat_ASV_04$Bourake_Rhabdastrella_globostellata_2019_1_V1V3,
                             dat_ASV_04$Bourake_Rhabdastrella_globostellata_2019_2_V1V3,
                             dat_ASV_04$Bourake_Rhabdastrella_globostellata_2019_3_V1V3)


colnames(dat_ASV_marker) <- c("ASV", "R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03")
head(dat_ASV_marker)
tail(dat_ASV_marker)

length(unique(dat_ASV_marker$ASV))

upset(dat_ASV_marker, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"
), sets.bar.color = c("chartreuse", "forestgreen", "coral3", "chocolate3", "chocolate1"),
order.by = "freq", empty.intersections = "on")



####################################################################################################
# by taxonomic drawer
####################################################################################################

head(dat_ASV)
#head(dat_ASV_marker_02_dupli)

link <- data.frame(dat_ASV$ASV, dat_ASV$kingdom, dat_ASV$phylum, dat_ASV$class, dat_ASV$order, dat_ASV$family, dat_ASV$genus, dat_ASV$species)
colnames(link) <- c("ASV", "kingdom", "phylum", "class", "order", "family", "genus", "species")
head(link)

link_02 <- link[!duplicated(link), ]
head(link_02)
tail(link_02)
length(link_02$ASV)

head(dat_ASV_marker)
tail(dat_ASV_marker)
length(unique(dat_ASV_marker$ASV))

dat_ASV_marker_03 <- merge(dat_ASV_marker, link_02, by="ASV",  all.x=T)
tail(dat_ASV_marker_03)
length(dat_ASV_marker_03)
length(unique(dat_ASV_marker_03$ASV))
length((dat_ASV_marker_03$ASV))

dat_ASV_marker_03_test <- dat_ASV_marker_03[!duplicated(dat_ASV_marker_03$ASV), ]
length(dat_ASV_marker_03_test$ASV)

dat_ASV_marker_03 <- dat_ASV_marker_03_test

head(dat_ASV_marker_03)
tail(dat_ASV_marker_03)
unique(dat_ASV_marker_03$phylum)

length(unique(dat_ASV_marker_03$ASV))

upset(dat_ASV_marker_03, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), 
      sets.bar.color = c("coral3", "chocolate3", "chartreuse", "chocolate1", "forestgreen"),
      order.by = "freq", empty.intersections = "on")




upset(dat_ASV_marker_03, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), sets.bar.color = c("chartreuse", "forestgreen", "coral3", "chocolate3", "chocolate1"),
      keep.order = T, empty.intersections = "on")


### CHLOROFLEXI 

dat_ASV_marker_03_chloroflexi <- dat_ASV_marker_03[ which(dat_ASV_marker_03$phylum=='Chloroflexi'), ]
head(dat_ASV_marker_03_chloroflexi)
tail(dat_ASV_marker_03_chloroflexi)
dat_ASV_marker_03_chloroflexi_02 <- dat_ASV_marker_03_chloroflexi[!duplicated(dat_ASV_marker_03_chloroflexi), ]
head(dat_ASV_marker_03_chloroflexi_02)
tail(dat_ASV_marker_03_chloroflexi_02)

upset(dat_ASV_marker_03_chloroflexi_02, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), sets.bar.color = c("chartreuse", "forestgreen", "coral3", "chocolate3", "chocolate1"),
      order.by = "freq", empty.intersections = "on")


upset(dat_ASV_marker_03_chloroflexi_02, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), 
      sets.bar.color = c("chartreuse","coral3", "forestgreen", "chocolate3",  "chocolate1"),
      order.by = "freq", empty.intersections = "on")


length(unique(dat_ASV_marker_03_chloroflexi_02$ASV))

### SAR202

head(dat_ASV_marker_03)

dat_ASV_marker_03_SAR202_clade <- dat_ASV_marker_03[ which(dat_ASV_marker_03$order=='SAR202_clade'), ]
head(dat_ASV_marker_03_SAR202_clade)
tail(dat_ASV_marker_03_SAR202_clade)
dat_ASV_marker_03_SAR202_clade_02 <- dat_ASV_marker_03_SAR202_clade[!duplicated(dat_ASV_marker_03_SAR202_clade), ]
head(dat_ASV_marker_03_SAR202_clade_02)
tail(dat_ASV_marker_03_SAR202_clade_02)

upset(dat_ASV_marker_03_SAR202_clade_02, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), sets.bar.color = c("chartreuse", "forestgreen", "coral3", "chocolate3", "chocolate1"),
      order.by = "freq", empty.intersections = "on")

################################################
## By phylum

head(dat_ASV_marker_03, n=100)
head(dat_ASV_marker_03)

dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$phylum), ]
head(dat_ASV_marker_03_phylum, n=100)
head(dat_ASV_marker_03_phylum)
length(unique(dat_ASV_marker_03_phylum$ASV))

upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), 
      sets.bar.color = c( "coral3", "chocolate3", "chartreuse", "chocolate1","forestgreen"),
      order.by = "freq", empty.intersections = "on")


################################################
## By class


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$class), ]


upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), sets.bar.color = c("chartreuse", "forestgreen", "coral3", "chocolate3", "chocolate1"),
      order.by = "freq", empty.intersections = "on")

upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), 
      sets.bar.color = c( "coral3", "chocolate3", "chartreuse", "chocolate1","forestgreen"),
      order.by = "freq", empty.intersections = "on")
################################################
## By order


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$order), ]


upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), sets.bar.color = c("chartreuse", "forestgreen", "coral3", "chocolate3", "chocolate1"),
      order.by = "freq", empty.intersections = "on")

upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), 
      sets.bar.color = c( "coral3", "chocolate3", "chartreuse", "chocolate1","forestgreen"),
      order.by = "freq", empty.intersections = "on")


################################################
## By family


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$family), ]


upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), sets.bar.color = c("chartreuse", "forestgreen", "coral3", "chocolate3", "chocolate1"),
      order.by = "freq", empty.intersections = "on")

upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), 
      sets.bar.color = c( "coral3", "chocolate3", "chartreuse", "chocolate1","forestgreen"),
      order.by = "freq", empty.intersections = "on")


################################################
## By genus


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$genus), ]


upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), sets.bar.color = c("chartreuse", "forestgreen", "coral3", "chocolate3", "chocolate1"),
      order.by = "freq", empty.intersections = "on")

upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), 
      sets.bar.color = c( "coral3", "chocolate3", "chartreuse", "chocolate1","forestgreen"),
      order.by = "freq", empty.intersections = "on")


################################################
## By species


dat_ASV_marker_03_phylum <- dat_ASV_marker_03[!is.na(dat_ASV_marker_03$species), ]


upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), sets.bar.color = c("chartreuse", "forestgreen", "coral3", "chocolate3", "chocolate1"),
      order.by = "freq", empty.intersections = "on")

upset(dat_ASV_marker_03_phylum, sets = c("R_globostellata_autochtone", "R_globostellata_transplant", "R_globostellata_01", "R_globostellata_02", "R_globostellata_03"), 
      sets.bar.color = c( "coral3", "chocolate3", "chartreuse", "chocolate1","forestgreen"),
      order.by = "freq", empty.intersections = "on")











####################################################################################################
# JUNK PART WARNING
####################################################################################################


movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=T, sep=";" )
mutations <- read.csv( system.file("extdata", "mutations.csv", package = "UpSetR"), header=T, sep = ",")

upset(movies,attribute.plots=list(gridrows=60,plots=list(list(plot=scatter_plot, x="ReleaseDate", y="AvgRating"),
                                                         list(plot=scatter_plot, x="ReleaseDate", y="Watches"),list(plot=scatter_plot, x="Watches", y="AvgRating"),
                                                         list(plot=histogram, x="ReleaseDate")), ncols = 2))

upset(movies,attribute.plots=list(gridrows=60,plots=list(list(plot=scatter_plot, x="ReleaseDate", y="AvgRating"))))



upset(mutations, sets = c("PTEN", "TP53", "EGFR", "PIK3R1", "RB1"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")




ECM_02 <- data.frame(dat_ASV$ASV, (ECM$MAO_106_1+
                                ECM$MAO_106_2+
                                ECM$MAO_106_3+
                                ECM$MAO_106_4+
                                ECM$MAO_308_1+
                                ECM$MAO_308_2+
                                ECM$MAO_308_3+
                                ECM$MAO_308_4+
                                ECM$MAO_346_1+
                                ECM$MAO_346_2+
                                ECM$MAO_346_3+
                                ECM$MAO_355_1+
                                ECM$MAO_355_2+
                                ECM$MAO_355_3+
                                ECM$MAO_406_1+
                                ECM$MAO_406_2+
                                ECM$MAO_406_3+
                                ECM$MAO_406_4+
                                ECM$MAO_69_2+
                                ECM$MAO_69_3+
                                ECM$MAO_69_4), (
                                  ECM$MD_225_1+
                                    ECM$MD_225_2+
                                    ECM$MD_225_3+
                                    ECM$MD_225_4+
                                    ECM$MD_263_1+
                                    ECM$MD_263_2+
                                    ECM$MD_263_3+
                                    ECM$MD_263_4+
                                    ECM$MD_347_1+
                                    ECM$MD_347_2+
                                    ECM$MD_347_3+
                                    ECM$MD_358_2+
                                    ECM$MD_358_3+
                                    ECM$MD_358_4+
                                    ECM$MD_60_1+
                                    ECM$MD_60_2+
                                    ECM$MD_60_3+
                                    ECM$MD_60_4+
                                    ECM$MD_91_1+
                                    ECM$MD_91_2+
                                    ECM$MD_91_3+
                                    ECM$MD_91_4), (
                                      ECM$MP_135_2+
                                        ECM$MP_135_3+
                                        ECM$MP_135_4+
                                        ECM$MP_141_1+
                                        ECM$MP_141_2+
                                        ECM$MP_141_3+
                                        ECM$MP_141_4+
                                        ECM$MP_152_1+
                                        ECM$MP_152_2+
                                        ECM$MP_152_3+
                                        ECM$MP_152_4+
                                        ECM$MP_159_1+
                                        ECM$MP_159_2+
                                        ECM$MP_159_3+
                                        ECM$MP_159_4+
                                        ECM$MP_197_1+
                                        ECM$MP_197_2+
                                        ECM$MP_197_3+
                                        ECM$MP_197_4+
                                        ECM$MP_93_1+
                                        ECM$MP_93_2+
                                        ECM$MP_93_3 ), ( 
                                          ECM$MPGd_178_1+
                                            ECM$MPGd_178_2+
                                            ECM$MPGd_178_3+
                                            ECM$MPGd_178_4+
                                            ECM$MPGd_182_1+
                                            ECM$MPGd_182_2+
                                            ECM$MPGd_182_3+
                                            ECM$MPGd_182_4+
                                            ECM$MPGd_184_1+
                                            ECM$MPGd_184_2+
                                            ECM$MPGd_184_3+
                                            ECM$MPGd_184_4+
                                            ECM$MPGd_209_1+
                                            ECM$MPGd_209_2+
                                            ECM$MPGd_209_3+
                                            ECM$MPGd_209_4+
                                            ECM$MPGd_212_1+
                                            ECM$MPGd_212_2+
                                            ECM$MPGd_212_3+
                                            ECM$MPGd_212_4+
                                            ECM$MPGd_221_1+
                                            ECM$MPGd_221_2+
                                            ECM$MPGd_221_3+
                                            ECM$MPGd_221_4), (
                                              ECM$ZR_A11_1+
                                                ECM$ZR_A11_2+
                                                ECM$ZR_A11_3+
                                                ECM$ZR_A11_4), (
                                                  ECM$ZR_A21_1+
                                                    ECM$ZR_A21_2+
                                                    ECM$ZR_A21_3+
                                                    ECM$ZR_A21_4+
                                                    ECM$ZR_A22_1+
                                                    ECM$ZR_A22_2+
                                                    ECM$ZR_A22_3+
                                                    ECM$ZR_A22_4), (
                                                      ECM$ZR_A31_1+
                                                        ECM$ZR_A31_2+
                                                        ECM$ZR_A31_3), (
                                                          ECM$ZR_B11_1+
                                                            ECM$ZR_B11_2+
                                                            ECM$ZR_B11_3+
                                                            ECM$ZR_B11_4), (
                                                              ECM$ZR_C11_1+
                                                                ECM$ZR_C11_2+
                                                                ECM$ZR_C11_3+
                                                                ECM$ZR_C11_4+
                                                                ECM$ZR_C12_1+
                                                                ECM$ZR_C12_2+
                                                                ECM$ZR_C12_3+
                                                                ECM$ZR_C12_4))

head(ECM_02)
colnames(ECM_02) <- c("sp", "MAO", "MD", "MP", "MPGd", "A1","A2","A3", "B", "C")

str(ECM_02)

upset(mutations, sets = c("PTEN", "TP53", "EGFR", "PIK3R1", "RB1"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")

ECM_02
ECM_03 <- ECM_02 %>% mutate(MAO = replace(MAO, MAO != 0, 1))
ECM_03 <- ECM_03 %>% mutate(MD = replace(MD, MD != 0, 1))
ECM_03 <- ECM_03 %>% mutate(MP = replace(MP, MP != 0, 1))
ECM_03 <- ECM_03 %>% mutate(MPGd = replace(MPGd, MPGd != 0, 1))
ECM_03 <- ECM_03 %>% mutate(A1 = replace(A1, A1 != 0, 1))
ECM_03 <- ECM_03 %>% mutate(A2 = replace(A2, A2 != 0, 1))
ECM_03 <- ECM_03 %>% mutate(A3 = replace(A3, A3 != 0, 1))
ECM_03 <- ECM_03 %>% mutate(B = replace(B, B != 0, 1))
ECM_03 <- ECM_03 %>% mutate(C = replace(C, C != 0, 1))


head(ECM_03)
colnames(ECM_03) <- c("sp", "OLM", "CLM", "TCM.Fs", "TCM.Gd", "A1","A2","A3", "B", "C")
head(ECM_03)

upset(ECM_03, sets = c("OLM", "CLM", "TCM.Fs", "TCM.Gd", "A1","A2","A3", "B", "C"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")

upset(ECM_03, sets = c("OLM", "CLM", "TCM.Fs", "TCM.Gd", "A1","A2","A3", "B", "C"), 
      sets.bar.color = c("forestgreen","darkred",  "darkblue", "darkorchid4", "darkorange1", "darkorchid2", "darkorchid1","darkorchid3","darkorchid"),
      order.by = "freq", empty.intersections = "on")

head(ECM_03)





####################################################################################################
# Now cleaned the data to be readed by ggplot2
####################################################################################################

unique(dat$kingdom)

## I'm doing this in order to delete the "did'nt concatenated black borders" under a same fill (here genus).
dat_summarise <- dat %>% 
  group_by(phylum, index) %>% 
  summarise(value = sum(value))

head(dat_summarise)

####################################################################################################
# ANBUNDANCE
####################################################################################################

# Run the plot


ggplot(dat_summarise, aes(x = index, y = value, fill = phylum)) + 
  ylab("Abundance") + 
  xlab("Individuals - markers") + 
  geom_bar(stat = "identity", position="fill", colour="black")  +
  theme_bw() + 
  scale_fill_discrete(name = "Phyla") + 
  ggtitle(expression("Bacterial phyla abundance in the"~italic("R. globostellata")~ "and"~italic("Mycale sp.")~ "sponge, in silt, sediment and sea water")) +
  theme(axis.text.x = element_text(angle = -45, hjust=0))

ggsave("bacterial_phyla_abundance.pdf", width = 30, height = 20, units = "cm")

## Create Table
pb_long_phyla <- spread(dat_summarise,
                              key=phylum,
                              value=value)
write.csv(as.data.frame(pb_long_phyla), "bacterial_phyla_abundance.csv")


#### Now with merged sites :

head(dat)
length(dat$index)
str(dat$index)
head(dat$index)
unique(dat$index)

index <- c(as.character(dat$index))
head(index)
str(index)
unique(index)
index <- gsub("Bourake_Rhabdastrella_globostellata_2018_1[^. ]*_V1V3", "Bourake_Rhabdastrella_globostellata_2018_1_V1V3", index)
index <- gsub("Bourake_Rhabdastrella_globostellata_2018_2[^. ]*_V1V3", "Bourake_Rhabdastrella_globostellata_2018_2_V1V3", index)
index <- gsub("Bourake_Rhabdastrella_globostellata_2018_3[^. ]*_V1V3", "Bourake_Rhabdastrella_globostellata_2018_3_V1V3", index)

index <- as.factor(index)
str(index)
length(index)

pb_long_merged_sites <- data.frame(dat$index, dat$phylum, dat$value, index)
head(pb_long_merged_sites)
colnames(pb_long_merged_sites) <- c("Sites", "Phylum", "value", "Sites_merged")
head(pb_long_merged_sites)

# I'm doing this in order to delete the "did'nt concatenated black borders" under a same fill (here genus).
pb_long_merged_sites <- pb_long_merged_sites %>% 
  group_by(Sites_merged, Phylum) %>% 
  summarise(value = sum(value))

head(pb_long_merged_sites)

## Plot :
ggplot(pb_long_merged_sites, aes(x = Sites_merged, y = value, fill = Phylum)) + 
  ylab("Abundance") + 
  xlab("Individuals - markers") + 
  geom_bar(stat = "identity", position="fill", colour="black")  +
  ggtitle(expression("Bacterial phyla abundance in the"~italic("R. globostellata")~ "and"~italic("Mycale sp.")~ "sponge, in silt, sediment and sea water - merged individuals 2018")) +
  theme_bw() + 
  scale_fill_discrete(name = "Phyla")  + 
  theme(axis.text.x = element_text(angle = -45, hjust=0))

ggsave("bacterial_phyla_abundance_merged_individuals.pdf", width = 30, height = 20, units = "cm")

## Create Table
pb_long_merged_sites_table <- spread(pb_long_merged_sites,
                                     key=Sites_merged,
                                     value=value)
write.csv(as.data.frame(pb_long_merged_sites_table), "bacterial_phyla_abundance_merged_individuals.csv")



####################################################################################################
####################################################################################################
####################################################################################################
# ALL_DATA without cyanobacteria and unknown 
####################################################################################################
####################################################################################################
####################################################################################################

head(dat)

dat_clean <- dat[complete.cases(dat[ , 5]),]
head(dat_clean)

unique(dat$phylum)

dat_clean_02 <- dat_clean[!grepl("Cyanobacteria", dat_clean$phylum),]
head(dat_clean_02)
unique(dat_clean_02$phylum)

dat <- dat_clean_02

####################################################################################################
# Now cleaned the data to be readed by ggplot2
####################################################################################################

## I'm doing this in order to delete the "did'nt concatenated black borders" under a same fill (here genus).
dat_summarise <- dat %>% 
  group_by(phylum, index) %>% 
  summarise(value = sum(value))

head(dat_summarise)

####################################################################################################
# ANBUNDANCE
####################################################################################################

# Run the plot


ggplot(dat_summarise, aes(x = index, y = value, fill = phylum)) + 
  ylab("Abundance") + 
  xlab("Individuals - markers") + 
  geom_bar(stat = "identity", position="fill", colour="black")  +
  theme_bw() + 
  scale_fill_discrete(name = "Phyla") + 
  ggtitle(expression("Bacterial phyla abundance in the"~italic("R. globostellata")~ "and"~italic("Mycale sp.")~ "sponge, in silt, sediment and sea water without cyanobacteria and unknown")) +
  theme(axis.text.x = element_text(angle = -45, hjust=0))

ggsave("bacterial_phyla_abundance_without_cyanobacteria_and_unknown.pdf", width = 30, height = 20, units = "cm")

## Create Table
pb_long_phyla <- spread(dat_summarise,
                        key=phylum,
                        value=value)
write.csv(as.data.frame(pb_long_phyla), "bacterial_phyla_abundance_without_cyanobacteria_and_unknown.csv")


#### Now with merged sites :

head(dat)
length(dat$index)
str(dat$index)
head(dat$index)
unique(dat$index)

index <- c(as.character(dat$index))
head(index)
str(index)
unique(index)
index <- gsub("Bourake_Rhabdastrella_globostellata_2018_1[^. ]*_V1V3", "Bourake_Rhabdastrella_globostellata_2018_1_V1V3", index)
index <- gsub("Bourake_Rhabdastrella_globostellata_2018_2[^. ]*_V1V3", "Bourake_Rhabdastrella_globostellata_2018_2_V1V3", index)
index <- gsub("Bourake_Rhabdastrella_globostellata_2018_3[^. ]*_V1V3", "Bourake_Rhabdastrella_globostellata_2018_3_V1V3", index)

index <- as.factor(index)
str(index)
length(index)

pb_long_merged_sites <- data.frame(dat$index, dat$phylum, dat$value, index)
head(pb_long_merged_sites)
colnames(pb_long_merged_sites) <- c("Sites", "Phylum", "value", "Sites_merged")
head(pb_long_merged_sites)

# I'm doing this in order to delete the "did'nt concatenated black borders" under a same fill (here genus).
pb_long_merged_sites <- pb_long_merged_sites %>% 
  group_by(Sites_merged, Phylum) %>% 
  summarise(value = sum(value))

head(pb_long_merged_sites)

## Plot :
ggplot(pb_long_merged_sites, aes(x = Sites_merged, y = value, fill = Phylum)) + 
  ylab("Abundance") + 
  xlab("Individuals - markers") + 
  geom_bar(stat = "identity", position="fill", colour="black")  +
  ggtitle(expression("Bacterial phyla abundance in the"~italic("R. globostellata")~ "and"~italic("Mycale sp.")~ "sponge, in silt, sediment and sea water - merged individuals 2018 without cyanobacteria and unknown")) +
  theme_bw() + 
  scale_fill_discrete(name = "Phyla")  + 
  theme(axis.text.x = element_text(angle = -45, hjust=0))

ggsave("bacterial_phyla_abundance_merged_individuals_without_cyanobacteria_and_unknown.pdf", width = 30, height = 20, units = "cm")

## Create Table
pb_long_merged_sites_table <- spread(pb_long_merged_sites,
                                     key=Sites_merged,
                                     value=value)
write.csv(as.data.frame(pb_long_merged_sites_table), "bacterial_phyla_abundance_merged_individuals_without_cyanobacteria_and_unknown.csv")





####################################################################################################
####################################################################################################
####################################################################################################
# Focus on chloroflexi
####################################################################################################
####################################################################################################
####################################################################################################

head(dat)
unique(dat$phylum)

chloroflexi <- dat[ which(dat$phylum=='Chloroflexi'), ]
head(chloroflexi)

####################################################################################################
# Now cleaned the data to be readed by ggplot2
####################################################################################################

## I'm doing this in order to delete the "did'nt concatenated black borders" under a same fill (here genus).
chloroflexi_summarise <- chloroflexi %>% 
  group_by(order, index) %>% 
  summarise(value = sum(value))

head(chloroflexi_summarise)

####################################################################################################
# ANBUNDANCE
####################################################################################################

# Run the plot


ggplot(chloroflexi_summarise, aes(x = index, y = value, fill = order)) + 
  ylab("Abundance") + 
  xlab("Individuals - markers") + 
  geom_bar(stat = "identity", position="fill", colour="black")  +
  theme_bw() + 
  scale_fill_discrete(name = "Order") + 
  ggtitle(expression("Chloroflexi order abundance in the"~italic("R. globostellata")~ "and"~italic("Mycale sp.")~ "sponge, in silt, sediment and sea water")) +
  theme(axis.text.x = element_text(angle = -45, hjust=0))

ggsave("chloroflexi_order_abundance.pdf", width = 30, height = 20, units = "cm")

## Create Table
pb_long_phyla <- spread(chloroflexi_summarise,
                        key=order,
                        value=value)
write.csv(as.data.frame(pb_long_phyla), "chloroflexi_order_abundance.csv")


#### Now with merged sites :

head(chloroflexi)
length(chloroflexi$index)
str(chloroflexi$index)
head(chloroflexi$index)
unique(chloroflexi$index)

index <- c(as.character(chloroflexi$index))
head(index)
str(index)
unique(index)
index <- gsub("Bourake_Rhabdastrella_globostellata_2018_1[^. ]*_V1V3", "Bourake_Rhabdastrella_globostellata_2018_1_V1V3", index)
index <- gsub("Bourake_Rhabdastrella_globostellata_2018_2[^. ]*_V1V3", "Bourake_Rhabdastrella_globostellata_2018_2_V1V3", index)
index <- gsub("Bourake_Rhabdastrella_globostellata_2018_3[^. ]*_V1V3", "Bourake_Rhabdastrella_globostellata_2018_3_V1V3", index)

index <- as.factor(index)
str(index)
length(index)

pb_long_merged_sites <- data.frame(chloroflexi$index, chloroflexi$order, chloroflexi$value, index)
head(pb_long_merged_sites)
colnames(pb_long_merged_sites) <- c("Sites", "Phylum", "value", "Sites_merged")
head(pb_long_merged_sites)

# I'm doing this in order to delete the "did'nt concatenated black borders" under a same fill (here genus).
pb_long_merged_sites <- pb_long_merged_sites %>% 
  group_by(Sites_merged, Phylum) %>% 
  summarise(value = sum(value))

head(pb_long_merged_sites)

## Plot :
ggplot(pb_long_merged_sites, aes(x = Sites_merged, y = value, fill = Phylum)) + 
  ylab("Abundance") + 
  xlab("Individuals - markers") + 
  geom_bar(stat = "identity", position="fill", colour="black")  +
  ggtitle(expression("Chloroflexi order abundance in the"~italic("R. globostellata")~ "and"~italic("Mycale sp.")~ "sponge, in silt, sediment and sea water - merged individuals 2018")) +
  theme_bw() + 
  scale_fill_discrete(name = "Order")  + 
  theme(axis.text.x = element_text(angle = -45, hjust=0))

ggsave("chloroflexi_order_abundance_merged_individuals.pdf", width = 30, height = 20, units = "cm")

## Create Table
pb_long_merged_sites_table <- spread(pb_long_merged_sites,
                                     key=Sites_merged,
                                     value=value)
write.csv(as.data.frame(pb_long_merged_sites_table), "chloroflexi_order_abundance_merged_individuals.csv")


####################################################################################################
####################################################################################################
# Without NAs
####################################################################################################
####################################################################################################

head(chloroflexi)
unique(chloroflexi$order)


chloroflexi_clean <- chloroflexi[complete.cases(chloroflexi[ , 7]),]
head(chloroflexi_clean)
unique(chloroflexi_clean$order)

chloroflexi <- chloroflexi_clean

####################################################################################################
# Now cleaned the data to be readed by ggplot2
####################################################################################################

## I'm doing this in order to delete the "did'nt concatenated black borders" under a same fill (here genus).
chloroflexi_summarise <- chloroflexi %>% 
  group_by(order, index) %>% 
  summarise(value = sum(value))

head(chloroflexi_summarise)

####################################################################################################
# ANBUNDANCE
####################################################################################################

# Run the plot


ggplot(chloroflexi_summarise, aes(x = index, y = value, fill = order)) + 
  ylab("Abundance") + 
  xlab("Individuals - markers") + 
  geom_bar(stat = "identity", position="fill", colour="black")  +
  theme_bw() + 
  scale_fill_discrete(name = "Order") + 
  ggtitle(expression("Chloroflexi order abundance in the"~italic("R. globostellata")~ "and"~italic("Mycale sp.")~ "sponge, in silt, sediment and sea water without unknown")) +
  theme(axis.text.x = element_text(angle = -45, hjust=0))

ggsave("chloroflexi_order_abundance_without_unknown.pdf", width = 30, height = 20, units = "cm")

## Create Table
pb_long_phyla <- spread(chloroflexi_summarise,
                        key=order,
                        value=value)
write.csv(as.data.frame(pb_long_phyla), "chloroflexi_order_abundance_without_unknown.csv")


#### Now with merged sites :

head(chloroflexi)
length(chloroflexi$index)
str(chloroflexi$index)
head(chloroflexi$index)
unique(chloroflexi$index)

index <- c(as.character(chloroflexi$index))
head(index)
str(index)
unique(index)
index <- gsub("Bourake_Rhabdastrella_globostellata_2018_1[^. ]*_V1V3", "Bourake_Rhabdastrella_globostellata_2018_1_V1V3", index)
index <- gsub("Bourake_Rhabdastrella_globostellata_2018_2[^. ]*_V1V3", "Bourake_Rhabdastrella_globostellata_2018_2_V1V3", index)
index <- gsub("Bourake_Rhabdastrella_globostellata_2018_3[^. ]*_V1V3", "Bourake_Rhabdastrella_globostellata_2018_3_V1V3", index)

index <- as.factor(index)
str(index)
length(index)

pb_long_merged_sites <- data.frame(chloroflexi$index, chloroflexi$order, chloroflexi$value, index)
head(pb_long_merged_sites)
colnames(pb_long_merged_sites) <- c("Sites", "Phylum", "value", "Sites_merged")
head(pb_long_merged_sites)

# I'm doing this in order to delete the "did'nt concatenated black borders" under a same fill (here genus).
pb_long_merged_sites <- pb_long_merged_sites %>% 
  group_by(Sites_merged, Phylum) %>% 
  summarise(value = sum(value))

head(pb_long_merged_sites)

## Plot :
ggplot(pb_long_merged_sites, aes(x = Sites_merged, y = value, fill = Phylum)) + 
  ylab("Abundance") + 
  xlab("Individuals - markers") + 
  geom_bar(stat = "identity", position="fill", colour="black")  +
  ggtitle(expression("Chloroflexi order abundance in the"~italic("R. globostellata")~ "and"~italic("Mycale sp.")~ "sponge, in silt, sediment and sea water - merged individuals 2018 without unknown")) +
  theme_bw() + 
  scale_fill_discrete(name = "Order")  + 
  theme(axis.text.x = element_text(angle = -45, hjust=0))

ggsave("chloroflexi_order_abundance_merged_individuals_without_unknown.pdf", width = 30, height = 20, units = "cm")

## Create Table
pb_long_merged_sites_table <- spread(pb_long_merged_sites,
                                     key=Sites_merged,
                                     value=value)
write.csv(as.data.frame(pb_long_merged_sites_table), "chloroflexi_order_abundance_merged_individuals_without_unknown.csv")


