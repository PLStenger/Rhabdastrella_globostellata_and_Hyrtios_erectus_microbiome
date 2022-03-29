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
library(pheatmap)
library(gplots)

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

setwd("~/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_01_Rhabdastrella_globostellata_and_Mycale_microbiome/Rhabdastrella_globostellata_and_Mycale_microbiome/08_barplots")

taxo <- read.table("taxonomy.tsv", sep="\t", header=T)
head(taxo)
str(taxo)

level <- read.csv("level-7.csv",  header=T)
head(level)
unique(level$FullName)

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

length(taxo_split$kingdom)
length(pb_long_02$FullName)

dat <- cbind(pb_long_02, taxo_split)
dat$index<-gsub("_V1V3","",as.character(dat$index))
dat$index<-gsub("Bourake_","",as.character(dat$index))
dat$index<-gsub("_2019","",as.character(dat$index))

head(dat)

unique(dat$index)

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
  ggtitle(expression("Bacterial phyla abundance in the"~italic("R. globostellata")~ "and"~italic("Mycale sp.")~ "sponge and in sea water")) +
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
  xlab("Individuals") + 
  geom_bar(stat = "identity", position="fill", colour="black")  +
  theme_bw() + 
  scale_fill_discrete(name = "Phyla") + 
  ggtitle(expression("Bacterial phyla abundance in the"~italic("R. globostellata")~ "and"~italic("Mycale sp.")~ "sponge and in sea water without cyanobacteria and unknown")) +
  theme(axis.text.x = element_text(angle = -45, hjust=0))

ggsave("bacterial_phyla_abundance_without_cyanobacteria_and_unknown.pdf", width = 30, height = 20, units = "cm")

## Create Table
pb_long_phyla <- spread(dat_summarise,
                        key=phylum,
                        value=value)
write.csv(as.data.frame(pb_long_phyla), "bacterial_phyla_abundance_without_cyanobacteria_and_unknown.csv")


####################################################################################################
####################################################################################################
# bar plot sans NA et cyanobactéries,  seulement V1V3, sans 2018, sans YATé, sans Cafard pour Phylum
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

head(dat)

unique(dat$index)
dat_clean_02 <- dat[!grepl("*_V4", dat$index),]
head(dat_clean_02)
unique(dat_clean_02$index)
dat <- dat_clean_02

head(dat)

unique(dat$index)
dat_clean_02 <- dat[!grepl("*Chimeric*", dat$index),]
head(dat_clean_02)
unique(dat_clean_02$index)
dat <- dat_clean_02

head(dat)

unique(dat$index)
dat_clean_02 <- dat[!grepl("Cafard_island*", dat$index),]
head(dat_clean_02)
unique(dat_clean_02$index)
dat <- dat_clean_02

head(dat)

unique(dat$index)
dat_clean_02 <- dat[!grepl("Yate_*", dat$index),]
head(dat_clean_02)
unique(dat_clean_02$index)
dat <- dat_clean_02

head(dat, n=30)
head(dat)

length(unique(dat$FullName))

####################################################################################################
# Now cleaned the data to be readed by ggplot2
####################################################################################################

## I'm doing this in order to delete the "did'nt concatenated black borders" under a same fill (here genus).
dat_summarise <- dat %>% 
  group_by(phylum, index) %>% 
  summarise(value = sum(value))

head(dat_summarise)

unique(dat$index)
new_labels <- c(expression(~italic("Mycale sp.")~ " individual 1"),
                expression(~italic("Mycale sp.")~ " individual 2"),
                expression(~italic("Mycale sp.")~ " individual 3"),
                expression(~italic("Rhabdastrella globostellata")~ " 1"),
                expression(~italic("Rhabdastrella globostellata")~ " 2"),
                expression(~italic("Rhabdastrella globostellata")~ " 3"),
                "Sea water downstream",
                "Sea water upstream",
                "Sediment 1",
                "Sediment 2",
                "Sediment 3")

ggplot(dat_summarise, aes(x = index, y = value, fill = phylum)) + 
  ylab("Abundance") + 
  xlab("Conditions") + 
  geom_bar(stat = "identity", position="fill", colour="black")  +
  theme_bw() + 
  scale_fill_discrete(name = "Phyla") + 
  ggtitle(expression("Bacterial phyla abundance in the"~italic("R. globostellata")~ "and"~italic("Mycale sp.")~ "sponge and in sea water")) +
  theme(axis.text.x = element_text(angle = -45, hjust=0)) + 
  scale_x_discrete(labels= new_labels)

ggsave("Figure_01_bacterial_phyla_abundance.pdf", width = 25, height = 15, units = "cm")


## Create Table
pb_long_phyla <- spread(dat_summarise,
                        key=phylum,
                        value=value)
write.csv(as.data.frame(pb_long_phyla), "Figure_01_output_table_of_bacterial_phyla_abundance.csv")


####################################################################################################
# heat map
####################################################################################################

head(dat)
length(dat$index)
str(dat$index)
head(dat$index)
unique(dat$index)
length(unique(dat$FullName))

head(dat)

sub_dat <- data.frame(dat$index, dat$FullName, dat$value)
colnames(sub_dat) <- c("index", "FullName", "value")
head(sub_dat)

sub_dat_spread <- spread(sub_dat,
                         key=index,
                         value=value)

head(sub_dat_spread)
row.names(sub_dat_spread) <- sub_dat_spread$FullName
sub_dat_spread$FullName <- NULL
str(sub_dat_spread)

sub_dat_spread_matrix <- data.matrix(sub_dat_spread)

heatmap(sub_dat_spread_matrix, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))
heatmap(sub_dat_spread_matrix, col = bluered(100))
pheatmap(sub_dat_spread_matrix)


heatmap(sub_dat_spread_matrix, col = bluered(100), Colv=NA, Rowv=NA)
heatmap(sub_dat_spread_matrix, col = bluered(100))


heatmap.2(sub_dat_spread_matrix, scale="none", col=colorRampPalette(c("blue", "white", "red"))(1000), Colv=NA, trace="none", dendrogram ="row", breaks=seq(-5, 5, 0.01)) 



lwid=c(0.2,5) #make column of dendrogram and key very small and other colum very big 
lhei=c(0.2,5) #make row of key and other dendrogram very small and other row big. 

png(filename='test.png', width=800, height=750)
par(mar=c(7,4,4,2)+0.1) 

test <- heatmap.2(sub_dat_spread_matrix,col=redgreen(75), scale="row",
          key=TRUE, symkey=FALSE, density.info="none",cexRow=1,cexCol=1,margins=c(15,15),trace="none",srtCol=45)
graphics.off()

vo <- plot.new()

library(ggplot2)
library(cowplot)
plot_grid(test, vo, labels = c('A', 'B'), label_size = 12)


pdf("heatmap.pdf",width=20,height=40)
heatmap(sub_dat_spread_matrix, col = bluered(100))
dev.off()





#Heatmap parameters
dist1 <- "euclidean"
clust <- "average"
#ann_colors = list(Color = c(Forest= "forestgreen",Long_fallow= "chartreuse3", Short_fallow= "chartreuse4"))


## Edit body of pheatmap:::draw_colnames, customizing it to your liking
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))


pdf("heatmap.pdf",width=15,height=13)

pheatmap(sub_dat_spread_matrix,
         clustering_distance_cols = dist1, 
clustering_method = clust, 
cluster_rows=TRUE,
cluster_cols=TRUE, 
#annotation_col=df,
show_rownames=TRUE,
show_colnames=TRUE,
fontsize_row=5,
fontsize_col=8,
fontsize=8)

dev.off()



unique(dat$genus)
head(dat)
newdata <- dat[ which(dat$genus=='TK10' & dat$index=='Bourake_Rhabdastrella_globostellata_2019_2_V1V3'),]
newdata <- dat[ which(dat$genus=='AqS1' & dat$index=='Bourake_Rhabdastrella_globostellata_2019_2_V1V3'),]

newdata <- dat[ which(dat$genus=='AqS1'),]
newdata <- dat[ which(dat$genus=='TK10'),]





####################################################################################################
####################################################################################################
# heat map phylum
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

head(dat)
sub_dat <- data.frame(dat$index, dat$phylum, dat$value)
colnames(sub_dat) <- c("index", "phylum", "value")
tail(sub_dat)



dat_summarise <- sub_dat %>% 
  group_by(phylum, index) %>% 
  summarise(value = sum(value))

head(dat_summarise)

sub_dat_spread <- spread(dat_summarise,
                         key=index,
                         value=value)

head(sub_dat_spread)
sub_dat_spread <- data.frame(sub_dat_spread)
row.names(sub_dat_spread) <- sub_dat_spread$phylum
sub_dat_spread$phylum <- NULL
str(sub_dat_spread)

sub_dat_spread_matrix <- data.matrix(sub_dat_spread)
head(sub_dat_spread_matrix)

####################################################################################################
# heat map phylum MYCALE
####################################################################################################

sub_dat_spread_matrix_mycale <- sub_dat_spread_matrix[,1:3]
head(sub_dat_spread_matrix_mycale)
#colnames(sub_dat_spread_matrix_mycale) <- c("Mycale sp. individual 1", "Mycale sp. individual 2", "Mycale sp. individual 3")
colnames(sub_dat_spread_matrix_mycale) <- c("My_1", "My_2", "My_3")


pdf("Figure_03_heatmap_mycale_phylum_count.pdf",width=7,height=10)
heatmap.2(sub_dat_spread_matrix_mycale, col = bluered(100), 
          cexRow=1, cexCol=1,
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          margins=c(12,8))
dev.off()
write.csv(sub_dat_spread_matrix_mycale, "Figure_03_heatmap_mycale_phylum_count_table.csv")


pdf("Figure_03_heatmap_mycale_phylum_z_score.pdf",width=7,height=10)
heatmap.2(scale(sub_dat_spread_matrix_mycale), col = bluered(100), 
          cexRow=1, cexCol=1,
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          margins=c(12,8))
dev.off()
write.csv(scale(sub_dat_spread_matrix_mycale), "Figure_03_heatmap_mycale_phylum_z_score_table.csv")


sub_dat_spread_matrix_mycale
log <- log10(sub_dat_spread_matrix_mycale)
replace(log, log==-Inf, 0)


pdf("Figure_03_heatmap_mycale_phylum_log10.pdf",width=7,height=10)
heatmap.2(replace(log, log==-Inf, 0), col = bluered(100), 
          cexRow=1, cexCol=1,
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          margins=c(12,8))
dev.off()

write.csv(replace(log, log==-Inf, 0), "Figure_03_heatmap_mycale_phylum_log10_table.csv")


####################################################################################################
# heat map phylum Rhabdastrella globostellata
####################################################################################################

sub_dat_spread_matrix_rhabdastrella <- sub_dat_spread_matrix[,c(4,5, 6)]
head(sub_dat_spread_matrix_rhabdastrella)
#colnames(sub_dat_spread_matrix_rhabdastrella) <- c("Rhabdastrella globostellata 1", "Rhabdastrella globostellata 2", "Rhabdastrella globostellata 3")
colnames(sub_dat_spread_matrix_rhabdastrella) <- c("RG_1", "RG_2", "RG_3")


pdf("Figure_03_heatmap_rhabdastrella_phylum_count.pdf",width=7,height=10)
heatmap.2(sub_dat_spread_matrix_rhabdastrella, col = bluered(100), 
          cexRow=1, cexCol=1,
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          margins=c(12,8))
dev.off()
write.csv(sub_dat_spread_matrix_mycale, "Figure_03_heatmap_rhabdastrella_phylum_count_table.csv")


pdf("Figure_03_heatmap_rhabdastrella_phylum_z_score.pdf",width=7,height=10)
heatmap.2(scale(sub_dat_spread_matrix_rhabdastrella), col = bluered(100), 
          cexRow=1, cexCol=1,
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          margins=c(12,8))
dev.off()
write.csv(scale(sub_dat_spread_matrix_mycale), "Figure_03_heatmap_rhabdastrella_phylum_z_score_table.csv")


sub_dat_spread_matrix_rhabdastrella
log <- log10(sub_dat_spread_matrix_rhabdastrella)
replace(log, log==-Inf, 0)


pdf("Figure_03_heatmap_rhabdastrella_phylum_log10.pdf",width=7,height=10)
heatmap.2(replace(log, log==-Inf, 0), col = bluered(100), 
          cexRow=1, cexCol=1,
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          margins=c(12,8))
dev.off()

write.csv(replace(log, log==-Inf, 0), "Figure_03_heatmap_rhabdastrella_phylum_log10_table.csv")

#####################################################################################################
## heat map chloro
#####################################################################################################
#
#head(dat)
#sub_dat <- data.frame(dat$index, dat$phylum, dat$genus, dat$value)
#colnames(sub_dat) <- c("index", "phylum", "genus", "value")
#newdata <- sub_dat[ which(sub_dat$phylum=='Chloroflexi'), ]
#
#head(newdata)
#unique(newdata$genus)
#tail(sub_dat)
#
#newdata$phylum <- NULL
#
#dat_summarise <- newdata %>% 
#  group_by(genus, index) %>% 
#  summarise(value = sum(value))
#
#head(dat_summarise)
#
#sub_dat_spread <- spread(dat_summarise,
#                         key=index,
#                         value=value)
#
#head(sub_dat_spread)
#sub_dat_spread <- data.frame(sub_dat_spread)
#sub_dat_spread <- sub_dat_spread[-c(8, 9), ]
#row.names(sub_dat_spread) <- sub_dat_spread$genus
#sub_dat_spread$genus <- NULL
#str(sub_dat_spread)
#head(sub_dat_spread)
#
#sub_dat_spread_matrix <- data.matrix(sub_dat_spread)
#head(sub_dat_spread_matrix)
#
#heatmap(sub_dat_spread_matrix, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))
#heatmap(sub_dat_spread_matrix, col = bluered(100))
#
#
#pdf("heatmap.pdf",width=40,height=20)
#heatmap(sub_dat_spread_matrix, col = bluered(100))
#dev.off()
#
#
#sub_dat_spread_matrix <- sub_dat_spread_matrix[,-13:-19]
#
#head(sub_dat_spread_matrix)
#
#
#pdf("heatmap.pdf",width=20,height=5)
#heatmap(sub_dat_spread_matrix, col = bluered(100))
#dev.off()
#

####################################################################################################
# heat map chloro ORDER PAS GENUS
####################################################################################################

head(dat)

sub_dat <- data.frame(dat$index, dat$phylum, dat$order, dat$value)
colnames(sub_dat) <- c("index", "phylum", "order", "value")
newdata <- sub_dat[ which(sub_dat$phylum=='Chloroflexi'), ]

head(newdata)
unique(newdata$genus)
tail(sub_dat)

newdata$phylum <- NULL

dat_summarise <- newdata %>% 
  group_by(order, index) %>% 
  summarise(value = sum(value))

head(dat_summarise)

sub_dat_spread <- spread(dat_summarise,
                         key=index,
                         value=value)
head(sub_dat_spread)

sub_dat_spread <- data.frame(sub_dat_spread)
sub_dat_spread <- sub_dat_spread[-c(11), ]
row.names(sub_dat_spread) <- sub_dat_spread$order
#sub_dat_spread$order <- NULL
str(sub_dat_spread)
head(sub_dat_spread)

sub_dat_spread_matrix <- data.matrix(sub_dat_spread)
head(sub_dat_spread_matrix)

#heatmap(sub_dat_spread_matrix, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))

sub_dat_spread_matrix_02 <- sub_dat_spread_matrix[,-1]
sub_dat_spread_matrix_mycale <- sub_dat_spread_matrix_02[,1:3]

#pdf("heatmap.pdf",width=20,height=5)
#heatmap(sub_dat_spread_matrix_mycale, col = bluered(100))
#dev.off()

colnames(sub_dat_spread_matrix_mycale) <- c("Mycale sp. individual 1", "Mycale sp. individual 2", "Mycale sp. individual 3")


pdf("Figure_03_heatmap_mycale_order_count.pdf",width=7,height=10)
heatmap.2(sub_dat_spread_matrix_mycale, col = bluered(100), 
          cexRow=1, cexCol=1,
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          margins=c(12,8))
dev.off()
write.csv(sub_dat_spread_matrix_mycale, "Figure_03_heatmap_mycale_order_count_table.csv")


pdf("Figure_03_heatmap_mycale_order_z_score.pdf",width=7,height=10)
heatmap.2(scale(sub_dat_spread_matrix_mycale), col = bluered(100), 
          cexRow=1, cexCol=1,
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          margins=c(12,8))
dev.off()
write.csv(scale(sub_dat_spread_matrix_mycale), "Figure_03_heatmap_mycale_order_z_score_table.csv")


sub_dat_spread_matrix_mycale
log <- log10(sub_dat_spread_matrix_mycale)
replace(log, log==-Inf, 0)


pdf("Figure_03_heatmap_mycale_order_log10.pdf",width=7,height=10)
heatmap.2(replace(log, log==-Inf, 0), col = bluered(100), 
          cexRow=1, cexCol=1,
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          margins=c(12,8))
dev.off()

write.csv(replace(log, log==-Inf, 0), "Figure_03_heatmap_mycale_order_log10_table.csv")



####################################################################################################
# heat map order Rhabdastrella globostellata
####################################################################################################

sub_dat_spread_matrix_rhabdastrella <- sub_dat_spread_matrix[,c(5,6,7)]
head(sub_dat_spread_matrix_rhabdastrella)
colnames(sub_dat_spread_matrix_rhabdastrella) <- c("Rhabdastrella globostellata 1", "Rhabdastrella globostellata 2", "Rhabdastrella globostellata 3")


pdf("Figure_03_heatmap_rhabdastrella_order_count.pdf",width=7,height=10)
heatmap.2(sub_dat_spread_matrix_rhabdastrella, col = bluered(100), 
          cexRow=1, cexCol=1,
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          margins=c(12,8))
dev.off()
write.csv(sub_dat_spread_matrix_mycale, "Figure_03_heatmap_rhabdastrella_order_count_table.csv")


pdf("Figure_03_heatmap_rhabdastrella_order_z_score.pdf",width=7,height=10)
heatmap.2(scale(sub_dat_spread_matrix_rhabdastrella), col = bluered(100), 
          cexRow=1, cexCol=1,
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          margins=c(12,8))
dev.off()
write.csv(scale(sub_dat_spread_matrix_mycale), "Figure_03_heatmap_rhabdastrella_order_z_score_table.csv")


sub_dat_spread_matrix_rhabdastrella
log <- log10(sub_dat_spread_matrix_rhabdastrella)
replace(log, log==-Inf, 0)


pdf("Figure_03_heatmap_rhabdastrella_order_log10.pdf",width=7,height=10)
heatmap.2(replace(log, log==-Inf, 0), col = bluered(100), 
          cexRow=1, cexCol=1,
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          margins=c(12,8))
dev.off()

write.csv(replace(log, log==-Inf, 0), "Figure_03_heatmap_rhabdastrella_order_log10_table.csv")















#### Now with merged sites :

head(dat)
length(dat$index)
str(dat$index)
head(dat$index)
unique(dat$index)
length(unique(dat$FullName))

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
  xlab("Individuals") + 
  geom_bar(stat = "identity", position="fill", colour="black")  +
  theme_bw() + 
  scale_fill_discrete(name = "Order") + 
  ggtitle(expression("Chloroflexi order abundance in the"~italic("R. globostellata")~ "and"~italic("Mycale sp.")~ "sponge and in sea water")) +
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
  xlab("Individuals") + 
  geom_bar(stat = "identity", position="fill", colour="black")  +
  theme_bw() + 
  scale_fill_discrete(name = "Order") + 
  ggtitle(expression("Chloroflexi order abundance in the"~italic("R. globostellata")~ "and"~italic("Mycale sp.")~ "sponge and in sea water without unknown")) +
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






####################################################################################################
####################################################################################################
# bar plot sans NA et cyanobactéries,  seulement V1V3, sans 2018, sans YATé, sans Cafard pour order
####################################################################################################
####################################################################################################

head(pb_long_merged_sites)

unique(pb_long_merged_sites$Sites_merged)
pb_long_merged_sites_clean_02 <- pb_long_merged_sites[!grepl("*_V4", pb_long_merged_sites$Sites_merged),]
head(pb_long_merged_sites_clean_02)
unique(pb_long_merged_sites_clean_02$Sites_merged)
pb_long_merged_sites <- pb_long_merged_sites_clean_02

head(pb_long_merged_sites)

unique(pb_long_merged_sites$Sites_merged)
pb_long_merged_sites_clean_02 <- pb_long_merged_sites[!grepl("*Chimeric*", pb_long_merged_sites$Sites_merged),]
head(pb_long_merged_sites_clean_02)
unique(pb_long_merged_sites_clean_02$Sites_merged)
pb_long_merged_sites <- pb_long_merged_sites_clean_02

head(pb_long_merged_sites)

unique(pb_long_merged_sites$Sites_merged)
pb_long_merged_sites_clean_02 <- pb_long_merged_sites[!grepl("Cafard_island*", pb_long_merged_sites$Sites_merged),]
head(pb_long_merged_sites_clean_02)
unique(pb_long_merged_sites_clean_02$Sites_merged)
pb_long_merged_sites <- pb_long_merged_sites_clean_02

head(pb_long_merged_sites)

unique(pb_long_merged_sites$Sites_merged)
pb_long_merged_sites_clean_02 <- pb_long_merged_sites[!grepl("Yate_*", pb_long_merged_sites$Sites_merged),]
head(pb_long_merged_sites_clean_02)
unique(pb_long_merged_sites_clean_02$Sites_merged)
pb_long_merged_sites <- pb_long_merged_sites_clean_02

head(pb_long_merged_sites)

new_labels <- c(expression(~italic("Mycale sp.")~ " individual 1"),
                expression(~italic("Mycale sp.")~ " individual 2"),
                expression(~italic("Mycale sp.")~ " individual 3"),
                expression(~italic("Rhabdastrella globostellata")~ " 1"),
                expression(~italic("Rhabdastrella globostellata")~ " 2"),
                expression(~italic("Rhabdastrella globostellata")~ " 3"),
                "Sea water downstream",
                "Sea water upstream",
                "Sediment 1",
                "Sediment 2",
                "Sediment 3")

## Plot :
ggplot(pb_long_merged_sites, aes(x = Sites_merged, y = value, fill = Phylum)) + 
  ylab("Abundance") + 
  xlab("Conditions") + 
  geom_bar(stat = "identity", position="fill", colour="black")  +
  ggtitle(expression("Chloroflexi order abundance in the"~italic("R. globostellata")~ "and"~italic("Mycale sp.")~ "sponge and in sea water")) +
  theme_bw() + 
  scale_fill_discrete(name = "Order")  + 
  theme(axis.text.x = element_text(angle = -45, hjust=0)) + 
  scale_x_discrete(labels= new_labels)

ggsave("Figure_02_chloroflexi_order_abundance.pdf", width = 25, height = 15, units = "cm")

## Create Table
pb_long_merged_sites_table <- spread(pb_long_merged_sites,
                                     key=Sites_merged,
                                     value=value)
write.csv(as.data.frame(pb_long_merged_sites_table), "Figure_02_output_table_chloroflexi_order_abundance.csv")




