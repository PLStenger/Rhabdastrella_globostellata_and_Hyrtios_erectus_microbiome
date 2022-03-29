setwd("~/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_01_Rhabdastrella_globostellata_and_Mycale_microbiome/Rhabdastrella_globostellata_and_Mycale_microbiome/10_ICP_dbRDA")

############################################################################################################
######### Upload files and charge libraries
############################################################################################################

library(ggplot2)
library(ade4)
library(factoextra)
library(magrittr)
library(tidyverse)
library(vegan)
library(ggordiplots)
library(multcomp)
library(ggpubr)
library(gridExtra)
library(grid)

#pc <- read.table("physico_chemical.txt", header=T)
#pc <- read.table("Godd_data_ICP_for_all.txt", header=T)
#pc <- read.table("New_Data_Decembre_2021.txt", header=T)
#pc <- read.table("New_Data_Decembre_2021_all.txt", header=T, sep="\t")
pc <- read.table("New_Data_Decembre_2021_all_DOM.txt", header=T, sep="\t")
head(pc)
str(pc)

#pc$index <- gsub("SeaWater_01", "SeaWater_downstream_0_2micro", pc$index)
#pc$index <- gsub("SeaWater_02", "SeaWater_upstream_0_2micro", pc$index)

bpa <- read.csv("Figure_01_output_table_of_bacterial_phyla_abundance.csv", sep=",")
head(bpa)
bpa$index <- gsub("Mycale", "Dy", bpa$index)
bpa$index <- gsub("Rhabdastrella_globostellata", "RG", bpa$index)
bpa$index <- gsub("SeaWater_downstream_0_2micro", "SW_1", bpa$index)
bpa$index <- gsub("SeaWater_upstream_0_2micro", "SW_3", bpa$index)
bpa$index <- gsub("Sediment", "Sed", bpa$index)
head(bpa)
bpa
# If physico_chemical.txt
#bpa_02 <- bpa[ c(1, 2, 3, 5,8,11), ]
#bpa_02 <- bpa[ c(1, 2, 3, 5,8,11, 13, 14), ]

bpa_02 <- bpa

head(bpa_02)

#a <- gsub("_1", "_01", pc$index)
#b <- gsub("_2", "_02", a)
#c <- gsub("_3", "_03", b)
#pc$index <- c

#a <- gsub("ata_1", "ata_01", bpa_02$index)
#b <- gsub("ata_2", "ata_02", a)
#c <- gsub("ata_3", "ata_03", b)
#d <- gsub("ale_1", "ale_01", c)
#e <- gsub("ale_2", "ale_02", d)
#f <- gsub("ale_3", "ale_03", e)

#bpa_02$index <- NULL
bpa_02$X <- NULL
bpa_02$X.NA. <- NULL
#bpa_02$index <- f

head(bpa_02)
head(pc)

bpa_02
pc

dat <- merge(pc, bpa_02, by.x="index", by.y="index")
head(dat)

str(dat)

# Too much integer...
dat[] <- lapply(dat, function(x) {
  if(is.integer(x)) as.numeric(as.character(x)) else x
})

sapply(dat, class)
str(dat) # ok now..

head(dat)

row.names(dat) <- c("Dy_1",  "Dy_2",  "Dy_3", 
                    "RG_1" , "RG_2"  ,"RG_3" , 
                    "Sed_1", "Sed_2", "Sed_3",
                    "SW_1"  ,"SW_3")

############################################################################################################
### PCA en fonction PHYSICO-CHIMIC
############################################################################################################

head(dat)

unique(dat$index)

#dat2 <- dat[, -1]
#dat2$Cyanobacteria <- NULL
#dat2 <- dat[, c(2:11) ]
dat2 <- dat[, c(2:13) ]
#dat2 <- dat[, c(2:20) ]
#dat2 <- dat[, c(11:20) ]

head(dat2)


#dat2 <- dat2[ -c(7:8), ]

res.pca <- dudi.pca(dat2,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5            # Nombre d'axes gardés
)

str(dat)

dat$index  <- factor(dat$index, levels = c("Dy_1",  "Dy_2",  "Dy_3", 
                                           "RG_1" , "RG_2"  ,"RG_3" , 
                                           "Sed_1", "Sed_2", "Sed_3",
                                           "SW_1"  ,"SW_3"))
#unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "all", 
                col.ind = c("Dy_1",  "Dy_2",  "Dy_3", 
                            "RG_1" , "RG_2"  ,"RG_3" , 
                            "Sed_1", "Sed_2", "Sed_3",
                            "SW_1"  ,"SW_3"),
                pointsize = 2,
                repel=T,
                #repel=TRUE,
                #xlim=c(-10, 10) ,
                #ylim=c(-10, 10) ,
                title = "PCA - Biplot for physico-chemical data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod", 
                            "burlywood4", "antiquewhite4", "bisque4",
                            "cadetblue3", "cadetblue3"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")


summary(res.pca)


fviz_pca_biplot(res.pca, 
                label = "ind", 
                invisible ="var",
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Dy_1",  "Dy_2",  "Dy_3", 
                            "RG_1" , "RG_2"  ,"RG_3" , 
                            "Sed_1", "Sed_2", "Sed_3",
                            "SW_1"  ,"SW_3"),
                pointsize = 2,
                #xlim=c(-10, 10) ,
                #ylim=c(-10, 10) ,
                title = "PCA for physico-chemical data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod", 
                            "burlywood4", "antiquewhite4", "bisque4",
                            "cadetblue3", "cadetblue3"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")


options(ggrepel.max.overlaps = Inf)


fviz_pca_biplot(res.pca, 
                label = "var", 
                invisible ="ind", 
                #geom.var = c("point", "text"),
                alpha.var=0.3,
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Dy_1",  "Dy_2",  "Dy_3", 
                            "RG_1" , "RG_2"  ,"RG_3" , 
                            "Sed_1", "Sed_2", "Sed_3",
                            "SW_1"  ,"SW_3"),
                pointsize = 2,
                #xlim=c(-4, 4) ,
                #ylim=c(-4, 4) ,
                title = "Contributive arrows for physico-chemical data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod", 
                            "burlywood4", "antiquewhite4", "bisque4",
                            "cadetblue3", "cadetblue3"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")



############################################################################################################
### PCA en fonction TAXO
############################################################################################################

head(dat)

#dat2 <- dat[, -1]
#dat2 <- dat[, -c(1:25) ]
dat2 <- dat[, -c(1:13) ]
dat2$Cyanobacteria <- NULL

head(dat2)


res.pca <- dudi.pca(dat2,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5            # Nombre d'axes gardés
)

str(dat)

#dat$index  <- factor(dat$index, levels = c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03"))
#unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "all", 
                col.ind = c("Dy_1",  "Dy_2",  "Dy_3", 
                            "RG_1" , "RG_2"  ,"RG_3" , 
                            "Sed_1", "Sed_2", "Sed_3",
                            "SW_1"  ,"SW_3"),
                pointsize = 2,
                #xlim=c(-5, 5) ,
                #ylim=c(-5, 5) ,
                title = "PCA - Biplot for phyla data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod", 
                            "burlywood4", "antiquewhite4", "bisque4",
                            "cadetblue3", "cadetblue3"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

fviz_pca_biplot(res.pca, 
                label = "ind", 
                invisible ="var", 
                # geom=c("text", "point"),
                col.ind = c("Dy_1",  "Dy_2",  "Dy_3", 
                            "RG_1" , "RG_2"  ,"RG_3" , 
                            "Sed_1", "Sed_2", "Sed_3",
                            "SW_1"  ,"SW_3"),
                pointsize = 2,
                #xlim=c(-5, 5) ,
                #ylim=c(-5, 5) ,
                title = "PCA for phyla data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod", 
                            "burlywood4", "antiquewhite4", "bisque4",
                            "cadetblue3", "cadetblue3"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

options(ggrepel.max.overlaps = Inf)



fviz_pca_biplot(res.pca, 
                label = "var", 
                invisible ="ind", 
                #geom.var = c("point", "text"),
                alpha.var=0.3,
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Dy_1",  "Dy_2",  "Dy_3", 
                            "RG_1" , "RG_2"  ,"RG_3" , 
                            "Sed_1", "Sed_2", "Sed_3",
                            "SW_1"  ,"SW_3"),
                pointsize = 2,
                #xlim=c(-4, 4) ,
                #ylim=c(-4, 4) ,
                title = "Contributive arrows for phyla data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod", 
                            "burlywood4", "antiquewhite4", "bisque4",
                            "cadetblue3", "cadetblue3"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

















############################################################################################################
######### db-RDA Par Phylum
############################################################################################################

# If data set without seawater
#taxa_bacteria <- dat[, -c(1:25) ]

taxa_bacteria <- dat[, -c(1:13) ]


head(dat)

# If data set without seawater
#physico_chimic <- dat[, c(2:25) ]
physico_chimic <- dat[, c(2:13) ]
#physico_chimic <- dat[, c(11:20) ]



head(physico_chimic)
head(taxa_bacteria)

######### db-RDA for bacteria
############################################################################################################

# Beta Diversity
# from https://www.rpubs.com/roalle/mres_2019

# Bray-Curtis Dissimilarity

# Calculate pairwise dissimilarity
taxa_bacteria.mdf <- as.matrix.data.frame(taxa_bacteria)
#rownames(taxa_bacteria.mdf) <- env.data$site

taxa_bacteria.bray <- vegdist(taxa_bacteria.mdf, method = "bray")
taxa_bacteria.bray

# PERMANOVA (permutational multivariate analysis of variance; Anderson 2001)
# Homogeneity of dispersion test
#permutest(betadisper(taxa_bacteria.bray, physico_chimic$ph_h2o))

# PERMANOVA analysis
# adonis(taxa_bacteria.bray ~ ph_h2o + ph_kcl + N_NO3_mg_kg + N_NH4_mg_kg + N_g_kg + C_g_kg + P_mg_kg + cations_change_Ca + cations_change_Mg + cations_change_Na + cations_change_K + change_capacity, data = physico_chimic, permutations = 999)

# Distance-Based Redundancy Analysis (dbRDA)
# from https://www.rpubs.com/roalle/mres_2019

# Distance-Based Redundancy Analysis (dbRDA; Legendre and Anderson, 1999) is an extension of Redundancy Analysis (RDA) which allows the use of 
# non-euclidean distance matrices as inputs (e.g. Bray-Curtis dissimilarity). The method works by first calculating a PCoA from the dissimilarity matrix, 
# and then subjecting the PCoA eigenvalues (which represent dissimilarities in euclidean space) to RDA. The method aims to detect linear relationships between 
# environmental variables and these dissimilarities. 
# dbRDA differs from PCoA in that it is a constrained analysis. 
# While PCoA axes are generated to explain maximum variation in the distance matrix, 
# dbRDA canonical axes are constructed as linear combinations of environmental variables. 
# Consequently these axes are constrained to the environmental variables in the model, and the ordination will be distinct from the PCoA. dbRDA allows us to 
# visualise how environmental variables contrain variation in community composition between our sites. 
# As dbRDA considers multiple environmental variables which are measured using different techniques, it is important to normalise environmental variables such 
# that they can be compared concurrently. This is achieved using a z-score transformation. Moreover, dbRDA is sensitive to multicolinearity 
# (i.e. high correlation between environmental variables). We can test whether any of our environmental variables are correlated by calculating variance 
# inflation factors (VIFs) and then removing selected terms from the model until all VIF scores are < 10.

# Z-Score Transformation
# Z-Score transform environmental data :
# The z-score transformation normalises environmental variables by their standard deviation from the mean 
env.data.z <- physico_chimic

#env.data.z$Ca <- (env.data.z$Ca - mean(env.data.z$Ca))/sd(env.data.z$Ca)
#env.data.z$Mg <- (env.data.z$Mg - mean(env.data.z$Mg))/sd(env.data.z$Mg)
#env.data.z$Na <- (env.data.z$Na - mean(env.data.z$Na))/sd(env.data.z$Na)
#env.data.z$K <- (env.data.z$K - mean(env.data.z$K))/sd(env.data.z$K)
#env.data.z$P <- (env.data.z$P - mean(env.data.z$P))/sd(env.data.z$P)
#env.data.z$Co <- (env.data.z$Co - mean(env.data.z$Co))/sd(env.data.z$Co)
env.data.z$Cr <- (env.data.z$Cr - mean(env.data.z$Cr))/sd(env.data.z$Cr)
env.data.z$Fe <- (env.data.z$Fe - mean(env.data.z$Fe))/sd(env.data.z$Fe)
env.data.z$Mn <- (env.data.z$Mn - mean(env.data.z$Mn))/sd(env.data.z$Mn)
env.data.z$Ni <- (env.data.z$Ni - mean(env.data.z$Ni))/sd(env.data.z$Ni)
#env.data.z$Al <- (env.data.z$Al - mean(env.data.z$Al))/sd(env.data.z$Al)
#env.data.z$Cu <- (env.data.z$Cu - mean(env.data.z$Cu))/sd(env.data.z$Cu)
env.data.z$Zn <- (env.data.z$Zn - mean(env.data.z$Zn))/sd(env.data.z$Zn)
#env.data.z$Mo <- (env.data.z$Mo - mean(env.data.z$Mo))/sd(env.data.z$Mo)
#env.data.z$pH <- (env.data.z$pH - mean(env.data.z$pH))/sd(env.data.z$pH)
env.data.z$SiOH <- (env.data.z$SiOH - mean(env.data.z$SiOH))/sd(env.data.z$SiOH)
env.data.z$Nox <- (env.data.z$Nox - mean(env.data.z$Nox))/sd(env.data.z$Nox)
env.data.z$PO4 <- (env.data.z$PO4 - mean(env.data.z$PO4))/sd(env.data.z$PO4)
env.data.z$NH4 <- (env.data.z$NH4 - mean(env.data.z$NH4))/sd(env.data.z$NH4)
#env.data.z$DIC <- (env.data.z$DIC - mean(env.data.z$DIC))/sd(env.data.z$DIC)
env.data.z$POC <- (env.data.z$POC - mean(env.data.z$POC))/sd(env.data.z$POC)
env.data.z$DOC <- (env.data.z$DOC - mean(env.data.z$DOC))/sd(env.data.z$DOC)
env.data.z$PON <- (env.data.z$PON - mean(env.data.z$PON))/sd(env.data.z$PON)
#env.data.z$PN <- (env.data.z$PN - mean(env.data.z$PN))/sd(env.data.z$PN)
#env.data.z$Temp <- (env.data.z$Temp - mean(env.data.z$Temp))/sd(env.data.z$Temp)
#env.data.z$DIN <- (env.data.z$DIN - mean(env.data.z$DIN))/sd(env.data.z$DIN)

# Perform dbRDA
# construct full model and calculate VIF

dbRDA.full <- capscale(taxa_bacteria.bray ~   Mn +
                         Ni +
                         Cr +
                         SiOH +
                         Nox +
                         PO4 +
                         NH4 +
                         POC +
                         DOC +
                         PON +
                         Fe +
                         Zn ,
                       env.data.z)
vif.cca(dbRDA.full)

dbRDA.mat<- capscale(taxa_bacteria.bray ~   Mn +
                       Ni +
                       Cr +
                       SiOH +
                       Nox +
                       PO4 +
                       NH4 +
                       POC +
                       DOC +
                       PON +
                       Fe +
                       Zn ,
                     env.data.z)

dbRDA.mat <- capscale(taxa_bacteria.bray ~  POC +
                        DOC,
                      env.data.z)
vif.cca(dbRDA.mat)
# VIF explaination, according http://jsnfri.fra.affrc.go.jp/gunshu/3dbrdaE.html
# Multiple-colinearlity between variables are examined by variance inflation factor (VIF). 
# When VIF > 10, problems of multiple-colinearlity become severe. So, remove variables of large VIF, one by one.

# construct reduced model and calculate VIF
dbRDA.mat <- capscale(taxa_bacteria.mdf ~   Mn +
                        Ni +
                        Cr +
                        Fe +
                        Zn ,
                      env.data.z)



# test overall significance of the analysis
anova(dbRDA.mat, permutations = 999)

# test significance of each environmental variable
anova(dbRDA.mat, by = "terms")

# summary of dbRDA model to extract total variance constrained and axis scores
summary(dbRDA.mat)

# Plot dbRDA

smry <- summary(dbRDA.mat)
scrs <- scores(dbRDA.mat)
df1  <- data.frame(smry$sites[,1:2]) # site scores for RDA1 and RDA2
df1$site <- rownames(df1)  #add site names
df2  <- data.frame(smry$biplot[,1:2])  # mapping environmental variables

site <- c("Dy", "Dy", "Dy", "RG", "RG", "RG", "SW", "SW",  "Sed", "Sed", "Sed")
ind <- c("Dy_1", "Dy_2", "Dy_3", "RG_1", "RG_2", "RG_3", "SW_1", "SW_3",  "Sed_1", "Sed_2", "Sed_3")
df_col <- data.frame(ind, site)
df_col$site <- factor(df_col$site, levels = c("Dy", "RG", "SW", "Sed"))
head(df_col)



legend_title <- "Conditions"

ggplot(df1, aes(x=CAP1, y=CAP2, colour = df_col$site, shape=df_col$site)) + 
  geom_point(size = 3) +
  ggtitle("Distance-Based Redundancy Analysis") +
  scale_colour_manual(legend_title,labels=c("Dy", "RG", "SW", "Sed"), values = c("coral3",  "darkgoldenrod4", "cadetblue3", "bisque4")) +
  scale_shape_manual(legend_title,labels=c("Dy", "RG", "SW", "Sed"), values = c(16, 16, 16, 16)) +
  scale_fill_manual(legend_title,labels=c("Dy", "RG", "SW", "Sed"), values = c("coral3",  "darkgoldenrod4", "cadetblue3", "bisque4")) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #stat_ellipse() +
  # #xlim(-50, 50) +
  # #ylim(-50, 50) +
  xlab("db-RDA1 (75.36%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("db-RDA2 (13.32%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  coord_fixed() +
  # ADD ggforce's ellipses
  ggforce::geom_mark_ellipse(aes(fill = df_col$site,
                                 color = df_col$site)) +
  # guides(fill=guide_legend("my awesome title")) +
  theme_bw() #+ theme(legend.position="none")









#install.packages("ggrepel")
library(ggrepel)

ggplot(df1, aes(x=CAP1, y=CAP2, colour = df_col$site)) + 
  geom_segment(data=df2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               # color="black", arrow=arrow(length=unit(0.01,"npc"))) +
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  ggrepel:: geom_text_repel(data=df2, aes(x=CAP1,y=CAP2,label=rownames(df2)), 
                            #    color="grey50", size=3) +
                            color="black", size=5) +
  #  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  # scale_colour_manual(values = c( "firebrick", "goldenrod", "olivedrab1", "green", "forestgreen")) +
  # geom_text(aes(label=rownames(df1),
  #                hjust=0,vjust=1.5), colour = "black",size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #  stat_ellipse() +
  #xlim(-2, 1.1) +
  # #ylim(-5, 6) +
  xlab("RDA1 (75.36%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("RDA2 (13.32%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  #coord_fixed() +
  theme_bw()



ggplot(df1, aes(x=CAP1, y=CAP2, colour = df_col$site))+#, shape=df_col$site)) + 
  geom_point(size = 3) +
 ggforce::geom_mark_ellipse(aes(fill = df_col$site,
                                 color = df_col$site)) +
  geom_segment(data=df2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  ggrepel:: geom_text_repel(data=df2, aes(x=CAP1,y=CAP2,label=rownames(df2)), 
                            #    color="grey50", size=3) +
                            color="black", size=5) 





ggplot(df1, aes(x=CAP1, y=CAP2, colour = df_col$site))+#, shape=df_col$site)) + 
  geom_point(size = 3) +
  ggforce::geom_mark_ellipse(aes(fill = df_col$site,
                                 color = df_col$site)) +
  geom_segment(data=df2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  ggrepel:: geom_text_repel(data=df2, aes(x=CAP1,y=CAP2,label=rownames(df2)), 
                            #    color="grey50", size=3) +
                            color="black", size=5) +
  ggtitle("Distance-Based Redundancy Analysis") +
  scale_colour_manual(legend_title,labels=c("Dy", "RG", "SW", "Sed"), values = c("coral3",  "darkgoldenrod4", "cadetblue3", "bisque4")) +
  scale_shape_manual(legend_title,labels=c("Dy", "RG", "SW", "Sed"), values = c(16, 16, 16, 16)) +
  scale_fill_manual(legend_title,labels=c("Dy", "RG", "SW", "Sed"), values = c("coral3",  "darkgoldenrod4", "cadetblue3", "bisque4")) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #stat_ellipse() +
  # #xlim(-50, 50) +
  # #ylim(-50, 50) +
  xlab("db-RDA1 (75.36%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("db-RDA2 (13.32%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  coord_fixed() +
  # ADD ggforce's ellipses
  ggforce::geom_mark_ellipse(aes(fill = df_col$site,
                                 color = df_col$site)) +
  # guides(fill=guide_legend("my awesome title")) +
  theme_bw() #+ theme(legend.position="none")










############################################################################################################
############################################################################################################
######### QUE CHLOROFLEXI
############################################################################################################
############################################################################################################

library(ggplot2)
library(ade4)
library(factoextra)
library(magrittr)
library(tidyverse)
library(vegan)
library(ggordiplots)
library(multcomp)
library(ggpubr)
library(gridExtra)
library(grid)

#pc <- read.table("physico_chemical.txt", header=T)
#pc <- read.table("physico_chemical_02.txt", header=T)
#pc <- read.table("Godd_data_ICP_for_all.txt", header=T)

head(pc)


#pc$index <- gsub("SeaWater_01", "SeaWater_downstream_0_2micro", pc$index)
#pc$index <- gsub("SeaWater_02", "SeaWater_upstream_0_2micro", pc$index)



bpa <- read.csv("Figure_02_output_table_chloroflexi_order_abundance.csv", sep=";")
head(bpa)
bpa$index <- gsub("Mycale", "Dy", bpa$index)
bpa$index <- gsub("Rhabdastrella_globostellata", "RG", bpa$index)
bpa$index <- gsub("SeaWater_downstream_0_2micro", "SW_1", bpa$index)
bpa$index <- gsub("SeaWater_upstream_0_2micro", "SW_3", bpa$index)
bpa$index <- gsub("Sediment", "Sed", bpa$index)
head(bpa)
bpa
# If physico_chemical.txt
#bpa_02 <- bpa[ c(1, 2, 3, 5,8,11), ]
bpa_02 <- bpa

head(bpa_02)

#a <- gsub("ata_1", "ata_01", bpa_02$index)
#b <- gsub("ata_2", "ata_02", a)
#c <- gsub("ata_3", "ata_03", b)
#d <- gsub("ale_1", "ale_01", c)
#e <- gsub("ale_2", "ale_02", d)
#f <- gsub("ale_3", "ale_03", e)

#bpa_02$index <- NULL
bpa_02$X <- NULL
bpa_02$X.NA. <- NULL
#bpa_02$index <- f

head(bpa_02)
head(pc)

dat <- merge(pc, bpa_02, by.x="index", by.y="index")
head(dat)

str(dat)

# Too much integer...
dat[] <- lapply(dat, function(x) {
  if(is.integer(x)) as.numeric(as.character(x)) else x
})

sapply(dat, class)
str(dat) # ok now..



row.names(dat) <- c("Dy_1", "Dy_2", "Dy_3", 
                    "RG_1", "RG_2", "RG_3", 
                    "Sed_1", "Sed_2", "Sed_3",
                    "SW_1", "SW_2")

############################################################################################################
### PCA en fonction PHYSICO-CHIMIC
############################################################################################################

#dat2 <- dat[, -1]
#dat2 <- dat[, c(2:25) ]
dat2 <- dat[, c(2:13) ]


res.pca <- dudi.pca(dat2,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5            # Nombre d'axes gardés
)

str(dat)

#dat$index  <- factor(dat$index, levels = c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03"))
#unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "all", 
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3", 
                            "Sed_1", "Sed_2", "Sed_3",
                            "SW_1", "SW_2"),                pointsize = 2,
                repel=T,
                #repel=TRUE,
                #xlim=c(-10, 10) ,
                #ylim=c(-10, 10) ,
                title = "PCA - Biplot for physico-chemical data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod", 
                            "burlywood4", "antiquewhite4", "bisque4",
                            "cadetblue3", "cadetblue3"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

summary(res.pca)

dat$index  <- factor(dat$index, levels = c("Dy_1", "Dy_2", "Dy_3", 
                                           "RG_1", "RG_2", "RG_3", 
                                           "Sed_1", "Sed_2", "Sed_3",
                                           "SW_1", "SW_2"))
unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "ind", 
                invisible ="var",
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3", 
                            "Sed_1", "Sed_2", "Sed_3",
                            "SW_1", "SW_2"),                pointsize = 2,
                #xlim=c(-10, 10) ,
                #ylim=c(-10, 10) ,
                title = "PCA for physico-chemical data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod", 
                            "burlywood4", "antiquewhite4", "bisque4",
                            "cadetblue3", "cadetblue3"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

options(ggrepel.max.overlaps = Inf)



fviz_pca_biplot(res.pca, 
                label = "var", 
                invisible ="ind", 
                #geom.var = c("point", "text"),
                alpha.var=0.3,
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3", 
                            "Sed_1", "Sed_2", "Sed_3",
                            "SW_1", "SW_2"),                pointsize = 2,
                #xlim=c(-4, 4) ,
                #ylim=c(-4, 4) ,
                title = "Contributive arrows for physico-chemical data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod", 
                            "burlywood4", "antiquewhite4", "bisque4",
                            "cadetblue3", "cadetblue3"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")



############################################################################################################
### PCA en fonction TAXO
############################################################################################################

#dat2 <- dat[, -1]
#dat2 <- dat[, -c(1:25) ]
dat2 <- dat[, -c(1:13) ]
dat2$Cyanobacteria <- NULL

#row.names(dat2) <- c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03", "SeaWater 01", "SeaWater 02", "Sediment 01", "Sediment 02", "Sediment 03")

res.pca <- dudi.pca(dat2,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5            # Nombre d'axes gardés
)

str(dat)

#dat$index  <- factor(dat$index, levels = c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03"))
#unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "all", 
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3", 
                            "Sed_1", "Sed_2", "Sed_3",
                            "SW_1", "SW_2"),                pointsize = 2,
                #xlim=c(-5, 5) ,
                #ylim=c(-5, 5) ,
                repel=T,
                title = "PCA - Biplot for phyla data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod", 
                            "burlywood4", "antiquewhite4", "bisque4",
                            "cadetblue3", "cadetblue3"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

fviz_pca_biplot(res.pca, 
                label = "ind", 
                invisible ="var", 
                # geom=c("text", "point"),
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3", 
                            "Sed_1", "Sed_2", "Sed_3",
                            "SW_1", "SW_2"),                pointsize = 2,
                repel=T,
                #xlim=c(-5, 5) ,
                #ylim=c(-5, 5) ,
                title = "PCA for phyla data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod", 
                            "burlywood4", "antiquewhite4", "bisque4",
                            "cadetblue3", "cadetblue3"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

options(ggrepel.max.overlaps = Inf)



fviz_pca_biplot(res.pca, 
                label = "var", 
                invisible ="ind", 
                #geom.var = c("point", "text"),
                alpha.var=0.3,
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3", 
                            "Sed_1", "Sed_2", "Sed_3",
                            "SW_1", "SW_2"),                pointsize = 2,
                #xlim=c(-4, 4) ,
                #ylim=c(-4, 4) ,
                title = "Contributive arrows for phyla data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod", 
                            "burlywood4", "antiquewhite4", "bisque4",
                            "cadetblue3", "cadetblue3"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

















############################################################################################################
######### db-RDA Par Phylum
############################################################################################################

# If data set without seawater
#taxa_bacteria <- dat[, -c(1:25) ]

taxa_bacteria <- dat[, -c(1:13) ]
#row.names(taxa_bacteria) <- c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03", "SeaWater 01", "SeaWater 02", "Sediment 01", "Sediment 02", "Sediment 03")

# If data set without seawater
#physico_chimic <- dat[, c(2:25) ]
physico_chimic <- dat[, c(2:13) ]
#row.names(physico_chimic) <- c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03", "SeaWater 01", "SeaWater 02", "Sediment 01", "Sediment 02", "Sediment 03")


head(physico_chimic)
head(taxa_bacteria)

######### db-RDA for bacteria
############################################################################################################

# Beta Diversity
# from https://www.rpubs.com/roalle/mres_2019

# Bray-Curtis Dissimilarity

# SUPPRESSION DES LIGNES OU IL N'Y A QUE DES ZEROS
taxa_bacteria_02 <- taxa_bacteria[-c(7:11),]
taxa_bacteria_02

# Calculate pairwise dissimilarity
taxa_bacteria.mdf <- as.matrix.data.frame(taxa_bacteria_02)
#rownames(taxa_bacteria.mdf) <- env.data$site

taxa_bacteria.bray <- vegdist(taxa_bacteria.mdf, method = "bray")
taxa_bacteria.bray

# PERMANOVA (permutational multivariate analysis of variance; Anderson 2001)
# Homogeneity of dispersion test
#permutest(betadisper(taxa_bacteria.bray, physico_chimic$ph_h2o))

# PERMANOVA analysis
# adonis(taxa_bacteria.bray ~ ph_h2o + ph_kcl + N_NO3_mg_kg + N_NH4_mg_kg + N_g_kg + C_g_kg + P_mg_kg + cations_change_Ca + cations_change_Mg + cations_change_Na + cations_change_K + change_capacity, data = physico_chimic, permutations = 999)

# Distance-Based Redundancy Analysis (dbRDA)
# from https://www.rpubs.com/roalle/mres_2019

# Distance-Based Redundancy Analysis (dbRDA; Legendre and Anderson, 1999) is an extension of Redundancy Analysis (RDA) which allows the use of 
# non-euclidean distance matrices as inputs (e.g. Bray-Curtis dissimilarity). The method works by first calculating a PCoA from the dissimilarity matrix, 
# and then subjecting the PCoA eigenvalues (which represent dissimilarities in euclidean space) to RDA. The method aims to detect linear relationships between 
# environmental variables and these dissimilarities. 
# dbRDA differs from PCoA in that it is a constrained analysis. 
# While PCoA axes are generated to explain maximum variation in the distance matrix, 
# dbRDA canonical axes are constructed as linear combinations of environmental variables. 
# Consequently these axes are constrained to the environmental variables in the model, and the ordination will be distinct from the PCoA. dbRDA allows us to 
# visualise how environmental variables contrain variation in community composition between our sites. 
# As dbRDA considers multiple environmental variables which are measured using different techniques, it is important to normalise environmental variables such 
# that they can be compared concurrently. This is achieved using a z-score transformation. Moreover, dbRDA is sensitive to multicolinearity 
# (i.e. high correlation between environmental variables). We can test whether any of our environmental variables are correlated by calculating variance 
# inflation factors (VIFs) and then removing selected terms from the model until all VIF scores are < 10.

# Z-Score Transformation
# Z-Score transform environmental data :
# The z-score transformation normalises environmental variables by their standard deviation from the mean 

physico_chimic

# IDEM : 
# SUPPRESSION DES LIGNES OU IL N'Y A QUE DES ZEROS
physico_chimic_02 <- physico_chimic[-c(7:11),]
physico_chimic_02
env.data.z <- physico_chimic_02

#env.data.z$Ca <- (env.data.z$Ca - mean(env.data.z$Ca))/sd(env.data.z$Ca)
#env.data.z$Mg <- (env.data.z$Mg - mean(env.data.z$Mg))/sd(env.data.z$Mg)
#env.data.z$Na <- (env.data.z$Na - mean(env.data.z$Na))/sd(env.data.z$Na)
#env.data.z$K <- (env.data.z$K - mean(env.data.z$K))/sd(env.data.z$K)
#env.data.z$P <- (env.data.z$P - mean(env.data.z$P))/sd(env.data.z$P)
#env.data.z$Co <- (env.data.z$Co - mean(env.data.z$Co))/sd(env.data.z$Co)
env.data.z$Cr <- (env.data.z$Cr - mean(env.data.z$Cr))/sd(env.data.z$Cr)
env.data.z$Fe <- (env.data.z$Fe - mean(env.data.z$Fe))/sd(env.data.z$Fe)
env.data.z$Mn <- (env.data.z$Mn - mean(env.data.z$Mn))/sd(env.data.z$Mn)
env.data.z$Ni <- (env.data.z$Ni - mean(env.data.z$Ni))/sd(env.data.z$Ni)
#env.data.z$Al <- (env.data.z$Al - mean(env.data.z$Al))/sd(env.data.z$Al)
#env.data.z$Cu <- (env.data.z$Cu - mean(env.data.z$Cu))/sd(env.data.z$Cu)
env.data.z$Zn <- (env.data.z$Zn - mean(env.data.z$Zn))/sd(env.data.z$Zn)
#env.data.z$Mo <- (env.data.z$Mo - mean(env.data.z$Mo))/sd(env.data.z$Mo)
#env.data.z$pH <- (env.data.z$pH - mean(env.data.z$pH))/sd(env.data.z$pH)
env.data.z$SiOH <- (env.data.z$SiOH - mean(env.data.z$SiOH))/sd(env.data.z$SiOH)
env.data.z$Nox <- (env.data.z$Nox - mean(env.data.z$Nox))/sd(env.data.z$Nox)
env.data.z$PO4 <- (env.data.z$PO4 - mean(env.data.z$PO4))/sd(env.data.z$PO4)
env.data.z$NH4 <- (env.data.z$NH4 - mean(env.data.z$NH4))/sd(env.data.z$NH4)
#env.data.z$DIC <- (env.data.z$DIC - mean(env.data.z$DIC))/sd(env.data.z$DIC)
env.data.z$POC <- (env.data.z$POC - mean(env.data.z$POC))/sd(env.data.z$POC)
env.data.z$DOC <- (env.data.z$DOC - mean(env.data.z$DOC))/sd(env.data.z$DOC)
env.data.z$PON <- (env.data.z$PON - mean(env.data.z$PON))/sd(env.data.z$PON)
#env.data.z$PN <- (env.data.z$PN - mean(env.data.z$PN))/sd(env.data.z$PN)
#env.data.z$Temp <- (env.data.z$Temp - mean(env.data.z$Temp))/sd(env.data.z$Temp)
#env.data.z$DIN <- (env.data.z$DIN - mean(env.data.z$DIN))/sd(env.data.z$DIN)


# Perform dbRDA
# construct full model and calculate VIF

dbRDA.full <- capscale(taxa_bacteria.bray ~ Mn +
                         Ni +
                         Cr +
                         SiOH +
                         Nox +
                         PO4 +
                         NH4 +
                         POC +
                         DOC +
                         PON +
                         Fe +
                         Zn,
                       env.data.z)
vif.cca(dbRDA.full)

dbRDA.full <- capscale(taxa_bacteria.mdf ~ Mn +
                         Ni +
                         Cr +
                         SiOH +
                         Nox +
                         PO4 +
                         NH4 +
                         POC +
                         DOC +
                         PON +
                         Fe +
                         Zn,
                       env.data.z)
vif.cca(dbRDA.full)

# VIF explaination, according http://jsnfri.fra.affrc.go.jp/gunshu/3dbrdaE.html
# Multiple-colinearlity between variables are examined by variance inflation factor (VIF). 
# When VIF > 10, problems of multiple-colinearlity become severe. So, remove variables of large VIF, one by one.

# construct reduced model and calculate VIF


dbRDA.mat <- capscale(taxa_bacteria.mdf ~ Mn +
                        Ni +
                        Cr +
                        SiOH +
                        Nox +
                        PO4 +
                        NH4 +
                        POC +
                        DOC +
                        PON +
                        Fe +
                        Zn,
                      env.data.z)


dbRDA.mat <- capscale(taxa_bacteria.mdf ~  Ni +
                        Zn,
                      env.data.z)



sort(vif.cca(dbRDA.mat))


# test overall significance of the analysis
anova(dbRDA.mat, permutations = 999)

# test significance of each environmental variable
anova(dbRDA.mat, by = "terms")

# summary of dbRDA model to extract total variance constrained and axis scores
summary(dbRDA.mat)

# Plot dbRDA

smry <- summary(dbRDA.mat)
scrs <- scores(dbRDA.mat)
df1  <- data.frame(smry$sites[,1:2]) # site scores for RDA1 and RDA2
df1$site <- rownames(df1)  #add site names
df2  <- data.frame(smry$biplot[,1:2])  # mapping environmental variables

site <- c("Dy", "Dy", "Dy", "RG", "RG", "RG")
ind <- c("Dy_1", "Dy_2", "Dy_3", "RG_1", "RG_2", "RG_3")
df_col <- data.frame(ind, site)
df_col$site <- factor(df_col$site, levels = c("Dy", "RG"))
head(df_col)



legend_title <- "Conditions"

ggplot(df1, aes(x=CAP1, y=CAP2, colour = df_col$site, shape=df_col$site)) + 
  geom_point(size = 3) +
  ggtitle("Distance-Based Redundancy Analysis") +
  scale_colour_manual(legend_title,labels=c("Dy", "RG"), values = c("coral3",  "darkgoldenrod4")) +
  scale_shape_manual(legend_title,labels=c("Dy", "RG"), values = c(16, 16)) +
  scale_fill_manual(legend_title,labels=c("Dy", "RG"), values = c("coral3",  "darkgoldenrod4")) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #stat_ellipse() +
  # #xlim(-50, 50) +
  # #ylim(-50, 50) +
  xlab("db-RDA1 (85.55%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("db-RDA2 (4.83%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  coord_fixed() +
  # ADD ggforce's ellipses
  ggforce::geom_mark_ellipse(aes(fill = df_col$site,
                                 color = df_col$site)) +
  # guides(fill=guide_legend("my awesome title")) +
  theme_bw() #+ theme(legend.position="none")









#install.packages("ggrepel")
library(ggrepel)

ggplot(df1, aes(x=CAP1, y=CAP2, colour = df_col$site)) + 
  geom_segment(data=df2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               # color="black", arrow=arrow(length=unit(0.01,"npc"))) +
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  ggrepel:: geom_text_repel(data=df2, aes(x=CAP1,y=CAP2,label=rownames(df2)), 
                            #    color="grey50", size=3) +
                            color="black", size=5) +
  #  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  # scale_colour_manual(values = c( "firebrick", "goldenrod", "olivedrab1", "green", "forestgreen")) +
  # geom_text(aes(label=rownames(df1),
  #                hjust=0,vjust=1.5), colour = "black",size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #  stat_ellipse() +
  #xlim(-2, 1.1) +
  # #ylim(-5, 6) +
  xlab("RDA1 (85.55%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("RDA2 (4.83%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  #coord_fixed() +
  theme_bw()







############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
######### QUE LES EPONGES
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################



#pc <- read.table("physico_chemical.txt", header=T)
pc <- read.table("Godd_data_ICP_for_all.txt", header=T)
pc <- read.table("New_Data_Decembre_2021.txt", header=T)
head(pc)
str(pc)
pc

#pc$index <- gsub("SeaWater_01", "SeaWater_downstream_0_2micro", pc$index)
#pc$index <- gsub("SeaWater_02", "SeaWater_upstream_0_2micro", pc$index)

bpa <- read.csv("Figure_01_output_table_of_bacterial_phyla_abundance.csv", sep=",")
head(bpa)
bpa$index <- gsub("Mycale", "Dy", bpa$index)
bpa$index <- gsub("Rhabdastrella_globostellata", "RG", bpa$index)
bpa$index <- gsub("SeaWater_downstream_0_2micro", "SW_1", bpa$index)
bpa$index <- gsub("SeaWater_upstream_0_2micro", "SW_3", bpa$index)
bpa$index <- gsub("Sediment", "Sed", bpa$index)
head(bpa)
bpa
# If physico_chemical.txt
#bpa_02 <- bpa[ c(1, 2, 3, 5,8,11), ]
#bpa_02 <- bpa[ c(1, 2, 3, 5,8,11, 13, 14), ]

bpa_02 <- bpa

head(bpa_02)

#a <- gsub("_1", "_01", pc$index)
#b <- gsub("_2", "_02", a)
#c <- gsub("_3", "_03", b)
#pc$index <- c

#a <- gsub("ata_1", "ata_01", bpa_02$index)
#b <- gsub("ata_2", "ata_02", a)
#c <- gsub("ata_3", "ata_03", b)
#d <- gsub("ale_1", "ale_01", c)
#e <- gsub("ale_2", "ale_02", d)
#f <- gsub("ale_3", "ale_03", e)

#bpa_02$index <- NULL
bpa_02$X <- NULL
bpa_02$X.NA. <- NULL
#bpa_02$index <- f

head(bpa_02)
head(pc)

bpa_02
pc

dat <- merge(pc, bpa_02, by.x="index", by.y="index")
head(dat)

str(dat)

# Too much integer...
dat[] <- lapply(dat, function(x) {
  if(is.integer(x)) as.numeric(as.character(x)) else x
})

sapply(dat, class)
str(dat) # ok now..



row.names(dat) <- c("Dy_1", "Dy_2", "Dy_3", 
                    "RG_1", "RG_2", "RG_3", 
                    "Sed_1", "Sed_2", "Sed_3",
                    "SW_1", "SW_3")



dat

dat <- dat[c(1:6) ,]

############################################################################################################
### PCA en fonction PHYSICO-CHIMIC
############################################################################################################

head(dat)

#dat2 <- dat[, -1]
#dat2$Cyanobacteria <- NULL
dat2 <- dat[, c(2:6) ]
#dat2 <- dat[, c(2:20) ]
#dat2 <- dat[, c(11:20) ]

head(dat2)


#dat2 <- dat2[ -c(7:8), ]

res.pca <- dudi.pca(dat2,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5            # Nombre d'axes gardés
)

str(dat)

#dat$index  <- factor(dat$index, levels = c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03"))
#unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "all", 
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3"),
                pointsize = 2,
                repel=T,
                #repel=TRUE,
                #xlim=c(-10, 10) ,
                #ylim=c(-10, 10) ,
                title = "PCA - Biplot for physico-chemical data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")


summary(res.pca)

dat$index  <- factor(dat$index, levels = c("Dy_1", "Dy_2", "Dy_3", 
                                           "RG_1", "RG_2", "RG_3"))
unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "ind", 
                invisible ="var",
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3"),
                pointsize = 2,
                #xlim=c(-10, 10) ,
                #ylim=c(-10, 10) ,
                title = "PCA for physico-chemical data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")


options(ggrepel.max.overlaps = Inf)


fviz_pca_biplot(res.pca, 
                label = "var", 
                invisible ="ind", 
                #geom.var = c("point", "text"),
                alpha.var=0.3,
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3"),
                pointsize = 2,
                #xlim=c(-4, 4) ,
                #ylim=c(-4, 4) ,
                title = "Contributive arrows for physico-chemical data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")



############################################################################################################
### PCA en fonction TAXO
############################################################################################################

head(dat)

#dat2 <- dat[, -1]
#dat2 <- dat[, -c(1:25) ]
dat2 <- dat[, -c(1:6) ]
dat2$Cyanobacteria <- NULL

head(dat2)


res.pca <- dudi.pca(dat2,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5            # Nombre d'axes gardés
)

str(dat)

#dat$index  <- factor(dat$index, levels = c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03"))
#unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "all", 
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3"),
                pointsize = 2,
                #xlim=c(-5, 5) ,
                #ylim=c(-5, 5) ,
                title = "PCA - Biplot for phyla data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

fviz_pca_biplot(res.pca, 
                label = "ind", 
                invisible ="var", 
                # geom=c("text", "point"),
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3"),
                pointsize = 2,
                #xlim=c(-5, 5) ,
                #ylim=c(-5, 5) ,
                title = "PCA for phyla data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

options(ggrepel.max.overlaps = Inf)



fviz_pca_biplot(res.pca, 
                label = "var", 
                invisible ="ind", 
                #geom.var = c("point", "text"),
                alpha.var=0.3,
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3"),
                pointsize = 2,
                #xlim=c(-4, 4) ,
                #ylim=c(-4, 4) ,
                title = "Contributive arrows for phyla data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

















############################################################################################################
######### db-RDA Par Phylum
############################################################################################################

# If data set without seawater
#taxa_bacteria <- dat[, -c(1:25) ]

taxa_bacteria <- dat[, -c(1:6) ]


head(dat)

# If data set without seawater
#physico_chimic <- dat[, c(2:25) ]
physico_chimic <- dat[, c(2:6) ]
#physico_chimic <- dat[, c(11:20) ]



head(physico_chimic)
head(taxa_bacteria)

######### db-RDA for bacteria
############################################################################################################

# Beta Diversity
# from https://www.rpubs.com/roalle/mres_2019

# Bray-Curtis Dissimilarity

# Calculate pairwise dissimilarity
taxa_bacteria.mdf <- as.matrix.data.frame(taxa_bacteria)
#rownames(taxa_bacteria.mdf) <- env.data$site

taxa_bacteria.bray <- vegdist(taxa_bacteria.mdf, method = "bray")
taxa_bacteria.bray

# PERMANOVA (permutational multivariate analysis of variance; Anderson 2001)
# Homogeneity of dispersion test
#permutest(betadisper(taxa_bacteria.bray, physico_chimic$ph_h2o))

# PERMANOVA analysis
# adonis(taxa_bacteria.bray ~ ph_h2o + ph_kcl + N_NO3_mg_kg + N_NH4_mg_kg + N_g_kg + C_g_kg + P_mg_kg + cations_change_Ca + cations_change_Mg + cations_change_Na + cations_change_K + change_capacity, data = physico_chimic, permutations = 999)

# Distance-Based Redundancy Analysis (dbRDA)
# from https://www.rpubs.com/roalle/mres_2019

# Distance-Based Redundancy Analysis (dbRDA; Legendre and Anderson, 1999) is an extension of Redundancy Analysis (RDA) which allows the use of 
# non-euclidean distance matrices as inputs (e.g. Bray-Curtis dissimilarity). The method works by first calculating a PCoA from the dissimilarity matrix, 
# and then subjecting the PCoA eigenvalues (which represent dissimilarities in euclidean space) to RDA. The method aims to detect linear relationships between 
# environmental variables and these dissimilarities. 
# dbRDA differs from PCoA in that it is a constrained analysis. 
# While PCoA axes are generated to explain maximum variation in the distance matrix, 
# dbRDA canonical axes are constructed as linear combinations of environmental variables. 
# Consequently these axes are constrained to the environmental variables in the model, and the ordination will be distinct from the PCoA. dbRDA allows us to 
# visualise how environmental variables contrain variation in community composition between our sites. 
# As dbRDA considers multiple environmental variables which are measured using different techniques, it is important to normalise environmental variables such 
# that they can be compared concurrently. This is achieved using a z-score transformation. Moreover, dbRDA is sensitive to multicolinearity 
# (i.e. high correlation between environmental variables). We can test whether any of our environmental variables are correlated by calculating variance 
# inflation factors (VIFs) and then removing selected terms from the model until all VIF scores are < 10.

# Z-Score Transformation
# Z-Score transform environmental data :
# The z-score transformation normalises environmental variables by their standard deviation from the mean 
env.data.z <- physico_chimic

#env.data.z$Ca <- (env.data.z$Ca - mean(env.data.z$Ca))/sd(env.data.z$Ca)
#env.data.z$Mg <- (env.data.z$Mg - mean(env.data.z$Mg))/sd(env.data.z$Mg)
#env.data.z$Na <- (env.data.z$Na - mean(env.data.z$Na))/sd(env.data.z$Na)
#env.data.z$K <- (env.data.z$K - mean(env.data.z$K))/sd(env.data.z$K)
#env.data.z$P <- (env.data.z$P - mean(env.data.z$P))/sd(env.data.z$P)
#env.data.z$Co <- (env.data.z$Co - mean(env.data.z$Co))/sd(env.data.z$Co)
env.data.z$Cr <- (env.data.z$Cr - mean(env.data.z$Cr))/sd(env.data.z$Cr)
env.data.z$Fe <- (env.data.z$Fe - mean(env.data.z$Fe))/sd(env.data.z$Fe)
env.data.z$Mn <- (env.data.z$Mn - mean(env.data.z$Mn))/sd(env.data.z$Mn)
env.data.z$Ni <- (env.data.z$Ni - mean(env.data.z$Ni))/sd(env.data.z$Ni)
#env.data.z$Al <- (env.data.z$Al - mean(env.data.z$Al))/sd(env.data.z$Al)
#env.data.z$Cu <- (env.data.z$Cu - mean(env.data.z$Cu))/sd(env.data.z$Cu)
env.data.z$Zn <- (env.data.z$Zn - mean(env.data.z$Zn))/sd(env.data.z$Zn)
#env.data.z$Mo <- (env.data.z$Mo - mean(env.data.z$Mo))/sd(env.data.z$Mo)
#env.data.z$pH <- (env.data.z$pH - mean(env.data.z$pH))/sd(env.data.z$pH)
#env.data.z$SiOH <- (env.data.z$SiOH - mean(env.data.z$SiOH))/sd(env.data.z$SiOH)
#env.data.z$Nox <- (env.data.z$Nox - mean(env.data.z$Nox))/sd(env.data.z$Nox)
#env.data.z$PO4 <- (env.data.z$PO4 - mean(env.data.z$PO4))/sd(env.data.z$PO4)
#env.data.z$NH4 <- (env.data.z$NH4 - mean(env.data.z$NH4))/sd(env.data.z$NH4)
#env.data.z$DIC <- (env.data.z$DIC - mean(env.data.z$DIC))/sd(env.data.z$DIC)
#env.data.z$POC <- (env.data.z$POC - mean(env.data.z$POC))/sd(env.data.z$POC)
#env.data.z$PN <- (env.data.z$PN - mean(env.data.z$PN))/sd(env.data.z$PN)
#env.data.z$Temp <- (env.data.z$Temp - mean(env.data.z$Temp))/sd(env.data.z$Temp)
#env.data.z$DIN <- (env.data.z$DIN - mean(env.data.z$DIN))/sd(env.data.z$DIN)

# Perform dbRDA
# construct full model and calculate VIF

dbRDA.full <- capscale(taxa_bacteria.bray ~ Mn +
                         Ni +
                         Cr +
                         Zn ,
                       env.data.z)
vif.cca(dbRDA.full)

dbRDA.mat <- capscale(taxa_bacteria.bray ~ Mn +
                        Ni +
                        Cr +
                        Zn ,
                      env.data.z)

# VIF explaination, according http://jsnfri.fra.affrc.go.jp/gunshu/3dbrdaE.html
# Multiple-colinearlity between variables are examined by variance inflation factor (VIF). 
# When VIF > 10, problems of multiple-colinearlity become severe. So, remove variables of large VIF, one by one.

# construct reduced model and calculate VIF
dbRDA.mat <- capscale(taxa_bacteria.mdf ~ Mn +
                        Ni +
                        Zn ,
                      env.data.z)


sort(vif.cca(dbRDA.mat))


# test overall significance of the analysis
anova(dbRDA.mat, permutations = 999)

# test significance of each environmental variable
anova(dbRDA.mat, by = "terms")

# summary of dbRDA model to extract total variance constrained and axis scores
summary(dbRDA.mat)

# Plot dbRDA

smry <- summary(dbRDA.mat)
scrs <- scores(dbRDA.mat)
df1  <- data.frame(smry$sites[,1:2]) # site scores for RDA1 and RDA2
df1$site <- rownames(df1)  #add site names
df2  <- data.frame(smry$biplot[,1:2])  # mapping environmental variables

site <- c("Dy", "Dy", "Dy", "RG", "RG", "RG")
ind <- c("Dy_1", "Dy_2", "Dy_3", "RG_1", "RG_2", "RG_3")
df_col <- data.frame(ind, site)
df_col$site <- factor(df_col$site, levels = c("Dy", "Dy", "Dy", "RG", "RG", "RG"))
head(df_col)



legend_title <- "Conditions"

ggplot(df1, aes(x=CAP1, y=CAP2, colour = df_col$site, shape=df_col$site)) + 
  geom_point(size = 3) +
  ggtitle("Distance-Based Redundancy Analysis") +
  scale_colour_manual(legend_title,labels=c("My", "RG"), values = c("coral3",  "darkgoldenrod4")) +
  scale_shape_manual(legend_title,labels=c("My", "RG"), values = c(16, 16, 16, 16)) +
  scale_fill_manual(legend_title,labels=c("My", "RG"), values = c("coral3",  "darkgoldenrod4")) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #stat_ellipse() +
  # #xlim(-50, 50) +
  # #ylim(-50, 50) +
  xlab("db-RDA1 (85.79%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("db-RDA2 (7.86%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  coord_fixed() +
  # ADD ggforce's ellipses
  ggforce::geom_mark_ellipse(aes(fill = df_col$site,
                                 color = df_col$site)) +
  # guides(fill=guide_legend("my awesome title")) +
  theme_bw() #+ theme(legend.position="none")









#install.packages("ggrepel")
library(ggrepel)

ggplot(df1, aes(x=CAP1, y=CAP2, colour = df_col$site)) + 
  geom_segment(data=df2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               # color="black", arrow=arrow(length=unit(0.01,"npc"))) +
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  ggrepel:: geom_text_repel(data=df2, aes(x=CAP1,y=CAP2,label=rownames(df2)), 
                            #    color="grey50", size=3) +
                            color="black", size=5) +
  #  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  # scale_colour_manual(values = c( "firebrick", "goldenrod", "olivedrab1", "green", "forestgreen")) +
  # geom_text(aes(label=rownames(df1),
  #                hjust=0,vjust=1.5), colour = "black",size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #  stat_ellipse() +
  #xlim(-2, 1.1) +
  # #ylim(-5, 6) +
  xlab("db-RDA1 (85.79%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("db-RDA2 (7.86%)")  + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  #coord_fixed() +
  theme_bw()







ggplot(df1, aes(x=CAP1, y=CAP2, colour = df_col$site)) + #, shape=df_col$site)) + 
  geom_point(size = 3) +
  geom_segment(data=df2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               # color="black", arrow=arrow(length=unit(0.01,"npc"))) +
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  ggrepel:: geom_text_repel(data=df2, aes(x=CAP1,y=CAP2,label=rownames(df2)), 
                            #    color="grey50", size=3) +
                            color="black", size=5) +
  ggtitle("Distance-Based Redundancy Analysis") +
  scale_colour_manual(legend_title,labels=c("My", "RG"), values = c("coral3",  "darkgoldenrod4")) +
  scale_shape_manual(legend_title,labels=c("My", "RG"), values = c(16, 16, 16, 16)) +
  scale_fill_manual(legend_title,labels=c("My", "RG"), values = c("coral3",  "darkgoldenrod4")) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #stat_ellipse() +
  # #xlim(-50, 50) +
  # #ylim(-50, 50) +
  xlab("db-RDA1 (85.79%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("db-RDA2 (7.86%)")  + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  coord_fixed() +
  # ADD ggforce's ellipses
  ggforce::geom_mark_ellipse(aes(fill = df_col$site,
                                 color = df_col$site)) +
  # guides(fill=guide_legend("my awesome title")) +
  theme_bw() #+ theme(legend.position="none")










############################################################################################################
############################################################################################################
######### QUE CHLOROFLEXI
############################################################################################################
############################################################################################################

library(ggplot2)
library(ade4)
library(factoextra)
library(magrittr)
library(tidyverse)
library(vegan)
library(ggordiplots)
library(multcomp)
library(ggpubr)
library(gridExtra)
library(grid)

#pc <- read.table("physico_chemical.txt", header=T)
#pc <- read.table("physico_chemical_02.txt", header=T)
#pc <- read.table("Godd_data_ICP_for_all.txt", header=T)

head(pc)


#pc$index <- gsub("SeaWater_01", "SeaWater_downstream_0_2micro", pc$index)
#pc$index <- gsub("SeaWater_02", "SeaWater_upstream_0_2micro", pc$index)



bpa <- read.csv("Figure_02_output_table_chloroflexi_order_abundance.csv", sep=";")
head(bpa)
bpa$index <- gsub("Mycale", "Dy", bpa$index)
bpa$index <- gsub("Rhabdastrella_globostellata", "RG", bpa$index)
bpa$index <- gsub("SeaWater_downstream_0_2micro", "SW_1", bpa$index)
bpa$index <- gsub("SeaWater_upstream_0_2micro", "SW_3", bpa$index)
bpa$index <- gsub("Sediment", "Sed", bpa$index)
head(bpa)
bpa
# If physico_chemical.txt
#bpa_02 <- bpa[ c(1, 2, 3, 5,8,11), ]
bpa_02 <- bpa

head(bpa_02)

#a <- gsub("ata_1", "ata_01", bpa_02$index)
#b <- gsub("ata_2", "ata_02", a)
#c <- gsub("ata_3", "ata_03", b)
#d <- gsub("ale_1", "ale_01", c)
#e <- gsub("ale_2", "ale_02", d)
#f <- gsub("ale_3", "ale_03", e)

#bpa_02$index <- NULL
bpa_02$X <- NULL
bpa_02$X.NA. <- NULL
#bpa_02$index <- f

head(bpa_02)
head(pc)

dat <- merge(pc, bpa_02, by.x="index", by.y="index")
head(dat)
dat

str(dat)

# Too much integer...
dat[] <- lapply(dat, function(x) {
  if(is.integer(x)) as.numeric(as.character(x)) else x
})

sapply(dat, class)
str(dat) # ok now..



row.names(dat) <- c("Dy_1", "Dy_2", "Dy_3", 
                    "RG_1", "RG_2", "RG_3", 
                    "SW_1", "SW_2",
                    "Sed_1", "Sed_2", "Sed_3")

dat <- dat[c(1:6) ,]

############################################################################################################
### PCA en fonction PHYSICO-CHIMIC
############################################################################################################

#dat2 <- dat[, -1]
#dat2 <- dat[, c(2:25) ]
dat2 <- dat[, c(2:6) ]


res.pca <- dudi.pca(dat2,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5            # Nombre d'axes gardés
)

str(dat)

#dat$index  <- factor(dat$index, levels = c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03"))
#unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "all", 
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3"),                pointsize = 2,
                repel=T,
                #repel=TRUE,
                #xlim=c(-10, 10) ,
                #ylim=c(-10, 10) ,
                title = "PCA - Biplot for physico-chemical data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

summary(res.pca)

dat$index  <- factor(dat$index, levels = c("Dy_1", "Dy_2", "Dy_3", 
                                           "RG_1", "RG_2", "RG_3"))
unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "ind", 
                invisible ="var",
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3"),                pointsize = 2,
                #xlim=c(-10, 10) ,
                #ylim=c(-10, 10) ,
                title = "PCA for physico-chemical data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

options(ggrepel.max.overlaps = Inf)



fviz_pca_biplot(res.pca, 
                label = "var", 
                invisible ="ind", 
                #geom.var = c("point", "text"),
                alpha.var=0.3,
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3"),                pointsize = 2,
                #xlim=c(-4, 4) ,
                #ylim=c(-4, 4) ,
                title = "Contributive arrows for physico-chemical data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")



############################################################################################################
### PCA en fonction TAXO
############################################################################################################

#dat2 <- dat[, -1]
#dat2 <- dat[, -c(1:25) ]
dat2 <- dat[, -c(1:6) ]
dat2$Cyanobacteria <- NULL

#row.names(dat2) <- c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03", "SeaWater 01", "SeaWater 02", "Sediment 01", "Sediment 02", "Sediment 03")

res.pca <- dudi.pca(dat2,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5            # Nombre d'axes gardés
)

str(dat)

#dat$index  <- factor(dat$index, levels = c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03"))
#unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "all", 
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3"),                pointsize = 2,
                #xlim=c(-5, 5) ,
                #ylim=c(-5, 5) ,
                repel=T,
                title = "PCA - Biplot for phyla data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

fviz_pca_biplot(res.pca, 
                label = "ind", 
                invisible ="var", 
                # geom=c("text", "point"),
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3"),                pointsize = 2,
                repel=T,
                #xlim=c(-5, 5) ,
                #ylim=c(-5, 5) ,
                title = "PCA for phyla data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

options(ggrepel.max.overlaps = Inf)



fviz_pca_biplot(res.pca, 
                label = "var", 
                invisible ="ind", 
                #geom.var = c("point", "text"),
                alpha.var=0.3,
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Dy_1", "Dy_2", "Dy_3", 
                            "RG_1", "RG_2", "RG_3"),                pointsize = 2,
                #xlim=c(-4, 4) ,
                #ylim=c(-4, 4) ,
                title = "Contributive arrows for phyla data",
                palette = c("coral3", "chocolate3", "chocolate1",  
                            "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

















############################################################################################################
######### db-RDA Par Phylum
############################################################################################################

# If data set without seawater
#taxa_bacteria <- dat[, -c(1:25) ]

taxa_bacteria <- dat[, -c(1:6) ]
#row.names(taxa_bacteria) <- c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03", "SeaWater 01", "SeaWater 02", "Sediment 01", "Sediment 02", "Sediment 03")

# If data set without seawater
#physico_chimic <- dat[, c(2:25) ]
physico_chimic <- dat[, c(2:6) ]
#row.names(physico_chimic) <- c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03", "SeaWater 01", "SeaWater 02", "Sediment 01", "Sediment 02", "Sediment 03")


head(physico_chimic)
head(taxa_bacteria)

######### db-RDA for bacteria
############################################################################################################

# Beta Diversity
# from https://www.rpubs.com/roalle/mres_2019

# Bray-Curtis Dissimilarity

# SUPPRESSION DES LIGNES OU IL N'Y A QUE DES ZEROS
#taxa_bacteria_02 <- taxa_bacteria[-c(7:11),]
#taxa_bacteria_02

# Calculate pairwise dissimilarity
taxa_bacteria.mdf <- as.matrix.data.frame(taxa_bacteria_02)
#rownames(taxa_bacteria.mdf) <- env.data$site

taxa_bacteria.bray <- vegdist(taxa_bacteria.mdf, method = "bray")
taxa_bacteria.bray

# PERMANOVA (permutational multivariate analysis of variance; Anderson 2001)
# Homogeneity of dispersion test
#permutest(betadisper(taxa_bacteria.bray, physico_chimic$ph_h2o))

# PERMANOVA analysis
# adonis(taxa_bacteria.bray ~ ph_h2o + ph_kcl + N_NO3_mg_kg + N_NH4_mg_kg + N_g_kg + C_g_kg + P_mg_kg + cations_change_Ca + cations_change_Mg + cations_change_Na + cations_change_K + change_capacity, data = physico_chimic, permutations = 999)

# Distance-Based Redundancy Analysis (dbRDA)
# from https://www.rpubs.com/roalle/mres_2019

# Distance-Based Redundancy Analysis (dbRDA; Legendre and Anderson, 1999) is an extension of Redundancy Analysis (RDA) which allows the use of 
# non-euclidean distance matrices as inputs (e.g. Bray-Curtis dissimilarity). The method works by first calculating a PCoA from the dissimilarity matrix, 
# and then subjecting the PCoA eigenvalues (which represent dissimilarities in euclidean space) to RDA. The method aims to detect linear relationships between 
# environmental variables and these dissimilarities. 
# dbRDA differs from PCoA in that it is a constrained analysis. 
# While PCoA axes are generated to explain maximum variation in the distance matrix, 
# dbRDA canonical axes are constructed as linear combinations of environmental variables. 
# Consequently these axes are constrained to the environmental variables in the model, and the ordination will be distinct from the PCoA. dbRDA allows us to 
# visualise how environmental variables contrain variation in community composition between our sites. 
# As dbRDA considers multiple environmental variables which are measured using different techniques, it is important to normalise environmental variables such 
# that they can be compared concurrently. This is achieved using a z-score transformation. Moreover, dbRDA is sensitive to multicolinearity 
# (i.e. high correlation between environmental variables). We can test whether any of our environmental variables are correlated by calculating variance 
# inflation factors (VIFs) and then removing selected terms from the model until all VIF scores are < 10.

# Z-Score Transformation
# Z-Score transform environmental data :
# The z-score transformation normalises environmental variables by their standard deviation from the mean 

physico_chimic

# IDEM : 
# SUPPRESSION DES LIGNES OU IL N'Y A QUE DES ZEROS
#physico_chimic_02 <- physico_chimic[-c(7:11),]
#physico_chimic_02
env.data.z <- physico_chimic

#env.data.z$Ca <- (env.data.z$Ca - mean(env.data.z$Ca))/sd(env.data.z$Ca)
#env.data.z$Mg <- (env.data.z$Mg - mean(env.data.z$Mg))/sd(env.data.z$Mg)
#env.data.z$Na <- (env.data.z$Na - mean(env.data.z$Na))/sd(env.data.z$Na)
#env.data.z$K <- (env.data.z$K - mean(env.data.z$K))/sd(env.data.z$K)
#env.data.z$P <- (env.data.z$P - mean(env.data.z$P))/sd(env.data.z$P)
#env.data.z$Co <- (env.data.z$Co - mean(env.data.z$Co))/sd(env.data.z$Co)
env.data.z$Cr <- (env.data.z$Cr - mean(env.data.z$Cr))/sd(env.data.z$Cr)
env.data.z$Fe <- (env.data.z$Fe - mean(env.data.z$Fe))/sd(env.data.z$Fe)
env.data.z$Mn <- (env.data.z$Mn - mean(env.data.z$Mn))/sd(env.data.z$Mn)
env.data.z$Ni <- (env.data.z$Ni - mean(env.data.z$Ni))/sd(env.data.z$Ni)
#env.data.z$Al <- (env.data.z$Al - mean(env.data.z$Al))/sd(env.data.z$Al)
#env.data.z$Cu <- (env.data.z$Cu - mean(env.data.z$Cu))/sd(env.data.z$Cu)
env.data.z$Zn <- (env.data.z$Zn - mean(env.data.z$Zn))/sd(env.data.z$Zn)
#env.data.z$Mo <- (env.data.z$Mo - mean(env.data.z$Mo))/sd(env.data.z$Mo)
#env.data.z$pH <- (env.data.z$pH - mean(env.data.z$pH))/sd(env.data.z$pH)
#env.data.z$SiOH <- (env.data.z$SiOH - mean(env.data.z$SiOH))/sd(env.data.z$SiOH)
#env.data.z$Nox <- (env.data.z$Nox - mean(env.data.z$Nox))/sd(env.data.z$Nox)
#env.data.z$PO4 <- (env.data.z$PO4 - mean(env.data.z$PO4))/sd(env.data.z$PO4)
#env.data.z$NH4 <- (env.data.z$NH4 - mean(env.data.z$NH4))/sd(env.data.z$NH4)
#env.data.z$DIC <- (env.data.z$DIC - mean(env.data.z$DIC))/sd(env.data.z$DIC)
#env.data.z$POC <- (env.data.z$POC - mean(env.data.z$POC))/sd(env.data.z$POC)
#env.data.z$PN <- (env.data.z$PN - mean(env.data.z$PN))/sd(env.data.z$PN)
#env.data.z$Temp <- (env.data.z$Temp - mean(env.data.z$Temp))/sd(env.data.z$Temp)
#env.data.z$DIN <- (env.data.z$DIN - mean(env.data.z$DIN))/sd(env.data.z$DIN)

# Perform dbRDA
# construct full model and calculate VIF

dbRDA.full <- capscale(taxa_bacteria.bray ~ Mn +
                         Ni +
                         Cr +
                         Zn,
                       env.data.z)
vif.cca(dbRDA.full)

# VIF explaination, according http://jsnfri.fra.affrc.go.jp/gunshu/3dbrdaE.html
# Multiple-colinearlity between variables are examined by variance inflation factor (VIF). 
# When VIF > 10, problems of multiple-colinearlity become severe. So, remove variables of large VIF, one by one.

# construct reduced model and calculate VIF


dbRDA.mat <- capscale(taxa_bacteria.mdf ~ Mn +
                        Ni +
                        Cr +
                        Fe +
                        Zn,
                      env.data.z)

dbRDA.mat <- capscale(taxa_bacteria.mdf ~ Cr +
                        Zn,
                      env.data.z)


sort(vif.cca(dbRDA.mat))


# test overall significance of the analysis
anova(dbRDA.mat, permutations = 999)

# test significance of each environmental variable
anova(dbRDA.mat, by = "terms")

# summary of dbRDA model to extract total variance constrained and axis scores
summary(dbRDA.mat)

# Plot dbRDA

smry <- summary(dbRDA.mat)
scrs <- scores(dbRDA.mat)
df1  <- data.frame(smry$sites[,1:2]) # site scores for RDA1 and RDA2
df1$site <- rownames(df1)  #add site names
df2  <- data.frame(smry$biplot[,1:2])  # mapping environmental variables

site <- c("Dy", "Dy", "Dy", "RG","RG","RG")
ind <- c("Dy_1", "Dy_2", "Dy_3", "RG_1","RG_2","RG_3")
df_col <- data.frame(ind, site)
df_col$site <- factor(df_col$site, levels = c("Dy", "RG"))
head(df_col)



legend_title <- "Conditions"

ggplot(df1, aes(x=CAP1, y=CAP2, colour = df_col$site, shape=df_col$site)) + 
  geom_point(size = 3) +
  ggtitle("Distance-Based Redundancy Analysis") +
  scale_colour_manual(legend_title,labels=c("Dy", "RG"), values = c("coral3",  "darkgoldenrod4")) +
  scale_shape_manual(legend_title,labels=c("Dy", "RG"), values = c(16, 16)) +
  scale_fill_manual(legend_title,labels=c("Dy", "RG"), values = c("coral3",  "darkgoldenrod4")) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #stat_ellipse() +
  # #xlim(-50, 50) +
  # #ylim(-50, 50) +
  xlab("db-RDA1 (85.33%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("db-RDA2 (52.96%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  coord_fixed() +
  # ADD ggforce's ellipses
  ggforce::geom_mark_ellipse(aes(fill = df_col$site,
                                 color = df_col$site)) +
  # guides(fill=guide_legend("my awesome title")) +
  theme_bw() #+ theme(legend.position="none")









#install.packages("ggrepel")
library(ggrepel)

ggplot(df1, aes(x=CAP1, y=CAP2, colour = df_col$site)) + 
  geom_segment(data=df2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               # color="black", arrow=arrow(length=unit(0.01,"npc"))) +
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  ggrepel:: geom_text_repel(data=df2, aes(x=CAP1,y=CAP2,label=rownames(df2)), 
                            #    color="grey50", size=3) +
                            color="black", size=5) +
  #  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  # scale_colour_manual(values = c( "firebrick", "goldenrod", "olivedrab1", "green", "forestgreen")) +
  # geom_text(aes(label=rownames(df1),
  #                hjust=0,vjust=1.5), colour = "black",size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #  stat_ellipse() +
  #xlim(-2, 1.1) +
  # #ylim(-5, 6) +
  xlab("RDA1 (85.33%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("RDA2 (52.96%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  #coord_fixed() +
  theme_bw()







############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
######### QUE LES EPONGES + ENV 2019
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################




#pc <- read.table("physico_chemical.txt", header=T)
pc <- read.table("Env_data_Federica_sans_sed_sans_SW_only_2019.txt", header=T)
head(pc)
str(pc)



bpa <- read.csv("Figure_01_output_table_of_bacterial_phyla_abundance.csv", sep=",")
head(bpa)
# If physico_chemical.txt
#bpa_02 <- bpa[ c(1, 2, 3, 5,8,11), ]
#bpa_02 <- bpa[ c(1, 2, 3, 5,8,11, 13, 14), ]

bpa_02 <- bpa

head(bpa_02)

#a <- gsub("_1", "_01", pc$index)
#b <- gsub("_2", "_02", a)
#c <- gsub("_3", "_03", b)
#pc$index <- c

#a <- gsub("ata_1", "ata_01", bpa_02$index)
#b <- gsub("ata_2", "ata_02", a)
#c <- gsub("ata_3", "ata_03", b)
#d <- gsub("ale_1", "ale_01", c)
#e <- gsub("ale_2", "ale_02", d)
#f <- gsub("ale_3", "ale_03", e)



############################################################################################################
### PCA en fonction PHYSICO-CHIMIC
############################################################################################################

head(dat)

#dat2 <- dat[, -1]
#dat2$Cyanobacteria <- NULL
dat2 <- dat[, c(2:7) ]
#dat2 <- dat[, c(2:20) ]
#dat2 <- dat[, c(11:20) ]

#dat2 <- dat2[c(2:3) ,]

head(dat2)

row.names(dat2) <- c("Outer_1",  "Inner")

#dat2 <- as.data.frame(dat2)

#dat2 <- dat2[ -c(7:8), ]

res.pca <- dudi.pca(dat2,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5            # Nombre d'axes gardés
)

str(dat2)

#dat$index  <- factor(dat$index, levels = c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03"))
#unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "all", 
                col.ind = c("Outer",  "Inner"),
                pointsize = 2,
                repel=T,
                #repel=TRUE,
                #xlim=c(-10, 10) ,
                #ylim=c(-10, 10) ,
                title = "PCA - Biplot for physico-chemical data",
                palette = c("lightgoldenrod4",  "mistyrose4"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

fviz_pca_biplot(res.pca, 
                label = "all", 
                col.ind = c("Outer (1)",  "Inner (3)"),
                pointsize = 2,
                repel=T,
                #repel=TRUE,
                xlim=c(-10, 10) ,
                ylim=c(-10, 10) ,
                title = "PCA - Biplot for physico-chemical data",
                palette = c("lightgoldenrod4",  "mistyrose4"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

summary(res.pca)

dat$index  <- factor(dat$index, levels = c("Outer (1)",  "Inner (3)"))
unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "ind", 
                invisible ="var",
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Outer (1)",  "Inner (3)"),
                pointsize = 2,
                #  xlim=c(-10, 10) ,
                #  ylim=c(-10, 10) ,
                title = "PCA for physico-chemical data",
                palette = c("lightgoldenrod4",  "mistyrose4"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

fviz_pca_biplot(res.pca, 
                label = "ind", 
                invisible ="var",
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("My_1", "My_2", "My_3", 
                            "RG_1", "RG_2", "RG_3", 
                            "SW_1", "SW_2",
                            "Sed_1", "Sed_2", "Sed_3"),
                pointsize = 2,
                xlim=c(-10, 10) ,
                ylim=c(-10, 10) ,
                title = "PCA for physico-chemical data",
                palette = c("coral3", "chocolate3", "chocolate1",  "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod", "cadetblue3", "cadetblue3", "burlywood4", "antiquewhite4", "bisque4"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

options(ggrepel.max.overlaps = Inf)



fviz_pca_biplot(res.pca, 
                label = "var", 
                invisible ="ind", 
                #geom.var = c("point", "text"),
                alpha.var=0.3,
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("My_1", "My_2", "My_3", 
                            "RG_1", "RG_2", "RG_3", 
                            "SW_1", "SW_2",
                            "Sed_1", "Sed_2", "Sed_3"),
                pointsize = 2,
                xlim=c(-4, 4) ,
                ylim=c(-4, 4) ,
                title = "Contributive arrows for physico-chemical data",
                palette = c("coral3", "chocolate3", "chocolate1",  "darkgoldenrod4", "darkgoldenrod3", "darkgoldenrod", "cadetblue3", "cadetblue3", "burlywood4", "antiquewhite4", "bisque4"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")



fviz_pca_biplot(res.pca, 
                label = "var", 
                invisible ="ind", 
                #geom.var = c("point", "text"),
                alpha.var=0.3,
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Outer (1)",  "Inner (3)"),
                pointsize = 2,
                xlim=c(-4, 4) ,
                ylim=c(-4, 4) ,
                title = "Contributive arrows for physico-chemical data",
                palette = c("lightgoldenrod4",  "mistyrose4"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")



############################################################################################################
### PCA en fonction TAXO
############################################################################################################

head(dat)

#dat2 <- dat[, -1]
#dat2 <- dat[, -c(1:25) ]
dat2 <- dat[, -c(1:12) ]
dat2$Cyanobacteria <- NULL

head(dat2)


dat2  <- rbind(rbind(dat2["My_1",] + dat2["RG_1",], dat2["My_2",] + dat2["RG_2",]), dat2["My_3",] + dat2["RG_3",])

row.names(dat2) <- c("Outer (1)",  "Inner (3)")

res.pca <- dudi.pca(dat2,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5            # Nombre d'axes gardés
)

str(dat)

dat$index  <- factor(dat$index, levels = c("Outer (1)",  "Inner (3)"))
unique(dat$index)



fviz_pca_biplot(res.pca, 
                label = "all", 
                col.ind = c("Outer (1)",  "Inner (3)"),
                pointsize = 2,
                xlim=c(-5, 5) ,
                ylim=c(-5, 5) ,
                title = "PCA - Biplot for phyla data",
                palette = c("lightgoldenrod4",  "mistyrose4"),
                
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

fviz_pca_biplot(res.pca, 
                label = "ind", 
                invisible ="var", 
                # geom=c("text", "point"),
                col.ind = c("Outer (1)",  "Inner (3)"),
                pointsize = 2,
                xlim=c(-5, 5) ,
                ylim=c(-5, 5) ,
                title = "PCA for phyla data",
                palette = c("lightgoldenrod4",  "mistyrose4"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

options(ggrepel.max.overlaps = Inf)



fviz_pca_biplot(res.pca, 
                label = "var", 
                invisible ="ind", 
                #geom.var = c("point", "text"),
                alpha.var=0.3,
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Outer (1)",  "Inner (3)"),
                pointsize = 2,
                xlim=c(-4, 4) ,
                ylim=c(-4, 4) ,
                title = "Contributive arrows for phyla data",
                palette = c("lightgoldenrod4",  "mistyrose4"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")




############################################################################################################
######### db-RDA Par Phylum
############################################################################################################

# If data set without seawater
#taxa_bacteria <- dat[, -c(1:25) ]

taxa_bacteria <- dat[, -c(1:12) ]
taxa_bacteria  <- rbind(rbind(taxa_bacteria[1,] + taxa_bacteria[4,], taxa_bacteria[2,] + taxa_bacteria[5,]), taxa_bacteria[3,] + taxa_bacteria[6,])
row.names(taxa_bacteria) <- c("Outer (1)",  "Inner (3)")

head(dat)

# If data set without seawater
#physico_chimic <- dat[, c(2:25) ]
physico_chimic <- dat[, c(2:12) ]
#physico_chimic <- dat[, c(11:20) ]
physico_chimic  <- rbind(rbind(physico_chimic[1,] + physico_chimic[4,], physico_chimic[2,] + physico_chimic[5,]), physico_chimic[3,] + physico_chimic[6,])
row.names(physico_chimic) <- c("Outer (1)",  "Inner (3)")


head(physico_chimic)
head(taxa_bacteria)

######### db-RDA for bacteria
############################################################################################################

# Beta Diversity
# from https://www.rpubs.com/roalle/mres_2019

# Bray-Curtis Dissimilarity

# Calculate pairwise dissimilarity
taxa_bacteria.mdf <- as.matrix.data.frame(taxa_bacteria)
#rownames(taxa_bacteria.mdf) <- env.data$site

taxa_bacteria.bray <- vegdist(taxa_bacteria.mdf, method = "bray")
taxa_bacteria.bray

# PERMANOVA (permutational multivariate analysis of variance; Anderson 2001)
# Homogeneity of dispersion test
#permutest(betadisper(taxa_bacteria.bray, physico_chimic$ph_h2o))

# PERMANOVA analysis
# adonis(taxa_bacteria.bray ~ ph_h2o + ph_kcl + N_NO3_mg_kg + N_NH4_mg_kg + N_g_kg + C_g_kg + P_mg_kg + cations_change_Ca + cations_change_Mg + cations_change_Na + cations_change_K + change_capacity, data = physico_chimic, permutations = 999)

# Distance-Based Redundancy Analysis (dbRDA)
# from https://www.rpubs.com/roalle/mres_2019

# Distance-Based Redundancy Analysis (dbRDA; Legendre and Anderson, 1999) is an extension of Redundancy Analysis (RDA) which allows the use of 
# non-euclidean distance matrices as inputs (e.g. Bray-Curtis dissimilarity). The method works by first calculating a PCoA from the dissimilarity matrix, 
# and then subjecting the PCoA eigenvalues (which represent dissimilarities in euclidean space) to RDA. The method aims to detect linear relationships between 
# environmental variables and these dissimilarities. 
# dbRDA differs from PCoA in that it is a constrained analysis. 
# While PCoA axes are generated to explain maximum variation in the distance matrix, 
# dbRDA canonical axes are constructed as linear combinations of environmental variables. 
# Consequently these axes are constrained to the environmental variables in the model, and the ordination will be distinct from the PCoA. dbRDA allows us to 
# visualise how environmental variables contrain variation in community composition between our sites. 
# As dbRDA considers multiple environmental variables which are measured using different techniques, it is important to normalise environmental variables such 
# that they can be compared concurrently. This is achieved using a z-score transformation. Moreover, dbRDA is sensitive to multicolinearity 
# (i.e. high correlation between environmental variables). We can test whether any of our environmental variables are correlated by calculating variance 
# inflation factors (VIFs) and then removing selected terms from the model until all VIF scores are < 10.

# Z-Score Transformation
# Z-Score transform environmental data :
# The z-score transformation normalises environmental variables by their standard deviation from the mean 
env.data.z <- physico_chimic

#env.data.z$Mo <- (env.data.z$Mo - mean(env.data.z$Mo))/sd(env.data.z$Mo)
env.data.z$pH <- (env.data.z$pH - mean(env.data.z$pH))/sd(env.data.z$pH)
env.data.z$SiOH <- (env.data.z$SiOH - mean(env.data.z$SiOH))/sd(env.data.z$SiOH)
env.data.z$Nox <- (env.data.z$Nox - mean(env.data.z$Nox))/sd(env.data.z$Nox)
env.data.z$PO4 <- (env.data.z$PO4 - mean(env.data.z$PO4))/sd(env.data.z$PO4)
env.data.z$NH4 <- (env.data.z$NH4 - mean(env.data.z$NH4))/sd(env.data.z$NH4)
env.data.z$DIC <- (env.data.z$DIC - mean(env.data.z$DIC))/sd(env.data.z$DIC)
env.data.z$DOC <- (env.data.z$DOC - mean(env.data.z$DOC))/sd(env.data.z$DOC)
env.data.z$POC <- (env.data.z$POC - mean(env.data.z$POC))/sd(env.data.z$POC)
env.data.z$PN <- (env.data.z$PN - mean(env.data.z$PN))/sd(env.data.z$PN)
env.data.z$Temp <- (env.data.z$Temp - mean(env.data.z$Temp))/sd(env.data.z$Temp)
env.data.z$DIN <- (env.data.z$DIN - mean(env.data.z$DIN))/sd(env.data.z$DIN)

# Perform dbRDA
# construct full model and calculate VIF


dbRDA.full <- capscale(taxa_bacteria.bray ~ pH +
                         SiOH +
                         Nox +
                         PO4 +
                         NH4 +
                         DOC +
                         DIC +
                         POC +
                         PN +
                         Temp +
                         DIN ,
                       env.data.z)
vif.cca(dbRDA.full)
vif.cca(dbRDA.mat)

dbRDA.mat <- capscale(taxa_bacteria.bray ~ SiOH +
                        DOC  ,
                      env.data.z)

# VIF explaination, according http://jsnfri.fra.affrc.go.jp/gunshu/3dbrdaE.html
# Multiple-colinearlity between variables are examined by variance inflation factor (VIF). 
# When VIF > 10, problems of multiple-colinearlity become severe. So, remove variables of large VIF, one by one.

# construct reduced model and calculate VIF

dbRDA.mat <- capscale(taxa_bacteria.bray ~ SiOH +
                        DOC  ,
                      env.data.z)

sort(vif.cca(dbRDA.mat))


# test overall significance of the analysis
anova(dbRDA.mat, permutations = 999)

# test significance of each environmental variable
anova(dbRDA.mat, by = "terms")

# summary of dbRDA model to extract total variance constrained and axis scores
summary(dbRDA.mat)

# Plot dbRDA

smry <- summary(dbRDA.mat)
scrs <- scores(dbRDA.mat)
df1  <- data.frame(smry$sites[,1:2]) # site scores for RDA1 and RDA2
df1$site <- rownames(df1)  #add site names
df2  <- data.frame(smry$biplot[,1:2])  # mapping environmental variables

site <- c("Mycale", "Mycale", "Mycale", "Rhabdastrella globostellata", "Rhabdastrella globostellata", "Rhabdastrella globostellata", "SeaWater", "SeaWater",  "Sediment", "Sediment", "Sediment")
ind <- c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03", "SeaWater 01", "SeaWater 02", "Sediment 01", "Sediment 02", "Sediment 03")
df_col <- data.frame(ind, site)
df_col$site <- factor(df_col$site, levels = c("Mycale", "Rhabdastrella globostellata", "SeaWater", "Sediment"))
head(df_col)



legend_title <- "Conditions"

ggplot(df1, aes(x=CAP1, y=CAP2, colour = df_col$site, shape=df_col$site)) + 
  geom_point(size = 3) +
  ggtitle("Distance-Based Redundancy Analysis") +
  scale_colour_manual(legend_title,labels=c("Mycale", "Rhabdastrella globostellata", "SeaWater", "Sediment"), values = c("coral3",  "darkgoldenrod4", "cadetblue3", "bisque4")) +
  scale_shape_manual(legend_title,labels=c("Mycale", "Rhabdastrella globostellata", "SeaWater", "Sediment"), values = c(16, 16, 16, 16)) +
  scale_fill_manual(legend_title,labels=c("Mycale", "Rhabdastrella globostellata", "SeaWater", "Sediment"), values = c("coral3",  "darkgoldenrod4", "cadetblue3", "bisque4")) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #stat_ellipse() +
  # xlim(-50, 50) +
  # ylim(-50, 50) +
  xlab("db-RDA1 (79.01%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("db-RDA2 (16.58%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  coord_fixed() +
  # ADD ggforce's ellipses
  ggforce::geom_mark_ellipse(aes(fill = df_col$site,
                                 color = df_col$site)) +
  # guides(fill=guide_legend("my awesome title")) +
  theme_bw() #+ theme(legend.position="none")









#install.packages("ggrepel")
library(ggrepel)

ggplot(df1, aes(x=CAP1, y=CAP2, colour = df_col$site)) + 
  geom_segment(data=df2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               # color="black", arrow=arrow(length=unit(0.01,"npc"))) +
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  ggrepel:: geom_text_repel(data=df2, aes(x=CAP1,y=CAP2,label=rownames(df2)), 
                            #    color="grey50", size=3) +
                            color="black", size=5) +
  #  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  # scale_colour_manual(values = c( "firebrick", "goldenrod", "olivedrab1", "green", "forestgreen")) +
  # geom_text(aes(label=rownames(df1),
  #                hjust=0,vjust=1.5), colour = "black",size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #  stat_ellipse() +
  xlim(-2, 1.1) +
  # ylim(-5, 6) +
  xlab("RDA1 (79.01%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("RDA2 (16.58%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  #coord_fixed() +
  theme_bw()












############################################################################################################
############################################################################################################
######### QUE CHLOROFLEXI
############################################################################################################
############################################################################################################

library(ggplot2)
library(ade4)
library(factoextra)
library(magrittr)
library(tidyverse)
library(vegan)
library(ggordiplots)
library(multcomp)
library(ggpubr)
library(gridExtra)
library(grid)


pc <- read.table("Env_data_Federica_sans_sed_sans_SW.txt", header=T)

head(pc)


#pc$index <- gsub("SeaWater_01", "SeaWater_downstream_0_2micro", pc$index)
#pc$index <- gsub("SeaWater_02", "SeaWater_upstream_0_2micro", pc$index)



bpa <- read.csv("Figure_02_output_table_chloroflexi_order_abundance.csv", sep=";")
head(bpa)
# If physico_chemical.txt
#bpa_02 <- bpa[ c(1, 2, 3, 5,8,11), ]
bpa_02 <- bpa

head(bpa_02)

#a <- gsub("ata_1", "ata_01", bpa_02$index)
#b <- gsub("ata_2", "ata_02", a)
#c <- gsub("ata_3", "ata_03", b)
#d <- gsub("ale_1", "ale_01", c)
#e <- gsub("ale_2", "ale_02", d)
#f <- gsub("ale_3", "ale_03", e)

#bpa_02$index <- NULL
bpa_02$X <- NULL
bpa_02$X.NA. <- NULL
#bpa_02$index <- f

head(bpa_02)
head(pc)

dat <- merge(pc, bpa_02, by.x="index", by.y="index")
head(dat)

str(dat)

# Too much integer...
dat[] <- lapply(dat, function(x) {
  if(is.integer(x)) as.numeric(as.character(x)) else x
})

sapply(dat, class)
str(dat) # ok now..

# Merge
dat$index <- NULL
dat  <- rbind(rbind(dat[1,] + dat[4,], dat[2,] + dat[5,]), dat[3,] + dat[6,])
row.names(dat)<- c("Outer (1)",  "Inner (3)")

############################################################################################################
### PCA en fonction PHYSICO-CHIMIC
############################################################################################################

#dat2 <- dat[, -1]
#dat2 <- dat[, c(2:25) ]
dat2 <- dat[, c(2:12) ]


res.pca <- dudi.pca(dat2,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5            # Nombre d'axes gardés
)

str(dat)

#dat$index  <- factor(dat$index, levels = c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03"))
#unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "all", 
                col.ind = c("Outer (1)",  "Inner (3)"),                
                pointsize = 2,
                repel=T,
                #repel=TRUE,
                #  xlim=c(-10, 10) ,
                #  ylim=c(-10, 10) ,
                title = "PCA - Biplot for physico-chemical data",
                palette = c("lightgoldenrod4",  "mistyrose4"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

summary(res.pca)

dat$index  <- factor(dat$index, levels = c("Outer (1)",  "Inner (3)"))
unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "ind", 
                invisible ="var",
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Outer (1)",  "Inner (3)"),               
                pointsize = 2,
                #  xlim=c(-10, 10) ,
                #  ylim=c(-10, 10) ,
                title = "PCA for physico-chemical data",
                palette = c("lightgoldenrod4",  "mistyrose4"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

options(ggrepel.max.overlaps = Inf)



fviz_pca_biplot(res.pca, 
                label = "var", 
                invisible ="ind", 
                #geom.var = c("point", "text"),
                alpha.var=0.3,
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Outer (1)",  "Inner (3)"),               
                ,                pointsize = 2,
                # xlim=c(-4, 4) ,
                #  ylim=c(-4, 4) ,
                title = "Contributive arrows for physico-chemical data",
                palette = c("lightgoldenrod4",  "mistyrose4"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")



############################################################################################################
### PCA en fonction TAXO
############################################################################################################

#dat2 <- dat[, -1]
#dat2 <- dat[, -c(1:25) ]
dat2 <- dat[, -c(1:12) ]
dat2$Cyanobacteria <- NULL

row.names(dat2) <- c("Outer (1)",  "Inner (3)")

res.pca <- dudi.pca(dat2,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5            # Nombre d'axes gardés
)

str(dat)

#dat$index  <- factor(dat$index, levels = c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03"))
#unique(dat$index)

fviz_pca_biplot(res.pca, 
                label = "all", 
                col.ind = c("Outer (1)",  "Inner (3)"),              
                pointsize = 2,
                #xlim=c(-5, 5) ,
                # ylim=c(-5, 5) ,
                repel=T,
                title = "PCA - Biplot for phyla data",
                palette = c("lightgoldenrod4",  "mistyrose4"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

fviz_pca_biplot(res.pca, 
                label = "ind", 
                invisible ="var", 
                # geom=c("text", "point"),
                col.ind = c("Outer (1)",  "Inner (3)"),                pointsize = 2,
                repel=T,
                # xlim=c(-5, 5) ,
                # ylim=c(-5, 5) ,
                title = "PCA for phyla data",
                palette = c("lightgoldenrod4",  "mistyrose4"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

options(ggrepel.max.overlaps = Inf)



fviz_pca_biplot(res.pca, 
                label = "var", 
                invisible ="ind", 
                #geom.var = c("point", "text"),
                alpha.var=0.3,
                repel=TRUE,
                # geom=c("text", "point"),
                col.ind = c("Outer (1)",  "Inner (3)"),                pointsize = 2,
                #xlim=c(-4, 4) ,
                # ylim=c(-4, 4) ,
                title = "Contributive arrows for phyla data",
                palette = c("lightgoldenrod4",  "mistyrose4"),
                ggtheme = theme_minimal()) + 
  scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16)) + 
  theme_bw() + 
  theme(legend.position = "none")

















############################################################################################################
######### db-RDA Par Phylum
############################################################################################################

# If data set without seawater
#taxa_bacteria <- dat[, -c(1:25) ]

taxa_bacteria <- dat[, -c(1:12) ]


# If data set without seawater
#physico_chimic <- dat[, c(2:25) ]
physico_chimic <- dat[, c(2:12) ]



head(physico_chimic)
head(taxa_bacteria)

######### db-RDA for bacteria
############################################################################################################

# Beta Diversity
# from https://www.rpubs.com/roalle/mres_2019

# Bray-Curtis Dissimilarity

# SUPPRESSION DES LIGNES OU IL N'Y A QUE DES ZEROS
taxa_bacteria_02 <- taxa_bacteria[-c(7:12),]
taxa_bacteria_02

# Calculate pairwise dissimilarity
taxa_bacteria.mdf <- as.matrix.data.frame(taxa_bacteria_02)
#rownames(taxa_bacteria.mdf) <- env.data$site

taxa_bacteria.bray <- vegdist(taxa_bacteria.mdf, method = "bray")
taxa_bacteria.bray

# PERMANOVA (permutational multivariate analysis of variance; Anderson 2001)
# Homogeneity of dispersion test
#permutest(betadisper(taxa_bacteria.bray, physico_chimic$ph_h2o))

# PERMANOVA analysis
# adonis(taxa_bacteria.bray ~ ph_h2o + ph_kcl + N_NO3_mg_kg + N_NH4_mg_kg + N_g_kg + C_g_kg + P_mg_kg + cations_change_Ca + cations_change_Mg + cations_change_Na + cations_change_K + change_capacity, data = physico_chimic, permutations = 999)

# Distance-Based Redundancy Analysis (dbRDA)
# from https://www.rpubs.com/roalle/mres_2019

# Distance-Based Redundancy Analysis (dbRDA; Legendre and Anderson, 1999) is an extension of Redundancy Analysis (RDA) which allows the use of 
# non-euclidean distance matrices as inputs (e.g. Bray-Curtis dissimilarity). The method works by first calculating a PCoA from the dissimilarity matrix, 
# and then subjecting the PCoA eigenvalues (which represent dissimilarities in euclidean space) to RDA. The method aims to detect linear relationships between 
# environmental variables and these dissimilarities. 
# dbRDA differs from PCoA in that it is a constrained analysis. 
# While PCoA axes are generated to explain maximum variation in the distance matrix, 
# dbRDA canonical axes are constructed as linear combinations of environmental variables. 
# Consequently these axes are constrained to the environmental variables in the model, and the ordination will be distinct from the PCoA. dbRDA allows us to 
# visualise how environmental variables contrain variation in community composition between our sites. 
# As dbRDA considers multiple environmental variables which are measured using different techniques, it is important to normalise environmental variables such 
# that they can be compared concurrently. This is achieved using a z-score transformation. Moreover, dbRDA is sensitive to multicolinearity 
# (i.e. high correlation between environmental variables). We can test whether any of our environmental variables are correlated by calculating variance 
# inflation factors (VIFs) and then removing selected terms from the model until all VIF scores are < 10.

# Z-Score Transformation
# Z-Score transform environmental data :
# The z-score transformation normalises environmental variables by their standard deviation from the mean 

physico_chimic

# IDEM : 
# SUPPRESSION DES LIGNES OU IL N'Y A QUE DES ZEROS
physico_chimic_02 <- physico_chimic[-c(7:12),]
physico_chimic_02
env.data.z <- physico_chimic_02

env.data.z$pH <- (env.data.z$pH - mean(env.data.z$pH))/sd(env.data.z$pH)
env.data.z$SiOH <- (env.data.z$SiOH - mean(env.data.z$SiOH))/sd(env.data.z$SiOH)
env.data.z$Nox <- (env.data.z$Nox - mean(env.data.z$Nox))/sd(env.data.z$Nox)
env.data.z$PO4 <- (env.data.z$PO4 - mean(env.data.z$PO4))/sd(env.data.z$PO4)
env.data.z$NH4 <- (env.data.z$NH4 - mean(env.data.z$NH4))/sd(env.data.z$NH4)
env.data.z$DIC <- (env.data.z$DIC - mean(env.data.z$DIC))/sd(env.data.z$DIC)
env.data.z$DOC <- (env.data.z$DOC - mean(env.data.z$DOC))/sd(env.data.z$DOC)
env.data.z$POC <- (env.data.z$POC - mean(env.data.z$POC))/sd(env.data.z$POC)
env.data.z$PN <- (env.data.z$PN - mean(env.data.z$PN))/sd(env.data.z$PN)
env.data.z$Temp <- (env.data.z$Temp - mean(env.data.z$Temp))/sd(env.data.z$Temp)
env.data.z$DIN <- (env.data.z$DIN - mean(env.data.z$DIN))/sd(env.data.z$DIN)

# Perform dbRDA
# construct full model and calculate VIF

dbRDA.full <- capscale(taxa_bacteria.bray ~ pH +
                         SiOH +
                         Nox +
                         PO4 +
                         NH4 +
                         DOC +
                         DIC +
                         POC +
                         PN +
                         Temp +
                         DIN,
                       env.data.z)
vif.cca(dbRDA.full)

dbRDA.full <- capscale(taxa_bacteria.mdf ~  DOC +
                         SiOH,
                       env.data.z)
vif.cca(dbRDA.full)

# VIF explaination, according http://jsnfri.fra.affrc.go.jp/gunshu/3dbrdaE.html
# Multiple-colinearlity between variables are examined by variance inflation factor (VIF). 
# When VIF > 10, problems of multiple-colinearlity become severe. So, remove variables of large VIF, one by one.

# construct reduced model and calculate VIF
dbRDA.mat <- capscale(taxa_bacteria.mdf ~ DOC +
                        SiOH,
                      env.data.z)

sort(vif.cca(dbRDA.mat))


# test overall significance of the analysis
anova(dbRDA.mat, permutations = 999)

# test significance of each environmental variable
anova(dbRDA.mat, by = "terms")

# summary of dbRDA model to extract total variance constrained and axis scores
summary(dbRDA.mat)

# Plot dbRDA

smry <- summary(dbRDA.mat)
scrs <- scores(dbRDA.mat)
df1  <- data.frame(smry$sites[,1:2]) # site scores for RDA1 and RDA2
df1$site <- rownames(df1)  #add site names
df2  <- data.frame(smry$biplot[,1:2])  # mapping environmental variables

site <- c("Mycale", "Mycale", "Mycale", "Rhabdastrella globostellata", "Rhabdastrella globostellata", "Rhabdastrella globostellata", "SeaWater", "SeaWater")
ind <- c("Mycale 01", "Mycale 02", "Mycale 03", "Rhabdastrella globostellata 01", "Rhabdastrella globostellata 02", "Rhabdastrella globostellata 03", "SeaWater 01", "SeaWater 02")
df_col <- data.frame(ind, site)
df_col$site <- factor(df_col$site, levels = c("Mycale", "Rhabdastrella globostellata", "SeaWater"))
head(df_col)



legend_title <- "Conditions"

ggplot(df1, aes(x=CAP1, y=CAP2, colour = df_col$site, shape=df_col$site)) + 
  geom_point(size = 3) +
  ggtitle("Distance-Based Redundancy Analysis") +
  scale_colour_manual(legend_title,labels=c("Mycale", "Rhabdastrella globostellata", "SeaWater"), values = c("coral3",  "darkgoldenrod4", "cadetblue3")) +
  scale_shape_manual(legend_title,labels=c("Mycale", "Rhabdastrella globostellata", "SeaWater"), values = c(16, 16, 16)) +
  scale_fill_manual(legend_title,labels=c("Mycale", "Rhabdastrella globostellata", "SeaWater"), values = c("coral3",  "darkgoldenrod4", "cadetblue3")) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #stat_ellipse() +
  # xlim(-50, 50) +
  # ylim(-50, 50) +
  xlab("db-RDA1 (92.97%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("db-RDA2 (1.78%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  coord_fixed() +
  # ADD ggforce's ellipses
  ggforce::geom_mark_ellipse(aes(fill = df_col$site,
                                 color = df_col$site)) +
  # guides(fill=guide_legend("my awesome title")) +
  theme_bw() #+ theme(legend.position="none")









#install.packages("ggrepel")
library(ggrepel)

ggplot(df1, aes(x=CAP1, y=CAP2, colour = df_col$site)) + 
  geom_segment(data=df2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               # color="black", arrow=arrow(length=unit(0.01,"npc"))) +
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  ggrepel:: geom_text_repel(data=df2, aes(x=CAP1,y=CAP2,label=rownames(df2)), 
                            #    color="grey50", size=3) +
                            color="black", size=5) +
  #  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  # scale_colour_manual(values = c( "firebrick", "goldenrod", "olivedrab1", "green", "forestgreen")) +
  # geom_text(aes(label=rownames(df1),
  #                hjust=0,vjust=1.5), colour = "black",size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #  stat_ellipse() +
  xlim(-2, 1.1) +
  # ylim(-5, 6) +
  xlab("RDA1 (92.97%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("RDA2 (1.78%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  #coord_fixed() +
  theme_bw()


