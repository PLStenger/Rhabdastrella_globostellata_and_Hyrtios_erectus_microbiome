# Change your stwd()
# The CSV/txt files come from the alpha-rarefaction.qzv (click on "Download CSV" when you are using QIIME2 view online tool))

setwd("~/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/07_01_Rhabdastrella_globostellata_and_Mycale_microbiome/Rhabdastrella_globostellata_and_Mycale_microbiome/06_rarefaction_curves")

library(ggplot2)
library(stringr)
library(plyr)



############################################################################################################################################
############################################################################################################################################
############ V4
############################################################################################################################################
############################################################################################################################################

######################################################################
############ SHANNON
######################################################################

data <- read.table("shannon.csv", sep=",", header=T)
rownames(data) <- data[,1]

p <- list()
for (i in 1:15) {
  
  sequencing_depth <- rownames(data.frame(t(data[i,])))
  shannon <- c(data.frame(t(data[i,]))[,1])
  data_02 <- data.frame(sequencing_depth, shannon)
  data_02 <- data_02[-c(1, 102, 103, 104, 105),]
  data_02$sequencing_depth <- sub("depth.", "", data_02$sequencing_depth)
  data_02$sequencing_depth <- str_split_fixed(data_02$sequencing_depth, "_iter", 2)[,1]
  #data_02$sequencing_depth <- as.factor(data_02$sequencing_depth)
  data_02$sequencing_depth <- as.numeric(data_02$sequencing_depth)
  data_02$shannon <- as.numeric(data_02$shannon)
  str(data_02)
  
  p[[i]] <- ddply(data_02, .(sequencing_depth), summarise, shannon = mean(shannon))
}


my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)


#pdf(file = "alpha-rarefaction__shannon.pdf",   
#    width = 10, 
#    height = 4)

png("alpha-rarefaction_shannon.png", 1000, 400)


ggplot() +
  
  geom_point(data = p[[1]], aes(as.factor(sequencing_depth), shannon), colour="deepskyblue4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[1]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "deepskyblue4",
            position = position_dodge(width=0.75)) + 
  geom_point(data = p[[2]], aes(as.factor(sequencing_depth), shannon), colour="cadetblue4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[2]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "cadetblue4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[3]], aes(as.factor(sequencing_depth), shannon), colour="chartreuse4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[3]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "chartreuse4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[4]], aes(as.factor(sequencing_depth), shannon), colour="chocolate4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[4]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "chocolate4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[5]], aes(as.factor(sequencing_depth), shannon), colour="darkgoldenrod4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[5]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "darkgoldenrod4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[6]], aes(as.factor(sequencing_depth), shannon), colour="firebrick4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[6]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "firebrick4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[7]], aes(as.factor(sequencing_depth), shannon), colour="forestgreen", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[7]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "forestgreen",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[8]], aes(as.factor(sequencing_depth), shannon), colour="deeppink3", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[8]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "deeppink3",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[9]], aes(as.factor(sequencing_depth), shannon), colour="gray4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[9]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "gray4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[10]], aes(as.factor(sequencing_depth), shannon), colour="chocolate1", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[10]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "chocolate1",
            position = position_dodge(width=0.75)) +
  # fake axes (x-axis stops at year 2009, y-axis stops at value 45)
  annotate("segment", x = 0, y = 0, xend = 15, yend = 0, colour = "white") +
  # annotate("segment", x = 1991, y = 5, xend = 1991, yend = 45) + 
  annotate("text", label = c(data.frame(t(data[1,]))[,1])[1], x = replace(min(which(is.na(p[[1 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[1 ]]$shannon)))+0.5), 10.5), y = my.max(p[[1]]$shannon), colour = "deepskyblue4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[2,]))[,1])[1], x = replace(min(which(is.na(p[[2 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[2 ]]$shannon)))+0.5), 10.5), y = my.max(p[[2]]$shannon), colour = "cadetblue4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[3,]))[,1])[1], x = replace(min(which(is.na(p[[3 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[3 ]]$shannon)))+0.5), 10.5), y = my.max(p[[3]]$shannon), colour = "chartreuse4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[4,]))[,1])[1], x = replace(min(which(is.na(p[[4 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[4 ]]$shannon)))+0.5), 10.5), y = my.max(p[[4]]$shannon), colour = "chocolate4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[5,]))[,1])[1], x = replace(min(which(is.na(p[[5 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[5 ]]$shannon)))+0.5), 10.5), y = my.max(p[[5]]$shannon), colour = "darkgoldenrod4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[6,]))[,1])[1], x = replace(min(which(is.na(p[[6 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[6 ]]$shannon)))+0.5), 10.5), y = my.max(p[[6]]$shannon), colour = "firebrick4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[7,]))[,1])[1], x = replace(min(which(is.na(p[[7 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[7 ]]$shannon)))+0.5), 10.5), y = my.max(p[[7]]$shannon), colour = "forestgreen", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[8,]))[,1])[1], x = replace(min(which(is.na(p[[8 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[8 ]]$shannon)))+0.5), 10.5), y = my.max(p[[8]]$shannon), colour = "deeppink3", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[9,]))[,1])[1], x = replace(min(which(is.na(p[[9 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[9 ]]$shannon)))+0.5), 10.5), y = my.max(p[[9]]$shannon), colour = "gray4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[10,]))[,1])[1], x = replace(min(which(is.na(p[[10]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[10]]$shannon)))+0.5), 10.5), y = my.max(p[[10]]$shannon), colour = "chocolate1", hjust = 0) +

  scale_x_discrete("sequencing depth") +
  scale_y_continuous("shannon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("alpha-rarefaction Shannon")

dev.off()



######################################################################
############ FAITH_PD
######################################################################

data <- read.table("faith_pd.csv", sep=",", header=T)
rownames(data) <- data[,1]

p <- list()
for (i in 1:15) {
  
  sequencing_depth <- rownames(data.frame(t(data[i,])))
  faith_pd <- c(data.frame(t(data[i,]))[,1])
  data_02 <- data.frame(sequencing_depth, faith_pd)
  data_02 <- data_02[-c(1, 102, 103, 104, 105),]
  data_02$sequencing_depth <- sub("depth.", "", data_02$sequencing_depth)
  data_02$sequencing_depth <- str_split_fixed(data_02$sequencing_depth, "_iter", 2)[,1]
  #data_02$sequencing_depth <- as.factor(data_02$sequencing_depth)
  data_02$sequencing_depth <- as.numeric(data_02$sequencing_depth)
  data_02$faith_pd <- as.numeric(data_02$faith_pd)
  str(data_02)
  
  p[[i]] <- ddply(data_02, .(sequencing_depth), summarise, faith_pd = mean(faith_pd))
}

#pdf(file = "alpha-rarefaction_faith_pd.pdf",   
#    width = 10, 
#    height = 4)

png("alpha-rarefaction_faith_pd.png", 1000, 400)

ggplot() +
  
  geom_point(data = p[[1]], aes(as.factor(sequencing_depth), faith_pd), colour="deepskyblue4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[1]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "deepskyblue4",
            position = position_dodge(width=0.75)) + 
  geom_point(data = p[[2]], aes(as.factor(sequencing_depth), faith_pd), colour="cadetblue4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[2]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "cadetblue4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[3]], aes(as.factor(sequencing_depth), faith_pd), colour="chartreuse4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[3]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "chartreuse4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[4]], aes(as.factor(sequencing_depth), faith_pd), colour="chocolate4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[4]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "chocolate4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[5]], aes(as.factor(sequencing_depth), faith_pd), colour="darkgoldenrod4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[5]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "darkgoldenrod4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[6]], aes(as.factor(sequencing_depth), faith_pd), colour="firebrick4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[6]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "firebrick4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[7]], aes(as.factor(sequencing_depth), faith_pd), colour="forestgreen", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[7]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "forestgreen",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[8]], aes(as.factor(sequencing_depth), faith_pd), colour="deeppink3", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[8]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "deeppink3",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[9]], aes(as.factor(sequencing_depth), faith_pd), colour="gray4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[9]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "gray4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[10]], aes(as.factor(sequencing_depth), faith_pd), colour="chocolate1", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[10]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "chocolate1",
            position = position_dodge(width=0.75)) +
  # fake axes (x-axis stops at year 2009, y-axis stops at value 45)
  annotate("segment", x = 0, y = 0, xend = 15, yend = 0, colour = "white") +
  # annotate("segment", x = 1991, y = 5, xend = 1991, yend = 45) + 
  annotate("text", label = c(data.frame(t(data[1,]))[,1])[1], x = replace(min(which(is.na(p[[1 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[1 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[1]]$faith_pd), colour = "deepskyblue4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[2,]))[,1])[1], x = replace(min(which(is.na(p[[2 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[2 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[2]]$faith_pd), colour = "cadetblue4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[3,]))[,1])[1], x = replace(min(which(is.na(p[[3 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[3 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[3]]$faith_pd), colour = "chartreuse4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[4,]))[,1])[1], x = replace(min(which(is.na(p[[4 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[4 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[4]]$faith_pd), colour = "chocolate4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[5,]))[,1])[1], x = replace(min(which(is.na(p[[5 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[5 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[5]]$faith_pd), colour = "darkgoldenrod4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[6,]))[,1])[1], x = replace(min(which(is.na(p[[6 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[6 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[6]]$faith_pd), colour = "firebrick4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[7,]))[,1])[1], x = replace(min(which(is.na(p[[7 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[7 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[7]]$faith_pd), colour = "forestgreen", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[8,]))[,1])[1], x = replace(min(which(is.na(p[[8 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[8 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[8]]$faith_pd), colour = "deeppink3", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[9,]))[,1])[1], x = replace(min(which(is.na(p[[9 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[9 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[9]]$faith_pd), colour = "gray4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[10,]))[,1])[1], x = replace(min(which(is.na(p[[10]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[10]]$faith_pd)))+0.5), 10.5), y = my.max(p[[10]]$faith_pd), colour = "chocolate1", hjust = 0) +

  scale_x_discrete("sequencing depth") +
  scale_y_continuous("faith_pd") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("alpha-rarefaction faith_pd")

dev.off()

######################################################################
############ OBSERVED OTU
######################################################################

data <- read.table("observed_features.csv", sep=",", header=T)
rownames(data) <- data[,1]

p <- list()
for (i in 1:15) {
  
  sequencing_depth <- rownames(data.frame(t(data[i,])))
  observed_otu <- c(data.frame(t(data[i,]))[,1])
  data_02 <- data.frame(sequencing_depth, observed_otu)
  data_02 <- data_02[-c(1, 102, 103, 104, 105),]
  data_02$sequencing_depth <- sub("depth.", "", data_02$sequencing_depth)
  data_02$sequencing_depth <- str_split_fixed(data_02$sequencing_depth, "_iter", 2)[,1]
  #data_02$sequencing_depth <- as.factor(data_02$sequencing_depth)
  data_02$sequencing_depth <- as.numeric(data_02$sequencing_depth)
  data_02$observed_otu <- as.numeric(data_02$observed_otu)
  str(data_02)
  
  p[[i]] <- ddply(data_02, .(sequencing_depth), summarise, observed_otu = mean(observed_otu))
}


#pdf(file = "alpha-rarefaction_observed_otu.pdf",   
#    width = 10, 
#    height = 4)

png("alpha-rarefaction_observed_otu.png", 1000, 400)

ggplot() +
  
  geom_point(data = p[[1]], aes(as.factor(sequencing_depth), observed_otu), colour="coral3", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[1]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "coral3",
            position = position_dodge(width=0.75)) + 
  geom_point(data = p[[2]], aes(as.factor(sequencing_depth), observed_otu), colour="chocolate3", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[2]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "chocolate3",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[3]], aes(as.factor(sequencing_depth), observed_otu), colour="chocolate1", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[3]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "chocolate1",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[4]], aes(as.factor(sequencing_depth), observed_otu), colour="darkgoldenrod4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[4]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "darkgoldenrod4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[5]], aes(as.factor(sequencing_depth), observed_otu), colour="darkgoldenrod3", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[5]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "darkgoldenrod3",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[6]], aes(as.factor(sequencing_depth), observed_otu), colour="darkgoldenrod", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[6]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "darkgoldenrod",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[7]], aes(as.factor(sequencing_depth), observed_otu), colour="cadetblue3", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[7]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "cadetblue3",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[8]], aes(as.factor(sequencing_depth), observed_otu), colour="cadetblue3", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[8]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "cadetblue3",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[9]], aes(as.factor(sequencing_depth), observed_otu), colour="burlywood4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[9]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "burlywood4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[10]], aes(as.factor(sequencing_depth), observed_otu), colour="antiquewhite4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[10]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "antiquewhite4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[11]], aes(as.factor(sequencing_depth), observed_otu), colour="bisque4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[11]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "bisque4",
            position = position_dodge(width=0.75)) +
  # fake axes (x-axis stops at year 2009, y-axis stops at value 45)
  annotate("segment", x = 0, y = 0, xend = 15, yend = 0, colour = "white") +
  # annotate("segment", x = 1991, y = 5, xend = 1991, yend = 45) + 
  annotate("text", label = c(data.frame(t(data[1,]))[,1])[1], x = replace(min(which(is.na(p[[1 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[1 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[1]]$observed_otu), colour = "coral3", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[2,]))[,1])[1], x = replace(min(which(is.na(p[[2 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[2 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[2]]$observed_otu), colour = "chocolate3", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[3,]))[,1])[1], x = replace(min(which(is.na(p[[3 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[3 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[3]]$observed_otu), colour = "chocolate1", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[4,]))[,1])[1], x = replace(min(which(is.na(p[[4 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[4 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[4]]$observed_otu), colour = "darkgoldenrod4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[5,]))[,1])[1], x = replace(min(which(is.na(p[[5 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[5 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[5]]$observed_otu), colour = "darkgoldenrod3", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[6,]))[,1])[1], x = replace(min(which(is.na(p[[6 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[6 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[6]]$observed_otu), colour = "darkgoldenrod", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[7,]))[,1])[1], x = replace(min(which(is.na(p[[7 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[7 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[7]]$observed_otu), colour = "cadetblue3", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[8,]))[,1])[1], x = replace(min(which(is.na(p[[8 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[8 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[8]]$observed_otu), colour = "cadetblue3", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[9,]))[,1])[1], x = replace(min(which(is.na(p[[9 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[9 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[9]]$observed_otu), colour = "burlywood4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[10,]))[,1])[1], x = replace(min(which(is.na(p[[10]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[10]]$observed_otu)))+0.5), 10.5), y = my.max(p[[10]]$observed_otu), colour = "antiquewhite4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[11,]))[,1])[1], x = replace(min(which(is.na(p[[11]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[11]]$observed_otu)))+0.5), 10.5), y = my.max(p[[11]]$observed_otu), colour = "bisque4", hjust = 0) +
  
  scale_x_discrete("Sequencing depth") +
  scale_y_continuous("Observed ASVs") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("alpha-rarefaction Observed ASVs")

dev.off()










############################################################################################################################################
############################################################################################################################################
############ Only V1V3
############################################################################################################################################
############################################################################################################################################

######################################################################
############ SHANNON
######################################################################

data <- read.table("alpha-rarefaction_V4/shannon_VIV3.csv", sep=",", header=T)
rownames(data) <- data[,1]

p <- list()
for (i in 1:15) {
  
  sequencing_depth <- rownames(data.frame(t(data[i,])))
  shannon <- c(data.frame(t(data[i,]))[,1])
  data_02 <- data.frame(sequencing_depth, shannon)
  data_02 <- data_02[-c(1, 102, 103, 104, 105),]
  data_02$sequencing_depth <- sub("depth.", "", data_02$sequencing_depth)
  data_02$sequencing_depth <- str_split_fixed(data_02$sequencing_depth, "_iter", 2)[,1]
  #data_02$sequencing_depth <- as.factor(data_02$sequencing_depth)
  data_02$sequencing_depth <- as.numeric(data_02$sequencing_depth)
  data_02$shannon <- as.numeric(data_02$shannon)
  str(data_02)
  
  p[[i]] <- ddply(data_02, .(sequencing_depth), summarise, shannon = mean(shannon))
}


my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)


#pdf(file = "alpha-rarefaction__shannon.pdf",   
#    width = 10, 
#    height = 4)

png("alpha-rarefaction_shannon_V1V3.png", 1000, 400)


ggplot() +
  
  geom_point(data = p[[1]], aes(as.factor(sequencing_depth), shannon), colour="deepskyblue4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[1]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "deepskyblue4",
            position = position_dodge(width=0.75)) + 
  geom_point(data = p[[2]], aes(as.factor(sequencing_depth), shannon), colour="cadetblue4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[2]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "cadetblue4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[3]], aes(as.factor(sequencing_depth), shannon), colour="chartreuse4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[3]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "chartreuse4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[4]], aes(as.factor(sequencing_depth), shannon), colour="chocolate4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[4]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "chocolate4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[5]], aes(as.factor(sequencing_depth), shannon), colour="darkgoldenrod4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[5]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "darkgoldenrod4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[6]], aes(as.factor(sequencing_depth), shannon), colour="firebrick4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[6]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "firebrick4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[7]], aes(as.factor(sequencing_depth), shannon), colour="forestgreen", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[7]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "forestgreen",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[8]], aes(as.factor(sequencing_depth), shannon), colour="deeppink3", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[8]], aes(as.factor(sequencing_depth), shannon, group = 1), colour = "deeppink3",
            position = position_dodge(width=0.75)) +
  # fake axes (x-axis stops at year 2009, y-axis stops at value 45)
  annotate("segment", x = 0, y = 0, xend = 15, yend = 0, colour = "white") +
  # annotate("segment", x = 1991, y = 5, xend = 1991, yend = 45) + 
  annotate("text", label = c(data.frame(t(data[1,]))[,1])[1], x = replace(min(which(is.na(p[[1 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[1 ]]$shannon)))+0.5), 10.5), y = my.max(p[[1]]$shannon), colour = "deepskyblue4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[2,]))[,1])[1], x = replace(min(which(is.na(p[[2 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[2 ]]$shannon)))+0.5), 10.5), y = my.max(p[[2]]$shannon), colour = "cadetblue4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[3,]))[,1])[1], x = replace(min(which(is.na(p[[3 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[3 ]]$shannon)))+0.5), 10.5), y = my.max(p[[3]]$shannon), colour = "chartreuse4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[4,]))[,1])[1], x = replace(min(which(is.na(p[[4 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[4 ]]$shannon)))+0.5), 10.5), y = my.max(p[[4]]$shannon), colour = "chocolate4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[5,]))[,1])[1], x = replace(min(which(is.na(p[[5 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[5 ]]$shannon)))+0.5), 10.5), y = my.max(p[[5]]$shannon), colour = "darkgoldenrod4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[6,]))[,1])[1], x = replace(min(which(is.na(p[[6 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[6 ]]$shannon)))+0.5), 10.5), y = my.max(p[[6]]$shannon), colour = "firebrick4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[7,]))[,1])[1], x = replace(min(which(is.na(p[[7 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[7 ]]$shannon)))+0.5), 10.5), y = my.max(p[[7]]$shannon), colour = "forestgreen", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[8,]))[,1])[1], x = replace(min(which(is.na(p[[8 ]]$shannon)))+0.5, is.infinite(min(which(is.na(p[[8 ]]$shannon)))+0.5), 10.5), y = my.max(p[[8]]$shannon), colour = "deeppink3", hjust = 0) +

  scale_x_discrete("sequencing depth") +
  scale_y_continuous("shannon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("alpha-rarefaction Shannon")

dev.off()



######################################################################
############ FAITH_PD
######################################################################

data <- read.table("alpha-rarefaction_V4/faith_pd.csv", sep=",", header=T)
rownames(data) <- data[,1]

p <- list()
for (i in 1:15) {
  
  sequencing_depth <- rownames(data.frame(t(data[i,])))
  faith_pd <- c(data.frame(t(data[i,]))[,1])
  data_02 <- data.frame(sequencing_depth, faith_pd)
  data_02 <- data_02[-c(1, 102, 103, 104, 105),]
  data_02$sequencing_depth <- sub("depth.", "", data_02$sequencing_depth)
  data_02$sequencing_depth <- str_split_fixed(data_02$sequencing_depth, "_iter", 2)[,1]
  #data_02$sequencing_depth <- as.factor(data_02$sequencing_depth)
  data_02$sequencing_depth <- as.numeric(data_02$sequencing_depth)
  data_02$faith_pd <- as.numeric(data_02$faith_pd)
  str(data_02)
  
  p[[i]] <- ddply(data_02, .(sequencing_depth), summarise, faith_pd = mean(faith_pd))
}

#pdf(file = "alpha-rarefaction_faith_pd.pdf",   
#    width = 10, 
#    height = 4)

png("alpha-rarefaction_faith_pd_V4.png", 1000, 400)

ggplot() +
  
  geom_point(data = p[[1]], aes(as.factor(sequencing_depth), faith_pd), colour="deepskyblue4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[1]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "deepskyblue4",
            position = position_dodge(width=0.75)) + 
  geom_point(data = p[[2]], aes(as.factor(sequencing_depth), faith_pd), colour="cadetblue4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[2]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "cadetblue4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[3]], aes(as.factor(sequencing_depth), faith_pd), colour="chartreuse4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[3]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "chartreuse4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[4]], aes(as.factor(sequencing_depth), faith_pd), colour="chocolate4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[4]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "chocolate4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[5]], aes(as.factor(sequencing_depth), faith_pd), colour="darkgoldenrod4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[5]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "darkgoldenrod4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[6]], aes(as.factor(sequencing_depth), faith_pd), colour="firebrick4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[6]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "firebrick4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[7]], aes(as.factor(sequencing_depth), faith_pd), colour="forestgreen", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[7]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "forestgreen",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[8]], aes(as.factor(sequencing_depth), faith_pd), colour="deeppink3", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[8]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "deeppink3",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[9]], aes(as.factor(sequencing_depth), faith_pd), colour="gray4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[9]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "gray4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[10]], aes(as.factor(sequencing_depth), faith_pd), colour="chocolate1", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[10]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "chocolate1",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[11]], aes(as.factor(sequencing_depth), faith_pd), colour="indianred2", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[11]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "indianred2",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[12]], aes(as.factor(sequencing_depth), faith_pd), colour="goldenrod4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[12]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "goldenrod4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[13]], aes(as.factor(sequencing_depth), faith_pd), colour="blue3", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[13]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "blue3",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[14]], aes(as.factor(sequencing_depth), faith_pd), colour="cornflowerblue", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[14]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "cornflowerblue",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[15]], aes(as.factor(sequencing_depth), faith_pd), colour="coral4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[15]], aes(as.factor(sequencing_depth), faith_pd, group = 1), colour = "coral4",
            position = position_dodge(width=0.75)) +
  # fake axes (x-axis stops at year 2009, y-axis stops at value 45)
  annotate("segment", x = 0, y = 0, xend = 15, yend = 0, colour = "white") +
  # annotate("segment", x = 1991, y = 5, xend = 1991, yend = 45) + 
  annotate("text", label = c(data.frame(t(data[1,]))[,1])[1], x = replace(min(which(is.na(p[[1 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[1 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[1]]$faith_pd), colour = "deepskyblue4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[2,]))[,1])[1], x = replace(min(which(is.na(p[[2 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[2 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[2]]$faith_pd), colour = "cadetblue4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[3,]))[,1])[1], x = replace(min(which(is.na(p[[3 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[3 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[3]]$faith_pd), colour = "chartreuse4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[4,]))[,1])[1], x = replace(min(which(is.na(p[[4 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[4 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[4]]$faith_pd), colour = "chocolate4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[5,]))[,1])[1], x = replace(min(which(is.na(p[[5 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[5 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[5]]$faith_pd), colour = "darkgoldenrod4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[6,]))[,1])[1], x = replace(min(which(is.na(p[[6 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[6 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[6]]$faith_pd), colour = "firebrick4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[7,]))[,1])[1], x = replace(min(which(is.na(p[[7 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[7 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[7]]$faith_pd), colour = "forestgreen", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[8,]))[,1])[1], x = replace(min(which(is.na(p[[8 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[8 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[8]]$faith_pd), colour = "deeppink3", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[9,]))[,1])[1], x = replace(min(which(is.na(p[[9 ]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[9 ]]$faith_pd)))+0.5), 10.5), y = my.max(p[[9]]$faith_pd), colour = "gray4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[10,]))[,1])[1], x = replace(min(which(is.na(p[[10]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[10]]$faith_pd)))+0.5), 10.5), y = my.max(p[[10]]$faith_pd), colour = "chocolate1", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[11,]))[,1])[1], x = replace(min(which(is.na(p[[11]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[11]]$faith_pd)))+0.5), 10.5), y = my.max(p[[11]]$faith_pd), colour = "indianred2", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[12,]))[,1])[1], x = replace(min(which(is.na(p[[12]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[12]]$faith_pd)))+0.5), 10.5), y = my.max(p[[12]]$faith_pd), colour = "goldenrod4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[13,]))[,1])[1], x = replace(min(which(is.na(p[[13]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[13]]$faith_pd)))+0.5), 10.5), y = my.max(p[[13]]$faith_pd), colour = "blue3", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[14,]))[,1])[1], x = replace(min(which(is.na(p[[14]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[14]]$faith_pd)))+0.5), 10.5), y = my.max(p[[14]]$faith_pd), colour = "cornflowerblue", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[15,]))[,1])[1], x = replace(min(which(is.na(p[[15]]$faith_pd)))+0.5, is.infinite(min(which(is.na(p[[15]]$faith_pd)))+0.5), 10.5), y = my.max(p[[15]]$faith_pd), colour = "coral4", hjust = 0) +
  
  scale_x_discrete("sequencing depth") +
  scale_y_continuous("faith_pd") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("alpha-rarefaction faith_pd")

dev.off()

######################################################################
############ OBSERVED OTU
######################################################################

data <- read.table("alpha-rarefaction_V4/observed_features.csv", sep=",", header=T)
rownames(data) <- data[,1]

p <- list()
for (i in 1:15) {
  
  sequencing_depth <- rownames(data.frame(t(data[i,])))
  observed_otu <- c(data.frame(t(data[i,]))[,1])
  data_02 <- data.frame(sequencing_depth, observed_otu)
  data_02 <- data_02[-c(1, 102, 103, 104, 105),]
  data_02$sequencing_depth <- sub("depth.", "", data_02$sequencing_depth)
  data_02$sequencing_depth <- str_split_fixed(data_02$sequencing_depth, "_iter", 2)[,1]
  #data_02$sequencing_depth <- as.factor(data_02$sequencing_depth)
  data_02$sequencing_depth <- as.numeric(data_02$sequencing_depth)
  data_02$observed_otu <- as.numeric(data_02$observed_otu)
  str(data_02)
  
  p[[i]] <- ddply(data_02, .(sequencing_depth), summarise, observed_otu = mean(observed_otu))
}


#pdf(file = "alpha-rarefaction_observed_otu.pdf",   
#    width = 10, 
#    height = 4)

png("alpha-rarefaction_observed_otu_V4.png", 1000, 400)

ggplot() +
  
  geom_point(data = p[[1]], aes(as.factor(sequencing_depth), observed_otu), colour="deepskyblue4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[1]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "deepskyblue4",
            position = position_dodge(width=0.75)) + 
  geom_point(data = p[[2]], aes(as.factor(sequencing_depth), observed_otu), colour="cadetblue4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[2]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "cadetblue4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[3]], aes(as.factor(sequencing_depth), observed_otu), colour="chartreuse4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[3]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "chartreuse4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[4]], aes(as.factor(sequencing_depth), observed_otu), colour="chocolate4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[4]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "chocolate4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[5]], aes(as.factor(sequencing_depth), observed_otu), colour="darkgoldenrod4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[5]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "darkgoldenrod4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[6]], aes(as.factor(sequencing_depth), observed_otu), colour="firebrick4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[6]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "firebrick4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[7]], aes(as.factor(sequencing_depth), observed_otu), colour="forestgreen", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[7]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "forestgreen",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[8]], aes(as.factor(sequencing_depth), observed_otu), colour="deeppink3", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[8]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "deeppink3",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[9]], aes(as.factor(sequencing_depth), observed_otu), colour="gray4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[9]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "gray4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[10]], aes(as.factor(sequencing_depth), observed_otu), colour="chocolate1", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[10]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "chocolate1",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[11]], aes(as.factor(sequencing_depth), observed_otu), colour="indianred2", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[11]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "indianred2",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[12]], aes(as.factor(sequencing_depth), observed_otu), colour="goldenrod4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[12]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "goldenrod4",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[13]], aes(as.factor(sequencing_depth), observed_otu), colour="blue3", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[13]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "blue3",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[14]], aes(as.factor(sequencing_depth), observed_otu), colour="cornflowerblue", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[14]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "cornflowerblue",
            position = position_dodge(width=0.75)) +
  geom_point(data = p[[15]], aes(as.factor(sequencing_depth), observed_otu), colour="coral4", 
             position = position_dodge(width=0.75)) +
  geom_line(data = p[[15]], aes(as.factor(sequencing_depth), observed_otu, group = 1), colour = "coral4",
            position = position_dodge(width=0.75)) +
  # fake axes (x-axis stops at year 2009, y-axis stops at value 45)
  annotate("segment", x = 0, y = 0, xend = 15, yend = 0, colour = "white") +
  # annotate("segment", x = 1991, y = 5, xend = 1991, yend = 45) + 
  annotate("text", label = c(data.frame(t(data[1,]))[,1])[1], x = replace(min(which(is.na(p[[1 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[1 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[1]]$observed_otu), colour = "deepskyblue4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[2,]))[,1])[1], x = replace(min(which(is.na(p[[2 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[2 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[2]]$observed_otu), colour = "cadetblue4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[3,]))[,1])[1], x = replace(min(which(is.na(p[[3 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[3 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[3]]$observed_otu), colour = "chartreuse4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[4,]))[,1])[1], x = replace(min(which(is.na(p[[4 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[4 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[4]]$observed_otu), colour = "chocolate4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[5,]))[,1])[1], x = replace(min(which(is.na(p[[5 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[5 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[5]]$observed_otu), colour = "darkgoldenrod4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[6,]))[,1])[1], x = replace(min(which(is.na(p[[6 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[6 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[6]]$observed_otu), colour = "firebrick4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[7,]))[,1])[1], x = replace(min(which(is.na(p[[7 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[7 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[7]]$observed_otu), colour = "forestgreen", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[8,]))[,1])[1], x = replace(min(which(is.na(p[[8 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[8 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[8]]$observed_otu), colour = "deeppink3", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[9,]))[,1])[1], x = replace(min(which(is.na(p[[9 ]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[9 ]]$observed_otu)))+0.5), 10.5), y = my.max(p[[9]]$observed_otu), colour = "gray4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[10,]))[,1])[1], x = replace(min(which(is.na(p[[10]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[10]]$observed_otu)))+0.5), 10.5), y = my.max(p[[10]]$observed_otu), colour = "chocolate1", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[11,]))[,1])[1], x = replace(min(which(is.na(p[[11]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[11]]$observed_otu)))+0.5), 10.5), y = my.max(p[[11]]$observed_otu), colour = "indianred2", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[12,]))[,1])[1], x = replace(min(which(is.na(p[[12]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[12]]$observed_otu)))+0.5), 10.5), y = my.max(p[[12]]$observed_otu), colour = "goldenrod4", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[13,]))[,1])[1], x = replace(min(which(is.na(p[[13]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[13]]$observed_otu)))+0.5), 10.5), y = my.max(p[[13]]$observed_otu), colour = "blue3", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[14,]))[,1])[1], x = replace(min(which(is.na(p[[14]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[14]]$observed_otu)))+0.5), 10.5), y = my.max(p[[14]]$observed_otu), colour = "cornflowerblue", hjust = 0) +
  annotate("text", label = c(data.frame(t(data[15,]))[,1])[1], x = replace(min(which(is.na(p[[15]]$observed_otu)))+0.5, is.infinite(min(which(is.na(p[[15]]$observed_otu)))+0.5), 10.5), y = my.max(p[[15]]$observed_otu), colour = "coral4", hjust = 0) +
  
  scale_x_discrete("sequencing depth") +
  scale_y_continuous("observed_otu") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("alpha-rarefaction observed_otu")

dev.off()


