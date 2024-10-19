#Author: Hannah Budroe
#Date: 10-16-24
#Purpose: This script is intended for figuring out some analyses/visualizations that 
#I find online and want to implement but have trouble running - to keep main scripts tidier

#Libraries
library(tidyverse)
library(vegan)
library(phyloseq)

#Easy data to use - phyloseq GlobalPatterns
data("GlobalPatterns")

#Set path
setwd("~/Desktop/UGA/Practice_plots_for_now_delete_later")

### Cordier Sci paper Plankton/benthic patterns in sediment (1c) ----
### Figure 1C: Accumulation curves sampling effort (planktonic pico nano size fractions) 

#Notes: and follow along
#this code starts from ASV tables that are subset according to various classifications
#When I do the specaccum curves, they all show the same pattern that sharply increases right after 0 samples
#This pattern holds for all the curves added with GP dataset, number of ASVs on 
#y axis jumps up to 3ish (even if this was bc of the 10^3 label on the yaxis there's 70k asvs in gp_feces)?
#this is far too low.. but idk what is going on

gp_soil <- subset_samples(GlobalPatterns, SampleType == "Soil")
gp_feces <- subset_samples(GlobalPatterns, SampleType == "Feces")
gp_skin <- subset_samples(GlobalPatterns, SampleType == "Skin")

#do the specaccum curves individually for each ASV table, method = random
sp_accum_gp_soil <- specaccum(otu_table(gp_soil), method = "random")
sp_accum_gp_feces <- specaccum(otu_table(gp_feces), method = "random")
sp_accum_gp_skin <- specaccum(otu_table(gp_skin), method = "random")

#plot the accumulation curves
pdf("Figure1/Figure_1C_ASVs_accumulation_curves_sampling_effort.pdf", width=3.5, height=4) #default from code
plot(sp_accum_gp_soil, col = "cadetblue1", ylab = "#ASVs (10³)", 
     xlab = "#Samples", ci.type="poly", lwd=2, ylim= c(1,length(otu_table(gp_soil))+1000), 
     ci.lty=0, ci.col=alpha("cadetblue1", alpha = 0.2), cex.axis=0.6) #have to change data name + ylim
lines(sp_accum_gp_feces, col = "darkgoldenrod1", ci.type="poly", lwd=2, ci.lty=0, ci.col=alpha("darkgoldenrod1", alpha = 0.2))
lines(sp_accum_gp_skin, col = "blue4", ci.type="poly", lwd=2, ci.lty=0, ci.col=alpha("blue4", alpha = 0.2))
dev.off()

#try again with vegan example
data(BCI)
sp1 <- specaccum(BCI, "random") 
#This works like it does in the example --> species as cols, samples as rows I think
#but the GP asv table has taxa as rows instead - try flipping orientation t()
sp_accum_gp_soil <- specaccum(t(otu_table(gp_soil)), method = "random")
sp_accum_gp_feces <- specaccum(t(otu_table(gp_feces)), method = "random")
sp_accum_gp_skin <- specaccum(t(otu_table(gp_skin)), method = "random")

#plot the accumulation curves
pdf("Practice_plots_for_now_delete_later", width=3.5, height=4) #default from code
plot(sp_accum_gp_soil, col = "cadetblue1", ylab = "#ASVs (10³)", 
     xlab = "#Samples", ci.type="poly", lwd=2, ylim= c(1,length(otu_table(gp_soil))+1000), 
     ci.lty=0, ci.col=alpha("cadetblue1", alpha = 0.2), cex.axis=0.6) #have to change data name + ylim
lines(sp_accum_gp_feces, col = "darkgoldenrod1", ci.type="poly", lwd=2, ci.lty=0, ci.col=alpha("darkgoldenrod1", alpha = 0.2))
lines(sp_accum_gp_skin, col = "blue4", ci.type="poly", lwd=2, ci.lty=0, ci.col=alpha("blue4", alpha = 0.2))
dev.off()
#this worked! The reason the original asv axis only went up to 3ish was that's the nsamples
#able to plot successfully if t() -- make sure taxa are cols, samples as rows
#figure out the pdf() -- first command is file path
#but it saved the file as the folder (since wouldn't let me put /name), and opening is text string???



#Compositional tree plots Cordier 2a ----

#library('doMC')
#registerDoMC(cores = 16)
library(ggplot2)
library(ggpubr)
library(treemapify)
library(seqinr)

#notes


pdf("Figure2/Figure_2_Treemap_richness_eupthotic.pdf", width = 8, height = 2) 
treemap::treemap(subset(df_rich_agg, df_rich_agg$variable == "Euphotic"), 
                 title = "Richness - Euphotic", 
                 algorithm = "pivotSize", border.lwds = c(1,0.5,0.1), 
                 border.col = c("black", "black", "black"),
                 mapping = c(0,0,0),
                 index=c("rank3","rank4", "rank5"),
                 vSize="value",
                 vColor="rank_col",
                 type="index",
                 fontsize.labels=c(18,14,10),             # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
                 fontcolor.labels=c("Grey", "orange","white"),    # Color of labels
                 fontface.labels=c(3,2,1),                # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
                 bg.labels=c("transparent"),              # Background color of labels
                 align.labels=list(
                   c("left", "top"), 
                   c("left", "bottom"),
                   c("center", "center")
                 ),                                   # Where to place labels in the rectangle?
                 overlap.labels=0.8,                  # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
                 inflate.labels=F,
                 force.print.labels = F) 
dev.off()


