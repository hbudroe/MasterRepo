#Author: Hannah Budroe
#Date: 10-16-24
#Purpose: This script is intended for figuring out some analyses/visualizations that 
#I find online and want to implement but have trouble running - to keep main scripts tidier


#Libraries
library(tidyverse)
library(vegan)
library(phyloseq)
#install.packages("devtools")
devtools::install_github("houyunhuang/ggcor")
library(ANCOMBC)
library(DT)
library(knitr)
library(kableExtra)

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


#mantel visualization ----
mantel02 <- fortify_mantel(varespec, varechem, 
                           spec.select = list(1:10, 5:14, 7:22, 9:32)) %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf), 
                 labels = c("<0.25", "0.25-0.5", ">=0.5"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">=0.05"),
                       right = FALSE))
quickcor(varechem, type = "upper") + geom_square() + 
  add_link(mantel02, mapping = aes(colour = p.value, size = r),
           diag.label = TRUE) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  add_diag_label() + remove_axis("x")
#> Warning: `add_diag_label()` is deprecated. Use `geom_diag_label()` instead.


###ANCOMBC2 ----
#(start with phyloseq object (taxa are rows = FALSE))

#follow tutotiral with their data
data(atlas1006, package = "microbiome")
tse = mia::convertFromPhyloseq(atlas1006)
# subset to baseline
tse = tse[, tse$time == 0]
# Re-code the bmi group
tse$bmi = recode(tse$bmi_group,
                 obese = "obese",
                 severeobese = "obese",
                 morbidobese = "obese")
# Subset to lean, overweight, and obese subjects
tse = tse[, tse$bmi %in% c("lean", "overweight", "obese")]
# Note that by default, levels of a categorical variable in R are sorted alphabetically. In this case, the reference level 
#for `bmi` will be `lean`. To manually change the reference level, for instance, setting `obese`as the reference level, use:
tse$bmi = factor(tse$bmi, levels = c("obese", "overweight", "lean"))
# You can verify the change by checking:
# levels(sample_data(tse)$bmi)
# Create the region variable
tse$region = recode(as.character(tse$nationality),
                    Scandinavia = "NE", UKIE = "NE", SouthEurope = "SE", 
                    CentralEurope = "CE", EasternEurope = "EE",
                    .missing = "unknown")
# Discard "EE" as it contains only 1 subject
# Discard subjects with missing values of region
tse = tse[, ! tse$region %in% c("EE", "unknown")]
print(tse)
out = ancombc(data = tse, assay_name = "counts", 
              tax_level = "Family", # phyloseq = NULL, 
              formula = "age + region + bmi", 
              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
              group = "bmi", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
              n_cl = 1, verbose = TRUE)
res = out$res
res_global = out$res_global
#primary result
tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "Age", "NE - CE", "SE - CE", 
             "US - CE", "Overweight - Obese", "Lean - Obese")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)
#visualize age
df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% #gets rid of taxon col for lfc and diff_abn
  mutate(taxon_id = res$diff_abn$taxon) %>% #taxon at end instead, renames taxon_id
  dplyr::select(taxon_id, everything()) #puts taxon back in front?
df_se = data.frame(res$se[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE") #renames Standard error columns to have "SE" on end

df_fig_age = df_lfc %>% 
  dplyr::left_join(df_se, by = "taxon_id") %>% #joins lfc with se dfs
  dplyr::transmute(taxon_id, age, ageSE) %>% #subsets just the age var (had bmi, region before)
  dplyr::filter(age != 0) %>% #subsets taxon with a lfc != 0 for age
  dplyr::arrange(desc(age)) %>% #arranges in decr lfc order for age
  dplyr::mutate(direct = ifelse(age > 0, "Positive LFC", "Negative LFC")) #adds qualifier col (pos/neg LFC)
df_fig_age$taxon_id = factor(df_fig_age$taxon_id, levels = df_fig_age$taxon_id) #make taxon a factor
df_fig_age$direct = factor(df_fig_age$direct, 
                           levels = c("Positive LFC", "Negative LFC")) #also makes 'direct' col (pos/neg LFC) factor
p_age = ggplot(data = df_fig_age, 
               aes(x = taxon_id, y = age, fill = direct, color = direct)) + #x factor(taxon), y LFC for age, color = pos/neg
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = age - ageSE, ymax = age + ageSE), width = 0.2,
                position = position_dodge(0.05), color = "black") + #error bars as ageSE
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of age") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))

#visualize bmi
df_fig_bmi = df_lfc %>% 
  filter(bmioverweight != 0 | bmilean != 0) %>% #filtering for any taxon with !=0 lfc for any bmi comparisons
  transmute(taxon_id, 
            `Overweight vs. Obese` = round(bmioverweight, 2),
            `Lean vs. Obese` = round(bmilean, 2)) %>% #subsetting only bmi cols, renaming col names, rounding to 2 dec
  pivot_longer(cols = `Overweight vs. Obese`:`Lean vs. Obese`, 
               names_to = "group", values_to = "value") %>% #for any taxa, separates separate comparisons of bmi into rows
  arrange(taxon_id)
lo = floor(min(df_fig_bmi$value)) 
up = ceiling(max(df_fig_bmi$value))
mid = (lo + up)/2 #rounding value(lfc) + averaging the up/down? 
p_bmi = df_fig_bmi %>%
  ggplot(aes(x = group, y = taxon_id, fill = value)) + #x comparison group, y taxon, fill = lfc value
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon_id, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

###try ancombc with own data, agglomerage at order level to start

mzg_raw %>% subset_taxa(Phylum!="Echinodermata") %>% subset_taxa(Phylum !="Bryozoa") %>% 
  subset_taxa(Phylum !="Sipuncula") %>% subset_taxa(Phylum !="Nemertea") %>%
  subset_taxa(Phylum !="Platyhelminthes") %>% subset_taxa(Phylum !="Rotifera") %>% 
  subset_taxa(Order!="Scleractinia") %>% filter_taxa(function(x) sum(x > 0) > 1, TRUE) -> mzg_filt_unnorm 
sample_data(mzg_filt_unnorm)$WaterType = factor(sample_data(mzg_filt_unnorm)$WaterType, levels = c("Subantarctic", "Subtropical"))
sample_data(mzg_filt_unnorm)$MinSize = factor(sample_data(mzg_filt_unnorm)$MinSize, levels = c("2", "5", "1"))

out <- ancombc(data = mzg_filt_unnorm, taxa_are_rows = FALSE, tax_level = "Order",
               formula = "WaterType + MinSize", group = "WaterType",
               p_adj_method = "fdr")
res = out$res

#make a table with results
knitr::kable(res$diff_abn) %>% kableExtra::kable_styling("striped") %>% 
  kableExtra::scroll_box(width = "100%")

#res_global = out$res_global
#primary result
tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "ST-SA", "05-02", "1-02")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#plot -- age plot
#visualize age
df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% #gets rid of taxon col for lfc and diff_abn
  mutate(taxon_id = res$diff_abn$taxon) %>% #taxon at end instead, renames taxon_id
  dplyr::select(taxon_id, everything()) #puts taxon back in front?
df_se = data.frame(res$se[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE") #renames Standard error columns to have "SE" on end

df_fig_age = df_lfc %>% 
  dplyr::left_join(df_se, by = "taxon_id") %>% #joins lfc with se dfs
  dplyr::transmute(taxon_id, WaterTypeSubtropical, WaterTypeSubtropicalSE) %>% #subsets just the age var (had bmi, region before)
  dplyr::filter(WaterTypeSubtropical != 0) %>% #subsets taxon with a lfc != 0 for age
  dplyr::arrange(desc(WaterTypeSubtropical)) %>% #arranges in decr lfc order for age
  dplyr::mutate(direct = ifelse(WaterTypeSubtropical > 0, "Positive LFC", "Negative LFC")) #adds qualifier col (pos/neg LFC)
df_fig_age$taxon_id = factor(df_fig_age$taxon_id, levels = df_fig_age$taxon_id) #make taxon a factor
df_fig_age$direct = factor(df_fig_age$direct, 
                           levels = c("Positive LFC", "Negative LFC")) #also makes 'direct' col (pos/neg LFC) factor
group_colors <- c("Negative LFC"= "#56B4E9", "Positive LFC" = "#E69F00")
ggplot(data = df_fig_age, 
               aes(x = taxon_id, y = WaterTypeSubtropical, fill = direct, color = direct)) + #x factor(taxon), y LFC for age, color = pos/neg
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = WaterTypeSubtropical - WaterTypeSubtropicalSE, ymax = WaterTypeSubtropical + WaterTypeSubtropicalSE), width = 0.2,
                position = position_dodge(0.05), color = "black") + #error bars as ageSE
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes in ST referenced to SA") + 
  scale_colour_manual(values=group_colors)+ #this adds colors according to group
  scale_fill_manual(values= group_colors)+
  #scale_fill_discrete(name = NULL) +
  #scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))+
  guides(fill = guide_legend(title = "LFC"), color = "none")

#visualize bmi
df_fig_bmi = df_lfc %>% 
  filter(WaterTypeSubtropical != 0) %>% #filtering for any taxon with !=0 lfc for any bmi comparisons
  transmute(taxon_id, 
            `ST vs. SA` = round(WaterTypeSubtropical, 2)) %>% #subsetting only bmi cols, renaming col names, rounding to 2 dec
  pivot_longer(cols = `ST vs. SA`, 
               names_to = "group", values_to = "value") %>% #for any taxa, separates separate comparisons of bmi into rows
  arrange(taxon_id)
lo = floor(min(df_fig_bmi$value)) 
up = ceiling(max(df_fig_bmi$value))
mid = (lo + up)/2 #rounding value(lfc) + averaging the up/down? 
df_fig_bmi %>%
  ggplot(aes(x = group, y = taxon_id, fill = value)) + #x comparison group, y taxon, fill = lfc value
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon_id, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to SA waters") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



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


