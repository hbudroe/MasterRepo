#####Short script to do basic analyses on SILVA database representation/phyla for Holly#####

##Libraries
library(readxl)

##Read in silva meiodb tax strings
path <- "~/Desktop/UGA/LabLocalFiles/CodingMe/MeioDB/meiodb_v5_updatedtaxstringsmeiofauna_09222025.xlsx"
read_xlsx(path) -> meiodb

#Create dataframe to store counts of each phyla &
meiodbtracker <- data.frame("Taxa" = c("Cnidaria", "Loricifera", "Kinorhyncha", 
                                          "Nematoda", "Gastrotricha", "Gnathostomulida", "Rotifera", 
                                          "Xenacoelomorpha", "Ostracod", "Harpacticoida", "Misophrioida", 
                                          "Platyhelminthes", "Polychaeta", "Oligochaeta", "Tardigrada",
                                          "Entoprocta", "Chordata", "Echinodermata", "Mollusca", 
                                       "Hemichordata", "Porifera", "Ctenophora", "Arthropoda_other", "Chaetognatha"))

#Count
targettaxa <- c("D_9__Cnidaria", "D_13__Loricifera", "D_13__Kinorhyncha", "D_13__Nematoda", 
                "D_13__Gastrotricha", "D_13__Gnathostomulida", "D_13__Rotifera", "D_13__Xenacoelomorpha", "D_14__Ostracoda", 
                "D_16__Harpacticoida", "D_16__Misophrioida", "D_13__Platyhelminthes", "D_14__Polychaeta", 
                "D_15__Oligochaeta", "D_13__Tardigrada", "D_12__Entoprocta", "D_11__Chordata", "D_12__Echinodermata",
                "D_12__Mollusca", "D_12__Hemichordata", "D_17__Porifera", "D_9__Ctenophora", "D_12__Arthropoda",
                "D_11__Chaetognatha") #Change to D_13_entoprocta once fixed
countvec <- c()
for(i in 1:length(targettaxa)){
  sum(grepl(targettaxa[i], as.matrix(meiodb), ignore.case = TRUE)) -> counter
  countvec <- c(countvec, counter)
}
meiodbtracker$Counts <- countvec

#Plot
ggplot(meiodbtracker, aes(x = Taxa, y = Counts))+
  geom_col()+
  theme(axis.text.x = element_text(angle = 45))+
  ggtitle("Unique Animal Sequences in SILVA")

#################################################

####Script to plot the %identity of all our sequences in ddt with silva/meiodb

#Read in the blastn alignment (from cluster blastn) of ddt rep-seqs against meiodb
path <- "~/Desktop/UGA/LabLocalFiles/CodingMe/ddt_qiime/ddt_blastn/ddt_meiodb_blastn_results.csv"
blastcolnames <- c("qseqid","sseqid", "pident", "length", 
                   "mismatch", "gapopen", "qstart", "qend", 
                   "sstart", "send", "evalue", "bitscore")
read_csv(path, col_names = blastcolnames) -> ddtblastres
#steal tax table with tax assignments & hashes
tax_table(ddt_phylo_norm_euks) %>% data.frame() %>% rownames_to_column(var = "qseqid") %>% inner_join(ddtblastres) -> ddt_blast_tax



###for loop to plot lots heheh
phynames = c("Nematoda", "Xenacoelomorpha", "Gastrotricha", "Kinorhyncha", "Rotifera", 
             "Arthropoda", "Platyhelminthes", "Annelida", "Entoprocta")
my_plots <- vector(mode = "list", length = length(phynames))
plot_list <- vector("list", length(phynames))
for (i in phynames){
  subset(ddt_blast_tax, Taxon14 == i) %>%
    ggplot(., aes(x = pident))+
    geom_density()+
    theme_minimal()+
    labs(x = "% Identity to MeioDB Sequences")+
    xlim(88,100)+
    ggtitle(i) -> p
  p -> my_plots[[i]]
}
ggarrange(plotlist = my_plots,  ncol = sqrt(length(phynames)), nrow = sqrt(length(phynames)))

###also plot non meiodb corrected phyla (so far)
phynames = c("Mollusca", "Hemichordata", "Arthropoda", "Echinodermata", "Nemertea")
my_plots2 <- vector(mode = "list", length = length(phynames))
for (i in phynames){
  subset(ddt_blast_tax, Taxon13 == i) %>%
    ggplot(., aes(x = pident))+
    geom_density()+
    theme_minimal()+
    labs(x = "% Identity to MeioDB Sequences")+
    xlim(88,100)+
    ggtitle(i) -> p
  print(i)
  p -> my_plots2[[i]]
}
ggarrange(plotlist = my_plots2,  ncol = 1, nrow = 5)

phynames = c("Cnidaria", "Ctenophora") 
my_plots3 <- vector(mode = "list", length = length(phynames))
for (i in phynames){
  subset(ddt_blast_tax, Taxon10 == i) %>%
    ggplot(., aes(x = pident))+
    geom_density()+
    theme_minimal()+
    labs(x = "% Identity to MeioDB Sequences")+
    xlim(88,100)+
    ggtitle(i) -> p
  print(i)
  p -> my_plots3[[i]]
}
ggarrange(plotlist = my_plots3,  ncol = 1, nrow = length(phynames))


subset(ddt_blast_tax, Taxon12 == "Chordata") %>%
    ggplot(., aes(x = pident))+
    geom_density()+
    theme_minimal()+
    labs(x = "% Identity to MeioDB Sequences")+
  xlim(88,100)+
    ggtitle("Chordata") 








#########sCRAP code 
###try purrr to store plots in list
phynames = c("Nematoda", "Xenacoelomorpha", "Gastrotricha")
my_plots <- map(
  .x = phynames,
  .f = ~subset(ddt_blast_tax, Taxon14 == i) %>%
    ggplot(., aes(x = pident))+
    geom_density()+
    ggtitle(i)
)
# print the ggplots (they are stored in a list)
ggarrange(plotlist = my_plots,  ncol = 1, nrow = 3)

