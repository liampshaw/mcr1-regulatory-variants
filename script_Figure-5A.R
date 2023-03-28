# Shen analysis: construction of panel A of figure 4
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggsankey)
library(ggalluvial)
theme_basic <- function () { 
  theme_bw(base_size=12) %+replace% 
    theme(
      axis.text=element_text(colour="black")
    ) %+replace% 
    theme(
      panel.grid=element_blank()
    )
}

# Colours
variant.pal <- RColorBrewer::brewer.pal(name="Paired", n=8)
names(variant.pal) <- c("PV1", "PV2", "PV3", "PV4",
                        "SDV1", "SDV2", "SDV3", "SDV4")
variant.pal <- c(variant.pal, "#000000", "#bdbdbd")
names(variant.pal)[9:10] <- c("Consensus","Other")



#d <- read.csv('~/Dropbox/_Projects/mcr-1/Lancet_Microbe_MCRPEC/my_assemblies/shen_table.tsv', 
 #        header=T,
  #       stringsAsFactors = F, sep='\t')

d <- read.csv('data_metadata.csv',       header=T,
              stringsAsFactors = F, row.names = 1)

#d$REGULATORY.VARIANT[is.na(d$REGULATORY.VARIANT)] <- "Unknown" 

# Omit those without mcr-1...might want to check this 
d <- d[which(!is.na(d$MCR1.CONTIG)),]


# Make plasmid variants slightly easier to understand
d$PLASMIDS.ON.MCR1.CONTIG.simplified <- gsub("_.*", "", gsub("_.*,", ",", d$PLASMIDS.ON.MCR1.CONTIG))
d$PLASMIDS.ON.MCR1.CONTIG.simplified[which(d$PLASMIDS.ON.MCR1.CONTIG.simplified=="")] <- "None"

# Make even easier
d$PLASMIDS.ON.MCR1.CONTIG.simplified <- sapply(d$PLASMIDS.ON.MCR1.CONTIG.simplified, 
                                               function(x) ifelse(x=="None", "No replicon on contig", 
                                                                  ifelse(x=="IncX4", "IncX4", 
                                                                         ifelse(x=="IncI2", "IncI2", "Other replicon"))))
d$PLASMIDS.ON.MCR1.CONTIG.simplified <- ordered(d$PLASMIDS.ON.MCR1.CONTIG.simplified,
                                                  levels=c("No replicon on contig", 
                                                           "IncX4", "IncI2", "Other replicon"))


alluvial.df <- d %>% group_by(REGULATORY.VARIANT, PLASMIDS.ON.MCR1.CONTIG.simplified) %>%
  summarise(n=length(REGULATORY.VARIANT))
alluvial.df$PLASMIDS.ON.MCR1.CONTIG.simplified[which(alluvial.df$PLASMIDS.ON.MCR1.CONTIG.simplified=="")] <- "No replicon on contig"
alluvial.df$PLASMIDS.ON.MCR1.CONTIG.simplified[which(is.na(alluvial.df$PLASMIDS.ON.MCR1.CONTIG.simplified))] <- "No replicon on contig"

# Order plasmids by abundance
alluvial.df$PLASMIDS.ON.MCR1.CONTIG.simplified <- ordered(alluvial.df$PLASMIDS.ON.MCR1.CONTIG.simplified,
                                                          levels=c("No replicon on contig", 
                                                                   "IncX4", "IncI2", "Other replicon"),
                                                          labels=c("No replicon\npresent on contig", 
                                                                   "IncX4", "IncI2", "Other replicon(s)"))
alluvial.df$REGULATORY.VARIANT[alluvial.df$REGULATORY.VARIANT=="Consensus"] <- "Wild-type"
alluvial.df$REGULATORY.VARIANT <- ordered(alluvial.df$REGULATORY.VARIANT,
                                        levels=c("Wild-type", "PV1", "PV2", "PV3", "PV4", 
                                                 "SDV1", "SDV2", "SDV3", "SDV4", "Other"))
                                 
variant.pal["Wild-type"] <- "#000000"
p.alluvial <- ggplot(data = alluvial.df[which(alluvial.df$REGULATORY.VARIANT!="Other"),],
       aes(axis1 = REGULATORY.VARIANT, axis2 = PLASMIDS.ON.MCR1.CONTIG.simplified, y = n)) +
  geom_alluvium(aes(fill = REGULATORY.VARIANT),
                curve_type="quintic", alpha=1) +
  geom_stratum(fill="white", alpha=1, size=0.1, width = 0.2) +
  scale_x_discrete(limits = c("Promoter variant", "Plasmid"),
                   expand = c(0.15, 0.05)) +
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values=variant.pal)+
  ggrepel::geom_text_repel(stat = "stratum", aes(label=after_stat(stratum)), direction = "y", 
                           size = c(rep(4, 8), rep(4, 4)), segment.color = 'black', 
                           nudge_x = c(rep(-0.3, 8), rep(1.5, 4)))+
  xlim(c(0.5,2.5))

pdf("figure_5A.pdf", width=8, height=8)
p.alluvial
dev.off()


alluvial.df$REGULATORY.VARIANT[which(is.na(alluvial.df$REGULATORY.VARIANT))] <- "Other"
p.alluvial.all <- ggplot(data = alluvial.df,
                     aes(axis1 = REGULATORY.VARIANT, axis2 = PLASMIDS.ON.MCR1.CONTIG.simplified, y = n)) +
  geom_alluvium(aes(fill = REGULATORY.VARIANT),
                curve_type="quintic", alpha=1) +
  geom_stratum(fill="white", alpha=1, size=0.1, width = 0.2) +
  scale_x_discrete(limits = c("Promoter variant", "Plasmid"),
                   expand = c(0.15, 0.05)) +
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values=variant.pal)+
  ggrepel::geom_text_repel(stat = "stratum", aes(label=after_stat(stratum)), direction = "y", 
                           size = c(rep(4, 9), rep(4, 4)), segment.color = 'black', 
                           nudge_x = c(rep(-0.3, 9), rep(1.5, 4)))+
  xlim(c(0.5,2.5))
pdf("figure_supplementary_sankey-with-other.pdf", width=8, height=8)
p.alluvial.all
dev.off()
