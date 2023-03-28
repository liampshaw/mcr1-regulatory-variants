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

# Need to also do the ISApl1 status
d$isapl1.downstream.length <- d$ISAPL1.DOWNSTREAM.REL.END-d$ISAPL1.DOWNSTREAM.REL.START
d$isapl1.upstream.length <- d$ISAPL1.UPSTREAM.REL.END-d$ISAPL1.UPSTREAM.REL.START

d$isapl1 <- ifelse(d$isapl1.downstream.length>1000 & d$isapl1.upstream.length>1000, "both", 
                   ifelse(d$isapl1.upstream.length>1000 & d$isapl1.downstream.length<1000, "upstream", ifelse(d$isapl1.upstream.length<1000 & d$isapl1.downstream.length>1000, "downstream", "neither")))
d$mcr1.contig.length <-as.numeric(gsub("_.*", "", gsub(".*length_", "", d$MCR1.CONTIG)))

# Test the hypothesis that extremely short contigs <3kb are strongly associated
# with the wildtype regulatory sequence
test.table <- table(d[which(d$REGULATORY.VARIANT!="Other"),]$mcr1.contig.length<3000, d[which(d$REGULATORY.VARIANT!="Other"),]$REGULATORY.VARIANT=="Consensus")
chisq.test(test.table)


# Do the year analysis
#metadata <- read.csv('lancet-metadata.csv', header=T, stringsAsFactors = F)
#rownames(metadata) <- metadata$Sample
#d$year <- metadata[rownames(d), "Epoch"]
#d$source <- metadata[rownames(d), "Source"]

# categories of IncX4/PV3 on same contig
d$PV3.and.IncX4 <- ifelse(d$REGULATORY.VARIANT=="PV3" & d$PLASMIDS.ON.MCR1.CONTIG=="IncX4_1", "IncX4 with PV3", 
                          ifelse(d$REGULATORY.VARIANT=="PV3" & d$PLASMIDS.ON.MCR1.CONTIG!="IncX4_1", "PV3", 
                                 ifelse(d$REGULATORY.VARIANT!="PV3" & d$PLASMIDS.ON.MCR1.CONTIG=="IncX4_1", "IncX4", "neither")))
# plot the results
pv3.incx4.df <- d %>% group_by(year, PV3.and.IncX4)%>%
  summarise(n=length(FILE))%>%
  mutate(prop=n/sum(n))
pv3.incx4.df$PV3.and.IncX4 <- ordered(pv3.incx4.df$PV3.and.IncX4,
                                      levels=c("IncX4 with PV3", "IncX4", "PV3", "neither"))

p.all.isolates <- ggplot(pv3.incx4.df, aes(year, prop, group=PV3.and.IncX4, fill=PV3.and.IncX4)) +
  geom_bar(stat="identity", colour="black")+
  theme_bw()+
  theme_basic()+
  scale_fill_manual(values=c("#BBDD93","#a1d99b", "#e5f5e0", "white"))+
  ylab("Proportion of isolates")+
  labs(fill="mcr-1 contig\ncarries")+
  xlab("")


# Check for just pigs
# plot the results
pigs.pv3.incx4.df <- d[which(d$source=="Pig"),] %>% group_by(year, PV3.and.IncX4)%>%
  summarise(n=length(FILE))%>%
  mutate(prop=n/sum(n))
pigs.pv3.incx4.df$PV3.and.IncX4 <- ordered(pigs.pv3.incx4.df$PV3.and.IncX4,
                                           levels=c("IncX4 with PV3", "IncX4", "PV3", "neither"))

p.pigs <- ggplot(pigs.pv3.incx4.df, aes(year, prop, group=PV3.and.IncX4, fill=PV3.and.IncX4)) +
  geom_bar(stat="identity", colour="black")+
  theme_bw()+
  theme_basic()+
  scale_fill_manual(values=c("#BBDD93","#a1d99b", "#e5f5e0", "white"), breaks=levels(pigs.pv3.incx4.df$PV3.and.IncX4))+
  ylab("Proportion of isolates")+
  labs(fill="mcr-1 contig\ncarries")+
  xlab("")+
  theme(legend.position="none")

pdf("figure_supplementary_increase-IncX4.pdf", width=8, height=4)
cowplot::plot_grid(p.all.isolates+ggtitle("A. All isolates"), 
                   p.pigs+ggtitle("B. Pig isolates"), 
                   rel_widths = c(1.2,0.8))
dev.off()
