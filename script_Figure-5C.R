# Shen analysis: construction of panel C of figure 4
library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)
library(cowplot)
library(dplyr)
library(vegan)
library(reshape2)
library(ggtreeExtra)
library(ggnewscale)


aln <- read.dna('data_snippy-core-results.aln', format='fasta')
dist.aln <- dist.dna(aln, model = "N")
nj.tree <- nj(dist.aln)
nj.tree.midpoint <- midpoint(nj.tree)
nj.tree.root <- root(nj.tree, outgroup = "Reference") # We root on the reference plasmid

metadata <- read.csv('data_metadata.csv', header=T, stringsAsFactors = F, row.names = 1)

metadata$FILE <- gsub("\\.fa", "", metadata$FILE)

# Look at the reference mapping statistics
mapping.stats <- read.csv('data_snippy-core-results.txt', header=T, sep='\t')

# Only those with 0 unaligned
all_aligned <- mapping.stats$ID[which(mapping.stats$UNALIGNED==0)]

nj.tree.root.all.aligned <- keep.tip(nj.tree.root, all_aligned)

incx4 <- metadata$ISOLATE[which( metadata$PLASMIDS.ON.MCR1.CONTIG=="IncX4_1")]

metadata <- metadata[incx4,]
metadata <- metadata %>% add_row(FILE="Reference", REGULATORY.VARIANT="PV3")


# Keep only overlap
#nj.tree.IncX4 <- keep.tip(new.tree, metadata$FILE)

# Use different tree
new.tree = read.tree('data_snippy-core-results.refined.tre')
# Refined from original FastTree with https://github.com/liampshaw/beta-lactamases/blob/main/scripts/pipeline/scripts/refine_tree.py). 

metadata <- metadata[which(metadata$FILE %in% new.tree$tip.label),]
new.tree = keep.tip(new.tree, metadata$FILE)
metadata <- metadata[which(metadata$FILE %in% new.tree$tip.label),]


# Add mlst
mlst <- read.csv('data_shen-mlst.tsv', sep='\t', header=F, stringsAsFactors = F)
rownames(mlst) <- gsub("\\.fa", "", mlst$V1)
mlst$ST <- mlst$V3
metadata$ST <- mlst[as.character(metadata$FILE), "ST"]
# Add shape
metadata$type <- "isolate"

# Add row for reference sequence

# Metadata aligned only
metadata.aligned <- metadata[which(metadata$FILE %in% all_aligned),]
metadata.aligned$id <- metadata.aligned$FILE
rownames(metadata.aligned) <- metadata.aligned$id
#metadata.aligned <- metadata.aligned[nj.tree.root.all.aligned$tip.label,]
# Prune
metadata.aligned = metadata.aligned[which(rownames(metadata.aligned) %in% new.tree$tip.label),]

# New tree
new.tree = keep.tip(new.tree, metadata.aligned$id )
# Order 
metadata.aligned = metadata.aligned[new.tree$tip.label,]
# Set promoter variant as "PV3" - which it is, rather than 

metadata.aligned <- metadata.aligned[,c("id", colnames(metadata.aligned)[1:ncol(metadata.aligned)-1])]

# Make promoter variant nicer
metadata.aligned$REGULATORY.VARIANT <- sapply(metadata.aligned$REGULATORY.VARIANT,
                                            function(x) ifelse(x=="consensus", "Wild-type",
                                                               ifelse(x=="Other", "Other", 
                                                                      ifelse(x=="PV2", "PV2",
                                                                             ifelse(x=="PV3", "PV3", 
                                                                                    ifelse(x=="REFERENCE", "PV3", x))))))

metadata.aligned$REGULATORY.VARIANT <- ordered(metadata.aligned$REGULATORY.VARIANT, 
                                             levels=c("Consensus", "PV2", "PV3", "Other", "unknown"))


# Make a matrix of dummy variables of ST to add as heatmap
si <- metadata.aligned[,c( "id", "ST")]
si$id = NULL
si$ID = rownames(si)
si$Pos = 1




# Try
metadata.aligned$Variant = metadata.aligned$REGULATORY.VARIANT

# We root the tree to the only isolate with the consensus regulatory variant
# SRR15731853
new.tree = root(new.tree, outgroup = "SRR15731853")

  
p.0 = ggtree(new.tree, layout = 'rectangular',size=0.75)+
  geom_tippoint(fill="white", colour="black", shape=21)+
  geom_treescale(width=1)#

p.0 =  ggtree(new.tree, layout = 'rectangular')+
  geom_treescale(width=1)
p.1 = p.0 %<+% metadata.aligned[,c("id", "Variant")] +   
  geom_tippoint(aes(fill=Variant), colour="black", shape=21)+
  geom_treescale(width=1)+
  scale_fill_manual(values=c("black", "#1F78B4", "#B2DF8A", "#5E5E5E"),
                    guide=guide_legend(keywidth=1, keyheight=1, order=3))


metadata.aligned$Variant2 = metadata.aligned$Variant
p.2 = p.1 +
  new_scale_fill()+
  geom_fruit(
  data=metadata.aligned[,c("id", "Variant2")],
  geom=geom_tile,
  mapping=aes(y=id,fill=Variant2),
  offset=0.12,   # The distance between external layers, default is 0.03 times of x range of tree.
  pwidth=0.5 # width of the external layer, default is 0.2 times of x range of tree.
)+
  scale_fill_manual(values=c("black", "#1F78B4", "#B2DF8A", "#5E5E5E"),
   guide=guide_legend(keywidth=2, keyheight=2, order=3))+
  labs(fill="Variant")

# ST.colours
ST.counts = sort(table(metadata.aligned$ST))
abundant.ST = names(ST.counts[ST.counts>3]) 
abundant.ST.colours = rainbow(length(abundant.ST))
abundant.ST.colours = c(abundant.ST.colours, "Grey")
names(abundant.ST.colours) = c(abundant.ST, "ST with <4 isolates")

metadata.aligned$ST.simple = sapply(metadata.aligned$ST,
                                    function(x) ifelse(x%in% abundant.ST, x, "ST with <4 isolates"))

p.3 = p.2 + 
  new_scale_fill() +
  geom_fruit(
  data=metadata.aligned[,c("id", "ST.simple")],
  geom=geom_tile,
  mapping=aes(y=id,fill=ST.simple),
  offset=0.22,   # The distance between external layers, default is 0.03 times of x range of tree.
  pwidth=0.25 # width of the external layer, default is 0.2 times of x range of tree.
)+
  scale_fill_manual(values=abundant.ST.colours,
                    guide=guide_legend(keywidth=1, keyheight=1, order=3))+
  labs(fill="ST")
  

pdf("figure_5C.pdf", width=6, height=13)
p.3
dev.off()

metadata.aligned$fake.year = metadata.aligned$year-2016

# Just those in 2016?
isolates.2016 = metadata.aligned$id[which(metadata.aligned$year=="2016" | metadata.aligned$id=="Reference")]

metadata.aligned.2016 = metadata.aligned[isolates.2016,]
new.tree.2016 = keep.tip(new.tree, isolates.2016)
p.0.2016 =  ggtree(new.tree.2016, layout = 'rectangular')+
  geom_treescale(width=1)
p.1.2016 = p.0.2016 %<+% metadata.aligned.2016[,c("id", "Variant")] +   
  geom_tippoint(aes(fill=Variant), colour="black", shape=21)+
  geom_treescale(width=1)+
  scale_fill_manual(values=c("black", "#1F78B4", "#B2DF8A", "#5E5E5E"),
                    guide=guide_legend(keywidth=1, keyheight=1, order=3))


p.2.2016 = p.1.2016 +
  new_scale_fill()+
  geom_fruit(
    data=metadata.aligned.2016[,c("id", "Variant2")],
    geom=geom_tile,
    mapping=aes(y=id,fill=Variant2),
    offset=0.12,   # The distance between external layers, default is 0.03 times of x range of tree.
    pwidth=0.5 # width of the external layer, default is 0.2 times of x range of tree.
  )+
  scale_fill_manual(values=c("black", "#1F78B4", "#B2DF8A", "#5E5E5E"),
                    guide=guide_legend(keywidth=2, keyheight=2, order=3))+
  labs(fill="Variant")

p.3.2016 = p.2.2016 + 
  new_scale_fill() +
  geom_fruit(
    data=metadata.aligned.2016[,c("id", "ST.simple")],
    geom=geom_tile,
    mapping=aes(y=id,fill=ST.simple),
    offset=0.22,   # The distance between external layers, default is 0.03 times of x range of tree.
    pwidth=0.25 # width of the external layer, default is 0.2 times of x range of tree.
  )+
  scale_fill_manual(values=abundant.ST.colours,
                    guide=guide_legend(keywidth=1, keyheight=1, order=3))+
  labs(fill="ST")


pdf("figure_supplementary_phylogeny-emergence-2016.pdf", width=6, height=7)
p.3.2016
dev.off()
