# Analysis of IncX4 plasmid variation
library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)
library(cowplot)
library(dplyr)
library(vegan)
library(reshape2)


aln <- read.dna('data_snippy-core-results.aln', format='fasta')
dist.aln <- dist.dna(aln, model = "N")
nj.tree <- nj(dist.aln)
nj.tree.midpoint <- midpoint(nj.tree)
nj.tree.root <- root(nj.tree, outgroup = "Reference") # We root on the reference plasmid

metadata <- read.csv('data_metadata.csv', header=T, stringsAsFactors = F, row.names =1)

metadata$FILE <- gsub("\\.fa", "", metadata$FILE)

# Look at the reference mapping statistics
mapping.stats <- read.csv('data_snippy-core-results.txt', header=T, sep='\t')

# Only those with 0 unaligned
all_aligned <- mapping.stats$ID[which(mapping.stats$UNALIGNED==0)]

nj.tree.root.all.aligned <- keep.tip(nj.tree.root, all_aligned)

incx4 <- metadata$ISOLATE[which(metadata$PLASMIDS.ON.MCR1.CONTIG=="IncX4_1")]

metadata <- metadata[incx4,]

# Keep only overlap
metadata <- metadata[which(metadata$FILE %in% nj.tree.root.all.aligned$tip.label),]
nj.tree.IncX4 <- keep.tip(nj.tree.root.all.aligned, metadata$FILE)


# Add mlst
mlst <- read.csv('data_shen-mlst.tsv', sep='\t', header=F, stringsAsFactors = F)
rownames(mlst) <- gsub("\\.fa", "", mlst$V1)
mlst$ST <- mlst$V3
metadata$ST <- mlst[as.character(metadata$FILE), "ST"]
# Add shape
metadata$type <- "isolate"

# Add row for reference sequence
# Set promoter variant as "REFERENCE" even though it is actually PV3
metadata <- metadata %>% add_row(FILE="Reference", REGULATORY.VARIANT="REFERENCE", type="reference")

# Metadata aligned only
metadata.aligned <- metadata[which(metadata$FILE %in% all_aligned),]
metadata.aligned$id <- metadata.aligned$FILE
rownames(metadata.aligned) <- metadata.aligned$id
metadata.aligned <- metadata.aligned[nj.tree.root.all.aligned$tip.label,]
metadata.aligned <- metadata.aligned[,c("id", colnames(metadata.aligned)[1:ncol(metadata.aligned)-1])]

# Make promoter variant nicer
metadata.aligned$REGULATORY.VARIANT <- sapply(metadata.aligned$REGULATORY.VARIANT,
                                            function(x) ifelse(x=="consensus", "Consensus",
                                                               ifelse(x=="Other", "Other", 
                                                                      ifelse(x=="PV2", "PV2",
                                                                             ifelse(x=="PV3", "PV3", 
                                                                                    ifelse(x=="REFERENCE", "PV3", "unknown"))))))

metadata.aligned$REGULATORY.VARIANT <- ordered(metadata.aligned$REGULATORY.VARIANT, 
                                             levels=c("Consensus", "PV2", "PV3", "Other", "unknown"))

# Plot only those that are aligned. 
p1 <- ggtree(nj.tree.IncX4, layout = 'rectangular') %<+% metadata.aligned + 
  geom_tippoint(aes(fill=REGULATORY.VARIANT), size=1.5, shape=21, colour="black")+
  scale_fill_manual(values=c("#1F78B4", "#B2DF8A", "#bdbdbd", "#000000"))+
  geom_treescale(width = 1, x=-2)+
  labs(fill="Promoter\nvariant")+
  theme(legend.position = "left")+
  guides(fill = guide_legend(override.aes = list(size=10)))

pdf('figure_supplementary_phylogeny-emergence.pdf', width=6, height=9); p1; dev.off()

