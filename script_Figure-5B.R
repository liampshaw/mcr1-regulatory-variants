# Shen analysis: construction of panel B of figure 4
library(ape)
library(circlize)
library(cowplot)
library(magick)
# Gu, Z. (2014) circlize implements and enhances circular visualization in R. Bioinformatics. DOI: 10.1093/bioinformatics/btu393

#aln <- read.dna('~/Dropbox/_Projects/mcr-1/IncX4/core.aln', format='fasta')

# Function for shannon
shannon <- function(a_vector){
  a_vector_norm <- a_vector/sum(a_vector)
  return (sum(-log(a_vector_norm)*a_vector_norm))
}
num_variants <- function(a_vector){
  return(sum(a_vector[2:length(a_vector)]))
}
# Number minority position
num_minority <- function(a_vector){
  return(sum(a_vector)-max(a_vector))
}

# Read in SNP positions 
snps <- read.csv('data_snippy-core-results.tab', sep='\t', header=T)

snps.df <- data.frame(pos=snps$POS,
                      entropy=apply(snps[,grepl("SRR", colnames(snps))], 
                                    MARGIN=1, FUN= function(x) shannon(table(as.character(x)))),
                      n.variants=apply(snps[,4:ncol(snps)], 
                                       MARGIN=1, FUN= function(x) num_variants(table(as.character(x)))),
                      n.minority=apply(snps[,4:ncol(snps)], 
                                       MARGIN=1, FUN= function(x) num_minority(table(as.character(x)))))


# gff3
gff <- read.csv('data_KU761327.1.gff3.tsv', sep='\t', header=T)
gff$gene <- rownames(gff)
gff$length <- gff$end-gff$start
gff$start.new <- sapply(1:nrow(gff), 
                        function(x) ifelse(gff[x,"strand"]=="+", gff[x,"start"], gff[x,"end"]))
gff$end.new <- sapply(1:nrow(gff), 
                      function(x) ifelse(gff[x,"strand"]=="+", gff[x,"end"], gff[x,"start"]))

gff$start.old <- gff$start
gff$end.old <- gff$end

# Sorting out gene group colours
gene.groups <- gff[-1,"group"]
#gene.groups.unique <- unique(gene.groups)
gene.group.colours <- c( "#fdae61", "grey50","#d7191c","#abd9e9", "#2c7bb6", "black", "yellow")
#gene.group.colours <- RColorBrewer::brewer.pal(n=length(gene.groups.unique), name="Set1")
#names(gene.group.colours) <- gene.groups.unique
names(gene.group.colours) <- c("toxin-antitoxin", "hypothetical", "mcr", 
                               "tra", "pil", "other", "transposase")

# Sorting out labels to print
gene.labels <- gff[-1,"product"] 
gene.labels[gene.labels=="hypothetical protein"] <- ""
gene.labels[gene.labels=="putative lipoprotein"] <- ""

#gene.labels[4] <- "PAP2-like ORF"
#gene.labels[8] <- "ParA"

# Plot v3
pdf("figure_5B.pdf", width=8, height=8)
circos.par(start.degree = 90, cell.padding =c(0.0,0, 0.0, 0))
circos.genomicInitialize(gff, major.by=2000, labels.cex = 1, axis.labels.cex = 1, sector.names = "")
circos.genomicTrack(gff[-1,], stack = TRUE, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicText(region, value, labels=gene.labels, y=0.2, cex = 0.7,niceFacing = TRUE)
                    }, track.height=0.05, bg.border=NA)
circos.track(ylim=c(1,2), bg.border=NA)
for (i in seq(2, nrow(gff))){
  if (gff[i,"strand"]=="+"){
    circos.arrow(gff[i,"start"], gff[i,"end"], y=1.55, 
                 width = 0.25, arrow.head.width = 0.5, 
                 arrow.head.length = 200, 
                 col = gene.group.colours[gff[i,"group"]])
  }
  if (gff[i,"strand"]=="-"){
    circos.arrow(gff[i,"start"], gff[i,"end"], y=1.75, 
                 width = 0.25, arrow.head.width = 0.5, 
                 arrow.head.length = 200,
                 arrow.position = "start", 
                 col=gene.group.colours[gff[i,"group"]])
  }
}
# Variant locations
circos.track(ylim=c(0,30), track.height=0.2 )
circos.yaxis(labels.cex = 0.5, tick = TRUE)
#circos.axis(h = "bottom",labels = FALSE)
for (x in 1:nrow(snps.df)){
  circos.trackPoints(sectors="chrKU761327.1", x=snps.df[x,"pos"], y=snps.df[x,"n.minority"], 
                     pch = 19,cex = 0.3)
}
# Masking
mask.df <- data.frame(chr=c("chrKU761327.1", "chrKU761327.1", "chrKU761327.1"), start=c(2337, 29959, 32039), end=c(4762, 30111, 32191), strand=c("+", "+", "+"))
circos.track(ylim=c(1,1.01), bg.border=NA, track.height=0.005)
circos.genomicTrack(mask.df, stack = TRUE, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, col="grey", border = NA, ...)
                    }, track.height=0.05, bg.border=NA)
dev.off()


# Look at the promoter variants in PV2
# All we need to do here is put the PV2 and PV3 which are 4807 and 4806 respectively 
# in a ggplot or something?
# Or maybe just in inkscape, although a good way of wasting time
table(as.character(snps[which(snps$POS==4807),grepl("SRR", colnames(snps))]))
table(as.character(snps[which(snps$POS==4806),grepl("SRR", colnames(snps))]))

