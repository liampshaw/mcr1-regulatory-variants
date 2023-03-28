# Combines the 3 panels for Figure 4
# N.B. the versions in the paper are edited and annotated before combining
# These make it look somewhat different - but the underlying data is the same. 
library(cowplot)
library(magick)

p.alluvial <- ggdraw() + draw_image(magick::image_read_pdf("figure_5A.pdf", density = 600))
p.circos <- ggdraw() + draw_image(magick::image_read_pdf("figure_5B.pdf", density = 600))
p.phylogeny <- ggdraw() + draw_image(magick::image_read_pdf("figure_5C.pdf", density = 600))


pdf('figure_5.pdf', width=14, height=4)
cowplot::plot_grid(p.alluvial, p.circos,
                   p.phylogeny, 
                   rel_widths = c(1,1, 1),
                   nrow=1)
dev.off()
 
