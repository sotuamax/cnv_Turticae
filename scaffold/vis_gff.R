
library("ggplot2")
suppressMessages(library("argparse"))

parse_arg <- function() {
    # accept arguments from terminal
    parser <- ArgumentParser(description = "Accept count file, and do differential analysis given compare-pair. Usage: Rscript ")
    parser$add_argument("-gff", "--gff", help = "parsed gff file ")
    parser$add_argument("-O", "--output", help = "the output name")
    parser$parse_args()
}

argv <- parse_arg()
gff <- argv$gff
output <- argv$output

gff <- read.table(gff, sep = "\t", header = T)
chr <- c(unique(gff["chromosome"]))

xline <- NULL
for (c in chr$chromosome) {
  m <- max(gff[gff["chromosome"] == c,"gorder"])
  xline <- c(xline, m)
}

p <- ggplot(gff, aes(x = gorder, y = torder, col = contig)) + geom_point(size = 0.3, shape = 1) + theme_classic() + theme(legend.position = "none") + xlab("ref") + ylab("query") + ggtitle(output)

for (x in xline) {
  p <- p + geom_vline(xintercept = x, linetype = "dashed")
}

tig <- c(unique(gff["contig"]))
yline <- NULL
for (t in tig$contig) {
  if (nrow(gff[gff["contig"] == t, ]) > 10) {
    gmin <- min(gff[gff["contig"] == t, "torder"])
    gmax <- max(gff[gff["contig"] == t, "torder"])
    yline <- c(yline, gmin, gmax)
  }
  
}

for (y in yline) {
  p <- p + geom_hline(yintercept = y, linetype = "dashed", size = 0.3)
}

ggsave(plot = p, file = paste0(output,".pdf"), width = 5, height = 5)
