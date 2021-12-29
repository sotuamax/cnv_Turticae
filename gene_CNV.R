# Update on Dec-28-2021
# Developed by Meiyuan Ji (for questions contact: meiyuan.ji@utah.edu)
# This program is based on the model parameters of the fitting distribution and estimate copy number for sites on gene region. 

suppressMessages(library("argparse"))
suppressMessages(library("stats4"))
suppressMessages(library("stringr"))
suppressMessages(library("dplyr"))

parse_arg <- function() {
    # accept arguments from terminal
    parser <- ArgumentParser(description = " ")
    parser$add_argument("-RD", "--RD", help = "Read-depth of site on gene region. ")
    parser$add_argument("-fit", "--fit", help = "parameters for the modeling distribution. ")
    parser$add_argument("-O", "--output", help = "the output folder")
    parser$parse_args()
}

argv <- parse_arg()
RD <- argv$RD
fit <- argv$fit
out <- argv$output

RD_df <- read.table(RD, header = T, sep = "\t")
fit_df <- read.table(fit, header = T, sep = "\t")
fit_mean <- fit_df[1, "mean"] # mean of single-copy coverage
fit_sd <- fit_df[1, "sd"] # standard deviation of single-copy coverage distribution
###
# RD_df$read_depth <- ifelse((RD_df$read_depth < 3) | (RD_df$read_depth <= fit_mean*0.1), 0, RD_df$read_depth) # filter noise mapping
RD_df$norm_depth <- round(RD_df$read_depth/fit_mean, 2)
RD_df$p <- dnorm(RD_df$read_depth, mean = fit_mean, sd = fit_sd)
RD_df$padj <- p.adjust(RD_df$p, method = "bonferroni")

# read_depth, norm_depth, p, padj

write.table(RD_df, file.path(out, "CNV.txt"), sep = "\t", quote = F, row.names = F)
