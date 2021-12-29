# Updated on Oct-04-2021
# Developed by Meiyuan Ji (for questions contact: meiyuan.ji@utah.edu)
# This program estimates the parameters (mean and standard devidation) for the distribution of single-copy read depth via normal mixture models with one component. 

suppressMessages(library("argparse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stats4"))
suppressMessages(library("mclust"))
# suppressMessages(library("bbmle"))

parse_arg <- function() {
    # accept arguments from terminal
    parser <- ArgumentParser(description = "Use the median read-depth value as representative observation for each gene region. Observed values was fitted using normal mixture model given one component. ")
    parser$add_argument("-RD", "--read_depth", help = "mean and median of read-depth of gene region, collected with stepsize n bp (n=100). ")
    parser$add_argument("-left_disperse", "--left_disperse", default = 0.45, help = "disperse range of read-depth on the left side of the peak-value. ")
    parser$add_argument("-right_disperse", "--right_disperse", default = 0.65, help = "disperse range of read-depth on the right side of the peak-value. ")
    parser$add_argument("-O", "--outdir", help = "output directory")
    parser$parse_args()
}

argv <- parse_arg()
rd <- argv$read_depth
outdir <- argv$outdir
ldis <- argv$left_disperse
rdis <- argv$right_disperse

rd_df <- read.table(rd, sep = "\t", header = T, row.names = 1)
med <- mean(rd_df$depth_median)
rd_df_filter <- rd_df[(rd_df$size > 2) & (rd_df$depth_median > med*0.1), ]

max_median <- ceiling(max(rd_df_filter$depth_median))
min_median <- floor(min(rd_df_filter$depth_median))

pdf(file.path(outdir, 'median_hist.pdf'), width = 4, height = 4)
rd_hist_freq <- hist(rd_df_filter$depth_median, breaks = (max_median-min_median)%/%10, freq = T) # %/% integer divistion; %% get modulo

hist_breaks <- rd_hist_freq$breaks 
hist_count <- rd_hist_freq$counts
hist_mids <- rd_hist_freq$mids

hist_freq <- data.frame(bar_mid = hist_mids, bar_count = hist_count)
total_count <- sum(hist_freq$bar_count)
hist_freq$bar_prop <- hist_freq$bar_count/total_count 
# observed number of frequency and its proportion
mode_mid <- hist_freq[hist_freq$bar_count == max(hist_freq$bar_count), ]$bar_mid # **** important (find mode)
abline(v = mode_mid, col = "red")

# normally, for inbreds, model fits on a broader range of read-depth (say, upper or lower boundary of 0.75 in respect to the mean); for non-inbreds, models fits on a narrower range of read-depth (say, upper or lower boundary of 0.35 in respect to the mean)
# a more stringent data set for modele fitting needs to cut more data and left data most around the mean for mean/sd optimization with normal distribution

if (!is.na(ldis) & !(is.na(rdis))) {
  mode_upper <- mode_mid + mode_mid*as.numeric(rdis)
  mode_lower <- mode_mid - mode_mid*as.numeric(ldis)
}

abline(v = mode_upper, col = "red")
abline(v = mode_lower, col = "red")
dev.off()

rd_df_around_mode <- rd_df_filter[rd_df_filter$depth_median > mode_lower & rd_df_filter$depth_median < mode_upper, ]

#sd <- sd(rd_df_around_mode$depth_median)
#data_init <- rd_df_around_mode$depth_median
#neg_log_gaussian <- function(p) {
#  r <- -sum(dnorm(data_init, mean = p[1], sd = p[2], log = T))
#  if (!is.finite(r)) {
#    r <- 1e+20
#  }
#  return(r)
#}

#gaussian_fit_init <- nlm(neg_log_gaussian, c(mode_mid, sd_median), hessian = T)

#gaussian_optimization_iteration <- function(df, mu, sd, times) {
#  upper <- mu + times*sd
#  lower <- mu - times*sd
#  new_df <- df[df$depth_median <= upper & df$depth_median >= lower,]
#  data <- new_df$depth_median
#  neg_log_gaussian <- function(p) {
#    r <- -sum(dnorm(data, mean = p[1], sd = p[2], log = T))
#    if (!is.finite(r)) {
#      r <- 1e+20
#    }
#    return(r)
#  }
#  gaussian_fit <- nlm(neg_log_gaussian, c(mu, sd), hessian = T)
#  return(gaussian_fit)
#}

# negative binomial distribution 
# neg_log_NB <- function(size, p) {
#   r <- -sum(dnbinom(as.integer(rd_df_filter$depth_median), size = size, prob = p, log = T))
#   if (!is.finite(r)) {
#     r <- 1e+20}
#   return(r)
# }

# NB_fit <- mle2(minuslogl = neg_log_NB, 
#           start = list(size = 5, p = 0.05), method = "L-BFGS-B")

# size <- NB_fit@coef[1]
# p <- NB_fit@coef[2]

#mu_init <- gaussian_fit_init$estimate[1]
#sd_init <- gaussian_fit_init$estimate[2]
#for (i in seq(iter)) {
  #print(paste0(mu_init, "-", sd_init))
#  gaussian_fit <- gaussian_optimization_iteration(rd_df, mu_init, sd_init, 2)
#  mu = gaussian_fit$estimate[1]
#  sd = gaussian_fit$estimate[2]
#  mu_init <- mu
#  sd_init <- sd
#}

# model-based clustering, classification, and density estimation based on finite normal mixture modelling. 
# parameter estimated via the EM algorithm for normal mixture models with a variety of covariance structures.

clust_M <- Mclust(rd_df_around_mode$depth_median, G = 1) 
mu <- clust_M$parameter$mean
sd <- clust_M$parameter$variance$sigmasq^0.5

hist_freq$exp_prop <- dnorm(hist_freq$bar_mid, mean = mu, sd = sd)
p <- ggplot(data = hist_freq) + geom_col(aes(x = bar_mid, y = bar_prop)) + geom_line(aes(x = bar_mid, y = exp_prop), col = "red") + geom_vline(xintercept = mu + 3*sd, col = "red") + geom_vline(xintercept = mu - 3*sd, col = "red") + xlim(c(mu-10*sd, mu+10*sd)) + theme_bw()
ggsave(filename = file.path(outdir, "depth_fit.pdf"), plot = p, width = 4, height = 4)

rd_model <- data.frame(mean = mu, sd = sd)
write.table(rd_model, file = file.path(outdir, "depth_fit.txt"), sep = "\t", row.names = F, quote = F)

