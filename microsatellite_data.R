#' ---
#' title: "Example of polyploid, clonal, microsatellite data in R"
#' author: "Zhian N. Kamvar, Jonah C. Brooks, and Niklaus J. Grünwald"
#' output:
#'  html_document:
#'    keep_md: true
#'    toc: true
#' ---
#'
#' The purpose of this analysis is to investigate the population structure of
#' *Phytophthora infestans* from North America and South America using the data
#' from Goss *et al.*, 2014.
#'
#' Data Setup and Import
#' ---------------------
library('poppr')
library('purrr')
library('dplyr')
library('ggtree')
library('ggrepel')
library('ggplot2')
library('RColorBrewer')
data(Pinf)
Pinf <- Pinf %>%
  recode_polyploids(newploidy = TRUE)
pinfreps <- c(Pi02 = 2, D13 = 2, Pi33 = 6, Pi04 = 2, Pi4B = 2, Pi16 = 2,
              G11 = 2, Pi56 = 2, Pi63 = 3, Pi70 = 3, Pi89 = 2)
pinfreps <- fix_replen(Pinf, pinfreps)
Pinf
ContinentPAL <- setNames(c("firebrick", "blue"), popNames(Pinf))
setPop(Pinf) <- ~Country
CountryPAL   <- setNames(RColorBrewer::brewer.pal(4, "Dark2"), popNames(Pinf))
#'
#' One of the first things to do in an analysis of microsatellite data is to
#' ensure that the we have provide enough information to accurately call
#' multilocus genotypes.
#+ gac, results = "hide", fig.show = "hide"
genotype_curve(Pinf, sample = 1e4)
gcp <- last_plot()
#'
gcp +
  theme_bw() +
  geom_smooth(aes(group = 1)) +
  theme(text = element_text(size = 18)) +
  theme(panel.grid.major.x = element_blank()) +
  ggtitle(expression(paste(italic("P. infestans"), " genotype accumulation curve")))
ggsave(filename = "microsatellite_data_files/gac.svg", scale = 1.2)
#' The data that we have is mixed diploid and polyploid microsatellite markers.
#' This means that we should use Bruvo's distance, which accounts for ploidy.
#' Here, I'm using the genome addition model to calculate the distance. First,
#' I'm going to filter the multilocus genotypes.
pinf.bd     <- bruvo.dist(Pinf, replen = pinfreps, add = TRUE, loss = FALSE)
pinf.filter <- filter_stats(Pinf, dist = pinf.bd, plot = TRUE)
rug(pinf.bd, col = "#4D4D4D80")
#'
#+ echo = FALSE
svglite::svglite(file="microsatellite_data_files/filter.svg", width = 10, height = 8, pointsize = 20)
pinf.filter <- filter_stats(Pinf, dist = pinf.bd, plot = TRUE)
rug(pinf.bd, col = "#4D4D4D80")
dev.off()
#'
pinf.cutoff <- pinf.filter %>%
  transpose() %>%        # Transpose the data
  as_data_frame() %>%    # into a data frame and then
  select(THRESHOLDS) %>% # get the thresholds,
  flatten() %>%          # flatten to a list of nearest, farthest, and average and
  map_dbl(cutoff_predictor)  # calculate the cutoff for

mlg.filter(Pinf, dist = bruvo.dist, replen = pinfreps, loss = FALSE) <- pinf.cutoff["farthest"]
Pinf
#'
#'
#' Bootstrap analysis
#' ------------------
#+ pinfboot, results = "hide", cache = TRUE
pboot <- bruvo.boot(Pinf, replen = pinfreps, sample = 100, loss = FALSE,
                    showtree = FALSE, cutoff = 75)
#'
#+ warning = FALSE, fig.height = 9, fig.width = 10
ptree   <- apeBoot(pboot, pboot$node.label)
pstrata <- data_frame(taxa = indNames(Pinf)) %>%
  bind_cols(strata(Pinf)) %>%
  as.data.frame()
gt <- ggtree(ptree, layout = "circular") +
  geom_label_repel(aes(label = bootstrap, size = bootstrap),
                   nudge_x = -0.015, nudge_y = -0.005) +
  scale_size(range = c(2, 4))
gt <- gt %<+% pstrata +
  geom_tippoint(aes(color = Country), size = 3) +
  theme(legend.position = "right") +
  scale_color_manual(values = CountryPAL) +
  theme(text = element_text(size = 18))
gt
ggsave(gt, filename = "microsatellite_data_files/tree.svg")
#'
#' Minimum Spanning Network
#' ------------------------
#'
#' I'm using the filtered genotypes for this analysis.
#+ fminspan, fig.width = 9, fig.height = 10
fmin_span_net <- bruvo.msn(Pinf, replen = pinfreps, add = TRUE, loss = FALSE,
                          showplot = FALSE,
                          include.ties = TRUE,
                          threshold = pinf.cutoff["farthest"],
                          clustering.algorithm = "farthest_neighbor")
set.seed(69)
plot_poppr_msn(Pinf,
               fmin_span_net,
               inds = "none",
               mlg = FALSE,
               gadj = 2,
               nodebase = 1.15,
               palette = CountryPAL,
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               layfun = igraph::layout_with_kk)
#' Here's the original network
#+ minspan, fig.width = 9, fig.height = 10
mll(Pinf) <- "original"
min_span_net <- bruvo.msn(Pinf, replen = pinfreps, add = TRUE, loss = FALSE,
                          showplot = FALSE,
                          include.ties = TRUE)
set.seed(69)
plot_poppr_msn(Pinf,
               min_span_net,
               inds = "none",
               mlg = FALSE,
               gadj = 2,
               nodebase = 1.15,
               palette = CountryPAL,
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               layfun = igraph::layout_with_kk)
#'
#+ echo = FALSE, results = "hide", fig.width = 9, fig.height = 10
svglite::svglite(file="microsatellite_data_files/fminspan.svg", width = 9, height = 10)
set.seed(69)
mll(Pinf) <- "contracted"
plot_poppr_msn(Pinf,
               fmin_span_net,
               inds = "none",
               mlg = FALSE,
               gadj = 2,
               nodebase = 1.15,
               palette = CountryPAL,
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               layfun = igraph::layout_with_kk)
dev.off()
svglite::svglite(file="microsatellite_data_files/minspan.svg", width = 9, height = 10)
set.seed(69)
mll(Pinf) <- "original"
plot_poppr_msn(Pinf,
               min_span_net,
               inds = "none",
               mlg = FALSE,
               gadj = 2,
               nodebase = 1.15,
               palette = CountryPAL,
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               layfun = igraph::layout_with_kk)
dev.off()
#' Session Information
#' ===================
#'
if (!interactive()) options(width = 100)
devtools::session_info()
