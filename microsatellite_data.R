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
library('igraph')
library('ggplot2')
library('RColorBrewer')
data(Pinf)
Pinf <- Pinf %>%
  recode_polyploids(newploidy = TRUE)
pinfreps <- c(Pi02 = 2, D13 = 2, Pi33 = 6, Pi04 = 2, Pi4B = 2, Pi16 = 2,
              G11 = 2, Pi56 = 2, Pi63 = 3, Pi70 = 3, Pi89 = 2)
pinfreps <- fix_replen(Pinf, pinfreps)
Pinf
#'
#'
#' The
pinf.bd <- bruvo.dist(Pinf, replen = pinfreps, add = TRUE, loss = FALSE)
pinf.filter <- filter_stats(Pinf, dist = pinf.bd, plot = TRUE)
rug(pinf.bd, col = "#4D4D4D80")

pinf.cutoff <- pinf.filter %>%
  transpose() %>%        # Transpose the data
  as_data_frame() %>%    # into a data frame and then
  select(THRESHOLDS) %>% # get the thresholds,
  flatten() %>%          # flatten to a list of nearest, farthest, and average and
  map_dbl(cutoff_predictor)  # calculate the cutoff for

setPop(Pinf) <- ~Country
mlg.filter(Pinf, dist = bruvo.dist, replen = pinfreps, loss = FALSE) <- pinf.cutoff["farthest"]
min_span_net <- bruvo.msn(Pinf, replen = pinfreps, add = TRUE, loss = FALSE,
                          showplot = FALSE,
                          include.ties = TRUE,
                          threshold = pinf.cutoff["farthest"],
                          clustering.algorithm = "farthest_neighbor")
set.seed(69)
plot_poppr_msn(Pinf,
               min_span_net,
               inds = "none",
               mlg = FALSE,
               gadj = 2,
               nodebase = 1.15,
               palette = RColorBrewer::brewer.pal(4, "Set1"),
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               layfun = layout_with_kk)
#' Session Information
#' ===================
#'
if (!interactive()) options(width = 100)
devtools::session_info()
