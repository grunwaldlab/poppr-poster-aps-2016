#' ---
#' title: "Example of genomic data in R"
#' author: "Zhian N. Kamvar, Jonah C. Brooks, and Niklaus J. Grünwald"
#' output:
#'  html_document:
#'    keep_md: true
#'    toc: true
#' ---
#'
#'
#' I've taken the example data from https://github.com/knausb/vcfR_class/ and
#' ran it through
#'
#' 1. vcf_import.Rmd
#' 2. depth_filtering.Rmd
#'
#' Now I have the file called `TASSEL_GBS0077_dp_filtered.vcf.gz`. This will be
#' my starting point.
#'
#' Data Setup and Import
#' ---------------------
#+ setup, message = FALSE
library('ape')     # Creating a tree
library('vcfR')    # Reading in VCF file and conversion to genlight
library('poppr')   # multilocus genotype and linkage analysis
library('readr')   # reading in tsv file
library('tidyr')   # creating tidy data
library('purrr')   # manipulating lists
library('dplyr')   # manipulating data frames
library('knitr')   # printing tables
library('ggtree')  # plotting trees
library('ggplot2') # plotting cool graphics
library('ggrepel') # repelling bootstrap labels
rubfrag <- read.vcfR('TASSEL_GBS0077_dp_filtered.vcf.gz', verbose = FALSE)
(rfstrata <- read_tsv("rub_frag_strata.txt"))
(rf.sc <- rubfrag %>% vcfR2genlight() %>% as.snpclone())
#'
#' Since there are 40 entries in the strata and 41 in the VCF file, we're going
#' to only include the samples that are present in the strata.
rf.sc <- rf.sc[indNames(rf.sc) %in% rfstrata$VCF_ID, ]
strata(rf.sc) <- rfstrata[match(indNames(rf.sc), rfstrata$VCF_ID), ]
setPop(rf.sc) <- ~State
rf.sc
PAL <- setNames(RColorBrewer::brewer.pal(nPop(rf.sc), "Set2"), popNames(rf.sc))
#'
#' ### Filtering
#'
#' Because GBS data can't be clone-corrected by default due to various factors,
#' this step is important. I'm filtering multilocus genotypes by degree of
#' similarity. The `filter_stats()` function allows me to do this over all three
#' algorithms across all distance thresholds. I use the output of that function
#' to pass to the cutoff predictor which will find the largest gap in the data
#' and create a cutoff within that gap.
#'
rf.filter <- filter_stats(rf.sc, plot = TRUE)
rug(bitwise.dist(rf.sc, percent = TRUE), col = "#4D4D4D80")
#+ echo = FALSE
svglite::svglite(file="genomic_data_files//filter.svg", width = 10, height = 6)
rf.filter <- filter_stats(rf.sc, plot = TRUE)
rug(bitwise.dist(rf.sc, percent = TRUE), col = "#4D4D4D80")
dev.off()
#'
# predict cutoff for each algorithm
rf.cutoff <- rf.filter %>%
  transpose() %>%        # Transpose the data
  as_data_frame() %>%    # into a data frame and then
  select(THRESHOLDS) %>% # get the thresholds,
  flatten() %>%          # flatten to a list of nearest, farthest, and average and
  map_dbl(cutoff_predictor, 0.75)  # calculate the cutoff for
rf.cutoff
#'
#' Now that we have the cutoff set, I can filter the data.
mlg.filter(rf.sc) <- rf.cutoff["farthest"]
rf.sc
#'
#' To get a sense of the distribution of the MLGs, we should create a table
rf.tab <- mlg.table(rf.sc)
#'
#' To avoid issues with other analyses, we'll stick to the OR, WA, and CA
#' populations since they have more than 2 individuals.
rf.cow <- popsub(rf.sc, sublist = c("OR", "WA", "CA"))
#'
#' Poor Man's Jackknife
#' --------------------
#'
#' Bootstrapping 43K snps can be done with the `aboot()` function, but here, we
#' are using a jacknife approach with 20 samples of 500 SNPs at a time. This is
#' a similar process that happens internally with `aboot()`, but instead of
#' rebuilding a matrix of 43K SNPs after randomly sampling with replacement, we
#' are sampling 500 SNPs without replacement and calculating a tree off of
#' those.
#'
#' This will give us a measure of the internal consistency of the data.
#'
#+ tree, cache = TRUE
rf.tree <- phangorn::upgma(bitwise.dist(rf.sc))
set.seed(20160719)
sample_fun   <- function(i) seploc(rf.sc, block = 500, random = TRUE) %>% sample(10)
rf.sc.trees <- lapply(1:20, sample_fun) %>%
  flatten() %>%
  lapply(bitwise.dist) %>%
  lapply(phangorn::upgma)
nodelabs <- prop.clades(rf.tree, rf.sc.trees, rooted = TRUE)/2
nodelabs[nodelabs < 75] <- NA
rf.tree$node.label <- nodelabs
#'
#'
#' Now, I'm going to plot the tree using the *ggtree* package from Bioconductor
#' and the *ggrepel* package for the bootstrap labels.
rf.treeb <- apeBoot(rf.tree, nodelabs)
rf.strata <- strata(rf.sc) %>% rename(taxa = VCF_ID)
gt <- ggtree(rf.treeb) +
  geom_label_repel(aes(label = bootstrap, size = bootstrap),
                   nudge_x = -0.005, nudge_y = -0.005) +
  scale_size(range = c(2, 4))
gt <- gt %<+% rf.strata +
  geom_tippoint(aes(color = State), size = 3) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = PAL) +
  theme(text = element_text(size = 18))
gt
ggsave(gt, filename = "genomic_data_files/tree.svg")
#'
#' Diveristy analysis
#' ------------------
#'
#' We'll do an analysis of genotypic diversity for our populations using the
#' function `diversity_ci()`. It should be noted that the confidence intervals
#' are adjusted.
#+ mlg_diversity, cache = TRUE
rf.div <- diversity_ci(rf.cow, n = 1e4,
                       parallel = "multicore", ncpus = 4L, raw = FALSE)
kable(rf.div, digits = 2)
#'
#' We can also do a rarefaction of these stats.
#'
#+ mlg_rarefy, cache = TRUE
rf.rare <- diversity_ci(rf.cow, n = 1e4, rarefy = TRUE,
                        parallel = "multicore", ncpus = 4L, raw = FALSE)
kable(rf.rare, digits = 2)
#'
#'
#' Index of Association ($\bar{r}_d$)
#' ----------------------------------
#'
#' The index of association can tell us how clonal these populations appear. We
#' have the option of either doing a sliding window or a random sampling of
#' SNPs. Since the sliding window assumes contiguous chromosomes, and we have
#' `r format(length(chromosome(rf.cow)), big.mark = ",")` contigs, it would be
#' best to conduct a sampling analysis.
#'
#' Here I've chosen to do **1000** replicates with **500** SNPs for each
#' population and the total data set.
#+ ia_analysis, cache = TRUE, results = "hide"
set.seed(20160718)
rf.ia <- seppop(rf.cow) %>% # separate each population
  c(Total = rf.cow) %>%     # add the total population
  lapply(samp.ia, threads = 0, n.snp = 500L, reps = 1000L) %>%
  data.frame %>%    # convert list to data frame w/ 1000 rows
  gather(state, value) # convert to long, tidy data
#'
#'
#+ iasave
head(rf.ia)
rf.ia$state <- factor(rf.ia$state, levels = c("OR", "WA", "CA", "Total"))
ggia <- ggplot(rf.ia, aes(x = state, y = value)) +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()) +
  theme(text = element_text(size = 18)) +
  ggtitle(expression(paste(bar(r)[d], " per population sampled over 500 SNPs")))
ggia
ggsave(ggia, filename = "genomic_data_files/ia.svg")
#'
#' Minimum Spanning Network
#' ------------------------
#'
rf.cow_dist <- bitwise.dist(rf.cow, percent = TRUE, mat = FALSE,
                            missing_match = TRUE, differences_only = FALSE,
                            threads = 0)
min_span_net <- poppr.msn(rf.cow, rf.cow_dist, showplot = FALSE,
                          include.ties = TRUE,
                          threshold = rf.cutoff["farthest"],
                          clustering.algorithm = "farthest")

set.seed(70)
plot_poppr_msn(rf.cow,
               min_span_net,
               inds = "none",
               mlg = TRUE,
               gadj = 6,
               nodebase = 1.15,
               palette = PAL,
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               vertex.label.font = 2)
#'
#+ msn_save, echo = FALSE
svglite::svglite(file="genomic_data_files/fminspan.svg", width = 9, height = 10)
set.seed(70)
plot_poppr_msn(rf.cow,
               min_span_net,
               inds = "none",
               mlg = TRUE,
               gadj = 6,
               nodebase = 1.15,
               palette = PAL,
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               vertex.label.font = 2)
dev.off()
#'
#' Discriminant Analysis of Principle Components
#' ---------------------------------------------
#'
#+ dapc, cache = TRUE
rf.dapc <- dapc(rf.cow[order(pop(rf.cow))], n.pca = 12, n.da = 2)
#'
scatter.dapc(rf.dapc, col = PAL)
compoplot(rf.dapc, col = PAL)
#' Session Information
#' ===================
#'
if (!interactive()) options(width = 100)
devtools::session_info()
