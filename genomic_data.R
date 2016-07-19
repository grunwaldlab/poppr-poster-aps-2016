#' Example of genomic data in R
#' ============================
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
#+ setup, cache = TRUE
library('vcfR')
library('poppr')
library('readr')
library('tidyr')
library('purrr')
library('ggplot2')
rubfrag <- read.vcfR('TASSEL_GBS0077_dp_filtered.vcf.gz', verbose = FALSE)
(rfstrata <- read_tsv("rub_frag_strata.txt"))
(rf.sc <- rubfrag %>% vcfR2genlight() %>% as.snpclone())
#'
#' Since there are 40 entries in the strata and 41 in the VCF file, we're going
#' to only include the samples that are present in the strata.
rf.sc <- rf.sc[indNames(rf.sc) %in% rfstrata$VCF_ID]
strata(rf.sc) <- rfstrata[match(indNames(rf.sc), rfstrata$VCF_ID), ]
setPop(rf.sc) <- ~State
rf.sc
#'
#' ### Filtering
#'
#' Because GBS data can't be clone-corrected by default due to various factors,
#' this step is important. I'm filtering multilocus genotypes by degree of
#' similarity. The `filter_stats()` function allows me to do this over all three
#' algorithms across all distance thresholds. I use the output of that function
#' to pass to the cutoff predictor which will find the largest gap in the data
#' and create a cutoff within that gap.
rf.filter <- filter_stats(rf.sc, plot = TRUE)
(rf.cutoff <- cutoff_predictor(rf.filter$farthest$THRESHOLDS))
#'
#' Now that we have the cutoff set, I can filter the data.
mlg.filter(rf.sc) <- rf.cutoff
rf.sc
#'
#' To get a sense of the distribution of the MLGs, we should create a table
rf.tab <- mlg.table(rf.sc)
#'
#' To avoid issues with other analyses, we'll stick to the OR, WA, and CA
#' populations since they have more than 2 individuals.
rf.cow <- popsub(rf.sc, sublist = c("OR", "WA", "CA"))
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
#+ ia_analysis, cache=TRUE, results = "hide"
set.seed(20160718)
rf.ia <- seppop(rf.cow) %>% # separate each population
	c(Total = rf.cow) %>%     # add the total population
	lapply(samp.ia, threads = 0, n.snp = 500L, reps = 1000L) %>%
	data.frame %>%    # convert list to data frame w/ 1000 rows
	gather(state, Ia) # convert to long, tidy data

ggplot(rf.ia, aes(x = state, y = Ia)) +
	geom_boxplot() +
	theme_bw() +
	ggtitle(expression(paste(bar(r)[d], " per population sampled over 500 SNPs"))) +
	ylab(expression(bar(r)[d]))

