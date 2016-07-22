# poppr-poster-aps-2016

source data for the APS 2016 poster on Tools for analysis of clonal population
genetic data in R

## Abstract

### [Tools for analysis of clonal population genetic data in R](http://www.apsnet.org/meetings/annual/abstracts/pages/abstractdetail.aspx?MID=816)

Knowledge of the population dynamics and evolution of plant pathogens allows
inferences on evolutionary processes involved in their adaptation to hosts,
pesticides, and other environmental pressures. These microbial pathogens
often require different tools for analysis due to the fact that they can be
clonal or partially clonal. With the advent of high-throughput sequencing
technologies, obtaining genome-wide data has become easier and cheaper than
ever before with techniques such as genotyping-by-sequencing. In 2013, we
created the widely used R package poppr for analysis of clonal populations.
In 2015, we published several additional extensions to poppr for use with
genomic data including the ability to define clones based on a genetic
distance threshold, minimum spanning networks with reticulation, and sliding
window analysis of the index of association. We present here an overview of
poppr as it pertains to traditional and high-throughput population genetic
data with select applications.

## Poster

The poster is located here: [poster/znk-aps-2016-poster.svg](poster/znk-aps-2016-poster.svg). This is an [Inkscape](https://inkscape.org/) SVG file.

The fonts utilized are: 

 - [Cooper Black](https://en.wikipedia.org/wiki/Cooper_Black)
 - [Nanum Gothic](https://en.wikipedia.org/wiki/Nanum_font)

The figures are derived from R scripts, detailed below. Minimal manipulation was
used on the figures in Inkscape (this includes rotation, resizing, and trimming).

## Scripts and Reports

The figures are generated from R scripts in this repository.

 - [genomic_data.R](genomic_data.R)
 - [microsatellite_data.R](microsatellite_data.R)

They are written in such a way that they generate reports, which you can find
here:

  - [genomic_data.md](genomic_data.md)
  - [microsatellite_data.md](microsatellite_data.md)

## Data

### Microsatellite data of *Phytophthora infestans*

These data are courtesy of Erica M. Goss from her 2014 paper:

Goss, Erica M., et al. "The Irish potato famine pathogen Phytophthora infestans
originated in central Mexico rather than the Andes." Proceedings of the National
Academy of Sciences 111.24 (2014): 8791-8796. [doi:
10.1073/pnas.1401884111](http://dx.doi.org/10.1073/pnas.1401884111)

It is now part of the *poppr* R package.

### GBS SNP data of *Phytophthora rubi*

These data are kindly provided by Javier F. Tabima and Niklaus J. Grünwald.
These data were first filtered with the package
[vcfR](https://github.com/knausb/vcfR#readme) before running the analysis. You
can find the scripts/tutorials for that at https://github.com/knausb/vcfR_class