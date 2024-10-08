% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plot_microbiota.R
\name{plot_microbiota}
\alias{plot_microbiota}
\title{Clustered stackbar coupled with DESeq2}
\usage{
plot_microbiota(
  ps_object = ps,
  exp_group = "group",
  subset_group = NULL,
  sample_name = "SampleID",
  main_level = "Phylum",
  sub_level = "Family",
  threshold = 1,
  n_phy = 4,
  mean_group = FALSE,
  hues = c("Oranges", "Greens", "Blues", "Purples"),
  color_bias = 2,
  n_row = 1,
  n_col = NULL,
  text_size = 9,
  legend_size = 7,
  x_axis_size = 8,
  differential_analysis = FALSE,
  mult_comp = FALSE,
  selected_comparisons = NULL,
  test = c("Wald", "LRT")[1],
  fdr_threshold = 0.05,
  sig_lab = FALSE,
  fitType = c("parametric", "local", "mean", "glmGamPoi")[1],
  sfType = c("ratio", "poscounts", "iterate")[1],
  betaPrior = FALSE,
  reduced = FALSE,
  quiet = TRUE,
  minReplicatesForReplace = 7,
  modelMatrixType = c("standard", "expanded")[1],
  useT = FALSE,
  minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5,
  parallel = FALSE
)
}
\arguments{
\item{ps_object}{A \linkS4class{phyloseq} object. Must have a \link{sample_data}, \link{otu_table}, and \link{tax_table} components.}

\item{exp_group}{A column in the \link{sample_data} containing the experimental group information.}

\item{subset_group}{Default NULL. Groups among the 'exp_group' column from the \link{sample_data} object to subset.}

\item{sample_name}{Name of the column in \link{sample_data} containing the unique sample identifier.}

\item{main_level}{Default 'Phylum'. Level present in the \link{tax_table} to which taxa will be clustered.}

\item{sub_level}{Default 'Family'. Level present in the \link{tax_table} to which taxa will be plotted and analyzed.}

\item{threshold}{Default 1, \% threshold to regroup taxa with lower relative abundance into the 'other' groups}

\item{n_phy}{Default 4, number of main_level to plot. Same number of colors must be given in the 'hues' parameter.}

\item{mean_group}{Default FALSE. Whether or not agglomerate samples belonging to the same exp_group.}

\item{hues}{Color used to represent main_level, should be the same number than n_phy parameter. See \link{colorRampPalette}.}

\item{color_bias}{Define the gradiant of shades among the colors. See \link{colorRampPalette}.}

\item{n_row}{Default 1. Define the number of row of the graph. See \link{facet_wrap}}

\item{n_col}{Default NULL. Define the number of column of the graph. See \link{facet_wrap}}

\item{text_size}{Control the text size of the graph. See \link{ggplot2}.}

\item{legend_size}{Control the legend text size of the graph. see ggplot2. See \link{ggplot2}.}

\item{x_axis_size}{Control the x axis text size of the graph. see ggplot2. See \link{ggplot2}.}

\item{differential_analysis}{Default TRUE. Whether or not use DESeq2 to do a differential abundance analysis. See \link{DESeq2-package}.}

\item{mult_comp}{Default False. Whether or not perform  differential abundance analysis pairwise comparisons.}

\item{selected_comparisons}{Default NULL. A list of vectors of length 2 to subset differential abundance analysis pairwise comparisons. Example: "list(c("0", "7"), c("0", "21"), c("0", "90"))"}

\item{test}{Default 'Wald'. See \link{DESeq2-package}.}

\item{fdr_threshold}{Default '0.05'. Threshold to which taxa bellow are considered significant.}

\item{sig_lab}{Default TRUE. Whether to add stars after taxa name reflecting statistical significance. Only available for 2 groups comparison.}

\item{fitType}{Default "parametric". See \link{DESeq2-package}.}

\item{sfType}{Default "ratio". See \link{DESeq2-package}.}

\item{betaPrior}{Default FALSE.See \link{DESeq2-package}.}

\item{reduced}{Default FALSE. See \link{DESeq2-package}. For example: "~1".}

\item{quiet}{Default TRUE. See \link{DESeq2-package}.}

\item{minReplicatesForReplace}{Default '7'. See \link{DESeq2-package}.}

\item{modelMatrixType}{Default "standard". See \link{DESeq2-package}.}

\item{useT}{Default FALSE. See \link{DESeq2-package}.}

\item{minmu}{See \link{DESeq2-package}.}

\item{parallel}{Default FALSE. See \link{DESeq2-package}.}
}
\value{
A \link{list} containing  \link{ggplot2} '$plot' object and, if differential_analysis = TRUE, \link{data.frame} from \code{\link{results}} DESeq2 function output of significant features '$significant_table_main' and '$significant_table_sub' .
}
\description{
Create a plot of taxa relative abundance, clustered by phylogeny with the possibility to highlight the differentially abundant features.
}
\examples{

\dontrun{

data(ps)

my_plot <- plot_microbiota(
 ps_object = ps,
 exp_group = 'timepoint',
 sample_name = 'SampleID',
 hues = c("Purples", "Blues", "Greens", "Oranges"),
 differential_analysis = T,
 sig_lab = T,
 fdr_threshold = 0.05
)

print(my_plot$plot)
print(my_plot$significant_table_main)
print(my_plot$significant_table_sub)


 }

}
\seealso{
\link{phyloseq}

\link{DESeq2-package}

\link{ggplot2}
}
\author{
Thibault Cuisiniere,  Manuela M. Santos
}
