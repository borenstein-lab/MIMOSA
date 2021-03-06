% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postProcessing.R
\name{plot_metabolite_counts}
\alias{plot_metabolite_counts}
\title{Make bar plot of MIMOSA outcome counts in different metabolite categories}
\usage{
plot_metabolite_counts(node_data, variable, threshold = "low",
  prediction_colors = c("#1B9E77", "lightgrey", "#F46D43"))
}
\arguments{
\item{node_data}{Processed core MIMOSA results table}

\item{variable}{Grouping variable for metabolites (must be a column name in node_data)}

\item{threshold}{Whether to use the lower or higher significance threshold}

\item{prediction_colors}{List of 3 colors to use as color scheme for Consistent, Contrasting, & Inconsistent metabolites}
}
\value{
plot object of MIMOSA outcome counts across specified metabolite categories
plot_metabolite_counts(node_data, "Category")
}
\description{
Make bar plot of MIMOSA outcome counts in different metabolite categories
}
