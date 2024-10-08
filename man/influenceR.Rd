% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/influence_function.R
\name{influenceR}
\alias{influenceR}
\title{Calculate and Plot Measures of Influence for Linear Models}
\usage{
influenceR(model, data, plot = "none")
}
\arguments{
\item{model}{a glm or lm object}

\item{data}{a dataframe}

\item{plot}{optional parameter to include plots in the output, default is "none", options include "none", "all", "Cooks", "DFFITS", "Hadi"}
}
\value{
Returns a dataframe with all influential measures. If `plot` parameter is set to anything other than "none", the function will return a list with dataframe and `ggplot2` plot objects.
}
\description{
This function calculates and evaluates key influence measures, including Cook's Distance, DFFITS, and Hadi's Influence, for linear or generalized linear models by computing diagnostic metrics to identify influential observations. The function has an optional setting that will generate plots to visualize these influence measures, aiding in the detection of data points that may impact the model's fit. This tool is helpful for model assessment and ensuring the reliability of statistical inferences.
}
\examples{
#import dataset mtcars
data1 <- mtcars
model1 <- lm(mpg ~ wt + hp, data = data1)

## run function to create influence measures + plots
sample1 <- influenceR(model1, data1, "all")

## show dataframe of measures
sample1$influence_measures

## show plots
sample1$plots

## show individual plots
sample1$plots$Cooks
sample1$plots$DFFITS
sample1$plots$Hadi

## run function to create dataframe only (no plots)
sample2 <- influenceR(model1, data1)
sample2

## view / manipulate individual influence measure
sample2$DFFITs
}
