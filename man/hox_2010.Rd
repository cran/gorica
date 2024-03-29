% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hox_2010.R
\docType{data}
\name{hox_2010}
\alias{hox_2010}
\title{Sesame Street data based on Hox (2010)}
\format{
A data frame with 2000 rows and 6 variables.
}
\usage{
data(hox_2010)
}
\description{
Synthetic data based on Hox (2010, p. 16). In the study, the outcome variable popular
represents the popularity score of pupils, ranging from 0 (very unpopular) to 10 (very
popular), for pupils nested in 100 classes of varying size. The popularity scores are
predicted by pupil level predictors gender (G) and pupil extraversion scores (PE) that
range from 1 (introversion) to 10 (extraversion), a class-level predictor teacher
experience (TE), and the cross-level interactions between G and TE as well as PE and TE.
Since standardization is recommended when the model contains interactions, we standardize
PS, PE and TE by means of grand mean centering. That is, we first substract the overall
means of the continuous variables PS, PE, and TE from each of their values, before dividing
these values by their standard deviations.
}
\details{
\tabular{lll}{
  \strong{ID} \tab \code{integer} \tab Pupil ID\cr
  \strong{class} \tab \code{integer} \tab Class ID\cr
  \strong{PE} \tab \code{numeric} \tab Pupil extraversion, standardized\cr
  \strong{G} \tab \code{factor} \tab Pupil sex\cr
  \strong{PS} \tab \code{numeric} \tab Popularity scores, standardized\cr
  \strong{TE} \tab \code{integer} \tab Teacher experience, standardized
}
}
\references{
Hox, J. J. (2010). Multilevel analysis: Techniques and
applications (2nd ed.). New York, NY: Routledge.
}
\keyword{datasets}
