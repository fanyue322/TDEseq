% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tdeseq.object.R
\docType{class}
\name{TDEseq-class}
\alias{TDEseq-class}
\alias{TDEseq}
\title{Each TDEseq object has a number of slots which store information. Key slots to access
are listed below.}
\description{
Each TDEseq object has a number of slots which store information. Key slots to access
are listed below.
}
\section{Slots}{

\describe{
\item{\code{assays}}{A list of assays for this project}

\item{\code{active.assay}}{Name of the active, or default, assay; settable using
\code{\link{DefaultAssay}}}

\item{\code{project.name}}{Name of the project}

\item{\code{misc}}{A list of miscellaneous information}

\item{\code{version}}{Version of TDEseq this object was built under}

\item{\code{species}}{Could be human, mouse, or rhe, etc}

\item{\code{commands}}{A list of logged commands run on this \code{TDEseq} object}

\item{\code{tools}}{A list of miscellaneous data generated by other tools, should be
filled by developers only using \code{\link{Tool}<-}}
}}

