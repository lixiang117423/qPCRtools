% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CalExp2dCt.R
\name{CalExp2dCt}
\alias{CalExp2dCt}
\title{Calculate expression using standard curve.}
\arguments{
\item{cq.table}{The data frame of the position and cq value.}

\item{design.table}{The data frame of the position and corresponding information.}

\item{ref.gene}{The name of reference gene.}
}
\value{
A list contain a table and a figure.
}
\description{
Calculate expression using standard curve.
}
\examples{
df1.path <- system.file("examples", "dct.cq.txt", package = "qPCRtools")
df2.path <- system.file("examples", "dct.design.txt", package = "qPCRtools")
cq.table <- read.table(df1.path, sep = ",", header = TRUE)
design.table <- read.table(df2.path, sep = ",", header = TRUE)
CalExp2dCt(cq.table,
           design.table,
           ref.gene = "Actin"
) -> res

}
\author{
Xiang LI <lixiang117423@gmail.com>
}
