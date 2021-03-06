% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/halton.R
\name{halton}
\alias{halton}
\title{Halton Sequences}
\usage{
halton(
  draw,
  dim,
  prime = NULL,
  scramble = F,
  randomise = F,
  normal = F,
  discard = NULL
)
}
\arguments{
\item{draw}{Number of draws for the sequence}

\item{dim}{Number of dimensions to draw}

\item{prime}{If \code{NULL}, the function will use prime numbers in ascending order. If specified, should be an atomic vector of length \code{dim} containing prime numbers to be used for each dimension.}

\item{scramble}{If \code{TRUE}, randomly permute the sequences}

\item{randomise}{If \code{TRUE}, randomise each draw}

\item{normal}{If \code{TRUE}, convert the sequence to a standard normal distribution}

\item{discard}{If provided, simulate an extra number of draws and discard the first few. For example, if 10, \code{draw} + 10 draws will be drawn, and the first 10 discarded.}
}
\value{

}
\description{
Halton Sequences
}
\references{
\insertRef{bhat2003}{rngbox}

\insertRef{train2003}{rngbox}

\insertRef{okten2012}{rngbox}
}
