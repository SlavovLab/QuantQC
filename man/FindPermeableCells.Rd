% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PermeabilityClassifyer.R
\name{FindPermeableCells}
\alias{FindPermeableCells}
\title{gives probability cell is permeable}
\usage{
FindPermeableCells(mat, species = "Human")
}
\arguments{
\item{mat}{an expression matrix (can have NAs) proteins rows cellID columns}

\item{species}{Human or Mouse}
}
\value{
A vector of probabilities length number columns of matrix
}
\description{
This function takes a matrix as an object with cell IDs on columns and uniprot IDs on rownames
}
\examples{
FindPermeableCells(protein_matrix, species = 'Mouse')
}
