% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/H5AD-class.R
\name{h5ad}
\alias{h5ad}
\alias{layer}
\alias{column_data}
\alias{row_data}
\alias{column_embedding}
\alias{row_embedding}
\title{H5AD ('anndata') file representation and access}
\usage{
h5ad(url)

layer(h5ad, layer = "X", i = NULL, j = NULL, index = 1:10)

column_data(h5ad, i = NULL, j = NULL)

row_data(h5ad, i = NULL, j = NULL)

column_embedding(
  h5ad,
  name = .column_embedding_names(h5ad),
  i = NULL,
  j = NULL,
  with = NULL
)

row_embedding(
  h5ad,
  name = .row_embedding_names(h5ad),
  i = NULL,
  j = NULL,
  with = NULL
)
}
\arguments{
\item{url}{\code{character(1)} path to an 'anndata' file}

\item{h5ad}{an object created by \code{h5ad()}, or a URL to an H5AD
file.}

\item{layer}{\code{character(1)} \code{"X"}, \code{"raw"} (if available), or other
layers named in the object and reported by h5ad().}

\item{i}{\code{logical()}, \code{numeric()} (coerced to
integer) or NULL row indexes to be used when extracting
elements (layers, data frames, embeddings) from the H5AD file.}

\item{j}{\code{character()} (for column and row data) or \verb{logical() or }numeric()` (coerced to integer, for layers and embeddings) or
NULL column indexes to be used when extracting elements
from the H5AD file.}

\item{index}{\code{numeric()} (coerced to integer) or NULL indexes into
the native data representation (e.g., \code{1:100} returns the first
100 non-zero elements, column-wise, of a 'csr' matrix). \code{NULL}
indicates that all data (!) should be input. If \code{i} or \code{j} is
defined, then \code{index} is treated as NULL.}

\item{with}{\code{character()} column names of row- or column data
frames corresponding the the row- or column embedding to be
added to the numerical embedding.}
}
\value{
\code{h5ad()} returns an object summarizing the H5AD file. The
object can be used in other calls documented on this help page.

\code{layer()} returns a tibble with row and column indexs and
associated data values.

\code{column_data()} returns a tibble of column (cell)
annotations.

\code{row_data()} returns a tibble of row (gene) annotations.

\code{column_embedding()} returns a tibble summarizing a column
embedding, containing the \code{j} dimensions of the embedding and
perhaps annotated using \verb{with=} columns from the column data.

\code{row_embedding()} returns a tibble summarizing a row
embedding, containing the \code{j} dimensions of the embedding and
perhaps annotated using \verb{with=} columns from the row data.
}
\description{
H5AD ('anndata') file representation and access
}
