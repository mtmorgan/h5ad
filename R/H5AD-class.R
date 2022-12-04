##
## H5AD-class
##

#' @importFrom methods setClass show is
.H5AD <- setClass(
    "H5AD",
    slots = c(
        url = "character",
        dim = "integer",
        has_raw = "logical",
        layers = "character",
        column_data_names = "character",
        row_data_names = "character",
        column_embedding_names = "character",
        row_embedding_names = "character"
    )
)

.url <- function(x) x@url

.dim <- function(x) x@dim

.has_raw <- function(x) x@has_raw

.layers <- function(x) x@layers

.column_data_names <-  function(x) x@column_data_names

.row_data_names <- function(x) x@row_data_names

.column_embedding_names <- function(x) x@column_embedding_names

.row_embedding_names <- function(x) x@row_embedding_names

setMethod("show", "H5AD", function(object) {
    cat(
        "class: ", class(object), "\n",
        "dim: ", .dim(object)[[1]], " x ", .dim(object)[[2]], "\n",
        .pretty(.layers(object), "layer names"), "\n",
        .pretty(.column_data_names(object), "column_data names"), "\n",
        .pretty(.row_data_names(object), "row_data names"), "\n",
        .pretty(.column_embedding_names(object), "column_embedding"), "\n",
        .pretty(.row_embedding_names(object), "row_embedding"), "\n",
        sep = ""
    )
})

#' @rdname H5AD-class
#'
#' @title H5AD ('anndata') file representation and access
#'
#' @param url `character(1)` path to an 'anndata' file
#'
#' @return `h5ad()` returns an object summarizing the H5AD file. The
#'     object can be used in other calls documented on this help page.
h5ad <-
    function(url)
{
    stopifnot(
        is.character(url), length(url) == 1L, file.exists(url)
    )

    dim <-
        h5readAttributes(url, "/X")[["shape"]] |>
        as.vector() |>
        rev()

    has_raw <- .h5ad_has_raw(url)

    layers <- .h5ad_layers(url)

    column_data_names <-
        h5readAttributes(url, "/obs")[["column-order"]] |>
        as.vector()

    row_data_names <-
        h5readAttributes(url, "/var")[["column-order"]] |>
        as.vector()

    column_embedding_names <-
        h5ad_ls(url) |>
        filter(startsWith(group, "/obsm")) |>
        pull(name)

    row_embedding_names <-
        h5ad_ls(url) |>
        filter(startsWith(group, "/varm")) |>
        pull(name)

    .H5AD(
        url = url,
        dim = dim,
        has_raw = has_raw,
        layers = layers,
        column_data_names = column_data_names,
        row_data_names = row_data_names,
        column_embedding_names = column_embedding_names,
        row_embedding_names = row_embedding_names
    )
}

setAs("character", "H5AD", function(from) h5ad(from))

#' @rdname H5AD-class
#'
#' @param h5ad an object created by `h5ad()`, or a URL to an H5AD
#'     file.
#'
#' @param layer `character(1)` `"X"`, `"raw"` (if available), or other
#'     layers named in the object and reported by h5ad().
#'
#' @param i `numeric()` (coerced to integer) or NULL row indexes to be
#'     used when extracting elements (layers, data frames, embeddings)
#'     from the H5AD file.
#'
#' @param j `numeric()` (coerced to integer) or NULL column indexes to
#'     be used when extracting elements (layers, data frames,
#'     embeddings) from the H5AD file.
#'
#' @param index `numeric()` (coerced to integer) or NULL indexes into
#'     the native data representation (e.g., `1:100` returns the first
#'     100 non-zero elements, column-wise, of a 'csr' matrix). `NULL`
#'     indicates that all data (!) should be input. If `i` or `j` is
#'     defined, then `index` is treated as NULL.
#'
#' @return `layer()` returns a tibble with row and column indexs and
#'     associated data values.
#' 
#' @export
layer <-
    function(h5ad, layer = "X", i = NULL, j = NULL, index = 1:10)
{
    h5ad <- as(h5ad, "H5AD")
    stopifnot(
        is.character(layer) &&
        length(layer) == 1L &&
        layer %in% .h5ad_layers(.url(h5ad)),
        is.null(i) || is.numeric(i),
        is.null(j) || is.numeric(j)
    )        

    if (identical(layer, "raw")) {
        layer <- "raw/X"
    } else if (!identical(layer, "X")) {
        layer <- paste0("layers/", layer)
    }

    .h5ad_matrix(.url(h5ad), layer, i, j, index = index)
}

#' @rdname H5AD-class
#'
#' @return `column_data()` returns a tibble of column (cell)
#'     annotations.
#'
#' @export
column_data <-
    function(h5ad, i = NULL, j = NULL)
{
    h5ad <- as(h5ad, "H5AD")
    stopifnot(
        is.null(i) || is.integer(i),
        is.null(j) || is.character(j)
    )

    .h5ad_dataframe_read(.url(h5ad), "obs", i, j)
}

#' @rdname H5AD-class
#'
#' @return `row_data()` returns a tibble of row (gene) annotations.
#'
#' @export
row_data <-
    function(h5ad, i = NULL, j = NULL)
{
    h5ad <- as(h5ad, "H5AD")
    stopifnot(
        is.null(i) || is.integer(i),
        is.null(j) || is.character(j)
    )

    .h5ad_dataframe_read(.url(h5ad), "var", i, j)
}

#' @rdname H5AD-class
#'
#' @param with `character()` column names of row- or column data
#'     frames corresponding the the row- or column embedding to be
#'     added to the numerical embedding.
#'
#' @return `column_embedding()` returns a tibble summarizing a column
#'     embedding, containing the `j` dimensions of the embedding and
#'     perhaps annotated using `with=` columns from the column data.
#'
#' @export
column_embedding <-
    function(
        h5ad,
        name = .column_embedding_names(h5ad), i = NULL, j = NULL,
        with = NULL)
{
    h5ad <- as(h5ad, "H5AD")
    stopifnot(
        is.null(with) || is.character(with)
    )
    name <- match.arg(name)
    if (is.na(name))
        stop("no column embeddings defined")

    embedding <- .h5ad_embedding_read(.url(h5ad), "obsm", name, list(i, j))

    if (!is.null(with)) {
        annotation <- column_data(h5ad, i = i, j = with)
        embedding <- bind_cols(embedding, annotation)
    }

    embedding
}

#' @rdname H5AD-class
#'
#' @return `row_embedding()` returns a tibble summarizing a row
#'     embedding, containing the `j` dimensions of the embedding and
#'     perhaps annotated using `with=` columns from the row data.
#'
#' @export
row_embedding <-
    function(
        h5ad,
        name = .row_embedding_names(h5ad), i = NULL, j = NULL,
        with = NULL)
{
    h5ad <- as(h5ad, "H5AD")
    stopifnot(
        is.null(with) || is.character(with)
    )
    name <- match.arg(name)
    if (is.na(name))
        stop("no row embeddings defined")

    .h5ad_embedding_read(.url(h5ad), "varm", name, list(i, j))
}
