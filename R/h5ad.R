#' @importFrom rhdf5 h5ls
#'
#' @importFrom dplyr .data as_tibble filter select everything
h5ad_ls <-
    function(path)
{
    h5ls(path) |>
        as_tibble() |>
        select("group", "name", "dim", everything())
}

.h5ad_has_raw <-
    function(path)
{
    nrow <-
        h5ad_ls(path) |>
        filter(startsWith(group, "/raw")) |>
        NROW()
    nrow > 0L
}

.h5ad_layers <-
    function(path)
{
    layers <- 
        h5ad_ls(path) |>
        filter(group == "/layers") |>
        pull("name")
    c("X", if (.h5ad_has_raw(path)) "raw", layers)
}

##
## matrix
##

#' @importFrom rlang .env
#'
#' @importFrom dplyr pull
#'
#' @importFrom rhdf5 h5readAttributes
.h5ad_matrix_attribute <-
    function(path, group = "/X")
{
    attributes <- h5readAttributes(path, group)
    non_zero_elements <-
        h5ad_ls(path) |>
        filter(.data$group == .env$group, .data$name == "data") |>
        pull(.data$dim) |>
        as.integer()
    structure(
        c(attributes, list(non_zero_elements = non_zero_elements)),
        class = c("h5ad_matrix_attribute", "h5ad")
    )
}

#' @export
print.h5ad_matrix_attribute <-
    function(x, ...)
{
    cat(
        class(x)[1], "\n",
        "encoding type: ", x[["encoding-type"]], "\n",
        "encoding version: ", x[["encoding-version"]], "\n",
        "shape: ", x[["shape"]][[1]], " x ", x[["shape"]][[2]], "\n",
        "non-zero elements: ", x[["non_zero_elements"]], "\n",
        sep = ""
    )
}

#' @importFrom rhdf5 h5read
#'
#' @importFrom dplyr tibble
#'
#' @importFrom pillar tbl_sum
#'
#' @importFrom dplyr dim_desc
#'
#' @importFrom S4Vectors Rle
.h5ad_csr_matrix <-
    function(path, group = "/X", i = NULL, j = NULL, index = 1:10)
{
    stopifnot(
        `'i=' not yet supported for a 'csr' matrix layer` = is.null(i),
        is.null(j) || is.numeric(j)
    )

    indptr <- h5read(path, paste0(group, "/indptr"), drop = TRUE)
    ## convert 'j' to index elements
    if (is.numeric(j)) {
        starts <- indptr[j] + 1L
        ends <- indptr[j + 1L]
        index <- unlist(Map(seq, starts, ends))
    }
    col <- Rle(
        lengths = diff(indptr),
        values = seq_len(length(indptr) - 1L)
    )[index] |> as.integer()

    row <- h5read(path, paste0(group, "/indices"), list(index), drop = TRUE)
    data <- h5read(path, paste0(group, "/data"), list(index), drop = TRUE)

    tbl <- tibble(index, row, col, data)
    class(tbl) <- c("h5ad_csr_matrix_tbl", class(tbl))
    attr(tbl, "non_zero_elements") <-
        .h5ad_matrix_attribute(path, group)$non_zero_elements

    tbl
}

.h5ad_matrix <-
    function(path, group = "X", i = NULL, j = NULL, index = 1:10)
{
    ## allow i and / or j, but not with index
    if (!is.null(i) || !is.null(j))
        index <- NULL
    group <- paste0("/", group)
    attributes <- .h5ad_matrix_attribute(path, group)
    encoding <- attributes[["encoding-type"]]
    if (identical(encoding, "csr_matrix")) {
        .h5ad_csr_matrix(path, group, i, j, index)
    } else {
        stop("unknown matrix encoding ", encoding)
    }
}

#' @export
tbl_sum.h5ad_csr_matrix_tbl <-
    function(x)
{
    msg <- paste0(
        attr(x, "non_zero_elements"), " non-zero elements"
    )
    c(`An h5ad csr matrix` = msg)
}

##
## data frame
##

.h5ad_dataframe_column_attribute <-
    function(path, group = c("obs", "var"), column)
{
    group <- paste0("/", match.arg(group))

    query <- paste0(group, "/", column)
    attributes <- withCallingHandlers({
        h5readAttributes(path, query)
    }, warning = function(w) {
        ## FIXME: In H5Aread(A, ...) : Reading attribute data of type
        ## 'ENUM' not yet implemented. Values replaced by NA's.
        invokeRestart("muffleWarning")
    })
    structure(
        attributes,
        class = c("h5ad_dataframe_column_attribute", "h5ad")
    )
}

#' @export
print.h5ad_dataframe_column_attribute <-
    function(x, ...)
{
    cat(
        class(x)[1], "\n",
        "encoding type: ", x[["encoding-type"]], "\n",
        "encoding version: ", x[["encoding-version"]], "\n",
        sep = ""
    )
}

.h5ad_dataframe_column_read_array <-
    function(path, group, i, j)
{
    query <- paste0("/", group, "/", j)
    values <- h5read(path, query, list(i), drop = TRUE)
    if (is.factor(values) && setequal(levels(values), c("TRUE", "FALSE")))
        values <- as.logical(values)
    values
}

.h5ad_dataframe_column_read_categorical <-
    function(path, group, i, j)
{
    ## categories
    query <- paste0("/", group, "/", j, "/", "categories")
    categories <- h5read(path, query, drop = TRUE)
    ## codes
    query <- paste0("/", group, "/", j, "/", "codes")
    codes <- h5read(path, query, list(i), drop = TRUE) + 1L

    factor(categories[codes], levels = categories)
}

.h5ad_dataframe_column_read <-
    function(path, group, i, j)
{
    attributes <- .h5ad_dataframe_column_attribute(path, group, j)
    switch(
        attributes[["encoding-type"]],
        "array" = 
            .h5ad_dataframe_column_read_array(path, group, i, j),
        "categorical" =
            .h5ad_dataframe_column_read_categorical(path, group, i, j)
    )
}

.h5ad_dataframe_attribute <-
    function(path, group = c("obs", "var"))
{
    group <- paste0("/", match.arg(group))
    attributes <- h5readAttributes(path, group)
    structure(
        attributes,
        class = c("h5ad_dataframe_attribute", "h5ad")
    )
}

#' @export
print.h5ad_dataframe_attribute <-
    function(x, ...)
{
    cat(
        class(x)[1], "\n",
        "encoding type: ", x[["encoding-type"]], "\n",
        "encoding version: ", x[["encoding-version"]], "\n",
        .pretty(x[["column-order"]], "column order"), "\n",
        sep = ""
    )
}

.h5ad_dataframe_read <-
    function(path, group = c("obs", "var"), i, j)
{
    column_names <- .h5ad_dataframe_attribute(path, group)[["column-order"]]
    if (is.null(j)) {
        j <- column_names
    } else {
        stopifnot(
            all(j %in% column_names)
        )
    }

    values <- lapply(
        j,
        \(name) .h5ad_dataframe_column_read(path, group, i, name)
    )
    names(values) <- j
    as_tibble(values)
}

##
## embedding
##

.h5ad_embedding_attribute <-
    function(path, group = c("obsm", "varm"), embedding)
{
    query <- paste0("/", match.arg(group), "/", embedding)
    attributes <- h5readAttributes(path, query)
}

.h5ad_embedding_read <-
    function(path, group = c("obsm", "varm"), embedding, index = NULL)
{
    query <- paste0("/", group, "/", embedding)
    values <- h5read(path, query, index = index, native = TRUE)
    if (is.null(index[[2]])) {
        colname_idx <- seq_len(ncol(values))
    } else {
        colname_idx <- index[[2]]
    }
    colnames(values) <- paste(embedding, colname_idx, sep = "_")
    as_tibble(values)
}
