# functions to create a voronoi shape paths

library(polyclip)

# xy2voronoi ----

xy2voronoi <- function (data, scales, bound = NULL, eps = 1e-09, max.radius = NULL,
          normalize = FALSE, asp.ratio = 1)
{
  data$group <- paste0(seq_len(nrow(data)), ":", data$group)
  if (any(duplicated(data[, c("x", "y")]))) {
    warning("stat_voronoi_tile: dropping duplicated points",
            call. = FALSE)
  }
  polybound <- NULL
  if (is.null(bound)) {
    if (!is.null(max.radius)) {
      bound <- c(range(data$x), range(data$y))
      bound[c(1, 3)] <- bound[c(1, 3)] - max.radius * 1.5
      bound[c(2, 4)] <- bound[c(2, 4)] + max.radius * 1.5
    }
  }
  else if (is.matrix(bound) || is.data.frame(bound)) {
    if (is.matrix(bound) && is.null(colnames(bound))) {
      colnames(bound) <- c("x", "y")
    }
    polybound <- as.data.frame(bound)
    bound <- c(range(polybound$x), range(polybound$y))
  }
  if (normalize) {
    x_range <- range(data$x, na.rm = TRUE, finite = TRUE)
    y_range <- range(data$y, na.rm = TRUE, finite = TRUE)
    data$x <- rescale(data$x, from = x_range) * asp.ratio
    data$y <- rescale(data$y, from = y_range)
    if (!is.null(bound)) {
      bound[1:2] <- rescale(bound[1:2], from = x_range) *
        asp.ratio
      bound[3:4] <- rescale(bound[3:4], from = y_range)
    }
    if (!is.null(polybound)) {
      polybound$x <- rescale(polybound$x, from = x_range) *
        asp.ratio
      polybound$y <- rescale(polybound$y, from = y_range)
    }
  }
  vor <- deldir::deldir(data$x, data$y, rw = bound, eps = eps,
                        suppressMsge = TRUE)
  tiles <- to_tile(vor)
  tiles$orig_x <- data$x[vor$ind.orig[tiles$group]]
  tiles$orig_y <- data$y[vor$ind.orig[tiles$group]]
  tiles$group <- data$group[vor$ind.orig[tiles$group]]
  tiles <- clip_tiles(tiles, max.radius, polybound)
  data$x <- NULL
  data$y <- NULL
  data <- merge(tiles, data, sort = FALSE, all.x = TRUE)
  if (normalize) {
    data$x <- rescale(data$x/asp.ratio, to = x_range, from = c(0,
                                                               1))
    data$y <- rescale(data$y, to = y_range, from = c(0, 1))
  }
  data
}

# clip_tiles ----

clip_tiles <- function (tiles, radius, bound)
{
  if (is.null(radius) && is.null(bound))
    return(tiles)
  p <- seq(0, 2 * pi, length.out = 361)[-361]
  circ <- list(x = cos(p) * radius, y = sin(p) * radius)
  dapply(tiles, "group", function(tile) {
    final_tile <- list(x = tile$x, y = tile$y)
    if (!is.null(radius)) {
      circ_temp <- list(x = circ$x + tile$orig_x[1], y = circ$y +
                          tile$orig_y[1])
      final_tile <- polyclip(final_tile, circ_temp, "intersection")
    }
    if (!is.null(bound)) {
      final_tile <- polyclip(final_tile, bound, "intersection")
    }
    if (length(final_tile) == 0)
      return(NULL)
    new_data_frame(list(x = final_tile[[1]]$x, y = final_tile[[1]]$y,
                        group = tile$group[1]))
  })
}

# dapply ----

dapply <- function (df, by, fun, ..., drop = TRUE)
{
  grouping_cols <- .subset(df, by)
  ids <- id(grouping_cols, drop = drop)
  group_rows <- split(seq_len(nrow(df)), ids)
  rbind_dfs(lapply(seq_along(group_rows), function(i) {
    cur_data <- df_rows(df, group_rows[[i]])
    res <- fun(cur_data, ...)
    if (is.null(res))
      return(res)
    if (length(res) == 0)
      return(new_data_frame())
    vars <- lapply(setNames(by, by), function(col) .subset2(cur_data,
                                                            col)[1])
    if (is.matrix(res))
      res <- split_matrix(res)
    if (is.null(names(res)))
      names(res) <- paste0("V", seq_along(res))
    new_data_frame(modify_list(unclass(vars), unclass(res)))
  }))
}

# df_rows ----

df_rows <- function (x, i)
{
  new_data_frame(lapply(x, `[`, i = i))
}

# id ----

id <- function (.variables, drop = FALSE)
{
  nrows <- NULL
  if (is.data.frame(.variables)) {
    nrows <- nrow(.variables)
    .variables <- unclass(.variables)
  }
  lengths <- vapply(.variables, length, integer(1))
  .variables <- .variables[lengths != 0]
  if (length(.variables) == 0) {
    n <- nrows %||% 0L
    id <- seq_len(n)
    attr(id, "n") <- n
    return(id)
  }
  if (length(.variables) == 1) {
    return(id_var(.variables[[1]], drop = drop))
  }
  ids <- rev(lapply(.variables, id_var, drop = drop))
  p <- length(ids)
  ndistinct <- vapply(ids, attr, "n", FUN.VALUE = numeric(1),
                      USE.NAMES = FALSE)
  n <- prod(ndistinct)
  if (n > 2^31) {
    char_id <- do.call("paste", c(ids, sep = "\r"))
    res <- match(char_id, unique(char_id))
  }
  else {
    combs <- c(1, cumprod(ndistinct[-p]))
    mat <- do.call("cbind", ids)
    res <- c((mat - 1L) %*% combs + 1L)
  }
  if (drop) {
    id_var(res, drop = TRUE)
  }
  else {
    res <- as.integer(res)
    attr(res, "n") <- n
    res
  }
}

# id_var ----

id_var <- function (x, drop = FALSE)
{
  if (length(x) == 0) {
    id <- integer()
    n <- 0L
  }
  else if (!is.null(attr(x, "n")) && !drop) {
    return(x)
  }
  else if (is.factor(x) && !drop) {
    x <- addNA(x, ifany = TRUE)
    id <- as.integer(x)
    n <- length(levels(x))
  }
  else {
    levels <- sort(unique(x), na.last = TRUE)
    id <- match(x, levels)
    n <- max(id)
  }
  attr(id, "n") <- n
  id
}

# modify_list ----

modify_list <- function (old, new)
{
  for (i in names(new)) old[[i]] <- new[[i]]
  old
}

# new_data_frame ----

new_data_frame <- function (x = list(), n = NULL)
{
  if (length(x) != 0 && is.null(names(x)))
    stop("Elements must be named", call. = FALSE)
  lengths <- vapply(x, length, integer(1))
  if (is.null(n)) {
    n <- if (length(x) == 0)
      0
    else max(lengths)
  }
  for (i in seq_along(x)) {
    if (lengths[i] == n)
      next
    if (lengths[i] != 1)
      stop("Elements must equal the number of rows or 1",
           call. = FALSE)
    x[[i]] <- rep(x[[i]], n)
  }
  class(x) <- "data.frame"
  attr(x, "row.names") <- .set_row_names(n)
  x
}

# rbind_dfs ----

rbind_dfs <- function (dfs)
{
  out <- list()
  columns <- unique(unlist(lapply(dfs, names)))
  nrows <- vapply(dfs, .row_names_info, integer(1), type = 2L)
  total <- sum(nrows)
  if (length(columns) == 0)
    return(new_data_frame(list(), total))
  allocated <- rep(FALSE, length(columns))
  names(allocated) <- columns
  col_levels <- list()
  for (df in dfs) {
    new_columns <- intersect(names(df), columns[!allocated])
    for (col in new_columns) {
      if (is.factor(df[[col]])) {
        all_factors <- all(vapply(dfs, function(df) {
          val <- .subset2(df, col)
          is.null(val) || is.factor(val)
        }, logical(1)))
        if (all_factors) {
          col_levels[[col]] <- unique(unlist(lapply(dfs,
                                                    function(df) levels(.subset2(df, col)))))
        }
        out[[col]] <- rep(NA_character_, total)
      }
      else {
        out[[col]] <- rep(.subset2(df, col)[1][NA],
                          total)
      }
    }
    allocated[new_columns] <- TRUE
    if (all(allocated))
      break
  }
  pos <- c(cumsum(nrows) - nrows + 1)
  for (i in seq_along(dfs)) {
    df <- dfs[[i]]
    rng <- seq(pos[i], length.out = nrows[i])
    for (col in names(df)) {
      if (inherits(df[[col]], "factor")) {
        out[[col]][rng] <- as.character(df[[col]])
      }
      else {
        out[[col]][rng] <- df[[col]]
      }
    }
  }
  for (col in names(col_levels)) {
    out[[col]] <- factor(out[[col]], levels = col_levels[[col]])
  }
  attributes(out) <- list(class = "data.frame", names = names(out),
                          row.names = .set_row_names(total))
  out
}

# to_tile ----

to_tile <- function (object)
{
  try_require("deldir", "to_tile")
  tiles <- rbind(structure(object$dirsgs[, c(1:2, 5)], names = c("x",
                                                                 "y", "group")), structure(object$dirsgs[, c(1:2, 6)],
                                                                                           names = c("x", "y", "group")), structure(object$dirsgs[,
                                                                                                                                                  c(3:5)], names = c("x", "y", "group")), structure(object$dirsgs[,
                                                                                                                                                                                                                  c(3:4, 6)], names = c("x", "y", "group")))
  tiles <- unique(tiles)
  tiles <- rbind(tiles, data.frame(x = object$rw[c(1, 2, 2,
                                                   1)], y = object$rw[c(3, 3, 4, 4)], group = deldir::get.cnrind(object$summary$x,
                                                                                                                 object$summary$y, object$rw)))
  tiles$theta <- atan2(tiles$y - object$summary$y[tiles$group],
                       tiles$x - object$summary$x[tiles$group])
  tiles$theta <- ifelse(tiles$theta > 0, tiles$theta, tiles$theta +
                          2 * pi)
  tiles[order(tiles$group, tiles$theta), ]
}

# try_require ----

try_require <- function (package, fun)
{
  if (requireNamespace(package, quietly = TRUE)) {
    return(invisible())
  }
  stop("Package `", package, "` required for `", fun, "`.\n",
       "Please install and try again.", call. = FALSE)
}



