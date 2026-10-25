setwd("C:/Users/adf44/source/r/frmtmb-wt-sarfcor/extensions/frmtmb.spatial")
rep1 <- function(p, old, new) {
  txt <- paste(readLines(p), collapse = "\n")
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop(p, " matches: ", n, " for ", substr(old, 1, 60))
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], p)
}
rep1("DESCRIPTION", "    rstan,\n    spdep,\n", "    rstan,\n")
rep1("R/sar.R", "  M <- lookup(spec[[\"M\"]], \"sar()\")
  if (inherits(M, c(\"listw\", \"nb\"))) {
    if (!requireNamespace(\"spdep\", quietly = TRUE)) {
      frm_stop(\"sar(): 'M' is a '\", class(M)[1L], \"' object, which the \",
               \"spdep package converts to a matrix, and spdep is not \",
               \"installed. Install it, or pass the weights matrix itself\",
               call. = FALSE)
    }
    M <- if (inherits(M, \"listw\")) {
      spdep::listw2mat(M)
    } else {
      spdep::nb2mat(M)
    }
  }",
"  M <- lookup(spec[[\"M\"]], \"sar()\")
  if (inherits(M, \"listw\")) {
    M <- sar_listw_matrix(M)
  } else if (inherits(M, \"nb\")) {
    M <- sar_nb_matrix(M)
  }")
rep1("R/sar.R", "#' rho from its unbounded parameter:",
"#' The weights matrix of an spdep `nb` object, as brms takes it:
#' `spdep::nb2mat()` with its defaults, the row-standardized binary
#' matrix (style \"W\"), which refuses a location without neighbors. An
#' `nb` object is a list holding each location's neighbor indices, with
#' `0L` for none, so no spdep call is needed to read it.
#'
#' @noRd
sar_nb_matrix <- function(nb) {
  n <- length(nb)
  empty <- vapply(nb, function(v) length(v) == 1L && v[1L] == 0L, TRUE)
  if (any(empty)) {
    frm_stop(\"sar(): the 'nb' object has \", sum(empty), \" location(s) \",
             \"without neighbors, which spdep::nb2mat(), the conversion \",
             \"brms uses, refuses (\\\"Empty neighbour sets found\\\")\",
             call. = FALSE)
  }
  M <- matrix(0, n, n)
  for (i in seq_len(n)) M[i, nb[[i]]] <- 1 / length(nb[[i]])
  M
}

#' The weights matrix of an spdep `listw` object, as
#' `spdep::listw2mat()` gives it: each location's weights at its
#' neighbors' columns.
#'
#' @noRd
sar_listw_matrix <- function(lw) {
  nb <- lw[[\"neighbours\"]]
  wt <- lw[[\"weights\"]]
  n <- length(nb)
  if (!is.list(nb) || !is.list(wt) || length(wt) != n) {
    frm_stop(\"sar(): the 'listw' object needs neighbours and weights, \",
             \"one entry per location\", call. = FALSE)
  }
  M <- matrix(0, n, n)
  for (i in seq_len(n)) {
    if (length(nb[[i]]) == 1L && nb[[i]][1L] == 0L) next
    M[i, nb[[i]]] <- wt[[i]]
  }
  M
}

#' rho from its unbounded parameter:")
rep1("R/sar.R", "#'     `Matrix`, or a `listw` or `nb` object of the spdep package, which
#'     is converted as brms converts it (`spdep::listw2mat()`,
#'     `spdep::nb2mat()`).",
"#'     `Matrix`, or a `listw` or `nb` object of the spdep package, which
#'     is read as brms reads it (`spdep::listw2mat()`, and
#'     `spdep::nb2mat()` with its row-standardized default) without
#'     needing spdep.")
