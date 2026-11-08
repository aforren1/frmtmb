# Hand resolutions of the merge conflicts of the round of 2026-11-06, one
# lane at a time. dev/rel080-merge.sh calls this after each lane's merge:
#
#   Rscript dev/rel080-resolve.R <lane>
#
# Each resolution names the file, the lane whose merge raised it and how many
# conflict blocks it expects, and stops if the file does not hold exactly
# those. NEWS files and the backlog are merged with --union and edited by
# dev/rel080-news.R; a man/*.Rd conflict takes the release side, since
# roxygen writes every Rd file again from R/.
args <- commandArgs(trailingOnly = TRUE)
lane <- args[1L]
REL <- "C:/Users/adf44/source/r/frmtmb-wt-rel080"
rd <- function(p) readLines(file.path(REL, p), warn = FALSE)
wr <- function(x, p) {
  con <- file(file.path(REL, p), open = "wb")
  on.exit(close(con))
  writeLines(x, con, sep = "\r\n")
}
blocks <- function(x) {
  st <- grep("^<<<<<<< ", x)
  mid <- grep("^=======$", x)
  en <- grep("^>>>>>>> ", x)
  stopifnot(length(st) == length(mid), length(mid) == length(en),
            all(st < mid), all(mid < en))
  lapply(seq_along(st), function(i) {
    list(st = st[i], en = en[i],
         a = if (mid[i] - st[i] > 1L) x[(st[i] + 1L):(mid[i] - 1L)] else character(),
         b = if (en[i] - mid[i] > 1L) x[(mid[i] + 1L):(en[i] - 1L)] else character())
  })
}
# replace block i by f(a, b), where a is the release side and b the lane's
resolve <- function(p, n, fs) {
  x <- rd(p)
  bl <- blocks(x)
  if (length(bl) != n) stop(p, ": expected ", n, " conflicts, found ", length(bl))
  out <- character()
  prev <- 1L
  for (i in seq_along(bl)) {
    if (bl[[i]]$st > prev) out <- c(out, x[prev:(bl[[i]]$st - 1L)])
    out <- c(out, fs[[i]](bl[[i]]$a, bl[[i]]$b))
    prev <- bl[[i]]$en + 1L
  }
  if (prev <= length(x)) out <- c(out, x[prev:length(x)])
  stopifnot(!any(grepl("^(<<<<<<<|=======$|>>>>>>>)", out)))
  if (grepl("[.]R$", p)) invisible(parse(text = out, keep.source = FALSE))
  wr(out, p)
  cat("RESOLVED", lane, p, n, "\n")
}
rel <- function(a, b) a
both <- function(a, b) c(a, b)
lane_first <- function(a, b) c(b, a)

switch(lane,
  mfxse = {
    # latperf's paragraphs on the sparse joint covariance, with mfxse's
    # sentence on the wider block a carried `vcov =` inverts
    resolve("R/autocor.R", 1L, list(function(a, b) {
      stopifnot(any(grepl("held as the sparse", a)),
                any(grepl("The same option bounds a wider block", b)))
      c(
"#' The joint covariance of the estimates and the latents is held as the sparse",
"#' Cholesky factor of their joint precision, which is banded along a chain, so",
"#' a standard error reads it at a cost that grows with the rows, not with",
"#' their square. [fitted()], [residuals()] and `frm_linpred(se.fit = TRUE)` on",
"#' the fitted rows read each row's variance through the selected inverse of",
"#' the factor: on one poisson series, 0.4 s at 8000 rows and 1.7 s at 32000.",
"#' A hurdle or zero-inflated mean reads the latent of `mu` and a level of",
"#' `hu ~ (1 | id)` or `zi ~ (1 | id)` on one row, a pair that no row's",
"#' likelihood couples; that pair is read from a column of the inverse solved",
"#' once per level, 3.3 s at 32000 rows. Rows in `newdata` read the estimates",
"#' and the levels of any other group-level block alone, [frm_curve()] and",
"#' [conditional_effects()] too.",
"#'",
"#' A dense block of the joint covariance is formed only where a function",
"#' returns one or must form one: [frm_joint_cov()], [frm_lp_basis()] on the",
"#' fitted rows, and the category probabilities of an ordinal family's fitted",
"#' rows. A block of more than `getOption(\"frmtmb.latres_dense_rows\", 20000L)^2`",
"#' numbers is refused, with its size in the message. The same option bounds a",
"#' wider block: on any fit with group-level effects, latent or not, a `vcov =`",
"#' that marginaleffects carries into those effects inverts the inner block, one",
"#' row per effect, densely, and past that many inner rows the effects are held",
"#' at their estimates with a warning ([frmtmb-marginaleffects]). Above 3000",
"#' rows the category probabilities' standard errors are `NA`, with a message of",
"#' class `frmtmb_latres_se_skipped`; `options(frmtmb.latres_se_rows = )` moves",
"#' that limit. A fit whose covariance of the estimates is held as a factor,",
"#' which an ill-conditioned design gives, keeps the dense joint covariance, and",
"#' both limits apply to every standard error of its fitted rows.")
    }))
    resolve("man/frmtmb-autocor.Rd", 1L, list(rel))
  },
  arstep = {
    # brief step 2: arstep's second descent ends the Laplace branch, before
    # nestquad's nested-quadrature branch
    resolve("R/fit.R", 1L, list(function(a, b) {
      stopifnot(grepl("else if \\(!is.null\\(qnest\\)\\)", a[1L]),
                grepl("gp_second_descent", b[1L]), length(b) == 2L)
      c(b, a)
    }))
  },
  fitedge = {
    # the same two gp_second_descent() lines, already in from arstep
    resolve("R/fit.R", 1L, list(function(a, b) {
      stopifnot(grepl("else if \\(!is.null\\(qnest\\)\\)", a[1L]), !length(b))
      a
    }))
    # rev-fitedge-src/merge2/resolve.R (the fitedge review's tested
    # resolution of diff3 -m -E fitedge base arstep), as written there
    p <- "R/gp-start.R"
    x <- rd(p)
    st <- grep("^<<<<<<< ", x); mid <- grep("^=======$", x); en <- grep("^>>>>>>> ", x)
    stopifnot(length(st) == 5L, length(mid) == 5L, length(en) == 5L)
    side <- function(i, which) {
      if (which == "f") x[(st[i] + 1L):(mid[i] - 1L)] else
        if (mid[i] + 1L <= en[i] - 1L) x[(mid[i] + 1L):(en[i] - 1L)] else character(0)
    }
    f <- lapply(1:5, side, "f"); a <- lapply(1:5, side, "a")
    res <- list()
    res[[1]] <- c(f[[1]], a[[1]][grepl("ard", a[[1]])])
    r2 <- f[[2]][seq_len(grep("@return", f[[2]]) - 1L)]
    res[[2]] <- c(r2,
      "#' @return `list(idx, theta, alts, body)`: the positions in the template's",
      "#'   `theta`, their values, the full `theta` vectors at each other local",
      "#'   maximum the working marginal reached (gp_marginal_alts()), and which of",
      "#'   the positions the body start set, or `NULL`.")
    res[[3]] <- c(a[[3]], f[[3]])
    stopifnot(grepl("^  if \\(length\\(idx\\)\\) list", f[[4]][1]))
    res[[4]] <- c(
      "  if (length(idx)) list(idx = idx, theta = theta, alts = alts,",
      "                        body = from_body)", f[[4]][-1L])
    tail5 <- c(
      "  out <- gp_start_scan_units(obj, template, frame, val, covered)",
      "  body <- intersect(covered, w[[\"body\"]])",
      "  if (!length(covered) || (!length(w$alts) && !length(body))) {",
      "    return(out)",
      "  }",
      "  saved <- obj_state_save(obj)",
      "  on.exit(obj_state_restore(obj, saved), add = TRUE)",
      "  alt <- list()",
      "  if (length(w$alts)) {",
      "    tpl <- template",
      "    tpl[[\"theta\"]][w$idx] <- w$alts[[1L]]",
      "    if (!is.na(val(tpl))) {",
      "      alt <- c(alt, list(outer_from_template(tpl, obj, frame)))",
      "    }",
      "  }",
      "  if (length(body)) {",
      "    keep <- setdiff(covered, body)",
      "    scanned[[\"theta\"]][keep] <- template[[\"theta\"]][keep]",
      "    p <- gp_start_scan_units(obj, scanned, frame, val, keep)$par",
      "    if (length(p) == length(out$par) && all(is.finite(p))) {",
      "      alt <- c(alt, list(p))",
      "    }",
      "  }",
      "  if (length(alt)) out$env$.frmtmb_gp_alt <- alt",
      "  out",
      "}")
    gsd <- grep("^#' A second descent from", f[[5]])
    res[[5]] <- c(tail5, "", f[[5]][gsd:length(f[[5]])],
                  a[[5]][(grep("^#' Whether the objective at the marginal", a[[5]]) - 2L):
                           length(a[[5]])])
    out <- character(0); prev <- 1L
    for (i in 1:5) {
      out <- c(out, x[prev:(st[i] - 1L)], res[[i]])
      prev <- en[i] + 1L
    }
    out <- c(out, x[prev:length(x)])
    k <- grep("^      from_body <- c\\(from_body, from\\)$", out)
    stopifnot(length(k) == 1L, grepl("theta <- c\\(theta, m\\$theta\\)", out[k - 1L]))
    j <- grep("^      idx <- c\\(idx, from\\)$", out)
    j <- j[j < k & j > k - 4L]
    stopifnot(length(j) == 1L)
    out <- append(out, "      alts <- lapply(alts, function(th) c(th, m$theta))", j - 1L)
    invisible(parse(text = out, keep.source = FALSE))
    stopifnot(sum(grepl("^gp_second_descent <- ", out)) == 1L,
              sum(grepl("^gp_marginal_agrees <- ", out)) == 1L,
              !any(grepl("^(<<<<<<<|=======|>>>>>>>)", out)))
    wr(out, p)
    cat("RESOLVED fitedge R/gp-start.R 5 (merge2/resolve.R)\n")
  },
  bigx = {
    # brief step 7: bigx's body-basis descent after arstep's second descent,
    # in the Laplace branch, before nestquad's
    resolve("R/fit.R", 1L, list(function(a, b) {
      stopifnot(grepl("gp_second_descent", a[1L]),
                grepl("else if \\(!is.null\\(qnest\\)\\)", a[3L]),
                any(grepl("nl_basis_descent", b)))
      c(a[1:2], b, a[-(1:2)])
    }))
  },
  nlstart = {
    # loglik's word for the value, nlstart's class on the message
    resolve("R/par-template.R", 1L, list(function(a, b) {
      stopifnot(grepl("^      llw, ", a[1L]), grepl("frmtmb_nl_restart", b[3L]))
      c(a[1L], sub(")$", ",", a[2L]), b[3L])
    }))
  },
  superassign = {
    # superassign review, MINOR 3: readers' me_node_value() (no memo, so no
    # `<<-` and no environment)
    resolve("R/me.R", 2L, list(rel, rel))
  },
  cat("nothing to resolve for", lane, "\n"))
