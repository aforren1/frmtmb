# The findings of the release review (dev/reviews/2026-11-04-release.md), each with its
# number. It runs last in dev/rel079-merge.sh, after the NEWS merge and before the line
# endings go back, and it was run once on the tree, so the two agree. Each edit matches
# once or the run stops, and a file keeps its line endings.
#   Rscript dev/rel079-release-fix2.R
rel <- Sys.getenv("REL079_ROOT", "C:/Users/adf44/source/r/frmtmb-wt-rel079")
setwd(rel)
latres_tree <- "C:/Users/adf44/source/r/frmtmb-wt-latres"
ed <- function(f, old, new, n = 1L) {
  raw <- readBin(f, "raw", file.size(f))
  crlf <- length(raw) > 1L && any(raw == as.raw(13L))
  s <- rawToChar(raw)
  Encoding(s) <- "bytes"
  s <- gsub("\r\n", "\n", s, fixed = TRUE, useBytes = TRUE)
  # a second run on the tree, after a stop, skips what the first applied
  if (grepl(new, s, fixed = TRUE, useBytes = TRUE)) {
    cat("already", f, substr(new, 1L, 40L), "\n")
    return(invisible())
  }
  k <- sum(gregexpr(old, s, fixed = TRUE, useBytes = TRUE)[[1L]] > 0L)
  if (k != n) stop(f, ": expected ", n, " match, found ", k, ": ", substr(old, 1L, 70L))
  s <- gsub(old, new, s, fixed = TRUE, useBytes = TRUE)
  if (crlf) s <- gsub("\n", "\r\n", s, fixed = TRUE, useBytes = TRUE)
  writeBin(charToRaw(s), f)
  cat("edited", f, "\n")
}

# ---- 1. fitted(newdata) beside another block ---------------------------------------
# The differences hold beta and theta at a fixed b, which is right on a latent alone
# and wrong beside any other block, by 1.12 to 3.7 times above the row limit. So beside
# a block the joint covariance is read at any size up to frmtmb.latres_dense_rows, and
# finding 13's handler reads NA, with a message, above that.
ed("R/predict.R",
   paste0("  if (!is.null(newdata) && !is.matrix(est) && scale == \"response\" &&\n",
          "        length(frame_latres(object$frame)) &&\n",
          "        !isTRUE(object$frame[[\"latres_off\"]])) {\n",
          "    fd <- as.vector(fit_fd_se(object, function(fit) {\n"),
   paste0("  if (!is.null(newdata) && !is.matrix(est) && scale == \"response\" &&\n",
          "        length(frame_latres(object$frame)) &&\n",
          "        !isTRUE(object$frame[[\"latres_off\"]]) &&\n",
          "        length(object$frame[[\"re_blocks\"]]) == 1L) {\n",
          "    fd <- as.vector(fit_fd_se(object, function(fit) {\n"))

# ---- 2. marginaleffects on a multivariate fit --------------------------------------
ed("R/interop.R",
   paste0("  if (type == \"response\") {\n",
          "    mm <- fitted_me_moments(model, newdata, dots[[\"re_formula\"]],\n",
          "                            dots[[\"resp\"]], dpar,\n",
          "                            isTRUE(dots[[\"allow_new_levels\"]]),\n",
          "                            mean_only = TRUE)\n",
          "    if (!is.null(mm)) {\n",
          "      attributes(mm$mean) <- attributes(p)\n",
          "      p <- mm$mean\n",
          "    } else {\n",
          "      p <- latres_fitted_mean(model, newdata, dots[[\"re_formula\"]],\n",
          "                              dots[[\"resp\"]], dpar,\n",
          "                              isTRUE(dots[[\"allow_new_levels\"]]), p)\n",
          "    }\n",
          "  }\n"),
   paste0("  if (type == \"response\") {\n",
          "    read_mean <- function(r, value) {\n",
          "      anl <- isTRUE(dots[[\"allow_new_levels\"]])\n",
          "      mm <- fitted_me_moments(model, newdata, dots[[\"re_formula\"]], r, dpar, anl,\n",
          "                              mean_only = TRUE)\n",
          "      if (is.null(mm)) {\n",
          "        return(latres_fitted_mean(model, newdata, dots[[\"re_formula\"]], r, dpar,\n",
          "                                  anl, value))\n",
          "      }\n",
          "      attributes(mm$mean) <- attributes(value)\n",
          "      mm$mean\n",
          "    }\n",
          "    # with no `resp` the frame holds every response, and each has its own load\n",
          "    rns <- names(model$spec$responses)\n",
          "    if (is.null(dots[[\"resp\"]]) && length(rns) > 1L &&\n",
          "          identical(colnames(p), rns)) {\n",
          "      for (r in rns) p[, r] <- read_mean(r, p[, r])\n",
          "    } else {\n",
          "      p <- read_mean(dots[[\"resp\"]], p)\n",
          "    }\n",
          "  }\n"))

# ---- 3. me() on a predictor other than mu beside a latent ---------------------------
ed("R/me.R",
   "  tot <- function(mp) Reduce(`+`, lapply(mp, function(e) (e$g * e$sd)^2))\n",
   paste0("  if (!is.null(sl) && !identical(names(parts)[1L], linpred_key(rn, \"mu\"))) {\n",
          "    # one node cannot read two independent draws, and the latent moves the mean\n",
          "    # by exp(s^2 / 2), where this me() value moves it by 0.1 percent\n",
          "    key <- linpred_key(rn, \"mu\")\n",
          "    s <- rep_len(sl, nrow(newdata))\n",
          "    return(list(key = key, s2 = s^2, by_key = stats::setNames(list(s), key)))\n",
          "  }\n",
          "  tot <- function(mp) Reduce(`+`, lapply(mp, function(e) (e$g * e$sd)^2))\n"))

# ---- 4. latres RF1's number, read on the final build ---------------------------------
ed("NEWS.md", "warned on 0.864 of the wrong fits", "warned on 0.865 of the wrong fits")
ed("R/autocor.R", "#' | all | this rule | 0.864 | 0.157 |", "#' | all | this rule | 0.865 | 0.157 |")

# ---- 6. the help the merge lost ----------------------------------------------------
lost <- sub("\r$", "", readLines(file.path(latres_tree, "R/autocor.R"), warn = FALSE))
a <- which(startsWith(lost, "#' [frm()] refuses `cov = FALSE`"))
b <- a[1L] + which(!startsWith(lost[(a[1L] + 1L):length(lost)], "#'"))[1L] - 1L
stopifnot(length(a) >= 1L, grepl("frmtmb_latres_dispersion|mi\\(\\)", paste(lost[a[1L]:b], collapse = " ")))
# what follows the first @section stays where it is
sec <- which(startsWith(lost[a[1L]:b], "#' @section "))
if (length(sec)) b <- a[1L] + sec[1L] - 2L
while (lost[b] == "#'") b <- b - 1L
cat("the help restored:", b - a[1L] + 1L, "lines\n")
ed("R/autocor.R",
   "#' fault is the Wald interval and not the approximation.\n",
   paste0("#' fault is the Wald interval and not the approximation.\n#'\n",
          paste(lost[a[1L]:b], collapse = "\n"), "\n"))

# ---- 7. frmtmb.latent's floor -------------------------------------------------------
ed("extensions/frmtmb.latent/DESCRIPTION", "Version: 0.7.1", "Version: 0.7.2")
ed("extensions/frmtmb.latent/DESCRIPTION", "frmtmb (>= 0.75.0)", "frmtmb (>= 0.79.0)")
ed("extensions/frmtmb.latent/NEWS.md", "# frmtmb.latent 0.7.1\n",
   paste0("# frmtmb.latent 0.7.2\n\nNeeds frmtmb 0.79.0. The package itself is unchanged, and one of its\n",
          "tests pins the interval of `fitted()` on a `mixture_mvn()` fit that\nfrmtmb 0.79.0 gives.\n\n",
          "# frmtmb.latent 0.7.1\n"))
ed("NEWS.md",
   "- Update the extensions together with frmtmb. frmtmb.sample 0.26.0 and\n  frmtmb.learn 0.7.4 need frmtmb 0.79.0.\n",
   "- Update the extensions together with frmtmb. frmtmb.sample 0.26.0,\n  frmtmb.learn 0.7.4 and frmtmb.latent 0.7.2 need frmtmb 0.79.0.\n")

# ---- 8. a category's chain step beside a latent ---------------------------------------
ed("R/autocor-latent.R",
   paste0("    linkfun <- fam[[\"links\"]][[\"mu\"]][[\"linkfun\"]]\n",
          "    function(fam, dpars, aterms) {\n",
          "      fit[[\"eta_shift\"]] <- list(key = linpred_key(rn, \"mu\"),\n",
          "                                 h = linkfun(dpars[[\"mu\"]]) - eta)\n"),
   paste0("    linkfun <- fam[[\"links\"]][[\"mu\"]][[\"linkfun\"]]\n",
          "    # eta holds fit_fd_se()'s chain step, and this list replaces it: without it\n",
          "    # d out / d eta reads 0, which drops a gp()'s kriging variance\n",
          "    key <- linpred_key(rn, \"mu\")\n",
          "    h0 <- if (identical(fit[[\"eta_shift\"]][[\"key\"]], key)) {\n",
          "      fit[[\"eta_shift\"]][[\"h\"]]\n",
          "    } else {\n",
          "      0\n",
          "    }\n",
          "    function(fam, dpars, aterms) {\n",
          "      fit[[\"eta_shift\"]] <- list(key = key,\n",
          "                                 h = linkfun(dpars[[\"mu\"]]) - eta + h0)\n"))

# ---- 9. the linear scale on a latent, and 13. a refused standard error says so --------
ed("R/predict.R",
   "  if (is.null(mm)) return(p$se.fit)\n",
   paste0("  if (is.null(mm)) {\n",
          "    # a new row's latent is in brms's linear predictor, as a me() draw is here\n",
          "    sl <- if (!is.null(newdata) && scale == \"linear\") {\n",
          "      latres_row_sd(object, resp %||% names(object$spec$responses)[1L], dpar)\n",
          "    }\n",
          "    return(if (is.null(sl)) p$se.fit else sqrt(p$se.fit^2 + sl^2))\n",
          "  }\n"))
ed("R/predict.R",
   paste0("        invokeRestart(\"muffleWarning\")\n",
          "      }\n",
          "    }),\n",
          "    error = keep_missing(NULL))\n",
          "  if (is.null(p)) return(NULL)\n"),
   paste0("        invokeRestart(\"muffleWarning\")\n",
          "      }\n",
          "    }),\n",
          "    # a refusal above the dense limit is said, not read as no standard error\n",
          "    error = function(e) {\n",
          "      if (inherits(e, \"frmtmb_latres_dense\")) {\n",
          "        frm_message(conditionMessage(e), class = \"frmtmb_latres_se_skipped\")\n",
          "      }\n",
          "      keep_missing(NULL)(e)\n",
          "    })\n",
          "  if (is.null(p)) return(NULL)\n"))
ed("R/predict.R",
   "\"and predict() do not read it\", call. = FALSE)\n",
   "\"and predict() do not read it\", call. = FALSE,\n               class = \"frmtmb_latres_dense\")\n")

# ---- 10. no mu, no latent residual ----------------------------------------------------
no_mu <- paste0("\"Refused as brms refuses it ('Explicit covariance terms can only be specified ",
                "on 'mu''): this family has no mu.\"")
for (fm in c("lba", "rdm")) {
  ed("extensions/frmtmb.eam/R/zzz.R",
     sprintf("  r(\"%s\", \"nl\", \"refused\",\n", fm),
     sprintf("  r(\"%s\", \"group:autocor\", \"refused\",\n    %s)\n  r(\"%s\", \"nl\", \"refused\",\n", fm, no_mu, fm))
}
ed("extensions/frmtmb.learn/R/zzz.R",
   "  r(nm, \"nl\", \"refused\",\n",
   sprintf("  r(nm, \"group:autocor\", \"refused\",\n    %s)\n  r(nm, \"nl\", \"refused\",\n", no_mu))
ed("extensions/frmtmb.eam/NEWS.md",
   "`nl = TRUE` needs a family with a single `mu`, and `nlf(v2 ~ ...)` is\n  the route.\n",
   paste0("`nl = TRUE` needs a family with a single `mu`, and `nlf(v2 ~ ...)` is\n  the route.\n\n",
          "* `frm_compat()`: `ar()`, `ma()`, `arma()`, `cosy()` and `unstr()` by `lba` and\n",
          "  `rdm` are \"refused\", as `frm()` and brms refuse them: these families have no\n",
          "  `mu`.\n"))
ed("extensions/frmtmb.learn/NEWS.md",
   "`nlf()` on one parameter is the\n  route.\n",
   paste0("`nlf()` on one parameter is the\n  route.\n\n",
          "* `frm_compat()`: `ar()`, `ma()`, `arma()`, `cosy()` and `unstr()` by these\n",
          "  families are \"refused\", as `frm()` and brms refuse them: they have no `mu`.\n"))

# ---- 11. the compat note --------------------------------------------------------------
ed("R/compat.R",
   paste0("The condition: frm() warns (class frmtmb_latres_laplace) where a row tells its latent too little for ",
          "Laplace: binary or ordinal rows, at most 3 trials, a mean share of the latent below 0.18, or a family ",
          "with an estimated dispersion that trades off against the latent."),
   paste0("The condition: frm() warns (class frmtmb_latres_laplace) on ordinal rows, and where a first-order ",
          "Laplace correction moves log sderr or ar[1] by more than 0.15 or cannot be read. Above 3000 rows it ",
          "reads the rows' share of the latent."),
   n = 1L)

# ---- 12. the regression, in NEWS --------------------------------------------------------
# the top section's own heading: older sections have one too
nw <- sub("\r$", "", readLines("NEWS.md", warn = FALSE))
if (!any(grepl("conditional_effects()` on a nonlinear body whose parameter reads", nw, fixed = TRUE))) {
  k <- which(nw == "## New features")[1L]
  stopifnot(k < which(startsWith(nw, "# frmtmb 0.78.0"))[1L])
  nw <- c(nw[seq_len(k - 1L)],
          "* `conditional_effects()` on a nonlinear body whose parameter reads `me()`",
          "  now stops, with \"frm_lp_basis(): the nonlinear body reaches 'a'\". It",
          "  answered on 0.78.0. A `gp()` in such a parameter stops on both, and a",
          "  smooth answers on both.", "",
          nw[k:length(nw)])
  raw <- readBin("NEWS.md", "raw", file.size("NEWS.md"))
  con <- file("NEWS.md", "wb")
  writeLines(nw, con, sep = if (any(raw == as.raw(13L))) "\r\n" else "\n")
  close(con)
  cat("edited NEWS.md: the regression\n")
}

# ---- 18. me(x) * me(z) -----------------------------------------------------------------
ed("R/predict.R",
   paste0("#' At a row whose `mi()` value left the objective, a reader's fit reads the\n"),
   paste0("#' The product of two `me()` values on one row is read at first order: its mean is\n",
          "#' low by 7 percent at 0.8 and 0.6 on a log link.\n#'\n",
          "#' At a row whose `mi()` value left the objective, a reader's fit reads the\n"))

# ---- 19. the nits -------------------------------------------------------------------------
ed("NEWS.md", "for `theta` is how to check a fit that matters. * `sigma` near",
   "for `theta` is how to check a fit that matters.\n\n* `sigma` near")
ed("NEWS.md",
   paste0("  `(Intercept)` at `y\n  * 1000` read flat and took every fitted standard error with it. ",
          "Other fits keep their Hessian and their verdicts bit for bit.\n"),
   paste0("  `(Intercept)` at `y * 1000` read flat and took every fitted standard\n",
          "  error with it. Other fits keep their Hessian and their verdicts bit\n  for bit.\n"))
ed("NEWS.md",
   "Peak\n  memory went down, from 1.1 to 0.72 GB at 500 positions.\n",
   paste0("Peak\n  memory for `frm()` and `summary()` went down, from 0.89 to 0.73 GB at\n",
          "  500 positions. A saved fit grows, from 25 to 36 MB there.\n"))
ub <- "dev/upstream-bugs.md"
ed(ub, "### stan-NN. `beta_binomial_lpmf()`", "### Stan-2. `beta_binomial_lpmf()`")
ed(ub, "- Seen: StanHeaders 2.39.1, 2026-11-02. Number to be settled at consolidation.",
   "- Seen: StanHeaders 2.39.1, 2026-11-02.")
ed(ub, "- Seen: 1.1.0, 2026-11-02. Number to be settled at consolidation.", "- Seen: 1.1.0, 2026-11-02.")
ub_x <- sub("\r$", "", readLines(ub, warn = FALSE))
h <- which(ub_x == "## rtdists")
r1 <- which(startsWith(ub_x, "### rtdists-1."))
if (length(h) == 1L && h < r1 - 2L) {
  stopifnot(ub_x[h + 1L] == "", startsWith(ub_x[h + 2L], "## "))
  ub_x <- c(ub_x[seq_len(r1 - 1L)], "## rtdists", "", ub_x[r1:length(ub_x)])[-c(h, h + 1L)]
  raw <- readBin(ub, "raw", file.size(ub))
  con <- file(ub, "wb")
  writeLines(ub_x, con, sep = if (any(raw == as.raw(13L))) "\r\n" else "\n")
  close(con)
  cat("edited", ub, ": rtdists heading\n")
}
ed("extensions/frmtmb.learn/NEWS.md",
   "The package code is unchanged but for the\n  note on its `frm_compat()` row.",
   "Apart from `rlddm()`'s `type`, the package code\n  is unchanged but for the note on its `frm_compat()` row.")
ed("R/compat.R",
   "Written as brms writes it, bf(..., nl = TRUE, autocor = ~ term) or + acformula(~ term), since a body is evaluated rather than read as terms and a term in the body is refused.",
   "A nonlinear shape or zi beside a linear mu keeps its latent residual, as brms builds it. Written as brms writes it, bf(..., nl = TRUE, autocor = ~ term) or + acformula(~ term), since a body is evaluated rather than read as terms and a term in the body is refused.")
ed("dev/equiv-specs.R",
   paste0("  \"exact gp(): the L1 gap is recorded, not scored; it is unexplained\",\n",
          "  \"(matching the nugget removes about 11 of 94 nats on p094, review\",\n",
          "  \"2026-10-28; dev/equiv-gp-gap.R cannot attribute the rest in double\",\n",
          "  \"precision)\")\n"),
   paste0("  \"exact gp(): the L1 gap is recorded, not scored; it is the nugget, absolute\",\n",
          "  \"1e-12 in brms and 1e-12 * sdgp^2 here: with brms's one line at the relative\",\n",
          "  \"size L1 reads under 1e-8 (release review 2026-11-04)\")\n"))
ed("dev/equiv-specs.R",
   paste0("# by tens to a thousand nats. The nugget explains only part of it\n",
          "# (review 2026-10-28, finding 3), so L1 is recorded and not scored. The\n"),
   "# by tens to a thousand nats. That is the nugget, so L1 is recorded and not scored. The\n")
ed("R/autocor.R",
   "#' `check_se` does not reach it. `frmtmb.sample::frm_sample()` samples the same model and agrees\n",
   "#' `check_se` does not reach it. `frmtmb.sample::frm_sample()` samples the same\n#' model and agrees\n")
ed("R/par-template.R",
   "#' another sign of `b1`, or a located prior, can show it. [frm_allfit()] cannot, as its optimizers\n",
   "#' another sign of `b1`, or a located prior, can show it. [frm_allfit()] cannot,\n#' as its optimizers\n")
ed("R/me.R", "  id\n}\n\n\n#' How one response's", "  id\n}\n\n#' How one response's")
ed("R/importance.R", "      centre <- uhat_g[[j]] - s\n", "      center <- uhat_g[[j]] - s\n")
ed("R/importance.R", "      draw <- centre[gidx] + rinv", "      draw <- center[gidx] + rinv")
ed("R/importance.R", "build the proposal at `at`, centre the weights", "build the proposal at `at`, center the weights")
ed("tests/testthat/test-nl-family.R", "# the first build of this lane, which carried", "# an earlier build, which carried")
ed("tests/testthat/test-nl-family.R", "# build of this lane.", "# an earlier build.")
ed("dev/.gitignore", "*.driverlog\n", "*.driverlog\n*.progress\n*.done\n")
# test-v07.R has no last newline, so it would commit as a line-ending diff
raw <- readBin("tests/testthat/test-v07.R", "raw", file.size("tests/testthat/test-v07.R"))
if (raw[length(raw)] != as.raw(10L)) {
  writeBin(c(raw, as.raw(10L)), "tests/testthat/test-v07.R")
  cat("edited tests/testthat/test-v07.R\n")
}
# backlog: the items the review moves, and what it files
bl <- "dev/test-backlog.md"
ed(bl,
   paste0("- A me() term on another predictor than `mu` beside a latent residual on `mu`:\n",
          "  one node cannot read two independent draws, and the latent's variance does not\n",
          "  join that node (found at this consolidation; `R/me.R`, `me_row_loads()`).\n"),
   paste0("- A `me()` term on another predictor than `mu` beside a latent residual on `mu`:\n",
          "  the `me()` draw on that predictor is not read there. The latent wins, since the\n",
          "  draw moves those means by 0.01 to 0.1 percent (`R/me.R`, `me_row_loads()`).\n",
          "- `conditional_effects()` on a body whose `nlpar` reads `me()` stops, a regression\n",
          "  from 0.78.0: the chain rule through a body is not written.\n",
          "- Above `frmtmb.latres_dense_rows` rows beside another block, `fitted(newdata)` standard\n",
          "  errors are `NA`, with a message. A sparse joint covariance would close it.\n",
          "- `frm_sample()`'s `posterior_epred(newdata, allow_new_levels = TRUE)` refuses a `me(gr =)`\n",
          "  new level (\"cannot predict a grouping level the fit did not see\"), where `fitted()` and brms\n",
          "  answer: 8.051 and 7.979 at one row, seed 7511. It is loud.\n",
          "- A saved fit after `fitted()` at 2000 rows on a latent is 30 MB, about the dense\n",
          "  covariance.\n",
          "- `lognormal`, `me()` on `mu` and `sigma`, sx .8: `Est.Error` 207 against Monte Carlo\n",
          "  3.56 at one row, with `Estimate` fine; not traced. A latent plus one `me()` on two\n",
          "  predictors: `hurdle` `Est.Error` .117 against .139, `lognormal` mean -2.5 percent.\n",
          "- `gp()` in a body: 1.12 nats short of the linear spelling at conv 0 on one of 4\n",
          "  designs. A body's `Est.Error` is `NA` with no message, a gap against brms.\n",
          "- `me(x, sx) * me(z, sz)` on a nonlinear link is read at first order; stated in\n",
          "  `?fitted`.\n"))
ed(bl, "The list in `dev/upstream-bugs.md` is renumbered: brms-45 to brms-50, RTMBdist-5\nand TMB-5 are not reported.",
   "The list in `dev/upstream-bugs.md` is renumbered: brms-45 to brms-50, RTMBdist-5,\nStan-2 and TMB-5 are not reported.")
cat("RELEASE FIXES 2 DONE\n")
# 1 and 13, in words
ed("NEWS.md",
   "  latent residual reads one mean and one variance over the two.\n",
   paste0("  latent residual reads one mean and one variance over the two. Beside\n",
          "  another group-level block it reads the joint covariance, up to\n",
          "  `frmtmb.latres_dense_rows` rows: 14 seconds and 6.4 GB at 20000 rows.\n",
          "  Above that its standard errors are `NA`, with a message of class\n",
          "  `frmtmb_latres_se_skipped`.\n"))
ed("R/autocor.R",
   "#' `options(frmtmb.latres_se_rows = )` moves the limit. Rows in `newdata`,\n#' [frm_curve()] and [conditional_effects()] read `beta` and the term's\n#' parameters alone, and keep their standard errors.\n",
   paste0("#' `options(frmtmb.latres_se_rows = )` moves the limit. Rows in `newdata`,\n",
          "#' [frm_curve()] and [conditional_effects()] keep their standard errors. Beside\n",
          "#' another group-level block, [fitted()] on new rows reads the joint covariance\n",
          "#' up to `frmtmb.latres_dense_rows` rows, and above that reads `NA` with a message.\n"))
cat("RELEASE FIXES 2, WORDS DONE\n")
