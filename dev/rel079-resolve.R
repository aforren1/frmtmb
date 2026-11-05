# Consolidation edits of the 0.79.0 release tree, applied after the
# textual merge (dev/rel079-merge.sh) and before the bin-2 generator and
# the NEWS merge (dev/rel079-news-merge.R). Each edit checks that its old
# text is there exactly once (or the count it names), so a rerun on a
# fresh merge either applies all of them or stops on the first that no
# longer fits. Every hand resolution of the round is here, each with its
# reason; the release fixes are in dev/rel079-release-fix.R, sourced at the
# end. Adapted from dev/rel078-resolve.R.
#
#   Rscript dev/rel079-resolve.R
rel <- Sys.getenv("REL079_ROOT", "C:/Users/adf44/source/r/frmtmb-wt-rel079")
git_dir <- Sys.getenv("REL079_GIT", rel)
Sys.setenv(GIT_OPTIONAL_LOCKS = "0")
setwd(rel)
rd <- function(f) sub("\r$", "", readLines(f, warn = FALSE))
wr <- function(x, f) {
  force(x)
  con <- file(f, "wb")
  writeLines(x, con, sep = "\n", useBytes = TRUE)
  close(con)
}
edit_file <- function(f, old, new, n = 1L) {
  s <- paste(rd(f), collapse = "\n")
  hits <- gregexpr(old, s, fixed = TRUE)[[1]]
  k <- sum(hits > 0)
  if (k != n) stop(f, ": expected ", n, " match, found ", k, ": ",
                   substr(old, 1, 70))
  s <- gsub(old, new, s, fixed = TRUE)
  wr(strsplit(s, "\n", fixed = TRUE)[[1]], f)
  cat("edited", f, "\n")
}
# A rerun starts each file no lane changed from its HEAD blob too, so that every
# edit below applies once however often this runs
chg <- system2("git", c("-C", git_dir, "diff", "--name-only", "HEAD"),
               stdout = TRUE)
for (f in setdiff(chg, readLines("dev/rel079-log/merge-work/all-files"))) {
  h <- system2("git", c("-C", git_dir, "show", paste0("HEAD:", f)), stdout = TRUE)
  wr(sub("\r$", "", h), f)
  cat("from HEAD:", f, "\n")
}
# one conflict block at a time; `pick(ours, theirs, n)` returns its lines
res <- function(f, pick, expect) {
  x <- rd(f)
  out <- character()
  i <- 1L
  n <- 0L
  while (i <= length(x)) {
    if (startsWith(x[i], "<<<<<<< ")) {
      rest <- (i + 1L):length(x)
      j2 <- rest[x[rest] == "======="][1L]
      j3 <- rest[startsWith(x[rest], ">>>>>>> ")][1L]
      ours <- if (j2 > i + 1L) x[(i + 1L):(j2 - 1L)] else character()
      theirs <- if (j3 > j2 + 1L) x[(j2 + 1L):(j3 - 1L)] else character()
      n <- n + 1L
      out <- c(out, pick(ours, theirs, n))
      i <- j3 + 1L
    } else {
      out <- c(out, x[i])
      i <- i + 1L
    }
  }
  if (n != expect) stop(f, ": expected ", expect, " hunks, found ", n)
  wr(out, f)
  cat(f, "hunks:", n, "\n")
}

# ---- 1. the textual conflicts ------------------------------------------
# frm(): the three checks of a latent residual and the one of a read mi()
# value are separate reports on the finished fit
res("R/fit.R", function(ours, theirs, n) c(theirs, ours), 1L)
# smooth's sds (a refusal lane) and latres's sderr each leave class sd
res("R/priors.R", function(ours, theirs, n) {
  if (n == 1L) {
    return(c(
      "    # a gp() term's sd is class \"sdgp\", a smooth's \"sds\" and a latent",
      "    # residual's \"sderr\", which the caller names",
      "    lapply(Filter(function(bk) {",
      "      !bk[[\"covstruct\"]] %in% c(\"gp\", \"hsgp\", \"smooth\", \"latres\")",
      "    }, frame[[\"re_blocks\"]]), function(bk) {"))
  }
  c("        # a gp() term's sd is brms's class \"sdgp\", not \"sd\", a smooth's is",
    "        # \"sds\" (sd_spec_reach() keeps those two apart) and a latent",
    "        # residual's is class \"sderr\"",
    "        if (bk[[\"covstruct\"]] %in% c(\"gp\", \"hsgp\", \"latres\")) next")
}, 2L)
res("extensions/frmtmb.sample/R/sample.R", function(ours, theirs, n) {
  if (n == 1L) {
    return(c("    # so is a latent residual's chain, in unit-scale coordinates",
             "    bk[[\"covstruct\"]] %in% c(\"rr\", \"gp\", \"latres\")"))
  }
  c("    if (identical(bk[[\"covstruct\"]], \"latres\")) next", ours)
}, 2L)
res("tests/testthat/helper-brms.R", function(ours, theirs, n) {
  if (n == 1L) return(c(ours, theirs))
  # the pattern of the inner parameters takes both lanes' spellings, and
  # both blocks of functions stay, the shared closing brace after theirs
  stopifnot(grepl("Ymi(_.+)?|", ours[1L], fixed = TRUE),
            grepl("zerr(_.+)?)$\")", theirs[2L], fixed = TRUE))
  c(ours[1L], theirs[2L], ours[-(1:2)], "}", theirs[-(1:2)])
}, 2L)
# the sampling API's export list: both lanes' names
res("R/sampling-api.R", function(ours, theirs, n) {
  c("#' @rawNamespace export(build_objective, gp_tape_scope, gp_tape_order,",
    "#'   gp_nugget_message, row_lpdf, with_cs_offsets,",
    "#'   cs_offsets_add, update_merge_prior, ce_by_response, resp_refuse_single,",
    "#'   ce_na_vars, autocor_for_newdata,")
}, 1L)
# S5. group:autocor x nl takes nlfix's "conditional"; the mixture row above it
# keeps latres's reading of brms's refusal. The ar() x predict row takes both
# sentences: latres's on a family with no residual, predread's otherwise.
res("R/compat.R", function(ours, theirs, n) {
  if (n == 1L) return(c(ours[1L], theirs[-1L]))
  stopifnot(grepl("Under cov = TRUE predict() draws", theirs, fixed = TRUE))
  sub("\"Under cov = TRUE predict() draws",
      paste0("\"On a family with no residual y - mu, predict() reads the ",
             "latent's conditional mode on fitted rows and draws new latents ",
             "on new rows under cov = TRUE. On gaussian and student, under ",
             "cov = TRUE predict() draws"),
      theirs, fixed = TRUE)
}, 2L)
# S4b: one mean. A row's me() draw and latent residual are one node.
res("R/interop.R", function(ours, theirs, n) {
  c("  # the expected response over a new row's me() value and latent residual,",
    "  # as fitted() and brms's posterior_epred() read it",
    "  if (type == \"response\") {",
    "    mm <- fitted_me_moments(model, newdata, dots[[\"re_formula\"]],",
    "                            dots[[\"resp\"]], dpar,",
    "                            isTRUE(dots[[\"allow_new_levels\"]]))",
    "    if (!is.null(mm)) {",
    "      attributes(mm$mean) <- attributes(p)",
    "      p <- mm$mean",
    "    } else {",
    "      p <- latres_fitted_mean(model, newdata, dots[[\"re_formula\"]],",
    "                              dots[[\"resp\"]], dpar,",
    "                              isTRUE(dots[[\"allow_new_levels\"]]), p)",
    "    }")
}, 1L)
# fitted_point_se(): latres's early returns, on no me() row, then predread's test
res("R/predict.R", function(ours, theirs, n) {
  c(ours[-length(ours)], theirs)
}, 1L)
# latres's shift, then predread's cov = TRUE frame on the shifted fit
res("extensions/frmtmb.sample/R/methods-draws.R", function(ours, theirs, n) {
  c(ours, theirs[-1L])
}, 1L)
# dev/upstream-bugs.md is renumbered below: both entries stay
res("dev/upstream-bugs.md", function(ours, theirs, n) c(ours, "", theirs), 1L)

# ---- 2. gpgap against ncgp ---------------------------------------------
# dev/rel079-merge.sh took the listed files whole from ncgp. What gpgap kept:
# brms_gp_order() (helpers, merged three-way), the workflows' step, the
# ?frmtmb-gp page, rewritten here for one tape.
page <- c(
  "#' Exact Gaussian processes and their nugget",
  "#'",
  "#' @description",
  "#' `gp(x)` in a formula fits an exact Gaussian process, and `gp(x, k = )` its",
  "#' Hilbert-space approximation, with the arguments of [brms::gp()]. This page",
  "#' states where the exact form is not brms's model.",
  "#'",
  "#' @section The nugget:",
  "#' A squared exponential kernel is near singular, so a nugget on its diagonal",
  "#' keeps it factorable. `frm()` and `frm_sample()` read one tape, the field",
  "#' `sdgp * L z`, with `L` the Cholesky factor of the kernel plus `1e-12` on the",
  "#' correlation: a nugget of `1e-12 * sdgp^2`. brms adds an absolute `1e-12`. At",
  "#' `sdgp` near 1 the two are one number. A relative nugget does not depend on",
  "#' the units of y. brms's does: at y times 1000 its kernel does not factor at",
  "#' the estimates (its run reads Rhat 3.77), and at y divided by 1000 its nugget",
  "#' floors `sigma`.",
  "#'",
  "#' The two nuggets give one model until the noise reaches the nugget, near",
  "#' `sigma = 1e-6 * sdgp`:",
  "#'",
  "#' * With one row per position and an estimated `sigma`, `sigma` can read low.",
  "#' * Where positions repeat, or `sigma` is known through `se()`, or the family",
  "#'   has no `sigma`, `sdgp` and the length scale can read low.",
  "#'",
  "#' @section The warning:",
  "#' `frm()` warns with class `frmtmb_gp_nugget`, on top of `frmtmb_warning`, when",
  "#' the 90th percentile over the positions of the nugget's variance over the",
  "#' noise at a position reaches `sqrt(2 / m)` for `m` positions, which is about",
  "#' one standard error of `sigma`. The noise at a position is read from the",
  "#' likelihood's curvature, so the rule is the same for every family and route.",
  "#' `frm()` does not check an `importance =` or `integrate =` fit.",
  "#' `frm_sample(fit)` on a fit with this warning starts at random, because chains",
  "#' started at that mode stay pinned.",
  "#'",
  "#' @section What to do about it:",
  "#' `frm_sample(laplace = TRUE)` integrates the field out, and it mixes fast at",
  "#' low noise, for the outer parameters only (`sigma`, `sdgp`, the length",
  "#' scales): it refuses `posterior_epred()`, `posterior_linpred()`,",
  "#' `posterior_predict()`, `fitted()`, `predict()`, `log_lik()`, `loo()` and",
  "#' `conditional_effects()`. The drawn route, the default, is the one for the",
  "#' field and for predictions, and at low noise it needs long chains (at sigma",
  "#' 1e-3 and 15 rows, Rhat on `sdgp` read 1.04 to 2.37 on short chains, and 1.03",
  "#' on `sigma` and 1.07 on `sdgp` after an hour, and not yet converged).",
  "#' `gp(x, k = )` has no nugget and reaches brms's numbers when `k` is large",
  "#' enough for the field. Too small a `k` is a different failure, with no",
  "#' warning and with `convergence` 0: `sin(8 * x)` on 0 to 6 with noise 1e-3 at",
  "#' `k = 30` has a basis that stops",
  "#' at frequency 6.5, so it reads `sigma` 720 times high; `k = 60` reads it.",
  "#' Compare several `k`. A smooth trend written as fixed effects also lowers what",
  "#' the field has to carry.",
  "#'",
  "#' @section The order of the field:",
  "#' `frm_sample()`'s `b[i]` columns and `fit$estimates$b` of an exact term hold",
  "#' z, in a farthest-first order of the positions. [ranef()], [coef()] and",
  "#' `getME(fit, \"b\")` read the field, at the sorted distinct positions, where",
  "#' brms keeps first appearance. None of these is in the order of the rows of",
  "#' the data. Read a field at the data through `posterior_epred()` or",
  "#' `frm_linpred()`, which are right row by row.",
  "#'",
  "#' @name frmtmb-gp",
  "#' @aliases frmtmb_gp_nugget",
  "#' @seealso [frm()], [frmtmb_control()], `vignette(\"brms-migration\")`",
  "NULL")
x <- rd("R/gp-by.R")
a <- which(x == "#' Exact Gaussian processes and their nugget")
b <- a - 1L + which(x[a:length(x)] == "NULL")[1L]
stopifnot(length(a) == 1L, x[b + 1L] == "" || b == length(x))
wr(c(x[seq_len(a - 1L)], page, x[-seq_len(b)]), "R/gp-by.R")
cat("rewrote the page frmtmb-gp in R/gp-by.R\n")
# the page's sentence in compatibility.Rmd is not the nugget route, so it stays
edit_file("vignettes/compatibility.Rmd",
          "`gp_pred()`, and `cs_pred()` to keep them apart.",
          paste0("`gp_pred()`, and `cs_pred()` to keep them apart. A row that ",
                 "marks the exact\n`gp()` as working means its model apart from ",
                 "the nugget, which `?\"frmtmb-gp\"`\nstates."))
edit_file("R/conditions.R",
          paste0("the nugget of an exact `gp()` term, and tells a boundary fit ",
                 "with class"),
          paste0("the nugget of an exact `gp()` term (see [`frmtmb-gp`]), and ",
                 "tells a boundary fit\n#' with class"))
# test-brms-lik-gp-exact.R: brms's own program now, at most 1e-12 apart
tf <- "tests/testthat/test-brms-lik-gp-exact.R"
x <- rd(tf)
a <- grep("^gpx_frm_nugget <- function", x)
b <- a - 1L + which(x[a:length(x)] == "}")[1L]
stopifnot(length(a) == 1L)
# its comment block sits right above it
ca <- a
while (ca > 1L && startsWith(x[ca - 1L], "#")) ca <- ca - 1L
cat("dropped from the test file:\n", paste0("  ", x[ca:b], "\n"), sep = "")
x <- x[-(ca:b)]
s <- paste(x, collapse = "\n")
s <- gsub("joint = TRUE,\n                  code_fix = gpx_frm_nugget)",
          "joint = TRUE)", s, fixed = TRUE)
s <- gsub("joint = TRUE,\n                code_fix = gpx_frm_nugget)",
          "joint = TRUE,\n                code_fix = gpx_rel_nugget)", s, fixed = TRUE)
s <- sub("exact gp() on unsorted positions, at frmtmb's nugget",
         "exact gp() on unsorted positions", s, fixed = TRUE)
stopifnot(!grepl("gpx_frm_nugget", s, fixed = TRUE))
x <- strsplit(s, "\n", fixed = TRUE)[[1]]
a <- grep("^test_that\\(\"exact gp\\(\\) under gr = FALSE keeps row order", x)
stopifnot(length(a) == 1L)
x <- c(x[seq_len(a - 1L)],
  "# The row below reads the order, on a kernel where the two nuggets' sizes differ by",
  "# 0.026 nats and 2.3 in the gradient. brms's one line is read at frmtmb's relative",
  "# size for that row alone; brms's own line is read in the other rows here and in",
  "# test-brms-likelihood.R.",
  "gpx_rel_nugget <- function(code) {",
  "  new <- gsub(\"cov[n, n] += 1e-12;\", \"cov[n, n] += 1e-12 * square(sdgp);\",",
  "              code, fixed = TRUE)",
  "  testthat::expect_false(identical(new, code))",
  "  new",
  "}",
  "",
  x[a:length(x)])
a <- grep("^test_that\\(\"the sampler's factor at its nugget", x)
stopifnot(length(a) == 1L)
nx <- which(startsWith(x, "test_that(") & seq_along(x) > a)
b <- if (length(nx)) nx[1L] - 1L else length(x)
x <- x[-(a:b)]
while (x[length(x)] == "") x <- x[-length(x)]
wr(x, tf)
cat("edited", tf, "\n")

# ---- 3. refuse x ncgp: an mi() retape of an exact gp() --------------------
# Without the scope the tape can hit TMB-5, a wrong gradient. Every other
# MakeADFun is listed in dev/round-20261104.md, "MakeADFun".
edit_file("R/frame.R",
          paste0("  fit$obj <- RTMB::MakeADFun(build_objective(fit$frame), fit$estimates,\n",
                 "                             random = random, map = fit$frame[[\"map\"]],\n",
                 "                             silent = TRUE)"),
          paste0("  fit$obj <- gp_tape_scope(fit$frame, RTMB::MakeADFun(\n",
                 "    build_objective(fit$frame), gp_tape_order(fit$frame, fit$estimates),\n",
                 "    random = random, map = fit$frame[[\"map\"]], silent = TRUE))"))

# ---- 4. predread + latres + nlfix ---------------------------------------
# S2: a dpar's own link in a body branch (a dpar with a link of its own reads
# wrong through the mu link)
edit_file("R/predict.R",
          "out <- if (type == \"response\") lp[[\"link\"]]$linkinv(eta) else eta\n      return(if (is.null(newdata)) napred(object, out) else out)",
          "out <- if (type == \"response\") {\n        dpar_response_link(rspec, lp, dpar)$linkinv(eta)\n      } else {\n        eta\n      }\n      return(if (is.null(newdata)) napred(object, out) else out)")
# S4a: one response_mean_gh(), latres's, which guards the tails. R/me.R
# collates after R/autocor-latent.R and would win silently.
x <- rd("R/me.R")
st <- which(startsWith(x, "response_mean_gh <- function"))
stopifnot(length(st) == 1L)
en <- st + which(x[(st + 1L):length(x)] == "}")[1L]
ca <- st
while (ca > 1L && startsWith(x[ca - 1L], "#'")) ca <- ca - 1L
cat("R/me.R: response_mean_gh removed, lines", ca, "to", en, "\n")
x <- x[-(ca:en)]
wr(x, "R/me.R")
# S4b and S4c, M3's rule on a latent residual: it is a draw of a new row as a
# me() value is, so the mean, Est.Error and Q all read the node. A latent
# residual alone gets a load too (S4c); beside me() its variance joins the node
# the me() call has made (S4b).
edit_file("R/me.R",
          "me_row_loads <- function(fit, rspec, newdata, use_re, anl, dpar = NULL) {\n  if (is.null(fit$frame[[\"me\"]]) || is.null(newdata)) return(NULL)\n  rn <- rspec$resp_name",
          "me_row_loads <- function(fit, rspec, newdata, use_re, anl, dpar = NULL) {\n  if (is.null(newdata)) return(NULL)\n  rn <- rspec$resp_name\n  sl <- latres_row_sd(fit, rn, dpar)\n  if (is.null(fit$frame[[\"me\"]])) {\n    if (is.null(sl)) return(NULL)\n    key <- linpred_key(rn, \"mu\")\n    s <- rep_len(sl, nrow(newdata))\n    return(list(key = key, s2 = s^2, by_key = stats::setNames(list(s), key)))\n  }")
edit_file("R/me.R",
          "  if (!length(parts)) return(NULL)\n  tot <- function(mp) Reduce(`+`, lapply(mp, function(e) (e$g * e$sd)^2))",
          "  if (!length(parts)) {\n    if (is.null(sl)) return(NULL)\n    key <- linpred_key(rn, \"mu\")\n    s <- rep_len(sl, nrow(newdata))\n    return(list(key = key, s2 = s^2, by_key = stats::setNames(list(s), key)))\n  }\n  tot <- function(mp) Reduce(`+`, lapply(mp, function(e) (e$g * e$sd)^2))")
edit_file("R/me.R",
          "  by_key[[names(parts)[1L]]] <- s\n  list(key = names(parts)[1L], s2 = s2, by_key = by_key)\n}",
          paste0(
            "  by_key[[names(parts)[1L]]] <- s\n",
            "  # a latent residual is independent of the me() draws and moves mu alone, so\n",
            "  # the lead's variance grows by it and the other predictors' share of the\n",
            "  # node falls\n",
            "  if (!is.null(sl) && identical(names(parts)[1L], linpred_key(rn, \"mu\"))) {\n",
            "    s_new <- sqrt(s2 + sl^2)\n",
            "    for (k in names(by_key)) {\n",
            "      by_key[[k]] <- ifelse(s_new > 0, by_key[[k]] * s / s_new, 0)\n",
            "    }\n",
            "    by_key[[names(parts)[1L]]] <- s_new\n",
            "    s2 <- s_new^2\n",
            "  }\n",
            "  list(key = names(parts)[1L], s2 = s2, by_key = by_key)\n}\n\n",
            "#' The marginal sd a latent residual adds to `mu` of a new row, or `NULL` where\n",
            "#' the read is not of `mu` or the fit has none.\n",
            "#'\n",
            "#' @noRd\n",
            "latres_row_sd <- function(fit, rn, dpar) {\n",
            "  if (!is.null(dpar) && !identical(dpar, \"mu\")) return(NULL)\n",
            "  bk <- frame_latres(fit$frame)[[rn]]\n",
            "  if (is.null(bk) || isTRUE(fit$frame[[\"latres_off\"]])) return(NULL)\n",
            "  latres_marginal_sd(fit$estimates[[\"theta\"]][bk[[\"theta_idx\"]]],\n",
            "                     bk[[\"aux\"]])\n",
            "}"))
edit_file("R/predict.R",
          "  if (is.null(newdata) || is.null(object$frame[[\"me\"]])) return(NULL)\n  resp <- resp_read(object, resp)",
          "  if (is.null(newdata)) return(NULL)\n  if (is.null(object$frame[[\"me\"]]) && !length(frame_latres(object$frame))) {\n    return(NULL)\n  }\n  resp <- resp_read(object, resp)")
edit_file("R/predict.R",
          "  if (!is.null(newdata) && scale == \"response\" &&\n        (is.null(dpar) || identical(dpar, \"mu\"))) {\n    out <- latres_fitted_mean(",
          "  if (is.null(mm) && !is.null(newdata) && scale == \"response\" &&\n        (is.null(dpar) || identical(dpar, \"mu\"))) {\n    out <- latres_fitted_mean(")

# S4c holds the mean's differences: frm_linpred(se.fit = TRUE) reads the joint
# covariance, 3 GB at 12000 rows where the differences take 0.47 GB, and above the
# dense limit it refuses, which the handler below would read as no standard error
edit_file("R/predict.R",
          paste0("    return(as.vector(fit_fd_se(object, function(fit) {\n",
                 "      fitted_point(fit, newdata, re_formula, scale, resp, dpar,\n",
                 "                   allow_new_levels)\n",
                 "    })))\n"),
          paste0("    fd <- as.vector(fit_fd_se(object, function(fit) {\n",
                 "      fitted_point(fit, newdata, re_formula, scale, resp, dpar,\n",
                 "                   allow_new_levels)\n",
                 "    }))\n",
                 "    # the draw's variance comes from the integral, beside the coefficients'\n",
                 "    return(if (is.null(mm)) fd else sqrt(fd^2 + mm$var))\n"))
# get_predict() reads the mean alone, and marginaleffects calls it once per step:
# the variance and the slope are 2 more integrals a call
edit_file("R/predict.R",
          "fitted_me_moments <- function(object, newdata, re_formula, resp, dpar,\n                              allow_new_levels) {",
          "fitted_me_moments <- function(object, newdata, re_formula, resp, dpar,\n                              allow_new_levels, mean_only = FALSE) {")
edit_file("R/predict.R",
          "  m1 <- over(val)\n  m2 <- over(function(z) val(z)^2)\n",
          "  m1 <- over(val)\n  if (mean_only) return(list(mean = drop(m1)))\n  m2 <- over(function(z) val(z)^2)\n")
edit_file("R/interop.R",
          "                            isTRUE(dots[[\"allow_new_levels\"]]))\n    if (!is.null(mm)) {",
          "                            isTRUE(dots[[\"allow_new_levels\"]]),\n                            mean_only = TRUE)\n    if (!is.null(mm)) {")

# ---- 5. predread x ncgp: the second maximum's numbers -------------------
# "0.46" was 0.78.0's second maximum (dev/reviews/2026-11-03-ncgp.md); the
# numbers are remeasured on the merged tree by dev/rel079-mfx-remeasure.R
edit_file("R/interop.R",
          "On a `gp()` the standard error is 0.46 on every row, where\n#'   [frm_linpred()] reads 0.089 to 0.880;",
          "On a `gp()` the standard error is 0.31 on every row, where\n#'   [frm_linpred()] reads 0.087 to 0.695;")
edit_file("R/interop.R",
          "0.46 on every row where `frm_linpred(se.fit = TRUE)` reads 0.089 to 0.880.",
          "0.31 on every row where `frm_linpred(se.fit = TRUE)` reads 0.087 to 0.695.")

# ---- 6. dev/upstream-bugs.md: one number per entry -----------------------
# latres keeps brms-45, which test-latres-brms.R and the decisions name.
# refuse's three follow it, then nlfix's and predread's brms-NN, and nlfix's
# RTMBdist-NN. The refuse entries move into the brms section. Nothing is posted.
ub <- "dev/upstream-bugs.md"
edit_file(ub, "### brms-NN. `predictor.bprepnl()`",
          "### brms-49. `predictor.bprepnl()`")
edit_file(ub, "- Seen: 2.23.0. Number to be settled at consolidation. Found 2026-11-02.",
          "- Seen: 2.23.0. Found 2026-11-02.")
edit_file(ub, "### brms-NN. `me()` on new data", "### brms-50. `me()` on new data")
edit_file(ub, "### RTMBdist-NN.", "### RTMBdist-5.")
x <- rd(ub)
heads <- grep("^### ", x)
a <- grep("^### brms-45 \\(refuse\\)", x)
b <- heads[heads > grep("^### brms-47 \\(refuse\\)", x)][1L] - 1L
blk <- x[a:b]
blk <- sub("^### brms-47 \\(refuse\\)", "### brms-48", blk)
blk <- sub("^### brms-46 \\(refuse\\)", "### brms-47", blk)
blk <- sub("^### brms-45 \\(refuse\\)", "### brms-46", blk)
x <- x[-(a:b)]
# after the last brms entry
heads <- grep("^### ", x)
last <- max(grep("^### brms-", x))
at <- heads[heads > last][1L] - 1L
sect <- grep("^## ", x)
if (any(sect > last & sect < at)) at <- sect[sect > last][1L] - 1L
x <- c(x[seq_len(at)], blk, x[-seq_len(at)])
wr(x, ub)
nn <- grep("brms-NN|RTMBdist-NN|\\(refuse\\)", x, value = TRUE)
if (length(nn)) stop("left in the upstream list: ", nn[1L])
ids <- sub("^### ([^. ]+).*", "\\1", grep("^### ", x, value = TRUE))
if (anyDuplicated(ids)) stop("duplicate entry: ", ids[duplicated(ids)][1L])
cat("upstream-bugs: ", length(ids), "entries, no number twice\n")

# ---- 7. versions and floors ---------------------------------------------
# dev/round-20261104.md, "Versions and floors"; the floors are measured by
# dev/rel079-floor.sh. Core: minor, breaking. frmtmb.sample calls core's new
# sampling-API exports. frmtmb.learn: its documented numbers and its check need
# the new core. eam and survival gain rows only, and call nothing new.
ver <- function(f, from, to) edit_file(f, paste0("Version: ", from),
                                       paste0("Version: ", to))
ver("DESCRIPTION", "0.78.0", "0.79.0")
edit_file("codemeta.json", '"version": "0.78.0",', '"version": "0.79.0",')
dsa <- "extensions/frmtmb.sample/DESCRIPTION"
ver(dsa, "0.25.0", "0.26.0")
edit_file(dsa, "frmtmb (>= 0.78.0)", "frmtmb (>= 0.79.0)")
dle <- "extensions/frmtmb.learn/DESCRIPTION"
ver(dle, "0.7.3", "0.7.4")
edit_file(dle, "frmtmb (>= 0.68.0)", "frmtmb (>= 0.79.0)")
ver("extensions/frmtmb.eam/DESCRIPTION", "0.13.0", "0.13.1")
ver("extensions/frmtmb.survival/DESCRIPTION", "0.2.0", "0.2.1")

# ---- 8. the release harness reads rellib-r17 -----------------------------
# a file no lane changed comes back from its HEAD blob, so a rerun applies once
merged <- readLines("dev/rel079-log/merge-work/files")
for (f in c("dev/brmsport-bin2-run.R", "dev/brmsport-guards.R",
            "dev/brmsport-run.R", "dev/ciharden-run1.R",
            "dev/release/run-check.ps1", "dev/release/run-docs.ps1",
            "dev/release/run-tests.R", "dev/warnleak-scan.R",
            "dev/release/run-session.R", "dev/brmsport-bin2-all.sh",
            "dev/brmsport-bin2-record.sh", "dev/ciharden-scan.sh")) {
  if (!f %in% merged) {
    h <- system2("git", c("-C", git_dir, "show", paste0("HEAD:", f)),
                 stdout = TRUE)
    wr(sub("\r$", "", h), f)
  }
  s <- paste(rd(f), collapse = "\n")
  k <- lengths(regmatches(s, gregexpr("rellib-r16", s, fixed = TRUE)))
  stopifnot(k >= 1L)
  wr(strsplit(gsub("rellib-r16", "rellib-r17", s, fixed = TRUE), "\n",
              fixed = TRUE)[[1]], f)
  cat("edited", f, "(", k, ")\n")
}

source("dev/rel079-adjudicate.R", local = FALSE)
source("dev/rel079-release-fix.R", local = FALSE)
cat("RESOLVE DONE\n")
