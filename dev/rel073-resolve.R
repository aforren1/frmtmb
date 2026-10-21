# Consolidation edits of the 0.73.0 release tree, applied after the
# textual merge (dev/rel073-merge.sh) and before the NEWS merge
# (dev/rel073-news-merge.R). Each edit checks that its old text is
# there exactly once, so a rerun on a fresh merge either applies all of
# them or stops on the first that no longer fits. A file no lane
# changed is restored from its HEAD blob (4854e425, 0.72.0) before it
# is edited, so every edit here is applied once however often this
# runs. Adapted from dev/rel072-resolve.R.
#
#   Rscript dev/rel073-resolve.R
rel <- Sys.getenv("REL073_ROOT", "C:/Users/adf44/source/r/frmtmb-wt-release")
# the checkout whose blobs are read (dev/rel073-merge.sh, REL073_GIT)
git_dir <- Sys.getenv("REL073_GIT", rel)
setwd(rel)
rd <- function(f) sub("\r$", "", readLines(f, warn = FALSE))
wr <- function(x, f) {
  force(x)
  con <- file(f, "wb")
  writeLines(x, con, sep = "\n")
  close(con)
}
from_head <- function(f) {
  x <- system2("git", c("-C", git_dir, "show", paste0("HEAD:", f)),
               stdout = TRUE)
  wr(sub("\r$", "", x), f)
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
merged <- readLines("dev/rel073-log/merge-work/files")
restore <- function(f) {
  if (f %in% merged) stop(f, " is a merged file; edit it, do not restore it")
  from_head(f)
}
copied <- function(lane) {
  readLines(file.path("dev/rel073-log", paste0("copied-", lane, ".lst")))
}

# ---- tests/testthat/helper-fuzz.R: cdfs against fuzzfix ----------------
# Both add one entry at the end of FUZZ_KNOWN_REFUSAL: cdfs's
# trunc-no-mean (fitted()'s refusal of a truncated mean brms does not
# have) and fuzzfix's quad-one-block-per-row (the quadrature refusal of
# a nested block). Different refusals on different specs: keep both
edit_file("tests/testthat/helper-fuzz.R",
'<<<<<<< release
  list(id = "trunc-no-mean",',
'  list(id = "trunc-no-mean",')
edit_file("tests/testthat/helper-fuzz.R",
'         grepl("has no truncated mean", f$detail))
=======
  list(id = "quad-one-block-per-row",',
'         grepl("has no truncated mean", f$detail)),
  list(id = "quad-one-block-per-row",')
edit_file("tests/testthat/helper-fuzz.R",
'               f$detail, fixed = TRUE))
>>>>>>> fuzzfix
)',
'               f$detail, fixed = TRUE))
)')

# ---- the upstream-bug numbers (brief item 2) ----------------------------
# Three lanes appended a "brms-30". In the brief's lane order: loodraws's
# ndraws-subsample entry stays brms-30; cdfs's four, in their order,
# become brms-31 to brms-34; nb2dm's rdirichlet entry becomes brms-35.
# One pass per text, so that 30 -> 31 is not then read as 31 -> 32
renum <- function(s, map) {
  m <- gregexpr("brms-[0-9]+", s)
  regmatches(s, m) <- lapply(regmatches(s, m), function(v) {
    k <- v %in% names(map)
    v[k] <- map[v[k]]
    v
  })
  s
}
map_cdfs <- c("brms-30" = "brms-31", "brms-31" = "brms-32",
              "brms-32" = "brms-33", "brms-33" = "brms-34")
map_nb2dm <- c("brms-30" = "brms-35")

ub <- "dev/upstream-bugs.md"
s <- paste(rd(ub), collapse = "\n")
# the brms block: release(loodraws) / cdfs, then that / nb2dm; each
# side ends where the shared "- Status: not reported." line follows
at <- function(pat, from = 1L) {
  p <- regexpr(pat, substring(s, from), fixed = TRUE)
  stopifnot(p > 0)
  c(from + p - 1L, from + p - 1L + attr(p, "match.length"))
}
pa <- at("<<<<<<< release\n<<<<<<< release\n### brms-30.")
pa[2] <- pa[1] + nchar("<<<<<<< release\n<<<<<<< release\n")
pb <- at("\n=======\n", pa[2])
pc <- at("\n>>>>>>> cdfs\n=======\n", pb[2])
pd <- at("\n>>>>>>> nb2dm\n- Status: not reported.\n", pc[2])
loo_part <- substr(s, pa[2], pb[1] - 1L)
cdfs_part <- renum(substr(s, pb[2], pc[1] - 1L), map_cdfs)
nb_part <- renum(substr(s, pc[2], pd[1] - 1L), map_nb2dm)
stopifnot(startsWith(loo_part, "### brms-30. `loo_epred(ndraws =)`"),
          startsWith(cdfs_part, "### brms-31. `zero_inflated_asym_laplace"),
          startsWith(nb_part, "### brms-35. `rdirichlet()`"))
st <- "\n- Status: not reported.\n"
s <- paste0(substr(s, 1L, pa[1] - 1L), loo_part, st, "\n", cdfs_part, st,
            "\n", nb_part, st, substring(s, pd[2]))
# rstan-3 (cdfs) and Stan-1 (nb2dm): both, each with its status line
s0 <- s
s <- sub(paste0("(?s)<<<<<<< release\n(### rstan-3[^\n]*\n.*?)=======\n",
                "(### Stan-1[^\n]*\n.*?)>>>>>>> nb2dm\n"),
         "\\1- Status: not reported.\n\n\\2", s, perl = TRUE)
stopifnot(!identical(s, s0))
stopifnot(!grepl("^(<<<<<<<|=======|>>>>>>>)", strsplit(s, "\n")[[1]]))
wr(strsplit(s, "\n", fixed = TRUE)[[1]], ub)
# a prefix can hold a digit (EMC2-1): the release review, finding 5
ids <- regmatches(s, gregexpr("(?m)^### [A-Za-z0-9]+-[0-9]+", s,
                              perl = TRUE))[[1]]
stopifnot(!anyDuplicated(ids))
cat("upstream-bugs.md: brms entries renumbered;", length(ids),
    "IDs, all unique\n")

# Every reference in the files the two lanes brought. The review of lane
# cdfs discusses the numbering itself in two passages, which keep the
# numbers they discussed and gain a note of the renumbering
renum_file <- function(f, map, keep = NULL) {
  x <- rd(f)
  skip <- logical(length(x))
  for (k in keep) {
    a <- which(startsWith(x, k[1]))
    stopifnot(length(a) == 1L)
    z <- a - 1L + which(grepl(k[2], x[a:length(x)], fixed = TRUE))[1]
    skip[a:z] <- TRUE
  }
  y <- x
  y[!skip] <- renum(x[!skip], map)
  if (!identical(x, y)) {
    wr(y, f)
    cat("renumbered", sum(x != y), "lines of", f, "\n")
  }
}
cdfs_files <- c(copied("cdfs"), "R/objective.R",
                "tests/testthat/test-cens-trunc.R")
cdfs_files <- cdfs_files[grepl("[.](R|md|sh|jobs)$", cdfs_files)]
rv <- "dev/reviews/2026-10-20-cdfs.md"
for (f in setdiff(cdfs_files, rv)) renum_file(f, map_cdfs)
renum_file(rv, map_cdfs, keep = list(
  c("2. **Upstream numbering.** Three lanes appended", "do not collide."),
  c("- **Upstream numbering.** The renumbering of the", "brms-33.")))
renote <- function(ind) paste0(
  ind, "[At the 0.73.0 consolidation, in the brief's lane order:\n",
  ind, "loodraws brms-30; cdfs's four brms-31 to brms-34, in their\n",
  ind, "order; nb2dm brms-35. The other references in this file use the\n",
  ind, "new numbers.]")
edit_file(rv,
  "   RTMBdist-3 and rstan-3 are this lane's alone and do not collide.",
  paste0("   RTMBdist-3 and rstan-3 are this lane's alone and do not ",
         "collide.\n", renote("   ")))
edit_file(rv,
  "  entries is still for consolidation. This lane now has brms-30 to
  brms-33.",
  paste0("  entries is still for consolidation. This lane now has ",
         "brms-30 to\n  brms-33.\n", renote("  ")))
for (f in c(copied("nb2dm"))) {
  if (grepl("[.](R|md|sh)$", f)) renum_file(f, map_nb2dm)
}

# ---- negbinomial2 in cdfs's CDF tables (brief item 1) -------------------
# nb2dm's negbinomial2 is negbinomial() underneath (fam_reparam()), so it
# inherits the CDF, cens(), trunc() and the truncated mean; the compat
# table names the families one by one, and brms writes neg_binomial_2
# for both (dev/reviews/2026-10-20-cdfs.md, consolidation note 1)
edit_file("R/compat.R",
'compat_cdf_discrete <- c("binomial", "beta_binomial", "negbinomial",
                         "geometric", "compois", "zero_inflated_poisson",',
'compat_cdf_discrete <- c("binomial", "beta_binomial", "negbinomial",
                         "negbinomial2", "geometric", "compois",
                         "zero_inflated_poisson",')
edit_file("R/compat.R",
'    stan <- switch(f, negbinomial = , geometric = "neg_binomial_2",',
'    stan <- switch(f, negbinomial = , negbinomial2 = ,
                   geometric = "neg_binomial_2",')
edit_file("R/compat.R",
'    tmean <- f %in% c("binomial", "negbinomial", "geometric")',
'    tmean <- f %in% c("binomial", "negbinomial", "negbinomial2",
                      "geometric")')

# ---- loodraws's last NITs (brief item 5) --------------------------------
# the review's check after fix round 2: the factor S belongs on the
# rounding term only, and log_lik() sums nothing over draws
tl <- "extensions/frmtmb.sample/tests/testthat/test-loo-predict.R"
edit_file(tl,
'  # relative difference of the prediction matrix. A leave-one-out
  # summary of a row is a sum over the S draws: each term\'s weight is
  # exp() of a log-likelihood, so its relative error is at most 2 * dl
  # once normalized, its value\'s at most dx, and the sum adds eps per
  # term. So a weighted mean errs by at most S * (2 dl + dx + eps) of
  # the largest value. The old r_eff = 1 rule on a subset moved',
'  # relative difference of the prediction matrix. A leave-one-out
  # summary of a row is a self-normalized weighted sum over the S
  # draws: each weight is exp() of a log-likelihood, so its relative
  # error is at most 2 * dl once normalized, and each value\'s is at
  # most dx. Those errors do not add up over the draws: a weighted mean
  # errs by about (2 dl + dx) of the largest value. Only the rounding
  # of the sum grows with S, by eps per term. The bounds below multiply
  # every term by S, so the dl and dx parts carry a headroom of S (1000
  # here). The old r_eff = 1 rule on a subset moved')
edit_file(tl,
'    # the inputs themselves agree to rounding: a guard on the guard
    expect_lte(dl, S * eps * max(abs(llb)))',
'    # the inputs themselves agree to rounding: a guard on the guard. A
    # log_lik() value sums nothing over the draws, so S is headroom over
    # its few eps of rounding, not a count of terms
    expect_lte(dl, S * eps * max(abs(llb)))')
edit_file(tl, "the slope between neighbouring", "the slope between neighboring")
# loo_moment_match() on a fit: moment matching transforms the draws and
# needs the sampler's log density; refitting is only brms's fallback.
# The helper stays one frm_stop() statement (test-arg-refusal.R)
lg <- "R/loo-generics.R"
edit_file(lg,
"#' The refusal of a refitting function on a fit. Sampling would not
#' help: frmtmb.sample refuses it on draws too. `loo_moment_match()`
#' is one: it moment-matches the draws toward each fold and refits the
#' folds that stay bad.
#'
#' @noRd
loo_refit_refusal <- function(what) {
  frm_stop(what, \" refits the sampler on modified data, and neither \",",
"#' The refusal of a function that needs the sampler itself, on a fit.
#' Sampling would not help: frmtmb.sample refuses it on draws too.
#' `how` says what the function does with the sampler: `kfold()` and
#' `reloo()` refit it, and `loo_moment_match()` transforms the draws
#' with its log density and refits only where that fails. `join` ends
#' that clause: a long one ends its sentence (the release review,
#' finding 7).
#'
#' @noRd
loo_refit_refusal <- function(what,
                              how = \" refits the sampler on modified data\",
                              join = \", and\") {
  frm_stop(what, how, join, \" neither \",")
edit_file(lg,
"loo_moment_match.frmtmb_fit <- function(x, ...) {
  loo_refit_refusal(\"loo_moment_match()\")
}",
"loo_moment_match.frmtmb_fit <- function(x, ...) {
  loo_refit_refusal(\"loo_moment_match()\", paste0(
    \" transforms the posterior draws toward each observation with a \",
    \"high Pareto k, which needs the sampler's log density of every \",
    \"draw, and it refits the sampler where that fails\"), \". Here\")
}")

# ---- the release review (dev/reviews/2026-10-21-release.md) -----------
# Finding 1: check_weights_unbounded() returned early on cens(), and
# lane cdfs made eight spike families censorable, so a censoring column
# with no row censored switched both refusals off. The check now runs
# on the uncensored rows, subtracts 1/2 per censored row under Laplace,
# and carries zero_inflated_asym_laplace(), whose runaway was measured
# (dev/rel073-wcens-probe.R, dev/rel073-zial-probe.R). The code and the
# tests are dev/rel073-review.patch against lane fuzzfix's files, which
# the merge copies untracked
rp <- system2("git", c("apply", "--whitespace=nowarn",
                         "dev/rel073-review.patch"), stdout = TRUE,
              stderr = TRUE)
if (!is.null(attr(rp, "status"))) stop("rel073-review.patch: ", rp)
# git apply writes CRLF under core.autocrlf; the tree keeps LF
wu <- "test-weights-unbounded.R"
for (f in c("R/weights-unbounded.R", file.path("tests/testthat", wu),
            file.path("extensions/frmtmb.sample/tests/testthat", wu))) {
  wr(rd(f), f)
}
cat("applied dev/rel073-review.patch\n")
# finding 7: the relative clause needs its "that"
edit_file("R/families.R",
"#' its own interval. A row rejection has not filled after `inv_after`",
"#' its own interval. A row that rejection has not filled after `inv_after`")

# ---- found at the consolidation: the core site index --------------------
# lane loodraws added core's man/loo_predict.Rd (the ten generics and
# their refusals on a fit) and no _pkgdown.yml entry, so
# pkgdown::check_pkgdown() stops on a topic missing from the index
# (dev/rel073-pkgdown.R, dev/rel073-log/pkgdown.txt), and so would the
# site workflow. It goes beside loo() and bayes_R2()
restore("_pkgdown.yml")
edit_file("_pkgdown.yml", "  - loo\n  - log_lik\n  - bayes_R2\n",
          "  - loo\n  - log_lik\n  - bayes_R2\n  - loo_predict\n")
# core documents these four generics and their refusals on a fit on
# loo_predict.html now; the *.frmtmb_draws pages stay on frmtmb.sample's
# refusals (the release review, section 9, decision 3)
for (g in c("kfold", "reloo", "loo_moment_match", "loo_subsample")) {
  edit_file("_pkgdown.yml",
            paste0("  - [\"reference/", g, ".html\", ",
                   "\"frmtmb.sample/reference/frmtmb-loo-refusals.html\"]"),
            paste0("  - [\"reference/", g, ".html\", ",
                   "\"reference/loo_predict.html\"]"))
}

# ---- found at the consolidation: a link roxygen resolved two ways -------
# `[loo_model_weights()]` in ?posterior_average became
# \link[loo:loo_model_weights] when frmtmb.sample was already installed
# in the library roxygen ran against, and \link[=loo_model_weights] when
# it was not (dev/rel073-log/snap-c.txt against snap-d.txt), so the Rd
# depended on the library. The link names the topic it means: the draws
# method this package documents
edit_file("extensions/frmtmb.sample/R/postfit-draws.R",
"#' * `\"stacking\"` (the default) and `\"pseudobma\"`: the weights of
#'   [loo_model_weights()].",
"#' * `\"stacking\"` (the default) and `\"pseudobma\"`: the weights of
#'   [loo_model_weights()][sample-loo_model_weights].")

# ---- the gated one-session LEAK lines (brief item 6) --------------------
# core's generated test-brms-suite-emmeans.R attached emmeans for every
# later file of a session; the generator now attaches it for that file
# only. withr::local_package() at a file's top level detaches at the end
# of that file; with teardown_env() it would stay attached for every
# later file, as library() did (dev/rel073-local-package.R, testthat
# 3.3.2, withr 3.0.3). The generated file is written the same way here,
# so the tree is right before dev/brmsport-gen.R runs, and the
# generator's output is compared with it after (dev/round-20261021.md)
restore("dev/brmsport-gen.R")
edit_file("dev/brmsport-gen.R",
'            "library(emmeans)",',
'            "withr::local_package(\\"emmeans\\")",')
restore("tests/testthat/test-brms-suite-emmeans.R")
edit_file("tests/testthat/test-brms-suite-emmeans.R",
'library(emmeans)',
'withr::local_package("emmeans")')
# rstan's Rcpp module classes: measured to come from rstan itself
# (dev/rel073-rcpp-globals.R), so run-session.R tags them as a load does
rs <- "dev/release/run-session.R"
restore(rs)
edit_file(rs,
'    by <- if (k %in% c("globals", "search path")) character(0) else
      if (k %in% c("options", "envvars")) ns_new else
        grep("^frmtmb", ns_new, value = TRUE)
    cat("LEAK ", f, " ", k, ": ", paste(what, collapse = " "),
        if (length(by)) {
          paste0(" [with load of ", paste(sort(by), collapse = " "), "]")
        }, "\\n", sep = "")',
'    by <- if (k %in% c("globals", "search path")) character(0) else
      if (k %in% c("options", "envvars")) ns_new else
        grep("^frmtmb", ns_new, value = TRUE)
    # rstan defines a Stan model\'s Rcpp module classes, and the S4
    # class cache entries they bring, in the global environment when it
    # instantiates the model, whoever calls it: rstan::sampling() of a
    # cached model inside a function, with nothing else loaded, leaves
    # the same names (dev/rel073-rcpp-globals.R). They are rstan\'s, so
    # they are tagged as a load, on a line of their own; every other
    # global stays the test\'s
    if (identical(k, "globals")) {
      st <- grepl(paste0("^[+](\\\\.__C__Rcpp_stan_fit4model|",
                         "\\\\.__C__C[+][+]Object$|\\\\.__global__$|",
                         "\\\\.requireCachedGenerics$)"), what)
      if (any(st)) {
        cat("LEAK ", f, " ", k, ": ", paste(what[st], collapse = " "),
            " [with load of a Stan model (rstan)]\\n", sep = "")
        what <- what[!st]
        if (!length(what) || all(what == "-")) next
      }
    }
    cat("LEAK ", f, " ", k, ": ", paste(what, collapse = " "),
        if (length(by)) {
          paste0(" [with load of ", paste(sort(by), collapse = " "), "]")
        }, "\\n", sep = "")')

# ---- versions and floors (brief item 8) ---------------------------------
ver <- function(f, from, to) edit_file(f, paste0("Version: ", from),
                                       paste0("Version: ", to))
ver("DESCRIPTION", "0.72.0", "0.73.0")
restore("codemeta.json")
edit_file("codemeta.json", '"version": "0.72.0",', '"version": "0.73.0",')
# lane cdfs added two Suggests; codemeta lists every one
cm_pkg <- function(p) paste0(
'    {
      "@type": "SoftwareApplication",
      "identifier": "', p, '",
      "name": "', p, '",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=', p, '"
    },')
edit_file("codemeta.json",
          paste0(cm_pkg("fmesher")),
          paste0(cm_pkg("extraDistr"), "\n", cm_pkg("fmesher")))
edit_file("codemeta.json",
          paste0(cm_pkg("rstan")),
          paste0(cm_pkg("Rmpfr"), "\n", cm_pkg("rstan")))
ds <- "extensions/frmtmb.sample/DESCRIPTION"
restore(ds)
ver(ds, "0.20.0", "0.21.0")
# it imports the generics frmtmb 0.73.0 defines and calls
# check_weights_unbounded() (measured: dev/rel073-floor.sh)
edit_file(ds, "frmtmb (>= 0.72.0)", "frmtmb (>= 0.73.0)")

# ---- bin 2: the verdict moves of this consolidation ---------------------
# The stale overrides the lanes named go: negbinomial2 is a frmtmb
# family (log_lik:288, posterior_predict:169), and so is
# dirichlet_multinomial (log_lik:510, posterior_predict:355, which
# dev/rel072-bin2-moves.R set). The rows that now hold, and the
# re-judged ones, are dev/rel073-bin2-moves.R, sourced after the
# 0.72.0 moves in both verdict scripts
v2p <- "dev/brmsport-bin2-verdicts-2p.R"
restore(v2p)
edit_file(v2p,
'for (id in c("log_lik:288", "posterior_predict:169")) {
  add(id, "cannot transfer", "absent", absent_family("negbinomial2",
    paste("brms\'s sigma = 1 / shape spelling of the negative binomial;",
          "frmtmb\'s negbinomial() is the same density in shape, and the",
          "block\'s negbinomial row (log_lik:285) holds")))
}
', '')
edit_file(v2p,
'W <- mv072(W, "values")
V <- mv072(V, "2p")',
'W <- mv072(W, "values")
V <- mv072(V, "2p")
# the verdict moves of the 0.73.0 consolidation
mv073 <- source("dev/rel073-bin2-moves.R")$value
W <- mv073(W, "values")
V <- mv073(V, "2p")')
m072 <- "dev/rel072-bin2-moves.R"
restore(m072)
edit_file(m072,
'    # the reason the 0.71.0 table gave these two is stale: dirichlet and
    # logistic_normal are frmtmb families, and their rows hold
    for (id in c("log_lik:510", "posterior_predict:355")) {
      X <- set(X, id, "cannot transfer", "absent", paste(
        "dirichlet_multinomial is not a frmtmb family (\'not a supported",
        "family\'): a gap, not a decision"))
    }
', '    # (two rows set to "absent" here, log_lik:510 and
    # posterior_predict:355, hold since 0.73.0: dev/rel073-bin2-moves.R)
')
v2t <- "dev/brmsport-bin2-verdicts-2t.R"
restore(v2t)
edit_file(v2t,
'V <- source("dev/rel071-bin2-moves.R")$value(V, "2t")',
'V <- source("dev/rel071-bin2-moves.R")$value(V, "2t")
# the verdict moves of the 0.73.0 consolidation
V <- source("dev/rel073-bin2-moves.R")$value(V, "2t")')

# ---- the release harness reads rellib-r11 -------------------------------
# rellib-r11 is this release's library, the base of the next round. The
# verification here sets every library explicitly (FRMTMB_RELLIB,
# FRMTMB_PORT_LIB, ROUND_BASE), so no run reads the default before
# rellib-r11 exists
for (f in c("dev/brmsport-bin2-run.R", "dev/brmsport-guards.R",
            "dev/brmsport-run.R", "dev/ciharden-run1.R",
            "dev/release/run-check.ps1", "dev/release/run-docs.ps1",
            "dev/release/run-tests.R", "dev/warnleak-scan.R",
            "dev/release/run-session.R")) {
  if (!f %in% merged && f != rs) restore(f)
  s <- paste(rd(f), collapse = "\n")
  k <- lengths(regmatches(s, gregexpr("rellib-r10", s, fixed = TRUE)))
  stopifnot(k >= 1L)
  wr(strsplit(gsub("C:/Users/adf44/source/r/rellib-r10",
                   "C:/Users/adf44/source/r/rellib-r11", s, fixed = TRUE),
              "\n", fixed = TRUE)[[1]], f)
  cat("edited", f, "(", k, ")\n")
}
for (f in c("dev/brmsport-bin2-all.sh", "dev/brmsport-bin2-record.sh",
            "dev/ciharden-scan.sh")) {
  if (!f %in% merged) restore(f)
  s <- paste(rd(f), collapse = "\n")
  wr(strsplit(gsub("rellib-r10", "rellib-r11", s, fixed = TRUE), "\n",
              fixed = TRUE)[[1]], f)
  cat("edited", f, "\n")
}

# ---- dev/test-backlog.md: closures in place, then this release's section
bl <- "dev/test-backlog.md"
restore(bl)
close_item <- function(old, note) edit_file(bl, old, paste0(old, note))
# each old text and note is given line by line, joined with newlines,
# so that no line of this script passes 80 columns
L <- function(...) paste(c(...), collapse = "\n")
# the release review, finding 6 (dev/rel073-rev-adhess19.R)
close_item(
  L(paste("- `y | trunc(lb = 1) ~ x` with poisson: the gradient at the",
          "optimum is"),
    "  not finite and `obj$fn()` at `last.par.best` is `Inf` (seed 19)."),
  L(" CLOSED",
    "  at 0.73.0 by lane cdfs's log-survivor normalizer: relative",
    "  convergence, no warning, gradient 7.6e-07 and 3.5e-04, and",
    "  `obj$fn(last.par.best)` 333.752789196, an independent",
    "  zero-truncated fit on `stats::dpois()` and `stats::ppois()` to 12",
    "  digits; SEs 0.35776564 and 0.00108982 against 0.35776574 and",
    "  0.00108982 (the release review, `dev/rel073-rev-adhess19.R`)."))
close_item(
  L("  a poisson `cens()` fit refuses with \"cens() is not supported for",
    "  discrete families yet (truncation is)\"."),
  L("", "  CLOSED: poisson censors since 0.54.0, and at 0.73.0 every",
    "  discrete family brms 2.23.0 censors does (lane cdfs)."))
close_item(
  L(paste("  holding 1e-40 of the mass before it, against Rmpfr, through",
          "the core"),
    "  path."),
  L("",
    "  CLOSED in core at 0.73.0 (lane cdfs): `frmtmb_family(cdf_log = )`",
    "  and `log_window()` form every censored row and truncation window on",
    "  the log scale. frmtmb.eam's `wiener()` declares no `cdf_log` and",
    "  keeps its workaround and its refusal."))
close_item(
  L(paste("- `test-brms-suite-emmeans.R` (generated by",
          "`dev/brmsport-gen.R`) calls"),
    "  `library(emmeans)` at file level, so a gated one-session run has",
    "  emmeans attached for every later file."),
  L("", "  CLOSED at 0.73.0: the generator attaches it with",
    "  `withr::local_package()` for that file only."))
close_item(
  L(paste("- The one-session runs are `NOT_CRAN=true` only; the gated tiers",
          "were"),
    "  not run in one session."),
  L(" CLOSED at 0.73.0: the release ran both",
    "  (`dev/round-20261021.md`)."))
close_item(
  L(paste("- The bin-2 rows `log_lik:510` and `posterior_predict:355` stay",
          "gaps:"),
    "  frmtmb has no `dirichlet_multinomial()`."),
  L(" CLOSED at 0.73.0 (lane nb2dm):", "  both hold."))
edit_file(bl, "\n## Reference\n",
          paste0("\n", paste(rd("dev/rel073-backlog-section.md"),
                             collapse = "\n"), "\n\n## Reference\n"))

# ---- dev/round-handoff.md: the 0.73.0 section ---------------------------
hf <- "dev/round-handoff.md"
restore(hf)
edit_file(hf,
"0.70.0 release, 2026-10-14 at the 0.71.0 release and 2026-10-20 at
the 0.72.0 release. Read this, then",
"0.70.0 release, 2026-10-14 at the 0.71.0 release, 2026-10-20 at the
0.72.0 release and 2026-10-21 at the 0.73.0 release. Read this, then")
edit_file(hf, "\n## The 0.72.0 round, 2026-10-14 to 2026-10-20\n",
          paste0("\n", paste(rd("dev/rel073-handoff-section.md"),
                             collapse = "\n"),
                 "\n\n## The 0.72.0 round, 2026-10-14 to 2026-10-20\n"))
cat("RESOLVE DONE\n")
