# The NEWS sections of frmtmb 0.80.0 and its extensions, from each lane's
# own development section (read from the lane trees), grouped by subject,
# with every lane's limits in one Limits section:
#
#   Rscript dev/rel080-news.R
#
# NEWS carries no lane names: the lane of each group is named only here.
R <- "C:/Users/adf44/source/r"
REL <- file.path(R, "frmtmb-wt-rel080")
rl <- function(p) readLines(p, warn = FALSE)
wr <- function(x, p) {
  con <- file(p, open = "wb")
  on.exit(close(con))
  writeLines(x, con, sep = "\r\n")
}
# a lane's development section, split into its sub-sections: a named list
# of character vectors, "main" for the bullets before any sub-heading
lane_section <- function(lane, file = "NEWS.md") {
  x <- rl(file.path(R, paste0("frmtmb-wt-", lane), file))
  stopifnot(grepl("(development version)", x[1L], fixed = TRUE))
  end <- grep("^# ", x)[2L] - 1L
  x <- x[2:end]
  hd <- grepl("^## ", x)
  h <- cumsum(hd)
  nm <- c("main", sub("^## ", "", x[hd]))
  out <- split(x[!hd], factor(nm[h[!hd] + 1L], levels = unique(nm)))
  out <- lapply(out, function(v) {
    # trim blank lines at both ends
    k <- which(nzchar(trimws(v)))
    if (length(k)) v[min(k):max(k)] else character()
  })
  out[lengths(out) > 0L]
}
fix <- function(v, old, new) {
  s <- paste(v, collapse = "\n")
  n <- lengths(regmatches(s, gregexpr(old, s, perl = TRUE)))
  if (n != 1L) stop("NEWS fix matched ", n, " times: ", old)
  strsplit(gsub(old, new, s, perl = TRUE), "\n", fixed = TRUE)[[1L]]
}
sec <- lapply(stats::setNames(nm = c("ci", "latperf", "mfxse", "loglik", "osa",
                                     "nestquad", "arstep", "fitedge", "bigx",
                                     "nlstart", "readers", "superassign")),
              lane_section)
# arstep review, nit on bullet 3 (the brief's release fix 15)
sec$arstep[["main"]] <- fix(sec$arstep[["main"]],
  "and 2\\.2 to 2\\.9 times on 600 rows of 50 positions, under 1\\n  second in each\\.",
  "and 2 to 5 times on 600 rows of 50 positions, under about 2\n  seconds in each.")
# the release review of 0.80.0 (dev/reviews/2026-11-07-release.md):
# item 3, numbers the final build reads (nlstart's from
# dev/rel080-log/rr-fix3/nlstart-freshsum.txt, after item 1)
nl_short <- 11L; nl_silent <- 5L; nl_wins <- 87L
sec$nlstart[["main"]] <- fix(sec$nlstart[["main"]],
  "now 11 are, 4 with no\\n    warning\\.",
  sprintf("now %d are, %d with no\n    warning.", nl_short, nl_silent))
sec$nlstart[["main"]] <- fix(sec$nlstart[["main"]],
  "on 1 of the 89\\n", sprintf("on 1 of the %d\n", nl_wins))
# the short check of the release fixes: the fifth silent row is named
sec$nlstart[["main"]] <- fix(sec$nlstart[["main"]],
  "warning\\. They are Michaelis-Menten bodies under `bernoulli\\(\\)` and\\n    `binomial\\(\\)`, a power law under `cumulative\\(\\)`, and fits that warn\\.",
  paste0("warning. They are Michaelis-Menten bodies under `bernoulli()` and\n",
         "    `binomial()`, a power law under `cumulative()`, a four-parameter\n",
         "    logistic under `bernoulli()` 2.8e-4 short, and fits that warn."))
sec$fitedge[["Bug fixes"]] <- fix(sec$fitedge[["Bug fixes"]],
  "Now 9 stop short, by up to 36\\.5 nats, also with no\\n  warning, and none reads below its old value\\. 5 of the 9 are low-count",
  "Now 8 stop short, by up to 36.5 nats, also with no\n  warning, and none reads below its old value. 5 of the 8 are low-count")
# item 4, BREAKING marks
sec$arstep[["main"]] <- fix(sec$arstep[["main"]],
  "^\\* An exact `gp\\(\\)` beside `ar\\(\\)` on a gaussian response now starts at the\\n  maximum of the model's own marginal\\. That marginal includes the `ar`\\n  coefficients\\. Before, the start read a marginal of the response alone\\.",
  "* BREAKING: an exact `gp()` beside `ar()` on a gaussian response now\n  starts at the maximum of the model's own marginal. That marginal\n  includes the `ar` coefficients. Before, the start read a marginal of the\n  response alone.")
sec$arstep[["main"]] <- fix(sec$arstep[["main"]],
  "\\n\\* Off the gaussian, where the start's marginal of the response has a\\n  second maximum",
  "\n* BREAKING: off the gaussian, where the start's marginal of the response\n  has a second maximum")
sec$readers[["Fixes"]] <- fix(sec$readers[["Fixes"]],
  "\\* `conditional_effects\\(\\)`'s band for a category of a row with a `me\\(\\)`\\n  value is `fitted\\(\\)`'s,",
  "* BREAKING: `conditional_effects()`'s band for a category of a row with a\n  `me()` value is `fitted()`'s,")
conv_classes <- c(
  "* The convergence warnings of a fit have their own classes, before",
  "  `frmtmb_warning`: `frmtmb_conv_code` (\"Optimizer did not report",
  "  convergence\"), `frmtmb_conv_gradient` (\"Large maximum absolute",
  "  gradient\") and `frmtmb_conv_gradient_nonfinite` (\"The gradient at the",
  "  reported optimum is not finite\"). The messages do not change. Code that",
  "  reads `class(w)[1]` sees the new class; `inherits(w, \"frmtmb_warning\")`",
  "  is still true.")
# item 2, the ordinal Limits bullet, after the bullet on a column far from
# zero against its spread
sec$bigx[["Limits"]] <- fix(sec$bigx[["Limits"]],
  "`\\(1 \\+ x \\| g\\)`, stops short of the maximum\\. Center the column\\.\\n",
  paste0("`(1 + x | g)`, stops short of the maximum. Center the column.\n",
         "* Under an ordinal family, whose thresholds take the intercept's place,\n",
         "  products and powers of columns shifted by 1e6 (`w1 * w2`,\n",
         "  `w1 + I(w1^2)`) stop 8.6 to 28.6 log-likelihood units short of the\n",
         "  model on the centered columns, each fit with a convergence warning.\n",
         "  Center the columns.\n"))
stopifnot(identical(names(sec$arstep), c("main", "Limits")),
          identical(names(sec$fitedge), c("Bug fixes", "Performance", "Limits")),
          identical(names(sec$bigx), c("main", "Limits")),
          identical(names(sec$readers), c("Breaking changes", "Fixes")),
          all(vapply(sec[c("ci", "latperf", "mfxse", "loglik", "osa", "nestquad",
                           "nlstart", "superassign")],
                     function(s) identical(names(s), "main"), NA)))
# measured at the consolidation (dev/rel080-nestshift.R): the nested
# rounds step in natural coordinates, and only the final tape takes the
# design's basis
nest_limit <- c(
  "* The rounds of nested quadrature also run in natural coordinates. On",
  "  a bernoulli `(1 | school/class)` fit with `x + 1e4` or `x + 1e6`, the",
  "  fit was 0.0007 to 0.039 below the centered fit and its slope up to",
  "  0.06 standard errors off, with a warning that the rule did not",
  "  settle and one of false convergence. Center the column.")
group <- function(title, parts) {
  c(paste("##", title), "", unlist(lapply(parts, function(p) c(p, ""))))
}
core <- c(
  "# frmtmb 0.80.0", "",
  group("Fitting", list(sec$nestquad[["main"]], sec$bigx[["main"]], sec$nlstart[["main"]],
                        sec$fitedge[["Bug fixes"]], sec$arstep[["main"]],
                        conv_classes)),
  group("Log-likelihood and priors", list(sec$loglik[["main"]])),
  group("Predictions and their standard errors",
        list(sec$readers[["Breaking changes"]], sec$readers[["Fixes"]],
             sec$latperf[["main"]], sec$mfxse[["main"]])),
  group("Residuals", list(sec$osa[["main"]])),
  group("Performance", list(sec$fitedge[["Performance"]])),
  group("Package code and continuous integration",
        list(sec$superassign[["main"]], sec$ci[["main"]])),
  group("Limits", list(sec$arstep[["Limits"]], sec$fitedge[["Limits"]],
                       sec$bigx[["Limits"]], nest_limit)))
x <- rl(file.path(REL, "NEWS.md"))
old <- grep("^# frmtmb 0\\.79\\.1$", x)
# a rerun replaces the section it wrote
stopifnot(length(old) == 1L,
          grepl("(development version)", x[1L], fixed = TRUE) ||
            identical(x[1L], "# frmtmb 0.80.0"))
wr(c(core, x[old:length(x)]), file.path(REL, "NEWS.md"))
cat("NEWS core:", length(core), "lines\n")

# extensions: one section each, the lanes' bullets, the floor named
ext <- list(
  frmtmb.sample = list(ver = "0.27.0", lanes = c("readers", "superassign")),
  frmtmb.latent = list(ver = "0.7.3", lanes = "loglik"),
  frmtmb.ode = list(ver = "0.7.3", lanes = "superassign"),
  frmtmb.survival = list(ver = "0.2.2", lanes = "superassign"))
for (pk in names(ext)) {
  f <- file.path("extensions", pk, "NEWS.md")
  b <- unlist(lapply(ext[[pk]]$lanes, function(l) c(lane_section(l, f)[["main"]], "")))
  s <- paste(b, collapse = "\n")
  s <- gsub("Needs the frmtmb release\\s+that carries this change\\.",
            "Needs frmtmb 0.80.0.", s)
  s <- gsub("Needs the frmtmb release that adds\\s+`logLik_penalized\\(\\)`\\.",
            "Needs frmtmb 0.80.0, which adds\n  `logLik_penalized()`.", s)
  b <- c(strsplit(s, "\n", fixed = TRUE)[[1L]], "")
  x <- rl(file.path(REL, f))
  old <- grep(paste0("^# ", gsub(".", "\\.", pk, fixed = TRUE), " [0-9]"), x)
  # a rerun replaces the section it wrote
  old <- setdiff(old, grep(paste0(" ", ext[[pk]]$ver, "$"), x))[1L]
  wr(c(paste0("# ", pk, " ", ext[[pk]]$ver), "", b, x[old:length(x)]),
     file.path(REL, f))
  cat("NEWS", pk, ext[[pk]]$ver, "\n")
}
