# The NEWS sections of frmtmb 0.81.0 and its changed extensions, from each
# lane's own development section (read from the lane trees), grouped by
# subject, with one Limits section:
#
#   Rscript dev/rel081-news.R
#
# NEWS carries no lane names: the lane of each group is named only here
# (compat: Compatibility table and the stage2() bug fix; code8:
# Convergence and standard errors; hollow: Hypothesis tests and the
# make_conditions() fix, frmtmb.sample's newdata bullets; support:
# Families, frmtmb.sample's exact log_lik(), frmtmb.survival).
R <- "C:/Users/adf44/source/r"
REL <- file.path(R, "frmtmb-wt-rel081")
rl <- function(p) sub("\r$", "", readLines(p, warn = FALSE))
wr <- function(x, p) {
  # an element written as one string over several lines is split, so
  # every line ends CRLF
  # (strsplit() of "" is empty, so a blank line is kept as it is)
  x <- unlist(lapply(x, function(s) {
    if (grepl("\n", s, fixed = TRUE)) strsplit(s, "\n", fixed = TRUE)[[1L]]
    else s
  }), use.names = FALSE)
  con <- file(p, open = "wb")
  on.exit(close(con))
  writeLines(x, con, sep = "\r\n")
}
# the released sections are read from the base commit, so a rerun never
# reads back its own output
base_news <- function(f) {
  sub("\r$", "", system2("git", c("-C", REL, "show", paste0("e0cce906:", f)),
                         stdout = TRUE))
}
# the first line of a NEWS file before this script ran, or after it ran
# once (the script then replaces its own section)
is_head <- function(l, ver) {
  grepl("(development version)", l, fixed = TRUE) || identical(l, ver)
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
    k <- which(nzchar(trimws(v)))
    if (length(k)) v[min(k):max(k)] else character()
  })
  out[lengths(out) > 0L]
}
# split a block of lines into its bullets (each starts with "* ")
bullets <- function(v) {
  v <- v[nzchar(trimws(v))]
  split(v, cumsum(grepl("^\\* ", v)))
}
fix <- function(v, old, new) {
  s <- paste(v, collapse = "\n")
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop("NEWS fix matched ", n, " times: ", old)
  strsplit(sub(old, new, s, fixed = TRUE), "\n", fixed = TRUE)[[1L]]
}
join <- function(bl) unlist(lapply(bl, function(b) c(b, "")), use.names = FALSE)

sec <- lapply(stats::setNames(nm = c("compat", "code8", "hollow", "support")),
              lane_section)

## Families (support): bullet 3 is a limit; bullet 4 names the von_mises()
## bound past pi (the support review's short check, note 1) and is marked
## BREAKING as the lane's record classes it (estimates change)
fam <- bullets(sec$support[["Families"]])
stopifnot(length(fam) == 4L)
fam_limit <- fam[[3L]]
fam4 <- fix(fam[[4L]],
"* `trunc(lb = )` at or below the lower end of the support is no
  truncation for `exponential()`, `weibull()`, `frechet()`,
  `shifted_lognormal()` and `von_mises()`, as for the other families.
  Before, `exponential()` added `lb / mu` to every row and gave wrong
  estimates without a warning (an intercept of 1.078 where the fit
  without `trunc()` gives 0.661), `von_mises()` with `lb = -4` read an
  arc of negative length and moved `logLik()` by 6.4, and `weibull()`,
  `frechet()` and `shifted_lognormal()` stopped at a `NaN` gradient.
  The shifted lognormal's floor at its shift had rounded to 0.",
"* BREAKING: A `trunc()` bound outside the support is no truncation,
  as for the other families: `lb` at or below the lower end for
  `exponential()`, `weibull()`, `frechet()` and `shifted_lognormal()`,
  and `lb` below -pi or `ub` above pi for `von_mises()`. Before,
  `exponential()` added `lb / mu` to every row and gave wrong estimates
  without a warning (an intercept of 1.078 where the fit without
  `trunc()` gives 0.661). `von_mises()` read an arc of negative length
  and moved `logLik()` by 6.4 with `lb = -4`, by 38.6 with `ub = 4`,
  and by 188 with `cens()` beside `lb = -4, ub = 4`, there with a
  parameter 69.1 away. `weibull()`, `frechet()` and
  `shifted_lognormal()` stopped at a `NaN` gradient. The shifted
  lognormal's floor at its shift had rounded to 0.")
families <- join(c(fam[1:2], list(fam4)))

## Bug fixes: compat's stage2() refusal, and make_conditions() of draws
## (hollow's core change, which its NEWS did not name)
bug <- c(sec$compat[["Bug fixes"]], "",
"* `make_conditions()` of draws from `frmtmb.sample::frm_sample()` reads
  the data of the fit they were drawn from, as it does for the fit. The
  draws now carry the fit's model frame as `$data`, which holds a
  transform's basis, such as `poly(x, 2)`, and not its variable, so
  `make_conditions(draws, \"x\")` stopped with \"`x` is not a column of
  the data\".", "")

limits <- join(list(
  fam_limit,
c("* A correlation row held at its edge in an unstructured block of three",
  "  or more loses the standard errors of all of its entries, although",
  "  the row's direction is determined: on `(1 + poly(x, 2) | g)` with a",
  "  tight prior on one standard deviation, `cor(Intercept, polyx22)` is",
  "  0.41 and has no standard error. Only the row's length is at the",
  "  boundary."),
c("* On a ridge that the standard-error check measures flat, the",
  "  convergence check asks the objective only along the flat directions.",
  "  Under a weak prior across the ridge (`normal(0.5, 1000)` on `k` of",
  "  `c0 + exp(a)^k`), two stops 1.9e-4 and 1.2e-3 short along it lost",
  "  the gradient warning that 0.80.0 gave. Stops short by 0.003 or more",
  "  warn, and none gets the message `frmtmb_converged_at_max`."),
c("* Some densities outside a support still read a finite floor, also for",
  "  a reader: an empty or unresolvable censoring or truncation window",
  "  (about -690.8), tied thresholds of an ordinal mixture, the tails of",
  "  OSA residuals (-1e4), and the floors of frmtmb.eam's and",
  "  frmtmb.survival's composed densities. A fit or a draw does not reach",
  "  the last of these: their links keep each fitted or sampled row away",
  "  from the floor."),
c("* `weights()`, `subset()` and `index()` with `disc`, `equidistant` and",
  "  `sum_to_zero` stay \"untested\" in `frm_compat()`: the ordinal",
  "  families take these terms, and a fit that runs is not a measured",
  "  agreement with brms."),
c("* `frmtmb_register_aterm()` does not refuse a term whose argument",
  "  keys collide with those of another term, such as `x` of three",
  "  arguments beside `x2`: there `y | x2(z) + x(u, v, w)` lets the",
  "  second argument of `x()` replace `x2(z)` without a word. No shipped",
  "  extension has such a pair.")))

## BREAKING marks the release review asked for (finding 3): code8's
## bullets 2 to 4 change standard errors and the gradient warning. The
## mark lengthens the first line, so the bullet is wrapped again at 72
## columns. compat's bullet has none: frm() already refused each of
## those models, and only rows of frm_compat() change (the review's
## short check).
breaking <- function(b) {
  stopifnot(grepl("^\\* ", b[1L]), !grepl("BREAKING", b[1L], fixed = TRUE))
  txt <- paste(trimws(b), collapse = " ")
  txt <- sub("^\\* ", "", txt)
  strwrap(paste("BREAKING:", txt), width = 73, initial = "* ",
          prefix = "  ")
}
conv <- bullets(sec$code8[["Convergence and standard errors"]])
stopifnot(length(conv) == 5L)
conv[2:4] <- lapply(conv[2:4], breaking)
comp <- bullets(sec$compat[["Compatibility table"]])
stopifnot(length(comp) == 1L)

new_core <- c("# frmtmb 0.81.0", "",
  "## Convergence and standard errors", "",
  join(conv),
  "## Families", "", families,
  "## Hypothesis tests", "",
  join(bullets(sec$hollow[["Hypothesis tests"]])),
  "## Compatibility table", "",
  join(comp),
  "## Bug fixes", "", join(bullets(bug)),
  "## Limits", "", limits)
new_core <- new_core[seq_len(max(which(nzchar(new_core))))]

p <- file.path(REL, "NEWS.md")
x <- c(rl(p)[1L], base_news("NEWS.md"))
stopifnot(is_head(x[1L], "# frmtmb 0.81.0"))
old_start <- grep("^# frmtmb 0\\.80\\.0$", x)
stopifnot(length(old_start) == 1L)
wr(c(new_core, "", x[old_start:length(x)]), p)
cat("NEWS.md: frmtmb 0.81.0,", sum(grepl("^\\* ", new_core)), "bullets\n")

## frmtmb.sample: hollow's three bullets, then support's; the floor named
s_h <- lane_section("hollow", "extensions/frmtmb.sample/NEWS.md")[["main"]]
s_s <- lane_section("support", "extensions/frmtmb.sample/NEWS.md")[["main"]]
# the section's first line names the floor; the bullet's own copy of it
# goes (release review, finding 4)
s_s <- fix(s_s, "\n  Needs the frmtmb that adds `row_lpdf(exact = )`.", "")
# the release's decision on a log_lik() cell that is not finite
# (dev/rel081-loo-inf.R; release review, finding 2, measured 0 of 44,000
# sampled draws outside in 22 fits that read the barrier)
s_add <- c(
  "* BREAKING: `loo()`, `waic()` and `psis()` of draws stop with an error",
  "  of class `frmtmb_loglik_not_finite` when a cell of `log_lik()` is not",
  "  finite, and give the number of cells, observations and draws. The",
  "  `loo_*()` functions and the `loo_*` types of `pp_check()` stop in the",
  "  same way. A `-Inf` cell is an observation that a draw puts outside",
  "  the support of the response's distribution. Before, `loo()` and",
  "  `psis()` stopped in the loo package with \"All input values must be",
  "  finite or -Inf\", and `waic()` gave an elpd of `NaN` with no error.",
  "  brms stops in all three, because its density there is `NaN`. Draws",
  "  from `frm_sample()` were not seen to make such a cell: in 22",
  "  `gen_extreme_value()` fits that sample the fitted objective (a",
  "  predicted `xi`, or `cens()`), 0 of 44,000 draws put a row outside",
  "  the support. The objective is continuous across the end of the",
  "  support and falls by 1e8 per unit of `t` outside it, so a draw",
  "  outside must lie in a band about 1e-8 wide. A draws object built or",
  "  edited by hand can have such a cell.")
new_s <- c("# frmtmb.sample 0.28.0", "",
           "Needs frmtmb 0.81.0, which adds `row_lpdf(exact = )`.", "",
           join(c(bullets(s_h), bullets(s_s), list(s_add))))
new_s <- new_s[seq_len(max(which(nzchar(new_s))))]
p <- file.path(REL, "extensions/frmtmb.sample/NEWS.md")
x <- c(rl(p)[1L], base_news("extensions/frmtmb.sample/NEWS.md"))
stopifnot(is_head(x[1L], "# frmtmb.sample 0.28.0"))
old_start <- grep("^# frmtmb\\.sample 0\\.27\\.0$", x)
stopifnot(length(old_start) == 1L)
wr(c(new_s, "", x[old_start:length(x)]), p)
cat("frmtmb.sample NEWS: 0.28.0,", sum(grepl("^\\* ", new_s)), "bullets\n")

## frmtmb.survival: support's bullet
s_v <- lane_section("support", "extensions/frmtmb.survival/NEWS.md")[["main"]]
new_v <- c("# frmtmb.survival 0.2.3", "", join(bullets(s_v)))
new_v <- new_v[seq_len(max(which(nzchar(new_v))))]
p <- file.path(REL, "extensions/frmtmb.survival/NEWS.md")
x <- c(rl(p)[1L], base_news("extensions/frmtmb.survival/NEWS.md"))
stopifnot(is_head(x[1L], "# frmtmb.survival 0.2.3"))
old_start <- grep("^# frmtmb\\.survival 0\\.2\\.2$", x)
stopifnot(length(old_start) == 1L)
wr(c(new_v, "", x[old_start:length(x)]), p)
cat("frmtmb.survival NEWS: 0.2.3,", sum(grepl("^\\* ", new_v)), "bullets\n")
cat("NEWS DONE\n")
