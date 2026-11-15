# The NEWS sections of frmtmb 0.83.0 and its changed extensions, built from
# each lane's own development section (read from the lane trees), grouped
# by subject. Rerunnable: the released sections are read from 8159568e,
# never from this script's output. Adapted from dev/rel082-news.R.
#
#   Rscript dev/rel083-news.R
#
# NEWS carries no lane names, so the lane of each group is named only here:
# selimit2, slopescale and quadseam2 ("Convergence"): Convergence and
# standard errors; quadseam2: Quadrature; nlstop2 and the consolidation's
# nl_starts bullet: Fitting; gpfield: Gaussian processes; priors2 and
# famend: Priors; priors2, silent2 and famend: Ordinal mixtures; famend:
# Families, Residuals; wcheck2: Weighted models with an observation-level
# effect, asym_laplace() with random effects; arspat2 and quadseam2:
# Residual correlation terms; refits: Response coding and refits; parity2:
# brms surface; predraw2: Predictions; silent2: Silent wrong answers fixed;
# hygiene: Compatibility table; hygiene and parity2: For extension
# packages; nlstop2 and hygiene: Documentation. Extensions: frmtmb.sample
# (priors2, silent2, parity2), frmtmb.spatial (arspat2), frmtmb.latent
# (quadseam2's and refits' tests).
R <- "C:/Users/adf44/source/r"
REL <- file.path(R, "frmtmb-wt-rel083")
BASE <- "ec8178d3"      # the lanes' base
RELBASE <- "8159568e"   # the release tree's base, 0.82.1
rl <- function(p) sub("\r$", "", readLines(p, warn = FALSE))
wr <- function(x, p) {
  con <- file(p, open = "wb")
  on.exit(close(con))
  writeLines(x, con, sep = "\r\n")
}
git_show <- function(rev, f) {
  sub("\r$", "", system2("git", c("-C", REL, "show", paste0(rev, ":", f)),
                         stdout = TRUE))
}
# a lane's development section: the lines above its first released heading
dev_part <- function(lane, f) {
  x <- rl(file.path(R, paste0("frmtmb-wt-", lane), f))
  rel <- grep("^# \\S+ [0-9]", x)[1L]
  b <- git_show(BASE, f)
  stopifnot(identical(x[rel:length(x)], b))
  x <- x[seq_len(rel - 1L)]
  x[-1L]  # the "# <pkg> (development version)" heading
}
blocks <- function(x) {
  out <- list(); cur <- character()
  for (l in x) {
    if (!nzchar(trimws(l)) || grepl("^\\* ", l)) {
      if (length(cur)) out[[length(out) + 1L]] <- cur
      cur <- if (nzchar(trimws(l))) l else character()
    } else cur <- c(cur, l)
  }
  if (length(cur)) out[[length(out) + 1L]] <- cur
  out
}
sections <- function(x) {
  h <- grep("^## ", x)
  starts <- c(1L, h); ends <- c(h - 1L, length(x))
  out <- list()
  for (i in seq_along(starts)) {
    if (starts[i] > ends[i]) next
    seg <- x[starts[i]:ends[i]]
    nm <- if (i == 1L) "" else sub("^## ", "", seg[1L])
    body <- if (i == 1L) seg else seg[-1L]
    bl <- blocks(body)
    if (length(bl)) out[[length(out) + 1L]] <- list(name = nm, blocks = bl)
  }
  out
}
emit <- function(heading, secs) {
  y <- c(heading, "")
  for (s in secs) {
    if (nzchar(s$name)) y <- c(y, paste("##", s$name), "")
    for (b in s$blocks) y <- c(y, b, "")
  }
  y
}
taken <- list()
take <- function(lane, f, name, pick = NULL, drop = NULL) {
  s <- sections(dev_part(lane, f))
  hit <- Filter(function(z) identical(z$name, name), s)
  if (length(hit) != 1L) stop(lane, ": section '", name, "' found ", length(hit))
  b <- hit[[1L]]$blocks
  first <- vapply(b, `[`, "", 1L)
  keep <- rep(TRUE, length(b))
  if (!is.null(pick)) keep <- grepl(pick, first, fixed = TRUE)
  if (!is.null(drop)) keep <- !grepl(drop, first, fixed = TRUE)
  if (!is.null(pick) || !is.null(drop)) {
    stopifnot(sum(keep) >= 1L, sum(!keep) >= 1L || is.null(drop))
  }
  key <- paste(lane, f)
  taken[[key]] <<- c(taken[[key]], first[keep])
  b[keep]
}
rewrap <- function(b) {
  txt <- paste(trimws(b), collapse = " ")
  strwrap(sub("^\\* ", "", txt), width = 72, initial = "* ", prefix = "  ")
}
# merge step 5, the one bullet the consolidation adds
nl_starts <- list(rewrap(c(
  "* `frmtmb_control(nl_starts = FALSE)` turns off the sign probe and the",
  "pole restart of a nonlinear fit: the fits from each sign pattern of",
  "the body's coefficients and the fits from beyond a pole of the body.",
  "The default, `TRUE`, keeps them. A fit can still restart from a saddle",
  "point (message class `frmtmb_nl_restart`); this switch does not change",
  "that."
)))
# priors2's bullet names a version as provenance (NEWS carries none)
no_version <- function(b) {
  i <- grep("0.82.0 ended at the bound with the objective", b, fixed = TRUE)
  if (!length(i)) return(b)
  b[i] <- sub("0.82.0 ended at the bound with the objective",
              "the fit ended at the bound before, with the objective", b[i],
              fixed = TRUE)
  rewrap(b)
}

## core
N <- "NEWS.md"
core <- list(
  list(name = "Convergence and standard errors",
       blocks = c(take("selimit2", N, "Convergence and standard errors"),
                  take("slopescale", N, "Convergence and standard errors"),
                  take("quadseam2", N, "Convergence"))),
  list(name = "Quadrature", blocks = take("quadseam2", N, "Quadrature")),
  list(name = "Fitting", blocks = c(take("nlstop2", N, "Fitting"), nl_starts)),
  list(name = "Gaussian processes",
       blocks = take("gpfield", N, "Gaussian processes")),
  list(name = "Priors",
       blocks = c(lapply(take("priors2", N, "Priors"), no_version),
                  take("famend", N, "Priors"))),
  list(name = "Ordinal mixtures",
       blocks = c(take("priors2", N, "Mixtures"),
                  take("silent2", N, "", pick = "ordinal-mixture"),
                  take("silent2", N, "", pick = "an ordinal mixture finds"),
                  take("famend", N, "Families", pick = "An ordinal mixture"))),
  list(name = "Families",
       blocks = take("famend", N, "Families", drop = "An ordinal mixture")),
  list(name = "Residuals", blocks = take("famend", N, "Residuals")),
  list(name = "Weighted models with an observation-level effect",
       blocks = take("wcheck2", N, "Weighted models with an observation-level effect")),
  list(name = "asym_laplace() with random effects",
       blocks = take("wcheck2", N, "asym_laplace() with random effects")),
  list(name = "Residual correlation terms",
       blocks = c(take("arspat2", N, "Residual correlation terms"),
                  take("quadseam2", N, "Residual correlation terms"))),
  list(name = "Response coding and refits",
       blocks = take("refits", N, "Response coding and refits")),
  list(name = "brms surface",
       blocks = take("parity2", N, "brms surface", drop = "The sampling API")),
  list(name = "Predictions", blocks = take("predraw2", N, "Predictions")),
  list(name = "Silent wrong answers fixed",
       blocks = take("silent2", N, "", pick = "BREAKING: ")),
  list(name = "Compatibility table",
       blocks = take("hygiene", N, "Compatibility table")),
  list(name = "For extension packages",
       blocks = c(take("hygiene", N, "For extension packages"),
                  take("parity2", N, "brms surface", pick = "The sampling API"))),
  list(name = "Documentation",
       blocks = c(take("nlstop2", N, "Documentation"),
                  take("hygiene", N, "Documentation")))
)
# silent2's "BREAKING: " pick also takes its two ordinal-mixture bullets;
# they belong under Ordinal mixtures only
om <- unlist(lapply(core[[6L]]$blocks, `[`, 1L))
k <- which(vapply(core, function(s) s$name, "") == "Silent wrong answers fixed")
core[[k]]$blocks <- Filter(function(b) !b[1L] %in% om, core[[k]]$blocks)
# every lane's development section is taken once and nothing is left out
lanes <- c("priors2", "slopescale", "quadseam2", "selimit2", "nlstop2",
           "gpfield", "silent2", "wcheck2", "arspat2", "parity2", "refits",
           "famend", "predraw2", "hygiene")
lane_firsts <- unlist(lapply(lanes, function(l) {
  unlist(lapply(sections(dev_part(l, N)), function(s) {
    vapply(s$blocks, `[`, "", 1L)
  }))
}))
out_firsts <- unlist(lapply(core, function(s) vapply(s$blocks, `[`, "", 1L)))
out_firsts <- setdiff(out_firsts, nl_starts[[1L]][1L])
stopifnot(length(out_firsts) == length(lane_firsts),
          setequal(out_firsts, lane_firsts), !anyDuplicated(out_firsts))
cat("core: lane blocks", length(lane_firsts), "| sections", length(core),
    "| BREAKING", sum(grepl("^\\* BREAKING", out_firsts)), "\n")
txt <- unlist(lapply(core, function(s) unlist(s$blocks)))
stopifnot(!any(grepl("0\\.8[0-9]\\.[0-9]", txt)), !any(grepl("dev/", txt)),
          !any(grepl("\u2014", txt)))
x <- rl(file.path(REL, N))
b <- git_show(RELBASE, N)
stopifnot(grepl("^# frmtmb 0.82.1", b[1L]))
wr(c(emit("# frmtmb 0.83.0", core), b), file.path(REL, N))

## frmtmb.sample
SN <- "extensions/frmtmb.sample/NEWS.md"
sb <- c(
  list("Needs frmtmb 0.83.0: it reads class `\"theta\"` as a mixture's weights and calls `ce_snl()` and `latres_err()` of the sampling API."),
  take("priors2", SN, "", pick = "* "),
  take("silent2", SN, ""),
  take("parity2", SN, "", pick = "* ")
)
sb[[1L]] <- strwrap(sb[[1L]], width = 72)
sb_n <- sum(vapply(c("priors2", "silent2", "parity2"), function(l) {
  sum(grepl("^\\* ", dev_part(l, SN)))
}, 1L))
stopifnot(sum(vapply(sb, function(b) grepl("^\\* ", b[1L]), NA)) == sb_n)
wr(c(emit("# frmtmb.sample 0.30.0", list(list(name = "", blocks = sb))),
     git_show(RELBASE, SN)), file.path(REL, SN))
cat("frmtmb.sample:", sb_n, "bullets\n")

## frmtmb.spatial
PN <- "extensions/frmtmb.spatial/NEWS.md"
pb <- c(list(strwrap("Needs frmtmb 0.83.0: `sar()` gives core the margin of the observed rows through `frmtmb_register_autocor(margin = )`.", width = 72)),
        take("arspat2", PN, "", pick = "* "))
wr(c(emit("# frmtmb.spatial 0.5.0", list(list(name = "", blocks = pb))),
     git_show(RELBASE, PN)), file.path(REL, PN))
cat("frmtmb.spatial:", length(pb) - 1L, "bullets\n")

## frmtmb.latent: tests only
LN <- "extensions/frmtmb.latent/NEWS.md"
lb <- list(strwrap(paste(
  "Needs frmtmb 0.83.0: a test reads the escape from an `hmm()` state",
  "family's stationary point, and a test reads `update()`, which takes",
  "`na.action` from the fit. No change to the package's code."), width = 72))
wr(c(emit("# frmtmb.latent 0.7.5", list(list(name = "", blocks = lb))),
     git_show(RELBASE, LN)), file.path(REL, LN))
cat("NEWS DONE\n")
