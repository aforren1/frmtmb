# The NEWS sections of frmtmb 0.82.0 and its changed extensions, built from
# each lane's own development section (read from the lane trees), grouped
# by subject, with one Limits section last. Rerunnable: the released
# sections are read from the base commit, never from this script's output.
#
#   Rscript dev/rel082-news.R
#
# NEWS carries no lane names, so the lane of each group is named only here:
# boundfam: Convergence and standard errors; quadseam: Quadrature; nlstop:
# Fitting (first part), Documentation; gpstop: Fitting (second part),
# Limits; ordstart: Ordinal mixtures; wcheck: Weighted likelihoods with no
# maximum; priortab: Prior tables; arspat: Residual correlation terms;
# ppdraw: Predictions; names: brms names and layouts; silent: Silent wrong
# answers fixed; sampfront: For extension packages (its one core bullet,
# which has no heading in its lane). Extensions: frmtmb.sample (priortab,
# names, sampfront), frmtmb.spatial (arspat), frmtmb.survival and
# frmtmb.latent (silent).
R <- "C:/Users/adf44/source/r"
REL <- file.path(R, "frmtmb-wt-rel082")
BASE <- "e8abed0a"
rl <- function(p) sub("\r$", "", readLines(p, warn = FALSE))
wr <- function(x, p) {
  con <- file(p, open = "wb")
  on.exit(close(con))
  writeLines(x, con, sep = "\r\n")
}
base_news <- function(f) {
  sub("\r$", "", system2("git", c("-C", REL, "show", paste0(BASE, ":", f)),
                         stdout = TRUE))
}
# a lane's development section: the lines above its first released heading
dev_part <- function(lane, f) {
  x <- rl(file.path(R, paste0("frmtmb-wt-", lane), f))
  rel <- grep("^# \\S+ [0-9]", x)[1L]
  b <- base_news(f)
  stopifnot(identical(x[rel:length(x)], b))
  x <- x[seq_len(rel - 1L)]
  x[-1L]  # the "# <pkg> (development version)" heading
}
# blocks: a run of lines between blank lines, also split where a new
# top-level bullet starts; a block keeps its own inner line breaks
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
# sections of a development part: list(heading = blocks); blocks before
# any heading go under ""
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
    else if (nzchar(nm)) out[[length(out) + 1L]] <- list(name = nm, blocks = list())
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
take <- function(lane, f, name) {
  s <- sections(dev_part(lane, f))
  hit <- Filter(function(z) identical(z$name, name), s)
  if (length(hit) != 1L) stop(lane, ": section '", name, "' found ", length(hit))
  hit[[1L]]$blocks
}

## core
core <- list(
  list(name = "Convergence and standard errors",
       blocks = take("boundfam", "NEWS.md", "Convergence and standard errors")),
  list(name = "Quadrature", blocks = take("quadseam", "NEWS.md", "Quadrature")),
  list(name = "Fitting", blocks = c(take("nlstop", "NEWS.md", "Fitting"),
                                    take("gpstop", "NEWS.md", "Fitting"))),
  list(name = "Ordinal mixtures",
       blocks = take("ordstart", "NEWS.md", "Ordinal mixtures")),
  list(name = "Weighted likelihoods with no maximum",
       blocks = take("wcheck", "NEWS.md", "Weighted likelihoods with no maximum")),
  list(name = "Prior tables and prior specifications",
       blocks = take("priortab", "NEWS.md", "Prior tables and prior specifications")),
  list(name = "Residual correlation terms",
       blocks = take("arspat", "NEWS.md", "Residual correlation terms")),
  list(name = "Predictions", blocks = take("ppdraw", "NEWS.md", "Predictions")),
  list(name = "brms names and layouts",
       blocks = take("names", "NEWS.md", "brms names and layouts")),
  list(name = "Silent wrong answers fixed",
       blocks = take("silent", "NEWS.md", "Silent wrong answers fixed")),
  list(name = "For extension packages",
       blocks = take("sampfront", "NEWS.md", "")),
  list(name = "Documentation", blocks = take("nlstop", "NEWS.md", "Documentation")),
  list(name = "Limits", blocks = take("gpstop", "NEWS.md", "Limits"))
)
# every lane's section is taken once and nothing is left behind
lanes <- c("boundfam", "quadseam", "nlstop", "gpstop", "ordstart", "wcheck",
           "priortab", "arspat", "ppdraw", "names", "silent", "sampfront")
nb_lane <- sum(vapply(lanes, function(l) {
  sum(vapply(sections(dev_part(l, "NEWS.md")), function(s) length(s$blocks), 1L))
}, 1L))
nb_out <- sum(vapply(core, function(s) length(s$blocks), 1L))
stopifnot(nb_lane == nb_out)
# wcheck's bullet 4: the word order the review asks for
fix_w <- function(b) {
  i <- grep("Under the log link the objective so$", b)
  if (length(i) == 1L) {
    b[i] <- sub("Under the log link the objective so$",
                "So under the log link, the objective falls like", b[i])
    stopifnot(grepl("^  falls like ", b[i + 1L]))
    b[i + 1L] <- sub("^  falls like ", "  ", b[i + 1L])
  }
  b
}
core[[5L]]$blocks <- lapply(core[[5L]]$blocks, fix_w)
stopifnot(any(grepl("So under the log link, the objective falls like",
                    unlist(core[[5L]]$blocks), fixed = TRUE)))
# The release review (minors 2, 3 and 7). A bullet is rewrapped at 72
# columns only when its text changes; each edit must find its bullet.
rewrap <- function(b) {
  txt <- paste(trimws(b), collapse = " ")
  strwrap(sub("^\\* ", "", txt), width = 72, initial = "* ", prefix = "  ")
}
edit_block <- function(starts, f) {
  hit <- 0L
  for (k in seq_along(core)) {
    for (j in seq_along(core[[k]]$blocks)) {
      b <- core[[k]]$blocks[[j]]
      if (startsWith(b[1L], starts)) {
        core[[k]]$blocks[[j]] <<- f(b)
        hit <- hit + 1L
      }
    }
  }
  if (hit != 1L) stop("bullet '", starts, "' found ", hit, " times")
}
# values or refusals change, so these are BREAKING (minor 2)
breaking <- c(
  "* A family's declared stationary point (the skew normal's alpha = 0) is",
  "* A prior on `sigma` that is finite at 0",
  "* A fit that stops on a bound with an objective that is not finite is",
  "* `variables()` of a fit with a `mo()` term lists the simplex,",
  "* `predict(summary = FALSE)` of one response has dimnames",
  "* `frmtmb_register_aterm()` refuses a name whose keys in `aterms` another",
  "* `conditional_effects()` refuses a `vint()`, `vreal()` or registered")
for (s in breaking) {
  edit_block(s, function(b) {
    b[1L] <- sub("^\\* ", "* BREAKING: ", b[1L])
    rewrap(b)
  })
}
# no version number in NEWS text (minor 3)
edit_block("* BREAKING: A nonlinear fit that no other check warns about", function(b) {
  txt <- paste(b, collapse = "\n")
  stopifnot(grepl("\"quieter lower maximum\" of 0.80.0.", txt, fixed = TRUE))
  rewrap(strsplit(sub("\"quieter lower maximum\" of 0.80.0.",
                      "\"quieter lower maximum\".", txt, fixed = TRUE),
                  "\n", fixed = TRUE)[[1L]])
})
# the seeded candidate of the autoscale choice takes the restart too
# (minor 7; the consolidation's cold-start rule)
edit_block("* BREAKING: When the likelihood has a free dispersion or scale", function(b) {
  txt <- paste(b, collapse = "\n")
  old <- "Only the user's own fit from the\n  default start takes the second fit;"
  stopifnot(grepl(old, txt, fixed = TRUE))
  txt <- sub(old, paste0("Only the user's own fit from the\n  default start ",
                         "takes the second fit, and under autoscale's choice ",
                         "between a fit\n  with and without its pre-fit, ",
                         "both of those take it;"), txt, fixed = TRUE)
  rewrap(strsplit(txt, "\n", fixed = TRUE)[[1L]])
})
out <- c(emit("# frmtmb 0.82.0", core), base_news("NEWS.md"))
wr(out, file.path(REL, "NEWS.md"))
cat("core: ", nb_out, " blocks in ", length(core), " sections, ",
    sum(grepl("^\\* BREAKING", unlist(lapply(core, `[[`, "blocks")))),
    " BREAKING bullets\n", sep = "")

## frmtmb.sample: one "Needs" line for the three lanes
f <- "extensions/frmtmb.sample/NEWS.md"
needs <- function(b) grepl("^Needs ", b[1L])
sb <- list()
for (l in c("names", "sampfront", "priortab")) {
  s <- sections(dev_part(l, f))
  stopifnot(length(s) == 1L, s[[1L]]$name == "")
  b <- s[[1L]]$blocks
  stopifnot(sum(vapply(b, needs, TRUE)) == 1L)
  sb <- c(sb, b[!vapply(b, needs, TRUE)])
}
need_line <- c(
  "Needs frmtmb 0.82.0, which adds `brms_variable_sort()` to the sampling",
  "API, exports `integrated_design_zero()`, and gives its prior tables the",
  "`lb` and `ub` columns and the support box below.")
wr(c("# frmtmb.sample 0.29.0", "", need_line, "",
     unlist(lapply(sb, function(b) c(b, ""))), base_news(f)),
   file.path(REL, f))
cat("frmtmb.sample: ", length(sb), " bullets\n", sep = "")

## one-lane extensions: the lane's text under the version heading
one <- function(lane, pkg, ver, need = NULL) {
  f <- file.path("extensions", pkg, "NEWS.md")
  x <- dev_part(lane, f)
  if (!is.null(need)) {
    i <- grep("^Needs the development version of frmtmb", x)
    stopifnot(length(i) == 1L)
    j <- i
    while (nzchar(x[j + 1L])) j <- j + 1L
    x <- c(x[seq_len(i - 1L)], need, x[-seq_len(j)])
  }
  while (length(x) && !nzchar(x[length(x)])) x <- x[-length(x)]
  wr(c(paste0("# ", pkg, " ", ver), x, "", base_news(f)), file.path(REL, f))
  cat(pkg, ver, "\n")
}
one("arspat", "frmtmb.spatial", "0.4.0",
    need = c("Needs frmtmb 0.82.0, which takes `mi()` beside a registered",
             "residual correlation term."))
one("silent", "frmtmb.survival", "0.3.0")
one("silent", "frmtmb.latent", "0.7.4")
cat("NEWS DONE\n")
