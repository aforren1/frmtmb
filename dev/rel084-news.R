# The NEWS sections of frmtmb 0.84.0 and its changed extensions, built from
# each lane's own development section (read from the lane trees), grouped
# by subject. Rerunnable: the released sections are read from b255350e,
# never from this script's output. Adapted from dev/rel083-news.R.
#
#   Rscript dev/rel084-news.R
#
# NEWS carries no lane names, so the lane of each group is named only here:
# squant: New family; cdfname and nbrobust: Families; alapse, slopetrig,
# sparsese and quadesc: Convergence and standard errors; nlband: Fitting;
# ordstarts: Ordinal mixtures; gpselimit: Gaussian processes; drawdir:
# Predictions. Extensions: frmtmb.eam, frmtmb.survival, frmtmb.coupling
# (cdfname). The consolidation adds quantile_elf() to asym_laplace()'s
# routes and the smoothing advice to quantile_elf()'s bullet (merge step 5),
# and names the floor 0.84.0 in the extensions' "Needs" bullets.
R <- "C:/Users/adf44/source/r"
REL <- file.path(R, "frmtmb-wt-rel084")
BASE <- "b255350e"
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
  stopifnot(grepl("\\(development version\\)$", x[1L]))
  stopifnot(identical(x[rel:length(x)], git_show(BASE, f)))
  x <- x[seq_len(rel - 1L)][-1L]
  # trim blank lines at both ends
  nz <- which(nzchar(trimws(x)))
  x[min(nz):max(nz)]
}
# "## Heading" -> lines under it (blank-trimmed)
sections <- function(x) {
  h <- grep("^## ", x)
  stopifnot(length(h) >= 1L, h[1L] == 1L)
  ends <- c(h[-1L] - 1L, length(x))
  out <- list()
  for (i in seq_along(h)) {
    body <- x[(h[i] + 1L):ends[i]]
    nz <- which(nzchar(trimws(body)))
    out[[sub("^## ", "", x[h[i]])]] <- body[min(nz):max(nz)]
  }
  out
}
# replace an exact run of lines once, or stop
swap <- function(x, old, new) {
  n <- length(old)
  hit <- which(vapply(seq_len(length(x) - n + 1L), function(i)
    identical(x[i:(i + n - 1L)], old), logical(1)))
  if (length(hit) != 1L) stop("swap: ", length(hit), " matches of ", old[1L])
  c(x[seq_len(hit - 1L)], new, x[-seq_len(hit + n - 1L)])
}

# ---- core ----
order <- list(
  "New family" = "squant",
  "Families" = c("cdfname", "nbrobust"),
  "Convergence and standard errors" = c("alapse", "slopetrig", "sparsese",
                                        "quadesc"),
  "Fitting" = "nlband",
  "Ordinal mixtures" = "ordstarts",
  "Gaussian processes" = "gpselimit",
  "Predictions" = "drawdir")
secs <- list()
lanes <- unique(unlist(order))
for (l in lanes) {
  s <- sections(dev_part(l, "NEWS.md"))
  for (nm in names(s)) {
    if (!any(vapply(order[[nm]] %||% character(), identical, logical(1), l))) {
      stop("lane ", l, " has a section not in the plan: ", nm)
    }
    secs[[nm]] <- c(secs[[nm]], if (length(secs[[nm]])) "", s[[nm]])
  }
}
# merge step 5: quantile_elf() first among asym_laplace()'s routes
secs[["Convergence and standard errors"]] <- swap(
  secs[["Convergence and standard errors"]],
  c("  regression inference\") and the warning name routes to intervals:",
    "  `frm_bootstrap()` and `frmtmb.sample::frm_sample()` read the same",
    "  likelihood, so they hold only where the data are asymmetric Laplace,",
    "  and the help no longer calls the bootstrap's intervals defensible."),
  c("  regression inference\") and the warning name routes to intervals:",
    "  first `quantile_elf()`, a smoothed asymmetric Laplace whose log",
    "  density has curvature, with `vcov_cluster()` for its standard",
    "  errors. `frm_bootstrap()` and `frmtmb.sample::frm_sample()` read the",
    "  same likelihood, so they hold only where the data are asymmetric",
    "  Laplace, and the help no longer calls the bootstrap's intervals",
    "  defensible."))
# the release review (dev/reviews/2026-11-18-rel084.md), F3 and F4 b: no
# version number as the "before"; "another link" is of mu; a bound keeps
# the fit local only where the maximum is inside it
cv <- secs[["Convergence and standard errors"]]
cv <- swap(cv,
  c("  does not cover (`cens()`, `trunc()`, another link, a nonlinear `mu`,",
    "  `me()`, a prior density on the coefficients) keeps its convergence"),
  c("  does not cover (`cens()`, `trunc()`, another link of `mu`, a",
    "  nonlinear `mu`, `me()`, a prior density on the coefficients) keeps",
    "  its convergence"))
cv <- swap(cv,
  c("  either (`set_prior(lb =, ub =)`) keeps an estimated quantile's fit",
    "  local, since the search does not read it. Such"),
  c("  either (`set_prior(lb =, ub =)`) keeps an estimated quantile's fit",
    "  local where the maximum is inside the bounds, since the search does",
    "  not read them; a fit stopped on a bound keeps its convergence",
    "  warnings. Such"))
secs[["Convergence and standard errors"]] <- cv
fm <- secs[["Families"]]
fm <- swap(fm,
  "  maximum. On 0.83.0, 5 of 6 such fits warned \"Large maximum absolute",
  "  maximum. Before, 5 of 6 such fits warned \"Large maximum absolute")
fm <- swap(fm,
  c("* The values of these densities are now correct at a large shape. In",
    "  0.83.0, `zero_inflated_negbinomial()` with a mean link that is not"),
  c("* The values of these densities are now correct at a large shape.",
    "  Before, `zero_inflated_negbinomial()` with a mean link that is not"))
secs[["Families"]] <- fm
secs[["New family"]] <- swap(
  secs[["New family"]],
  c("  cluster per row for the quantile estimate. `?quantile_elf` has a",
    "  section on how to choose between the two families. brms has no such",
    "  family."),
  c("  cluster per row for the quantile estimate. With a group-level term,",
    "  a small `smoothing` can give the Laplace objective local maxima; use",
    "  `smoothing = 1` there. `?quantile_elf` has a section on how to choose",
    "  between the two families, and the refusal of `asym_laplace()`'s",
    "  Hessian standard errors names it first. brms has no such family."))
out <- c("# frmtmb 0.84.0", "")
for (nm in names(order)) {
  out <- c(out, paste("##", nm), "", secs[[nm]], "")
}
wr(c(out, git_show(BASE, "NEWS.md")), file.path(REL, "NEWS.md"))
cat("core:", sum(grepl("^\\* ", out)), "bullets,",
    sum(grepl("^\\* BREAKING", out)), "BREAKING,", length(order), "sections\n")

# ---- extensions (lane cdfname) ----
ext <- c("frmtmb.eam" = "0.14.0", "frmtmb.survival" = "0.4.0",
         "frmtmb.coupling" = "0.6.2")
for (p in names(ext)) {
  f <- file.path("extensions", p, "NEWS.md")
  x <- dev_part("cdfname", f)
  x <- sub("Needs the frmtmb version that renames the CDF slot of",
           "Needs frmtmb 0.84.0, which renames the CDF slot of", x, fixed = TRUE)
  wr(c(paste("#", p, ext[[p]]), "", x, "", git_show(BASE, f)),
     file.path(REL, f))
  cat(p, ext[[p]], ":", sum(grepl("^\\* ", x)), "bullets\n")
}
cat("NEWS DONE\n")
