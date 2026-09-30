# Reviewer, punch round 1: frm_bootstrap() and confint(method = "boot")
# on models with and without smooth / gp() / hsgp() / fs terms, base vs
# lane at the same seed. Saves dev/postfit2-rev-log/p1-boot-<arm>.rds;
# `compare` diffs them.
#   Rscript dev/postfit2-rev-p1-boot.R <base|lane|compare>
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "lane"
out <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2/dev/postfit2-rev-log"
say <- function(...) cat(sprintf(...), "\n", sep = "")
if (arm == "compare") {
  a <- readRDS(file.path(out, "p1-boot-base.rds"))
  b <- readRDS(file.path(out, "p1-boot-lane.rds"))
  for (nm in names(a)) {
    x <- a[[nm]]; y <- b[[nm]]
    if (is.character(x) || is.character(y)) {
      say("%-36s base %s | lane %s", nm, substr(paste(x), 1, 60),
          substr(paste(y), 1, 60))
      next
    }
    same <- identical(x, y)
    gap <- if (same) 0 else max(abs(x - y), na.rm = TRUE) /
      max(abs(x), na.rm = TRUE)
    rat <- if (is.matrix(x) && nrow(x) > 2) {
      sx <- apply(x, 2, sd, na.rm = TRUE); sy <- apply(y, 2, sd, na.rm = TRUE)
      paste(sprintf("%.3f", sy / sx), collapse = " ")
    } else ""
    say("%-36s identical %s; max rel gap %.3g%s", nm, same, gap,
        if (nzchar(rat)) paste0("; replicate sd lane/base ", rat) else "")
  }
  quit(save = "no")
}
libs <- c("C:/Users/adf44/source/r/wt-postfit2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages(library(frmtmb))
say("ARM %s from %s", arm, find.package("frmtmb"))
res <- list()
keep <- function(nm, expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  res[[nm]] <<- r
}
set.seed(51)
d <- data.frame(x = runif(150), z = rnorm(150), g = factor(rep(1:10, 15)))
d$y <- sin(2 * pi * d$x) + 0.5 * d$z + rnorm(10, 0, 0.5)[d$g] +
  rnorm(150, 0, 0.3)
d$cnt <- rpois(150, exp(0.3 + 0.4 * d$z))
m_lm <- frm(bf(y ~ z), family = gaussian(), data = d)
m_re <- frm(bf(y ~ z + (1 | g)), family = gaussian(), data = d)
m_po <- frm(bf(cnt ~ z + (1 | g)), family = poisson(), data = d)
m_s <- frm(bf(y ~ z + s(x)), family = gaussian(), data = d)
m_sre <- frm(bf(y ~ z + s(x) + (1 | g)), family = gaussian(), data = d)
m_fs <- frm(bf(y ~ z + s(x, g, bs = "fs", k = 5)), family = gaussian(),
            data = d)
m_gp <- tryCatch(frm(bf(y ~ z + gp(x)), family = gaussian(), data = d),
                 error = function(e) e)
m_hs <- tryCatch(frm(bf(y ~ z + gp(x, k = 10, c = 5/4)), family = gaussian(),
                     data = d), error = function(e) e)
fe <- function(f) fixef(f, flatten = TRUE)
for (nm in c("m_lm", "m_re", "m_po", "m_s", "m_sre", "m_fs", "m_gp", "m_hs")) {
  m <- get(nm)
  if (inherits(m, "error")) {
    res[[paste(nm, "fit")]] <- paste("ERROR:", conditionMessage(m))
    next
  }
  keep(paste(nm, "frm_bootstrap NA"),
       frm_bootstrap(m, nsim = 25, seed = 7)$t)
  keep(paste(nm, "frm_bootstrap NULL"),
       frm_bootstrap(m, nsim = 25, seed = 7, re_formula = NULL)$t)
  keep(paste(nm, "confint boot"),
       unclass(confint(m, method = "boot", nsim = 25, seed = 7)))
}
keep("m_re influence", unclass(influence(m_re)$fixed %||%
                                 influence(m_re)[[1]]))
saveRDS(res, file.path(out, paste0("p1-boot-", arm, ".rds")))
say("saved %d", length(res))
