# Review of lane slopescale, attack items 1 and 5: where the default
# engages. (a) The plan's spread is now the root mean square, not the sd:
# an uncentered slope column x = 5 + s * z has RMS near 5, so the default
# no longer engages on it, where 0.82.0 engaged on y ~ x + (1 + x | g)
# by sd(x) = s. Is the plain fit then short? (b) Ordinary small-unit
# designs (a time in years, proportions) near the 0.05 threshold: does
# the default engage, is the plain fit short there, and does the note
# fire on a fit that is not short?
# Reference: lme4 (lmer REML = FALSE) on the same data and, for frmtmb,
# the best logLik over all arms. One process, libraries swapped.
# Usage: Rscript slopescale-review-rms.R
LANE <- "C:/Users/adf44/source/r/wt-slopescale-lib"
BASE <- "C:/Users/adf44/source/r/rellib-r21"
.libPaths(c(BASE, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
use <- function(lib) {
  if ("frmtmb" %in% loadedNamespaces()) unloadNamespace("frmtmb")
  suppressPackageStartupMessages(library(frmtmb, lib.loc = lib))
}
suppressPackageStartupMessages(library(lme4))
data("sleepstudy", package = "lme4")
mk <- function(seed, s, c0, fam = "gaussian") {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  z <- rnorm(300)
  eta <- 0.5 + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, 0.5)[g] * z + 0.3 * z
  y <- if (fam == "gaussian") eta + rnorm(300) else rpois(300, exp(eta))
  data.frame(g = g, x = c0 + s * z, y = y)
}
cases <- list()
add <- function(nm, fo, d, fam = "gaussian", grp = "a") {
  cases[[length(cases) + 1L]] <<- list(nm = nm, fo = fo, d = d, fam = fam,
                                       grp = grp)
}
for (s in c(0.03, 0.01, 1e-3)) for (seed in 1:3) {
  add(sprintf("x5_xx_s%g_%d", s, seed), y ~ x + (1 + x | g), mk(seed, s, 5))
  add(sprintf("x5_1x_s%g_%d", s, seed), y ~ 1 + (1 + x | g), mk(seed, s, 5))
}
for (seed in 1:2) {
  add(sprintf("x5p_xx_s0.01_%d", seed), y ~ x + (1 + x | g),
      mk(seed, 0.01, 5, "poisson"), "poisson")
}
for (dv in c(250, 100, 365)) {
  d <- sleepstudy
  d$t <- d$Days / dv
  add(sprintf("sleep_t%g", dv), Reaction ~ t + (t | Subject), d, grp = "b")
  add(sprintf("sleep_1t%g", dv), Reaction ~ 1 + (t | Subject), d, grp = "b")
  add(sprintf("sleep_hom_t%g", dv), Reaction ~ t + homdiag(1 + t | Subject),
      d, grp = "b")
}
for (top in c(0.08, 0.1, 0.15)) for (seed in 1:2) {
  set.seed(100 + seed)
  g <- factor(rep(1:20, each = 15))
  p <- runif(300, 0, top)
  y <- 1 + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, 5)[g] * p + rnorm(300, 0, 0.3)
  d <- data.frame(g = g, p = p, y = y)
  add(sprintf("prop%g_%d", top, seed), y ~ p + (1 + p | g), d, grp = "b")
  add(sprintf("prop_hom%g_%d", top, seed), y ~ p + homdiag(1 + p | g), d,
      grp = "b")
}
fam_of <- function(f) if (f == "gaussian") gaussian() else poisson()
run1 <- function(cs, ctl) {
  cond <- character(0)
  f <- withCallingHandlers(
    tryCatch(frm(cs$fo, family = fam_of(cs$fam), data = cs$d, control = ctl),
             error = function(e) e),
    warning = function(w) {
      cond <<- c(cond, paste0("W:", class(w)[1]))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cond <<- c(cond, paste0("M:", class(m)[1], ":",
                              substr(conditionMessage(m), 1, 60)))
      invokeRestart("muffleMessage")
    })
  if (inherits(f, "error")) return(list(ll = NA, cond = cond, as = NA))
  list(ll = as.numeric(logLik(f)), cond = cond, as = isTRUE(f$autoscaled))
}
res <- list()
use(BASE)
res$base <- lapply(cases, run1, ctl = frmtmb_control())
use(LANE)
res$lane <- lapply(cases, run1, ctl = frmtmb_control())
res$laneF <- lapply(cases, run1, ctl = frmtmb_control(autoscale = FALSE))
res$laneT <- lapply(cases, run1, ctl = frmtmb_control(autoscale = TRUE))
plan_info <- lapply(cases, function(cs) {
  fr <- frm(cs$fo, family = fam_of(cs$fam), data = cs$d, dry_run = "frame")
  pl <- frmtmb:::autoscale_plan(fr)
  zs <- unlist(lapply(pl, function(p) vapply(p$z, `[[`, 0, "scale")))
  xs <- unlist(lapply(pl, `[[`, "scale"))
  dec <- !is.null(frmtmb:::autoscale_decide(fr, frmtmb_control(), FALSE))
  c(zrms = if (length(zs)) min(zs) else NA, xsd = if (length(xs)) min(xs) else NA,
    engage = dec)
})
lme <- vapply(cases, function(cs) {
  if (cs$fam != "gaussian" || grepl("homdiag", deparse(cs$fo))) return(NA_real_)
  m <- tryCatch(suppressMessages(suppressWarnings(
    lmer(cs$fo, data = cs$d, REML = FALSE))), error = function(e) NULL)
  if (is.null(m)) NA_real_ else as.numeric(logLik(m))
}, 0)
tab <- do.call(rbind, lapply(seq_along(cases), function(i) {
  ll <- c(base = res$base[[i]]$ll, lane = res$lane[[i]]$ll,
          laneF = res$laneF[[i]]$ll, laneT = res$laneT[[i]]$ll, lme4 = lme[i])
  best <- max(ll, na.rm = TRUE)
  data.frame(case = cases[[i]]$nm, grp = cases[[i]]$grp,
             zrms = signif(plan_info[[i]]["zrms"], 3),
             xsd = signif(plan_info[[i]]["xsd"], 3),
             engage = as.logical(plan_info[[i]]["engage"]),
             as_base = res$base[[i]]$as, as_lane = res$lane[[i]]$as,
             short_base = signif(best - ll[["base"]], 3),
             short_lane = signif(best - ll[["lane"]], 3),
             short_F = signif(best - ll[["laneF"]], 3),
             short_T = signif(best - ll[["laneT"]], 3),
             short_lme4 = signif(best - ll[["lme4"]], 3),
             note_lane = any(grepl("did not rescale", res$lane[[i]]$cond)),
             warn_lane = sum(startsWith(res$lane[[i]]$cond, "W:")),
             warn_base = sum(startsWith(res$base[[i]]$cond, "W:")))
}))
options(width = 250)
print(tab, row.names = FALSE)
a <- tab[tab$grp == "a", ]
cat("\n(a) uncentered x = 5 + s z: base short > 1e-4:", sum(a$short_base > 1e-4),
    "; lane short > 1e-4:", sum(a$short_lane > 1e-4), "of", nrow(a),
    "; max base", max(a$short_base), "; max lane", max(a$short_lane), "\n")
b <- tab[tab$grp == "b", ]
cat("(b) small-unit designs: lane note fired:", sum(b$note_lane),
    "; of which laneF short > 1e-4:", sum(b$note_lane & b$short_F > 1e-4),
    "; lane short > 1e-4:", sum(b$short_lane > 1e-4), "of", nrow(b), "\n")
