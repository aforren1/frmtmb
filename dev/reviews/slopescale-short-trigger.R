# Short check of lane slopescale's fix round, item 1: attack the trigger
# `spread` (sd over the nonzero rows of a slope column) < 0.05. Designs
# where the nonzero rows have small sd but the column is well scaled
# (dose 0 or near 1), two distinct nonzero values plus zero, an
# indicator times a constant with ulp noise, and the absent cases (two
# distinct values, a 0/1 dummy, a mostly-zero N(0, 1) column). Arms:
# base default, lane default, lane FALSE, lane TRUE, in ONE process
# (namespace unloaded and reloaded; rtmb-pitfalls item 21).
# Usage: Rscript slopescale-short-trigger.R
LANE <- "C:/Users/adf44/source/r/wt-slopescale-lib"
BASE <- "C:/Users/adf44/source/r/rellib-r21"
.libPaths(c(BASE, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
use <- function(lib) {
  if ("frmtmb" %in% loadedNamespaces()) unloadNamespace("frmtmb")
  suppressPackageStartupMessages(library(frmtmb, lib.loc = lib))
}
mk <- function(seed, kind) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  n <- 300
  trt <- rep(c(TRUE, FALSE), length.out = n)
  z <- rnorm(n)
  x <- switch(kind,
    dose01 = ifelse(trt, 1 + 0.01 * z, 0),
    twonz = ifelse(trt, sample(c(1, 1.01), n, TRUE), 0),
    ulpind = ifelse(trt, sample(c(0.3, 0.1 * 3), n, TRUE), 0),
    indconst = ifelse(trt, 3, 0),
    mostzero = ifelse(runif(n) < 0.9, 0, z),
    stop(kind))
  eta <- 0.5 + rnorm(20, 0, 0.5)[g] + (0.3 + rnorm(20, 0, 0.5)[g]) * x
  data.frame(g = g, x = x, y = eta + rnorm(n))
}
kinds <- c("dose01", "twonz", "ulpind", "indconst", "mostzero")
forms <- c("y ~ x + (1 + x | g)", "y ~ 1 + (1 + x | g)")
cases <- list()
for (k in kinds) for (fo in forms) for (seed in 1:5) {
  d <- mk(seed, k)
  nz <- d$x[d$x != 0]
  cases[[length(cases) + 1L]] <- list(
    kind = k, fo = stats::as.formula(fo), fos = fo, seed = seed, d = d,
    ndist = length(unique(d$x)), spread = if (length(nz) > 1) sd(nz) else NA,
    rms = sqrt(mean(nz^2)), sdall = sd(d$x))
}
run1 <- function(cs, ctl) {
  cond <- character(0)
  f <- withCallingHandlers(
    tryCatch(frm(cs$fo, data = cs$d, control = ctl), error = function(e) e),
    warning = function(w) {
      cond <<- c(cond, paste0("W:", conditionMessage(w)))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cond <<- c(cond, paste0("M:", conditionMessage(m)))
      invokeRestart("muffleMessage")
    })
  if (inherits(f, "error")) {
    return(list(ll = NA, err = conditionMessage(f), cond = cond, as = NA,
                par = NULL))
  }
  list(ll = as.numeric(logLik(f)), par = f$opt$par, cond = cond,
       as = isTRUE(f$autoscaled), err = NA)
}
use(BASE)
rb <- lapply(cases, run1, ctl = frmtmb_control())
use(LANE)
cat("lane decide reads spread:",
    any(grepl("spread", deparse(frmtmb:::autoscale_decide))), "\n")
rl <- lapply(cases, run1, ctl = frmtmb_control())
rf <- lapply(cases, run1, ctl = frmtmb_control(autoscale = FALSE))
rt <- lapply(cases, run1, ctl = frmtmb_control(autoscale = TRUE))
nw <- function(r) sum(startsWith(r$cond, "W:"))
tab <- do.call(rbind, lapply(seq_along(cases), function(i) {
  cs <- cases[[i]]
  lls <- c(rb[[i]]$ll, rl[[i]]$ll, rf[[i]]$ll, rt[[i]]$ll)
  best <- max(lls, na.rm = TRUE)
  data.frame(kind = cs$kind, form = sub(" \\+ \\(1 \\+ x \\| g\\)", "",
                                        cs$fos),
             seed = cs$seed, ndist = cs$ndist,
             spread = signif(cs$spread, 3), rms = signif(cs$rms, 3),
             sdall = signif(cs$sdall, 3),
             as_b = rb[[i]]$as, as_l = rl[[i]]$as,
             sh_b = signif(best - lls[1], 3), sh_l = signif(best - lls[2], 3),
             sh_F = signif(best - lls[3], 3), sh_T = signif(best - lls[4], 3),
             w_b = nw(rb[[i]]), w_l = nw(rl[[i]]),
             l_eq_b = identical(rl[[i]]$par, rb[[i]]$par) &&
               identical(rl[[i]]$cond, rb[[i]]$cond),
             l_eq_F = identical(rl[[i]]$par, rf[[i]]$par),
             err = !is.na(rl[[i]]$err))
}))
options(width = 220)
print(tab, row.names = FALSE)
k <- 1e-3
cat("\nfits per arm:", nrow(tab),
    "; lane errors:", sum(tab$err),
    "; lane default engaged:", sum(tab$as_l, na.rm = TRUE),
    "; base engaged:", sum(tab$as_b, na.rm = TRUE),
    "; engaged on lane not base:", sum(tab$as_l & !tab$as_b, na.rm = TRUE), "\n")
cat("short >", k, ": base", sum(tab$sh_b > k, na.rm = TRUE),
    "; lane default", sum(tab$sh_l > k, na.rm = TRUE),
    "; lane FALSE", sum(tab$sh_F > k, na.rm = TRUE),
    "; lane TRUE", sum(tab$sh_T > k, na.rm = TRUE),
    "; lane short where base not:",
    sum(tab$sh_l > k & tab$sh_b <= k, na.rm = TRUE),
    "; max lane short", max(tab$sh_l, na.rm = TRUE), "\n")
cat("not engaged on lane: identical to base (par and conditions)",
    sum(tab$l_eq_b[!tab$as_l]), "of", sum(!tab$as_l), "\n")
for (kd in kinds) {
  s <- tab[tab$kind == kd, ]
  cat(sprintf("%-9s ndist %s spread %s..%s engaged lane %d base %d ; max short lane %.3g base %.3g FALSE %.3g ; warn lane %d base %d\n",
              kd, paste(range(s$ndist), collapse = "-"),
              signif(min(s$spread, na.rm = TRUE), 3),
              signif(max(s$spread, na.rm = TRUE), 3),
              sum(s$as_l), sum(s$as_b), max(s$sh_l), max(s$sh_b),
              max(s$sh_F), sum(s$w_l), sum(s$w_b)))
}
dif <- which(tab$as_l & !vapply(seq_along(cases), function(i)
  identical(rl[[i]]$cond, rb[[i]]$cond), NA))
for (i in dif) {
  cat("\n== conditions differ:", cases[[i]]$kind, cases[[i]]$fos,
      cases[[i]]$seed, "\n base:", substr(rb[[i]]$cond, 1, 200),
      "\n lane:", substr(rl[[i]]$cond, 1, 200), "\n")
}
