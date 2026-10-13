# Reviewer of lane optrap, consolidation check: newton_polish()'s gate,
# a step kept only if it gains at least 0.1 of its prediction. The gate
# is patched in this process to G (0.1 as built, 0.01, 0 = any
# decrease), and every attempted step's actual gain over its prediction
# is recorded. Designs, seeds 1 to 20: raw polynomials of degree 5 to 8
# on x in [1, 2] (dev/optrap-p1-polystudy.R's data), without and with
# (1 | g), y x 1 and 16; the y ~ big + (1 | g) design of
# dev/optrap-rev4-big.R. Short: more than 1e-4 below lm() or lmer().
#   Rscript dev/optrap-rev4-gate.R <G> <label>
a <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
G <- as.numeric(a[1])
lab <- a[2]
ns <- asNamespace("frmtmb")
f0 <- get("newton_polish", ns)
src <- deparse(f0)
i <- grep("v - vq < 0.1 * gain", src, fixed = TRUE)
stopifnot(length(i) == 1L)
src[i] <- sub("v - vq < 0.1 * gain",
              "{ .rev4_ratio(v, vq, gain); v - vq < .rev4_G * gain }",
              src[i], fixed = TRUE)
ratios <- numeric(0)
pe <- new.env(parent = ns)
pe$.rev4_G <- G
pe$.rev4_ratio <- function(v, vq, gain) {
  if (is.finite(vq)) ratios <<- c(ratios, (v - vq) / gain)
}
f1 <- eval(parse(text = src))
environment(f1) <- pe
utils::assignInNamespace("newton_polish", f1, "frmtmb")
cat(lab, "gate", G, "lib", dirname(find.package("frmtmb")), "\n")
X <- list()
one <- function(design, s, fo, d, ref) {
  ratios <<- numeric(0)
  w <- character()
  f <- withCallingHandlers(suppressMessages(frm(fo, data = d)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  X[[length(X) + 1L]] <<- data.frame(
    design = design, seed = s, code = f$opt$convergence,
    gap = as.numeric(logLik(f)) - ref, nwarn = length(w),
    steps = f$opt$polish[["steps"]] %||% 0,
    ratios = paste(signif(ratios, 3), collapse = ","))
}
for (deg in 5:8) for (re in c(FALSE, TRUE)) for (yk in c(1, 16)) {
  for (s in 1:20) {
    set.seed(1000 * deg + s)
    d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
    d$y <- (sin(3 * d$x) + (if (re) rnorm(20, 0, 0.3)[d$g] else 0) +
              rnorm(200, 0, 0.2)) * yk
    for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
    rhs <- paste0("x", 1:deg, collapse = " + ")
    fo <- as.formula(paste("y ~", rhs, if (re) "+ (1 | g)"))
    # lmer fails on some degree-8 data ("Downdated VtV"); no reference
    ref <- if (re) {
      tryCatch(as.numeric(logLik(suppressMessages(suppressWarnings(
        lmer(fo, data = d, REML = FALSE))))), error = function(e) NA_real_)
    } else as.numeric(logLik(lm(fo, data = d)))
    one(sprintf("deg%d%s_y%g", deg, if (re) "re" else "", yk), s, fo, d, ref)
  }
}
for (s in 1:20) {
  set.seed(s)
  n <- 400
  d <- data.frame(g = factor(rep(1:40, each = 10)), big = 1e5 + runif(n))
  d$y2 <- 0.5 * (d$big - 1e5) + rnorm(40, 0, 0.5)[d$g] + rnorm(n)
  fo <- y2 ~ big + (1 | g)
  ref <- as.numeric(logLik(suppressMessages(suppressWarnings(
    lmer(fo, data = d, REML = FALSE)))))
  one("big", s, fo, d, ref)
}
X <- do.call(rbind, X)
write.table(X, sprintf("dev/optrap-rev4-log/gate-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
for (ds in unique(X$design)) {
  x <- X[X$design == ds, ]
  sh <- !is.na(x$gap) & x$gap < -1e-4
  cat(sprintf("%-14s | short %2d (silent %2d, min %9.3g) | stepped %2d\n",
              ds, sum(sh), sum(sh & x$code == 0 & x$nwarn == 0),
              suppressWarnings(min(x$gap, na.rm = TRUE)), sum(x$steps > 0)))
}
r <- as.numeric(unlist(strsplit(X$ratios[nzchar(X$ratios)], ",")))
cat("attempted steps:", length(r), "| gain / prediction below 0.1:",
    sum(r < 0.1), "| in [0, 0.1):", sum(r >= 0 & r < 0.1), "\n")
print(round(quantile(r, c(0, 0.01, 0.05, 0.5, 0.95, 1)), 3))
