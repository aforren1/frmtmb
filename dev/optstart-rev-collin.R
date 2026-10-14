# Reviewer of lane optstart: near-collinear designs below the rule's
# kappa bound of 1e5, with a group term, where newton_polish() does not
# run (it takes RE fits on only above predictor_design_kappa() 1e5).
#   Rscript dev/optstart-rev-collin.R <lib> <label>
# y = 2 + x1 + b_g + e, x2 = x1 + delta * e2, (1 | g), 15 groups of 20.
# Per seed and delta: the default fit, a fit from 0.70.0's start given
# through `start`, and lme4 (ML).
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
cat(a[2], "lib", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
one <- function(fo, d, start = NULL) {
  w <- character()
  f <- tryCatch(withCallingHandlers(suppressMessages(
    frm(fo, data = d, start = start)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }), error = function(e) NULL)
  list(ll = if (is.null(f)) NA else as.numeric(logLik(f)),
       code = if (is.null(f)) NA else f$opt$convergence,
       nw = length(w), w = substr(paste(unique(substr(w, 1, 30)),
                                       collapse = "|"), 1, 60))
}
rows <- list()
fo <- y ~ x1 + x2 + (1 | g)
for (s in 1:10) for (delta in c(1e-1, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4,
                                5e-5, 3e-5, 2.1e-5)) {
  set.seed(s)
  n <- 300
  d <- data.frame(x1 = rnorm(n), e = rnorm(n), g = factor(rep(1:15, 20)))
  d$y <- 2 + d$x1 + rnorm(15, 0, 1)[d$g] + rnorm(n)
  d$x2 <- d$x1 + delta * d$e
  fr <- frm(fo, data = d, dry_run = "frame")
  kap <- ns$predictor_design_kappa(fr)
  l4 <- suppressMessages(suppressWarnings(lmer(fo, data = d, REML = FALSE)))
  ll4 <- as.numeric(logLik(l4))
  df <- one(fo, d)
  od <- one(fo, d, list(beta = c(mean(d$y), 0, 0),
                        betad = log(sd(d$y))))
  rows[[length(rows) + 1L]] <- data.frame(seed = s, delta = delta,
    kappa = signif(kap, 3), gap_default = df$ll - ll4, code_d = df$code,
    nw_d = df$nw, gap_old = od$ll - ll4, code_o = od$code, nw_o = od$nw,
    w_d = df$w)
}
X <- do.call(rbind, rows)
options(width = 200)
print(X, digits = 4, row.names = FALSE)
cat("\nshort by more than 1e-4 (and silent: code 0, no warning):\n")
print(aggregate(cbind(def_short = gap_default < -1e-4,
                      def_silent = gap_default < -1e-4 & code_d == 0 &
                        nw_d == 0,
                      old_short = gap_old < -1e-4,
                      old_silent = gap_old < -1e-4 & code_o == 0 &
                        nw_o == 0) ~ delta, data = X, FUN = sum))
