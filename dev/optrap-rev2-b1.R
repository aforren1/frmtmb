# Reviewer of lane optrap, re-check after punch round 1: designs where
# least squares on the fixed design is a poor proxy for the residual
# scale. Per fit: code, warnings, logLik against lme4 (where it fits the
# same model) and the unit the lane chose.
#   Rscript dev/optrap-rev2-b1.R <lib or none> <case> <label> [seeds]
# cases: re10 re100 re1000 (group sd 10/100/1000, residual sd 1),
#   wide (n = 101000, 50 columns, y = k x1 + ..., k = 1e4, so that
#   n * p passes predictor_resid_sd()'s 5e6 cap), narrow (the same with
#   5 columns), bnd_lo / bnd_hi (residual sd on a grid across 0.71 and
#   1.41), t1 / t2 (student data, t with 1 and 2 df, at y x 1 and 1e3)
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
case <- a[2]
lab <- a[3]
seeds <- if (length(a) > 3) eval(parse(text = a[4])) else 1:20
cat(lab, case, "lib", dirname(find.package("frmtmb")), "\n")
units_of <- function(f) {
  pu <- tryCatch(get("predictor_units", asNamespace("frmtmb")),
                 error = function(e) NULL)
  if (is.null(pu)) return("none")
  u <- pu(f$frame)
  if (!length(u)) "1" else paste(signif(u, 3), collapse = ",")
}
one <- function(form, d, family = gaussian(), ref = TRUE) {
  w <- character()
  t0 <- proc.time()[[3]]
  f <- withCallingHandlers(suppressMessages(frm(form, family = family,
                                                data = d)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  el <- proc.time()[[3]] - t0
  l4 <- if (ref) suppressMessages(suppressWarnings(lmer(form, data = d,
    REML = FALSE, control = lmerControl(check.conv.singular = "ignore"))))
  data.frame(code = f$opt$convergence, ll = as.numeric(logLik(f)),
             dll = if (ref) as.numeric(logLik(f)) - as.numeric(logLik(l4))
               else NA,
             nwarn = length(w), unit = units_of(f), secs = el,
             warn = substr(gsub("[\t\n]", " ", paste(unique(substr(w, 1, 50)),
                                                     collapse = " | ")),
                           1, 150))
}
X <- list()
add <- function(s, setting, r) {
  X[[length(X) + 1L]] <<- cbind(case = case, seed = s, setting = setting, r)
}
for (s in seeds) {
  set.seed(s)
  if (case %in% c("re10", "re100", "re1000")) {
    S <- as.numeric(sub("re", "", case))
    d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
    d$y <- 1 + 0.5 * d$x + rnorm(20, 0, S)[d$g] + rnorm(100)
    add(s, S, one(y ~ x + (1 | g), d))
  } else if (case %in% c("wide", "narrow")) {
    n <- 101000
    p <- if (case == "wide") 49 else 4
    Z <- matrix(rnorm(n * p), n, p)
    colnames(Z) <- paste0("z", seq_len(p))
    d <- data.frame(Z, g = factor(rep(seq_len(1000), length.out = n)))
    d$y <- 1 + 1e4 * d$z1 + rnorm(1000, 0, 0.5)[d$g] + rnorm(n)
    fo <- as.formula(paste("y ~", paste(colnames(Z), collapse = " + "),
                           "+ (1 | g)"))
    add(s, p + 1, one(fo, d))
  } else if (case %in% c("bnd_lo", "bnd_hi")) {
    grid <- if (case == "bnd_lo") seq(0.66, 0.76, by = 0.01) else
      seq(1.36, 1.46, by = 0.01)
    for (sd_e in grid) {
      set.seed(s)
      d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
      d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.5 * sd_e)[d$g] +
        rnorm(100, 0, sd_e)
      add(s, sd_e, one(y ~ x + (1 | g), d))
    }
  } else if (case %in% c("t1", "t2")) {
    df <- as.numeric(sub("t", "", case))
    d0 <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
    d0$y <- 1 + 0.5 * d0$x + rnorm(20, 0, 0.5)[d0$g] + rt(100, df)
    for (sc in c(1, 1e3)) {
      d <- d0
      d$y <- d$y * sc
      r <- one(y ~ x + (1 | g), d, family = student(), ref = FALSE)
      r$ll <- r$ll + 100 * log(sc)
      add(s, sc, r)
    }
  }
}
X <- do.call(rbind, X)
write.table(X, sprintf("dev/optrap-rev2-log/b1-%s-%s.tsv", case, lab),
            sep = "\t", quote = FALSE, row.names = FALSE)
for (st in unique(X$setting)) {
  x <- X[X$setting == st, ]
  cat(sprintf(paste0("%s %s setting %s | units %s | warned %d | code != 0 ",
                     "%d | short of lme4 > 1e-4 %s | secs median %.2f\n"),
              lab, case, st, paste(unique(x$unit), collapse = " "),
              sum(x$nwarn > 0), sum(x$code != 0),
              if (all(is.na(x$dll))) "-" else
                paste(sum(x$dll < -1e-4), sprintf("(min %.3g)", min(x$dll))),
              median(x$secs)))
}
w <- X$warn[nzchar(X$warn)]
if (length(w)) print(table(substr(w, 1, 60)))
