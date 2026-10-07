# Reviewer of lane setier, re-check: can se_curvature_real()'s loss-ratio
# window (3 to 5.5 for a doubled step) lose an IDENTIFIED direction whose
# likelihood is skewed or asymmetric? Designs whose tier 3 is reached
# through a unit-diagonal eigenvalue at or below 1e-9 of the largest, so
# the probe decides, against glm()/lm() on the same data.
#   Rscript dev/setier-rev2-window.R lane|base [seeds]
# Designs (seed 1..S):
#   pois   poisson, y ~ x1 + x2 with cor(x1, x2) about 1 - 1e-10, mean 0.4
#   bern   bernoulli, rare events (about 5 percent), the same pair
#   ppoly  poisson log link, raw polynomial of degree 5 on [1, 2]
#   bpoly  bernoulli logit, raw polynomial of degree 5 on [1, 2]
#   ubnear degree-5 polynomial, an INACTIVE ub on x5 one tenth of an SE
#          above its estimate (the probe's step crosses it)
#   expb   y ~ a * exp(b * x), x on [1, 1 + 1e-5] (identified, SE of b
#          about 1e4)
args <- commandArgs(TRUE)
arm <- args[1]
S <- if (length(args) > 1) as.integer(args[2]) else 20L
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", arm, find.package("frmtmb"), "\n")
ev_ratio <- function(f) {
  H <- tryCatch(f$obj$he(f$opt$par), error = function(e) NULL)
  if (is.null(H) || !all(is.finite(H))) return(NA)
  D <- sqrt(abs(diag(H)))
  if (!all(D > 0)) return(NA)
  ev <- eigen(H / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
  min(ev) / max(abs(ev))
}
rows <- list()
for (des in c("pois", "bern", "ppoly", "bpoly", "ubnear", "expb")) {
  for (s in seq_len(S)) {
    set.seed(s)
    n <- 400
    ref <- NULL
    if (des %in% c("pois", "bern")) {
      d <- data.frame(x1 = rnorm(n))
      d$x2 <- d$x1 + rnorm(n) * 1e-5
      if (des == "pois") {
        d$y <- rpois(n, exp(-1 + 0.3 * d$x1))
        fam <- poisson(); gfam <- poisson()
      } else {
        d$y <- rbinom(n, 1, plogis(-3 + 0.5 * d$x1))
        fam <- bernoulli(); gfam <- binomial()
      }
      fo <- y ~ x1 + x2
    } else if (des %in% c("ppoly", "bpoly", "ubnear")) {
      d <- data.frame(x = runif(n, 1, 2))
      for (k in 1:5) d[[paste0("x", k)]] <- d$x^k
      if (des == "ppoly") {
        d$y <- rpois(n, exp(0.5 + sin(3 * d$x)))
        fam <- poisson(); gfam <- poisson()
      } else if (des == "bpoly") {
        d$y <- rbinom(n, 1, plogis(-1 + sin(3 * d$x)))
        fam <- bernoulli(); gfam <- binomial()
      } else {
        d$y <- sin(3 * d$x) + rnorm(n, 0, 0.2)
        fam <- gaussian(); gfam <- gaussian()
      }
      fo <- y ~ x1 + x2 + x3 + x4 + x5
    } else {
      d <- data.frame(x = 1 + runif(60) * 1e-5)
      d$y <- 2 * exp(0.5 * d$x) + rnorm(60, 0, 1)
    }
    w <- character()
    catchw <- function(expr) withCallingHandlers(expr, warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(x) invokeRestart("muffleMessage"))
    f <- tryCatch(catchw(
      if (des == "expb") {
        frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), data = d,
            start = list(beta = c(2, 0.5)))
      } else if (des == "ubnear") {
        f0 <- suppressWarnings(frm(fo, data = d))
        b5 <- fixef(f0)["x5", ]
        frm(fo, data = d, prior = set_prior("", class = "b", coef = "x5",
                                            ub = b5[["Estimate"]] +
                                              0.1 * b5[["Est.Error"]]))
      } else frm(fo, family = fam, data = d)), error = function(e) e)
    if (inherits(f, "error")) {
      cat(des, s, "ERROR", conditionMessage(f), "\n"); next
    }
    se <- suppressWarnings(fixef(f))[, "Est.Error"]
    if (des != "expb") {
      g <- glm(fo, family = gfam, data = d)
      rs <- sqrt(diag(vcov(g)))
      if (des %in% c("ubnear")) rs <- rs * sqrt((n - 6) / n)
      rel <- max(abs(unname(se) / unname(rs) - 1))
    } else {
      H <- f$obj$he(f$opt$par)
      rs <- tryCatch(sqrt(diag(solve(H)))[1:2], error = function(e) NA)
      rel <- max(abs(unname(se) / unname(rs) - 1))
    }
    lost <- ns$sdr_of(f)$se_lost
    rows[[length(rows) + 1L]] <- data.frame(
      des = des, seed = s, code = f$opt$convergence, ev = ev_ratio(f),
      lost = if (length(lost)) paste(names(lost), lost, sep = ":",
                                     collapse = ",") else "",
      nlost = length(lost), nfin = sum(is.finite(se)), ncoef = length(se),
      rel = rel, warn = length(w), stringsAsFactors = FALSE)
  }
}
X <- do.call(rbind, rows)
write.table(X, paste0("dev/setier-rev2-log/window-", arm, ".tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
for (des in unique(X$des)) {
  x <- X[X$des == des, ]
  cat(sprintf(paste0("%-6s fits %d | code!=0 %d | ev ratio <= 1e-9: %d | ",
                     "fits losing an SE %d | max |SE/ref - 1| among fits ",
                     "losing none %.3g\n"),
              des, nrow(x), sum(x$code != 0), sum(x$ev <= 1e-9, na.rm = TRUE),
              sum(x$nlost > 0),
              suppressWarnings(max(x$rel[x$nlost == 0], na.rm = TRUE))))
  bad <- x[x$nlost > 0, ]
  if (nrow(bad)) {
    for (i in seq_len(min(nrow(bad), 6))) {
      cat(sprintf("   seed %d code %d ev %.2g lost %s\n", bad$seed[i],
                  bad$code[i], bad$ev[i], bad$lost[i]))
    }
  }
}
