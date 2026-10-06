# Reviewer: lane fixes' flat-case table (dev/fixes-p2-table.R in the
# fixes worktree, body copied below this header by
# dev/nanse-rev-fixes-build.sh), run on the trial merge: per model, every
# warning frm() raises, classified, so a fit warned twice shows up.
#   Rscript dev/nanse-rev-fixes-table.R merge|release
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  release = "C:/Users/adf44/source/r/rellib-r6",
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
kind <- function(m) {
  if (grepl("^Standard errors are not available", m)) return("SE")
  if (grepl("are not identified: at the optimum", m)) return("NLFLAT")
  if (grepl("^Optimizer did not report|gradient", m)) return("CONV")
  if (grepl("^Some standard errors are not finite", m)) return("VCOV")
  if (grepl("^Hessian is not positive", m)) return("PDHESS")
  "OTHER"
}
run <- function(lab, expr, truth) {
  w <- character()
  r <- tryCatch({
    fit <- suppressMessages(withCallingHandlers(expr, warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }))
    wv <- character()
    se <- withCallingHandlers(fixef(fit)[, "Est.Error"],
      warning = function(x) {
        wv <<- c(wv, conditionMessage(x)); invokeRestart("muffleWarning")
      })
    k <- vapply(w, kind, "")
    kv <- vapply(wv, kind, "")
    sprintf("%-4s conv %d  frm(): %-14s fixef(): %-6s finite SE %d/%d",
            truth, fit$opt$convergence,
            if (length(k)) paste(k, collapse = "+") else "none",
            if (length(kv)) paste(kv, collapse = "+") else "none",
            sum(is.finite(se)), length(se))
  }, error = function(e) paste(truth, "ERROR:",
                               substr(conditionMessage(e), 1, 100)))
  cat(sprintf("%-46s %s\n", lab, r))
}
n <- 60
dd <- data.frame(x = rnorm(n), z = rnorm(n))
dd$y <- 3 + 0.5 * dd$x + 0.3 * dd$z + rnorm(n, 0, 0.4)
cat("== unidentified\n")
run("a + b, shared 1 + x", frm(bf(y ~ a + b, a ~ 1 + x, b ~ 1 + x,
                                  nl = TRUE), data = dd), "UNID")
run("a - b", frm(bf(y ~ a - b, a ~ 1 + x, b ~ 1, nl = TRUE), data = dd),
    "UNID")
run("a + 2 * b", frm(bf(y ~ a + 2 * b, a ~ 1 + x, b ~ 1, nl = TRUE),
                     data = dd), "UNID")
run("partial overlap", frm(bf(y ~ a + b, a ~ 1 + x, b ~ 0 + x + z,
                              nl = TRUE), data = dd), "UNID")
run("a * b intercepts", frm(bf(y ~ a * b, a ~ 1, b ~ 1, nl = TRUE),
                            data = dd, start = list(beta = c(1, 1))), "UNID")
run("a + exp(b), b ~ 1", frm(bf(y ~ a + exp(b), a ~ 1 + x, b ~ 1, nl = TRUE),
                             data = dd, start = list(beta = c(1, 0, -1))),
    "UNID")
set.seed(77)
G <- 8
dn <- data.frame(g = factor(rep(seq_len(G), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(G, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
run("c0 + exp(a)^k, a ~ 1 + (1 | g)",
    frm(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1,
           nl = TRUE), data = dn), "UNID")
k <- matrix(c(1,1,1,1,0,0,1,1,0,1,0,0,1,0,0,1,0,1,0,0,
              0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
              0,0,1,1,0,1,1,0,0,0,0,1,0,0,0,0,0,1,0,0,
              0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,
              1,0,1,1,0,1,1,1,0,1,0,0,1,0,0,0,0,1,0,0,
              1,1,0,1,0,0,0,1,0,1,0,1,1,0,0,1,0,1,0,0,
              0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,1,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,1,
              1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0),
            nrow = 10, byrow = TRUE)
db <- data.frame(k = as.vector(k), person = factor(rep(1:10, times = 20)),
                 question = factor(rep(1:20, each = 10)))
run("BCM plogis(lp) * plogis(lq)",
    frm(bf(k ~ plogis(lp) * plogis(lq), lp ~ 0 + person, lq ~ 0 + question,
           nl = TRUE), family = bernoulli(link = "identity"), data = db),
    "UNID")
cat("== identified\n")
run("a + b, disjoint", frm(bf(y ~ a + b, a ~ 1 + x, b ~ 0 + z, nl = TRUE),
                           data = dd), "ID")
run("a * exp(b * x)", frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE),
                          data = dd, start = list(beta = c(3, 0))), "ID")
run("a + exp(b), b ~ 1 + x", frm(bf(y ~ a + exp(b), a ~ 1, b ~ 1 + x,
                                    nl = TRUE), data = dd,
                                 start = list(beta = c(1, 0.5, 0.1))), "ID")
set.seed(7)
nn <- 3000
pt <- runif(nn, 0, 0.6)
f1 <- pnorm(pt, 0.2, 0.045); f2 <- pnorm(pt, 0.32, 0.045)
pp <- f2 * 0.95 + f1 * (1 - f2) * 0.2 + (1 - f1) * (1 - f2) * 0.5
d1 <- data.frame(pt = pt, y = rbinom(nn, 1, pp))
run("B1 response preparation",
    frm(bf(y ~ pnorm(pt, m2, exp(ls)) * 0.95 +
             pnorm(pt, m1, exp(ls)) * (1 - pnorm(pt, m2, exp(ls))) * 0.2 +
             (1 - pnorm(pt, m1, exp(ls))) * (1 - pnorm(pt, m2, exp(ls))) * 0.5,
           m1 ~ 1, m2 ~ 1, ls ~ 1, nl = TRUE),
        family = bernoulli(link = "identity"), data = d1,
        start = list(beta = c(0.15, 0.4, -3))), "ID")
set.seed(11)
x <- runif(800, -2, 2)
d2 <- data.frame(x = x, y = rbinom(800, 1, pnorm(x, 0.3, 0.7)))
run("B1 psychometric", frm(bf(y ~ pnorm(x, m0, exp(ls)), m0 ~ 1, ls ~ 1,
                              nl = TRUE), family = bernoulli(link = "identity"),
                           data = d2, start = list(beta = c(0, 0))), "ID")
set.seed(12)
x <- runif(300, -3, 3)
d3 <- data.frame(x = x, y = 2 * dnorm(x, 0.5, 1.2) + rnorm(300, 0, 0.05))
run("B1 bump", frm(bf(y ~ h * dnorm(x, c0, exp(ls)), h ~ 1, c0 ~ 1, ls ~ 1,
                      nl = TRUE), data = d3, start = list(beta = c(1, 0, 0))),
    "ID")
set.seed(22)
yr <- seq(1900, 2000, length.out = 120)
d4 <- data.frame(yr = yr, y = 50 / (1 + exp((1950 - yr) / 12)) +
                   rnorm(120, 0, 1.5))
run("B2 logistic growth, years", frm(bf(y ~ Asym / (1 + exp((xmid - yr) /
                                                           exp(lscal))),
                                        Asym ~ 1, xmid ~ 1, lscal ~ 1,
                                        nl = TRUE), data = d4,
                                     start = list(beta = c(50, 1950,
                                                           log(12)))), "ID")
set.seed(21)
x <- runif(800, 200, 400)
d5 <- data.frame(x = x, y = rbinom(800, 1, 0.04 * 0.5 + 0.96 *
                                     pnorm(x, 300, 25)))
run("B2 lapse psychometric, x 200..400",
    frm(bf(y ~ lapse * 0.5 + (1 - lapse) * pnorm(x, m0, exp(ls)),
           lapse ~ 1, m0 ~ 1, ls ~ 1, nl = TRUE),
        family = bernoulli(link = "identity"), data = d5,
        start = list(beta = c(0.05, 300, log(25)))), "ID")
set.seed(21)
x <- runif(600, 50, 150)
d6 <- data.frame(x = x, y = rbinom(600, 1, pnorm(x, 100, 12)))
run("psychometric x 50..150", frm(bf(y ~ pnorm(x, m0, exp(ls)), m0 ~ 1,
                                     ls ~ 1, nl = TRUE),
                                  family = bernoulli(link = "identity"),
                                  data = d6,
                                  start = list(beta = c(100, log(12)))), "ID")
set.seed(23)
tt <- runif(150, 1000, 5000)
d7 <- data.frame(tt = tt, y = 10 * exp(-5e-4 * tt) + rnorm(150, 0, 0.2))
run("decay tt 1000..5000", frm(bf(y ~ y0 * exp(-exp(lk) * tt), y0 ~ 1,
                                  lk ~ 1, nl = TRUE), data = d7,
                               start = list(beta = c(10, log(5e-4)))), "ID")
set.seed(25)
dose <- rep(c(100, 200, 400, 600, 800, 1000), each = 20)
d8 <- data.frame(dose = dose, y = 2 + 10 * dose^3 / (400^3 + dose^3) +
                   rnorm(120, 0, 0.8))
run("Emax dose 100..1000", frm(bf(y ~ e0 + emax * dose^3 / (exp(led50)^3 +
                                                             dose^3),
                                  e0 ~ 1, emax ~ 1, led50 ~ 1, nl = TRUE),
                               data = d8, start = list(beta = c(2, 10,
                                                                log(400)))),
    "ID")
set.seed(24)
d9 <- data.frame(x = runif(100, 1e3, 2e3), z = rnorm(100))
d9$y <- 1 + 0.002 * d9$x + 0.5 * d9$z + rnorm(100, 0, 0.3)
run("a + b, covariate 1e3..2e3", frm(bf(y ~ a + b, a ~ 1 + x, b ~ 0 + z,
                                        nl = TRUE), data = d9), "ID")
d10 <- data.frame(x = c(1, 2, 3, 4), y = c(2.1, 3.9, 8.2, 15.8))
run("n = 4, a * exp(b * x)", frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1,
                                    nl = TRUE), data = d10,
                                 start = list(beta = c(1, 0.7))), "ID")
set.seed(31)
d11 <- data.frame(x = rnorm(200))
d11$x2 <- d11$x + rnorm(200, 0, 0.03)
d11$y <- 1 + d11$x - 0.5 * d11$x2 + rnorm(200, 0, 0.5)
run("near-collinear a ~ 1 + x, b ~ 0 + x2 (cor 0.9995)",
    frm(bf(y ~ a + b, a ~ 1 + x, b ~ 0 + x2, nl = TRUE), data = d11), "ID")
set.seed(505)
ng <- 40
g <- factor(rep(seq_len(ng), each = 10))
wg <- rnorm(ng)
d12 <- data.frame(g = g, x = rnorm(400), w = wg[g])
d12$y <- 1 + 0.5 * d12$x + rnorm(ng, 0, exp(-0.5 + 0.5 * wg))[g] +
  rnorm(400, 0, 0.7)
run("drmTMB b0 + exp(lsd) * zz", frm(bf(y ~ b0 + exp(lsd) * zz, b0 ~ x,
                                        lsd ~ 0 + w, zz ~ 0 + (1 | g),
                                        nl = TRUE), data = d12), "ID")
set.seed(26)
d13 <- data.frame(f = factor(sample(c("p", "q", "r"), 90, TRUE)),
                  x = runif(90))
d13$y <- c(p = 1, q = 2, r = 3)[as.character(d13$f)] * exp(0.4 * d13$x) +
  rnorm(90, 0, 0.1)
d13 <- d13[d13$f != "r", ]
run("absent level r, a ~ 0 + f, a * exp(b * x)",
    frm(bf(y ~ a * exp(b * x), a ~ 0 + f, b ~ 1, nl = TRUE), data = d13,
        start = list(beta = c(1, 1, 0))), "ID")
