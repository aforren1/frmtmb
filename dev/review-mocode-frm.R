## Reviewer probe: frmtmb on the cases of dev/review-mocode-brms.R.
##   Rscript dev/review-mocode-frm.R lane|base
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), arm, "lib:",
    dirname(dirname(system.file(package = "frmtmb"))), "\n")
qm <- function(e) suppressMessages(e)
run <- function(lab, expr) {
  w <- character(0)
  r <- withCallingHandlers(tryCatch(qm(expr), error = function(e) e),
                           warning = function(c) {
                             w <<- c(w, conditionMessage(c))
                             invokeRestart("muffleWarning")
                           })
  out <- if (inherits(r, "error")) {
    paste("ERROR:", gsub("\n", " ", conditionMessage(r)))
  } else if (is.list(r) && !is.data.frame(r)) {
    paste(vapply(names(r), function(n) paste0(n, "=", paste(
      if (is.numeric(r[[n]])) signif(r[[n]], 7) else r[[n]],
      collapse = " ")), ""), collapse = "; ")
  } else paste(signif(as.numeric(as.matrix(r)), 7), collapse = " ")
  cat(sprintf("%-34s %s\n", lab, out))
  if (length(w)) cat(sprintf("%-34s   warning: %s\n", "",
                             substr(unique(w), 1, 160)))
}
mo_entries <- function(fit) {
  lps <- fit$frame$linpreds
  out <- list()
  for (nm in names(lps)) for (mi in lps[[nm]]$mo) {
    out[[paste0(nm, ":", deparse1(mi$expr))]] <- c(
      D = mi$D, min = mi$min %||% NA, codes_head = head(mi$codes, 6))
  }
  out
}
`%||%` <- function(a, b) if (is.null(a)) b else a
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$k0 <- d$k - 1L
d$kn <- d$k - 3L
d$kd <- as.numeric(d$k)
d$ko <- factor(d$k, levels = 1:4, ordered = TRUE)
d$z <- rnorm(n)
d$g <- factor(sample(1:8, n, replace = TRUE))
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
d$y2 <- d$y + rnorm(n)
d$k2 <- d$k + 1L
d$kg <- c(1L, 2L, 5L)[(d$k %% 3L) + 1L]
cat("brms's codes of the first 6 rows of k: ",
    paste(head(d$k - 1L, 6), collapse = " "), "\n")

cat("\n== newdata on y ~ mo(k), k in 1..4\n")
fit <- qm(frm(y ~ mo(k), data = d))
nd <- function(lab, x) run(lab, fitted(fit, newdata = data.frame(k = x))[, 1])
nd("k = 1 (one row, the minimum)", 1L)
nd("k = 4 (one row, the maximum)", 4L)
nd("k = 1:4", 1:4)
nd("k = NA", NA_integer_)
nd("k = c(2, NA)", c(2L, NA))
nd("k = -1", -1L)
nd("k = 2.5", 2.5)
nd("k = 2.4", 2.4)
nd("k = 1.6", 1.6)
nd("k = 0.6 (below the minimum)", 0.6)
nd("k = 4.4 (above the maximum)", 4.4)
nd("k = 2.0 (double)", 2)
nd("k = ordered factor 2", factor(2, levels = 1:4, ordered = TRUE))
nd("k = character '2'", "2")
run("predict(k = 2.5)", predict(fit, newdata = data.frame(k = 2.5))[, 1])
run("simulate(newdata k = 2.5)",
    simulate(fit, nsim = 1, seed = 1, newdata = data.frame(k = 2.5))[[1]])

cat("\n== double fit, integer newdata\n")
fitd <- qm(frm(y ~ mo(kd), data = d))
run("kd = 2L", fitted(fitd, newdata = data.frame(kd = 2L))[, 1])

cat("\n== ordered-factor fit, integer newdata\n")
fito <- qm(frm(y ~ mo(ko), data = d))
run("ko = 2L", fitted(fito, newdata = data.frame(ko = 2L))[, 1])
run("ko = '2'", fitted(fito, newdata = data.frame(ko = "2"))[, 1])
run("ko fitted 1:4", fitted(fito, newdata = data.frame(ko = 1:4))[, 1])

cat("\n== a gap: kg in 1, 2, 5\n")
fitg <- NULL
run("fit", {fitg <- frm(y ~ mo(kg), data = d); mo_entries(fitg)})
run("fitted kg = 1:5 Estimate",
    fitted(fitg, newdata = data.frame(kg = 1:5))[, "Estimate"])
run("fitted kg = 1:5 Est.Error",
    fitted(fitg, newdata = data.frame(kg = 1:5))[, "Est.Error"])
run("fixef", fixef(fitg)[, 1:2])

cat("\n== negative: kn in -2..1\n")
fitn <- qm(frm(y ~ mo(kn), data = d))
run("kn = -2", fitted(fitn, newdata = data.frame(kn = -2L))[, 1])
run("kn = -3", fitted(fitn, newdata = data.frame(kn = -3L))[, 1])

cat("\n== the minimum on a row dropped for a missing response\n")
dm <- d
dm$y[dm$k == 1L] <- NA
fitm <- NULL
run("y NA where k == 1", {fitm <- frm(y ~ mo(k), data = dm); mo_entries(fitm)})
run("  newdata k = 1", fitted(fitm, newdata = data.frame(k = 1L))[, 1])
run("  newdata k = 2", fitted(fitm, newdata = data.frame(k = 2L))[, 1])
run("  fitted k = 2..4 vs mo(k0) on the same rows",
    fitted(fitm, newdata = data.frame(k = 2:4))[, 1] -
      fitted(qm(frm(y ~ mo(k0), data = dm)),
             newdata = data.frame(k0 = 1:3))[, 1])

cat("\n== multivariate: two responses, different minima\n")
fitmv <- NULL
run("y ~ mo(k), y2 ~ mo(k2)", {
  fitmv <- frm(bf(y ~ mo(k)) + bf(y2 ~ mo(k2)) + set_rescor(FALSE),
               data = d)
  mo_entries(fitmv)
})
run("  newdata k = 1, k2 = 2",
    fitted(fitmv, newdata = data.frame(k = 1L, k2 = 2L)))
run("  newdata k = 1, k2 = 1",
    fitted(fitmv, newdata = data.frame(k = 1L, k2 = 1L)))
run("  logLik mv", logLik(fitmv))
run("  logLik y ~ mo(k) + y2 ~ mo(k)",
    logLik(qm(frm(bf(y ~ mo(k)) + bf(y2 ~ mo(k)) + set_rescor(FALSE),
                  data = d))))

cat("\n== subset(): the minimum of the subset rows\n")
ds <- d
ds$s1 <- ds$k > 1L
ds$s2 <- TRUE
run("y | subset(s1) ~ mo(k)", mo_entries(qm(frm(
  bf(y | subset(s1) ~ mo(k)) + bf(y2 | subset(s2) ~ mo(k)) +
    set_rescor(FALSE), data = ds))))

cat("\n== dpar and nonlinear formulas\n")
fs <- NULL
run("sigma ~ mo(k2)", {fs <- frm(bf(y ~ mo(k), sigma ~ mo(k2)), data = d)
  mo_entries(fs)})
run("  logLik", logLik(fs))
run("  logLik on k0, k (shift of 0 and -1)",
    logLik(qm(frm(bf(y ~ mo(k0), sigma ~ mo(k)), data = d))))
run("  fitted sigma k2 = 2:5, k = 1:4",
    fitted(fs, newdata = data.frame(k = 1:4, k2 = 2:5), dpar = "sigma")[, 1])
run("  fitted sigma k2 = 1", fitted(fs, newdata = data.frame(k = 1L,
                                                             k2 = 1L),
                                     dpar = "sigma")[, 1])
fnl <- NULL
run("nl: a ~ mo(k)", {
  fnl <- frm(bf(y ~ a + exp(b) * z, a ~ mo(k), b ~ 1, nl = TRUE), data = d)
  mo_entries(fnl)
})
run("  nl logLik vs a ~ mo(k0)",
    as.numeric(logLik(fnl)) - as.numeric(logLik(qm(frm(
      bf(y ~ a + exp(b) * z, a ~ mo(k0), b ~ 1, nl = TRUE), data = d)))))
run("  nl fitted k = 1:4", fitted(fnl, newdata = data.frame(k = 1:4,
                                                            z = 0))[, 1])
run("  nl fitted k = 0", fitted(fnl, newdata = data.frame(k = 0L, z = 0)))
fi <- NULL
run("mo(k):z", {fi <- frm(y ~ mo(k):z, data = d); mo_entries(fi)})
run("  fitted k = 1:4, z = 1",
    fitted(fi, newdata = data.frame(k = 1:4, z = 1))[, 1])
run("  vs mo(k0):z",
    fitted(qm(frm(y ~ mo(k0):z, data = d)),
           newdata = data.frame(k0 = 0:3, z = 1))[, 1])
fre <- NULL
run("(mo(k) | g)", {fre <- frm(y ~ mo(k) + (mo(k) | g), data = d)
  mo_entries(fre)})
run("  fitted k = 1:4, g = 1",
    fitted(fre, newdata = data.frame(k = 1:4, g = "1"))[, 1])
run("  vs (mo(k0) | g)",
    fitted(qm(frm(y ~ mo(k0) + (mo(k0) | g), data = d)),
           newdata = data.frame(k0 = 0:3, g = "1"))[, 1])
run("  fitted k = 0, g = 1",
    fitted(fre, newdata = data.frame(k = 0L, g = "1"))[, 1])

cat("\n== mi() inside mo()\n")
di <- d
di$k[c(3, 7, 11)] <- NA
run("mo(mi(k))", frm(bf(y ~ mo(mi(k))) + bf(k | mi() ~ 1) +
                       set_rescor(FALSE), data = di))

cat("\n== attr(x, 'min') (brms-40)\n")
da <- d
attr(da$k, "min") <- 0L
run("fit with attr min 0", mo_entries(qm(frm(y ~ mo(k), data = da))))

cat("\n== two and one categories\n")
run("k in 1, 2", frm(y ~ mo(x), data = data.frame(x = rep(1:2, 10),
                                                  y = rnorm(20))))
run("k in 0, 1", frm(y ~ mo(x), data = data.frame(x = rep(0:1, 10),
                                                  y = rnorm(20))))
run("k all 3", frm(y ~ mo(x), data = data.frame(x = rep(3L, 10),
                                                y = rnorm(10))))
run("k logical", frm(y ~ mo(x), data = data.frame(
  x = rep(c(TRUE, FALSE), 10), y = rnorm(20))))
run("k with Inf", frm(y ~ mo(x), data = data.frame(
  x = c(Inf, rep(1:3, 6)), y = rnorm(19))))
run("k with NA (na.omit)", mo_entries(qm(frm(y ~ mo(x), data = data.frame(
  x = c(NA, 0L, rep(1:3, 6)), y = rnorm(20))))))

cat("\n== conditional_effects: the value a mo() variable is held at\n")
fz <- qm(frm(y ~ mo(k) + z, data = d))
cat("mean(k) =", mean(d$k), " round(median(k)) (brms) =",
    round(stats::median(d$k)), "\n")
ce <- qm(conditional_effects(fz, effects = "z"))[[1L]]
cat("ce held k:", format(unique(ce$k)), "\n")
run("ce estimate__ at the first z",
    ce$estimate__[1L])
run("fitted at that z with k = round(median)",
    fitted(fz, newdata = data.frame(z = ce$z[1L],
                                    k = round(stats::median(d$k))))[, 1])
run("fitted at that z with k = mean",
    fitted(fz, newdata = data.frame(z = ce$z[1L], k = mean(d$k)))[, 1])
ce2 <- qm(conditional_effects(fz, effects = "z:k"))[[1L]]
cat("ce z:k second values of k:", format(unique(ce2$k)), "\n")
cat("brms rule round((-1:1) * mad + median):",
    round((-1:1) * stats::mad(d$k) + stats::median(d$k)), "\n")
