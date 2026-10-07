# Reviewer of lane setier, final check of RB2's upward check
# (se_sd_gain_up()).
#   Rscript dev/setier-rev3-short.R <lib> [seeds]
# Part 1, false "short": genuinely singular designs (group sd 0), seeds
#   1..S each, against lme4's isSingular():
#   bin60   binomial(10 trials) y ~ x + (1 | g), 60 groups of 6
#   pois60  poisson y ~ x + (1 | g), 60 groups of 6
#   both0   gaussian y ~ x + (1 + x | g), 30 groups of 6, both sds 0
#   slope0  gaussian y ~ x + (1 + x | g), intercept sd 0.7, slope sd 0
# Part 2, trapped fits: the start puts the group sd at exp(-12) where the
#   true sd is 0.8 (seeds 1..10): gaussian, binomial and poisson with the
#   Laplace approximation, binomial and poisson with quadrature = TRUE.
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
S <- if (length(args) > 1) as.integer(args[2]) else 40L
suppressMessages({library(frmtmb); library(lme4)})
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), "\n")
told <- function(w, m, lost) {
  k <- c(if (any(grepl("^Optimizer did not report|gradient", w))) "CONV",
         if (any(grepl("^Standard errors are not available", w))) "SEWARN",
         if (any(grepl("^Boundary", m))) "BOUNDARY",
         if (any(lost == "short")) "short")
  if (length(k)) paste(k, collapse = "+") else "nothing"
}
run <- function(expr) {
  w <- character(); m <- character()
  f <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(x) {
    m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
  })
  list(f = f, w = w, m = m, lost = ns$sdr_of(f)$se_lost)
}
gen <- function(des, s, sd = 0) {
  set.seed(s)
  if (des %in% c("bin60", "pois60")) {
    d <- data.frame(g = factor(rep(1:60, each = 6)), x = rnorm(360))
    eta <- 0.2 + 0.5 * d$x + rnorm(60, 0, sd)[d$g]
    if (des == "bin60") {
      d$k <- rbinom(360, 10, plogis(eta)); d$nk <- 10 - d$k
      list(d = d, ff = cbind(k, nk) ~ x + (1 | g), fam = binomial())
    } else {
      d$y <- rpois(360, exp(eta))
      list(d = d, ff = y ~ x + (1 | g), fam = poisson())
    }
  } else if (des %in% c("both0", "slope0")) {
    d <- data.frame(g = factor(rep(1:30, each = 6)), x = rnorm(180))
    s0 <- if (des == "both0") 0 else 0.7
    d$y <- 1 + 0.5 * d$x + rnorm(30, 0, s0)[d$g] + rnorm(180)
    list(d = d, ff = y ~ x + (1 + x | g), fam = gaussian())
  } else {
    d <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(120))
    d$y <- 1 + 0.5 * d$x + rnorm(20, 0, sd)[d$g] + rnorm(120)
    list(d = d, ff = y ~ x + (1 | g), fam = gaussian())
  }
}
lme <- function(g) suppressMessages(suppressWarnings(
  if (g$fam$family == "gaussian") lmer(g$ff, data = g$d, REML = FALSE)
  else glmer(g$ff, data = g$d, family = g$fam)))
cat("\n== part 1: false 'short' on singular designs\n")
for (des in c("bin60", "pois60", "both0", "slope0")) {
  tab <- character(); nshort <- 0; sing <- 0; bm <- 0
  for (s in seq_len(S)) {
    g <- gen(des, s)
    r <- run(frm(g$ff, family = g$fam, data = g$d))
    l4 <- lme(g)
    isg <- isSingular(l4)
    sing <- sing + isg
    sh <- any(r$lost == "short")
    nshort <- nshort + (sh & isg)
    bm <- bm + (isg & any(grepl("^Boundary", r$m)))
    if (sh) {
      cat(sprintf("   %s seed %d: short; lme4 singular %s; logLik - lme4 %.3g\n",
                  des, s, isg, as.numeric(logLik(r$f)) -
                    as.numeric(logLik(l4))))
    }
    tab <- c(tab, told(r$w, r$m, r$lost))
  }
  tb <- table(tab)
  cat(sprintf("%-7s lme4 singular %d | boundary msg on them %d | short on them %d | %s\n",
              des, sing, bm, nshort, paste(names(tb), tb, collapse = ", ")))
}
cat("\n== part 2: fits started with the sd at exp(-12), true sd 0.8\n")
for (cs in list(list("gaus", FALSE), list("bin60", FALSE),
                list("pois60", FALSE), list("bin60", TRUE),
                list("pois60", TRUE))) {
  tab <- character(); short <- 0; dl <- numeric()
  for (s in 1:10) {
    g <- gen(cs[[1]], s, sd = 0.8)
    nb <- 2L
    st <- list(theta = -12)
    r <- tryCatch(run(frm(g$ff, family = g$fam, data = g$d, start = st,
                          quadrature = cs[[2]])), error = function(e) e)
    if (inherits(r, "error")) {
      tab <- c(tab, paste("ERROR", substr(conditionMessage(r), 1, 40)))
      next
    }
    l4 <- lme(g)
    dl <- c(dl, as.numeric(logLik(r$f)) - as.numeric(logLik(l4)))
    tab <- c(tab, told(r$w, r$m, r$lost))
  }
  tb <- table(tab)
  cat(sprintf("%-6s quadrature %-5s | logLik - lme4 min %.3g median %.3g | %s\n",
              cs[[1]], cs[[2]], min(dl), stats::median(dl),
              paste(names(tb), tb, collapse = ", ")))
}
