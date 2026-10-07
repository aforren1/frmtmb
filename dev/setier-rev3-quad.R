# Reviewer of lane setier, final check: quadrature fits and the boundary
# message, which se_sd_gain_up() does not check under quadrature.
#   Rscript dev/setier-rev3-quad.R <lib> [seeds]
# A  default start, binomial (10 trials) and poisson y ~ x + (1 | g), 60
#    groups of 6, group sd 0.3 and 0 (seeds 1..S): boundary messages
#    against glmer()'s isSingular() and logLik (Laplace glmer and
#    nAGQ = 25), quadrature = TRUE.
# B  the same with the sd started at exp(-12) and true sd 0.8: what the
#    user is told.
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
S <- if (length(args) > 1) as.integer(args[2]) else 20L
suppressMessages({library(frmtmb); library(lme4)})
cat("lib", find.package("frmtmb"), "\n")
gen <- function(fam, s, sd) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:60, each = 6)), x = rnorm(360))
  eta <- 0.2 + 0.5 * d$x + rnorm(60, 0, sd)[d$g]
  if (fam == "bin") {
    d$k <- rbinom(360, 10, plogis(eta)); d$nk <- 10 - d$k
    list(d = d, ff = cbind(k, nk) ~ x + (1 | g), fam = binomial())
  } else {
    d$y <- rpois(360, exp(eta))
    list(d = d, ff = y ~ x + (1 | g), fam = poisson())
  }
}
told <- function(w, m) {
  k <- c(if (any(grepl("^Optimizer did not report|gradient", w))) "CONV",
         if (any(grepl("^Standard errors are not available", w))) "SEWARN",
         if (any(grepl("^Boundary", m))) "BOUNDARY")
  if (length(k)) paste(k, collapse = "+") else "nothing"
}
for (part in c("A", "B")) for (fam in c("bin", "pois")) {
  for (sd in if (part == "A") c(0.3, 0) else 0.8) {
    tab <- character(); short <- 0; fals <- 0
    for (s in seq_len(if (part == "A") S else 5)) {
      g <- gen(fam, s, sd)
      w <- character(); m <- character()
      f <- withCallingHandlers(
        frm(g$ff, family = g$fam, data = g$d, quadrature = TRUE,
            start = if (part == "B") list(theta = -12)),
        warning = function(x) {w <<- c(w, conditionMessage(x))
          invokeRestart("muffleWarning")},
        message = function(x) {m <<- c(m, conditionMessage(x))
          invokeRestart("muffleMessage")})
      l25 <- suppressMessages(suppressWarnings(glmer(g$ff, data = g$d,
                                                      family = g$fam,
                                                      nAGQ = 25)))
      dl <- as.numeric(logLik(f)) - as.numeric(logLik(l25))
      bm <- any(grepl("^Boundary", m))
      short <- short + (dl < -1e-3)
      fals <- fals + (bm & !isSingular(l25))
      tab <- c(tab, told(w, m))
    }
    tb <- table(tab)
    cat(sprintf(paste0("%s %-4s sd %.1f | fits more than 1e-3 below ",
                       "glmer(nAGQ = 25): %d | boundary msg where glmer is ",
                       "not singular: %d | %s\n"), part, fam, sd, short, fals,
                paste(names(tb), tb, collapse = ", ")))
  }
}
