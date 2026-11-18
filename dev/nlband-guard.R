# Lane nlband: guard-absent and neighbouring cases of the hole start. Each
# ROW prints the fit at full precision (logLik, estimates, standard
# errors), its warnings and messages, and with NLBAND_ARM=after the number
# of hole starts the escape ran (HOLES). Run once per arm and compare with
# dev/nlband-cmp.R.
# NLBAND_ARM=before|after Rscript dev/nlband-guard.R
source("dev/nlband-common.R")
holes_run <- 0L
if (arm == "after") {
  trace("nl_pole_escape", exit = quote({
    if (exists("checks", inherits = FALSE)) {
      holes_run <<- holes_run +
        sum(vapply(checks, function(ck) identical(ck$kind, "hole"), NA))
    }
  }), print = FALSE, where = ns)
}
row <- function(tag, expr) {
  holes_run <<- 0L
  r <- tt(expr)
  if (iserr(r)) {
    cat("ROW", tag, "| ERROR", substr(r$msg, 1, 160), "\n")
    return(invisible(NULL))
  }
  se <- tryCatch(sqrt(diag(suppressWarnings(vcov(r)))),
                 error = function(e) NA)
  cat("ROW", tag, "| ll", format(llv(r), digits = 17), "code",
      r$opt$convergence, "beta",
      paste(format(r$estimates$beta, digits = 17), collapse = " "),
      "se", paste(format(se, digits = 12), collapse = " "),
      "warn", wsum(r), "msg", msum(r), "| HOLES", holes_run, "\n")
  for (m in attr(r, "m")) if (grepl("pole", m)) cat("   MSG", m, "\n")
  for (w in attr(r, "w")) if (grepl("pole", w)) cat("   WARN", w, "\n")
  invisible(r)
}
# 1. no hole: the band data's design with every row kept, true pole
# beyond the data (b = 0.5) and inside it (b = -2)
for (s in 1:5) {
  set.seed(s)
  d <- data.frame(x = runif(400, 0, 4))
  d$y <- 1 / (d$x + 0.5) + rnorm(400, 0, 0.2)
  row(paste("nohole clear", s), frm(band_form, data = d))
}
# 2. a hole where the inside start is worse: true pole beyond the data
for (s in 1:5) {
  d <- band_data(s)
  d$y <- 1 / (d$x + 0.5) + rnorm(nrow(d), 0, 0.2)
  row(paste("hole worse", s), frm(band_form, data = d))
}
# 3. the band data from a given start clear of the poles: kept, warned
for (s in 1:2) {
  d <- band_data(s)
  row(paste("given clear", s), frm(band_form, data = d,
                                   start = list(beta = c(-0.5, 0.4))))
  row(paste("given inside", s), frm(band_form, data = d,
                                    start = list(beta = c(1, -1.9))))
}
# 4. binomial and bernoulli band data: a pole near a row has a supremum
for (s in 1:5) {
  d <- band_data(s)
  d$nt <- 10L
  d$yb <- rbinom(nrow(d), 10, plogis(1 / (d$x - 2)))
  d$y1 <- rbinom(nrow(d), 1, plogis(2 / (d$x - 2)))
  row(paste("binomial band", s),
      frm(bf(yb | trials(nt) ~ a / (x + b), a + b ~ 1, nl = TRUE),
          family = binomial(), data = d))
  row(paste("bernoulli band", s),
      frm(bf(y1 ~ a / (x + b), a + b ~ 1, nl = TRUE),
          family = bernoulli(), data = d))
}
# 5. a saturating body with the pole in the hole: K = -2
for (s in 1:3) {
  d <- band_data(s)
  d$y <- 2 * d$x / (d$x - 2) + rnorm(nrow(d), 0, 0.2)
  row(paste("Vm*x/(K+x)", s),
      frm(bf(y ~ Vm * x / (K + x), Vm + K ~ 1, nl = TRUE), data = d))
}
# 6. two holes, the pole in the second (x without |x - 1| < 0.25 and
# without |x - 3| < 0.25, pole at 3)
for (s in 1:3) {
  set.seed(s)
  x <- runif(500, 0, 4)
  x <- x[abs(x - 1) > 0.25 & abs(x - 3) > 0.25]
  d <- data.frame(x = x, y = 1 / (x - 3) + rnorm(length(x), 0, 0.2))
  row(paste("two holes", s), frm(band_form, data = d))
}
# 7. nl_starts = FALSE turns the hole start off
d <- band_data(1)
row("nl_starts FALSE", frm(band_form, data = d,
                           control = frmtmb_control(nl_starts = FALSE)))
# 8. a K + xp body on two clusters of xp (a hole between them)
for (s in 1:3) {
  set.seed(s)
  xp <- c(runif(150, 0.1, 1.5), runif(150, 3.5, 5))
  d <- data.frame(xp = xp, y = 3 * xp / (1 + xp) + rnorm(300, 0, 0.5))
  row(paste("mm clusters", s),
      frm(bf(y ~ Vm * xp / (K + xp), Vm + K ~ 1, nl = TRUE), data = d))
}
cat("DONE\n")
