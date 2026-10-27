## Reviewer probe: brms 2.23.0's mo() coding on edge cases the lane did
## not run. standata() and brm(empty = TRUE) only; nothing compiles.
##   Rscript dev/review-mocode-brms.R > dev/review-mocode-log/brms.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
mo_fields <- function(s) {
  nm <- grep("^(Xmo|Jmo)", names(s), value = TRUE)
  paste(vapply(nm, function(n) paste0(n, "=", paste(s[[n]], collapse = " ")),
               ""), collapse = "; ")
}
sd_show <- function(lab, expr) {
  r <- tryCatch(q(expr), error = function(e) e)
  cat(sprintf("%-34s %s\n", lab, if (inherits(r, "error")) {
    paste("ERROR:", gsub("\n", " ", conditionMessage(r)))
  } else mo_fields(r)))
}
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

cat("\n== newdata on an empty fit of y ~ mo(k), k in 1..4\n")
fit <- q(brm(y ~ mo(k), data = d, empty = TRUE))
nd <- function(lab, x) {
  sd_show(lab, standata(fit, newdata = data.frame(k = x, y = 0)))
}
nd("k = 1 (one row, the minimum)", 1L)
nd("k = 4 (one row, the maximum)", 4L)
nd("k = NA", NA_integer_)
nd("k = c(2, NA)", c(2L, NA))
nd("k = -1", -1L)
nd("k = 2.5", 2.5)
nd("k = 2.0 (double)", 2)
nd("k = ordered factor 2", factor(2, levels = 1:4, ordered = TRUE))
nd("k = character '2'", "2")

cat("\n== integer fit on a double column, integer newdata\n")
fitd <- q(brm(y ~ mo(kd), data = d, empty = TRUE))
sd_show("kd = 2L", standata(fitd, newdata = data.frame(kd = 2L, y = 0)))

cat("\n== ordered-factor fit, integer newdata\n")
fito <- q(brm(y ~ mo(ko), data = d, empty = TRUE))
sd_show("ko = 2L", standata(fito, newdata = data.frame(ko = 2L, y = 0)))
sd_show("ko = '2'", standata(fito, newdata = data.frame(ko = "2", y = 0)))

cat("\n== a gap: kg in 1, 2, 5; newdata at the empty categories\n")
fitg <- q(brm(y ~ mo(kg), data = d, empty = TRUE))
sd_show("fit", standata(fitg))
sd_show("kg = 3, 4", standata(fitg, newdata = data.frame(kg = 3:4, y = 0)))

cat("\n== negative: kn in -2..1\n")
fitn <- q(brm(y ~ mo(kn), data = d, empty = TRUE))
sd_show("kn = -2", standata(fitn, newdata = data.frame(kn = -2L, y = 0)))
sd_show("kn = -3", standata(fitn, newdata = data.frame(kn = -3L, y = 0)))

cat("\n== the minimum on a row dropped for a missing response\n")
dm <- d
dm$y[dm$k == 1L] <- NA
sd_show("y NA where k == 1", standata(y ~ mo(k), data = dm))
fitm <- q(brm(y ~ mo(k), data = dm, empty = TRUE))
sd_show("  newdata k = 1", standata(fitm, newdata = data.frame(k = 1L, y = 0)))
sd_show("  newdata k = 2", standata(fitm, newdata = data.frame(k = 2L, y = 0)))

cat("\n== multivariate: two responses, different minima\n")
sd_show("y ~ mo(k), y2 ~ mo(k2)",
        standata(bf(y ~ mo(k)) + bf(y2 ~ mo(k2)) + set_rescor(FALSE),
                 data = d))
fitmv <- q(brm(bf(y ~ mo(k)) + bf(y2 ~ mo(k2)) + set_rescor(FALSE),
               data = d, empty = TRUE))
sd_show("  newdata k = 1, k2 = 2",
        standata(fitmv, newdata = data.frame(k = 1L, k2 = 2L, y = 0,
                                             y2 = 0)))
sd_show("  newdata k = 1, k2 = 1",
        standata(fitmv, newdata = data.frame(k = 1L, k2 = 1L, y = 0,
                                             y2 = 0)))

cat("\n== subset(): the minimum of the subset or of all rows?\n")
ds <- d
ds$s1 <- ds$k > 1L
ds$s2 <- TRUE
sd_show("y | subset(s1) ~ mo(k)",
        standata(bf(y | subset(s1) ~ mo(k)) + bf(y2 | subset(s2) ~ mo(k)) +
                   set_rescor(FALSE), data = ds))

cat("\n== dpar and nonlinear formulas\n")
sd_show("sigma ~ mo(k2)", standata(bf(y ~ mo(k), sigma ~ mo(k2)), data = d))
sd_show("nl: a ~ mo(k)",
        standata(bf(y ~ a + exp(b) * z, a ~ mo(k), b ~ 1, nl = TRUE),
                 data = d, prior = prior(normal(0, 1), nlpar = "a") +
                   prior(normal(0, 1), nlpar = "b")))
sd_show("mo(k):z", standata(y ~ mo(k):z, data = d))
sd_show("(mo(k) | g)", standata(y ~ mo(k) + (mo(k) | g), data = d))

cat("\n== mi() inside mo()\n")
di <- d
di$k[c(3, 7, 11)] <- NA
sd_show("mo(mi(k))", standata(bf(y ~ mo(mi(k))) + bf(k | mi() ~ 1) +
                                set_rescor(FALSE), data = di))

cat("\n== brms-40: attr(x, 'min')\n")
set.seed(1)
da <- data.frame(x = c(1, 2, 3, 4, 1, 2), y = rnorm(6))
attr(da$x, "min") <- 0
fa <- q(brm(y ~ mo(x), data = da, empty = TRUE))
sd_show("fit with attr min 0", standata(fa))
sd_show("newdata x = 1:4",
        standata(fa, newdata = data.frame(x = 1:4, y = 0)))
cat("attr(fa$data$x, 'min') =", format(attr(fa$data$x, "min")), "\n")
db <- data.frame(x = c(1, 2, 3, 4, 1, 2), y = rnorm(6))
attr(db$x, "min") <- 2
sd_show("fit with attr min 2 (above data)", standata(y ~ mo(x), data = db))

cat("\n== two and one categories\n")
sd_show("k in 1, 2", standata(y ~ mo(x), data = data.frame(
  x = rep(1:2, 10), y = rnorm(20))))
sd_show("default_prior k in 1, 2", {
  p <- default_prior(y ~ mo(x), data = data.frame(x = rep(1:2, 10),
                                                    y = rnorm(20)))
  list(Jmo = paste(p$class, p$coef, p$prior, collapse = " | "))
})
