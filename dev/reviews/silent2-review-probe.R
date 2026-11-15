# Review of lane silent2: (a) the probit fit whose threshold sits at -37
# with latent SE 38.5 (silent); (b) fixed-only factor(cyl) under helmert.
.libPaths(c("C:/Users/adf44/source/r/wt-silent2-lib", "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
gen <- function(seed, n = 300, K = 4, b1 = 1.5, b2 = -0.8, c1 = 1, c2 = -1,
                p = 0.4, noise = "logit") {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, p)
  e <- switch(noise, logit = rlogis(n), probit = rnorm(n))
  lat <- ifelse(cls == 1, b1 * d$x + c1, b2 * d$x + c2) + e
  cuts <- stats::quantile(lat, seq_len(K - 1) / K)
  d$y <- 1L + rowSums(outer(lat, cuts, ">"))
  d
}
d <- gen(1, n = 80, K = 4, noise = "probit")
w <- character()
f <- withCallingHandlers(frm(bf(y ~ x), family = mixture(cumulative("probit"), cumulative("probit")), data = d),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") },
  message = function(m) invokeRestart("muffleMessage"))
cat("warnings:", length(w), "\n")
print(suppressWarnings(fixef(f)))
print(ns$mixture_ord_degeneracy(f, "y"))
cat("latent SE:", ns$mixture_ord_thres_se(f, "y"), "\n")
# category 1's probability in component 2 at every row
dp <- ns$eval_dpars(f)$y
eta2 <- rep_len(as.numeric(dp$mu2), nrow(d))
fx <- suppressWarnings(fixef(f))
t1 <- fx["mu2_Intercept[1]", "Estimate"]
cat("component 2, P(y = 1) over rows: max", format(max(pnorm(t1 - eta2)), digits = 3), "\n")
cat("rows with y = 1:", sum(d$y == 1), "\n")

# (b) fixed-only inline factor() under helmert at predict
set.seed(7)
dc <- data.frame(g = factor(rep(1:10, each = 9)), cyl = sample(c(4, 6, 8), 90, TRUE))
dc$y <- dc$cyl / 4 + rnorm(90)
fs <- withr::with_options(list(contrasts = c("contr.sum", "contr.poly")), {
  ff <- suppressMessages(frm(y ~ factor(cyl), data = dc))
  list(f = ff, p = fitted(ff, newdata = dc[1:5, ])[, "Estimate"])
})
ph <- withr::with_options(list(contrasts = c("contr.helmert", "contr.poly")),
                          fitted(fs$f, newdata = dc[1:5, ])[, "Estimate"])
cat("fixed-only factor(cyl) sum -> helmert, max |diff|:", max(abs(ph - fs$p)), "\n")
fr <- withr::with_options(list(contrasts = c("contr.sum", "contr.poly")),
  suppressWarnings(suppressMessages(frm(y ~ 1 + (1 + factor(cyl) | g), data = dc))))
p0 <- withr::with_options(list(contrasts = c("contr.sum", "contr.poly")),
                          fitted(fr, newdata = dc[1:5, ])[, "Estimate"])
p1 <- withr::with_options(list(contrasts = c("contr.helmert", "contr.poly")),
                          fitted(fr, newdata = dc[1:5, ])[, "Estimate"])
cat("RE-only (1 + factor(cyl) | g) sum -> helmert, max |diff|:", max(abs(p1 - p0)),
    "; in-sample vs fit-time newdata:", max(abs(unname(fitted(fr)[1:5, "Estimate"]) - unname(p0))), "\n")
