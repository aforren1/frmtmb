# Reviewer of lane optstart: edge cases of the least-squares start.
#   Rscript dev/optstart-rev-edge.R <lib> <label>
# Per case: the cold start (par_template()), the fit's log likelihood,
# code and warnings, and a reference where one exists. Run on the lane
# and on rellib-r8 and compare.
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
cat(a[2], "lib", dirname(find.package("frmtmb")), "\n")
options(width = 200)
cases <- list()
add <- function(label, fo, d, fam = gaussian(), ref = NULL) {
  cases[[length(cases) + 1L]] <<- list(label = label, fo = fo, d = d,
                                       fam = fam, ref = ref)
}
set.seed(1)
n <- 200
d0 <- data.frame(x = rnorm(n), g = factor(rep(1:20, 10)),
                 w = runif(n, 0.2, 3), o = runif(n, 0, 50))
d0$y <- 10 + 100 * d0$x + rnorm(20, 0, 2)[d0$g] + rnorm(n, 0, 1.5)
# a scale with a non-default link
add("gaussian link_sigma identity", bf(y ~ x + (1 | g)), d0,
    brmsfamily("gaussian", link_sigma = "identity"))
add("gaussian link_sigma softplus", bf(y ~ x + (1 | g)), d0,
    brmsfamily("gaussian", link_sigma = "softplus"))
add("student link_sigma identity", bf(y ~ x + (1 | g)), d0,
    brmsfamily("student", link_sigma = "identity"))
dg <- d0
dg$yg <- rgamma(n, shape = 3, rate = 3 / exp(1 + 0.8 * d0$x))
add("Gamma(log) link_shape identity", bf(yg ~ x + (1 | g)), dg,
    brmsfamily("Gamma", link = "log", link_shape = "identity"))
# sigma with its own predictor, and one without an intercept
add("sigma ~ x", bf(y ~ x + (1 | g), sigma ~ x), d0)
add("sigma ~ 0 + x", bf(y ~ x + (1 | g), sigma ~ 0 + x), d0)
# weights: unequal, with zeros, all equal to 2
add("weights unequal", bf(y | weights(w) ~ x + (1 | g)), d0)
dz <- d0
dz$w[1:20] <- 0
dz$y[1:20] <- 1e6
add("weights with zero rows holding 1e6", bf(y | weights(w) ~ x + (1 | g)),
    dz)
d2 <- d0
d2$w <- 2
add("weights all 2", bf(y | weights(w) ~ x + (1 | g)), d2)
# offsets
d3 <- d0
d3$y2 <- d0$y + d0$o * 20
add("offset term", bf(y2 ~ x + offset(20 * o) + (1 | g)), d3)
add("intercept and offset only", bf(y2 ~ 1 + offset(20 * o)), d3)
add("offset only, no intercept", bf(y2 ~ 0 + offset(20 * o)), d3)
# one row per group
d4 <- data.frame(x = rnorm(60), g = factor(1:60))
d4$y <- 3 + 2 * d4$x + rnorm(60, 0, 2) + rnorm(60)
add("one row per group", bf(y ~ x + (1 | g)), d4)
# a perfect fit, and a nearly perfect one
d5 <- data.frame(x = seq(0, 1, length.out = 50), g = factor(rep(1:5, 10)))
d5$y <- 1 + 2 * d5$x
add("perfect fit", bf(y ~ x), d5)
d5$yn <- 1 + 2 * d5$x + rnorm(50, 0, 1e-7)
add("near-perfect fit (sd 1e-7)", bf(yn ~ x), d5)
add("near-perfect fit with (1 | g)", bf(yn ~ x + (1 | g)), d5)
# n smaller than p, and n = p + 1
d6 <- as.data.frame(matrix(rnorm(10 * 14), 10, 14))
d6$y <- rnorm(10)
add("n 10 < p 15", bf(as.formula(paste("y ~", paste(names(d6)[1:14],
                                                   collapse = "+")))), d6)
d7 <- as.data.frame(matrix(rnorm(12 * 10), 12, 10))
d7$y <- d7$V1 + rnorm(12)
d7$g <- factor(rep(1:3, 4))
add("n 12, p 11, (1 | g)", bf(as.formula(paste("y ~",
    paste(names(d7)[1:10], collapse = "+"), "+ (1 | g)"))), d7)
# a response with NA rows under mi(), and weights beside it
d8 <- d0
d8$y[sample(n, 50)] <- NA
add("mi() response with NA rows", bf(y | mi() ~ x + (1 | g)), d8)
add("mi() and weights", bf(y | mi() + weights(w) ~ x + (1 | g)), d8)
# a factor with one row per level, so MAD of the residuals is 0
d9 <- data.frame(f = factor(c(1:30, rep(31, 30))), x = rnorm(60))
d9$y <- 5 + 3 * d9$x + rnorm(60)
add("student, half the residuals exactly 0", bf(y ~ f + x), d9, student())
add("huber, half the residuals exactly 0", bf(y ~ f + x), d9, huber())
# a negative binomial with nearly all zeros
d10 <- data.frame(x = rnorm(200))
d10$c <- rbinom(200, 1, 0.02) * rpois(200, 3)
add("negbinomial, 98% zeros", bf(c ~ x), d10, negbinomial())
# a multivariate model with subset(): the excluded rows of y1 hold junk
d11 <- d0
d11$s1 <- rep(c(TRUE, FALSE), length.out = n)
d11$y1 <- ifelse(d11$s1, d0$y, 1e6)
d11$y2 <- d0$y + rnorm(n)
add("subset(): excluded rows hold 1e6",
    mvbf(bf(y1 | subset(s1) ~ x), bf(y2 ~ x)) + set_rescor(FALSE), d11)
d11b <- d11
d11b$y1[!d11b$s1] <- d0$y[!d11b$s1]
add("subset(): excluded rows hold real values",
    mvbf(bf(y1 | subset(s1) ~ x), bf(y2 ~ x)) + set_rescor(FALSE), d11b)
# beta near its bounds
d12 <- data.frame(x = rnorm(300))
d12$yb <- pmin(pmax(rbeta(300, plogis(-4 + 2 * d12$x) * 5,
                          (1 - plogis(-4 + 2 * d12$x)) * 5), 1e-8),
               1 - 1e-8)
add("beta near 0", bf(yb ~ x), d12, Beta())
res <- list()
for (cs in cases) {
  st <- tryCatch(unclass(par_template(cs$fo, data = cs$d, family = cs$fam)),
                 error = function(e) e)
  stx <- if (inherits(st, "error")) paste("ERROR", conditionMessage(st)) else
    paste(signif(c(st$beta[1:min(2, length(st$beta))], st$betad), 5),
          collapse = " ")
  w <- character()
  f <- tryCatch(withCallingHandlers(suppressMessages(
    frm(cs$fo, data = cs$d, family = cs$fam)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  res[[length(res) + 1L]] <- data.frame(case = cs$label,
    start = substr(stx, 1, 60),
    ll = if (inherits(f, "error")) NA else as.numeric(logLik(f)),
    code = if (inherits(f, "error")) NA else f$opt$convergence,
    nw = length(w),
    note = substr(if (inherits(f, "error")) {
      paste("ERROR", conditionMessage(f))
    } else paste(unique(substr(w, 1, 45)), collapse = " | "), 1, 90))
}
R <- do.call(rbind, res)
print(R, digits = 12, row.names = FALSE, right = FALSE)
# references
cat("\nreferences:\n")
cat(sprintf("lme4 (1|g) ML, d0: %.10f\n",
    as.numeric(logLik(lmer(y ~ x + (1 | g), d0, REML = FALSE)))))
cat(sprintf("lme4 weighted, d0 (lme4 weights scale the residual variance): %.10f\n",
    as.numeric(logLik(lmer(y ~ x + (1 | g), d0, weights = w,
                           REML = FALSE)))))
cat(sprintf("lme4 offset, d3: %.10f\n",
    as.numeric(logLik(lmer(y2 ~ x + offset(20 * o) + (1 | g), d3,
                           REML = FALSE)))))
cat(sprintf("lm one row per group (sigma^2 = total), d4: %.10f\n",
    as.numeric(logLik(lm(y ~ x, d4)))))
ls8 <- lm(y ~ x, d8)
cat(sprintf("lm on the observed rows of d8 (no group): %.10f\n",
            as.numeric(logLik(ls8))))
cat(sprintf("lmer on observed rows of d8: %.10f\n",
    as.numeric(logLik(lmer(y ~ x + (1 | g), d8[!is.na(d8$y), ],
                           REML = FALSE)))))
