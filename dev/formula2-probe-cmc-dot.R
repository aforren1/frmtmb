# cmc and dot expansion against brms 2.23.0's standata(), and the
# family list against brms's default_prior(). No seeds: deterministic
# data except where set.seed() is called.
args <- commandArgs(TRUE)
lib <- if (length(args) && args[1] == "before") character() else
  "C:/Users/adf44/source/r/wt-formula2-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
cat("frmtmb from", find.package("frmtmb"), "\n")
tr <- function(label, expr) {
  cat("==", label, "\n")
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  print(r)
  invisible(r)
}
strip <- function(X) {
  X <- as.matrix(X)
  attr(X, "assign") <- NULL
  attr(X, "contrasts") <- NULL
  colnames(X)[colnames(X) == "(Intercept)"] <- "Intercept"
  X
}
df <- data.frame(y = 1:10, g = rep(c("a", "b"), 5))
# standata:732 and :742
# frmtmb's cumulative() has no disc, so brms's case is run on sigma
tr("disc in cumulative", tryCatch(frm(bf(y ~ g) + lf(disc ~ 0 + g, cmc = FALSE) +
  cumulative(), data = df, dry_run = "frame"), error = conditionMessage))
b1 <- brms::bf(y ~ g) + brms::lf(sigma ~ 0 + g + (0 + g | y), cmc = FALSE)
s1 <- brms::standata(b1, df)
f1 <- frm(bf(y ~ g) + lf(sigma ~ 0 + g + (0 + g | y), cmc = FALSE),
          data = df, dry_run = "frame")
lpd <- f1$linpreds[["y.sigma"]]
tr("X_sigma equal", all.equal(strip(lpd$X), unname(s1$X_sigma),
                             check.attributes = FALSE))
tr("X_sigma colnames", list(frm = colnames(lpd$X),
                           brms = colnames(s1$X_sigma)))
cps <- Filter(function(c) c$lp_key == "y.sigma",
              unlist(lapply(f1$re_blocks, `[[`, "components"),
                     recursive = FALSE))
tr("Z disc cnms", lapply(cps, `[[`, "cnms"))
Zd <- as.matrix(f1$linpreds[["y.sigma"]]$Z)
tr("Z_1_sigma_1 brms", s1$Z_1_sigma_1)
tr("Z_1_sigma_2 in brms", "Z_1_sigma_2" %in% names(s1))
tr("Z disc row sums", rowSums(Zd))
b2 <- brms::bf(y ~ 0 + g + (1 | y), cmc = FALSE)
s2 <- brms::standata(b2, df)
f2 <- frm(bf(y ~ 0 + g + (1 | y), cmc = FALSE), data = df,
          dry_run = "frame")
tr("X equal", all.equal(strip(f2$linpreds[["y.mu"]]$X), unname(s2$X),
                        check.attributes = FALSE))
tr("X colnames", list(frm = colnames(f2$linpreds[["y.mu"]]$X),
                      brms = colnames(s2$X)))
# cmc = TRUE (default) unchanged
f3 <- frm(bf(y ~ 0 + g), data = df, dry_run = "frame")
tr("default cmc colnames", colnames(f3$linpreds[["y.mu"]]$X))
# three levels, RE with cmc = FALSE, numerically against brms
set.seed(5)
d3 <- data.frame(g = factor(rep(c("a", "b", "c"), 20)),
                 h = factor(rep(1:6, each = 10)))
d3$y <- rnorm(60, as.numeric(d3$g))
b4 <- brms::bf(y ~ 0 + g + (0 + g | h), cmc = FALSE)
s4 <- brms::standata(b4, d3)
f4 <- frm(bf(y ~ 0 + g + (0 + g | h), cmc = FALSE), data = d3,
          dry_run = "frame")
tr("3-level X", all.equal(strip(f4$linpreds[["y.mu"]]$X), unname(s4$X),
                          check.attributes = FALSE))
cp4 <- f4$re_blocks[[1]]$components[[1]]
tr("3-level RE cnms", cp4$cnms)
# the lp's Z is level-major: coefficient k of level l is column
# (l - 1) * 2 + k
Zl <- as.matrix(f4$linpreds[["y.mu"]]$Z)
tr("Z dims", dim(Zl))
tr("Z_1_1 equal", all.equal(rowSums(Zl[, seq(1, ncol(Zl), by = 2)]),
                            as.numeric(s4$Z_1_1), check.attributes = FALSE))
tr("Z_1_2 equal", all.equal(rowSums(Zl[, seq(2, ncol(Zl), by = 2)]),
                            as.numeric(s4$Z_1_2), check.attributes = FALSE))
# a fit, and a prediction on new data, with cmc = FALSE
fit4 <- frm(bf(y ~ 0 + g + (0 + g | h), cmc = FALSE), data = d3)
tr("fit4 fixef", fixef(fit4))
tr("fit4 predict newdata", predict(fit4, newdata = d3[1:4, ]))
tr("fit4 fitted equals newdata linpred",
   all.equal(frm_linpred(fit4), frm_linpred(fit4, newdata = d3)))
# same likelihood as the full-rank cmc = TRUE parameterization with
# intercept: y ~ g + (0 + g | h) spans the same X but not the same RE
fitA <- frm(bf(y ~ 0 + g, cmc = FALSE), data = d3)
fitB <- frm(bf(y ~ g), data = d3)
tr("0 + g cmc FALSE vs g: logLik", c(logLik(fitA), logLik(fitB)))

# dot expansion
dat <- data.frame(y = 1:10, x1 = 1:10, x2 = (1:10)^2, g = rep(1:2, 5))
for (f in list(y ~ ., y ~ . - x2, y ~ 0 + ., y ~ . + (1 | g),
               y ~ x1 * .)) {
  sb <- brms::standata(f, dat)
  fr <- frm(f, data = dat, dry_run = "frame")
  tr(paste("dot", deparse(f)),
     list(frm = colnames(fr$linpreds[["y.mu"]]$X),
          brms = colnames(sb$X),
          frm_formula = deparse(fr$spec$responses$y$dpars$mu$fixed),
          brms_formula = deparse(brms:::validate_formula(
            f, data = dat)$formula)))
}
sb <- brms::standata(brms::bf(y ~ x1, sigma ~ .), dat)
fr <- frm(bf(y ~ x1, sigma ~ .), data = dat, dry_run = "frame")
tr("dot sigma", list(frm = colnames(fr$linpreds[["y.sigma"]]$X),
                     brms = colnames(sb$X_sigma)))
dw <- data.frame(y = 1:10, x1 = rnorm(10), w = runif(10))
sb <- brms::standata(y | weights(w) ~ ., dw)
fr <- frm(y | weights(w) ~ ., data = dw, dry_run = "frame")
tr("dot weights", list(frm = colnames(fr$linpreds[["y.mu"]]$X),
                       brms = colnames(sb$X)))

# family list
d2 <- data.frame(y1 = rnorm(10), y2 = c(1, rep(1:3, 3)), x = rnorm(10),
                 g = rep(1:2, 5))
fl <- list(gaussian, poisson())
pb <- brms::default_prior(brms::bf(y1 ~ x) + brms::bf(y2 ~ 1), data = d2,
                          family = fl)
pf <- get_prior(bf(y1 ~ x) + bf(y2 ~ 1), data = d2, family = fl)
tr("family list brms", pb[, c("class", "coef", "resp", "dpar")])
tr("family list frm", pf[, c("class", "coef", "resp", "dpar")])
tr("family list fit", fixef(frm(bf(y1 ~ x) + bf(y2 ~ 1), data = d2,
                                family = fl)))
tr("family list length", frm(bf(y1 ~ x) + bf(y2 ~ 1), data = d2,
                             family = list(gaussian())))
tr("family list univariate", frm(bf(y1 ~ x), data = d2,
                                 family = list(gaussian())))
tr("family list keeps bf family",
   vapply(frm(bf(y1 ~ x) + bf(y2 ~ 1) + poisson(), data = d2,
              family = list(gaussian(), gaussian()),
              dry_run = "spec")$responses,
          function(r) r$family$family, ""))
