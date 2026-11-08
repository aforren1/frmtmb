# Probe of items 2, 3, 5, 7 and 8 on one build (READERS_ARM). Seeds inline.
source("dev/readers-setup.R")
q <- function(e) tryCatch(suppressWarnings(e), error = function(err) {
  cat("ERROR:", conditionMessage(err), "\n"); NULL
})
set.seed(7971)
n <- 400
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.4)
d$xo <- d$lat + rnorm(n, 0, 0.4)
eta <- 0.4 + 0.7 * d$lat
d$cnt <- rpois(n, exp(eta))
d$pos <- exp(eta + rnorm(n, 0, 0.5))
d$hp <- ifelse(d$cnt == 0, 1L, d$cnt)
d$hp[runif(n) < 0.3] <- 0L
d$cat <- as.integer(cut(eta + rlogis(n), c(-Inf, 0, 1.5, Inf), labels = FALSE))
d$x <- d$lat

cat("\n== item 3: se__ on plain fits, x at 0\n")
for (sp in list(list("cnt", poisson()), list("pos", lognormal()),
                list("hp", hurdle_poisson()))) {
  f <- frm(bf(as.formula(paste(sp[[1]], "~ x"))), family = sp[[2]], data = d)
  ce <- q(conditional_effects(f, effects = "x",
                              int_conditions = list(x = c(-1, 0, 1))))
  fv <- fitted(f, newdata = data.frame(x = c(-1, 0, 1)))
  if (!is.null(ce)) {
    print(cbind(ce[[1]][, c("x", "estimate__", "se__", "lower__", "upper__")],
                fitted_Est = fv[, "Estimate"], fitted_EstErr = fv[, "Est.Error"]))
  }
}

cat("\n== item 2: categorical under me(), band on category 3\n")
fc <- frm(bf(cat ~ me(xo, sx)), family = categorical(), data = d)
nd <- data.frame(xo = 1, sx = 1.2)
fv <- fitted(fc, newdata = nd)
print(fv[1, , ])
ce <- q(conditional_effects(fc, effects = "xo", categorical = TRUE,
                            int_conditions = list(xo = 1),
                            conditions = data.frame(sx = 1.2)))
if (!is.null(ce)) print(ce[[1]][, c("xo", "cats__", "estimate__", "se__",
                                    "lower__", "upper__")])

cat("\n== item 5: conditional_effects() on a body whose nlpar reads me()\n")
set.seed(5101)
d5 <- data.frame(xo = rnorm(200), sx = 0.3, z = runif(200))
d5$y <- 2 * exp(0.5 * (d5$xo + rnorm(200, 0, 0.3))) * d5$z + rnorm(200, 0, 0.5)
f5 <- q(frm(bf(y ~ a * exp(b * z), a ~ me(xo, sx), b ~ 1, nl = TRUE),
            family = gaussian(), data = d5))
if (!is.null(f5)) {
  ce5 <- q(conditional_effects(f5, effects = "z"))
  if (!is.null(ce5)) print(head(ce5[[1]][, c("z", "estimate__", "se__",
                                              "lower__", "upper__")], 3))
  cat("fitted on new rows:\n")
  print(q(fitted(f5, newdata = data.frame(xo = 0.5, sx = 0.3, z = 0.5))))
}

cat("\n== item 8: a body's Est.Error\n")
set.seed(8101)
d8 <- data.frame(x = runif(200, 0, 3))
d8$y <- 3 * (1 - exp(-0.8 * d8$x)) + rnorm(200, 0, 0.2)
f8 <- frm(bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE),
          family = gaussian(), data = d8,
          prior = c(prior(normal(3, 2), nlpar = "b1"),
                    prior(normal(1, 1), nlpar = "b2")))
print(q(fitted(f8, newdata = data.frame(x = c(0.5, 2)))))
print(q(fitted(f8)[1:2, ]))
d8$yp <- rpois(200, 3 * (1 - exp(-0.8 * d8$x)))
f8p <- frm(bf(yp ~ log(b1 * (1 - exp(-b2 * x))), b1 + b2 ~ 1, nl = TRUE),
           family = poisson(), data = d8)
print(q(fitted(f8p, newdata = data.frame(x = c(0.5, 2)))))
print(q(predict(f8, newdata = data.frame(x = c(0.5, 2)))))
