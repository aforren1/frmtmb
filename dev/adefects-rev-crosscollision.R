.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("--- brms:::combine_groups ---\n")
print(brms:::combine_groups)
set.seed(1)
d <- data.frame(t = c(1:6, 7:12, 1:6, 1:6),
                a = c(rep("x_1", 6), rep("x", 6), rep("p", 6),
                      rep("r", 6)),
                b = c(rep("2", 6), rep("1_2", 6), rep("q", 6),
                      rep("s", 6)))
d$x <- rnorm(24); d$y <- rnorm(24)
sd <- tryCatch(brms::make_standata(y ~ x + ar(t, a:b, cov = TRUE), d,
        family = brms::brmsfamily("gaussian")), error = function(e) e)
if (inherits(sd, "condition")) {
  cat("brms a:b on the collision data ERR:", conditionMessage(sd), "\n")
} else {
  cat("brms a:b on the collision data: N =", sd$N,
      " N_tg =", sd$N_tg, " nobs_tg =", paste(sd$nobs_tg, collapse = ","),
      "\n")
}
# control: the same data with no collision
d2 <- d; d2$a[1:6] <- "y"
sd2 <- tryCatch(brms::make_standata(y ~ x + ar(t, a:b, cov = TRUE), d2,
        family = brms::brmsfamily("gaussian")), error = function(e) e)
if (inherits(sd2, "condition")) {
  cat("control ERR:", conditionMessage(sd2), "\n")
} else {
  cat("control (no collision): N_tg =", sd2$N_tg,
      " nobs_tg =", paste(sd2$nobs_tg, collapse = ","), "\n")
}
