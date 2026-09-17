## Reviewer recheck, USER decision: frm_simulate(newparams =) takes brms's
## names only. Does any bare or old name still map silently?
##   Rscript dev/brmsnames-rev2-simulate.R
## Data seed 81.
.libPaths(c("C:/Users/adf44/source/r/brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
try1 <- function(e) {
  tryCatch(e, error = function(err) paste("ERROR:", conditionMessage(err)))
}
q(library(frmtmb))
set.seed(81)
n <- 400
d <- data.frame(x = rnorm(n), b_x = rnorm(n), g = factor(rep(1:10, 40)),
                y = 0)
run <- function(label, f, np, fam = gaussian()) {
  s <- try1(q(frm_simulate(f, data = d, family = fam, newparams = np,
                           seed = 1)))
  out <- if (is.data.frame(s)) {
    sprintf("ran; sd(y) %.3f, cor(y, x) %.3f, cor(y, b_x) %.3f",
            stats::sd(s[[1]]), stats::cor(s[[1]], d$x),
            stats::cor(s[[1]], d$b_x))
  } else substr(gsub("\n", " ", s), 1, 170)
  cat(sprintf("%-44s %s\n", label, out))
}
run("covariate b_x only, b_x = 3", bf(y ~ b_x),
    list(b_Intercept = 0, b_x = 3, sigma = 1))
run("covariate b_x only, b_b_x = 3", bf(y ~ b_x),
    list(b_Intercept = 0, b_b_x = 3, sigma = 1))
run("x and b_x, b_x = 3 b_b_x = 0", bf(y ~ x + b_x),
    list(b_Intercept = 0, b_x = 3, b_b_x = 0, sigma = 1))
run("old name sigma_Intercept = 0", bf(y ~ x),
    list(b_Intercept = 0, b_x = 1, sigma_Intercept = 0))
run("sigma = 2 (natural)", bf(y ~ x),
    list(b_Intercept = 0, b_x = 0, sigma = 2))
run("sigma ~ 1 written, b_sigma_Intercept = log 2", bf(y ~ x, sigma ~ 1),
    list(b_Intercept = 0, b_x = 0, b_sigma_Intercept = log(2)))
run("sigma ~ 1 written, sigma = 2", bf(y ~ x, sigma ~ 1),
    list(b_Intercept = 0, b_x = 0, sigma = 2))
run("bare x = 1", bf(y ~ x), list(b_Intercept = 0, x = 1, sigma = 1))
run("sd_g__(Intercept) old", bf(y ~ x + (1 | g)),
    list(b_Intercept = 0, b_x = 0, sigma = 1, `sd_g__(Intercept)` = 2))
run("sd_g__Intercept = 2", bf(y ~ x + (1 | g)),
    list(b_Intercept = 0, b_x = 0, sigma = 1, sd_g__Intercept = 2))
d$yc <- 0L
run("negbinomial shape_Intercept old", bf(yc ~ x),
    list(b_Intercept = 1, b_x = 0, shape_Intercept = 1), negbinomial())
run("negbinomial shape = 1", bf(yc ~ x),
    list(b_Intercept = 1, b_x = 0, shape = 1), negbinomial())
run("nonlinear a_Intercept old", bf(y ~ a * x, a ~ 1, nl = TRUE),
    list(a_Intercept = 2, sigma = 1))
run("nonlinear b_a_Intercept", bf(y ~ a * x, a ~ 1, nl = TRUE),
    list(b_a_Intercept = 2, sigma = 1))
