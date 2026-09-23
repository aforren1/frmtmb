# Punch round 1: the measurements behind each fix, on the lane build.
#   Rscript dev/correct-punch1.R
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-prelude.R")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-data.R")
show <- function(label, expr) {
  r <- tryCatch({
    v <- withCallingHandlers(force(expr), warning = function(w) {
      cat(sprintf("%-42s WARNING: %s\n", label,
                  substr(conditionMessage(w), 1, 90)))
      invokeRestart("muffleWarning")
    })
    if (inherits(v, "ggplot")) invisible(ggplot2::ggplot_build(v))
    "OK"
  }, error = function(e) {
    paste0("ERROR [", class(e)[1], "]: ",
           substr(gsub("\n", " ", conditionMessage(e)), 1, 150))
  })
  cat(sprintf("%-42s %s\n", label, r))
}

cat("== minor 6: the refusal quotes the caller's type ==\n")
d <- correct_data_ord()
fo <- frm(bf(ord ~ x) + cumulative(), data = d)
for (ty in c("response", "ordinary", "pearson")) {
  cat(sprintf("%-12s %s\n", ty,
              tryCatch(residuals(fo, type = ty),
                       error = function(e) substr(conditionMessage(e), 1, 70))))
}

cat("\n== minor 2: resp on a one-response fit, and minor 8: the warning ==\n")
dg <- correct_data_gauss()
fg <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dg)
set.seed(3)
p1 <- suppressMessages(pp_check(fg, type = "dens_overlay", ndraws = 5))
set.seed(3)
p2 <- suppressMessages(pp_check(fg, type = "dens_overlay", ndraws = 5,
                                resp = "w"))
cat("resp = \"w\" plots the same data:", identical(p1$data, p2$data), "\n")
show("dens_overlay, group = g", suppressMessages(
  pp_check(fg, type = "dens_overlay", group = "g", ndraws = 5)))
show("dens_overlay, group = nosuch", suppressMessages(
  pp_check(fg, type = "dens_overlay", group = "nosuch", ndraws = 5)))
show("stat_grouped, group = nosuch", suppressMessages(
  pp_check(fg, type = "stat_grouped", group = "nosuch", ndraws = 5)))

cat("\n== minor 3: error_binned on ordinal draws and on the fit ==\n")
show("fit, error_binned", pp_check(fo, type = "error_binned", ndraws = 5))
ds <- suppressWarnings(suppressMessages(
  frm_sample(fo, chains = 1, iter = 200, refresh = 0, seed = 3)))
show("draws, error_binned", suppressMessages(
  pp_check(ds, type = "error_binned", ndraws = 5)))

cat("\n== minor 4: variables() order on a mi() interaction ==\n")
dmi <- correct_data_mi()
f2 <- frm(bf(y ~ mi(xm) * z) + bf(xm | mi() ~ z) + set_rescor(FALSE) +
            gaussian(), data = dmi)
cat("variables:", paste(variables(f2), collapse = " "), "\n")
cat("fixef rows:", paste(rownames(fixef(f2)), collapse = " "), "\n")

cat("\n== minor 5: class sd on a location-scale smooth ==\n")
set.seed(44)
n <- 400
dd <- data.frame(x = runif(n), z = runif(n))
dd$y <- rnorm(n, sin(3 * dd$x), exp(0.8 * cos(5 * dd$z) - 0.5))
fs <- frm(bf(y ~ s(x), sigma ~ s(z)) + gaussian(), data = dd)
reach <- function(pl) {
  e <- frmtmb:::resolve_prior_input(fs, pl)$entries
  as.numeric(sort(unlist(lapply(e, function(z) {
    if (identical(z$comp, "theta")) z$idx
  }))))
}
blk <- vapply(fs$frame$re_blocks, function(b) {
  paste0(b$dpar, "/", paste(b$theta_idx, collapse = ","))
}, "")
cat("smooth blocks:", blk, "\n")
cat("class sd reaches theta:", reach(set_prior("normal(0, 5)",
                                               class = "sd")), "\n")
cat("class sd dpar sigma reaches theta:",
    reach(set_prior("normal(0, 5)", class = "sd", dpar = "sigma")), "\n")
tab <- default_prior(bf(y ~ s(x), sigma ~ s(z)) + gaussian(), data = dd)
print(as.data.frame(tab)[tab$class == "sd", c("prior", "class", "coef",
                                              "group", "dpar")])
ts <- default_prior(bf(y ~ s(x), sigma ~ s(z)) + gaussian(), data = dd,
                    route = "sample")
print(as.data.frame(ts)[ts$class == "sd", c("prior", "class", "coef",
                                            "group", "dpar")])

cat("\n== minor 2, the draws half: resp on a one-response draws object ==\n")
dsg <- suppressWarnings(suppressMessages(
  frm_sample(fg, chains = 1, iter = 200, refresh = 0, seed = 5)))
set.seed(3)
q1 <- suppressMessages(pp_check(dsg, type = "dens_overlay", ndraws = 5))
set.seed(3)
q2 <- suppressMessages(pp_check(dsg, type = "dens_overlay", ndraws = 5,
                                resp = "w"))
cat("draws, resp = \"w\" plots the same data:",
    identical(q1$data, q2$data), "\n")

cat("\n== minor 1: what a multinomial fit answers, against brms ==\n")
dm <- correct_data_multinom()
fm <- frm(bf(Y | trials(n) ~ x) + multinomial(K = 3), data = dm)
show("multinomial, error_binned", pp_check(fm, type = "error_binned",
                                           ndraws = 5))
show("multinomial, error_hist", pp_check(fm, type = "error_hist",
                                         ndraws = 5))
