# Lane priortab: frm_sample() side of the changes, before and after.
# PRIORTAB_LIB unset: rellib-r20; set: the lane library.
# Rscript dev/priortab-probe2.R > dev/priortab-out/probe2-<arm>.txt
lib <- Sys.getenv("PRIORTAB_LIB")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", dirname(find.package("frmtmb")), "\n")
options(width = 150)
show <- function(label, expr) {
  cat("\n=== ", label, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(force(expr),
         warning = function(w) {
           cat("WARN: ", substr(conditionMessage(w), 1, 160), "\n")
           invokeRestart("muffleWarning")
         },
         message = function(m) invokeRestart("muffleMessage")),
       error = function(e) {
         cat("ERROR: ", conditionMessage(e), "\n")
         NULL
       })
  if (!is.null(r)) print(r)
  invisible(r)
}
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$y <- 1 + d$z * 0.3 + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$y2 <- d$y + 0.5 * d$z
d$Y_a <- d$y
d$Y_b <- d$y2 + rnorm(n)

# 4925 on the sampling route: are the draws inside the support?
show("frm_sample uniform(0, 0.2) on z, no bounds", {
  dr <- frm_sample(frm(y2 ~ z, data = d), chains = 1, iter = 400,
                   warmup = 200, seed = 3, refresh = 0,
                   prior = set_prior("uniform(0, 0.2)", class = "b",
                                     coef = "z"))
  zz <- as.matrix(dr)[, "b_z"]
  c(n = length(zz), min = min(zz), max = max(zz))
})

# 4266: coerced sampling tables
ps <- default_prior(y ~ z + (1 | g), data = d, route = "sample")
f1 <- frm(y ~ z + (1 | g), data = d)
run <- function(pr) {
  dr <- frm_sample(f1, prior = pr, chains = 1, iter = 40, warmup = 20,
                   seed = 3, refresh = 0)
  "sampled"
}
pd <- as.data.frame(ps)
show("4266a coerced sample table + set_prior('', sigma)",
     run(as.brmsprior(pd[pd$class != "sigma", ]) +
           set_prior("", class = "sigma")))
show("4266c coerced sample table edited to '' on sd", {
  pp <- as.data.frame(ps)
  pp$prior[pp$class == "sd" & pp$group == ""] <- ""
  run(as.brmsprior(pp))
})
show("4266d hand-built (flat) on sigma (still refused)",
     run(as.brmsprior(data.frame(prior = "(flat)", class = "sigma"))))
show("4266 guard: coerced FIT table (still refused)",
     run(as.brmsprior(as.data.frame(default_prior(y ~ z + (1 | g),
                                                  data = d)))))

# Y_a and Ya as two slot keys
fm <- bf(mvbind(Y_a, Y_b) ~ z) + set_rescor(FALSE)
show("Y_a and Ya on one slot", {
  f <- frm(fm, data = d, prior = set_prior("normal(0, 1)", class = "b",
                                           resp = "Y_a") +
             set_prior("normal(5, 0.01)", class = "b", resp = "Ya"))
  fixef(f)[, 1]
})
cat("\ndone\n")
