# Reviewer of lane priortab: Stan-name resp in a sampling table passed
# back to frm_sample() (resolution only, no Stan).
# PRIORTAB_LIB unset = lane library; PRIORTAB_LIB=none = rellib-r20.
lib <- Sys.getenv("PRIORTAB_LIB", "C:/Users/adf44/source/r/wt-priortab-lib")
if (identical(lib, "none")) lib <- ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", dirname(find.package("frmtmb")), "\n")
options(width = 150, digits = 10)
show <- function(label, expr) {
  cat("\n=== ", label, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(force(expr),
         warning = function(w) {
           cat("WARN: ", substr(gsub("\n", " ", conditionMessage(w)), 1, 200),
               "\n")
           invokeRestart("muffleWarning")
         },
         message = function(m) invokeRestart("muffleMessage")),
       error = function(e) {
         cat("ERROR: ", substr(conditionMessage(e), 1, 400), "\n")
         NULL
       })
  if (!is.null(r)) print(r)
  invisible(r)
}
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$Y_a <- 1 + d$z * 0.3 + rnorm(n, sd = 0.3)
d$Y_b <- d$Y_a + rnorm(n)
f <- bf(mvbind(Y_a, Y_b) ~ z) + set_rescor(FALSE)
fit <- frm(f, data = d, dry_run = "objective")
res <- function(pr) {
  pl <- frmtmb:::as_priorlist(pr)
  frmtmb:::check_prior_slots(pl)
  rp <- frmtmb.sample:::sample_resolve_priors(fit, pl)
  x <- as.data.frame(rp$effective)
  x <- x[x$class %in% c("sigma", "Intercept", "b"),
         intersect(c("prior", "class", "coef", "resp", "lb"), names(x))]
  list(eff = x, n_entries = length(rp$ri$entries),
       entries = vapply(rp$ri$entries, function(e)
         paste(e$comp, paste(e$idx, collapse = ","), e$dist$kind), ""))
}
tab <- default_prior(f, data = d, route = "sample")
show("T0 sampling table resp column", unique(as.data.frame(tab)$resp))
show("T1 sampling table passed back unmodified", res(tab))
show("T2 sampling table, sigma Ya/Y_a row edited to normal(0, 1)", {
  t2 <- tab
  t2$prior[t2$class == "sigma" & t2$resp %in% c("Ya", "Y_a")] <-
    "normal(0, 1)"
  res(t2)
})
show("T3 set_prior sigma resp = 'Ya'", res(set_prior("normal(0, 1)",
                                                    class = "sigma",
                                                    resp = "Ya")))
show("T4 set_prior sigma resp = 'Y_a'", res(set_prior("normal(0, 1)",
                                                     class = "sigma",
                                                     resp = "Y_a")))
show("T5 defaults only", res(NULL))
cat("\ndone\n")
