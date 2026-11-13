# Fix round of lane priortab: the reviewer's cases on the lane build.
# Rscript dev/priortab-fixprobe.R > dev/priortab-out/fixprobe-<arm>.txt
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
           cat("WARN: ", substr(gsub("\n", " ", conditionMessage(w)), 1, 150),
               "\n")
           invokeRestart("muffleWarning")
         },
         message = function(m) {
           cat("MSG: ", conditionMessage(m))
           invokeRestart("muffleMessage")
         }),
       error = function(e) {
         cat("ERROR: ", substr(conditionMessage(e), 1, 300), "\n")
         NULL
       })
  if (!is.null(r)) print(r)
  invisible(r)
}
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$y <- 1 + d$z * 0.3 + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$yneg <- d$y - 0.8 * d$z
fz <- function(f) c(conv = f$opt$convergence, obj = f$opt$objective,
                    z = unname(fixef(f)["z", "Estimate"]))
for (lbv in list(NULL, 0, 1e-8, 1e-10)) {
  show(paste("beta(2, 2) yneg, lb =", format(lbv)), fz(frm(yneg ~ z,
    data = d, prior = set_prior("beta(2, 2)", class = "b", coef = "z",
                                lb = lbv %||% NA, ub = if (!is.null(lbv)) 1
                                else NA))))
}
show("beta(2, 2) y (interior, no bound)", fz(frm(y ~ z, data = d,
  prior = set_prior("beta(2, 2)", class = "b", coef = "z"))))
show("exponential(1) yneg, no bound", fz(frm(yneg ~ z, data = d,
  prior = set_prior("exponential(1)", class = "b", coef = "z"))))
show("uniform(0.6, 1) on sigma, no bound", {
  a <- frm(y ~ z, data = d, prior = set_prior("uniform(0.6, 1)",
                                              class = "sigma"))
  b <- frm(y ~ z, data = d, prior = set_prior("uniform(0.6, 1)",
                                              class = "sigma", lb = 0.6,
                                              ub = 1))
  list(identical = identical(a$opt$par, b$opt$par), sigma = sigma(a),
       conv = a$opt$convergence)
})
show("uniform(0, 0.5) on an uncentered Intercept", {
  a <- frm(bf(y ~ 0 + Intercept + z), data = d,
           prior = set_prior("uniform(0, 0.5)", class = "b",
                             coef = "Intercept"))
  b <- frm(bf(y ~ 0 + Intercept + z), data = d,
           prior = set_prior("uniform(0, 0.5)", class = "b",
                             coef = "Intercept", lb = 0, ub = 0.5))
  list(identical = identical(a$opt$par, b$opt$par), est = fixef(a)[, 1])
})
show("uniform(0, 0.5) on the centered Intercept", fz(frm(y ~ z, data = d,
  prior = set_prior("uniform(0, 0.5)", class = "Intercept"))))
for (lk in c("softplus", "squareplus")) {
  show(paste("sigma lb under", lk), {
    tab <- as.data.frame(default_prior(bf(y ~ z), data = d,
                                       family = gaussian(link_sigma = lk)))
    tab[tab$class == "sigma", c("lb", "ub")]
  })
}
fml <- frm(y ~ z + (1 | g), data = d)
eff <- function(pl, base = NULL, fit = fml) {
  rp <- frmtmb.sample:::sample_resolve_priors(fit, pl, base = base)
  x <- as.data.frame(rp$effective)
  x[, c("prior", "class", "coef", "group", "lb")]
}
for (cl in c("sd", "sigma")) {
  show(paste("MAJOR 1:", cl, "lb = 0 via fit and direct"), {
    sp <- set_prior("", class = cl, lb = 0)
    f <- frm(y ~ z + (1 | g), data = d, prior = sp)
    list(prior_null = is.null(f$prior),
         via_fit = eff(NULL, base = f$prior, fit = f),
         direct = eff(sp))
  })
}
show("MAJOR 1: table passed back still drops", {
  f <- frm(y ~ z + (1 | g), data = d,
           prior = default_prior(y ~ z + (1 | g), data = d))
  is.null(f$prior)
})
show("H1 coerced fit table, sigma set to brms's default", {
  tab <- as.data.frame(default_prior(y ~ z + (1 | g), data = d))
  tab$prior[tab$class == "sigma"] <- "student_t(3, 0, 2.5)"
  eff(as.brmsprior(tab))
})
show("H3 coerced sample table, untouched", {
  eff(as.brmsprior(as.data.frame(default_prior(y ~ z + (1 | g), data = d,
                                               route = "sample"))))
})
show("H4 coerced sample table, Intercept to (flat)", {
  tab <- as.data.frame(default_prior(y ~ z + (1 | g), data = d,
                                     route = "sample"))
  tab$prior[tab$class == "Intercept"] <- "(flat)"
  eff(as.brmsprior(tab))
})
show("4266 coerced sample table, sd to ''", {
  tab <- as.data.frame(default_prior(y ~ z + (1 | g), data = d,
                                     route = "sample"))
  tab$prior[tab$class == "sd" & tab$group == ""] <- ""
  eff(as.brmsprior(tab))
})
line <- function(pr, f) {
  m <- character()
  withCallingHandlers(frmtmb.sample:::sample_resolve_priors(f, pr),
                      message = function(e) {
                        m <<- c(m, conditionMessage(e))
                        invokeRestart("muffleMessage")
                      })
  grep("^  b ", strsplit(paste(m, collapse = "\n"), "\n")[[1]], value = TRUE)
}
f2 <- frm(y ~ z + x, data = d, dry_run = "objective")
show("slopes: class b on all", line(set_prior("normal(0, 1)", class = "b"), f2))
show("slopes: coef z only", line(set_prior("normal(0, 1)", class = "b",
                                           coef = "z"), f2))
show("slopes: no prior", line(NULL, f2))
d$Y_a <- d$y
d$Y_b <- d$y + rnorm(n)
fm <- bf(mvbind(Y_a, Y_b) ~ z) + set_rescor(FALSE)
show("MINOR 4: Y_a table passed back, effective rows", {
  fmv <- frm(fm, data = d)
  tab <- default_prior(fm, data = d, route = "sample")
  nrow(eff(tab, fit = fmv))
})
cat("\ndone\n")
