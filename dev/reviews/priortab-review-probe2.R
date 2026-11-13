# Reviewer of lane priortab: support box, drop_unbounding_slots(), 4266.
# PRIORTAB_LIB unset = lane library; PRIORTAB_LIB=none = rellib-r20.
# Rscript dev/reviews/priortab-review-probe2.R > dev/reviews/priortab-review-probe2-<arm>.txt
lib <- Sys.getenv("PRIORTAB_LIB", "C:/Users/adf44/source/r/wt-priortab-lib")
if (identical(lib, "none")) lib <- ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", dirname(find.package("frmtmb")), "\n")
cat("frmtmb.sample", dirname(find.package("frmtmb.sample")), "\n")
options(width = 150, digits = 10)
show <- function(label, expr) {
  cat("\n=== ", label, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(force(expr),
         warning = function(w) {
           cat("WARN: ", substr(gsub("\n", " ", conditionMessage(w)), 1, 140),
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
d$y <- 1 + d$z * 0.3 + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$yneg <- d$y - 0.8 * d$z
d$o <- cut(d$y, 4, labels = FALSE)
d$cat <- factor(sample(c("a", "b", "c"), n, TRUE))
d$ys <- 1 + rnorm(n, sd = exp(0.2 + 0.3 * d$x))
summ <- function(f, coef = NULL) {
  out <- list(conv = f$opt$convergence, obj = f$opt$objective,
              msg = f$opt$message)
  if (!is.null(coef)) out$par <- f$opt$par[coef]
  out$prior_null <- is.null(f$prior)
  out
}

## --- support box on half-bounded and other placements
show("B1 lognormal(0,1) on b z, ML slope -0.5", {
  f <- frm(yneg ~ z, data = d, prior = set_prior("lognormal(0, 1)",
                                                class = "b", coef = "z"))
  c(summ(f), z = unname(fixef(f)["z", "Estimate"]))
})
show("B1e lognormal(0,1) on b z, lb = 0 given", {
  f <- frm(yneg ~ z, data = d, prior = set_prior("lognormal(0, 1)",
                                                class = "b", coef = "z", lb = 0))
  c(summ(f), z = unname(fixef(f)["z", "Estimate"]))
})
show("B2 gamma(0.5, 1) on b z, ML slope -0.5", {
  f <- frm(yneg ~ z, data = d, prior = set_prior("gamma(0.5, 1)",
                                                class = "b", coef = "z"))
  c(summ(f), z = unname(fixef(f)["z", "Estimate"]))
})
show("B2e gamma(0.5, 1) on b z, lb = 0 given", {
  f <- frm(yneg ~ z, data = d, prior = set_prior("gamma(0.5, 1)",
                                                class = "b", coef = "z", lb = 0))
  c(summ(f), z = unname(fixef(f)["z", "Estimate"]))
})
show("B3 exponential(1) on b z, ML slope -0.5", {
  f <- frm(yneg ~ z, data = d, prior = set_prior("exponential(1)",
                                                class = "b", coef = "z"))
  c(summ(f), z = unname(fixef(f)["z", "Estimate"]))
})
show("B4 beta(2, 2) on b z", {
  f <- frm(y ~ z, data = d, prior = set_prior("beta(2, 2)",
                                             class = "b", coef = "z"))
  c(summ(f), z = unname(fixef(f)["z", "Estimate"]))
})
show("B5 uniform(0, 0.2) class-wide b, two slopes", {
  f <- frm(y ~ z + x, data = d, prior = set_prior("uniform(-0.1, 0.2)",
                                                 class = "b"))
  c(summ(f), b = fixef(f)[c("z", "x"), "Estimate"])
})
show("B6 uniform(-1, 0.1) on sigma's slope x (sigma ~ x)", {
  f <- frm(bf(ys ~ 1, sigma ~ x), data = d,
           prior = set_prior("uniform(-1, 0.1)", class = "b", coef = "x",
                             dpar = "sigma"))
  summ(f)
})
show("B7 uniform(0, 0.5) on sigma's Intercept (sigma ~ x)", {
  f <- frm(bf(ys ~ 1, sigma ~ x), data = d,
           prior = set_prior("uniform(0, 0.1)", class = "Intercept",
                             dpar = "sigma"))
  summ(f)
})
show("B8 uniform(0, 0.2) on nlpar b", {
  f <- frm(bf(y ~ a + b * z, a ~ 1, b ~ 1, nl = TRUE), data = d,
           prior = set_prior("uniform(0, 0.2)", nlpar = "b"))
  summ(f)
})
show("B8e uniform(0, 0.2) on nlpar b with lb/ub", {
  f <- frm(bf(y ~ a + b * z, a ~ 1, b ~ 1, nl = TRUE), data = d,
           prior = set_prior("uniform(0, 0.2)", nlpar = "b", lb = 0, ub = 0.2))
  summ(f)
})
show("B9 uniform(0.6, 1) on class sigma (natural, not boxed)", {
  f <- frm(y ~ z, data = d, prior = set_prior("uniform(0.6, 1)",
                                             class = "sigma"))
  summ(f)
})
show("B10 uniform(0, 0.5) on centered mu Intercept", {
  f <- frm(y ~ z, data = d, prior = set_prior("uniform(0, 0.5)",
                                             class = "Intercept"))
  summ(f)
})
show("B11 acat cs(z) uniform(0, 0.1) on b z", {
  f <- frm(o ~ cs(z), data = d, family = acat(),
           prior = set_prior("uniform(0, 0.1)", class = "b", coef = "z"))
  summ(f)
})
show("B12 categorical uniform(0, 0.1) on b z dpar mub", {
  f <- frm(cat ~ z, data = d, family = categorical(),
           prior = set_prior("uniform(0, 0.1)", class = "b", coef = "z",
                             dpar = "mub"))
  summ(f)
})
show("B13 uniform(0, 0.2) on b z + lb = -1 (user lb outside support)", {
  f <- frm(y ~ z, data = d, prior = set_prior("uniform(0, 0.2)", class = "b",
                                             coef = "z", lb = -1))
  summ(f)
})
show("B14 class b normal + coef z uniform(0, 0.1)", {
  f <- frm(y ~ z + x, data = d, prior = set_prior("normal(0, 1)", class = "b") +
             set_prior("uniform(0, 0.1)", class = "b", coef = "z"))
  c(summ(f), b = fixef(f)[c("z", "x"), "Estimate"])
})
show("B15 mo(): uniform(0, 0.1) on b moom", {
  d$om <- factor(sample(1:4, n, TRUE), ordered = TRUE)
  f <- frm(y ~ mo(om), data = d, prior = set_prior("uniform(0, 0.1)",
                                                  class = "b", coef = "moom"))
  summ(f)
})

## --- drop_unbounding_slots and frm_sample(fit)
fml <- frm(y ~ z + (1 | g), data = d)
eff_rows <- function(rp) {
  if (is.null(rp$effective)) return("effective NULL")
  x <- as.data.frame(rp$effective)
  x[x$class %in% c("sd", "sigma"), intersect(c("prior", "class", "group",
                                                "lb", "ub"), names(x))]
}
for (nm in c("sd", "sigma")) {
  sp <- set_prior("", class = nm, lb = 0)
  show(paste("D1", nm, "lb=0: frm fit$prior"), {
    f <- frm(y ~ z + (1 | g), data = d, prior = sp)
    list(prior_null = is.null(f$prior), obj = f$opt$objective)
  })
  show(paste("D2", nm, "lb=0: frm_sample(fit) effective"), {
    f <- frm(y ~ z + (1 | g), data = d, prior = sp)
    eff_rows(frmtmb.sample:::sample_resolve_priors(f, NULL, base = f$prior,
                                                   defaults = TRUE))
  })
  show(paste("D3", nm, "lb=0: frm_sample(ml fit, prior = sp) effective"), {
    eff_rows(frmtmb.sample:::sample_resolve_priors(
      fml, frmtmb:::as_priorlist(sp), base = NULL, defaults = TRUE))
  })
}
show("D4 frm(prior = table + set_prior same slot)", {
  tab <- default_prior(y ~ z + x, data = d)
  sp <- set_prior("normal(0, 1)", class = "b", coef = "z")
  pr <- tryCatch(tab + sp, error = function(e) c(tab, sp))
  f <- frm(y ~ z + x, data = d, prior = pr)
  as.data.frame(f$prior)
})
show("D5 frm(prior = sigma lb=0 twice)", {
  f <- frm(y ~ z, data = d, prior = c(set_prior("", class = "sigma", lb = 0),
                                      set_prior("", class = "sigma", lb = 0)))
  is.null(f$prior)
})
show("D6 frm(prior = table) with an sd row edited to lb = 0.05", {
  tab <- default_prior(y ~ z + (1 | g), data = d)
  tab$lb[tab$class == "sd" & tab$group == ""] <- 0.05
  f <- frm(y ~ z + (1 | g), data = d, prior = tab)
  list(prior_null = is.null(f$prior), obj = f$opt$objective,
       sd = exp(f$opt$par[grep("theta", names(f$opt$par))]))
})

## --- 4266: coerced fit-route table with one row equal to a brms default
show("H1 coerced fit table, sigma row set to brms's default", {
  tab <- as.data.frame(default_prior(y ~ z + (1 | g), data = d))
  tab$prior[tab$class == "sigma"] <- "student_t(3, 0, 2.5)"
  pl <- frmtmb:::as_priorlist(tab)
  rp <- frmtmb.sample:::sample_resolve_priors(fml, pl, base = NULL,
                                              defaults = TRUE)
  x <- as.data.frame(rp$effective)
  x[, intersect(c("prior", "class", "coef", "group", "lb"), names(x))]
})
show("H2 coerced fit table, untouched", {
  tab <- as.data.frame(default_prior(y ~ z + (1 | g), data = d))
  pl <- frmtmb:::as_priorlist(tab)
  rp <- frmtmb.sample:::sample_resolve_priors(fml, pl, base = NULL,
                                              defaults = TRUE)
  "accepted"
})
show("H3 coerced sample table, untouched", {
  tab <- as.data.frame(default_prior(y ~ z + (1 | g), data = d,
                                     route = "sample"))
  pl <- frmtmb:::as_priorlist(tab)
  rp <- frmtmb.sample:::sample_resolve_priors(fml, pl, base = NULL,
                                              defaults = TRUE)
  x <- as.data.frame(rp$effective)
  x[, intersect(c("prior", "class", "coef", "group", "lb"), names(x))]
})
cat("\ndone\n")
