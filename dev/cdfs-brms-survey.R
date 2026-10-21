# Which families does brms 2.23.0 censor and truncate, and with which
# Stan calls? stancode() only; no compile. Output: dev/cdfs-out/survey.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(20261020)
n <- 40
x <- rnorm(n)
cc <- rep(c(-1, 0, 1, 2), length.out = n)
pos <- rgamma(n, 3, 1)
cnt <- rpois(n, 3)
cnt0 <- cnt; cnt0[1:8] <- 0
unit <- rbeta(n, 2, 3)
unit0 <- unit; unit0[1:5] <- 0
pos0 <- pos; pos0[1:5] <- 0
ang <- runif(n, -pi, pi)
tr <- rep(10, n)
bin <- rbinom(n, 10, 0.4)
real <- rnorm(n, 1)
cases <- list(
  student = list(y = real, fam = student()),
  skew_normal = list(y = real, fam = skew_normal()),
  exgaussian = list(y = pos, fam = exgaussian()),
  shifted_lognormal = list(y = pos + 0.5, fam = shifted_lognormal()),
  asym_laplace = list(y = real, fam = asym_laplace()),
  zero_inflated_asym_laplace = list(y = c(0, real[-1]),
                                    fam = brmsfamily("zero_inflated_asym_laplace")),
  von_mises = list(y = ang, fam = von_mises()),
  Gamma = list(y = pos, fam = Gamma(link = "log")),
  Beta = list(y = unit, fam = Beta()),
  xbeta = list(y = unit, fam = brmsfamily("xbeta")),
  binomial = list(y = bin, fam = binomial(), trials = TRUE),
  bernoulli = list(y = as.integer(bin > 4), fam = bernoulli()),
  beta_binomial = list(y = bin, fam = beta_binomial(), trials = TRUE),
  negbinomial = list(y = cnt, fam = negbinomial()),
  negbinomial2 = list(y = cnt, fam = brmsfamily("negbinomial2")),
  geometric = list(y = cnt, fam = geometric()),
  com_poisson = list(y = cnt, fam = brmsfamily("com_poisson")),
  poisson = list(y = cnt, fam = poisson()),
  inverse.gaussian = list(y = pos, fam = inverse.gaussian()),
  zero_inflated_poisson = list(y = cnt0, fam = zero_inflated_poisson()),
  zero_inflated_negbinomial = list(y = cnt0,
                                   fam = zero_inflated_negbinomial()),
  zero_inflated_binomial = list(y = bin, fam = zero_inflated_binomial(),
                                trials = TRUE),
  zero_inflated_beta_binomial = list(y = bin,
                                     fam = zero_inflated_beta_binomial(),
                                     trials = TRUE),
  zero_inflated_beta = list(y = unit0, fam = zero_inflated_beta()),
  zero_one_inflated_beta = list(y = unit0, fam = zero_one_inflated_beta()),
  hurdle_poisson = list(y = cnt0, fam = hurdle_poisson()),
  hurdle_negbinomial = list(y = cnt0, fam = hurdle_negbinomial()),
  hurdle_gamma = list(y = pos0, fam = hurdle_gamma()),
  hurdle_lognormal = list(y = pos0, fam = hurdle_lognormal())
)
grab <- function(code, fam) {
  m <- strsplit(code, "\n")[[1]]
  st <- grep("^model \\{", m)
  model <- m[st:length(m)]
  hits <- grep("_lcdf|_lccdf|_lcdf|log_diff_exp|log_sum_exp", model,
               value = TRUE)
  hits <- unique(trimws(hits))
  fun <- m[seq_len(st)]
  defs <- grep("real [a-z_0-9]+_(lcdf|lccdf)\\(", fun, value = TRUE)
  list(model = hits, defs = unique(trimws(defs)))
}
run1 <- function(nm, cs, form_rhs) {
  d <- data.frame(y = cs$y, x = x, cc = cc, tr = tr)
  d$y2 <- d$y + ifelse(cs$fam$family %in% c("von_mises"), 0.3,
                       ifelse(is.integer(d$y) || all(d$y == round(d$y)),
                              2, 0.5))
  if (cs$fam$family %in% c("beta", "zero_inflated_beta",
                           "zero_one_inflated_beta", "xbeta")) {
    d$y2 <- pmin(d$y + 0.1, 0.999)
  }
  lhs <- if (isTRUE(cs$trials)) "y | trials(tr) +" else "y |"
  f <- as.formula(paste(lhs, form_rhs, "~ x"))
  r <- tryCatch({
    code <- stancode(f, data = d, family = cs$fam)
    g <- grab(code, cs$fam)
    paste0("OK\n    model: ", paste(g$model, collapse = "\n           "),
           "\n    defs : ", paste(g$defs, collapse = "\n           "))
  }, error = function(e) paste0("ERROR: ", conditionMessage(e)))
  cat("==", nm, "|", form_rhs, "\n   ", r, "\n")
}
for (nm in names(cases)) {
  cs <- cases[[nm]]
  run1(nm, cs, "cens(cc, y2)")
  lb <- if (cs$fam$family %in% c("von_mises")) -2 else
    if (cs$fam$family %in% c("student", "skew_normal", "asym_laplace",
                             "zero_inflated_asym_laplace")) -1 else
      if (cs$fam$family %in% c("beta", "zero_inflated_beta",
                               "zero_one_inflated_beta", "xbeta")) 0.05 else 1
  ub <- if (cs$fam$family %in% c("von_mises")) 2 else
    if (cs$fam$family %in% c("beta", "zero_inflated_beta",
                             "zero_one_inflated_beta", "xbeta")) 0.95 else 30
  run1(nm, cs, sprintf("trunc(lb = %s, ub = %s)", lb, ub))
  run1(nm, cs, sprintf("trunc(lb = %s)", lb))
}
