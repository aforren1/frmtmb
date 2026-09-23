# Lane wt-correct: what brms 2.23.0 does for each of the seven items,
# measured on the data the frmtmb probe (dev/correct-probe.R) uses.
#
# Usage: Rscript dev/correct-brms-probe.R <part>
#   part = resid | ppcheck | priors | mi | theta | all
# Output goes to stdout; the caller tees it into dev/correct-log/.
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-prelude.R")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), " bayesplot",
    format(packageVersion("bayesplot")), "\n")
part <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(part)) part <- "all"
FITS <- "C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-brmsfits"
dir.create(FITS, showWarnings = FALSE)
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-data.R")

bfit <- function(name, ...) {
  brm(..., chains = 1, iter = 600, warmup = 300, refresh = 0,
      seed = 1, backend = "rstan", silent = 2,
      file = file.path(FITS, name))
}
try_msg <- function(expr) {
  r <- tryCatch({
    v <- force(expr)
    if (inherits(v, "ggplot")) invisible(ggplot2::ggplot_build(v))
    "OK"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  gsub("\n", " ", r)
}

if (part %in% c("resid", "all")) {
  cat("\n== item 1: residuals() by family ==\n")
  d <- correct_data_ord()
  fits <- list(
    cumulative = bfit("ord-cum", ord ~ x, data = d, family = cumulative()),
    sratio = bfit("ord-sratio", ord ~ x, data = d, family = sratio()),
    categorical = bfit("cat", cat ~ x, data = d, family = categorical()),
    gaussian = bfit("gauss-ord", y ~ x, data = d, family = gaussian())
  )
  dm <- correct_data_multinom()
  fits$multinomial <- bfit("multinom", bf(Y | trials(n) ~ x), data = dm,
                           family = multinomial())
  for (nm in names(fits)) {
    f <- fits[[nm]]
    for (ty in c("ordinary", "pearson")) {
      cat(sprintf("%-12s residuals(type = %-9s) %s\n", nm,
                  paste0("\"", ty, "\""),
                  try_msg(suppressWarnings(residuals(f, type = ty)))))
    }
    for (m in c("posterior_predict", "posterior_epred")) {
      cat(sprintf("%-12s predictive_error(method = %s) %s\n", nm, m,
                  try_msg(predictive_error(f, method = m))))
    }
    cat(sprintf("%-12s pp_check(type = \"error_binned\") %s\n", nm,
                try_msg(suppressMessages(pp_check(f, type = "error_binned",
                                                  ndraws = 5)))))
    cat(sprintf("%-12s pp_check(type = \"error_hist\") %s\n", nm,
                try_msg(suppressMessages(pp_check(f, type = "error_hist",
                                                  ndraws = 5)))))
  }
}

if (part %in% c("ppcheck", "all")) {
  cat("\n== item 2: pp_check(), every type bayesplot exposes ==\n")
  d <- correct_data_gauss()
  f <- bfit("gauss-g", y ~ x + (1 | g), data = d, family = gaussian())
  ty <- sub("^ppc_", "", as.character(bayesplot::available_ppc("")))
  ns <- asNamespace("bayesplot")
  for (t in ty) {
    fa <- names(formals(get(paste0("ppc_", t), ns)))
    args <- list(f, type = t, ndraws = 20)
    if ("group" %in% fa) args$group <- "g"
    if ("x" %in% fa) args$x <- "x"
    cat(sprintf("%-28s %s\n", t,
                try_msg(suppressWarnings(suppressMessages(
                  do.call(pp_check, args))))))
  }
  cat("\n-- refusals --\n")
  cat("violin:", try_msg(pp_check(f, type = "violin")), "\n")
  cat("stat_grouped, no group:",
      try_msg(pp_check(f, type = "stat_grouped")), "\n")
  cat("stat_grouped, group = nosuch:",
      try_msg(pp_check(f, type = "stat_grouped", group = "nosuch")), "\n")
  cat("stat_grouped, group = z (in data, not in model):",
      try_msg(pp_check(f, type = "stat_grouped", group = "z")), "\n")
  cat("intervals, x = z (in data, not in model):",
      try_msg(pp_check(f, type = "intervals", x = "z")), "\n")
  cat("stat_grouped, group = c('g','g'):",
      try_msg(pp_check(f, type = "stat_grouped", group = c("g", "g"))), "\n")
  cat("dens_overlay, group = g (ignored?):",
      try_msg(suppressMessages(pp_check(f, type = "dens_overlay",
                                        group = "g", ndraws = 5))), "\n")
  nd <- d[1:12, ]
  cat("stat_grouped, newdata 12 rows:",
      try_msg(suppressMessages(pp_check(f, type = "stat_grouped",
                                        group = "g", newdata = nd))), "\n")
  p <- suppressMessages(pp_check(f, type = "stat_grouped", group = "g",
                                 ndraws = 5))
  cat("stat_grouped group column in plot data equals d$g:",
      identical(as.character(p$data$group %||% NA), NA_character_), "\n")
  print(utils::head(ggplot2::ggplot_build(p)$data[[1]], 3))
}

if (part %in% c("ppcheck2", "all")) {
  cat("\n== item 2: pp_check() on discrete responses ==\n")
  d <- correct_data_counts()
  fits <- list(
    poisson = bfit("pois-g", cnt ~ x + (1 | g), data = d, family = poisson()),
    bernoulli = bfit("bern-g", bin ~ x + (1 | g), data = d,
                     family = bernoulli()))
  ty <- sub("^ppc_", "", as.character(bayesplot::available_ppc("")))
  ns <- asNamespace("bayesplot")
  for (nm in names(fits)) {
    for (t in ty) {
      fa <- names(formals(get(paste0("ppc_", t), ns)))
      args <- list(fits[[nm]], type = t, ndraws = 20)
      if ("group" %in% fa) args$group <- "g"
      if ("x" %in% fa) args$x <- "x"
      cat(sprintf("%-10s %-28s %s\n", nm, t,
                  try_msg(suppressWarnings(suppressMessages(
                    do.call(pp_check, args))))))
    }
  }
  cat("\n-- prefix ppd --\n")
  for (t in c("dens_overlay", "stat_grouped", "intervals", "violin_grouped")) {
    fa <- names(formals(get(paste0("ppd_", t), ns)))
    args <- list(fits$poisson, type = t, ndraws = 20, prefix = "ppd")
    if ("group" %in% fa) args$group <- "g"
    if ("x" %in% fa) args$x <- "x"
    cat(sprintf("ppd %-24s %s\n", t, try_msg(suppressWarnings(
      suppressMessages(do.call(pp_check, args))))))
  }
  cat("ppd violin:", try_msg(pp_check(fits$poisson, type = "violin",
                                      prefix = "ppd")), "\n")
}

if (part %in% c("priors2", "all")) {
  cat("\n== item 3: sd prior scope, nonlinear and cor ==\n")
  set.seed(111)
  n <- 120
  dnl <- data.frame(x = runif(n, 0, 3), g = factor(rep(1:8, length.out = n)))
  dnl$y <- (2 + rnorm(8, 0, 0.3)[dnl$g]) * exp(-0.5 * dnl$x) +
    rnorm(n, 0, 0.1)
  fnl <- bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE)
  pnl <- set_prior("normal(0, 1)", nlpar = "a") +
    set_prior("normal(0, 1)", nlpar = "b")
  for (p in list(set_prior("normal(0, 5)", class = "sd"),
                 set_prior("normal(0, 5)", class = "sd", group = "g"),
                 set_prior("normal(0, 5)", class = "sd", nlpar = "a"))) {
    cat(sprintf("nl class sd group %-2s nlpar %-2s: ", p$group, p$nlpar))
    r <- tryCatch({
      sc <- stancode(fnl, data = dnl, prior = pnl + p)
      l <- grep("lprior \\+= .*sd_", strsplit(sc, "\n")[[1]], value = TRUE)
      paste(trimws(l), collapse = " | ")
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(gsub("\n", " ", r), "\n")
  }
  print(as.data.frame(default_prior(fnl, dnl))[, c("prior", "class", "coef",
                                                   "group", "dpar",
                                                   "nlpar")])
  d <- correct_data_beta()
  d$x2 <- rnorm(nrow(d))
  fb <- bf(yb ~ x + (1 + x | g), phi ~ (1 + x2 | g))
  for (p in list(set_prior("lkj(2)", class = "cor"),
                 set_prior("lkj(2)", class = "cor", group = "g"),
                 set_prior("lkj(2)", class = "cor", dpar = "phi"))) {
    cat(sprintf("Beta class cor group %-2s dpar %-3s: ", p$group, p$dpar))
    r <- tryCatch({
      sc <- stancode(fb, data = d, family = Beta(), prior = p)
      l <- grep("lkj", strsplit(sc, "\n")[[1]], value = TRUE)
      paste(trimws(l), collapse = " | ")
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(gsub("\n", " ", r), "\n")
  }
  print(as.data.frame(default_prior(fb, d, family = Beta()))[
    , c("prior", "class", "coef", "group", "dpar")])
  cat("\n-- (1 | q | g) across mu and sigma --\n")
  dg <- correct_data_gauss()
  fq <- bf(y ~ x + (1 | q | g), sigma ~ (1 | q | g))
  for (p in list(set_prior("normal(0, 5)", class = "sd"),
                 set_prior("normal(0, 5)", class = "sd", group = "g"),
                 set_prior("normal(0, 5)", class = "sd", dpar = "sigma"))) {
    cat(sprintf("q class sd group %-2s dpar %-5s: ", p$group, p$dpar))
    r <- tryCatch({
      sc <- stancode(fq, data = dg, prior = p)
      l <- grep("lprior \\+= .*sd_", strsplit(sc, "\n")[[1]], value = TRUE)
      paste(trimws(l), collapse = " | ")
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(gsub("\n", " ", r), "\n")
  }
  print(as.data.frame(default_prior(fq, dg))[, c("prior", "class", "coef",
                                                 "group", "dpar")])
}

if (part %in% c("priors", "all")) {
  cat("\n== item 3: sd prior scope ==\n")
  d <- correct_data_beta()
  pr <- set_prior("normal(0, 5)", class = "sd", group = "g")
  sc <- stancode(bf(yb ~ x + (1 | g), phi ~ (1 | g)), data = d,
                 family = Beta(), prior = pr)
  cat(grep("sd_|lprior", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
  vp <- validate_prior(pr, bf(yb ~ x + (1 | g), phi ~ (1 | g)), data = d,
                       family = Beta())
  print(vp[vp$class == "sd", c("prior", "class", "coef", "group", "dpar",
                               "source")])
  pr2 <- set_prior("normal(0, 5)", class = "sd")
  sc2 <- stancode(bf(yb ~ x + (1 | g), phi ~ (1 | g)), data = d,
                  family = Beta(), prior = pr2)
  cat("-- class sd, no group, no dpar --\n")
  cat(grep("target|lprior", strsplit(sc2, "\n")[[1]], value = TRUE),
      sep = "\n")
  pr3 <- set_prior("normal(0, 5)", class = "sd", coef = "Intercept",
                   group = "g")
  sc3 <- stancode(bf(yb ~ x + (1 | g), phi ~ (1 | g)), data = d,
                  family = Beta(), prior = pr3)
  cat("-- class sd, group g, coef Intercept, no dpar --\n")
  cat(grep("lprior", strsplit(sc3, "\n")[[1]], value = TRUE), sep = "\n")

  cat("\n-- multivariate: sd group = g, no resp --\n")
  dmv <- correct_data_mv()
  f_mv <- bf(mvbind(y1, y2) ~ x + (1 | g)) + set_rescor(TRUE)
  for (p in list(set_prior("normal(0, 5)", class = "sd", group = "g"),
                 set_prior("normal(0, 5)", class = "sd"),
                 set_prior("normal(0, 5)", class = "sd", group = "g",
                           resp = "y1"),
                 set_prior("normal(0, 5)", class = "b"),
                 set_prior("normal(0, 5)", class = "b", coef = "x"),
                 set_prior("normal(0, 5)", class = "Intercept"),
                 set_prior("normal(0, 5)", class = "sigma"))) {
    cat(sprintf("class %-9s coef %-2s group %-2s resp %-3s: ", p$class,
                p$coef, p$group, p$resp))
    r <- tryCatch({
      sc <- stancode(f_mv, data = dmv, family = gaussian(), prior = p)
      l <- grep("normal_lpdf|normal_lupdf", strsplit(sc, "\n")[[1]],
                value = TRUE)
      paste(trimws(l), collapse = " | ")
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(gsub("\n", " ", r), "\n")
  }

  cat("\n== item 4: default priors ==\n")
  dp15 <- correct_data_p15()
  b15 <- bf(mvbind(y1, y2) ~ x + (x | ID1 | g)) + set_rescor(TRUE)
  t15 <- default_prior(b15, dp15, family = gaussian())
  print(as.data.frame(t15)[, c("prior", "class", "coef", "group", "resp",
                               "dpar", "nlpar", "lb", "ub", "source")])
  dp27 <- data.frame(y = rep(c(1, 3), each = 5), off = 10)
  t27 <- default_prior(y ~ 1 + offset(off), dp27)
  print(as.data.frame(t27)[, c("prior", "class", "coef", "resp", "dpar")])
  do <- correct_data_offset()
  for (fm in list(y ~ x + offset(o), bf(y ~ x + offset(o), sigma ~ z))) {
    t <- default_prior(fm, do)
    print(as.data.frame(t)[, c("prior", "class", "coef", "dpar")])
  }
  dpc <- correct_data_offset_pois()
  t <- default_prior(cnt ~ x + offset(log(expo)), dpc, family = poisson())
  print(as.data.frame(t)[, c("prior", "class", "coef", "dpar")])
  t <- default_prior(bf(y ~ x + offset(o), sigma ~ z + offset(o2)), do)
  print(as.data.frame(t)[, c("prior", "class", "coef", "dpar")])
}

if (part %in% c("mi", "all")) {
  cat("\n== item 5: mi() names ==\n")
  dmi <- correct_data_mi()
  bmi <- bf(y ~ mi(xm) + z) + bf(xm | mi() ~ z) + set_rescor(FALSE)
  print(as.data.frame(default_prior(bmi, dmi))[, c("prior", "class",
                                                   "coef", "resp")])
  f <- bfit("mi", bmi, data = dmi)
  cat("variables:", paste(variables(f)[!grepl("^Ymi|^lp__|^lprior",
                                            variables(f))],
                          collapse = " "), "\n")
  print(fixef(f))
  print(prior_summary(f))
  sm <- summary(f)
  print(rownames(sm$fixed))
  bmi2 <- bf(y ~ mi(xm) * z) + bf(xm | mi() ~ z) + set_rescor(FALSE)
  f2 <- bfit("mi2", bmi2, data = dmi)
  cat("variables (interaction):",
      paste(grep("^b", variables(f2), value = TRUE), collapse = " "), "\n")
  print(fixef(f2))
  print(as.data.frame(default_prior(bmi2, dmi))[, c("prior", "class",
                                                    "coef", "resp")])
  cat("hypothesis on bsp name:",
      try_msg(hypothesis(f, "y_mixm > 0")), "\n")
  cat("posterior_summary rows:",
      paste(rownames(posterior_summary(f, variable = "^b",
                                       regex = TRUE)), collapse = " "), "\n")
}

if (part %in% c("theta", "all")) {
  cat("\n== item 6: a written theta formula's reference component ==\n")
  dth <- correct_data_mix()
  mix2 <- mixture(gaussian(), gaussian())
  mix3 <- mixture(gaussian(), gaussian(), gaussian())
  show_theta <- function(sc) {
    l <- strsplit(sc, "\n")[[1]]
    cat(grep("theta", l, value = TRUE), sep = "\n")
  }
  for (spec in list(list("K2 theta1 ~ x", bf(y ~ 1, theta1 ~ x), mix2),
                    list("K2 theta2 ~ x", bf(y ~ 1, theta2 ~ x), mix2),
                    list("K3 theta1 ~ x", bf(y ~ 1, theta1 ~ x), mix3),
                    list("K3 theta2,3 ~ x", bf(y ~ 1, theta2 ~ x,
                                               theta3 ~ x), mix3),
                    list("K3 theta1,2 ~ x", bf(y ~ 1, theta1 ~ x,
                                               theta2 ~ x), mix3),
                    list("K2 theta1 ~ x, theta2 ~ x", bf(y ~ 1, theta1 ~ x,
                                                         theta2 ~ x), mix2))) {
    cat("--", spec[[1]], "--\n")
    r <- tryCatch({
      show_theta(stancode(spec[[2]], data = dth, family = spec[[3]]))
      tab <- default_prior(spec[[2]], data = dth, family = spec[[3]])
      print(as.data.frame(tab)[grepl("theta", tab$dpar) |
                                 grepl("theta", tab$class),
                               c("prior", "class", "coef", "dpar")])
      "OK"
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(r, "\n")
  }
  f <- bfit("mix-theta2", bf(y ~ 1, theta2 ~ x), data = dth, family = mix2,
            init = function() list(Intercept_mu1 = -1, Intercept_mu2 = 2))
  print(fixef(f))
  cat("theta2 at x = 0, 1:",
      format(colMeans(posterior_epred(f, dpar = "theta2",
                                      newdata = data.frame(x = c(0, 1))))),
      "\n")
  cat("theta1 at x = 0, 1:",
      format(colMeans(posterior_epred(f, dpar = "theta1",
                                      newdata = data.frame(x = c(0, 1))))),
      "\n")
}
