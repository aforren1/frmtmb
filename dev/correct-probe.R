# Lane wt-correct: what frmtmb does for each of the seven items, on the
# data dev/correct-brms-probe.R measures brms on.
#
# Usage: Rscript dev/correct-probe.R <arm> <part>
#   arm  = base (rellib-r3 only) | lane (correct-lib first)
#   part = resid | ppcheck | priors | mi | theta | all
arm <- commandArgs(trailingOnly = TRUE)[1]
part <- commandArgs(trailingOnly = TRUE)[2]
if (is.na(part)) part <- "all"
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-prelude.R")
if (identical(arm, "base")) .libPaths(.libPaths()[-1L])
suppressMessages(library(frmtmb))
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-data.R")

try_msg <- function(expr) {
  r <- tryCatch({
    v <- force(expr)
    if (inherits(v, "ggplot")) invisible(ggplot2::ggplot_build(v))
    "OK"
  }, error = function(e) {
    paste0("ERROR [", paste(setdiff(class(e), c("error", "condition")),
                            collapse = ","), "]: ", conditionMessage(e))
  })
  substr(gsub("\n", " ", r), 1, 300)
}

if (part %in% c("resid", "all")) {
  cat("\n== item 1: residuals() by family ==\n")
  d <- correct_data_ord()
  fits <- list(
    cumulative = frm(bf(ord ~ x) + cumulative(), data = d),
    sratio = frm(bf(ord ~ x) + sratio(), data = d),
    categorical = frm(bf(cat ~ x) + categorical(), data = d),
    gaussian = frm(bf(y ~ x) + gaussian(), data = d)
  )
  dm <- correct_data_multinom()
  fits$multinomial <- frm(bf(Y | trials(n) ~ x) + multinomial(K = 3), data = dm)
  for (nm in names(fits)) {
    f <- fits[[nm]]
    for (ty in c("response", "ordinary", "pearson", "deviance", "osa")) {
      cat(sprintf("%-12s residuals(type = %-10s) %s\n", nm,
                  paste0("\"", ty, "\""),
                  try_msg(suppressWarnings(residuals(f, type = ty)))))
    }
    cat(sprintf("%-12s pp_check(type = \"error_binned\") %s\n", nm,
                try_msg(pp_check(f, type = "error_binned", ndraws = 5))))
    cat(sprintf("%-12s pp_check(type = \"error_hist\") %s\n", nm,
                try_msg(pp_check(f, type = "error_hist", ndraws = 5))))
    cat(sprintf("%-12s plot() %s\n", nm, try_msg({
      grDevices::pdf(NULL); on.exit(grDevices::dev.off()); plot(f)
    })))
  }
}

if (part %in% c("ppcheck", "all")) {
  cat("\n== item 2: pp_check(), every type bayesplot exposes ==\n")
  d <- correct_data_gauss()
  f <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d)
  ty <- sub("^ppc_", "", as.character(bayesplot::available_ppc("")))
  ns <- asNamespace("bayesplot")
  for (t in ty) {
    fa <- names(formals(get(paste0("ppc_", t), ns)))
    args <- list(f, type = t, ndraws = 20)
    if ("group" %in% fa) args$group <- "g"
    if ("x" %in% fa) args$x <- "x"
    set.seed(1)
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
  cat("dens_overlay, group = g:",
      try_msg(pp_check(f, type = "dens_overlay", group = "g", ndraws = 5)),
      "\n")
}

if (part %in% c("priors", "all")) {
  cat("\n== item 3: sd prior scope ==\n")
  d <- correct_data_beta()
  fit <- suppressWarnings(frm(bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta(),
                              data = d))
  show <- function(label, pl, on = fit) {
    r <- tryCatch({
      e <- frmtmb:::resolve_prior_input(on, pl)$entries
      paste(vapply(e, function(z) {
        paste0(z$comp, z$idx, "=", z$dist$kind, "(",
               paste(unlist(z$dist[-1]), collapse = ","), ")")
      }, ""), collapse = "  ")
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(sprintf("%-40s %s\n", label, r))
  }
  cat("blocks (theta idx):", vapply(fit$frame$re_blocks, function(b) {
    paste0(b$group_name, "/", paste(b$theta_idx, collapse = ","))
  }, ""), "\n")
  show("sd group g", set_prior("normal(0, 5)", class = "sd", group = "g"))
  show("sd", set_prior("normal(0, 5)", class = "sd"))
  show("sd group g coef Intercept",
       set_prior("normal(0, 5)", class = "sd", group = "g",
                 coef = "Intercept"))
  show("sd dpar phi", set_prior("normal(0, 5)", class = "sd", dpar = "phi"))
  vp <- validate_prior(set_prior("normal(0, 5)", class = "sd", group = "g"),
                       bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta(), data = d)
  print(as.data.frame(vp)[vp$class == "sd", c("prior", "class", "coef",
                                             "group", "dpar", "source")])

  cat("\n-- multivariate: sd group = g, no resp --\n")
  dmv <- correct_data_mv()
  f_mv <- bf(mvbind(y1, y2) ~ x + (1 | g)) + set_rescor(TRUE)
  fmv <- frm(f_mv + gaussian(), data = dmv)
  for (p in list(set_prior("normal(0, 5)", class = "sd", group = "g"),
                 set_prior("normal(0, 5)", class = "sd"),
                 set_prior("normal(0, 5)", class = "sd", group = "g",
                           resp = "y1"),
                 set_prior("normal(0, 5)", class = "b"),
                 set_prior("normal(0, 5)", class = "b", coef = "x"),
                 set_prior("normal(0, 5)", class = "Intercept"),
                 set_prior("normal(0, 5)", class = "sigma"))) {
    s <- unclass(p)[[1]]
    cl <- if (isTRUE(s$natural)) s$dpar else s$class
    show(sprintf("mv class %s coef %s group %s resp %s", cl, s$coef,
                 s$group, s$resp %||% ""), p, on = fmv)
  }

  cat("\n== item 4: default priors, route = \"sample\" ==\n")
  suppressMessages(library(frmtmb.sample))
  dp15 <- correct_data_p15()
  b15 <- bf(mvbind(y1, y2) ~ x + (x | ID1 | g)) + set_rescor(TRUE)
  t15 <- default_prior(b15, dp15, family = gaussian(), route = "sample")
  print(as.data.frame(t15)[, c("prior", "class", "coef", "group", "resp",
                               "dpar", "nlpar")])
  dp27 <- data.frame(y = rep(c(1, 3), each = 5), off = 10)
  t27 <- default_prior(y ~ 1 + offset(off), dp27, route = "sample")
  print(as.data.frame(t27)[, c("prior", "class", "coef", "resp", "dpar")])
  do <- correct_data_offset()
  for (fm in list(y ~ x + offset(o), bf(y ~ x + offset(o), sigma ~ z))) {
    t <- default_prior(fm, do, route = "sample")
    print(as.data.frame(t)[, c("prior", "class", "coef", "dpar")])
  }
  dpc <- correct_data_offset_pois()
  t <- default_prior(bf(cnt ~ x + offset(log(expo))) + poisson(), dpc,
                     route = "sample")
  print(as.data.frame(t)[, c("prior", "class", "coef", "dpar")])
  t <- default_prior(bf(y ~ x + offset(o), sigma ~ z + offset(o2)), do,
                     route = "sample")
  print(as.data.frame(t)[, c("prior", "class", "coef", "dpar")])
}

if (part %in% c("mi", "all")) {
  cat("\n== item 5: mi() names ==\n")
  dmi <- correct_data_mi()
  bmi <- bf(y ~ mi(xm) + z) + bf(xm | mi() ~ z) + set_rescor(FALSE)
  print(as.data.frame(default_prior(bmi, dmi))[, c("prior", "class",
                                                   "coef", "resp")])
  f <- frm(bmi + gaussian(), data = dmi)
  v <- variables(f)
  cat("variables:", paste(v, collapse = " "), "\n")
  print(fixef(f))
  sm <- summary(f)
  print(rownames(sm$fixed))
  print(rownames(vcov(f)))
  bmi2 <- bf(y ~ mi(xm) * z) + bf(xm | mi() ~ z) + set_rescor(FALSE)
  f2 <- frm(bmi2 + gaussian(), data = dmi)
  cat("variables (interaction):", paste(variables(f2), collapse = " "), "\n")
  print(fixef(f2))
  print(as.data.frame(default_prior(bmi2, dmi))[, c("prior", "class",
                                                    "coef", "resp")])
  cat("hypothesis y_mixm > 0:", try_msg(hypothesis(f, "y_mixm > 0")), "\n")
  cat("prior class b coef mixm resp y:",
      try_msg(frm(bmi + gaussian(), data = dmi,
                  prior = set_prior("normal(0, 1)", class = "b",
                                    coef = "mixm", resp = "y"))), "\n")
}

if (part %in% c("theta", "all")) {
  cat("\n== item 6: a written theta formula's reference component ==\n")
  dth <- correct_data_mix()
  mix2 <- mixture(gaussian(), gaussian())
  mix3 <- mixture(gaussian(), gaussian(), gaussian())
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
      tab <- default_prior(spec[[2]], data = dth, family = spec[[3]])
      print(as.data.frame(tab)[grepl("theta", tab$dpar) |
                                 grepl("theta", tab$class),
                               c("prior", "class", "coef", "dpar")])
      "OK"
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(r, "\n")
  }
  f <- tryCatch(frm(bf(y ~ 1, theta2 ~ x) + mix2, data = dth),
                error = function(e) e)
  if (inherits(f, "error")) {
    cat("fit K2 theta2 ~ x: ERROR:", conditionMessage(f), "\n")
  } else {
    print(fixef(f))
    nd <- data.frame(x = c(0, 1))
    for (dp in c("theta1", "theta2")) {
      cat(dp, "at x = 0, 1:", try_msg(print(frm_linpred(
        f, newdata = nd, dpar = dp, type = "response"))), "\n")
    }
  }
  f1 <- frm(bf(y ~ 1, theta1 ~ x) + mix2, data = dth)
  print(fixef(f1))
}
