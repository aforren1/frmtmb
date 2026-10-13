# Lane optrap: suite fits that moved on the lane build, traced.
#   Rscript dev/optrap-diag4.R <lib or ""> <case>
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
show <- function(f) {
  cat(sprintf("code %d (%s) obj %.12f iter %s evals %s\n",
              f$opt$convergence, f$opt$message, f$opt$objective,
              f$opt$iterations, paste(f$opt$evaluations, collapse = "/")))
  cat("par", sprintf("%.10g", f$opt$par), "\n")
  cat("grad*u", signif(f$obj$gr(f$opt$par) * (f$par_units %||% 1), 3), "\n")
}
if (a[2] == "iq") {
  x <- matrix(c(90, 95, 100, 105, 110, 115, 150, 155, 160), 3, 3,
              byrow = TRUE)
  d <- data.frame(y = as.vector(t(x)), person = factor(rep(1:3, each = 3)))
  f <- frm(y ~ 0 + person, data = d, verbose = 2)
  show(f)
  cat("exact sigma", sqrt(mean((d$y - rep(c(95, 110, 155), each = 3))^2)),
      "fit", exp(f$opt$par[4]), "\n")
}
if (a[2] == "w0") {
  set.seed(1)
  gm <- rnorm(10, sd = 5) + c(rep(0, 8), -50, 50)
  dw <- data.frame(group = factor(rep(LETTERS[1:10], each = 3)),
                   x = rnorm(30))
  dw$y <- gm[as.integer(dw$group)] + 10 * dw$x + rnorm(30)
  dw$w <- as.integer(dw$group %in% LETTERS[1:8])
  sub <- frm(bf(y ~ x + (1 | group)) + gaussian(),
             data = subset(dw, w == 1), verbose = 2)
  show(sub)
  wt <- frm(bf(y | weights(w) ~ x + (1 | group)) + gaussian(), data = dw,
            verbose = 2)
  show(wt)
}
if (a[2] == "sleep") {
  data(sleepstudy, package = "lme4")
  f <- frm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
           verbose = 2)
  show(f)
  l4 <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy,
                   REML = FALSE)
  cat("logLik - lme4", as.numeric(logLik(f)) - as.numeric(logLik(l4)), "\n")
}
if (a[2] == "mo") {
  src <- readLines("tests/testthat/test-optima.R")
  eval(parse(text = src[grep("^optima_mo_data <- function", src):
                          (grep("^test_that", src)[1] - 1L)]))
  for (s in as.integer(strsplit(a[3], ",")[[1]])) {
    d <- optima_mo_data(s)
    f <- frm(ls ~ mo(income) * age, data = d, verbose = 1)
    show(f)
    cat("seed", s, "gap to exact", optima_mo_exact(d) - as.numeric(logLik(f)),
        "\n")
  }
}
if (a[2] == "sleep2") {
  data(sleepstudy, package = "lme4")
  f <- frm(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  l4 <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy,
                   REML = FALSE)
  g <- glmmTMB::glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy,
                        REML = FALSE)
  show(f)
  cat("logLik - lme4", as.numeric(logLik(f)) - as.numeric(logLik(l4)),
      "| - glmmTMB", as.numeric(logLik(f)) - as.numeric(logLik(g)), "\n")
  cat("fixef - lme4", fixef_by_dpar(f)$mu - lme4::fixef(l4), "\n")
  cat("glmmTMB fixef - lme4", glmmTMB::fixef(g)$cond - lme4::fixef(l4), "\n")
  se_f <- summary(f)$coefficients$mu[, "Std. Error"]
  se_g <- summary(g)$coefficients$cond[, "Std. Error"]
  cat("se - glmmTMB", se_f - se_g, "\n")
}
if (a[2] == "m3") {
  set.seed(101)
  d3 <- data.frame(x = rnorm(300), g1 = factor("only"))
  d3$y <- rnorm(300, 1 + 2 * d3$x, 1)
  if (identical(a[3], "half")) d3$y <- d3$y / 2
  w <- character()
  f <- withCallingHandlers(frm(bf(y ~ x + (1 | g1)), family = gaussian(),
    data = d3, prior = set_prior("", class = "b", ub = 0.1), verbose = 1),
    warning = function(x) {w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")})
  show(f)
  cat("units", f$par_units, "\n")
  for (x in w) cat("W:", substr(x, 1, 200), "\n")
  print(frmtmb:::sdr_of(f)$se_lost)
}
if (a[2] == "m3b") {
  set.seed(101)
  d3 <- data.frame(x = rnorm(300), g1 = factor("only"))
  d3$y <- rnorm(300, 1 + 2 * d3$x, 1)
  for (th in c(0, log(2))) for (pr in c(TRUE, FALSE)) {
    f <- suppressWarnings(frm(bf(y ~ x + (1 | g1)), family = gaussian(),
      data = d3, prior = if (pr) set_prior("", class = "b", ub = 0.1),
      start = list(theta = th)))
    cat(sprintf("theta0 %.3f bound %s: code %d evals %s obj %.6f theta %.3f\n",
                th, pr, f$opt$convergence,
                paste(f$opt$evaluations, collapse = "/"), f$opt$objective,
                f$opt$par[length(f$opt$par)]))
  }
}
if (a[2] == "mo7") {
  src <- readLines("tests/testthat/test-optima.R")
  eval(parse(text = src[grep("^optima_mo_data <- function", src):
                          (grep("^test_that", src)[1] - 1L)]))
  for (s in as.integer(strsplit(a[3], ",")[[1]])) {
    d <- optima_mo_data(s)
    f <- frm(ls ~ mo(income) * age, data = d, verbose = 1)
    show(f)
    cat("seed", s, "gap to exact", optima_mo_exact(d) - as.numeric(logLik(f)),
        "\n")
    print(summary(f)[["mo"]])
  }
}
if (a[2] == "moface") {
  src <- readLines("tests/testthat/test-optima.R")
  eval(parse(text = src[grep("^optima_mo_data <- function", src):
                          (grep("^test_that", src)[1] - 1L)]))
  ns <- asNamespace("frmtmb")
  s <- as.integer(a[3])
  d <- optima_mo_data(s)
  f <- suppressWarnings(frm(ls ~ mo(income) * age, data = d))
  terms <- ns$mo_search_terms(f$obj, f$frame)
  f0 <- f$obj$fn(f$opt$par)
  for (tm in terms) {
    w <- ns$mo_simplex(f$opt$par[tm$z])
    cat("term z", tm$z, "b", tm$b, "w", signif(w, 4), "b value",
        f$opt$par[tm$b], "\n")
    for (j in seq_len(tm$D)) for (aa in c(1e-6, 1e-4, 1e-2, 0.1, 0.5)) {
      p <- f$opt$par
      p[tm$z] <- ns$mo_coords((1 - aa) * w + aa * replace(numeric(tm$D), j, 1))
      cat(sprintf("  j %d a %g: f0 - f %.3g\n", j, aa, f0 - f$obj$fn(p)))
    }
  }
  cat("gap", optima_mo_exact(d) - as.numeric(logLik(f)), "\n")
}
if (a[2] == "mocone") {
  src <- readLines("tests/testthat/test-optima.R")
  eval(parse(text = src[grep("^optima_mo_data <- function", src):
                          (grep("^test_that", src)[1] - 1L)]))
  ns <- asNamespace("frmtmb")
  s <- as.integer(a[3])
  d <- optima_mo_data(s)
  f <- suppressWarnings(frm(ls ~ mo(income) * age, data = d))
  ex <- optima_mo_exact(d)
  tm <- ns$mo_search_terms(f$obj, f$frame)[[2]]
  np <- length(f$opt$par)
  lo <- rep(-Inf, np); up <- rep(Inf, np)
  b <- f$opt$par[tm$b]
  if (b < 0) lo[tm$b] <- 0 else up[tm$b] <- 0
  ctl <- f$control; ctl$restarts <- 0L
  for (k in 1:3) for (mix in c(0.1, 0.01)) for (bs in c(-b, -3 * b, -10 * b)) {
    p <- f$opt$par
    p[tm$b] <- bs
    p[tm$z] <- ns$mo_coords((1 - mix) * replace(numeric(3), k, 1) + mix / 3)
    o <- ns$optimize_obj(f$obj, ctl, list(lower = lo, upper = up),
                         f$par_units, start_par = pmin(pmax(p, lo), up))
    cat(sprintf("vertex %d mix %g b0 %.4f: code %d obj %.6f (fit %.6f, exact %.6f) b %.4f w %s\n",
                k, mix, bs, o$convergence, o$objective, f$opt$objective, -ex,
                o$par[tm$b], paste(signif(ns$mo_simplex(o$par[tm$z]), 3), collapse = ",")))
  }
}
if (a[2] == "motrace") {
  src <- readLines("tests/testthat/test-optima.R")
  eval(parse(text = src[grep("^optima_mo_data <- function", src):
                          (grep("^test_that", src)[1] - 1L)]))
  ns <- asNamespace("frmtmb")
  trace("optimize_obj", where = ns, print = FALSE,
        tracer = quote(cat("  optimize_obj from b", start_par[3:4], "z",
                           signif(start_par[6:9], 3), "lower b",
                           rep_len(bounds$lower, 9)[3:4], "upper b",
                           rep_len(bounds$upper, 9)[3:4], "\n")),
        exit = quote(cat("    -> obj", returnValue()$objective, "code",
                         returnValue()$convergence, "b",
                         returnValue()$par[3:4], "\n")))
  s <- as.integer(a[3])
  d <- optima_mo_data(s)
  f <- suppressWarnings(frm(ls ~ mo(income) * age, data = d, verbose = 1))
  cat("gap", optima_mo_exact(d) - as.numeric(logLik(f)), "\n")
}
if (a[2] == "monewton") {
  src <- readLines("tests/testthat/test-optima.R")
  eval(parse(text = src[grep("^optima_mo_data <- function", src):
                          (grep("^test_that", src)[1] - 1L)]))
  ns <- asNamespace("frmtmb")
  s0 <- as.integer(a[3])
  d <- optima_mo_data(s0)
  f <- suppressWarnings(frm(ls ~ mo(income) * age, data = d,
                            control = frmtmb_control(mo_search = FALSE)))
  tm <- ns$mo_search_terms(f$obj, f$frame)[[2]]
  b <- f$opt$par[tm$b]; s <- sign(b)
  cat("units", f$par_units, "\n")
  for (k in 1:3) {
    p <- f$opt$par; p[tm$b] <- 0
    p[tm$z] <- ns$mo_coords(replace(numeric(3), k, 1))
    g0 <- f$obj$gr(p)[tm$b]
    for (db in c(1e-3, 1e-2, 0.016, 0.1)) {
      p[tm$b] <- -s * db
      h <- (f$obj$gr(p)[tm$b] - g0) / (-s * db)
      cat(sprintf("vertex %d db %g: g0 %.4g h %.4g b0 %.4g\n", k, db, g0, h,
                  -g0 / h))
    }
  }
}
