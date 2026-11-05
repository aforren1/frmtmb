# Sweep of frm_compat()'s refused rules: for each rule with a brms
# spelling, does brms 2.23.0 build AND compile the model (stancode() plus a
# cached compile, as dev/equiv-flow.R does), and does frm() refuse it?
# Seed 4040. One row per spec to dev/refuse-sweep-out/<id>.tsv.
#   Rscript dev/refuse-sweep.R <id>
# With no id it lists the specs.
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true")
S <- function(id, rule, formula, family = "gaussian()", args = "",
              brms_formula = formula, brms_family = family,
              frm_extra = "") {
  list(id = id, rule = rule, formula = formula, family = family, args = args,
       brms_formula = brms_formula, brms_family = brms_family,
       frm_extra = frm_extra)
}
mv <- function(a, b, rescor = FALSE) {
  list(frm = sprintf("mvbf(%s, %s, rescor = %s)", a, b, rescor),
       brm = sprintf("%s + %s + set_rescor(%s)", a, b, rescor))
}
SPECS <- list(
  S("cens_mixture", "cens() x mixture", "bf(yp | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x)",
    "mixture(gaussian(), gaussian())"),
  S("trunc_mixture", "trunc() x mixture",
    "bf(yp | trunc(lb = 0) ~ 1, mu1 ~ x, mu2 ~ x)",
    "mixture(gaussian(), gaussian())"),
  S("cens_simplex", "cens() x group:simplex", "bf(cbind(p1, p2, p3) | cens(ce) ~ x)",
    "dirichlet()"),
  S("trunc_simplex", "trunc() x group:simplex",
    "bf(cbind(p1, p2, p3) | trunc(lb = 0) ~ x)", "dirichlet()"),
  S("cens_categorical", "cens() x categorical", "bf(yo | cens(ce) ~ x)",
    "categorical()"),
  S("trunc_categorical", "trunc() x categorical",
    "bf(yo | trunc(ub = 5) ~ x)", "categorical()"),
  S("cens_ordinal", "cens() x group:ordinal", "bf(yo | cens(ce) ~ x)",
    "cumulative()"),
  S("trunc_ordinal", "trunc() x group:ordinal", "bf(yo | trunc(ub = 5) ~ x)",
    "cumulative()"),
  S("se_family", "se() x kind:family (skew_normal)", "bf(y | se(s) ~ x)",
    "skew_normal()"),
  S("se_lognormal", "se() x kind:family (lognormal)", "bf(yp | se(s) ~ x)",
    "lognormal()"),
  S("mi_cens", "mi() x cens()", "bf(ym | mi() + cens(ce) ~ x)"),
  S("mi_trunc", "mi() x trunc()", "bf(ym | mi() + trunc(lb = -50) ~ x)"),
  S("mi_se", "mi() x se()", "bf(ym | mi() + se(s) ~ x)"),
  S("mi_sdy_se", "mi(sdy) x se()", "bf(ym | mi(s) + se(s) ~ x)"),
  S("mi_rescor", "mi() x rescor", mv("bf(ym | mi() ~ x)", "bf(y ~ x)", TRUE)$frm,
    brms_formula = mv("bf(ym | mi() ~ x)", "bf(y ~ x)", TRUE)$brm),
  S("thres_mvbf", "thres() x mvbf",
    mv("bf(yo | thres(3) ~ x, family = cumulative())", "bf(y ~ x)")$frm, "",
    brms_formula = mv("bf(yo | thres(3) ~ x, family = cumulative())",
                      "bf(y ~ x)")$brm, brms_family = ""),
  S("ordinal_mvbf", "ordinal x mvbf",
    mv("bf(yo ~ x, family = cumulative())", "bf(y ~ x)")$frm, "",
    brms_formula = mv("bf(yo ~ x, family = cumulative())", "bf(y ~ x)")$brm,
    brms_family = ""),
  S("subset_rescor", "subset() x rescor",
    mv("bf(y | subset(sb) ~ x)", "bf(y2 ~ x)", TRUE)$frm,
    brms_formula = mv("bf(y | subset(sb) ~ x)", "bf(y2 ~ x)", TRUE)$brm),
  S("subset_me", "subset() x me()",
    mv("bf(y | subset(sb) ~ me(x, s))", "bf(y2 ~ x)")$frm,
    brms_formula = mv("bf(y | subset(sb) ~ me(x, s))", "bf(y2 ~ x)")$brm),
  S("rescor_student_mix", "rescor x kind:family (poisson)",
    mv("bf(yc ~ x, family = poisson())", "bf(y ~ x)", TRUE)$frm, "",
    brms_formula = mv("bf(yc ~ x, family = poisson())", "bf(y ~ x)", TRUE)$brm,
    brms_family = ""),
  S("rescor_cens", "rescor x cens()",
    mv("bf(y | cens(ce) ~ x)", "bf(y2 ~ x)", TRUE)$frm,
    brms_formula = mv("bf(y | cens(ce) ~ x)", "bf(y2 ~ x)", TRUE)$brm),
  S("rescor_trunc", "rescor x trunc()",
    mv("bf(y | trunc(lb = -50) ~ x)", "bf(y2 ~ x)", TRUE)$frm,
    brms_formula = mv("bf(y | trunc(lb = -50) ~ x)", "bf(y2 ~ x)", TRUE)$brm),
  S("rescor_se", "rescor x se()",
    mv("bf(y | se(s) ~ x)", "bf(y2 ~ x)", TRUE)$frm,
    brms_formula = mv("bf(y | se(s) ~ x)", "bf(y2 ~ x)", TRUE)$brm),
  S("rescor_weights", "rescor x weights()",
    mv("bf(y | weights(wt) ~ x)", "bf(y2 ~ x)", TRUE)$frm,
    brms_formula = mv("bf(y | weights(wt) ~ x)", "bf(y2 ~ x)", TRUE)$brm),
  S("id_diag", "|ID| x kind:covstruct (gr(id, cor = FALSE))",
    "bf(y ~ x + (1 + w | gr(g, id = 'p', cor = FALSE)), sigma ~ (1 | gr(g, id = 'p', cor = FALSE)))"),
  S("id_rr_na", "|ID| x group:student_blocks",
    "bf(y ~ x + (1 | p | gr(g, dist = 'student')), sigma ~ (1 | p | gr(g, dist = 'student')))"),
  S("mm_student", "mm() x group:student_blocks",
    "bf(y ~ x + (1 | mm(g1, g2, dist = 'student')))"),
  S("mm_pw", "mm(pw = )", "bf(y ~ x + (1 | mm(g1, g2, pw = cbind(wt, wt))))"),
  S("student_grcov", "group:student_blocks x gr_cov",
    "bf(y ~ x + (1 | gr(g, cov = A10, dist = 'student')))"),
  S("mixture_rescor", "mixture x rescor",
    mv("bf(y ~ x, family = mixture(gaussian(), gaussian()))", "bf(y2 ~ x)", TRUE)$frm, "",
    brms_formula = mv("bf(y ~ x, family = mixture(gaussian(), gaussian()))",
                      "bf(y2 ~ x)", TRUE)$brm, brms_family = ""),
  S("mixture_mi", "mixture x mi()", "bf(ym | mi() ~ 1, mu1 ~ x, mu2 ~ x)",
    "mixture(gaussian(), gaussian())"),
  S("mixture_se", "mixture x se()", "bf(y | se(s) ~ 1, mu1 ~ x, mu2 ~ x)",
    "mixture(gaussian(), gaussian())"),
  S("mo_mi", "mo() x mi_pred()",
    mv("bf(y ~ mo(ord) : mi(ym))", "bf(ym | mi() ~ x)")$frm,
    brms_formula = mv("bf(y ~ mo(ord) : mi(ym))", "bf(ym | mi() ~ x)")$brm),
  S("bar_crossing", "bar_crossing x *", "bf(y ~ x * (1 | g))"),
  S("autocor_aterm_cens", "group:autocor x kind:aterm (cens)",
    "bf(y | cens(ce) ~ x + ar(p = 1, cov = TRUE))"),
  S("autocor_aterm_weights", "group:autocor x kind:aterm (weights)",
    "bf(y | weights(wt) ~ x + ar(p = 1, cov = TRUE))"),
  S("autocor_aterm_trunc", "group:autocor x kind:aterm (trunc)",
    "bf(y | trunc(lb = -50) ~ x + ar(p = 1, cov = TRUE))"),
  S("autocor_se_cov", "ar() x se() (cov = TRUE)",
    "bf(y | se(s) ~ x + ar(p = 1, cov = TRUE))"),
  S("autocor_se_nocov", "ar() x se() (cov = FALSE)",
    "bf(y | se(s) ~ x + ar(p = 1))"),
  S("autocor_rescor", "group:autocor x rescor",
    mv("bf(y ~ x + ar(p = 1, cov = TRUE))", "bf(y2 ~ x + ar(p = 1, cov = TRUE))", TRUE)$frm,
    brms_formula = mv("bf(y ~ x + ar(p = 1, cov = TRUE))",
                      "bf(y2 ~ x + ar(p = 1, cov = TRUE))", TRUE)$brm),
  S("autocor_mixture", "group:autocor x mixture",
    "bf(y ~ 1, mu1 ~ x + ar(p = 1, cov = TRUE), mu2 ~ x)",
    "mixture(gaussian(), gaussian())"),
  S("autocor_two", "group:autocor x kind:autocor",
    "bf(y ~ x + ar(p = 1) + cosy())"),
  S("grby_grcov", "gr_by x gr_cov",
    "bf(y ~ x + (1 | gr(g, by = f10, cov = A10)))"),
  S("mmby_cov", "mm(by = , cov = )",
    "bf(y ~ x + (1 | mm(g1, g2, by = cbind(fm1, fm2), cov = A)))"),
  S("grcov_diag", "gr_cov with || and several columns",
    "bf(y ~ x + (1 + w || gr(g, cov = A10)))"),
  S("cs_family", "cs_pred() x kind:family", "bf(y ~ cs(x))"),
  S("zero_intercept_ord", "0 + Intercept x group:ordinal",
    "bf(yo ~ 0 + Intercept + x)", "cumulative()"),
  S("cens_rate", "rate() x kind:family (gaussian)", "bf(y | rate(wt) ~ x)"),
  S("disc_gaussian", "disc x kind:family", "bf(y ~ x, disc ~ 1)"),
  S("thres_gaussian", "thres() x kind:family", "bf(y | thres(2) ~ x)"),
  S("cens_disc", "p010 cens() x disc", "bf(yo | cens(ce) ~ x, disc ~ 1)",
    "cumulative()"),
  S("cens_equidistant", "p011 cens() x equidistant", "bf(yo | cens(ce) ~ x)",
    "cumulative(threshold = 'equidistant')"),
  S("vint_stz", "p059 vint() x sum_to_zero", "bf(yo | vint(ce) ~ x)",
    "cumulative(threshold = 'sum_to_zero')"),
  S("trials_cs", "p146 trials() x cs_pred()", "bf(yb | trials(5) ~ cs(x))",
    "binomial()"),
  S("trials_cs2", "p146 trials() x cs_pred() on ordinal", "bf(yo | trials(5) ~ cs(x))",
    "cumulative()"),
  S("mi_pred_gaussian_student_x", "mi_pred under poisson",
    mv("bf(yc ~ mi(ym) + x, family = poisson())", "bf(ym | mi() ~ z)")$frm, "",
    brms_formula = mv("bf(yc ~ mi(ym) + x, family = poisson())", "bf(ym | mi() ~ z)")$brm,
    brms_family = "")
)
ids <- vapply(SPECS, `[[`, "", "id")
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  cat(ids, sep = "\n")
  quit(save = "no")
}
spec <- SPECS[[match(args[1L], ids)]]
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
attach(list(me = brms::me, mmc = brms::mmc, mm = brms::mm, gr = brms::gr,
            mo = brms::mo, mi = brms::mi),
       name = "refuse_brms_specials", warn.conflicts = FALSE)
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
set.seed(4040)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                g = factor(rep(sprintf("g%02d", 1:10), each = 20)),
                t = rep(1:20, 10))
d$y <- 1 + d$x + rnorm(10)[d$g] + rnorm(n, 0, 0.5)
d$y2 <- d$x + rnorm(n, 0, 0.5)
d$yp <- exp(0.3 * d$x + rnorm(n, 0, 0.3))
d$yc <- rpois(n, exp(0.5 + 0.3 * d$x))
d$yb <- rbinom(n, 5, 0.4)
d$yo <- findInterval(d$x + rnorm(n), c(-1, 0, 1)) + 1L
d$ym <- d$y
d$ym[sample(n, 15)] <- NA
d$ce <- sample(c(0L, 0L, 0L, 1L), n, TRUE)
d$s <- runif(n, 0.2, 0.5)
d$wt <- runif(n, 0.5, 2)
d$sb <- rep(c(TRUE, TRUE, TRUE, FALSE), length.out = n)
d$ord <- sample(1:4, n, TRUE)
d$g1 <- factor(sample(letters[1:8], n, TRUE), levels = letters[1:8])
d$g2 <- factor(sample(letters[1:8], n, TRUE), levels = letters[1:8])
d$fm1 <- ifelse(d$g1 %in% letters[1:4], "u", "v")
d$fm2 <- ifelse(d$g2 %in% letters[1:4], "u", "v")
d$f10 <- factor(ifelse(as.integer(d$g) <= 5, "a", "b"))
p <- MASS::mvrnorm(n, rep(0, 3), diag(3))
p <- exp(p) / rowSums(exp(p))
d$p1 <- p[, 1]; d$p2 <- p[, 2]; d$p3 <- p[, 3]
A <- 0.5^abs(outer(1:8, 1:8, "-"))
dimnames(A) <- list(letters[1:8], letters[1:8])
A10 <- 0.5^abs(outer(1:10, 1:10, "-"))
lv <- levels(d$g)
dimnames(A10) <- list(lv, lv)
d2 <- list(A = A, A10 = A10)
one_line <- function(x) substr(gsub("[\r\n\t]+", " ", x), 1, 260)
fr <- tryCatch({
  fm <- eval(parse(text = spec$formula), asNamespace("frmtmb"))
  fa <- if (nzchar(spec$family)) {
    list(family = eval(parse(text = spec$family), asNamespace("frmtmb")))
  } else {
    list()
  }
  fit <- suppressWarnings(do.call(frm, c(list(fm), fa,
                                         list(data = d, data2 = d2))))
  paste0("fits, logLik ", format(as.numeric(logLik(fit)), digits = 8))
}, error = function(e) paste0("REFUSED: ", one_line(conditionMessage(e))))
br <- tryCatch({
  bf_env <- asNamespace("brms")
  fm <- eval(parse(text = spec$brms_formula), bf_env)
  fa <- if (nzchar(spec$brms_family)) {
    list(family = eval(parse(text = spec$brms_family), bf_env))
  } else {
    list()
  }
  code <- suppressWarnings(do.call(brms::stancode,
                                   c(list(fm), fa, list(data = d, data2 = d2))))
  cc <- tryCatch({
    suppressWarnings(E$brms_stan_model(code))
    "builds and compiles"
  }, error = function(e) paste0("builds, COMPILE FAILS: ",
                                one_line(conditionMessage(e))))
  cc
}, error = function(e) paste0("brms refuses: ", one_line(conditionMessage(e))))
out <- data.frame(id = spec$id, rule = spec$rule, frm = fr, brms = br)
dir.create("dev/refuse-sweep-out", showWarnings = FALSE)
utils::write.table(out, file.path("dev/refuse-sweep-out", paste0(spec$id, ".tsv")),
                   sep = "\t", row.names = FALSE)
cat(spec$id, "\n  frm: ", fr, "\n  brms:", br, "\n")
