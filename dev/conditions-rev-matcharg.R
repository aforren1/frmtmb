# Reviewer, lane wt-conditions (recheck 1): frm_match_arg() against
# match.arg() on the same input, in the shapes its sites use and in the
# frames that could mislead it; then the frmtmb_package tag through
# serialization, update(), modification and comparison.
#   Rscript dev/conditions-rev-matcharg.R
.libPaths(c("C:/Users/adf44/source/r/conditions-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)
  library(frmtmb.latent); library(frmtmb.spline)})
res <- function(expr) tryCatch(list(ok = TRUE, v = expr),
                               error = function(e) list(ok = FALSE,
                                                        v = class(e)[1L]))
same <- function(a, b) identical(a$ok, b$ok) && (!a$ok || identical(a$v, b$v))

cat("== 1. direct: one formal with a default vector\n")
fb <- function(type = c("link", "response", "terms")) match.arg(type)
fl <- function(type = c("link", "response", "terms")) frm_match_arg(type)
fbs <- function(type = c("link", "response", "terms")) match.arg(type, several.ok = TRUE)
fls <- function(type = c("link", "response", "terms")) frm_match_arg(type, several.ok = TRUE)
fbc <- function(m = "a") match.arg(m, c("alpha", "beta", "al"))
flc <- function(m = "a") frm_match_arg(m, c("alpha", "beta", "al"))
inputs <- list("<missing>", NULL, "link", "response", "resp", "r", "l",
               "te", "", "LINK", c("link", "response", "terms"),
               c("response", "link"), c("link", "zz"), character(0), 1,
               NA_character_, NA, c("terms", "link", "response"),
               factor("link"), list("link"))
n <- 0; bad <- 0
for (i in seq_along(inputs)) {
  miss <- i == 1L
  a <- if (miss) "<missing>" else inputs[[i]]
  pairs <- list(c("fb", "fl"), c("fbs", "fls"), c("fbc", "flc"))
  for (p in pairs) {
    call_b <- if (miss) call(p[1]) else call(p[1], a)
    call_l <- if (miss) call(p[2]) else call(p[2], a)
    rb <- res(eval(call_b)); rl <- res(eval(call_l))
    n <- n + 1
    if (!same(rb, rl)) {
      bad <- bad + 1
      cat(sprintf("  DIFF %s(%s): base %s %s | lane %s %s\n", p[2],
                  paste(deparse(a), collapse = ""), rb$ok,
                  paste(format(rb$v), collapse = ","), rl$ok,
                  paste(format(rl$v), collapse = ",")))
    }
  }
}
cat("  compared", n, "differ", bad, "\n")

cat("\n== 2. frames: S3 method, nested function, do.call, argument reassigned\n")
gen <- function(x, ...) UseMethod("gen")
gen.default <- function(x, type = c("aa", "bb"), ...) frm_match_arg(type)
gen.foo <- function(x, type = c("cc", "dd"), ...) NextMethod()
cat("  S3 default          ", format(res(gen(1, "b"))$v), "\n")
cat("  NextMethod from foo ", format(res(gen(structure(1, class = "foo")))$v),
    " (match.arg:", format(res((function(x) { g2 <- function(x, ...) UseMethod("g2"); g2.default <- function(x, type = c("aa", "bb"), ...) match.arg(type); g2.foo <- function(x, type = c("cc", "dd"), ...) NextMethod(); g2(structure(1, class = "foo")) })(1))$v), ")\n")
outer_fn <- function(type = c("outer1", "outer2")) {
  inner <- function(type = c("inner1", "inner2")) frm_match_arg(type)
  c(inner(), frm_match_arg(type))
}
cat("  nested              ", outer_fn(), "\n")
cat("  do.call             ", do.call(fl, list("resp")), "\n")
reassign_b <- function(type = c("x1", "x2")) { type <- "x2"; match.arg(type) }
reassign_l <- function(type = c("x1", "x2")) { type <- "x2"; frm_match_arg(type) }
cat("  reassigned          base", reassign_b(), " lane", reassign_l(), "\n")
lapply_l <- function(type = c("p", "q")) vapply("q", function(z) frm_match_arg(type), "")
cat("  inside lapply FUN   ", format(res(lapply_l())$v), " (match.arg:",
    format(res((function(type = c("p", "q")) vapply("q", function(z) match.arg(type), ""))())$v), ")\n")
wrapped <- function(type = c("w1", "w2")) {
  helper <- function(t2) frm_match_arg(t2, c("w1", "w2"))
  helper(type)
}
cat("  explicit choices in helper, full vector passed:", format(res(wrapped())$v), "\n")

cat("\n== 3. the real sites, valid values through each\n")
set.seed(1)
d <- data.frame(y = rnorm(40), x = rnorm(40), g = gl(8, 5))
fit <- frm(y ~ x + (1 | g), data = d)
chk <- list(
  predict_default = quote(predict(fit)),
  predict_resp_partial = quote(predict(fit, type = "resp")),
  fitted_linear = quote(fitted(fit, scale = "linear")),
  residuals_pearson = quote(residuals(fit, type = "pearson")),
  confint_Wald = quote(confint(fit, method = "Wald")),
  confint_profile = quote(confint(fit, parm = "x", method = "profile")),
  hypothesis_scope_ranef = quote(hypothesis(fit, "Intercept = 0", scope = "ranef")),
  ce_method_predict = quote(conditional_effects(fit, method = "predict")),
  control_ignore = quote(frmtmb_control(check_nlev_1 = "ignore", check_olre = "stop")),
  get_prior_default = quote(get_prior(y ~ x, data = d)),
  vcov_cluster_CR1 = quote(vcov_cluster(fit, ~ g, type = "CR1")),
  periodogram = quote(frm_periodogram(rnorm(64), taper = "hann", detrend = "linear")),
  royston = quote(royston_parmar(scale = "odds")),
  hmm_init = quote(hmm(2, init = "uniform")),
  gddm_control = quote(gddm_control(tridiagonal = "atomic")))
for (nm in names(chk)) cat(sprintf("  %-24s %s\n", nm,
  tryCatch({eval(chk[[nm]]); "OK"}, error = function(e) paste("ERROR", conditionMessage(e)))))

cat("\n== 4. frmtmb_package tag\n")
w <- wiener(max_ndt = 0.3)
h <- hmm(K = 2, gaussian(), time = t, group = id)
g <- gaussian()
show_tag <- function(lbl, f) cat(sprintf("  %-40s %s\n", lbl,
  paste(format(frm_family_package(f)), collapse = ",")))
show_tag("wiener()", w)
show_tag("hmm()", h)
show_tag("gaussian()", g)
tf <- tempfile(fileext = ".rds"); saveRDS(w, tf)
show_tag("wiener() after saveRDS/readRDS", readRDS(tf))
w2 <- w; w2$link <- "log"; show_tag("wiener() after $<-", w2)
show_tag("wiener() after modifyList", utils::modifyList(w, list(zz = 1)))
show_tag("wiener() after [ subset", structure(unclass(w)[names(w)], class = class(w)))
show_tag("mixture(wiener(), wiener())", tryCatch(mixture(w, w), error = function(e) "mixture error"))
set.seed(2)
dd <- ddm_simulate(200, mu = 1, bs = 1.5, ndt = 0.25)
fw <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5), family = wiener(), data = dd)
show_tag("fit$family of a wiener fit", fw$family)
fam_in_spec <- fw$spec$responses[[1]]$family
show_tag("spec family of a wiener fit", fam_in_spec)
tf2 <- tempfile(fileext = ".rds"); saveRDS(fw, tf2)
show_tag("spec family after saveRDS of fit", readRDS(tf2)$spec$responses[[1]]$family)
fw_u <- update(fw, data = dd[1:150, ])
show_tag("spec family after update()", fw_u$spec$responses[[1]]$family)
e <- tryCatch(residuals(fw, type = "deviance"), error = identity)
cat("  refusal after fit: ", paste(class(e), collapse = "/"), "\n")
e2 <- tryCatch(residuals(readRDS(tf2), type = "deviance"), error = identity)
cat("  same after readRDS:", paste(class(e2), collapse = "/"), "\n")
cat("  identical(wiener(), wiener()):", identical(wiener(), wiener()),
    " all.equal:", isTRUE(all.equal(wiener(), wiener())), "\n")
cat("  identical(gaussian(), gaussian()):", identical(gaussian(), gaussian()), "\n")
uf <- frmtmb_family(family = "myfam", dpars = "mu",
                    links = list(mu = "identity"),
                    lpdf = function(y, mu) -0.5 * (y - mu)^2)
show_tag("user family at top level", uf)
cat("  user family attribute present:", !is.null(attr(uf, "frmtmb_package")), "\n")
