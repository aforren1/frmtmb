# Area 4: spellings of bf(nl = TRUE, autocor =) against brms, refits and readers.
source("dev/review-nlfix-common.R")
set.seed(4401)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 20)), t = rep(1:20, 10))
d$y <- 1 + 0.5 * d$x + ave(rnorm(n, 0, 0.5), d$g, FUN = function(z) as.numeric(stats::filter(z, 0.5, "recursive")))
d$yc <- rpois(n, exp(0.5 + 0.2 * d$x))
d$yo <- cut(d$y, c(-Inf, 0, 1, 2, Inf), labels = 1:4, ordered_result = TRUE)
d$cat <- factor(sample(c("a", "b", "c"), n, TRUE))
d$pos <- exp(d$y / 3)
fb <- function(txt, fam = "gaussian()") {
  r <- tryCatch(suppressWarnings(suppressMessages(brms::stancode(eval(parse(text = txt), envir = asNamespace("brms")), data = d, family = eval(parse(text = fam), envir = asNamespace("brms"))))), error = function(e) e)
  if (inherits(r, "error")) return(paste("REFUSES:", substr(gsub("\n", " ", conditionMessage(r)), 1, 220)))
  ln <- strsplit(as.character(r), "\n")[[1]]
  paste("builds;", paste(grep("Err\\[n, 1|cholesky_cor|normal_time|ar\\[|cosy|mu\\[n\\] =", ln, value = TRUE)[1:3], collapse = " | "))
}
ff <- function(txt, fam = "gaussian()", ...) {
  r <- tt(frm(eval(parse(text = txt)), family = eval(parse(text = fam)), data = d, ...))
  if (iserr(r)) return(paste("REFUSES:", substr(r$msg, 1, 300)))
  ac <- r$frame[["autocor"]]
  paste0("ll=", format(as.numeric(logLik(r)), digits = 10), " ac=", paste(names(ac), collapse = ","), " label=", paste(vapply(ac, function(a) a[["label"]] %||% "?", ""), collapse = ","),
         if (length(attr(r, "w"))) paste0(" WARN[", substr(attr(r, "w")[1], 1, 100), "]"))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
st <- list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 0.5))
sp <- list(
  arg = c("bf(y ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE))", "gaussian()"),
  plus = c("bf(y ~ a + b * x, a + b ~ 1, nl = TRUE) + acformula(~ ar(time = t, gr = g, cov = TRUE))", "gaussian()"),
  plus2 = c("bf(y ~ a + b * x, a + b ~ 1, nl = TRUE) + acformula(~ ar(time = t, gr = g, cov = TRUE)) + acformula(~ cosy(time = t, gr = g))", "gaussian()"),
  arg_plus = c("bf(y ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE)) + acformula(~ cosy(time = t, gr = g))", "gaussian()"),
  two_terms = c("bf(y ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE) + cosy(time = t, gr = g))", "gaussian()"),
  linear_plus2 = c("bf(y ~ x) + acformula(~ ar(time = t, gr = g, cov = TRUE)) + acformula(~ cosy(time = t, gr = g))", "gaussian()"),
  nlf_mu = c("bf(y ~ 1, autocor = ~ ar(time = t, gr = g, cov = TRUE)) + nlf(mu ~ a + b * x) + lf(a ~ 1) + lf(b ~ 1)", "gaussian()"),
  nlf_mu_plus = c("bf(y ~ 1) + nlf(mu ~ a + b * x) + lf(a ~ 1) + lf(b ~ 1) + acformula(~ ar(time = t, gr = g, cov = TRUE))", "gaussian()"),
  term_on_rhs_first = c("bf(y ~ 1 + ar(time = t, gr = g, cov = TRUE)) + nlf(mu ~ a + b * x) + lf(a ~ 1) + lf(b ~ 1)", "gaussian()"),
  nl_later = c("bf(y ~ a + b * x, autocor = ~ ar(time = t, gr = g, cov = TRUE)) + lf(a ~ 1, b ~ 1, nl = TRUE)", "gaussian()"),
  poisson_cov = c("bf(yc ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE))", "poisson()"),
  poisson_nocov = c("bf(yc ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g))", "poisson()"),
  cumul_cov = c("bf(yo ~ b * x, b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE))", "cumulative()"),
  gamma_cov = c("bf(pos ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE))", "Gamma('log')"),
  categ = c("bf(cat ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE))", "categorical()"),
  mixture = c("bf(y ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE))", "mixture(gaussian(), gaussian())"),
  sigma_ac = c("bf(y ~ a + b * x, a + b ~ 1, sigma ~ z, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = FALSE))", "gaussian()"),
  re_ac = c("bf(y ~ a + b * x, a ~ 1 + (1 | g), b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE))", "gaussian()"),
  cond_var = c("bf(y ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE, p = 2))", "gaussian()")
)
for (nm in names(sp)) {
  cat("\n==", nm, ":", sp[[nm]][1], "|", sp[[nm]][2], "\n")
  cat("  brms:", fb(sp[[nm]][1], sp[[nm]][2]), "\n")
  cat("  frm :", ff(sp[[nm]][1], sp[[nm]][2]), "\n")
}
cat("\n--- bf_update path\n")
b0 <- bf(y ~ a + b * x, a + b ~ 1, nl = TRUE)
b1 <- bf(b0, autocor = ~ ar(time = t, gr = g, cov = TRUE)); print(b1)
b2 <- bf(b1, autocor = ~ cosy(time = t, gr = g)); print(b2)
bl <- bf(bf(y ~ x, autocor = ~ ar(time = t, gr = g, cov = TRUE)), autocor = ~ cosy(time = t, gr = g)); print(bl)
br1 <- brms::bf(y ~ a + b * x, a + b ~ 1, nl = TRUE)
br2 <- brms::bf(br1, autocor = ~ ar(time = t, gr = g, cov = TRUE))
br3 <- brms::bf(br2, autocor = ~ cosy(time = t, gr = g)); print(br3$autocor)
cat("\n--- readers and refits on one autocor body\n")
for (ter in c("ar(time = t, gr = g, cov = TRUE)", "ar(time = t, gr = g, cov = FALSE)", "ma(time = t, gr = g, cov = TRUE)", "cosy(time = t, gr = g)")) {
  cat("\n term", ter, "\n")
  fn <- suppressWarnings(frm(eval(parse(text = sprintf("bf(y ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ %s)", ter))), data = d, start = st))
  fl <- suppressWarnings(frm(eval(parse(text = sprintf("bf(y ~ x, autocor = ~ %s)", ter))), data = d))
  cat(" logLik", as.numeric(logLik(fn)), as.numeric(logLik(fl)), "\n")
  nd <- d[1:30, ]
  chk <- function(lab, f) {
    a <- tt(f(fn)); b <- tt(f(fl))
    if (iserr(a) || iserr(b)) cat(sprintf("  %-22s nl=%s | lin=%s\n", lab, if (iserr(a)) paste("ERR", substr(a$msg, 1, 170)) else "ok", if (iserr(b)) paste("ERR", substr(b$msg, 1, 170)) else "ok"))
    else { a <- as.numeric(unlist(a)); b <- as.numeric(unlist(b)); cat(sprintf("  %-22s %s\n", lab, if (length(a) == length(b)) sprintf("max diff %.3e", max(abs(a - b))) else sprintf("lengths %d %d", length(a), length(b)))) }
  }
  chk("fitted", function(f) fitted(f))
  chk("fitted newdata", function(f) fitted(f, newdata = nd))
  chk("linpred", function(f) frm_linpred(f))
  chk("linpred newdata", function(f) frm_linpred(f, newdata = nd))
  chk("predict", function(f) { set.seed(1); predict(f) })
  chk("predict newdata", function(f) { set.seed(1); predict(f, newdata = nd) })
  chk("simulate", function(f) { set.seed(1); as.numeric(as.matrix(simulate(f, nsim = 2))) })
  chk("simulate newdata", function(f) { set.seed(1); as.numeric(as.matrix(simulate(f, nsim = 2, newdata = nd))) })
  chk("residuals", function(f) residuals(f))
  chk("resid pearson", function(f) residuals(f, type = "pearson"))
  chk("osa", function(f) residuals(f, type = "osa"))
  chk("update data", function(f) as.numeric(logLik(update(f, data = d[d$g != 1, ]))))
  chk("refit", function(f) as.numeric(logLik(update(f))))
  chk("confint_varcorr", function(f) as.data.frame(confint_varcorr(f))[, 2:3])
  chk("autocor_matrix", function(f) autocor_matrix(f))
  chk("frm_curve", function(f) as.data.frame(frm_curve(f, newdata = nd))$.estimate)
  chk("emmeans", function(f) as.data.frame(emmeans::emmeans(f, ~ x, at = list(x = 0:1)))$emmean)
  chk("predictions", function(f) marginaleffects::predictions(f, newdata = nd)$estimate)
  chk("cond eff", function(f) nrow(conditional_effects(f, effects = "x")[[1]]))
  chk("get_prior", function(f) nrow(get_prior(f)))
  x <- tt(print(formula(fn))); cat("  formula ok\n")
}
