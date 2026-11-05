# Item 6: the sampler route. D1 and D7 on frm_sample(). Arm "before" drops
# the lane's library (NLFIX_LIB unset below). Seeds 49, 158.
lib <- Sys.getenv("NLFIX_LIB",
                  "C:/Users/adf44/source/r/wt-nlfix-lib2;C:/Users/adf44/source/r/rellib-r16")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("core:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\nsample:", find.package("frmtmb.sample"), "\n")
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("dev/equiv-specs.R", envir = E)
d <- E$eq_dat("ordinal", seed = 49)
run <- function(lab, expr) {
  cat("\n==", lab, "\n")
  r <- tryCatch(suppressWarnings(suppressMessages(expr)), error = function(e) e)
  if (inherits(r, "error")) {
    cat("  ERROR:", substr(conditionMessage(r), 1, 240), "\n")
    return(invisible(NULL))
  }
  invisible(r)
}
s <- run("D1 nl cumulative", frm_sample(
  bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), data = d, family = cumulative(),
  chains = 1, iter = 500, warmup = 250, seed = 49, refresh = 0))
if (!is.null(s)) {
  cat("variables:", paste(head(variables(s), 12), collapse = ", "), "\n")
  m <- as.matrix(s)
  cat("posterior means:", format(colMeans(m)[seq_len(min(5, ncol(m)))],
                                 digits = 4), "\n")
  cat("names:", names(colMeans(m))[seq_len(min(5, ncol(m)))], "\n")
  e <- tryCatch(posterior_epred(s), error = function(e) conditionMessage(e))
  cat("posterior_epred dim:", if (is.character(e)) e else dim(e), "\n")
  p <- tryCatch(posterior_predict(s), error = function(e) conditionMessage(e))
  cat("posterior_predict dim:", if (is.character(p)) p else dim(p), "\n")
  l <- tryCatch(log_lik(s), error = function(e) conditionMessage(e))
  cat("log_lik dim:", if (is.character(l)) l else dim(l), "\n")
  print(prior_summary(s))
}
pb <- brms::default_prior(brms::bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), data = d,
                          family = brms::cumulative())
print(pb[, c("prior", "class", "coef", "nlpar")])
print(default_prior(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), data = d,
                    family = cumulative()))

d2 <- E$eq_dat(seed = 158)
d2$y <- 2 * exp(0.3 * d2$x) + rnorm(nrow(d2), 0, 0.5)
M <- diag(nrow(d2))
s2 <- run("D7 nl fcor", frm_sample(
  bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE, autocor = ~ fcor(M)),
  data = d2, data2 = list(M = M), family = gaussian(), chains = 1, iter = 500,
  warmup = 250, seed = 158, refresh = 0,
  prior = c(set_prior("normal(1, 2)", nlpar = "b1"),
            set_prior("normal(0, 1)", nlpar = "b2"))))
if (!is.null(s2)) {
  cat("variables:", paste(head(variables(s2), 12), collapse = ", "), "\n")
  print(colMeans(as.matrix(s2))[1:3])
}
s3 <- run("D7 nl ar(cov = TRUE)", frm_sample(
  bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE,
     autocor = ~ ar(time = t, gr = g, cov = TRUE)),
  data = d2, family = gaussian(), chains = 1, iter = 500,
  warmup = 250, seed = 158, refresh = 0,
  prior = c(set_prior("normal(1, 2)", nlpar = "b1"),
            set_prior("normal(0, 1)", nlpar = "b2"))))
if (!is.null(s3)) {
  cat("variables:", paste(head(variables(s3), 12), collapse = ", "), "\n")
  print(colMeans(as.matrix(s3))[1:4])
  p <- tryCatch(posterior_predict(s3), error = function(e) conditionMessage(e))
  cat("posterior_predict dim:", if (is.character(p)) p else dim(p), "\n")
}
s4 <- run("multi-mu nl categorical", {
  set.seed(5001)
  dc <- data.frame(x = rnorm(150))
  dc$cat <- factor(sample(c("a", "b", "c"), 150, TRUE))
  frm_sample(bf(cat ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), data = dc,
             family = categorical(), chains = 1, iter = 400, warmup = 200,
             seed = 1, refresh = 0,
             prior = c(set_prior("normal(0, 2)", nlpar = "b0"),
                       set_prior("normal(0, 2)", nlpar = "b1")))
})
if (!is.null(s4)) {
  cat("variables:", paste(head(variables(s4), 12), collapse = ", "), "\n")
  e <- tryCatch(posterior_epred(s4), error = function(e) conditionMessage(e))
  cat("posterior_epred dim:", if (is.character(e)) e else dim(e), "\n")
}
