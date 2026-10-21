# Reviewer, fix round 2: edges of prior_cellmeans_verdict(). Data of
# dev/fuzzfix-rev-r1-extra.R (seed 20261021, n = 80, sum(w) - n = 19.80).
#   Rscript dev/fuzzfix-rev-r2-edges.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a
run <- function(name, expr) {
  ws <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  if (inherits(r, "error")) {
    m <- conditionMessage(r)
    tag <- if (grepl("has no maximum|is improper", m)) "REFUSED" else "ERROR"
    cat(sprintf("ARM %-36s %s | %s\n", name, tag, substr(m, 1, 120)))
    invisible(m)
  } else {
    cat(sprintf("ARM %-36s ok logLik %.6g | betad %s | conv %s | warnings %d%s\n",
                name, as.numeric(logLik(r)),
                paste(sprintf("%.4g", r$estimates$betad), collapse = " "),
                r$opt$convergence %||% NA, length(ws),
                if (length(ws)) paste0(" (", substr(ws[1], 1, 50), ")") else ""))
    invisible(r)
  }
}
set.seed(20261021)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                fac = factor(rep(c("a", "b"), 40)),
                w = runif(n, 0.5, 2), g = factor(rep(1:20, each = 4)))
eta <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)
d$y <- eta + rnorm(n, 0, 0.8)
d$yg <- rgamma(n, shape = 3, rate = 3 / exp(eta))
cat(sprintf("excess: all %.3f | fac a %.3f | fac b %.3f\n", sum(d$w) - n,
            sum(d$w[d$fac == "a"] - 1), sum(d$w[d$fac == "b"] - 1)))
G <- gaussian()
P <- function(...) set_prior(...)
fo <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 0 + fac)
run("cellmeans, no prior", frm(fo, family = G, data = d))
run("cellmeans, normal(0,1) on b", frm(fo, family = G, data = d,
    prior = P("normal(0, 1)", class = "b", dpar = "sigma")))
run("cellmeans, normal on faca only", frm(fo, family = G, data = d,
    prior = P("normal(0, 1)", class = "b", dpar = "sigma", coef = "faca")))
run("cellmeans, normal on facb only", frm(fo, family = G, data = d,
    prior = P("normal(0, 1)", class = "b", dpar = "sigma", coef = "facb")))
run("cellmeans, flat b prior", frm(fo, family = G, data = d,
    prior = P("", class = "b", dpar = "sigma")))
run("cellmeans, student_t on b", frm(fo, family = G, data = d,
    prior = P("student_t(3, 0, 2.5)", class = "b", dpar = "sigma")))
cat("--- sigma ~ 0 + fac + x\n")
fox <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 0 + fac + x)
run("fac + x, no prior", frm(fox, family = G, data = d))
run("fac + x, normal(0,1) on b", frm(fox, family = G, data = d,
    prior = P("normal(0, 1)", class = "b", dpar = "sigma")))
run("fac + x, normal on faca + facb", frm(fox, family = G, data = d,
    prior = c(P("normal(0, 1)", class = "b", dpar = "sigma", coef = "faca"),
              P("normal(0, 1)", class = "b", dpar = "sigma", coef = "facb"))))
cat("--- Gamma, shape ~ 0 + fac\n")
fg <- bf(yg | weights(w) ~ x + (1 | obs), shape ~ 0 + fac)
GA <- Gamma("log")
run("Gamma cellmeans, no prior", frm(fg, family = GA, data = d))
run("Gamma cellmeans, gamma(2,1) on b", frm(fg, family = GA, data = d,
    prior = P("gamma(2, 1)", class = "b", dpar = "shape")))
run("Gamma cellmeans, normal(0,1) on b", frm(fg, family = GA, data = d,
    prior = P("normal(0, 1)", class = "b", dpar = "shape")))
cat("--- the message, gaussian (1 | obs), no prior\n")
m <- tryCatch(frm(bf(y | weights(w) ~ x + (1 | obs)), family = G, data = d),
              error = function(e) conditionMessage(e))
cat(m, "\n")
cat("DONE\n")
