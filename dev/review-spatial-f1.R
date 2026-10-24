# Reviewer, short check of fix round 1: the lpdf heuristic on correct
# structures, reserved names, replace =, compat rows, wrappers, gr(g).
.libPaths(c("C:/Users/adf44/source/r/wt-spatial-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
try1 <- function(label, expr) {
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  [warning]", substr(conditionMessage(w), 1, 160), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) paste("ERROR:", substr(gsub("\n", " ",
                                                      conditionMessage(e)), 1, 260)))
  if (is.character(r) && length(r) == 1L && startsWith(r, "ERROR")) {
    cat(label, "->", r, "\n")
  } else if (inherits(r, "frmtmb_fit")) {
    cat(label, "-> FIT logLik", format(as.numeric(logLik(r)), digits = 10),
        "\n")
  } else cat(label, "-> ok\n")
  invisible(r)
}
local_frmtmb <- if (exists("frmtmb_local_registries")) "exported" else "absent"
cat("frmtmb_local_registries:", local_frmtmb, "\n")

# 1. An improper RW1 field (rank n - 1), the local-level model
#    y_t = b_t + e_t with no intercept: a correct model whose density is
#    flat along the constant vector, so lpdf(1000 * 1) == lpdf(0).
n <- 30
D <- diff(diag(n))
Q1 <- crossprod(D)
dimnames(Q1) <- list(as.character(1:n), as.character(1:n))
rw_entry <- function(rank_def) list(
  lpdf = function(b, theta, blk) {
    Q <- blk[["aux"]][["Q"]]
    k <- nrow(Q) - rank_def
    -k * theta[1] - 0.5 * exp(-2 * theta[1]) * sum(b * (Q %*% b)) -
      0.5 * k * log(2 * pi)
  },
  vcov = function(theta, blk) matrix(exp(2 * theta[1]), 1, 1,
                                     dimnames = list("sd", "sd")),
  npar = 1L, sd_idx = 1L,
  natural = function(theta, blk) list(term = "sd(rw)", type = "sd",
                                      est = theta[1], grad = matrix(1, 1, 1)),
  levels_coupled = TRUE,
  term = list(args = c("Q", "gr"), usage = "rw(Q, gr = t)",
              parse = function(args, env) list(Q = args[["Q"]]),
              build = function(spec, gr, lookup) {
                Q <- lookup(spec[["Q"]], "rw()")
                lv <- rownames(Q)
                list(levels = lv, index = match(as.character(gr), lv),
                     aux = list(Q = Q))
              }))
e1 <- rw_entry(1L)
try1("register rw1 (improper)", do.call(frmtmb_register_covstruct,
                                        c(list("rwone"), e1)))
set.seed(7)
d <- data.frame(t = 1:n, x = rnorm(n))
d$y <- 0.5 * d$x + cumsum(rnorm(n, 0, 0.7)) + rnorm(n, 0, 0.5)
try1("fit y ~ 0 + x + rwone(Q1, gr = t), improper RW1",
     frm(y ~ 0 + x + rwone(Q1, gr = t), data = d, data2 = list(Q1 = Q1)))
# the same model through core's gr(prec = Q)? (a singular precision)
try1("same model, core (1 | gr(t, prec = Q1))",
     frm(y ~ 0 + x + (1 | gr(t, prec = Q1)), data = d,
         data2 = list(Q1 = Q1)))
# the proper version: Q + 1e-6 I, as a control
Q1p <- Q1 + diag(1e-6, n)
dimnames(Q1p) <- dimnames(Q1)
try1("fit proper RW1 (Q + 1e-6 I) through the seam",
     frm(y ~ 0 + x + rwone(Q1p, gr = t), data = d, data2 = list(Q1 = Q1p)))
# what the lane round accepted: the improper fit, by direct evaluation
# 2. A bar structure with a nonzero fixed mean (contrived, valid density)
for (m in c(100, 600)) {
  nm <- paste0("shift", m)
  try1(paste("register", nm), frmtmb_register_covstruct(
    nm,
    lpdf = local({mm <- m; function(b, theta, blk) {
      sum(RTMB::dnorm(b, mm, exp(theta[1]), log = TRUE))}}),
    vcov = function(theta, blk) {
      matrix(exp(2 * theta[1]), 1, 1, dimnames = list(blk$cnms, blk$cnms))
    }, npar = 1L, sd_idx = 1L))
}
set.seed(1)
dg <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
dg$y <- dg$x + rnorm(20)[dg$g] + rnorm(100)
try1("fit shift100(1 | g)", frm(y ~ x + shift100(1 | g), data = dg))
try1("fit shift600(1 | g)", frm(y ~ x + shift600(1 | g), data = dg))
# 3. A heavy-tailed density (Cauchy), a correct structure
try1("register cauchy", frmtmb_register_covstruct(
  "cauchyre",
  lpdf = function(b, theta, blk) {
    s <- exp(theta[1])
    sum(-log(pi) - log(s) - log1p((b / s)^2))
  },
  vcov = function(theta, blk) matrix(exp(2 * theta[1]), 1, 1,
                                     dimnames = list(blk$cnms, blk$cnms)),
  npar = 1L, sd_idx = 1L))
try1("fit cauchyre(1 | g)", frm(y ~ x + cauchyre(1 | g), data = dg))
# 4. Reserved names: plausible structure names
lp <- function(b, theta, blk) sum(RTMB::dnorm(b, 0, exp(theta[1]), log = TRUE))
vc <- function(theta, blk) matrix(exp(2 * theta[1]), 1, 1,
                                  dimnames = list(blk$cnms, blk$cnms))
for (nm in c("toeplitz", "spline", "kernel", "cov", "cor", "dist", "filter",
             "matern", "rw1", "rw2", "iid", "besag", "bym", "ar2",
             "exch", "kronecker", "nested", "circulant", "seasonal",
             "exponential", "gaussian", "spherical", "power", "sparse")) {
  r <- tryCatch({frmtmb_register_covstruct(nm, lpdf = lp, vcov = vc,
                                           npar = 1L, sd_idx = 1L); "ok"},
                error = function(e) paste("refused:", substr(conditionMessage(e), 1, 110)))
  cat(sprintf("  name %-11s %s\n", nm, r))
}
# 5. replace =
f5 <- function() frmtmb_register_covstruct("rw1", lpdf = lp, vcov = vc,
                                           npar = 1L, sd_idx = 1L)
try1("re-register rw1 from the same (global) caller", f5())
ns <- new.env(); environment(f5) <- ns
try1("re-register rw1 from another env", f5())
try1("re-register with replace = TRUE", frmtmb_register_covstruct(
  "rw1", lpdf = lp, vcov = vc, npar = 1L, sd_idx = 1L, replace = TRUE))
# 6. compat rows
print(as.data.frame(frm_compat("cauchyre"))[
  frm_compat("cauchyre")$feature_b %in% c("quadrature", "importance",
                                          "gr_by", "|ID|", "mm()"),
  c("feature_a", "feature_b", "status")])
# 7. wrappers
try1("foo(bar(1 | g))", frm(y ~ x + foo(bar(1 | g)), data = dg))
f0 <- frm(y ~ x + (1 | g), data = dg)
try1("predict re_formula = ~ foo(1 | g)", predict(f0, re_formula = ~ foo(1 | g)))
try1("predict re_formula = ~ (1 | g)", predict(f0, re_formula = ~ (1 | g)))
try1("predict re_formula = NA", predict(f0, re_formula = NA))
# 8. B1
fa <- frm(y ~ x + (1 | gr(g)), data = dg)
cat("B1: (1 | gr(g)) logLik", format(as.numeric(logLik(fa)), digits = 12),
    " (1 | g)", format(as.numeric(logLik(f0)), digits = 12),
    " variables identical:", identical(variables(fa), variables(f0)), "\n")
try1("(1 + x | gr(g))", frm(y ~ x + (1 + x | gr(g)), data = dg))
try1("diag(1 | gr(g))", frm(y ~ x + diag(1 | gr(g)), data = dg))
try1("(1 | gr(g, by = NULL))", frm(y ~ x + (1 | gr(g, by = NULL)), data = dg))
