# Reviewer: fit_fd_se()'s chain-rule route against the per-coefficient
# route on models where the chain rule's premise (an output row depends
# on the block's effects only through its own eta, linearly) is under
# stress: a gp() inside an nlpar, a gp() in disc (log link), mo() and
# me() beside a gp(), a categorical family's per-category predictor, a
# new grouping level. Both routes run on the lane build in one process;
# fitted()'s public Est.Error is printed for the arm given.
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
set.seed(11)
n <- 120
d <- data.frame(x = round(runif(n, 0, 5), 1), z = rnorm(n),
                g = factor(sample(letters[1:8], n, TRUE)),
                m = factor(sample(1:4, n, TRUE), ordered = TRUE))
d$sdx <- 0.2
d$xo <- d$z + rnorm(n, 0, 0.2)
lat <- sin(d$x) + 0.5 * d$z + rnorm(8)[as.integer(d$g)] * 0.5 +
  0.3 * as.integer(d$m) + rlogis(n)
d$y <- cut(lat, quantile(lat, c(0, 0.3, 0.6, 1)), include.lowest = TRUE,
           labels = FALSE)
d$y <- factor(d$y, ordered = TRUE)
d$yc <- factor(sample(c("A", "B", "C"), n, TRUE))
nd_new <- data.frame(x = c(5.6, 6.3, 2.55), z = c(0, 1, -1),
                     g = factor(c("a", "b", "c"), levels = levels(d$g)),
                     m = factor(c(1, 2, 3), levels = 1:4, ordered = TRUE),
                     xo = c(0, 1, -1), sdx = 0.2)
nd_lev <- nd_new
nd_lev$g <- factor(c("zz", "zz", "yy"))

route <- function(fit, newdata, re_formula = NULL, allow_new_levels = FALSE,
                  chain = TRUE, use_extra = TRUE) {
  est <- fitted(fit, newdata = newdata, re_formula = re_formula,
                allow_new_levels = allow_new_levels)
  f <- function(fx) {
    ns$fitted_point(fx, newdata, re_formula, "response", NULL, NULL,
                    allow_new_levels)
  }
  want <- ns$smooth_b_idx(fit)
  if (ns$re_form_keeps(re_formula)) {
    want <- sort(unique(c(want, ns$re_governed_b(fit))))
  }
  b_idx <- if (!length(want)) NULL else {
    used <- ns$re_used_b(fit, newdata, NULL, allow_new_levels)
    if (is.null(used)) want else intersect(want, used)
  }
  bt <- if (length(b_idx)) {
    ns$re_b_batches(fit, newdata, NULL, allow_new_levels, b_idx)
  }
  ch <- if (chain && exists("fd_eta_chain", ns)) {
    ns$fd_eta_chain(fit, newdata, NULL, allow_new_levels,
                    ns$re_form_keeps(re_formula), b_idx)
  } else list(chain = NULL, extra = NULL)
  n_eval <- 0L
  f2 <- function(fx) { n_eval <<- n_eval + 1L; f(fx) }
  args <- list(fit, f2, b_idx = b_idx, b_batch = bt)
  if ("b_chain" %in% names(formals(ns$fit_fd_se))) {
    args$b_chain <- ch$chain
    if (use_extra) args$extra <- ch$extra
  }
  se <- do.call(ns$fit_fd_se, args)
  list(se = se, n = n_eval, nchain = length(ch$chain),
       extra = length(ch$extra), pub = est)
}
cases <- list(
  list("cumulative gp(x) exact, new pos", y ~ gp(x), "cumulative", nd_new,
       NULL, FALSE),
  list("cumulative gp(x) exact, in sample", y ~ gp(x), "cumulative", NULL,
       NULL, FALSE),
  list("cumulative mo(m) + gp(x, k = 8)", y ~ mo(m) + gp(x, k = 8),
       "cumulative", nd_new, NULL, FALSE),
  list("cumulative me(xo, sdx) + gp(x)", y ~ me(xo, sdx) + gp(x),
       "cumulative", nd_new, NULL, FALSE),
  list("cumulative gp(x) + (1|g) new level", y ~ gp(x) + (1 | g),
       "cumulative", nd_lev, NULL, TRUE),
  list("cumulative s(x) + (1|g) in sample", y ~ s(x, k = 6) + (1 | g),
       "cumulative", NULL, NULL, FALSE),
  list("categorical muB ~ gp(x, k = 6)",
       bf(yc ~ 1, muB ~ gp(x, k = 6), muC ~ z), "categorical", nd_new,
       NULL, FALSE),
  list("cumulative disc ~ gp(x, k = 6)",
       bf(y ~ z, disc ~ 0 + gp(x, k = 6)), "cumulative", nd_new, NULL,
       FALSE),
  list("cumulative nl: a ~ gp(x, k = 6)",
       bf(y ~ a * z + b, a ~ gp(x, k = 6), b ~ 0 + x, nl = TRUE),
       "cumulative", nd_new, NULL, FALSE),
  list("cumulative nl exp(a): a ~ gp(x)",
       bf(y ~ exp(a) * z, a ~ gp(x), nl = TRUE), "cumulative", nd_new, NULL,
       FALSE)
)
for (cs in cases) {
  lab <- cs[[1]]
  fam <- switch(cs[[3]], cumulative = cumulative(),
                categorical = categorical())
  fo <- if (inherits(cs[[2]], "formula")) bf(cs[[2]]) else cs[[2]]
  fit <- tryCatch(suppressWarnings(frm(fo, data = d, family = fam)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("CASE %-38s | FIT REFUSED: %s\n", lab,
                substr(conditionMessage(fit), 1, 160)))
    next
  }
  r1 <- tryCatch(route(fit, cs[[4]], cs[[5]], cs[[6]], chain = TRUE),
                 error = function(e) e)
  r0 <- tryCatch(route(fit, cs[[4]], cs[[5]], cs[[6]], chain = FALSE),
                 error = function(e) e)
  if (inherits(r1, "error") || inherits(r0, "error")) {
    cat(sprintf("CASE %-38s | ROUTE ERROR %s %s\n", lab,
                if (inherits(r1, "error")) conditionMessage(r1) else "",
                if (inherits(r0, "error")) conditionMessage(r0) else ""))
    next
  }
  pub <- r1$pub[, "Est.Error", , drop = TRUE]
  if (arm == "lane") {
    # the chain route without its extra term against the per-coefficient
    # route: the coefficient part must agree
    rx <- route(fit, cs[[4]], cs[[5]], cs[[6]], chain = TRUE, use_extra = FALSE)
    ratio <- rx$se / r0$se; rtot <- r1$se / r0$se
    cat(sprintf(paste0("CASE %-38s | chain blocks %d extra %d | evals ",
                       "chain %d vs per-coef %d | se chain/per-coef ratio ",
                       "range [%.10f, %.10f] | public == chain %s\n"),
                lab, r1$nchain, r1$extra, r1$n, r0$n, min(ratio),
                max(ratio), min(rtot), max(rtot), isTRUE(all.equal(unname(as.vector(pub)),
                                             unname(as.vector(r1$se)),
                                             tolerance = 1e-12))))
  } else {
    cat(sprintf("CASE %-38s | base public Est.Error %s\n", lab,
                paste(sprintf("%.10f", as.vector(pub)[1:min(6, length(pub))]),
                      collapse = " ")))
  }
  if (arm == "lane") {
    cat(sprintf("     lane public Est.Error %s\n",
                paste(sprintf("%.10f", as.vector(pub)[1:min(6, length(pub))]),
                      collapse = " ")))
  }
}
cat("DONE\n")
