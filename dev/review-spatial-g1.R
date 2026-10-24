# Reviewer: car() and spde() spellings with and without frmtmb.spatial.
# Usage: Rscript review-spatial-g1.R <libs> <load spatial: 0|1>
a <- commandArgs(TRUE)
.libPaths(c(strsplit(a[1], ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
if (a[2] == "1") suppressMessages(library(frmtmb.spatial))
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), " spatial loaded:",
    isNamespaceLoaded("frmtmb.spatial"), "\n")
lattice_W <- function(r, c) {
  g <- expand.grid(r = seq_len(r), c = seq_len(c))
  n <- nrow(g)
  W <- 1 * (abs(outer(g$r, g$r, "-")) + abs(outer(g$c, g$c, "-")) == 1)
  dimnames(W) <- list(paste0("L", seq_len(n)), paste0("L", seq_len(n)))
  W
}
set.seed(3)
W <- lattice_W(4, 4)
d <- data.frame(loc = factor(rep(rownames(W), each = 4), levels = rownames(W)),
                x = rnorm(64), car = rnorm(64), node = rep(1:16, each = 4))
d$y <- d$x + rnorm(16)[as.integer(d$loc)] + rnorm(64, 0, 0.5)
d$y2 <- rnorm(64)
# a 1-D chain mesh as the fem triple
n <- 16
M0 <- diag(1, n); M1 <- diag(2, n); M1[cbind(1:(n - 1), 2:n)] <- -1
M1[cbind(2:n, 1:(n - 1))] <- -1; M1[1, 1] <- M1[n, n] <- 1
M2 <- M1 %*% M1
fem <- list(M0 = M0, M1 = M1, M2 = M2)
cases <- list(
  car_bare = quote(frm(y ~ x + car(W, gr = loc), data = d, data2 = list(W = W))),
  car_bf = quote(frm(bf(y ~ x + car(W, gr = loc, type = "icar")), data = d,
                     data2 = list(W = W))),
  car_bar = quote(frm(y ~ x + (1 | car(W, gr = loc)), data = d,
                      data2 = list(W = W))),
  car_sigma = quote(frm(bf(y ~ x, sigma ~ car(W, gr = loc)), data = d,
                        data2 = list(W = W))),
  car_nl = quote(frm(bf(y ~ a + b * x, a ~ 1 + car(W, gr = loc), b ~ 1,
                        nl = TRUE), data = d, data2 = list(W = W),
                     prior = c(set_prior("normal(0, 5)", nlpar = "a"),
                               set_prior("normal(0, 5)", nlpar = "b")))),
  car_mv = quote(frm(bf(y ~ x + car(W, gr = loc)) + bf(y2 ~ x) +
                       set_rescor(FALSE), data = d, data2 = list(W = W))),
  car_ns = quote(frm(y ~ x + frmtmb.spatial::car(W, gr = loc), data = d,
                     data2 = list(W = W))),
  car_column = quote(frm(y ~ x + car, data = d)),
  car_brms_ns = quote(frm(y ~ x + brms::car(W, gr = loc), data = d,
                          data2 = list(W = W))),
  spde_bare = quote(frm(y ~ x + spde(fem, gr = node), data = d,
                        data2 = list(fem = fem))),
  spde_bar = quote(frm(y ~ x + (1 | spde(fem, gr = node)), data = d,
                       data2 = list(fem = fem))),
  car_dryrun = quote(frm(bf(y ~ x + car(M, gr = g)), data = NULL,
                         dry_run = "spec")),
  car_inside_I = quote(frm(y ~ x + I(car(W, gr = loc)), data = d,
                           data2 = list(W = W))),
  car_reformula = quote({
    f <- frm(y ~ x + (1 | loc), data = d)
    predict(f, re_formula = ~ car(W, gr = loc))
  }))
for (nm in names(cases)) {
  r <- tryCatch({
    v <- suppressWarnings(suppressMessages(eval(cases[[nm]])))
    if (inherits(v, "frmtmb_fit")) {
      paste("FIT", paste(vapply(v$frame$re_blocks, `[[`, "", "covstruct"),
                         collapse = ","), "logLik",
            format(as.numeric(logLik(v)), digits = 10))
    } else "OK"
  }, error = function(e) paste("ERR", substr(gsub("\n", " ",
                                                  conditionMessage(e)), 1, 220)))
  cat(sprintf("%-14s -> %s\n", nm, r))
}
