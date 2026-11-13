# Review of lane sampfront: attack the subset-newdata route of log_lik().
#   Rscript dev/review-sampfront-subset.R [base]
# Draws are built from frm() fits as in dev/sampfront-repro.R's lap_pair()
# (full draws, no Stan). For each case: REFUSED (with the reason's tail),
# or the columns returned, compared with the fitted column of the row the
# user meant ("meant" is the row of x$data the subset came from).
args <- commandArgs(TRUE)
base <- length(args) && identical(args[1], "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb); library(frmtmb.sample)
})
cat("arm", if (base) "base" else "lane", find.package("frmtmb.sample"), "\n")
full_draws <- function(fit, n = 20L, seed = 1L) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit),
            class = "frmtmb_draws")
}
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
# meant: the x$data rows the subset is meant to be; NULL = should refuse
try_sub <- function(lab, x, nd, meant, fun = log_lik) {
  cat("\n-- ", lab, "\n", sep = "")
  ll <- tryCatch(suppressWarnings(fun(x, newdata = NULL)),
                 error = function(e) NULL)
  r <- tryCatch(suppressWarnings(fun(x, newdata = nd)), error = function(e) e)
  if (inherits(r, "error")) {
    m <- gsub("\\s+", " ", conditionMessage(r))
    tail <- if (grepl("This newdata", m)) sub(".*(This newdata)", "\\1", m) else
      substr(m, 1, 160)
    cat("  REFUSED:", tail, "\n")
    if (!is.null(meant)) cat("  (a correct subset was refused)\n")
    return(invisible())
  }
  if (is.null(meant)) {
    cat("  ACCEPTED, dim", dim(r), " (expected a refusal)\n")
    return(invisible(r))
  }
  d <- max(abs(unclass(r)[, , drop = FALSE] -
                 unclass(ll)[, meant, drop = FALSE]))
  cat("  ACCEPTED, dim", dim(r), "; max |col - fitted col of meant row| =",
      format(d), if (d > 0) " SILENT WRONG" else "", "\n")
  invisible(r)
}

set.seed(11)
d <- data.frame(g = factor(rep(1:6, each = 8)), x = rnorm(48),
                t = rep(1:8, 6))
d$y <- 1 + 0.5 * d$x + rnorm(6)[d$g] * 0.7 + rnorm(48)
d$w <- runif(48, 0.5, 2)
d$o <- log(runif(48, 1, 3))
d$cnt <- rpois(48, exp(0.3 + 0.2 * d$x))
d$hi <- d$y + runif(48, 0.1, 1)
d$cc <- ifelse(seq_len(48) %% 5 == 0, "interval", "none")

f0 <- q(bf(y ~ x + (1 | g)), data = d)
x0 <- full_draws(f0)
try_sub("control: d[1:10, ]", x0, d[1:10, ], 1:10)
try_sub("reordered d[c(20, 5, 3), ]", x0, d[c(20, 5, 3), ], c(20, 5, 3))
try_sub("duplicated d[c(3, 3, 3), ]", x0, d[c(3, 3, 3), ], c(3, 3, 3))
try_sub("fit$data[4:6, ] (the frame)", x0, x0$fit$data[4:6, ], 4:6)
dr <- d[11:20, ]; rownames(dr) <- NULL
try_sub("names reset, d[11:20, ] named 1:10", x0, dr, NULL)
# rows that coincide in value: copy rows 1:3 to 41:43 in the data, refit
dco <- d; dco[41:43, ] <- dco[1:3, ]
fco <- q(bf(y ~ x + (1 | g)), data = dco)
xco <- full_draws(fco)
sub <- dco[41:43, ]; rownames(sub) <- c("1", "2", "3")
try_sub("rows 41:43 equal to rows 1:3 in every variable (g too), named 1:3",
        xco, sub, 41:43)
dco2 <- d; dco2[41:43, ] <- dco2[1:3, ]; dco2$g[41:43] <- dco2$g[1:3]
fco2 <- q(bf(y ~ x + (1 | g)), data = dco2)
xco2 <- full_draws(fco2)
sub2 <- dco2[41:43, ]; rownames(sub2) <- c("1", "2", "3")
try_sub("rows 41:43 equal in every variable incl. g to rows 1:3, named 1:3",
        xco2, sub2, 41:43)
dv <- d[1:5, ]; dv$x[2] <- dv$x[2] + 1e-12
try_sub("value changed by 1e-12, name kept", x0, dv, NULL)
dy <- d[1:5, ]; dy$y[2] <- dy$y[2] + 1
try_sub("response changed, name kept", x0, dy, NULL)
dl <- droplevels(d[1:10, ])
try_sub("droplevels() subset (g has 2 of 6 levels)", x0, dl, 1:10)
dfl <- d[1:10, ]; dfl$g <- factor(as.character(dfl$g), levels = c(2, 1, 3:6))
try_sub("factor with reordered levels", x0, dfl, 1:10)
dch <- d[1:10, ]; dch$g <- as.character(dch$g)
try_sub("g as character", x0, dch, 1:10)
dnum <- d[1:10, ]; dnum$g <- as.integer(as.character(dnum$g))
try_sub("g as integer", x0, dnum, 1:10)
dextra <- d[1:10, ]; dextra$junk <- 1
try_sub("extra column", x0, dextra, 1:10)
if (requireNamespace("tibble", quietly = TRUE)) {
  tb <- tibble::as_tibble(d)
  try_sub("tibble[1:10, ]", x0, tb[1:10, ], 1:10)
  try_sub("tibble[11:20, ] (names 1:10)", x0, tb[11:20, ], 11:20)
}
# custom row names in the data
dn <- d; rownames(dn) <- paste0("s", 1:48)
fn <- q(bf(y ~ x + (1 | g)), data = dn)
xn <- full_draws(fn)
cat("\n  rownames(fit$data)[1:3]:", rownames(xn$fit$data)[1:3], "\n")
try_sub("custom row names s1..s48, dn[5:8, ]", xn, dn[5:8, ], 5:8)
# a variable with NA: na.omit drops row 3
dna <- d; dna$x[3] <- NA
fna <- q(bf(y ~ x + (1 | g)), data = dna)
xna <- full_draws(fna)
cat("\n  rownames(fit$data)[1:4]:", rownames(xna$fit$data)[1:4], "\n")
try_sub("NA dropped row 3; dna[c(1, 2, 4, 5), ] (fitted rows 1,2,3,4)", xna,
        dna[c(1, 2, 4, 5), ], 1:4)
try_sub("NA dropped row 3; dna[1:5, ] includes the dropped row", xna,
        dna[1:5, ], NULL)
# row names "3.1" in the fitted data itself, from rbind of duplicates
drb <- rbind(d[1:20, ], d[3, ])
cat("\n  rbind names:", tail(rownames(drb), 2), "\n")
frb <- q(bf(y ~ x + (1 | g)), data = drb)
xrb <- full_draws(frb)
try_sub("data has row '31' (copy of 3); drb[c('3', '31'), ]", xrb,
        drb[c("3", "31"), ], c(3, 21))

# addition terms and other model pieces
fw <- q(bf(y | weights(w) ~ x + (1 | g)), data = d)
try_sub("weights(w): d[1:10, ]", full_draws(fw), d[1:10, ], 1:10)
dw <- d[1:10, ]; dw$w[2] <- dw$w[2] * 2
try_sub("weights(w): weight moved", full_draws(fw), dw, NULL)
fo <- q(bf(cnt ~ x + offset(o) + (1 | g)), family = poisson(), data = d)
try_sub("offset(o): d[1:10, ]", full_draws(fo), d[1:10, ], 1:10)
do <- d[1:10, ]; do$o[2] <- do$o[2] + 1
try_sub("offset(o): offset moved", full_draws(fo), do, NULL)
fce <- q(bf(y | cens(cc, hi) ~ x), data = d)
try_sub("cens(cc, hi): d[1:10, ]", full_draws(fce), d[1:10, ], 1:10)
dce <- d[1:10, ]; dce$hi[5] <- dce$hi[5] + 1
try_sub("cens(cc, hi): y2 moved on an interval row", full_draws(fce), dce,
        NULL)
ftr <- q(bf(y | trunc(lb = -10) ~ x), data = d)
try_sub("trunc(lb = -10): d[1:10, ]", full_draws(ftr), d[1:10, ], 1:10)
fpo <- q(bf(y ~ poly(x, 2) + (1 | g)), data = d)
try_sub("poly(x, 2): d[1:10, ]", full_draws(fpo), d[1:10, ], 1:10)
fsc <- q(bf(y ~ scale(x) + (1 | g)), data = d)
try_sub("scale(x): d[1:10, ]", full_draws(fsc), d[1:10, ], 1:10)

# structures where a row's column depends on other rows
for (sp in list(
  list("ar(time = t, gr = g)", bf(y ~ x + ar(time = t, gr = g))),
  list("ar(cov = TRUE)", bf(y ~ x + ar(time = t, gr = g, cov = TRUE))),
  list("arma", bf(y ~ x + arma(time = t, gr = g))),
  list("cosy(gr = g)", bf(y ~ x + cosy(gr = g))),
  list("unstr(time = t, gr = g)", bf(y ~ x + unstr(time = t, gr = g))))) {
  fs <- tryCatch(q(sp[[2]], data = d), error = function(e) e)
  if (inherits(fs, "error")) {
    cat("\n-- ", sp[[1]], ": fit failed: ", conditionMessage(fs), "\n")
    next
  }
  try_sub(paste0(sp[[1]], ": d[1:10, ]"), full_draws(fs), d[1:10, ], NULL)
}
M <- diag(48) * 0.5 + 0.5
fco_ <- tryCatch(q(bf(y ~ x + fcor(M)), data = d, data2 = list(M = M)),
                 error = function(e) e)
if (!inherits(fco_, "error")) {
  try_sub("fcor(M): d[1:10, ]", full_draws(fco_), d[1:10, ], NULL)
} else cat("\n-- fcor fit failed:", conditionMessage(fco_), "\n")
fgp <- tryCatch(q(bf(y ~ gp(x)), data = d), error = function(e) e)
if (!inherits(fgp, "error")) {
  try_sub("gp(x): d[1:10, ]", full_draws(fgp), d[1:10, ], 1:10)
} else cat("\n-- gp fit failed:", conditionMessage(fgp), "\n")
fs_ <- tryCatch(q(bf(y ~ s(x)), data = d), error = function(e) e)
if (!inherits(fs_, "error")) {
  try_sub("s(x): d[1:10, ]", full_draws(fs_), d[1:10, ], 1:10)
}
d$sdx <- 0.2
fme <- tryCatch(q(bf(y ~ me(x, sdx)), data = d), error = function(e) e)
if (!inherits(fme, "error")) {
  try_sub("me(x, sdx): d[1:10, ]", full_draws(fme), d[1:10, ], NULL)
}

# multivariate
d$z <- -1 + 0.3 * d$x + rnorm(48)
fmv <- q(mvbf(bf(y ~ x), bf(z ~ x), rescor = FALSE), data = d)
try_sub("mvbf no rescor: d[1:10, ]", full_draws(fmv), d[1:10, ], 1:10)
fmr <- q(mvbf(bf(y ~ x), bf(z ~ x), rescor = TRUE), data = d)
try_sub("mvbf rescor: d[1:10, ]", full_draws(fmr), d[1:10, ], 1:10)
d$sel <- rep(c(TRUE, FALSE), 24)
fsub <- tryCatch(q(mvbf(bf(y ~ x), bf(z | subset(sel) ~ x), rescor = FALSE),
                   data = d), error = function(e) e)
if (!inherits(fsub, "error")) {
  try_sub("mvbf subset(): d[1:10, ]", full_draws(fsub), d[1:10, ], NULL)
  try_sub("mvbf subset(): resp = 'y', d[1:10, ]", full_draws(fsub),
          d[1:10, ], NULL,
          fun = function(x, newdata) log_lik(x, newdata = newdata, resp = "y"))
}
# variables read only by a distributional or nonlinear formula, trials()
d$zz <- rnorm(48)
fsg <- q(bf(y ~ x, sigma ~ zz), data = d)
dz <- d[1:10, ]; dz$zz[3] <- dz$zz[3] + 1
try_sub("sigma ~ zz: zz moved", full_draws(fsg), dz, NULL)
try_sub("sigma ~ zz: d[1:10, ]", full_draws(fsg), d[1:10, ], 1:10)
d$nt <- 10L
d$k <- rbinom(48, 10, plogis(0.2 * d$x))
fbi <- q(bf(k | trials(nt) ~ x), family = binomial(), data = d)
dt <- d[1:10, ]; dt$nt[3] <- 12L
try_sub("trials(nt): nt moved", full_draws(fbi), dt, NULL)
fnl <- tryCatch(q(bf(y ~ a + b * zz, a ~ 1, b ~ 1, nl = TRUE), data = d,
                  prior = c(set_prior("normal(0, 2)", nlpar = "a"),
                            set_prior("normal(0, 2)", nlpar = "b"))),
                error = function(e) e)
if (!inherits(fnl, "error")) {
  try_sub("nl: zz moved", full_draws(fnl), dz, NULL)
}
# loo() and psis() through the same route
try_sub("loo(): names reset", x0, dr, NULL,
        fun = function(x, newdata) loo(x, newdata = newdata))
try_sub("psis(): names reset", x0, dr, NULL,
        fun = function(x, newdata) psis(x, newdata))
try_sub("loo_model_weights(): names reset", x0, dr, NULL,
        fun = function(x, newdata) loo_model_weights(x, x0, newdata = newdata))
cat("\nDONE\n")
