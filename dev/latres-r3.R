# R3: statistics for a unit-free M1 condition, on three sets of designs, each fit
# the one the exact MLE was measured on:
#   lane  dev/latres-m1.R,  seeds 1 to 40 (exact in dev/latres-out/m1)
#   m1b   dev/review-latres-m1b.R, seeds 3001 to 3040 (exact in its tsv)
#   m4b   dev/review-latres-m4b.R, seeds 2001 to 2040, cov = TRUE (exact in m4b2)
# Usage: Rscript dev/latres-r3.R <set> <design> <from> <to> <outfile>
args <- commandArgs(TRUE)
set <- args[[1]]; design <- args[[2]]
s1 <- as.integer(args[[3]]); s2 <- as.integer(args[[4]]); outf <- args[[5]]
LIB <- "C:/Users/adf44/source/r/wt-latres-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-latres")
suppressPackageStartupMessages(library(frmtmb))
options(warn = -1)
ns <- asNamespace("frmtmb")
`%||%` <- function(a, b) if (is.null(a)) b else a
assignInNamespace("latres_refuse_margin", function(fit) invisible(NULL), "frmtmb")
load_gen <- function(file, stop_at, argv) {
  l <- readLines(file, warn = FALSE)
  l <- l[seq_len(grep(stop_at, l)[1L] - 1L)]
  e <- new.env(parent = globalenv())
  e$commandArgs <- function(...) argv
  eval(parse(text = l), e)
  e
}
gen <- switch(set,
  lane = load_gen("dev/latres-m1.R", "^row_info <-|^# curvature of each row", c(LIB, design, s1, s2, "x")),
  m1b = load_gen("dev/review-latres-m1b.R", "^extra_par <-", c(design, s1, s2, "x")),
  m4b = load_gen("dev/review-latres-m4b.R", "^exact_T <-", c(design, s1, s2, "x")))
options(warn = -1)
fam <- switch(set, lane = gen$famf, gen$fam)
lhs <- switch(set,
  lane = if (gen$D$fam == "binomial") "y | trials(tr)" else "y",
  m1b = if (gen$D$dist %in% c("binom", "zib")) "y | trials(tr)" else "y",
  m4b = if (gen$D$dist %in% c("binom", "bb")) "y | trials(tr)" else "y")
form <- as.formula(paste(lhs, "~ x + ar(t, g, cov = TRUE)"))
derivs <- function(fit, bk) {
  r <- bk$aux$ac$resp
  f <- fit$spec$responses[[r]]$family
  dp <- ns$eval_dpars(fit)[[r]]
  lk <- fit$frame$linpreds[[paste0(r, ".mu")]]$link
  eta <- lk$linkfun(as.numeric(dp$mu))
  y <- fit$frame$y[[r]]
  av <- fit$frame$aterm_values[[r]]
  ex <- ns$fit_extras(fit)
  ll <- function(e) {
    dp$mu <- lk$linkinv(e)
    as.numeric(ns$row_lpdf(f, y, y, dp, av, ex))
  }
  h <- 1e-3 * max(1, stats::median(abs(eta)))
  l0 <- ll(eta); lp <- ll(eta + h); ln <- ll(eta - h)
  lpp <- ll(eta + 2 * h); lnn <- ll(eta - 2 * h)
  list(d2 = -(lp - 2 * l0 + ln) / h^2,
       d3 = (lpp - 2 * lp + 2 * ln - lnn) / (2 * h^3),
       d4 = (lpp - 4 * lp + 6 * l0 - 4 * ln + lnn) / h^4)
}
first <- TRUE
for (s in s1:s2) {
  d <- gen$simd(s)
  w <- character(0)
  f <- withCallingHandlers(
    tryCatch(frm(form, data = d, family = fam), error = function(e) e),
    warning = function(c) {
      if (inherits(c, "frmtmb_latres_laplace")) w <<- conditionMessage(c)
      invokeRestart("muffleWarning")
    },
    message = function(c) invokeRestart("muffleMessage"))
  o <- data.frame(set = set, design = design, seed = s, err = "", sderr = NA_real_,
                  ar = NA_real_, v = NA_real_, info_med = NA_real_, info_mean = NA_real_,
                  kappa_med = NA_real_, lambda = NA_real_, g3 = NA_real_, g4 = NA_real_,
                  disp = NA_real_, ldisp = NA_real_, se_ls = NA_real_,
                  warned = NA, reason = "")
  if (inherits(f, "error")) {
    o$err <- substr(conditionMessage(f), 1, 60)
  } else {
    bk <- ns$frame_latres(f$frame)[[1L]]
    th <- f$estimates$theta[bk$theta_idx]
    v <- ns$latres_marginal_sd(th, bk$aux)^2
    g <- tryCatch(derivs(f, bk), error = function(e) NULL)
    o$sderr <- exp(th[1]); o$ar <- th[2] / sqrt(1 + th[2]^2); o$v <- v
    o$warned <- length(w) > 0L
    o$reason <- if (length(w)) substr(sub(".*this response: ", "", w), 1, 24) else ""
    if (!is.null(g)) {
      i2 <- pmax(g$d2, 0)
      o$info_med <- median(i2); o$info_mean <- mean(i2[is.finite(i2)])
      o$kappa_med <- median(i2) * v
      o$lambda <- mean(i2 * v / (1 + i2 * v), na.rm = TRUE)
      H <- i2 + 1 / v
      o$g3 <- mean(abs(g$d3) / H^1.5, na.rm = TRUE)
      o$g4 <- mean(abs(g$d4) / H^2, na.rm = TRUE)
    }
    links <- f$spec$responses[[1L]]$family$links
    dd <- setdiff(names(links), "mu")
    dd <- dd[vapply(links[dd], function(l) identical(l$name, "log"), NA)]
    o$disp <- length(dd)
    par <- f$opt$par
    if (any(names(par) == "betad")) o$ldisp <- par[names(par) == "betad"][1]
    vc <- tryCatch(as.matrix(vcov(f, full = TRUE)), error = function(e) NULL)
    o$se_ls <- tryCatch(sqrt(vc["theta_1", "theta_1"]), error = function(e) NA)
  }
  write.table(o, outf, sep = "\t", row.names = FALSE, quote = FALSE,
              col.names = first, append = !first)
  first <- FALSE
}
