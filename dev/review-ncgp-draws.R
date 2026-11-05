# Item 4, draws side. Each labelled column of frm_sample()'s draws is matched to the
# stanfit column it is a monotone function of (Spearman 1): the map must read
# b_* <- beta, sigma <- betad, sd*/lscale <- theta, r_* <- b. Then each natural
# column against frm()'s estimate, in posterior sds. Chains are short: the map is
# deterministic on any chain that moves. Args: model names, or none for all.
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
suppressPackageStartupMessages(library(frmtmb.sample))
cat("lib sample:", dirname(system.file(package = "frmtmb.sample")), "\n")
set.seed(4242)
n <- 70
d <- data.frame(x = runif(n, 0, 6), z = rnorm(n), g = factor(sample(letters[1:7], n, TRUE)),
                o = ordered(sample(1:4, n, TRUE)), f = factor(sample(c("u", "v"), n, TRUE)))
re <- rnorm(7, sd = 0.6)
d$y <- 1.7 + sin(d$x) + 0.5 * d$z + re[as.integer(d$g)] + as.integer(d$o) * 0.3 + rnorm(n, sd = 0.25)
d$cnt <- rpois(n, exp(1 + 0.5 * sin(d$x) + 0.2 * d$z))
d$bin <- rbinom(n, 1, plogis(sin(d$x)))
d$y2 <- 0.5 + cos(d$x) + rnorm(n, sd = 0.3)
d$zm <- d$z; d$zm[sample(n, 6)] <- NA
d$zobs <- d$z + rnorm(n, sd = 0.1); d$sdz <- 0.1
d$ygev <- 2 + sin(d$x) + 0.4 * (((-log(runif(n)))^(-0.1) - 1) / 0.1)
models <- list(
  gp_re = list(f = y ~ z + gp(x) + (1 | g)),
  gp_s = list(f = y ~ gp(x) + s(z)),
  gp_mo = list(f = y ~ z + gp(x) + mo(o)),
  gp_me = list(f = y ~ me(zobs, sdz) + gp(x)),
  gp_mi = list(f = bf(y ~ mi(zm) + gp(x)) + bf(zm | mi() ~ 1) + set_rescor(FALSE)),
  mvbf = list(f = bf(y ~ gp(x)) + bf(y2 ~ z + gp(x)) + set_rescor(FALSE)),
  mvbf_res = list(f = bf(y ~ gp(x)) + bf(y2 ~ z) + set_rescor(TRUE)),
  sigma_gp = list(f = bf(y ~ z, sigma ~ gp(x))),
  mix = list(f = bf(y ~ gp(x)), family = mixture(gaussian(), gaussian())),
  pois = list(f = cnt ~ z + gp(x) + (1 | g), family = poisson()),
  bern = list(f = bin ~ gp(x), family = bernoulli()),
  ord = list(f = o ~ gp(x), family = cumulative()),
  by = list(f = y ~ gp(x, by = f)),
  student = list(f = y ~ gp(x), family = student()),
  gev = list(f = ygev ~ gp(x), family = gen_extreme_value()),
  noncent = list(f = y ~ gp(x) + (1 | g), sample = list(reparameterize = TRUE)),
  laplace = list(f = y ~ z + gp(x) + (1 | g), sample = list(laplace = TRUE)),
  reml_no = list(f = y ~ gp(x) + (z | g))
)
pick <- if (length(args)) args else names(models)
mapcols <- function(m, raw) {
  # m: draws matrix (chain-major). raw: iter x chain x par array of the stanfit
  r <- do.call(rbind, lapply(seq_len(dim(raw)[2]), function(ch) raw[, ch, ]))
  colnames(r) <- dimnames(raw)[[3]]
  ra <- apply(r, 2, rank, ties.method = "average")
  out <- character(ncol(m))
  names(out) <- colnames(m)
  for (j in seq_len(ncol(m))) {
    if (sd(m[, j]) == 0) { out[j] <- "const"; next }
    cj <- abs(as.vector(cor(rank(m[, j]), ra)))
    i <- which.max(cj)
    out[j] <- if (isTRUE(cj[i] > 0.999999)) colnames(r)[i] else sprintf("~%s(%.3f)", colnames(r)[i], cj[i])
  }
  out
}
for (mn in pick) {
  s <- models[[mn]]
  a <- s; a$f <- NULL; a$sample <- NULL
  t0 <- Sys.time()
  fit <- tryCatch(suppressWarnings(suppressMessages(do.call(frm, c(list(s$f, data = d), a)))), error = function(e) e)
  if (inherits(fit, "error")) { cat("RD", mn, "FIT ERROR", substr(conditionMessage(fit), 1, 150), "\n"); next }
  dr <- tryCatch(suppressWarnings(suppressMessages(do.call(frm_sample, c(list(fit, chains = 1, iter = 500, warmup = 250, seed = 11, refresh = 0,
                 control = list(max_treedepth = 8)), s$sample)))), error = function(e) e)
  if (inherits(dr, "error")) { cat("RD", mn, "SAMPLE ERROR", substr(conditionMessage(dr), 1, 200), "\n"); next }
  m <- dr$draws
  raw <- rstan::extract(dr$stanfit, permuted = FALSE)
  map <- mapcols(m, raw)
  keep <- !grepl("^r_|^b\\[|^zgp|^s_|^Ymi|^Xme|^lp__", names(map))
  cat(sprintf("RD %s rawnames %s ... labels %d sec %.0f\n", mn, paste(head(dimnames(raw)[[3]], 3), collapse = ","), length(map), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  for (k in names(map)[keep]) {
    cat(sprintf("RD   %-30s <- %-14s median %.4f sd %.4f\n", k, map[[k]], median(m[, k]), sd(m[, k])))
  }
  bcols <- names(map)[grepl("^b\\[", names(map))]
  ok <- sum(map[bcols] == bcols)
  cat(sprintf("RD   b[i] columns reading their own raw column: %d of %d; first wrong: %s\n", ok, length(bcols),
              paste(head(paste0(bcols[map[bcols] != bcols], "<-", map[bcols][map[bcols] != bcols]), 3), collapse = " ")))
  rcols <- names(map)[grepl("^r_|^s_", names(map))]
  if (length(rcols)) cat("RD   r_/s_ map:", paste(head(paste0(rcols, "<-", map[rcols]), 6), collapse = " "), "\n")
  # natural reading: frm()'s fixed effects and epred
  fe <- tryCatch(fixef(fit, flatten = TRUE), error = function(e) NULL)
  if (!is.null(fe)) for (nm in names(fe)) {
    col <- intersect(c(paste0("b_", sub("\\(Intercept\\)", "Intercept", nm))), colnames(m))
    if (length(col)) cat(sprintf("RD   z %-24s frm %.4f draws %.4f z %.2f\n", col, fe[[nm]], median(m[, col]), (median(m[, col]) - fe[[nm]]) / sd(m[, col])))
  }
  ep <- tryCatch(suppressWarnings(fitted(dr, summary = TRUE)), error = function(e) e)
  fp <- tryCatch(suppressWarnings(fitted(fit)), error = function(e) e)
  if (!inherits(ep, "error") && !inherits(fp, "error") && length(dim(ep)) == 2L && identical(dim(ep)[1], dim(fp)[1])) {
    z <- (ep[, 1] - fp[, 1]) / pmax(ep[, 2], 1e-8)
    cat(sprintf("RD   epred against frm fitted: max |z| %.2f median |z| %.2f\n", max(abs(z)), median(abs(z))))
  } else cat("RD   epred:", if (inherits(ep, "error")) substr(conditionMessage(ep), 1, 120) else "shape differs", "\n")
  # warmup, same columns
  w <- tryCatch(as.array(dr, inc_warmup = TRUE), error = function(e) e)
  cat("RD   warmup array:", if (inherits(w, "error")) substr(conditionMessage(w), 1, 100) else paste(dim(w), collapse = "x"), "\n")
  flush.console()
}
cat("DONE\n")
