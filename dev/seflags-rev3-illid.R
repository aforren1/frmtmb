# Reviewer of lane seflags, final check of B4. On "illcond" fits (raw
# degree-8 polynomial on [1, 2], with and without (1 | g)) plus two
# well-determined terms, z ~ N(0, 1) and a two-level factor w: which
# predictions keep an SE, and are they right? The reference is the same
# likelihood with poly(x, 8), which is well conditioned (the column
# space is the same, so b_z, b_wb, sigma, the sd and every fitted value
# are the same parameters). Also every row's share as jc_nonest() reads
# it, (g . d)^2 / lambda over g' V g, against the true missing share,
# reference variance over g' V g minus 1, on a grid of x.
#   Rscript dev/seflags-rev3-illid.R <lane|op3|m9|base> <out.tsv>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev3-common.R")
out <- commandArgs(trailingOnly = TRUE)[2]
rec <- list()
if (exists("se_ill_share", ns)) {
  txt <- deparse(ns$jc_nonest)
  i <- grep("bad <- bad | (vi > se_ill_share * vk) %in% TRUE", txt,
            fixed = TRUE)
  stopifnot(length(i) == 1L)
  txt <- append(txt, c("s_ <- fl$dscale[il$rows]; if (length(s_) != length(il$rows)) s_ <- rep(1, length(il$rows)); s_[!is.finite(s_) | !(s_ > 0)] <- 1", "Gu_ <- sweep(Gr, 2L, s_, \"/\"); du_ <- il$d * s_", "cs_ <- abs(Gu_ %*% du_) / outer(sqrt(rowSums(Gu_^2)), sqrt(colSums(du_^2)))", "rec[[length(rec) + 1L]] <<- cbind(vi = vi, vk = vk, cosmax = apply(cs_, 1, max))"),
                after = i - 1L)
  f2 <- eval(parse(text = txt))
  environment(f2) <- ns
  utils::assignInNamespace("jc_nonest", f2, ns = "frmtmb")
}
take <- function(expr) {
  rec <<- list()
  r <- cap(expr)
  r$rec <- if (length(rec)) rec[[length(rec)]] else NULL
  r
}
rs <- function(r) if (is.null(r)) "-" else
  paste(signif(r[, "vi"] / r[, "vk"], 3), collapse = " ")
rc <- function(r) if (is.null(r)) "-" else
  paste(signif(r[, "cosmax"], 3), collapse = " ")
cat("seed\tre\tx\tse\tref\tvi_vk\ttrue_share\twarn\tcosmax\n", file = out)
for (cs in list(c(3, 1), c(8, 1), c(10, 1), c(4, 0), c(5, 0))) {
  s <- cs[1]; re <- cs[2] == 1
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)),
                  z = stats::rnorm(200),
                  w = factor(sample(c("a", "b"), 200, TRUE)))
  d$y <- sin(3 * d$x) + 0.3 * d$z + 0.5 * (d$w == "b") +
    stats::rnorm(20, 0, 0.3)[d$g] + stats::rnorm(200, 0, 0.2)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  rhs <- paste(paste0("x", 1:8, collapse = " + "), "+ z + w")
  reterm <- if (re) "+ (1 | g)" else ""
  f <- cap(frm(stats::as.formula(paste("y ~", rhs, reterm)), data = d))$v
  fr <- suppressWarnings(frm(stats::as.formula(paste("y ~ poly(x, 8) + z + w",
                                                     reterm)), data = d))
  cat(sprintf("\n== seed %d re %s | code %d (ref code %d) | lost %s\n", s, re,
              f$opt$convergence, fr$opt$convergence, lost_of(f)))
  # well-determined terms: coefficients, sigma, sd
  for (h in c("b_z = 0", "b_wb = 0", "b_z - b_wb = 0")) {
    a <- take(hypothesis(f, h, class = NULL))
    b <- suppressWarnings(hypothesis(fr, h, class = NULL))
    cat(sprintf("  hyp %-15s se %-10s ref %-10s warn %d | share %s | cos %s\n", h,
                format(a$v$hypothesis$Est.Error, digits = 5),
                format(b$hypothesis$Est.Error, digits = 5), length(a$w),
                rs(a$rec), rc(a$rec)))
  }
  sm <- function(m) {
    t <- capture.output(print(summary(m)))
    t[grepl("^sigma |sd\\(Intercept\\)", trimws(t))]
  }
  cat("  summary lane:", paste(trimws(sm(f)), collapse = " | "), "\n")
  cat("  summary ref :", paste(trimws(sm(fr)), collapse = " | "), "\n")
  nd <- data.frame(x = c(1, 1.5, 2), z = c(-1, 0, 1), w = c("a", "b", "a"),
                   g = d$g[1])
  for (k in 1:8) nd[[paste0("x", k)]] <- nd$x^k
  rf <- if (re) NA else NULL
  a <- take(fitted(f, newdata = nd, re_formula = rf, dpar = "sigma"))
  b <- suppressWarnings(fitted(fr, newdata = nd, re_formula = rf,
                               dpar = "sigma"))
  cat("  fitted(dpar = sigma) se", format(a$v[, 2], digits = 4), "ref",
      format(b[, 2], digits = 4), "warn", length(a$w), "| share", rs(a$rec),
      "| cos", rc(a$rec), "\n")
  # the grid: share read against the true one
  ng <- data.frame(x = seq(0.9, 2.2, by = 0.05), z = 0, w = "a", g = d$g[1])
  for (k in 1:8) ng[[paste0("x", k)]] <- ng$x^k
  a <- take(fitted(f, newdata = ng, re_formula = rf))
  b <- suppressWarnings(fitted(fr, newdata = ng, re_formula = rf))
  vr <- a$rec
  sh <- if (is.null(vr)) rep(NA, nrow(ng)) else vr[, "vi"] / vr[, "vk"]
  ts <- if (is.null(vr)) rep(NA, nrow(ng)) else b[, 2]^2 / vr[, "vk"] - 1
  for (i in seq_len(nrow(ng))) {
    cat(s, re, ng$x[i], a$v[i, 2], b[i, 2], sh[i], ts[i], length(a$w),
        if (is.null(vr)) NA else vr[i, "cosmax"], sep = "\t", file = out, append = TRUE)
    cat("\n", file = out, append = TRUE)
  }
  cat(sprintf("  grid: %d rows, finite %d, warn %d | est share range %s | true share range %s\n",
              nrow(ng), sum(is.finite(a$v[, 2])), length(a$w),
              paste(signif(range(sh), 3), collapse = " to "),
              paste(signif(range(ts), 3), collapse = " to ")))
}
