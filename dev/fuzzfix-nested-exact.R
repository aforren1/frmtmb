# What does quadrature = TRUE compute over a nested (1 | ga/gb) block?
# Compares, at the quadrature optimum's outer parameters:
#   - the quadrature objective (fit$obj$fn),
#   - the Laplace objective,
#   - the exact iterated integral by Gauss-Hermite in R over frmtmb's
#     own joint (no random effects) objective.
# and lists which random effects the TMB transform left to Laplace.
#   Rscript dev/fuzzfix-nested-exact.R <lib|base> <case>
# case: cum (the fuzz spec data), bern (bernoulli nested, seed 11)
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
case <- a[2]
if (case == "cum") {
  d <- readRDS("dev/fuzzfix-out/data-cum.rds")
  bf0 <- bf(y | weights(w) ~ 1 + x + (1 | ga/gb)) + cumulative()
} else {
  set.seed(4)
  ng <- 20; nt <- 5; n <- ng * nt
  d <- data.frame(g = factor(rep(seq_len(ng), each = nt)),
                  gb = factor(rep(c("i", "ii"), length.out = n)),
                  x = rnorm(n))
  d$ga <- d$g
  eta <- 0.5 + 0.4 * d$x + rnorm(ng, 0, 0.8)[d$g] +
    rnorm(ng * 2, 0, 0.8)[as.integer(d$ga) * 2 + as.integer(d$gb) - 2]
  set.seed(11)
  d$y <- rbinom(n, 1, plogis(eta))
  bf0 <- bf(y ~ 1 + x + (1 | ga/gb)) + bernoulli()
}
fq <- suppressWarnings(frm(bf0, data = d, quadrature = TRUE))
fl <- suppressWarnings(frm(bf0, data = d))
fr <- fq$frame
for (bk in fr$re_blocks) {
  cat("block", bk$term_label, "group", bk$group_name, "levels",
      bk$n_levels, "b_idx", range(bk$b_idx), "\n")
}
o <- fq$obj
nb <- length(fq$estimates$b)
cat("b total:", nb, " left random in quadrature tape:",
    length(o$env$random), "\n")
# which b entries survived: the transform drops integrated ones from
# the tape's parameter vector, and the survivors keep the calibration
# template's values, which are the Laplace modes there
bl <- fr$re_blocks
p <- fq$opt$par
cat("quadrature objective at its optimum:", format(o$fn(p), digits = 12),
    "\n")
lo <- fl$obj
cat("laplace objective at quadrature optimum:",
    format(lo$fn(p), digits = 12), "\n")
cat("laplace fit objective:", format(fl$opt$objective, digits = 12), "\n")

# exact iterated integral at outer p
tpl <- fq$estimates
J <- RTMB::MakeADFun(frmtmb:::build_objective(fr), tpl, map = fr$map,
                     silent = TRUE)
jp <- J$par
jp[names(jp) != "b"] <- p
ib <- which(names(jp) == "b")
gh <- function(k) {
  i <- seq_len(k - 1)
  e <- eigen(matrix(0, k, k) + diag(0, k) +
               {m <- matrix(0, k, k); m[cbind(i, i + 1)] <- sqrt(i / 2);
                m[cbind(i + 1, i)] <- sqrt(i / 2); m}, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}
nodes <- gh(80)
# block sds at p, from the fit's own reporting
sds <- sapply(bl, function(bk) {
  th <- fq$estimates$theta[bk$theta_idx]
  sqrt(frmtmb:::covstruct_registry[[bk$covstruct]]$vcov(th, bk)[1, 1])
})
cat("block sds:", sds, "\n")
# group structure from the data: outer level of each row and the inner
# b index of each row, read off the Z matrix columns
Z <- Filter(Negate(is.null), lapply(fr$linpreds, `[[`, "Z"))[[1]]
stopifnot(!is.null(Z))
outer_blk <- which(vapply(bl, function(bk) bk$group_name == "ga", TRUE))
inner_blk <- setdiff(seq_along(bl), outer_blk)
J0 <- J$fn(jp)
fb <- function(set) { q <- jp; q[ib[set$i]] <- set$v; J$fn(q) - J0 }
total <- -J0
Zm <- as.matrix(Z)
for (u_i in bl[[outer_blk]]$b_idx) {
  rows <- which(Zm[, u_i] != 0)
  inner <- bl[[inner_blk]]$b_idx
  inner <- inner[colSums(Zm[rows, inner, drop = FALSE] != 0) > 0]
  cu <- sqrt(2) * sds[outer_blk]; cv <- sqrt(2) * sds[inner_blk]
  lu <- vapply(nodes$x, function(xu) {
    du <- fb(list(i = u_i, v = cu * xu))
    acc <- -du
    for (v_j in inner) {
      lv <- vapply(nodes$x, function(xv) {
        -(fb(list(i = c(u_i, v_j), v = c(cu * xu, cv * xv))) - du)
      }, 1)
      m <- max(lv + nodes$x^2)
      acc <- acc + log(cv) + m +
        log(sum(nodes$w * exp(lv + nodes$x^2 - m)))
    }
    acc
  }, 1)
  m <- max(lu + nodes$x^2)
  total <- total + log(cu) + m + log(sum(nodes$w * exp(lu + nodes$x^2 - m)))
}
cat("exact iterated integral objective (-loglik):",
    format(-total, digits = 12), "\n")
cat("quadrature - exact:", format(o$fn(p) + total, digits = 4), "\n")
cat("laplace - exact:", format(lo$fn(p) + total, digits = 4), "\n")

# Two passes: the inner block first, each of its effects a 1-D integral
# given the outer one, then the outer block over the inner integrals.
gkc <- function(which) structure(list(dim = 1, adaptive = FALSE,
                                      debug = FALSE,
                                      method = "marginal_gk",
                                      which = which), class = "GK")
nll <- frmtmb:::build_objective(fr)
tpl2 <- fq$estimates[names(fr$par_template)]
ib_in <- bl[[inner_blk]]$b_idx
ib_out <- bl[[outer_blk]]$b_idx
two <- tryCatch(RTMB::MakeADFun(nll, tpl2, random = "b", map = fr$map,
  integrate = list(list(b = gkc(ib_in)), list(b = gkc(seq_along(ib_out)))),
  silent = TRUE), error = function(e) e)
if (inherits(two, "error")) {
  cat("two-pass tape error:", conditionMessage(two), "\n")
} else {
  cat("two-pass left random:", length(two$env$random), "\n")
  f2 <- two$fn(p)
  cat("two-pass objective:", format(f2, digits = 12), "\n")
  cat("two-pass - exact:", format(f2 + total, digits = 4), "\n")
  g2 <- two$gr(p)
  cat("two-pass gradient finite:", all(is.finite(g2)), "\n")
}
