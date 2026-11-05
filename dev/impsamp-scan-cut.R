# Two scans for test designs: where is the likelihood's curvature negative at
# the mode; and the objective off the anchor against integrate(), coordinate
# by coordinate. Log dev/impsamp-out/scan-cut.txt
.libPaths(c("C:/Users/adf44/source/r/wt-impsamp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
sq <- function(e) suppressWarnings(suppressMessages(e))
parts <- function(f, fam, d) {
  fr <- frm(f, family = fam, data = d, dry_run = "frame")
  nll <- ns$build_objective(fr)
  tpl <- ns$make_start(fr, NULL, NULL)
  lap <- RTMB::MakeADFun(nll, tpl, random = "b", map = fr[["map"]], silent = TRUE)
  opt <- sq(nlminb(lap$par, lap$fn, lap$gr))
  lay <- ns$imp_layout(fr[["re_blocks"]])
  list(frame = fr, nll = nll, tpl = tpl, lap = lap, opt = opt, lay = lay,
       gmap = ns$imp_group_map(fr, lay))
}
for (shift in c(6, 12, 25)) for (s in 1:6) {
  set.seed(s)
  g <- factor(rep(1:15, each = 2L))
  x <- rnorm(30L)
  y <- 0.3 + 0.5 * x + rnorm(15L, 0, 0.7)[g] + 0.3 * rt(30L, 3L)
  # one row of each of three groups far from its partner
  y[c(1L, 7L, 13L)] <- y[c(1L, 7L, 13L)] + shift
  d <- data.frame(y = y, x = x, g = g)
  p <- tryCatch(parts(bf(y ~ x + (1 | g), nu = 2), student(), d),
                error = function(e) NULL)
  if (is.null(p)) next
  plan <- ns$imp_plan(p$lap, p$frame, p$lay, p$opt$par, 100L, 1L)
  fz <- plan$hlik[1, 1, ]
  cat(sprintf("student shift %g seed %d: cut groups %d of 15, min hlik %.3f, prec %.3f\n",
              shift, s, sum(fz == 0), min(fz), plan$prec[1, 1]))
}
set.seed(124)
n <- 200L
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                g = factor(sprintf("g%02d", sort(rep(1:10, length.out = n)))))
sample(c("a", "b", "c"), n, TRUE)
sample(0:3, n, TRUE)
u <- rnorm(10L, 0, 0.5)
v <- rnorm(10L, 0, 0.3)
d$eta <- 0.3 + 0.5 * d$x - 0.3 * d$z + u[as.integer(d$g)] + v[as.integer(d$g)] * d$x
d$ntr <- 10L
d$y <- rbinom(n, 10L, plogis(d$eta))
f <- bf(y | trials(ntr) ~ x + (1 | g))
p <- parts(f, binomial(), d)
lap <- frm(f, family = binomial(), data = d)
se <- sqrt(diag(ns$sdr_of(lap)$cov.fixed))
exact_nll <- function(par) {
  s1 <- exp(par[[3L]])
  -sum(vapply(split(seq_len(nrow(d)), d$g), function(r) {
    ff <- function(u) vapply(u, function(ui) {
      exp(sum(dbinom(d$y[r], 10L, plogis(par[[1L]] + par[[2L]] * d$x[r] + ui),
                     log = TRUE))) * dnorm(ui, 0, s1)
    }, 0)
    log(integrate(ff, -8 * s1, 8 * s1, subdivisions = 400L, rel.tol = 1e-9)$value)
  }, 0))
}
anchor <- p$opt$par
for (N in c(200L, 2000L)) for (s in 1:3) {
  plan <- ns$imp_plan(p$lap, p$frame, p$lay, anchor, N, s)
  io <- ns$build_importance_objective(p$frame, p$lay, p$gmap, plan)
  a0 <- as.matrix(io$amat(ns$imp_par_list(p$tpl, anchor)))
  io$set_offset(apply(a0, 1L, max))
  obj <- RTMB::MakeADFun(io$fn, ns$imp_template(p$tpl, "b"), silent = TRUE)
  for (sgn in list(c(0, 0, 0), c(2, 2, 2), c(-2, 2, -2), c(0, 0, 2), c(0, 0, -2),
                   c(2, 0, 0), c(0, 2, 0))) {
    at <- anchor + sgn * se
    es <- ns$imp_ess(as.matrix(io$amat(ns$imp_par_list(p$tpl, at))), N)
    cat(sprintf("N=%d seed %d sgn %s: err %+.5f mcse %.5f ess %.3f\n", N, s,
                paste(sgn, collapse = ","), obj$fn(at) - exact_nll(at), es$mcse,
                min(es$ess)))
  }
}
