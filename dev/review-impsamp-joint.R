# smallpois seed 7, the far "joint" point: mean ratio 7 with a small error in
# one group. Is exact_group_1d() wrong there, or the estimator?
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs()
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-designs.R")
de <- design_list$smallpois(7L)
ns <- asNamespace("frmtmb")
fr <- frm(de$formula, family = de$family, data = de$data, dry_run = "frame")
nll <- ns$build_objective(fr)
tpl <- fr[["par_template"]]
lap <- RTMB::MakeADFun(nll, tpl, random = "b", map = fr[["map"]], silent = TRUE)
opt <- suppressWarnings(nlminb(lap$par, lap$fn, lap$gr))
lay <- ns$imp_layout(fr[["re_blocks"]])
gmap <- ns$imp_group_map(fr, lay)
anchor <- setNames(opt$par, names(lap$par))
se <- sqrt(diag(solve(numDeriv::jacobian(lap$gr, opt$par))))
p <- anchor + 4 * se * rep_len(c(1, -1), 3)
cat("anchor", anchor, "p", p, "sd", exp(p[3]), "\n")
gr <- split(seq_len(nrow(de$data)), de$data$g)
dl <- lapply(gr, function(i) de$data[i, , drop = FALSE])
piece <- function(d, par) {
  sd <- exp(par[3])
  f <- function(u) exp(de$ll(unname(par), d, u) + dnorm(u, 0, sd, log = TRUE))
  br <- seq(-10 * sd - 8, 10 * sd + 8, length.out = 201L)
  v <- 0
  for (i in seq_len(length(br) - 1L)) {
    v <- v + integrate(f, br[i], br[i + 1L], rel.tol = 1e-10, abs.tol = 0,
                       stop.on.error = FALSE)$value
  }
  log(v)
}
eg <- vapply(dl, function(d) {
  exact_group_1d(function(u) de$ll(unname(p), d, u), exp(p[3]))
}, 0)
eg2 <- vapply(dl, piece, 0, par = p)
print(round(cbind(mine1 = eg, piece = eg2, diff = eg - eg2), 4)[
  order(-abs(eg - eg2))[1:8], ])
S <- 60L
N <- 200L
r <- matrix(0, length(eg), S)
for (s in seq_len(S)) {
  plan <- ns$imp_plan(lap, fr, lay, anchor, N, s)
  io <- ns$build_importance_objective(fr, lay, gmap, plan)
  pl <- ns$imp_par_list(tpl, p, fr[["map"]])
  a <- as.matrix(io$amat(pl))
  pr <- environment(io$fn)$proposal(pl)
  lm <- apply(a, 1L, log_mean_exp) + 0.5 * log(2 * pi) + as.numeric(pr$logdet)
  r[, s] <- exp(lm - eg2)
}
mr <- rowMeans(r)
se_r <- apply(r, 1, sd) / sqrt(S)
print(round(cbind(n_y0 = vapply(dl, function(d) sum(d$y), 0),
                  mean_ratio = mr, se = se_r, med = apply(r, 1, median),
                  rate = plan$hlik_rate[, 1], rate3 = plan$hlik_rate[, 3],
                  hlik = plan$hlik[1, 1, ]), 4)[order(-abs(mr - 1))[1:12], ])
cat("total log-estimate - exact: mean", mean(colSums(log(r))), "sd",
    sd(colSums(log(r))), "\n")
