# Estimates, SEs, stationarity and the objective away from the optimum, all
# against an independent exact marginal, on designs not on the lane's grid.
# Usage: Rscript dev/review-impsamp-accuracy.R <design> <seed> [outdir]
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs()
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-designs.R")
args <- commandArgs(TRUE)
nm <- args[[1]]
seed <- as.integer(args[[2]])
out <- if (length(args) >= 3) args[[3]] else
  "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-out/accuracy"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
de <- design_list[[nm]](seed)

gr <- split(seq_len(nrow(de$data)), de$data$g)
dl <- lapply(gr, function(i) de$data[i, , drop = FALSE])

# AGH at order 61 per group, with integrate() at the reference optimum as a
# cross-check of the rule itself.
gh <- gh_rule(61L)
exact_nll <- function(par) {
  sd <- de$sdof(par)
  if (!is.finite(sd) || sd <= 0) return(Inf)
  tot <- 0
  for (d in dl) {
    fn <- function(u) de$ll(par, d, u) + dnorm(u, 0, sd, log = TRUE)
    o <- optimize(function(u) -fn(u), c(-8 * sd - 15, 8 * sd + 15), tol = 1e-7)
    o2 <- nlminb(o$minimum, function(u) -fn(u))
    m <- if (o2$objective < o$objective) o2$par else o$minimum
    h <- 1e-3
    cur <- -(fn(m + h) - 2 * fn(m) + fn(m - h)) / h^2
    s <- if (is.finite(cur) && cur > 1e-8) 1 / sqrt(cur) else sd
    nodes <- m + sqrt(2) * s * gh$x
    a <- fn(nodes) + gh$x^2 + log(gh$w)
    tot <- tot - (max(a) + log(sum(exp(a - max(a)))) + log(sqrt(2) * s))
  }
  tot
}
exact_nll_int <- function(par) {
  tot <- 0
  for (d in dl) {
    tot <- tot - exact_group_1d(function(u) de$ll(par, d, u), de$sdof(par))
  }
  tot
}
num_grad <- function(f, x, h = 1e-4) {
  vapply(seq_along(x), function(i) {
    e <- replace(numeric(length(x)), i, h)
    (f(x + e) - f(x - e)) / (2 * h)
  }, 0)
}
num_hess <- function(f, x, h = 1e-3) {
  k <- length(x)
  H <- matrix(0, k, k)
  for (i in seq_len(k)) {
    ei <- replace(numeric(k), i, h)
    H[i, ] <- num_grad(f, x + ei, 1e-3) - num_grad(f, x - ei, 1e-3)
    H[i, ] <- H[i, ] / (2 * h)
  }
  0.5 * (H + t(H))
}

warn <- character(0)
t0 <- Sys.time()
fit <- withCallingHandlers(
  tryCatch(frm(de$formula, family = de$family, data = de$data,
               importance = de$draws),
           error = function(e) e),
  warning = function(w) {
    warn <<- c(warn, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
tfit <- as.numeric(Sys.time() - t0, units = "secs")
if (inherits(fit, "error")) {
  cat(sprintf("%s\t%d\tERROR\t%s\n", nm, seed,
              substr(conditionMessage(fit), 1, 300)),
      file = file.path(out, sprintf("%s-%d.tsv", nm, seed)))
  cat("ERROR:", conditionMessage(fit), "\n")
  quit("no")
}
par <- unname(fit$opt$par)
cat("names:", names(fit$opt$par), "\n")
cat("par:", format(par, digits = 6), "\n")
# reference optimum from a start at the fit; nlminb on the exact nll with
# central differences
ref <- nlminb(par, exact_nll,
              gradient = function(p) num_grad(exact_nll, p),
              control = list(iter.max = 200, rel.tol = 1e-10))
pe <- ref$par
He <- num_hess(exact_nll, pe)
Ve <- tryCatch(solve(He), error = function(e) matrix(NA, length(pe), length(pe)))
se_exact <- sqrt(diag(Ve))
se_fit <- tryCatch(sqrt(diag(as.matrix(vcov(fit, full = TRUE)))),
                   error = function(e) NULL)
if (is.null(se_fit) || length(se_fit) != length(pe)) {
  sdr <- fit$cache$sdr %||% tryCatch(frmtmb:::sdr_of(fit), error = function(e) NULL)
  se_fit <- if (!is.null(sdr)) sqrt(diag(sdr$cov.fixed)) else rep(NA, length(pe))
}
g_at_fit <- num_grad(exact_nll, par)
newton <- sqrt(max(0, sum(g_at_fit * (Ve %*% g_at_fit))))
ll_exact <- -ref$objective
int_check <- exact_nll_int(pe)
ll_imp <- -fit$opt$objective
lap <- suppressWarnings(tryCatch(
  frm(de$formula, family = de$family, data = de$data), error = function(e) e))
lap_shift <- if (inherits(lap, "error")) rep(NA, length(pe)) else
  (unname(lap$opt$par) - pe) / se_exact
res <- data.frame(
  design = nm, seed = seed, draws = de$draws, rounds = fit$importance$rounds,
  grad = fit$importance$grad, ess_min = fit$importance$ess_min,
  mcse = fit$importance$mcse, ll_imp = ll_imp, ll_exact = ll_exact,
  int_minus_agh = -int_check - ll_exact,
  ll_at_fit_exact = -exact_nll(par),
  max_z_est = max(abs(par - pe) / se_exact),
  max_z_lap = max(abs(lap_shift)),
  newton = newton,
  se_ratio_min = min(se_fit / se_exact), se_ratio_max = max(se_fit / se_exact),
  n_warn = length(warn), tfit = tfit,
  warn = substr(paste(warn, collapse = " | "), 1, 200))
print(t(res))
cat("est-ref/se:", format((par - pe) / se_exact, digits = 3), "\n")
cat("se_fit/se_exact:", format(se_fit / se_exact, digits = 3), "\n")
write.table(res, file.path(out, sprintf("%s-%d.tsv", nm, seed)), sep = "\t",
            quote = FALSE, row.names = FALSE)

# Away from the optimum: the objective against the exact, the AD gradient
# against a central difference of the fit's own objective, and against the
# exact gradient. The points are along each coordinate and jointly.
pts <- list(opt = par)
for (i in seq_along(par)) {
  for (k in c(-3, 3, -8, 8)) {
    pts[[paste0("c", i, "_", k)]] <- replace(par, i, par[i] + k * se_exact[i])
  }
}
pts$sd_x3 <- replace(par, length(par), par[length(par)] + log(3))
pts$sd_d3 <- replace(par, length(par), par[length(par)] - log(3))
pts$joint <- par + 3 * se_exact * rep_len(c(1, -1), length(par))
ns <- asNamespace("frmtmb")
prop <- ns$imp_frozen_proposal(fit)
rows <- lapply(names(pts), function(pn) {
  p <- pts[[pn]]
  named <- setNames(p, names(fit$opt$par))
  v <- fit$obj$fn(named)
  e <- exact_nll(p)
  eg <- ns$imp_ess_at(prop, named)
  ad <- as.numeric(fit$obj$gr(named))
  fd <- num_grad(function(q) fit$obj$fn(setNames(q, names(fit$opt$par))), p,
                 h = 1e-5)
  ge <- num_grad(exact_nll, p)
  data.frame(design = nm, seed = seed, pt = pn, fn = v, exact = e,
             diff = v - e, mcse_at = eg$mcse, ess_min_at = min(eg$ess),
             z = (v - e) / eg$mcse,
             ad_vs_fd = max(abs(ad - fd)) / max(1e-6, max(abs(fd))),
             ad_vs_exact = max(abs(ad - ge)) / max(1e-3, max(abs(ge))),
             stringsAsFactors = FALSE)
})
rows <- do.call(rbind, rows)
print(format(rows[, -(1:2)], digits = 4))
write.table(rows, file.path(out, sprintf("%s-%d-away.tsv", nm, seed)),
            sep = "\t", quote = FALSE, row.names = FALSE)
