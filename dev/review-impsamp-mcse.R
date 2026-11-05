# Item 5: is the reported MCSE ever too LOW? Full fits over importance_seed on
# one dataset, error against an independent exact value at the fit's own
# estimates. Usage: Rscript dev/review-impsamp-mcse.R <design> <dseed> <draws> [S]
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs()
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-designs.R")
args <- commandArgs(TRUE)
nm <- args[[1]]
dseed <- as.integer(args[[2]])
N <- as.integer(args[[3]])
S <- if (length(args) >= 4) as.integer(args[[4]]) else 100L
out <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-out/mcse"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

extra <- list(
  # integrands symmetric about their mode: mean 0.5 groups, intercept near 0.
  # The weight is then an even function of z, a pair shares it, and a pair is
  # one draw, so the iid formula should read low by sqrt(2). Half the groups
  # sit at 10 of 20, the others at 3 or 17 of 20 to keep the sd away from 0.
  sym20 = function(seed) {
    set.seed(seed)
    G <- 24
    y <- rep(c(10, 10, 3, 17), length.out = G)
    d <- data.frame(y = y, n = 20L, g = factor(seq_len(G)))
    d <- d[rep(seq_len(G), 1), ]
    list(data = d, formula = bf(y | trials(n) ~ 1 + (1 | g)),
         family = binomial(), draws = N,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             sum(dbinom(d$y, d$n, plogis(par[1] + ui), log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[2]))
  },
  # symmetric and heavy: student t with one row per group at mean 0 is even
  # around each mode only for y = 0; two opposite rows per group are even
  # about the group's midpoint
  symt = function(seed) {
    set.seed(seed)
    G <- 20
    a <- rnorm(G, 0, 1)
    d <- data.frame(y = c(rbind(a - 1.5, a + 1.5)),
                    g = factor(rep(seq_len(G), each = 2)))
    list(data = d, formula = bf(y ~ 1 + (1 | g), nu = 4), family = student(),
         draws = N,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             s <- exp(par[2])
             sum(dt((d$y - par[1] - ui) / s, 4, log = TRUE) - log(s))
           }, 0)
         },
         sdof = function(par) exp(par[3]))
  }
)
de <- c(design_list, extra)[[nm]](dseed)
de$draws <- N
dl <- split(de$data, de$data$g)
gh <- gh_rule(81L)
exact_nll <- function(par) {
  sd <- de$sdof(par)
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
rows <- list()
pars <- list()
for (s in seq_len(S)) {
  wn <- 0L
  fit <- withCallingHandlers(
    tryCatch(frm(de$formula, family = de$family, data = de$data,
                 importance = N,
                 control = frmtmb_control(importance_seed = s)),
             error = function(e) e),
    warning = function(w) {
      wn <<- wn + 1L
      invokeRestart("muffleWarning")
    })
  if (inherits(fit, "error")) {
    rows[[s]] <- data.frame(seed = s, err = NA, mcse = NA, ess = NA,
                            rounds = NA, nwarn = wn,
                            msg = substr(conditionMessage(fit), 1, 80))
    next
  }
  p <- unname(fit$opt$par)
  pars[[s]] <- p
  rows[[s]] <- data.frame(seed = s, err = -fit$opt$objective + exact_nll(p),
                          mcse = fit$importance$mcse,
                          ess = fit$importance$ess_min,
                          rounds = fit$importance$rounds, nwarn = wn, msg = "")
}
r <- do.call(rbind, rows)
ok <- is.finite(r$err)
P <- do.call(rbind, pars)
res <- data.frame(
  design = nm, dseed = dseed, draws = N, S = S, nok = sum(ok),
  mean_err = mean(r$err[ok]), sd_err = sd(r$err[ok]),
  med_mcse = median(r$mcse[ok]), sd_over_mcse = sd(r$err[ok]) / median(r$mcse[ok]),
  rms_over_rmcse = sqrt(mean(r$err[ok]^2)) / sqrt(mean(r$mcse[ok]^2)),
  cover95 = mean(abs(r$err[ok]) <= 1.96 * r$mcse[ok]),
  cover_mean_z = mean(r$err[ok] / r$mcse[ok]),
  max_abs_z = max(abs(r$err[ok] / r$mcse[ok])),
  ess_min = min(r$ess[ok]), ess_med = median(r$ess[ok]),
  warned = mean(r$nwarn[ok] > 0), max_rounds = max(r$rounds[ok]))
print(t(res))
cat("err quantiles", format(quantile(r$err[ok], c(0, .05, .5, .95, 1)), digits = 3), "\n")
cat("par sd over seeds:", format(apply(P, 2, sd), digits = 3), "\n")
write.table(res, file.path(out, sprintf("%s-%d-%d.tsv", nm, dseed, N)),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(r, file.path(out, sprintf("%s-%d-%d.rows.txt", nm, dseed, N)),
            sep = "\t", quote = FALSE, row.names = FALSE)
