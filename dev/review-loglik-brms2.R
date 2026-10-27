# Reviewer of lane loglik: frmtmb.sample's log_lik() against brms
# 2.23.0's log_lik() at the same parameter vectors, on cases the lane
# did not run. Seed 71. Output: dev/review-loglik-log/brms2.txt
source("C:/Users/adf44/source/r/frmtmb-wt-loglik/dev/loglik-brms-lib.R")
LIB <- Sys.getenv("RV_LIB", "C:/Users/adf44/source/r/wt-loglik-lib")
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.spatial)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"),
    find.package("frmtmb.spatial"), "\n")

draws_at <- function(fit, vals) {
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  k <- max(lengths(vals))
  m <- matrix(0, k, length(lab), dimnames = list(NULL, lab))
  for (nm in names(vals)) {
    if (!nm %in% lab) stop("no label ", nm, " in ", paste(lab, collapse = " "))
    m[, nm] <- vals[[nm]]
  }
  structure(list(stanfit = NULL, draws = m, fit = fit),
            class = "frmtmb_draws")
}
pinv <- function(r) r / sqrt(1 - r^2)
three <- function(f) lapply(1:3, f)
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
cmp <- function(nm, llf, llb) {
  cat(sprintf("%-40s %d x %d  rel = %.3e\n", nm, nrow(llf), ncol(llf),
              rel(llf, llb)))
}
try_case <- function(nm, expr) {
  r <- tryCatch(expr, error = function(e) {
    cat(sprintf("%-40s ERROR: %s\n", nm, conditionMessage(e)))
    NULL
  })
  invisible(r)
}

set.seed(71)
# unequal groups, factor levels NOT in sorted order, rows shuffled
sizes <- c(6, 2, 5, 1, 7, 4)
labs <- c("d", "b", "f", "a", "e", "c")
d <- do.call(rbind, lapply(seq_along(sizes), function(k) {
  data.frame(g = labs[k], t = seq_len(sizes[k]))
}))
d$g <- factor(d$g, levels = c("e", "a", "c", "f", "b", "d"))
d$x <- rnorm(nrow(d))
d$y <- 1 + 0.5 * d$x + rnorm(nrow(d)) + rnorm(6)[as.integer(d$g)]
d <- d[sample(nrow(d)), ]
rownames(d) <- NULL
n <- nrow(d)
cat("n =", n, " group sizes:", paste(table(d$g), collapse = " "), "\n")
b0 <- c(1.1, 0.9, 1.3)
b1 <- c(0.4, 0.6, 0.2)
sg <- c(0.8, 1.2, 1.0)
# 2. cosy(gr =) with one-row groups
for (fam in c("gaussian", "student")) {
  try_case(paste("cosy one-row groups", fam), {
    nus <- c(1.5, 4, 30)
    fit <- frm(y ~ x + cosy(t, g), family = get(fam)(), data = d,
               dry_run = "objective")
    dd <- fit$frame$autocor$y$d
    a <- 1 / (dd - 1)
    rho <- c(0.4, 0.05, 0.9)
    vals <- list(b_Intercept = b0, b_x = b1, sigma = sg,
                 thetaac_1 = qlogis((rho + a) / (1 + a)))
    if (fam == "student") vals$nu <- nus
    llf <- log_lik(draws_at(fit, vals))
    inits <- three(function(s) {
      p <- list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(d$x),
                sigma = sg[s], cosy = rho[s])
      if (fam == "student") p$nu <- nus[s]
      p
    })
    llb <- log_lik(ll_fixed_fit(brms::bf(y ~ x + cosy(t, g)),
                                brms::brmsfamily(fam), d, inits))
    cmp(paste("cosy one-row groups", fam), llf, llb)
  })
}

# 3. unstr() with missing times in some groups (and one-row group)

# 6. sar() with an ASYMMETRIC row-standardized W (a K / K' slip is
# invisible on a symmetric W), both types, nu small, rho near its bound
set.seed(74)
P <- matrix(0, n, n)
for (i in 1:n) {
  nb <- sample(setdiff(1:n, i), sample(1:2, 1))
  P[i, nb] <- 1
}
P <- pmax(P, t(P)); W <- P / rowSums(P)
cat("W asymmetric:", !isSymmetric(W), "\n")
for (type in c("lag", "error")) {
  for (fam in c("gaussian", "student")) {
    nm <- paste("sar asym W", type, fam)
    try_case(nm, {
      nus <- c(1.2, 3, 10)
      fo <- if (type == "lag") y ~ x + sar(W) else
        y ~ x + sar(W, type = "error")
      fit <- frm(fo, family = get(fam)(), data = d, data2 = list(W = W),
                 dry_run = "objective")
      aux <- fit$frame$autocor$y$aux
      rho <- c(0.4, -0.3, 0.99 * aux$hi)
      vals <- list(b_Intercept = b0, b_x = b1, sigma = sg,
                   thetaac_1 = qlogis((rho - aux$lo) / (aux$hi - aux$lo)))
      if (fam == "student") vals$nu <- nus
      llf <- log_lik(draws_at(fit, vals))
      inits <- three(function(s) {
        p <- list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(d$x),
                  sigma = sg[s])
        if (type == "lag") p$lagsar <- rho[s] else p$errorsar <- rho[s]
        if (fam == "student") p$nu <- nus[s]
        p
      })
      llb <- log_lik(ll_fixed_fit(brms::bf(fo), brms::brmsfamily(fam), d,
                                  inits, data2 = list(W = W)))
      cmp(nm, llf, llb)
    })
  }
}
cat("DONE\n")
