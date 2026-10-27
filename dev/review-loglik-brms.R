# Reviewer of lane loglik: frmtmb.sample's log_lik() against brms
# 2.23.0's log_lik() at the same parameter vectors, on cases the lane
# did not run. Seed 71. Output: dev/review-loglik-log/brms.txt
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

# 1. ar(p = 1) with gr =, ragged, shuffled; nu small (1.2)
for (fam in c("gaussian", "student")) {
  try_case(paste("ar gr shuffled", fam), {
    a <- c(0.5, -0.3, 0.95)
    nus <- c(1.2, 2.5, 1.05)
    fit <- frm(y ~ x + ar(t, g, cov = TRUE), family = get(fam)(), data = d,
               dry_run = "objective")
    vals <- list(b_Intercept = b0, b_x = b1, sigma = sg, thetaac_1 = pinv(a))
    if (fam == "student") vals$nu <- nus
    llf <- log_lik(draws_at(fit, vals))
    inits <- three(function(s) {
      p <- list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(d$x),
                sigma = sg[s] * sqrt(1 - a[s]^2), ar = array(a[s], 1))
      if (fam == "student") p$nu <- nus[s]
      p
    })
    bfit <- ll_fixed_fit(brms::bf(y ~ x + ar(t, g, cov = TRUE)),
                         brms::brmsfamily(fam), d, inits)
    llb <- log_lik(bfit)
    cmp(paste("ar gr shuffled", fam), llf, llb)
    # row order: brms returns the data's order; a sorted reading would
    # show here as a large difference
    o <- order(d$g, d$t)
    cat(sprintf("%-40s rel vs brms columns sorted by (g, t) = %.3e\n",
                "  (order probe)", rel(llf, llb[, o])))
  })
}

# 2. cosy(gr =) with one-row groups
for (fam in c("gaussian", "student")) {
  try_case(paste("cosy one-row groups", fam), {
    nus <- c(1.5, 4, 30)
    fit <- frm(y ~ x + cosy(t, g), family = get(fam)(), data = d,
               dry_run = "objective")
    dd <- fit$frame$autocor$y$d
    a <- 1 / (dd - 1)
    rho <- c(0.4, -0.1, 0.9)
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
dgap <- d[!(d$g == "e" & d$t %in% c(2, 5)) & !(d$g == "d" & d$t == 1), ]
rownames(dgap) <- NULL
for (fam in c("gaussian", "student")) {
  try_case(paste("unstr gaps one-row", fam), {
    nus <- c(1.5, 4, 30)
    fit <- frm(y ~ x + unstr(t, g), family = get(fam)(), data = dgap,
               dry_run = "objective")
    ac <- fit$frame$autocor$y
    set.seed(72)
    ths <- three(function(s) rnorm(ac$npar, 0, 0.4))
    vals <- list(b_Intercept = b0, b_x = b1, sigma = sg)
    for (k in seq_len(ac$npar)) {
      vals[[paste0("thetaac_", k)]] <- vapply(ths, `[`, 0, k)
    }
    if (fam == "student") vals$nu <- nus
    llf <- log_lik(draws_at(fit, vals))
    inits <- three(function(s) {
      p <- list(b = array(b1[s], 1),
                Intercept = b0[s] + b1[s] * mean(dgap$x), sigma = sg[s],
                Lcortime = as.matrix(frmtmb:::us_chol_L(ths[[s]], ac$d)))
      if (fam == "student") p$nu <- nus[s]
      p
    })
    llb <- log_lik(ll_fixed_fit(brms::bf(y ~ x + unstr(t, g)),
                                brms::brmsfamily(fam), dgap, inits))
    cmp(paste("unstr gaps one-row", fam), llf, llb)
  })
}

# 4. ar() on an interior gap: the documented difference
try_case("ar interior gap", {
  a <- c(0.5, -0.3, 0.8)
  fit <- frm(y ~ x + ar(t, g, cov = TRUE), data = dgap,
             dry_run = "objective")
  vals <- list(b_Intercept = b0, b_x = b1, sigma = sg, thetaac_1 = pinv(a))
  llf <- log_lik(draws_at(fit, vals))
  inits <- three(function(s) {
    list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(dgap$x),
         sigma = sg[s] * sqrt(1 - a[s]^2), ar = array(a[s], 1))
  })
  llb <- log_lik(ll_fixed_fit(brms::bf(y ~ x + ar(t, g, cov = TRUE)),
                              brms::brmsfamily("gaussian"), dgap, inits))
  cmp("ar interior gap (expected to differ)", llf, llb)
  ok <- !(dgap$g %in% c("e"))   # group d lost t = 1, not interior
  cat(sprintf("%-40s rel on groups with no interior gap = %.3e\n",
              "  ", rel(llf[, ok], llb[, ok])))
})

# 5. fcor() with a dense random SPD M (not Toeplitz), sigma ~ x, small nu
set.seed(73)
A <- matrix(rnorm(n * n), n, n)
M <- crossprod(A) / n + diag(0.2, n)
M <- (M + t(M)) / 2
for (fam in c("gaussian", "student")) {
  try_case(paste("fcor dense M", fam), {
    nus <- c(1.1, 2, 7)
    fit <- frm(y ~ x + fcor(M), family = get(fam)(), data = d,
               data2 = list(M = M), dry_run = "objective")
    vals <- list(b_Intercept = b0, b_x = b1, sigma = sg)
    if (fam == "student") vals$nu <- nus
    llf <- log_lik(draws_at(fit, vals))
    inits <- three(function(s) {
      p <- list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(d$x),
                sigma = sg[s])
      if (fam == "student") p$nu <- nus[s]
      p
    })
    llb <- log_lik(ll_fixed_fit(brms::bf(y ~ x + fcor(M)),
                                brms::brmsfamily(fam), d, inits,
                                data2 = list(M = M)))
    cmp(paste("fcor dense M", fam), llf, llb)
  })
}
try_case("fcor dense M sigma~x", {
  bs0 <- c(-0.2, 0.1, 0)
  bs1 <- c(0.3, -0.2, 0.1)
  fit <- frm(bf(y ~ x + fcor(M), sigma ~ x), data = d,
             data2 = list(M = M), dry_run = "objective")
  vals <- list(b_Intercept = b0, b_x = b1, b_sigma_Intercept = bs0,
               b_sigma_x = bs1)
  llf <- log_lik(draws_at(fit, vals))
  inits <- three(function(s) {
    list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(d$x),
         b_sigma = array(bs1[s], 1),
         Intercept_sigma = bs0[s] + bs1[s] * mean(d$x))
  })
  llb <- log_lik(ll_fixed_fit(brms::bf(y ~ x + fcor(M), sigma ~ x),
                              brms::brmsfamily("gaussian"), d, inits,
                              data2 = list(M = M)))
  cmp("fcor dense M sigma~x", llf, llb)
})

# 6. sar() with an ASYMMETRIC row-standardized W (a K / K' slip is
# invisible on a symmetric W), both types, nu small, rho near its bound
set.seed(74)
P <- matrix(0, n, n)
for (i in 1:n) {
  nb <- sample(setdiff(1:n, i), sample(1:4, 1))
  P[i, nb] <- 1
}
W <- P / rowSums(P)
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

# 7. multivariate, autocor in one response
d$y2 <- 0.3 * d$x + rnorm(n)
try_case("mv ar in y only", {
  a <- c(0.5, -0.3, 0.8)
  fit <- frm(bf(y ~ x + ar(t, g, cov = TRUE)) + bf(y2 ~ x) +
               set_rescor(FALSE), data = d, dry_run = "objective")
  lab <- frmtmb::brms_par_labels(fit)
  cat("mv labels:", paste(lab, collapse = " "), "\n")
  th <- grep("^thetaac", lab, value = TRUE)
  vals <- list(b_y_Intercept = b0, b_y_x = b1, sigma_y = sg,
               b_y2_Intercept = c(0, 0.1, -0.1), b_y2_x = c(0.3, 0.2, 0.4),
               sigma_y2 = c(1, 0.9, 1.1))
  vals[[th]] <- pinv(a)
  ds <- draws_at(fit, vals)
  inits <- three(function(s) {
    list(b_y = array(b1[s], 1), Intercept_y = b0[s] + b1[s] * mean(d$x),
         sigma_y = sg[s] * sqrt(1 - a[s]^2), ar_y = array(a[s], 1),
         b_y2 = array(c(0.3, 0.2, 0.4)[s], 1),
         Intercept_y2 = c(0, 0.1, -0.1)[s] + c(0.3, 0.2, 0.4)[s] * mean(d$x),
         sigma_y2 = c(1, 0.9, 1.1)[s])
  })
  bfit <- ll_fixed_fit(brms::bf(y ~ x + ar(t, g, cov = TRUE)) +
                         brms::bf(y2 ~ x) + brms::set_rescor(FALSE),
                       brms::brmsfamily("gaussian"), d, inits)
  cmp("mv resp = y", log_lik(ds, resp = "y"), log_lik(bfit, resp = "y"))
  cmp("mv resp = y2", log_lik(ds, resp = "y2"), log_lik(bfit, resp = "y2"))
  cmp("mv combined", log_lik(ds), log_lik(bfit))
})
cat("DONE\n")
