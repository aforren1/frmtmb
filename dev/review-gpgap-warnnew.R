# Short check of fix round 1, item 3: the warning on designs its calibration
# did not see: new seeds (9001 up), new fields, clustered and unbalanced x, new
# families and structures. "Moved" as the lane defines it: any outer parameter
# at 1 SE of the reference, or 0.5 nats either way. The reference is the same
# tape non-centered at brms's line (relative 1e-12 away from y's unit), best of
# 3 starts.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-warnnew.R <log> <chunk> <chunks>
a <- as.integer(commandArgs(trailingOnly = TRUE))
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
orig <- get("gp_corr", envir = ns)
put <- function(v) {
  unlockBinding("gp_corr", ns); assign("gp_corr", v, envir = ns); lockBinding("gp_corr", ns)
}
corr_at <- function(kind) function(theta, blk) {
  Q <- 0
  for (j in seq_along(blk[["aux_D2"]])) {
    rho <- if (isTRUE(blk[["gp_iso"]])) exp(theta[2]) else exp(theta[1 + j])
    Q <- Q + blk[["aux_D2"]][[j]] / (2 * rho^2)
  }
  I <- diag(1, nrow(blk[["aux_D2"]][[1]]))
  if (kind == "abs") exp(-Q) + I * (1e-12 * exp(-2 * theta[1])) else exp(-Q) + I * 1e-12
}
fields <- list(
  tanh = function(x) 2 * tanh(2 * (x - 3)),
  two = function(x) sin(1.5 * x) + 0.5 * cos(4 * x),
  damp = function(x) exp(-0.4 * x) * sin(3 * x) * 3,
  bump = function(x) 3 * exp(-(x - 2)^2) - 2 * exp(-(x - 4.5)^2 / 0.5))
xs <- list(
  unif = function(n) round(runif(n, 0, 6), 4),
  clust = function(n) round(c(runif(n %/% 2, 0.5, 1.5), runif(n - n %/% 2, 4, 5.5)), 4),
  skew = function(n) round(6 * rbeta(n, 1, 3), 4))
des <- list()
add <- function(tag, fn, ref = "abs") des[[length(des) + 1L]] <<- list(tag = tag, fn = fn, ref = ref)
seed <- 9000L
for (fl in names(fields)) for (xk in names(xs)) for (r in c(0.1, 0.02, 0.008, 0.004, 0.002, 0.001, 3e-4)) {
  seed <- seed + 1L
  local({
    s0 <- seed; fl <- fl; xk <- xk; r <- r
    add(sprintf("gauss %s %s r%g s%d", fl, xk, r, s0), function() {
      set.seed(s0)
      d <- data.frame(x = xs[[xk]](50))
      d$y <- fields[[fl]](d$x) + r * rnorm(50)
      frm(bf(y ~ gp(x)), data = d)
    })
  })
}
mk <- function(tag, code, ref = "abs") {
  seed <<- seed + 1L
  s0 <- seed
  add(tag, function() {
    set.seed(s0)
    n <- 60
    d <- data.frame(x = round(runif(n, 0, 6), 3), z = round(runif(n, 0, 6), 3),
                    f = factor(rep(c("a", "b"), length.out = n)),
                    g = factor(rep(1:6, length.out = n)))
    d$e <- rnorm(n)
    code(d)
  }, ref)
}
for (sh in c(1e2, 1e4, 1e6, 1e8)) local({
  sh <- sh
  mk(sprintf("gamma shape %g", sh), function(d) {
    d$y <- rgamma(nrow(d), shape = sh, rate = sh / exp(sin(2 * d$x)))
    frm(bf(y ~ gp(x)), family = Gamma(link = "log"), data = d)
  })
  mk(sprintf("beta phi %g", sh), function(d) {
    m <- plogis(sin(2 * d$x))
    d$y <- rbeta(nrow(d), m * sh, (1 - m) * sh)
    frm(bf(y ~ gp(x)), family = brms::Beta(), data = d)
  })
  mk(sprintf("negbin size %g mean e^12", sh), function(d) {
    d$y <- rnbinom(nrow(d), size = sh, mu = exp(12 + sin(2 * d$x)))
    frm(bf(y ~ gp(x)), family = brms::negbinomial(), data = d)
  })
})
for (lm in c(7, 9, 10, 11, 13, 15)) local({
  lm <- lm
  mk(sprintf("poisson log mean %g", lm), function(d) {
    d$y <- rpois(nrow(d), exp(lm + sin(1.5 * d$x)))
    frm(bf(y ~ gp(x)), family = poisson(), data = d)
  })
})
for (tr in c(3e3, 3e4, 1e5, 3e5, 3e6)) local({
  tr <- tr
  for (k in 1:3) mk(sprintf("binomial trials %g rep %d", tr, k), function(d) {
    d$tr <- tr
    d$y <- rbinom(nrow(d), d$tr, plogis(sin(1.5 * d$x)))
    frm(bf(y | trials(tr) ~ gp(x)), family = binomial(), data = d)
  })
})
for (r in c(0.03, 0.01, 0.004, 0.001)) local({
  r <- r
  mk(sprintf("2-D iso r%g", r), function(d) {
    d$y <- sin(d$x) * cos(0.7 * d$z) + r * d$e
    frm(bf(y ~ gp(x, z)), data = d)
  })
  mk(sprintf("2-D aniso r%g", r), function(d) {
    d$y <- sin(d$x) * cos(0.7 * d$z) + r * d$e
    frm(bf(y ~ gp(x, z, iso = FALSE)), data = d)
  })
  mk(sprintf("by f r%g", r), function(d) {
    d$y <- ifelse(d$f == "a", sin(2 * d$x), cos(1.5 * d$x)) + r * d$e
    frm(bf(y ~ gp(x, by = f)), data = d)
  })
  mk(sprintf("two terms r%g", r), function(d) {
    d$y <- sin(2 * d$x) + 0.3 * cos(d$z) + r * d$e
    frm(bf(y ~ gp(x) + gp(z)), data = d)
  })
  mk(sprintf("mvbf r%g", r), function(d) {
    d$y <- sin(2 * d$x) + r * d$e
    d$y2 <- 0.5 * d$z + rnorm(nrow(d), 0, 0.3)
    frm(mvbf(bf(y ~ gp(x)), bf(y2 ~ z), rescor = FALSE), data = d)
  })
  mk(sprintf("gp plus (1|g) r%g", r), function(d) {
    d$y <- sin(2 * d$x) + r * d$e
    frm(bf(y ~ gp(x) + (1 | g)), data = d)
  })
  mk(sprintf("skew_normal r%g", r), function(d) {
    d$y <- sin(2 * d$x) + r * d$e
    frm(bf(y ~ gp(x)), family = brms::skew_normal(), data = d)
  })
  # unbalanced replicates: one position holds a third of the rows
  mk(sprintf("unbalanced reps r%g", r), function(d) {
    d$x[1:20] <- d$x[1]
    d$y <- sin(2 * d$x) + r * d$e
    frm(bf(y ~ gp(x)), data = d)
  })
})
des <- des[seq_along(des) %% a[2] == a[1] %% a[2]]
for (de in des) {
  warned <- FALSE; err <- NULL
  fr <- withCallingHandlers(
    tryCatch(de$fn(), error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }),
    warning = function(w) {
      if (inherits(w, "frmtmb_gp_nugget")) warned <<- TRUE
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  if (is.null(fr)) {
    cat(sprintf("WN %s | ERROR %s\n", de$tag, substr(err, 1, 90)))
    next
  }
  bks <- which(vapply(fr$frame$re_blocks, function(b) identical(b$covstruct, "gp"), NA))
  rd <- tryCatch({
    vapply(bks, function(i) {
      bk <- fr$frame$re_blocks[[i]]
      inf <- frmtmb:::gp_position_info(fr, bk)
      nv <- 1e-6 * exp(2 * fr$estimates$theta[bk$theta_idx][[1]])
      stats::quantile(nv * inf, 0.9, names = FALSE) / sqrt(2 / length(inf))
    }, 0)
  }, error = function(e) NA_real_)
  frn <- fr$frame
  frn$ncp_blocks <- bks
  frn$ncp_blocks <- bks
  put(corr_at(de$ref))
  best <- NULL
  th0 <- fr$estimates$theta
  for (st in 1:3) {
    r <- tryCatch({
      fn <- frmtmb::build_objective(frn)
      tp <- frn$par_template
      tp$b[] <- 0
      if (st == 2L) tp$theta[] <- th0
      if (st == 3L) tp$theta[] <- th0 + 0.5
      o <- RTMB::MakeADFun(fn, tp, random = "b", map = frn$map, silent = TRUE)
      op <- suppressWarnings(nlminb(o$par, o$fn, o$gr,
                                    control = list(iter.max = 800, eval.max = 1200)))
      list(o = o, op = op)
    }, error = function(e) NULL)
    if (!is.null(r) && is.finite(r$op$objective) &&
          (is.null(best) || r$op$objective < best$op$objective)) best <- r
  }
  line <- if (is.null(best)) {
    "REF FAILED"
  } else {
    sdr <- tryCatch(suppressWarnings(RTMB::sdreport(best$o, best$op$par)), error = function(e) NULL)
    se <- if (!is.null(sdr)) sqrt(diag(sdr$cov.fixed)) else rep(NA_real_, length(best$op$par))
    dz <- abs(fr$opt$par - best$op$par) / se
    dz <- dz[is.finite(dz)]
    gap <- fr$opt$objective - best$op$objective
    moved <- isTRUE(max(c(dz, 0)) >= 1) || abs(gap) >= 0.5
    sprintf("gap %.3f dzmax %.2f moved %s refconv %d refgr %.1e", gap, max(c(dz, 0)),
            moved, best$op$convergence, max(abs(best$o$gr(best$op$par))))
  }
  put(orig)
  cat(sprintf("WN %s | m %s conv %d warned %s r90/thr %s | %s\n", de$tag,
              paste(vapply(bks, function(i) nrow(fr$frame$re_blocks[[i]]$aux_D2[[1]]), 0L),
                    collapse = ","),
              fr$opt$convergence, warned, paste(signif(rd, 3), collapse = ","), line))
}
