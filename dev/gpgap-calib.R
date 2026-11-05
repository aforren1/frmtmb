# Calibrates the nugget warning: per design, the nugget's variance against
# each position's working noise (the inverse of the likelihood's curvature in
# the field, read off the tape), and how far frm() sits from the same model
# non-centered at brms's line (at a relative 1e-12 where y has been rescaled,
# since brms's absolute line fails there: review 2026-11-02, M7).
#
#   Rscript dev/gpgap-calib.R <chunk> <n chunks>     # one CALIB line per row
#
# The designs are the review's (floor, rep, silent, trends) and this lane's.
# Seeds are in the design list. Base build by default: GPGAP_LIB.
args <- as.integer(commandArgs(trailingOnly = TRUE))
chunk <- args[1]
nch <- args[2]
.libPaths(c(strsplit(Sys.getenv("GPGAP_LIB",
                                "C:/Users/adf44/source/r/rellib-r16"), ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
orig <- get("gp_corr", envir = ns)
put <- function(val) {
  unlockBinding("gp_corr", ns)
  assign("gp_corr", val, envir = ns)
  lockBinding("gp_corr", ns)
}
corr_at <- function(kind) {
  function(theta, blk) {
    Q <- 0
    for (j in seq_along(blk[["aux_D2"]])) {
      rho <- if (isTRUE(blk[["gp_iso"]])) exp(theta[2]) else exp(theta[1 + j])
      Q <- Q + blk[["aux_D2"]][[j]] / (2 * rho^2)
    }
    I <- diag(1, nrow(blk[["aux_D2"]][[1]]))
    if (kind == "abs") exp(-Q) + I * (1e-12 * exp(-2 * theta[1])) else
      exp(-Q) + I * 1e-12
  }
}

# ------------------------------------------------------------ the designs
dd <- function(n, xs = 6, seed = 1, rep = 1L) {
  set.seed(seed)
  x <- round(runif(n %/% rep, 0, xs), 6)
  x <- rep(x, each = rep)
  x <- x[sample(length(x))]
  list(x = x, e = rnorm(length(x)), n = length(x))
}
des <- list()
add <- function(tag, fn, ref = "abs") {
  des[[length(des) + 1L]] <<- list(tag = tag, fn = fn, ref = ref)
}
gauss <- function(tag, n, ratio, seed, A = 1, xs = 6, rep = 1L, ref = "abs",
                  fam = "gaussian") {
  force(n); force(ratio); force(seed); force(A); force(xs); force(rep)
  add(tag, function() {
    g <- dd(n, xs, seed, rep)
    d <- data.frame(x = g$x)
    d$y <- A * sin(2 * d$x * 6 / xs) + ratio * A * g$e
    if (fam == "student") {
      frm(bf(y ~ gp(x)), family = student(), data = d)
    } else {
      frm(bf(y ~ gp(x)), data = d)
    }
  }, ref)
}
# unreplicated ladders: the review's floor designs, several n and scales
for (r in c(1e-1, 3e-2, 2e-2, 1.2e-2, 1e-2, 7e-3, 5e-3, 3.5e-3, 2.5e-3, 2e-3,
            1.5e-3, 1e-3, 3e-4, 1e-4)) {
  for (s in c(101L, 111L)) {
    for (n in c(30L, 60L, 100L, 200L)) {
      gauss(sprintf("floor n%d r%g s%d", n, r, s), n, r, s + n)
    }
  }
  gauss(sprintf("floor n300 r%g", r), 300L, r, 131L)
  gauss(sprintf("scale A1000 n60 r%g", r), 60L, r, 121L, A = 1000, ref = "rel")
  gauss(sprintf("scale A1e-3 n60 r%g", r), 60L, r, 122L, A = 1e-3, ref = "rel")
  gauss(sprintf("xs10 n100 r%g", r), 100L, r, 141L, xs = 10)
}
# replicated positions: the review's M5
for (r in c(1e-1, 5e-2, 3e-2, 2e-2, 1.4e-2, 1e-2, 5e-3, 2e-3, 1e-3, 3e-4)) {
  gauss(sprintf("rep 20x30 r%g", r), 600L, r, 131L, rep = 30L)
  gauss(sprintf("rep 60x4 r%g", r), 240L, r, 132L, rep = 4L)
  gauss(sprintf("rep 30x10 r%g", r), 300L, r, 133L, rep = 10L)
}
# the review's M2 rows, on its seed 5601
sil <- function() {
  set.seed(5601)
  n <- 60
  x <- round(runif(n, 0, 6), 3)
  e <- rnorm(n)
  list(d = data.frame(x = x, z = round(runif(n, 0, 6), 3),
                      g = factor(rep(1:6, length.out = n))), e = e, n = n)
}
for (lm in c(2, 4, 6, 8, 10, 11, 12, 13, 14, 16)) {
  add(sprintf("poisson log mean %g", lm), local({
    lm <- lm
    function() {
      s <- sil()
      d <- s$d
      d$cnt <- rpois(s$n, exp(lm + sin(2 * d$x)))
      frm(bf(cnt ~ gp(x)), family = poisson(), data = d)
    }
  }))
}
for (tr in c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7)) {
  add(sprintf("binomial trials %g", tr), local({
    tr <- tr
    function() {
      s <- sil()
      d <- s$d
      d$tr <- tr
      d$suc <- rbinom(s$n, d$tr, plogis(sin(2 * d$x)))
      frm(bf(suc | trials(tr) ~ gp(x)), family = binomial(), data = d)
    }
  }))
}
for (se in c(1e-1, 1e-2, 5e-3, 3e-3, 2e-3, 1e-3, 3e-4)) {
  add(sprintf("se(%g) known", se), local({
    se <- se
    function() {
      s <- sil()
      d <- s$d
      d$se <- se
      d$ys <- sin(2 * d$x) + se * s$e
      frm(bf(ys | se(se) ~ gp(x)), data = d)
    }
  }))
}
for (no in c(1e-1, 1e-2, 3e-3, 1e-3, 3e-4)) {
  mk <- function(fun) {
    no <- no
    force(fun)
    function() {
      s <- sil()
      d <- s$d
      d$y <- sin(2 * d$x) + no * s$e
      d$yy <- d$y + 0.1 * d$z
      d$no <- no
      fun(d)
    }
  }
  add(sprintf("sigma~(1|g) noise %g", no),
      mk(function(d) frm(bf(y ~ gp(x), sigma ~ (1 | g)), data = d)))
  add(sprintf("sigma~s(z) noise %g", no),
      mk(function(d) frm(bf(y ~ gp(x), sigma ~ s(z, k = 5)), data = d)))
  add(sprintf("mixture noise %g", no),
      mk(function(d) frm(bf(y ~ 1, mu1 ~ gp(x), mu2 ~ 1),
                         family = mixture(gaussian(), gaussian()), data = d)))
  add(sprintf("nl b1~gp(x) noise %g", no),
      mk(function(d) frm(bf(yy ~ b1 + b2 * z, b1 ~ gp(x), b2 ~ 1, nl = TRUE),
                         data = d)))
  add(sprintf("student noise %g", no),
      mk(function(d) frm(bf(y ~ gp(x)), family = student(), data = d)))
  add(sprintf("lognormal noise %g", no),
      mk(function(d) {
        d$yl <- exp(d$y)
        frm(bf(yl ~ gp(x)), family = lognormal(), data = d)
      }))
  add(sprintf("sigma offset noise %g", no),
      mk(function(d) {
        d$wsm <- d$no
        frm(bf(y ~ gp(x), sigma ~ 0 + offset(log(wsm))), data = d)
      }))
  add(sprintf("sigma ~ z noise %g", no),
      mk(function(d) frm(bf(y ~ gp(x), sigma ~ z), data = d)))
}
# sigma low on a share of the rows only
for (share in c(1 / 60, 0.1, 0.25, 0.5)) {
  add(sprintf("sigma low on share %.3f", share), local({
    share <- share
    function() {
      s <- sil()
      d <- s$d
      k <- seq_len(ceiling(share * s$n))
      d$low <- factor(ifelse(seq_len(s$n) %in% k, "a", "b"))
      sg <- ifelse(d$low == "a", 1e-3, 0.3)
      d$y <- sin(2 * d$x) + sg * s$e
      frm(bf(y ~ gp(x), sigma ~ low), data = d)
    }
  }))
}
# trends at sigma 0.3: the review's M4, its seeds
for (seed in c(5701L, 5702L)) {
  for (lab in c("0.2x", "x", "5x", "0.1x^2", "x^2", "exp(0.5x)", "exp(0.8x)",
                "logis")) {
    add(sprintf("trend %s s%d", lab, seed), local({
      lab <- lab
      seed <- seed
      function() {
        set.seed(seed)
        x <- round(runif(100, 0, 10), 3)
        e <- rnorm(100)
        f <- switch(lab, "0.2x" = 0.2 * x, x = x, "5x" = 5 * x,
                    "0.1x^2" = 0.1 * x^2, "x^2" = x^2,
                    "exp(0.5x)" = exp(0.5 * x), "exp(0.8x)" = exp(0.8 * x),
                    logis = 100 * plogis(x - 5))
        frm(bf(y ~ gp(x)), data = data.frame(x = x, y = f + 0.3 * e))
      }
    }), ref = "rel")
  }
}

# ------------------------------------------------------------ the reading
pos_info <- function(fit, bk) {
  o <- fit$obj
  env <- o$env
  par <- env$last.par.best
  H <- env$spHess(par, random = TRUE)
  ri <- which(names(par)[env$random] == "b")[bk$b_idx]
  dH <- Matrix::diag(H)[ri]
  th <- fit$estimates$theta[bk$theta_idx]
  V <- frmtmb:::cs_entry("gp")$vcov(th, bk)
  pmax(dH - diag(chol2inv(chol(V))), 0)
}

ids <- seq_along(des)
ids <- ids[(ids %% nch) == (chunk %% nch)]
tag <- Sys.getenv("GPGAP_CALIB_TAG")
if (nzchar(tag)) {
  ids <- which(vapply(des, function(d) grepl(tag, d$tag), NA))
}
for (i in ids) {
  t0 <- proc.time()[[3]]
  de <- des[[i]]
  warned <- FALSE
  fr <- tryCatch(withCallingHandlers(de$fn(), warning = function(w) {
    if (inherits(w, "frmtmb_gp_nugget")) warned <<- TRUE
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) conditionMessage(e))
  if (is.character(fr)) {
    cat(sprintf("CALIB %s | ERROR %s\n", de$tag, substr(fr, 1, 90)))
    next
  }
  bks <- Filter(function(b) identical(b$covstruct, "gp"), fr$frame$re_blocks)
  bk <- bks[[1]]
  sd_f <- exp(fr$estimates$theta[bk$theta_idx][[1]])
  info <- tryCatch(pos_info(fr, bk), error = function(e) rep(NA_real_, bk$dim))
  rr <- 1e-6 * sd_f^2 * info
  srt <- sort(rr, decreasing = TRUE)
  m <- length(rr)
  # the outer point of frm() and the reference, on the reference's SEs
  frn <- fr$frame
  frn$ncp_blocks <- which(vapply(frn$re_blocks,
                                 function(b) identical(b$covstruct, "gp"), NA))
  put(corr_at(de$ref))
  ref <- tryCatch({
    fn <- frmtmb::build_objective(frn)
    best <- NULL
    tpl0 <- frn$par_template
    tpl0$b[] <- 0
    starts <- list(tpl0)
    # also from frm()'s own point, so a stop short of the reference's
    # maximum shows rather than posing as an effect
    tpl1 <- tpl0
    for (nm in setdiff(names(tpl1), "b")) tpl1[[nm]] <- fr$estimates[[nm]]
    starts[[2]] <- tpl1
    tpl2 <- tpl1
    tpl2$theta[bk$theta_idx[1]] <- tpl2$theta[bk$theta_idx[1]] + 1
    starts[[3]] <- tpl2
    for (tp in starts) {
      # a start may fail where another does not: brms's absolute line does
      # not factor at some of them (review of 2026-11-02, M7)
      tryCatch({
        o <- RTMB::MakeADFun(fn, tp, random = "b", map = frn$map, silent = TRUE)
        op <- suppressWarnings(stats::nlminb(
          o$par, o$fn, o$gr, control = list(iter.max = 800, eval.max = 1200)))
        if (is.finite(op$objective) &&
              (is.null(best) || op$objective < best$op$objective)) {
          best <- list(op = op, o = o)
        }
      }, error = function(e) NULL)
    }
    if (is.null(best)) stop("no start of the reference finished")
    sdr <- tryCatch(RTMB::sdreport(best$o, best$op$par), error = function(e) NULL)
    se <- if (!is.null(sdr)) sqrt(pmax(diag(sdr$cov.fixed), 0)) else NULL
    list(ll = -best$op$objective, par = best$op$par, se = se,
         conv = best$op$convergence, gr = max(abs(best$o$gr(best$op$par))))
  }, error = function(e) conditionMessage(e))
  put(orig)
  if (!is.list(ref)) {
    cat(sprintf("CALIB %s | REF FAILED %s\n", de$tag, substr(ref, 1, 80)))
    next
  }
  pf <- fr$opt$par
  dz <- if (!is.null(ref$se) && length(ref$se) == length(pf)) {
    d <- abs(pf - ref$par) / ref$se
    d[!is.finite(d)] <- NA
    d
  } else {
    rep(NA_real_, length(pf))
  }
  th_i <- which(names(pf) == "theta")[bk$theta_idx]
  sg_i <- which(!(seq_along(pf) %in% th_i) & names(pf) != "theta" &
                  grepl("sigma|beta", names(pf)))
  # the reference's own one-parameter 95 percent boundary is 1.92 nats
  cat(sprintf(paste0(
    "CALIB %s | m %d rows %d conv %d warned %s | rmed %.4g rmean %.4g rq90 %.4g",
    " rmax %.4g scan10 %.4g | sd_f %.5g sd_r %.5g | ll_f %.4f ll_r %.4f gap %.4f",
    " | dzmax %.3f dzsd %.3f | refconv %d refgr %.1e | sec %.0f\n"),
    de$tag, m, fr$frame$n_obs %||% NA_integer_, fr$opt$convergence, warned,
    stats::median(rr), mean(rr), stats::quantile(rr, 0.9), max(rr),
    max(vapply(seq_len(m)[seq_len(m) >= max(3, m %/% 10)],
               function(k) srt[k] / sqrt(2 / k), 0)),
    sd_f, exp(ref$par[th_i[1]]), -fr$opt$objective, ref$ll,
    ref$ll + fr$opt$objective,
    if (all(is.na(dz))) NA_real_ else max(dz, na.rm = TRUE),
    dz[th_i[1]], ref$conv, ref$gr, proc.time()[[3]] - t0))
}
