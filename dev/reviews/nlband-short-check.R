# Short check of lane nlband's fix round 1 (round of 2026-11-16).
# Sections: repro (separated bernoulli seeds 1 to 10, student seeds 1 to
# 3), band (seeds 1 to 3, 201 to 205 against the exact profile maximum),
# gate (nl_hole_located() on narrow-hole, noisy, n = 60 designs, and a
# supremum in a coefficient other than the hole's), cost (evaluations
# and the cost of nl_hole_located() against an emulation of the first
# build: nl_hole_located() TRUE and nl_cold_par() NULL).
# NLBAND_ARM=before|after Rscript dev/reviews/nlband-short-check.R <secs>
args <- commandArgs(trailingOnly = TRUE)
secs <- if (length(args)) strsplit(args[1], ",")[[1]] else
  c("repro", "band", "gate", "cost")
source("dev/nlband-common.R")
`%||%` <- function(a, b) if (is.null(a)) b else a
if (arm == "after") {
  # the installed build is the worktree's source
  src <- new.env()
  sys.source("R/fit-restart.R", envir = src, keep.source = FALSE)
  for (f in c("nl_hole_located", "nl_cold_par", "nl_pole_escape")) {
    cat("SRC", f, "identical to installed:",
        identical(deparse(body(src[[f]])), deparse(body(ns[[f]]))), "\n")
  }
}
q <- function(expr) suppressWarnings(suppressMessages(expr))
prof1 <- function(x, y, b) {
  vapply(b, function(bb) {
    z <- 1 / (x + bb)
    a <- sum(y * z) / sum(z * z)
    rss <- sum((y - a * z)^2)
    n <- length(y)
    -n / 2 * (log(2 * pi * rss / n) + 1)
  }, 1)
}
# exact profile maximum over b of y ~ a / (x + b) + N(0, s): 15 points in
# every gap between walls, 60 beyond each end, then optimize()
glob1 <- function(x, y) {
  w <- sort(unique(-x))
  r <- diff(range(w))
  u <- seq(0, 1, length.out = 17)[2:16]
  b <- c(as.vector(outer(u, diff(w)) + rep(w[-length(w)], each = 15)),
         w[1] - r * 10^seq(-4, 2, length.out = 60),
         w[length(w)] + r * 10^seq(-4, 2, length.out = 60))
  b <- sort(b)
  p <- prof1(x, y, b)
  p[!is.finite(p)] <- -Inf
  k <- which.max(p)
  o <- optimize(function(bb) prof1(x, y, bb), c(b[max(1, k - 1)],
                b[min(length(b), k + 1)]), maximum = TRUE, tol = 1e-12)
  c(max(o$objective, p[k]), if (o$objective >= p[k]) o$maximum else b[k])
}
cls <- function(r, k) {
  v <- attr(r, k)
  if (length(v)) paste(unique(v), collapse = "+") else "-"
}

if ("repro" %in% secs) {
  for (s in 1:10) {
    d <- band_data(s)
    d$y <- as.integer(d$x > 2)
    r <- tt(frm(band_form, data = d, family = bernoulli()))
    cat(sprintf("SEP seed %d | ll %.8g code %d beta %s | warn %s msg %s\n",
                s, llv(r), r$opt$convergence,
                paste(sprintf("%.6g", r$estimates$beta), collapse = " "),
                cls(r, "wc"), cls(r, "mc")))
  }
  for (s in 1:3) {
    set.seed(s)
    x <- band_data(s)$x
    d <- data.frame(x = x, y = 1 / (x - 2) + 0.2 * rt(length(x), 3))
    r <- tt(frm(band_form, data = d, family = student()))
    g <- tt(frm(band_form, data = d, family = student(),
                start = list(beta = c(1, -2)),
                control = frmtmb_control(nl_starts = FALSE)))
    # the maximum: the best of the start in the hole and 0.83.0's own
    cat(sprintf(paste0("STUD seed %d | default ll %.6f code %d par %s | ",
                       "warn %s msg %s | from (1,-2) ll %.6f | short %.6f\n"),
                s, llv(r), r$opt$convergence,
                paste(sprintf("%.6g", r$opt$par), collapse = " "),
                cls(r, "wc"), cls(r, "mc"), llv(g), llv(g) - llv(r)))
    for (m in attr(r, "m")) cat("   MSG", gsub("\n", " ", m), "\n")
    for (w in attr(r, "w")) cat("   WARN", substr(w, 1, 300), "\n")
  }
}

if ("band" %in% secs) {
  for (s in c(1:3, 201:205)) {
    d <- band_data(s)
    r <- tt(frm(band_form, data = d))
    g <- glob1(d$x, d$y)
    cat(sprintf(paste0("BAND seed %d | ll %.6f code %d b %.5f | warn %s msg %s",
                       " | GLOB %.6f at %.5f | short %.2e\n"),
                s, llv(r), r$opt$convergence, r$estimates$beta[2],
                cls(r, "wc"), cls(r, "mc"), g[1], g[2], g[1] - llv(r)))
  }
}

# a record of every nl_hole_located() call: the SE, the half width, the
# verdict, the refit's logLik
loc <- new.env()
loc$rows <- list()
trace_loc <- function() {
  suppressMessages(trace("nl_hole_located", where = ns, print = FALSE,
    tracer = quote({
      .t0 <- proc.time()[["elapsed"]]
      .cached <- is.environment(fit$cache) && !is.null(fit$cache$se_fit_an)
    }),
    exit = quote({
      .row <- list(se = if (exists("s", inherits = FALSE)) s else NA,
                   half = (ck$hi - ck$lo) / 2, lo = ck$lo, hi = ck$hi,
                   ok = returnValue(),
                   ll = tryCatch(as.numeric(logLik(fit)), error = function(e) NA),
                   b = fit$estimates$beta[2], code = fit$opt$convergence,
                   cached = .cached,
                   t = proc.time()[["elapsed"]] - .t0)
      assign("rows", c(get("rows", envir = loc), list(.row)), envir = loc)
    })))
}
show_loc <- function(tag) {
  for (r in loc$rows) {
    cat(sprintf(paste0("   LOC %s | se %.4g half %.4g hole (%.4f, %.4f) -> %s",
                       " | refit ll %.6f b %.5f code %d | cached %s t %.3f\n"),
                tag, r$se, r$half, r$lo, r$hi, r$ok, r$ll, r$b, r$code,
                r$cached, r$t))
  }
  loc$rows <- list()
}

if ("gate" %in% secs && arm == "after") {
  trace_loc()
  n <- 60
  # the narrowest hole about the pole at 2 that nl_pole_holes() flags on
  # n rows of U(0, 4)
  narrow <- function(seed) {
    for (h in seq(0.2, 1.6, by = 0.02)) {
      set.seed(seed)
      x <- runif(2000, 0, 4)
      x <- x[abs(x - 2) > h][seq_len(n)]
      hs <- ns$nl_pole_holes(sort(unique(-x)))
      if (!is.null(hs) && any(hs[, "lo"] < -2 & hs[, "hi"] > -2)) {
        return(list(x = x, h = h))
      }
    }
    NULL
  }
  for (amp in c(1, 0.5, 0.3, 0.15)) for (sdv in c(1, 0.2)) for (s in 1:4) {
    nw <- narrow(s)
    if (is.null(nw)) { cat("GATE seed", s, "no hole flagged\n"); next }
    set.seed(500 + s)
    d <- data.frame(x = nw$x)
    d$y <- amp / (d$x - 2) + rnorm(n, 0, sdv)
    g <- glob1(d$x, d$y)
    r <- tt(frm(band_form, data = d))
    ingap <- g[2] > -2 - nw$h && g[2] < -2 + nw$h
    cat(sprintf(paste0("GATE amp %.2f sd %.1f seed %d h %.2f | ll %.6f code %d",
                       " b %.5f | warn %s msg %s | GLOB %.6f at %.5f in hole %s",
                       " | short %.4f\n"),
                amp, sdv, s, nw$h, llv(r), r$opt$convergence,
                r$estimates$beta[2], cls(r, "wc"), cls(r, "mc"), g[1], g[2],
                ingap, g[1] - llv(r)))
    show_loc(sprintf("amp %.2f sd %.1f seed %d", amp, sdv, s))
  }
  # a supremum in another coefficient: bernoulli, z separates (y = 1
  # wherever z = 1), the hole's b located by the rows with z = 0
  form2 <- bf(y ~ a / (x + b) + c * z, a + b + c ~ 1, nl = TRUE)
  for (s in 1:5) {
    d <- band_data(s)
    set.seed(900 + s)
    d$z <- rbinom(nrow(d), 1, 0.2)
    d$y <- rbinom(nrow(d), 1, plogis(1 / (d$x - 2)))
    d$y[d$z == 1] <- 1L
    r <- tt(frm(form2, data = d, family = bernoulli()))
    cat(sprintf(paste0("SUPZ seed %d | ll %.6f code %d beta %s | warn %s ",
                       "msg %s\n"), s, llv(r), r$opt$convergence,
                paste(sprintf("%.5g", r$estimates$beta), collapse = " "),
                cls(r, "wc"), cls(r, "mc")))
    for (w in attr(r, "w")) cat("   WARN", substr(w, 1, 200), "\n")
    show_loc(sprintf("supz seed %d", s))
  }
  untrace("nl_hole_located", where = ns)
}

if ("cost" %in% secs && arm == "after") {
  orig <- list(nl_hole_located = ns$nl_hole_located,
               nl_cold_par = ns$nl_cold_par)
  put <- function(nm, f) {
    unlockBinding(nm, ns)
    assign(nm, f, envir = ns)
    lockBinding(nm, ns)
  }
  first_build <- function(on) {
    if (on) {
      f1 <- function(fit, ck) TRUE
      f2 <- function(frame, pn) NULL
      environment(f1) <- ns
      environment(f2) <- ns
      put("nl_hole_located", f1)
      put("nl_cold_par", f2)
    } else {
      put("nl_hole_located", orig$nl_hole_located)
      put("nl_cold_par", orig$nl_cold_par)
    }
  }
  n_ev <- 0
  suppressMessages(trace("optimize_obj",
                         exit = quote(n_ev <<- n_ev + returnValue()$evals),
                         print = FALSE, where = ns))
  n_he <- 0
  suppressMessages(trace("obj_he_cached",
                         tracer = quote(n_he <<- n_he + 1),
                         print = FALSE, where = ns))
  worse <- function(s) {
    d <- band_data(s)
    d$y <- 1 / (d$x + 0.5) + rnorm(nrow(d), 0, 0.2)
    d
  }
  re <- function(s) {
    d <- band_data(s)
    d$g <- factor(rep_len(1:8, nrow(d)))
    set.seed(300 + s)
    u <- rnorm(8, 0, 0.1)
    d$y <- (1 + u[d$g]) / (d$x - 2) + rnorm(nrow(d), 0, 0.2)
    d
  }
  sets <- list(
    list("hole wins (band data)", band_form, lapply(1:3, band_data)),
    list("hole loses (truth b 0.5)", band_form, lapply(1:3, worse)),
    list("hole wins, a ~ 1 + (1 | g)",
         bf(y ~ a / (x + b), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
         lapply(1:2, re)))
  for (st in sets) for (i in seq_along(st[[3]])) {
    d <- st[[3]][[i]]
    out <- list()
    for (a in c("fixed", "first")) {
      first_build(a == "first")
      if (a == "fixed") trace_loc()
      n_ev <- 0; n_he <- 0
      t0 <- proc.time()[["elapsed"]]
      r <- q(frm(st[[2]], data = d))
      el <- proc.time()[["elapsed"]] - t0
      out[[a]] <- c(ev = n_ev, he = n_he, t = el, ll = llv(r))
      if (a == "fixed") {
        lr <- loc$rows
        untrace("nl_hole_located", where = ns)
      }
    }
    first_build(FALSE)
    cat(sprintf(paste0("COST %s row %d | fixed ev %d he %d t %.2f ll %.6f | ",
                       "first ev %d he %d t %.2f ll %.6f | ev ratio %.3f\n"),
                st[[1]], i, as.integer(out$fixed[["ev"]]),
                as.integer(out$fixed[["he"]]), out$fixed[["t"]],
                out$fixed[["ll"]], as.integer(out$first[["ev"]]),
                as.integer(out$first[["he"]]), out$first[["t"]],
                out$first[["ll"]], out$fixed[["ev"]] / out$first[["ev"]]))
    loc$rows <- lr
    show_loc(sprintf("%s row %d", st[[1]], i))
  }
}
cat("DONE\n")
