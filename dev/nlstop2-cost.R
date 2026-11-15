# Lane nlstop2: the cost of the lane's three changes, in one process,
# arms interleaved, minimum of 3 rounds, each with a control arm built
# from the same code that must read 1.
#   pole:  nl_pole_escape() on (the lane) against replaced by NULL (off)
#          and off again (off2, the control). Blocks: "clear", the
#          fresh mm rows whose K is clear of its poles (gaussian, poisson,
#          gamma, cumulative, seeds 11 to 55); "pressed", binomial mm 11,
#          33 and 44, where the escape runs.
#   limit: nl_limit_message() on against the 0.82.0 copy (one move, no
#          step across) on the fits of the scan's 45 bodies at seed 101,
#          and the copy again (control); the call alone, repeated 100 times.
#   huber: y ~ x + (1 | g), the default (quadrature) against
#          quadrature = FALSE, and FALSE again (control), on 100 rows
#          (20 groups, seeds 1 to 10) and on 2000 rows (100 groups, seeds
#          1 to 3).
# Usage: Rscript dev/nlstop2-cost.R <pole|limit|huber>
source("dev/nlstop2-common.R")
fresh_defs()
ns <- asNamespace("frmtmb")
what <- commandArgs(TRUE)[1]
tm <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}
q <- function(expr) suppressWarnings(suppressMessages(expr))
rounds <- function(arms, block, n = 3) {
  res <- matrix(NA, n, length(arms), dimnames = list(NULL, names(arms)))
  for (i in seq_len(n)) {
    ord <- if (i %% 2) seq_along(arms) else rev(seq_along(arms))
    for (k in ord) {
      arms[[k]]()
      res[i, k] <- tm(block())
    }
  }
  print(res)
  m <- apply(res, 2, min)
  m
}
set_fun <- function(nm, f) {
  unlockBinding(nm, ns)
  assign(nm, f, envir = ns)
  lockBinding(nm, ns)
}
if (what == "pole") {
  orig <- ns$nl_pole_escape
  arms <- list(on = function() set_fun("nl_pole_escape", orig),
               off = function() set_fun("nl_pole_escape", function(...) NULL),
               off2 = function() set_fun("nl_pole_escape",
                                         function(...) NULL))
  clear <- list()
  for (fam in c("gaussian", "poisson", "gamma", "cumulative")) {
    for (seed in c(11, 22, 33, 44, 55)) {
      clear[[length(clear) + 1L]] <- list(
        fam = fam, d = gen(fam, "mm", seed * 1000 + 1), f = form(fam, "mm"))
    }
  }
  pressed <- lapply(c(11, 33, 44), function(seed) {
    list(fam = "binomial", d = gen("binomial", "mm", seed * 1000 + 1),
         f = form("binomial", "mm"))
  })
  run_block <- function(rows, reps) function() {
    for (r in seq_len(reps)) for (x in rows) {
      q(frm(x$f, family = F[[x$fam]], data = x$d))
    }
  }
  m <- rounds(arms, run_block(clear, 2))
  cat(sprintf("pole clear: min s on %.2f off %.2f off2 %.2f | on/off %.3f control %.3f\n",
              m["on"], m["off"], m["off2"], m["on"] / m["off"],
              m["off2"] / m["off"]))
  m <- rounds(arms, run_block(pressed, 2))
  cat(sprintf("pole pressed: min s on %.2f off %.2f off2 %.2f | on/off %.3f control %.3f\n",
              m["on"], m["off"], m["off2"], m["on"] / m["off"],
              m["off2"] / m["off"]))
  set_fun("nl_pole_escape", orig)
}
if (what == "polecount") {
  # load-independent: the objective and gradient evaluations of every
  # optimize_obj() run of a fit (opt$evals), escape on and off
  n_ev <- 0
  trace("optimize_obj", exit = quote(n_ev <<- n_ev + returnValue()$evals),
        print = FALSE, where = ns)
  orig <- ns$nl_pole_escape
  rows <- c(lapply(c(11, 33, 44), function(s) c("binomial", s)),
            lapply(c(11, 22, 33, 44, 55), function(s) c("gaussian", s)),
            lapply(c(11, 22, 33, 44, 55), function(s) c("poisson", s)))
  for (x in rows) {
    d <- gen(x[1], "mm", as.numeric(x[2]) * 1000 + 1)
    ev <- c(on = NA, off = NA)
    ts <- ev
    for (arm in names(ev)) {
      set_fun("nl_pole_escape",
              if (arm == "on") orig else function(...) NULL)
      n_ev <- 0
      ts[arm] <- tm(q(frm(form(x[1], "mm"), family = F[[x[1]]], data = d)))
      ev[arm] <- n_ev
    }
    cat(sprintf(paste0("polecount %s mm %s: evaluations on %d off %d ",
                       "ratio %.3f | s on %.2f off %.2f\n"),
                x[1], x[2], ev[["on"]], ev[["off"]], ev[["on"]] / ev[["off"]],
                ts[["on"]], ts[["off"]]))
  }
  set_fun("nl_pole_escape", orig)
}
if (what == "limit") {
  new <- ns$nl_limit_message
  old <- local({
    e <- new.env(parent = ns)
    Sys.setenv(GIT_OPTIONAL_LOCKS = "0")
    src <- system(paste("git -C C:/Users/adf44/source/r/frmtmb-wt-nlstop2",
                        "show ec8178d3:R/fit-restart.R"), intern = TRUE)
    eval(parse(text = src), e)
    e$nl_limit_message
  })
  bodies <- list(
    exp = quote(b1 * exp(b2 * x)), sat = quote(b1 * (1 - exp(-b2 * x))),
    logi = quote(b1 / (1 + exp(b2 * x))), prod = quote(a + b1 * b2 * x),
    two = quote(b1 * exp(b2 * x) + b3 * z))
  fits <- list()
  tfit <- 0
  for (fam in c("gaussian", "poisson", "bernoulli")) {
    for (bn in names(bodies)) {
      set.seed(101)
      d <- data.frame(x = round(rnorm(300), 3), z = round(rnorm(300), 3))
      lat <- 1.5 * exp(0.6 * d$x) + 0.3 * d$z
      d$y <- switch(fam, gaussian = lat + rnorm(300, 0, 0.5),
                    poisson = rpois(300, lat),
                    bernoulli = rbinom(300, 1, plogis(lat - 2)))
      p <- setdiff(all.vars(bodies[[bn]]), c("x", "z"))
      fm <- do.call(bf, c(list(as.formula(paste("y ~",
                                                deparse1(bodies[[bn]])))),
                          lapply(p, function(v) as.formula(paste(v, "~ 1"))),
                          list(nl = TRUE)))
      t0 <- proc.time()[["elapsed"]]
      fit <- tryCatch(q(frm(fm, family = get(fam)(), data = d)),
                      error = function(e) NULL)
      tfit <- tfit + proc.time()[["elapsed"]] - t0
      if (!is.null(fit)) fits[[length(fits) + 1L]] <- fit
    }
  }
  cat("fits:", length(fits), "\n")
  cur <- new
  arms <- list(on = function() cur <<- new, off = function() cur <<- old,
               off2 = function() cur <<- old)
  block <- function() for (r in 1:100) for (f in fits) q(cur(f))
  m <- rounds(arms, block)
  cat(sprintf("limit: min s on %.2f off %.2f off2 %.2f | on/off %.3f control %.3f\n",
              m["on"], m["off"], m["off2"], m["on"] / m["off"],
              m["off2"] / m["off"]))
  cat(sprintf("limit: per call on %.4f s off %.4f s against one fit's %.4f s\n",
              m["on"] / (100 * length(fits)), m["off"] / (100 * length(fits)),
              tfit / length(fits)))
}
if (what == "huber") {
  mk <- function(seed, ng, per) {
    set.seed(seed)
    d <- data.frame(g = factor(rep(seq_len(ng), each = per)),
                    x = rnorm(ng * per))
    d$y <- 1 + 0.5 * d$x + 0.5 * rnorm(ng)[d$g] + rt(ng * per, 3)
    d
  }
  for (sz in list(c(20, 5, 10), c(100, 20, 3))) {
    ds <- lapply(seq_len(sz[3]), mk, ng = sz[1], per = sz[2])
    mode <- NULL
    arms <- list(on = function() mode <<- "default",
                 off = function() mode <<- "laplace",
                 off2 = function() mode <<- "laplace")
    block <- function() for (d in ds) {
      if (mode == "default") {
        q(frm(y ~ x + (1 | g), family = huber(), data = d))
      } else {
        q(frm(y ~ x + (1 | g), family = huber(), data = d,
              quadrature = FALSE))
      }
    }
    m <- rounds(arms, block)
    cat(sprintf("huber %d rows x %d seeds: min s default %.2f laplace %.2f laplace2 %.2f | default/laplace %.2f control %.3f\n",
                sz[1] * sz[2], sz[3], m["on"], m["off"], m["off2"],
                m["on"] / m["off"], m["off2"] / m["off"]))
  }
}
cat("DONE\n")
