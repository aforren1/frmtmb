# Lane ordstarts fix round (backlog 6136, review item 6): one refit of
# an ordinal mixture from its default fit A with the components'
# parameters permuted, for every assignment but the identity (K! - 1).
# Also the default fit's warnings at summary(), which the main run did
# not read (the standard-error check of a fit with random effects often
# waits for the first use).
#
# The permutation sigma gives new component k the parameters of old
# component sigma(k): its coefficients (mu<k>_*), its random-effect
# blocks (theta entries and modes, block by block), and its thresholds.
# The mixing weights are the multinomial logit against the last
# component, eta_K = 0, so the new eta_k is eta_sigma(k) - eta_sigma(K).
# Two threshold variants:
#   value: the thresholds' values move, and each component's raw values
#          are rebuilt from them (an ordered component, as cumulative is
#          here, stores its first threshold and the logs of the gaps;
#          values that are not increasing are sorted and given a gap of
#          at least 0.1 first). sratio stores the values themselves.
#   raw:   the raw vectors move by position, as the review's
#          ordstarts-review-swap.R did (between cumulative and sratio
#          that puts values into log-gap slots).
# dfn: the Laplace objective at the permuted start less A's, before any
# optimization: 0 up to rounding when the permutation is the same model.
# Each refit is refit(A, start = ) with starts = 0.
#
# Designs: the K = 2 designs of dev/ordstarts-run.R (same data), and
# k3: mixture(cumulative(), sratio(), acat()) on three-class data
# (class probabilities 0.4, 0.35, 0.25; latent slopes 1.5, -0.8, 0.3 on
# x and shifts 1, -1, 0; logistic noise; the same group effect and cut
# points as the K = 2 data), for which the default fit with
# starts = 10 is fitted too (arm B), since the main run has no K = 3.
#   Rscript dev/ordstarts-swap.R <family> <re> <n> <G> <seeds> <out.tsv>
#     [lean]
#     family: cum_cum | cum_sratio | k3. "lean" (used for k3): no arm B
#     and the value variant only, because on the pilot (seed 1) the
#     starts fit took 294.5 s and each permuted refit about 18 s, so the
#     full k3 design would not fit in an hour on two processes.
args <- commandArgs(trailingOnly = TRUE)
famname <- args[1]
re <- args[2]
nobs <- as.integer(args[3])
G <- as.integer(args[4])
seeds <- eval(parse(text = args[5]))
out <- args[6]
lean <- identical(args[7], "lean")
# k3 class separation: the class slopes and shifts times `sep` (8th
# argument, default 1; the converging K = 3 design of
# dev/ordstarts-k3probe.R is n = 3000, sep = 2, re = none)
sep <- if (length(args) > 7) as.numeric(args[8]) else 1
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
# the data of dev/ordstarts-run.R, read from it so the two cannot differ
run <- readLines("dev/ordstarts-run.R")
eval(parse(text = run[grep("^mkg <- ", run):(grep("^form <- ", run) - 1L)]))
mk3 <- function(seed) {
  set.seed(seed)
  x <- rnorm(nobs)
  z <- rnorm(nobs)
  cls <- sample.int(3L, nobs, replace = TRUE, prob = c(0.4, 0.35, 0.25))
  lat <- sep * (c(1.5, -0.8, 0.3)[cls] * x + c(1, -1, 0)[cls]) +
    rlogis(nobs)
  g <- rep_len(seq_len(G), nobs)
  if (re != "none") lat <- lat + rnorm(G, 0, 0.5)[g]
  if (re %in% c("slope", "slopeu")) lat <- lat + rnorm(G, 0, 0.3)[g] * x
  y <- 1L + (lat > -1 + 0.1 * z) + (lat > 0) + (lat > 1 - 0.1 * z)
  data.frame(x = x, z = z, y = y, g = factor(g))
}
# re = "none": no random effects, the data of dev/ordstart-starts.R
# (mk(seed, n, mix = TRUE), no group effect). Then A is the fit with
# starts = 0 and B today's default (10 starts), so the swap is measured
# against the default that ships without random effects.
none <- identical(re, "none")
form <- switch(re,
  none = bf(y ~ x + z),
  int = bf(y ~ x + z + (1 | g)),
  slope = bf(y ~ x + z + (1 + x | g)))
fam <- function() {
  switch(famname,
    cum_cum = mixture(cumulative(), cumulative()),
    cum_sratio = mixture(cumulative(), sratio()),
    k3 = mixture(cumulative(), sratio(), acat()))
}
q <- function(e) suppressWarnings(suppressMessages(e))
grab <- function(expr) {
  w <- character()
  m <- character()
  t0 <- proc.time()[["elapsed"]]
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage")
    }), error = function(e) e)
  list(v = v, w = w, m = m, el = proc.time()[["elapsed"]] - t0)
}
clean <- function(s) paste(substr(gsub("[\t\n\r]", " ", s), 1, 100),
                           collapse = " | ")
perms <- function(K) {
  if (K == 1L) return(list(1L))
  out <- list()
  for (i in seq_len(K)) {
    for (r in perms(K - 1L)) out[[length(out) + 1L]] <- c(i, setdiff(
      seq_len(K), i)[r])
  }
  out
}
# the permuted start for sigma, as frm()'s `start` list
perm_start <- function(f, sigma, variant) {
  resp <- names(f$spec$responses)[1]
  K <- length(sigma)
  p <- f$opt$par
  nm <- ns$outer_par_map(f)$names
  stopifnot(length(nm) == length(p))
  s <- p
  # coefficients
  for (k in seq_len(K)) {
    to <- grep(paste0("^mu", k, "_"), nm)
    from <- match(sub(paste0("^mu", k, "_"), paste0("mu", sigma[k], "_"),
                      nm[to]), nm)
    s[to] <- p[from]
  }
  # weights against the last component
  eta <- c(vapply(seq_len(K - 1L), function(k) {
    p[match(paste0("theta", k, "_(Intercept)"), nm)]
  }, 0), 0)
  for (k in seq_len(K - 1L)) {
    s[match(paste0("theta", k, "_(Intercept)"), nm)] <-
      eta[sigma[k]] - eta[sigma[K]]
  }
  # thresholds
  obs <- lapply(seq_len(K), function(k) {
    lp <- f$frame$linpreds[[ns$linpred_key(resp, paste0("mu", k))]]
    ns$ord_lp_block(f$frame, f$spec, lp)
  })
  rawk <- lapply(obs, function(ob) as.numeric(f$estimates[[ob$comp]]))
  valk <- lapply(seq_len(K), function(k) {
    as.numeric(ns$ord_threshold_values(obs[[k]]$fam, rawk[[k]]))
  })
  for (k in seq_len(K)) {
    pos <- which(names(p) == obs[[k]]$comp)
    src <- sigma[k]
    newraw <- if (variant == "raw") rawk[[src]] else {
      v <- valk[[src]]
      if (isTRUE(obs[[k]]$ordered)) {
        v <- sort(v)
        for (j in seq_along(v)[-1L]) v[j] <- max(v[j], v[j - 1L] + 0.1)
        c(v[1L], log(diff(v)))
      } else v
    }
    stopifnot(length(pos) == length(newraw))
    s[pos] <- newraw
  }
  # random-effect blocks: theta entries and modes, block for block
  thpos <- which(names(p) == "theta")
  blks <- f$frame$re_blocks
  bd <- vapply(blks, function(b) b$dpar %||% "", "")
  bnew <- f$estimates$b
  for (i in seq_along(blks)) {
    k <- as.integer(sub("^mu", "", bd[i]))
    if (is.na(k)) next
    j <- which(bd == paste0("mu", sigma[k]) &
                 vapply(blks, function(b) {
                   identical(b$cnms, blks[[i]]$cnms) &&
                     identical(b$group_name, blks[[i]]$group_name)
                 }, NA))
    stopifnot(length(j) == 1L)
    s[thpos[blks[[i]]$theta_idx]] <- p[thpos[blks[[j]]$theta_idx]]
    bnew[blks[[i]]$b_idx] <- f$estimates$b[blks[[j]]$b_idx]
  }
  dfn <- tryCatch(q(as.numeric(f$obj$fn(s))) - f$opt$objective,
                  error = function(e) NA_real_)
  st <- f$obj$env$parList(s)
  st[["b"]][] <- bnew
  list(start = st, dfn = dfn)
}
row_of <- function(r, what, variant = "", dfn = NA) {
  f <- r$v
  base <- data.frame(family = famname, re = re, n = nobs, G = G,
                     seed = seed, what = what, variant = variant,
                     elapsed = r$el, dfn = dfn)
  if (inherits(f, "error")) {
    base$err <- clean(conditionMessage(f))
    return(base)
  }
  resp <- names(f$spec$responses)[1]
  base$err <- ""
  base$ll <- as.numeric(logLik(f))
  base$code <- f$opt$convergence
  base$ms_starts <- f$opt$mixture_starts[["starts"]] %||% 0
  cl <- tryCatch(ns$mixture_ord_start_class(f, resp), error = function(e) e)
  base$cls <- if (inherits(cl, "error")) "error" else
    if (length(cl$step)) "step" else if (length(cl$boundary)) "flat" else
      "proper"
  base$degW <- any(grepl("ended at a degenerate boundary", r$w))
  base$lostW <- any(grepl("Standard errors are not available", r$w))
  base$sdzero <- any(grepl("Boundary [(]singular[)] fit", c(r$w, r$m)))
  base$nwarn <- length(r$w)
  base$deferred <- is.null(f$cache$se_fit_an) && is.null(f$cache$sdr)
  q(ns$se_fit_analysis(f, force = TRUE))
  tse <- tryCatch(q(ns$mixture_ord_thres_se(f, resp)),
                  error = function(e) NA)
  base$tse <- suppressWarnings(max(tse, na.rm = TRUE))
  base$warns <- clean(r$w)
  base$msgs <- clean(r$m)
  base
}
cols <- c("family", "re", "n", "G", "seed", "what", "variant", "elapsed",
          "dfn", "ll", "code", "ms_starts", "cls", "degW", "lostW", "sdzero",
          "nwarn", "deferred", "tse", "s_nwarn", "s_lostW", "s_degW",
          "s_sdzero", "err", "warns", "msgs", "s_warns")
first <- !file.exists(out)
emit <- function(row) {
  for (cn in setdiff(cols, names(row))) row[[cn]] <- NA
  utils::write.table(row[cols], out, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = first, append = !first)
  first <<- FALSE
}
K <- if (famname == "k3") 3L else 2L
P <- Filter(function(s) !identical(s, seq_len(K)), perms(K))
for (seed in seeds) {
  d <- if (K == 3L) mk3(seed) else if (none) {
    mk(seed, n = nobs, mix = TRUE)
  } else mkg(seed)
  a <- grab(frm(form, family = fam(), data = d,
                control = if (none) frmtmb_control(starts = 0) else
                  frmtmb_control()))
  row <- row_of(a, "A")
  if (!inherits(a$v, "error")) {
    A <- a$v
    sm <- grab(summary(A))
    row$s_nwarn <- length(sm$w)
    row$s_lostW <- any(grepl("Standard errors are not available", sm$w))
    row$s_degW <- any(grepl("ended at a degenerate boundary", sm$w))
    row$s_sdzero <- any(grepl("Boundary [(]singular[)] fit",
                              c(sm$w, sm$m)))
    row$s_warns <- clean(sm$w)
  }
  emit(row)
  if ((K == 3L && !lean) || none) {
    b <- grab(frm(form, family = fam(), data = d,
                  control = frmtmb_control(starts = 10)))
    emit(row_of(b, "B"))
  }
  if (inherits(a$v, "error")) next
  A0 <- A
  A0$control$starts <- 0L
  for (variant in if (lean) "value" else c("value", "raw")) {
    for (sg in P) {
      lab <- paste(sg, collapse = "")
      ps <- tryCatch(perm_start(A, sg, variant), error = function(e) e)
      if (inherits(ps, "error")) {
        emit(data.frame(family = famname, re = re, n = nobs, G = G,
                        seed = seed, what = paste0("perm", lab),
                        variant = variant, elapsed = NA, dfn = NA,
                        err = clean(conditionMessage(ps))))
        next
      }
      r <- grab(refit(A0, d$y, start = ps$start))
      emit(row_of(r, paste0("perm", lab), variant, ps$dfn))
    }
  }
  cat(famname, re, nobs, "seed", seed, "done\n")
}
cat("DONE swap\n")
