# Short check 2, items 1 and 3. Per dataset, on the final build: the fit, the
# margin (with the post-fit refusal off in-process, so the fit is on the record),
# the warning, and the first-order Laplace correction read directly: the
# Tierney-Kadane term C of the latents, and the shift it implies in the outer
# parameters, V grad(C).
#   log Z = Laplace + C,  with h = -log density of a row in eta,
#   C = -1/8 sum h4_i v_ii^2 + 1/8 sum h3_i h3_j v_ii v_ij v_jj + 1/12 sum h3_i h3_j v_ij^3
# v = Cov(err | y) at the inner mode = A H^-1 A'.
# Usage: Rscript review-latres-tk.R <set> <design> <s1> <s2> <outfile>
root <- "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/"
args0 <- commandArgs(TRUE)
set <- args0[[1]]; design <- args0[[2]]; s1 <- as.integer(args0[[3]]); s2 <- as.integer(args0[[4]]); outf <- args0[[5]]
source(paste0(root, "review-latres-common.R"))
options(warn = -1)
commandArgs <- function(trailingOnly = FALSE) c(lane_lib, design, as.character(s1), as.character(s2), "x")
envd <- new.env(parent = globalenv())
file_of <- switch(set, lane = "latres-m1.R", m1b = "review-latres-m1b.R", fresh = "review-latres-fresh.R", m4b = "review-latres-m4b.R")
src <- readLines(paste0(root, file_of))
cut_at <- switch(set, lane = grep("^row_dens <- ", src), grep("^one <- ", src))
if (set != "lane") {
  commandArgs <- function(trailingOnly = FALSE) c(design, as.character(s1), as.character(s2), "x")
}
eval(parse(text = src[seq_len(cut_at - 1L)]), envd)
if (set == "lane") {
  famx <- envd$famf
  lhs <- if (envd$D$fam == "binomial") "y | trials(tr)" else "y"
  is_cum <- envd$D$fam == "cumulative"
  str_true <- envd$D$s
  ex_file <- paste0(root, "latres-out/m1/", design, ".tsv")
  ex <- read.table(ex_file, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE)
  ex$x_ls <- log(ex$x_sderr)
} else {
  famx <- envd$fam
  dist <- envd$D$dist
  lhs <- if (dist %in% c("binom", "bb", "zib")) "y | trials(tr)" else "y"
  is_cum <- identical(dist, "cum")
  str_true <- envd$D$s
  dir <- switch(set, m1b = "m1b", fresh = "fresh", m4b = "m4b2")
  fs <- list.files(paste0(root, "review-latres-out/", dir), paste0("^", design, "(-[0-9]+)?[.]tsv$"), full.names = TRUE)
  ex <- do.call(rbind, lapply(fs, function(f) read.table(f, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE)))
}
ns <- asNamespace("frmtmb")
real_margin <- get("latres_margin", ns)
assignInNamespace("latres_refuse_margin", function(fit) invisible(NULL), ns = "frmtmb")
form <- as.formula(paste(lhs, "~ x + ar(t, g, cov = TRUE)"))

row_d <- function(fit, bk) {
  a <- bk$aux
  r <- a$ac$resp
  fam <- fit$spec$responses[[r]]$family
  dp <- fi("eval_dpars")(fit)[[r]]
  lk <- fit$frame$linpreds[[fi("linpred_key")(r, "mu")]]$link
  eta <- lk$linkfun(as.numeric(dp$mu))
  y <- fit$frame$y[[r]]
  av <- fit$frame$aterm_values[[r]]
  exx <- fi("fit_extras")(fit)
  ll <- function(e) {
    dp$mu <- lk$linkinv(e)
    as.numeric(fi("row_lpdf")(fam, y, y, dp, av, exx))
  }
  sds <- fi("latres_marginal_sd")(fit$estimates$theta[bk$theta_idx], a)
  h <- 0.02 * sds
  l0 <- ll(eta); lp <- ll(eta + h); ln <- ll(eta - h); lp2 <- ll(eta + 2 * h); ln2 <- ll(eta - 2 * h)
  list(h2 = -(lp - 2 * l0 + ln) / h^2,
       h3 = -(lp2 - 2 * lp + 2 * ln - ln2) / (2 * h^3),
       h4 = -(lp2 - 4 * lp + 6 * l0 - 4 * ln + ln2) / h^4)
}
map_of <- function(fit) fi("outer_par_map")(fit)
tk_at <- function(fit, par, bk, map) {
  env <- fit$obj$env
  v <- suppressWarnings(fit$obj$fn(par))
  if (!is.finite(v)) return(NA_real_)
  full <- env$last.par
  f2 <- fi("fit_set_outer")(fit, par, map)
  f2$estimates$b <- as.numeric(full[env$random])
  d <- row_d(f2, bk)
  H <- as.matrix(env$spHess(full, random = TRUE))
  pj <- fi("latres_jacobian")(bk, f2$estimates$theta[bk$theta_idx])
  A <- matrix(0, bk$aux$n, length(env$random))
  A[cbind(pj$i, pj$j)] <- pj$x
  Vv <- A %*% solve(H, t(A))
  di <- diag(Vv)
  -sum(d$h4 * di^2) / 8 + sum(outer(d$h3 * di, d$h3 * di) * Vv) / 8 + sum(outer(d$h3, d$h3) * Vv^3) / 12
}
one <- function(seed) {
  d <- envd$simd(seed)
  w <- character(0); wc <- character(0)
  f <- withCallingHandlers(tryCatch(frm(form, data = d, family = famx), error = function(e) e),
    warning = function(c) { w <<- c(w, conditionMessage(c)); wc <<- c(wc, class(c)[1]); invokeRestart("muffleWarning") },
    message = function(c) invokeRestart("muffleMessage"))
  out <- list(set = set, design = design, seed = seed)
  if (inherits(f, "error")) { out$err <- substr(conditionMessage(f), 1, 80); return(out) }
  out$err <- ""
  bk <- fi("frame_latres")(f$frame)[[1]]
  th <- f$estimates$theta[bk$theta_idx]
  par <- f$opt$par
  nm <- names(par)
  out$ll <- as.numeric(logLik(f)); out$code <- f$opt$convergence
  out$ls <- th[1]; out$t <- th[2]; out$ar <- th[2] / sqrt(1 + th[2]^2)
  out$warned <- as.integer(any(grepl("Laplace approximation of", w, fixed = TRUE)))
  out$why <- substr(paste(w[grepl("Laplace approximation of", w, fixed = TRUE)], collapse = ""), 60, 120)
  out$boxed <- as.integer(any(grepl("ran to log", w, fixed = TRUE)))
  m <- tryCatch(real_margin(f, bk), error = function(e) NULL)
  out$margin <- if (is.null(m)) NA else m
  V <- tryCatch(as.matrix(vcov(f, full = TRUE)), error = function(e) NULL)
  if (!is.null(V)) {
    out$se_ls <- sqrt(V[nrow(V) - 1L, nrow(V) - 1L])
    out$se_t <- sqrt(V[nrow(V), nrow(V)])
  }
  map <- map_of(f)
  t0 <- proc.time()[[3]]
  C0 <- tryCatch(tk_at(f, par, bk, map), error = function(e) NA_real_)
  out$C <- C0
  out$n <- bk$aux$n
  # the correction away from the fit: at twice and half sderr, and across the
  # fit's own interval, with the others held
  i_ls <- length(par) - 1L
  at_ls <- function(dl) {
    pp <- par
    pp[i_ls] <- pp[i_ls] + dl
    tryCatch(tk_at(f, pp, bk, map), error = function(e) NA_real_)
  }
  out$C_up <- at_ls(log(2))
  out$C_dn <- at_ls(-log(2))
  se0 <- out$se_ls
  dl <- if (is.finite(se0 %||% NA)) min(1.96 * se0, 2) else log(2)
  out$C_chi <- at_ls(dl)
  out$C_clo <- at_ls(-dl)
  g <- rep(NA_real_, length(par))
  # theta and the dispersions alone, under REV_TK_RESTRICT: beta's gradient is small
  # and each entry costs two inner solves
  ks <- if (nzchar(Sys.getenv("REV_TK_RESTRICT"))) which(nm != "beta") else seq_along(par)
  for (k in ks) {
    hh <- 0.05
    pp <- par; pp[k] <- pp[k] + hh
    pm <- par; pm[k] <- pm[k] - hh
    g[k] <- tryCatch((tk_at(f, pp, bk, map) - tk_at(f, pm, bk, map)) / (2 * hh), error = function(e) NA_real_)
  }
  out$t_tk <- proc.time()[[3]] - t0
  it <- length(par) - 1L
  if (!is.null(V) && nrow(V) == length(par) && all(is.finite(V[it:length(par), ]))) {
    ok <- is.finite(g)
    if (nzchar(Sys.getenv("REV_TK_RESTRICT"))) {
      g[nm == "beta"] <- 0
      ok <- is.finite(g)
    }
    g0 <- g
    g0[!ok] <- 0
    shift <- as.numeric(V %*% g0)
    shift[!is.finite(shift)] <- NA
    out$b_ls <- shift[it]; out$b_t <- shift[it + 1L]
    t2 <- th[2] + shift[it + 1L]
    out$b_ar <- t2 / sqrt(1 + t2^2) - out$ar
    out$gnan <- sum(!ok)
  }
  # the exact value at this fit's point, where the reference takes it
  xr <- ex[ex$seed == seed, ]
  if (nrow(xr)) {
    out$x_ls <- xr$x_ls[1]; out$x_ar <- xr$x_ar[1]; out$x_ll <- xr$x_ll[1]
  }
  if (set == "m4b") {
    pf <- c(f$estimates$beta[1:2], th[1], pmin(pmax(th[2], -8), 8))
    ib <- which(nm == "betad")
    if (envd$npar > 4) pf <- c(pf, min(par[ib[1]], log(1e5)))
    if (envd$D$dist == "zoib") {
      bd <- par[ib]
      envd$D$zoi <- plogis(bd[2]); envd$D$coi <- plogis(bd[3])
    }
    out$x_at <- -envd$exact_nll(pf, d)
  } else if (set == "fresh" && envd$D$dist %in% c("pois", "binom", "bern", "geom", "expo")) {
    out$x_at <- -envd$exact(c(f$estimates$beta[1:2], th[1], pmin(pmax(th[2], -8), 8)), d)
  } else if (set == "lane" && envd$D$fam %in% c("poisson", "binomial", "bernoulli")) {
    pf <- c(f$estimates$beta[2], th[1], pmin(pmax(th[2], -8), 8), f$estimates$beta[1])
    out$x_at <- tryCatch(-get("exact_nll", envir = lane_exact)(pf, d), error = function(e) NA)
  }
  out$b1 <- f$estimates$beta[if (is_cum) 1 else 2]
  out$s_true <- str_true
  out
}
if (set == "lane") {
  lane_exact <- new.env(parent = envd)
  j <- grep("^row_dens <- ", src)
  k <- grep("^exact_nll <- ", src)
  e_end <- k + which(src[(k + 1):length(src)] == "}")[1]
  eval(parse(text = src[j:e_end]), lane_exact)
}
cols <- c("set","design","seed","err","ll","code","ls","t","ar","warned","boxed","margin","se_ls","se_t","C","C_up","C_dn","C_chi","C_clo","n","t_tk","b_ls","b_t","b_ar","gnan","x_ls","x_ar","x_ll","x_at","b1","s_true","why")
first <- TRUE
for (sd in s1:s2) {
  o <- tryCatch(one(sd), error = function(e) list(set = set, design = design, seed = sd, err = substr(conditionMessage(e), 1, 90)))
  for (cn in setdiff(cols, names(o))) o[[cn]] <- NA
  write.table(as.data.frame(o[cols]), outf, sep = "\t", row.names = FALSE, quote = TRUE, col.names = first, append = !first)
  first <- FALSE
}
