# Follow-ups: frm_allfit, profile intervals and their warning classes on
# a shifted design, and a search for a natural_finish() that relabels a
# stop short of the maximum as converged with no warning.
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
only <- Sys.getenv("REV_ONLY", "")
run_part <- function(p) !nzchar(only) || grepl(only, p)
gd <- function(seed, m, n = 150L, G = 15L) {
  set.seed(seed)
  g <- factor(rep(seq_len(G), length.out = n))
  z <- rnorm(n)
  eta <- 0.5 + 0.4 * z + rnorm(G, 0, 0.7)[g]
  data.frame(y = rpois(n, exp(eta)), z = z, w = z + m, g = g)
}
cls <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(x) {
      w <<- c(w, paste0("[", paste(setdiff(class(x), c("warning", "condition")),
                                    collapse = "/"), "] ",
                        substr(conditionMessage(x), 1, 90)))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  list(v = v, w = unique(w))
}

if (run_part("allfit")) {
  d <- gd(2, 1e6)
  fz <- quiet(frm(bf(y ~ z + (1 | g)), family = poisson(), data = d))
  fw <- quiet(frm(bf(y ~ w + (1 | g)), family = poisson(), data = d))
  bz <- quiet(fixef(fz))["z", 1:2]
  for (sfm in c(TRUE, FALSE)) {
    af <- cls(frm_allfit(fw, start_from_mle = sfm))
    cat("allfit start_from_mle =", sfm, "warnings:", paste(af$w, collapse = " | "), "\n")
    if (inherits(af$v, "error")) { cat("  error", conditionMessage(af$v), "\n"); next }
    for (nm in names(af$v$fits)) {
      o <- af$v$fits[[nm]]
      if (inherits(o, "frmtmb_fit")) {
        cat(sprintf("  %-28s basis %s slope shift %+.2e SE dll %+.2e conv %s\n",
                    nm, basis_on(o), (quiet(fixef(o))["w", 1] - bz[1]) / bz[2],
                    ll_of(o) - ll_of(fz), o$opt$convergence))
      } else cat("  ", nm, "not a fit:", class(o)[1], "\n")
    }
    if (length(af$v$errors)) print(af$v$errors)
  }
}

if (run_part("profile")) {
  for (s in 1:3) {
    d <- gd(s, 1e6)
    fz <- quiet(frm(bf(y ~ z + (1 | g)), family = poisson(), data = d))
    fw <- quiet(frm(bf(y ~ w + (1 | g)), family = poisson(), data = d))
    pz <- cls(confint(fz, parm = "z", method = "profile"))
    pw <- cls(confint(fw, parm = "w", method = "profile"))
    pi <- cls(confint(fw, parm = "theta_1", method = "profile"))
    piz <- cls(confint(fz, parm = "theta_1", method = "profile"))
    show <- function(p) if (inherits(p$v, "error")) paste("ERROR", conditionMessage(p$v)) else
      paste(format(unlist(p$v[1, 1:2]), digits = 6), collapse = " ")
    cat(sprintf("seed %d profile slope: w %s | z %s\n   w warnings: %s\n", s,
                show(pw), show(pz), paste(pw$w, collapse = " | ")))
    cat(sprintf("   profile theta: w %s | z %s\n   w warnings: %s\n",
                show(pi), show(piz), paste(pi$w, collapse = " | ")))
  }
}

if (run_part("natfin")) {
  cat("== iteration caps 1 to 4, seeds 1 to 25, m = 1e6 and 1e7\n")
  tab <- NULL
  for (m in c(1e6, 1e7)) for (s in 1:25) {
    d <- gd(s, m)
    fz <- quiet(frm(bf(y ~ z + (1 | g)), family = poisson(), data = d))
    for (it in 1:4) {
      ctl <- frmtmb_control(optCtrl = list(iter.max = it, eval.max = 400))
      r <- rfit(bf(y ~ w + (1 | g)), family = poisson(), data = d,
                control = ctl)
      if (r$err) next
      nf <- !is.null(r$f$opt$natural_finish)
      dll <- ll_of(r$f) - ll_of(fz)
      tab <- rbind(tab, data.frame(m = m, seed = s, it = it, nf = nf,
                                   conv = r$f$opt$convergence, dll = dll,
                                   warned = length(r$w) > 0,
                                   w = substr(short_w(r), 1, 90)))
    }
  }
  cat("runs", nrow(tab), "natural_finish kept", sum(tab$nf), "\n")
  print(tab[tab$nf, ], row.names = FALSE)
  cat("silent and short by more than 1e-4:\n")
  print(tab[!tab$warned & tab$dll < -1e-4, ], row.names = FALSE)
}
