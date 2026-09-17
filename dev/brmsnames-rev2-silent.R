## Reviewer RECHECK (round 1 build), copy of dev/brmsnames-rev-silent.R
## writing its own cache files. Claim 6: same numbers under the renaming. Run once per arm;
## each arm saves what its accessors return, then `compare` pairs them by
## POSITION and VALUE, never by name.
##   Rscript dev/brmsnames-rev2-silent.R base
##   Rscript dev/brmsnames-rev2-silent.R lane
##   Rscript dev/brmsnames-rev2-silent.R compare
## Data seed 31, sampler seed 5.
arm <- commandArgs(trailingOnly = TRUE)[1L]
libs <- list(
  base = c("C:/Users/adf44/source/r/rellib-r3"),
  lane = c("C:/Users/adf44/source/r/brmsnames-lib",
           "C:/Users/adf44/source/r/rellib-r3"),
  compare = c("C:/Users/adf44/source/r/rellib-r3"))
.libPaths(c(libs[[arm]], "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache"))
q <- function(e) suppressWarnings(suppressMessages(e))
outf <- function(a) sprintf("dev/stan-cache/brmsnames-rev2-silent-%s.rds", a)

if (arm != "compare") {
  q(library(frmtmb)); q(library(frmtmb.sample))
  cat("arm", arm, "frmtmb from", dirname(system.file(package = "frmtmb")), "\n")
  set.seed(31)
  G <- 10; n <- 200
  d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:G, 20)))
  u <- matrix(rnorm(G * 3, 0, 0.5), G)
  d$y <- 1 + 0.5 * d$x + u[d$g, 1] + u[d$g, 2] * d$x +
    rnorm(n, 0, exp(0.2 * d$z + u[d$g, 3]))
  d$y1 <- d$y; d$y2 <- 0.4 * d$x + u[d$g, 1] + rnorm(n)
  d$az <- abs(d$z)
  d$yn <- 2 * exp(-0.7 * d$az) * exp(u[d$g, 1] / 5) + rnorm(n, 0, 0.1)
  d$yo <- cut(d$y, c(-Inf, 0, 1, 2, Inf), labels = c("a", "b", "c", "d"),
              ordered_result = TRUE)
  ms <- list(
    dist = list(bf(y ~ x + (1 + x | g), sigma ~ z + (1 | g)), gaussian()),
    mv = list(mvbf(bf(y1 ~ x + (1 | p | g)), bf(y2 ~ x + (1 | p | g)),
                   rescor = FALSE), gaussian()),
    nl = list(bf(yn ~ a * exp(-b * az), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
              gaussian()),
    ord = list(bf(yo ~ x + (1 | g)), cumulative()),
    smooth = list(bf(y ~ s(x) + (1 | g)), gaussian()),
    gp = list(bf(y ~ gp(x)), gaussian())
  )
  res <- list(formals = list(
    fixef = names(formals(getS3method("fixef", "frmtmb_fit"))),
    ranef = names(formals(getS3method("ranef", "frmtmb_fit"))),
    coef = names(formals(getS3method("coef", "frmtmb_fit"))),
    VarCorr = names(formals(getS3method("VarCorr", "frmtmb_fit"))),
    hypothesis = names(formals(getS3method("hypothesis", "frmtmb_fit")))))
  for (nm in names(ms)) {
    fit <- q(frm(ms[[nm]][[1]], family = ms[[nm]][[2]], data = d,
                 start = if (nm == "nl") list(beta = c(2, 0.7))))
    vo <- frmtmb:::hyp_vals_only(fit)
    env <- frmtmb:::hyp_env_vals(fit, vo$vals, vo$comp)
    vcm <- if (arm == "lane") varcorr_matrices(fit) else VarCorr(fit)
    try1 <- function(e) tryCatch(e, error = function(err) paste("ERROR:", conditionMessage(err)))
    ds <- q(frm_sample(fit, chains = 1, iter = 120, refresh = 0, seed = 5))
    r <- list(
      fixef = fixef(fit), fixef_flat = fixef(fit, flatten = TRUE),
      ranef = lapply(unclass(ranef(fit)), function(m) { attributes(m)[c("term")] <- NULL; m }),
      coef = try1(coef(fit)), vcm = lapply(unclass(vcm), unclass),
      confint = try1(confint(fit)), env = unlist(env),
      pos_fixef_TRUE = try1(fixef(fit, TRUE)),
      pos_ranef_TRUE = try1(ranef(fit, TRUE)),
      draws = unname(ds$draws), draws_names = colnames(ds$draws),
      variables_fit = try1(variables(fit)),
      variables_ds = try1(variables(ds))
    )
    if (arm == "lane") {
      r$VarCorr <- try1(VarCorr(fit))
      r$ranef_ds <- try1(ranef(ds, summary = FALSE))
      r$coef_ds <- try1(coef(ds, summary = FALSE))
      r$fixef_ds <- try1(fixef(ds, summary = FALSE))
      r$VarCorr_ds <- try1(VarCorr(ds, summary = FALSE))
      # independent per-draw reference: the ML accessors at each draw
      idx <- frmtmb.sample:::draws_par_index(fit)
      r$at_draw <- lapply(c(1, 30, 60), function(i) {
        sh <- frmtmb.sample:::draws_fit_at(ds, i, idx)
        list(i = i, ranef = lapply(unclass(ranef(sh)), unclass),
             vcm = lapply(unclass(varcorr_matrices(sh)), unclass),
             fixef = fixef(sh))
      })
    }
    res[[nm]] <- r
    cat(nm, "done\n")
  }
  saveRDS(res, outf(arm))
  cat("saved\n")
} else {
  b <- readRDS(outf("base")); l <- readRDS(outf("lane"))
  cat("== formals ==\n")
  for (k in names(b$formals)) cat(sprintf("%-10s base: %s\n%-10s lane: %s\n", k,
      paste(b$formals[[k]], collapse = ","), "", paste(l$formals[[k]], collapse = ",")))
  for (nm in setdiff(names(b), "formals")) {
    B <- b[[nm]]; L <- l[[nm]]
    cat("\n== model", nm, "==\n")
    # round 1 stores an unmodeled dpar on its natural scale: count the
    # draws columns that are identical, and check each other one against
    # exp() of base, the only link these six models use there
    sameB <- vapply(seq_len(ncol(B$draws)), function(j)
      identical(B$draws[, j], L$draws[, j]), TRUE)
    cat("  draws columns identical:", sum(sameB), "of", length(sameB), "\n")
    for (j in which(!sameB)) {
      cat(sprintf("    %s -> %s: max rel |lane - exp(base)| %.3g\n",
                  B$draws_names[j], L$draws_names[j],
                  max(abs(L$draws[, j] - exp(B$draws[, j])) /
                        exp(B$draws[, j]))))
    }
    for (k in c("fixef", "fixef_flat", "ranef", "coef", "vcm", "confint", "draws")) {
      cat(sprintf("  %-12s identical: %s\n", k, identical(B[[k]], L[[k]])))
      if (!identical(B[[k]], L[[k]])) {
        ae <- all.equal(B[[k]], L[[k]])
        cat("     all.equal:", head(ae, 4), sep = "\n     ")
      }
    }
    cat("  positional fixef(fit, TRUE): base class", class(B$pos_fixef_TRUE),
        "lane class", class(L$pos_fixef_TRUE), " identical:",
        identical(B$pos_fixef_TRUE, L$pos_fixef_TRUE), "\n")
    rb <- B$pos_ranef_TRUE; rl <- L$pos_ranef_TRUE
    cat("  positional ranef(fit, TRUE): base condSD attr present:",
        !is.null(attr(unclass(rb)[[1]], "condSD")), " lane:",
        if (is.character(rl)) rl else !is.null(attr(unclass(rl)[[1]], "condSD")), "\n")
    # env values: pair by value
    eb <- B$env; el <- L$env
    cat("  env: base", length(eb), "names, lane", length(el), "names\n")
    for (i in seq_along(el)) {
      hit <- names(eb)[eb == el[i]]
      cat(sprintf("     lane %-40s = base %s\n", names(el)[i],
                  if (length(hit)) paste(hit, collapse = " | ") else "<NO BASE VALUE>"))
    }
    miss <- names(eb)[!eb %in% el]
    if (length(miss)) cat("     base names whose value no lane name carries:", miss, "\n")
    cat("  draws names by position (base -> lane), first differing 12:\n")
    dn <- which(B$draws_names != L$draws_names)
    for (i in head(dn, 12)) cat("     ", B$draws_names[i], "->", L$draws_names[i], "\n")
    # r_ columns against the ML ranef() at the same draw
    bad <- 0; checked <- 0
    for (ad in L$at_draw) {
      row <- L$draws[ad$i, ]
      names(row) <- L$draws_names
      for (gk in seq_along(ad$ranef)) {
        M <- ad$ranef[[gk]]
        g <- names(ad$ranef)[gk]
        for (lv in rownames(M)) for (cn in colnames(M)) {
          pat <- paste0("^r_", g, "(__[^[]+)?\\[", lv, ",", gsub("[()]", "", cn), "\\]$")
          cols <- grep(pat, names(row))
          if (length(cols) == 0) next
          checked <- checked + 1
          if (!any(row[cols] == M[lv, cn])) bad <- bad + 1
        }
      }
    }
    cat("  r_ columns checked against ranef() at draws 1,30,60:", checked,
        " mismatched:", bad, "\n")
  }
}
