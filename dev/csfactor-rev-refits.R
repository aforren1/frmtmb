# Claim 1: what the new identifiability check does to every in-package
# refit when a cs() factor level, or an ordinal response category, is
# absent from the rows a refit sees.
#   Rscript dev/csfactor-rev-refits.R <lib> <tag>
# Run with the lane library and with rellib-r3 to compare.
args <- commandArgs(TRUE)
LIB <- args[[1L]]; TAG <- args[[2L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("== build", TAG, ":", as.character(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), "\n")
options(warn = 1)

say <- function(...) cat(..., "\n", sep = "")
hdr <- function(x) cat("\n---- ", x, " ----\n", sep = "")
short <- function(e) {
  m <- conditionMessage(e)
  substr(gsub("[\r\n]+", " ", m), 1L, 220L)
}
try_fit <- function(expr) {
  tryCatch(suppressWarnings(eval.parent(substitute(expr))),
           error = function(e) structure(list(msg = short(e)),
                                         class = "revfail"))
}

# ---------------------------------------------------------------
# The data: 3-level cs() factor, one level with a SINGLE row, and an
# ordinal response with 3 categories, one of them rare.
set.seed(90291)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n))
# level "c" gets exactly one row; "b" gets 40
d$f <- factor(c(rep("a", n - 41L), rep("b", 40L), "c"),
              levels = c("a", "b", "c"))
eta <- 0.6 * d$x + c(a = 0, b = 0.8, c = -0.5)[as.character(d$f)]
pr <- function(e) {
  p1 <- plogis(-0.6 - e); p2 <- plogis(0.9 - e)
  cbind(p1, p2 - p1, 1 - p2)
}
P <- pr(eta)
d$yo <- apply(P, 1L, function(p) sample.int(3L, 1L, prob = pmax(p, 1e-9)))
say("table(f) = ", paste(table(d$f), collapse = " "))
say("table(yo) = ", paste(table(d$yo), collapse = " "))

fit <- try_fit(frm(bf(yo ~ x + cs(f)), family = sratio(), data = d))
if (inherits(fit, "revfail")) {
  say("BASE FIT FAILED: ", fit$msg); quit(save = "no")
}
say("base logLik = ", sprintf("%.6f", as.numeric(logLik(fit))),
    "  npar = ", length(fit$opt$par))
say("coef names: ", paste(rownames(fixef(fit)), collapse = ", "))

# ---------------------------------------------------------------
hdr("1a influence() observation-wise (deletes the only row of level c)")
inf <- try_fit(influence(fit, force = TRUE))
if (inherits(inf, "revfail")) {
  say("influence() ERRORED: ", inf$msg)
} else {
  fe <- inf$fixed
  naro <- which(apply(fe, 1L, function(r) all(is.na(r))))
  say("influence rows = ", nrow(fe), "  all-NA rows = ", length(naro))
  say("all-NA row ids = ", paste(head(names(naro), 10), collapse = ","))
  say("row of the single 'c' observation = ", which(d$f == "c"))
  cd <- try_fit(cooks.distance(inf))
  if (inherits(cd, "revfail")) {
    say("cooks.distance ERRORED: ", cd$msg)
  } else {
    say("cooks.distance NA count = ", sum(is.na(cd)), " of ", length(cd))
    say("max cooks (na.rm) = ", sprintf("%.6g", max(cd, na.rm = TRUE)))
  }
}

hdr("1a' the refusal message a deleted singleton level produces")
sub <- d[-which(d$f == "c"), , drop = FALSE]
say("levels kept in the subset: ", paste(levels(sub$f), collapse = ","),
    "   counts ", paste(table(sub$f), collapse = " "))
r <- tryCatch({
  fr <- frmtmb:::assemble_frame(fit$spec, sub, sparse_x = FALSE,
                                data2 = list())
  paste0("assembled, cs entries = ", length(fr$lp[[1L]][["cs"]]))
}, error = function(e) paste0("REFUSED: ", short(e)))
say(r)

hdr("1b frm_bootstrap() nsim=200 on the cs(f) model")
t0 <- proc.time()[["elapsed"]]
bs <- try_fit(frm_bootstrap(fit, nsim = 200, seed = 11))
if (inherits(bs, "revfail")) {
  say("frm_bootstrap ERRORED: ", bs$msg)
} else {
  nar <- sum(apply(bs$t, 1L, function(r) all(is.na(r))))
  say("nsim = ", bs$nsim, "  all-NA replicate rows = ", nar,
      "  not converged = ", sum(!bs$converged))
  say("t0 names: ", paste(names(bs$t0), collapse = ","))
  say("elapsed = ", sprintf("%.1f", proc.time()[["elapsed"]] - t0), "s")
}

hdr("1b' a bootstrap replicate whose response lacks a category")
# force it: refit() to a response that never takes the value 3
yr <- d$yo; yr[yr == 3L] <- 2L
say("forced response table = ", paste(table(yr), collapse = " "))
rf <- try_fit(refit(fit, yr))
if (inherits(rf, "revfail")) {
  say("refit on a category-short response ERRORED: ", rf$msg)
} else {
  say("refit logLik = ", sprintf("%.6f", as.numeric(logLik(rf))),
      "  npar = ", length(rf$opt$par))
}

hdr("1c frm_allfit()")
af <- try_fit(frm_allfit(fit))
if (inherits(af, "revfail")) {
  say("frm_allfit ERRORED: ", af$msg)
} else {
  say("allfit class = ", paste(class(af), collapse = "/"))
  print(af)
}

hdr("1d anova(refit = TRUE)")
fit2 <- try_fit(frm(bf(yo ~ x + z + cs(f)), family = sratio(), data = d))
if (inherits(fit2, "revfail")) {
  say("second fit FAILED: ", fit2$msg)
} else {
  a <- try_fit(anova(fit, fit2, refit = TRUE))
  if (inherits(a, "revfail")) say("anova ERRORED: ", a$msg) else print(a)
}
hdr("1d' drop1()")
d1 <- try_fit(drop1(fit2, test = "Chisq"))
if (inherits(d1, "revfail")) say("drop1 ERRORED: ", d1$msg) else print(d1)

hdr("1e conditional_effects(band = 'boot')")
ce <- try_fit(conditional_effects(fit, effects = "f", band = "boot",
                                  boot = 30, seed = 5))
if (inherits(ce, "revfail")) {
  say("conditional_effects boot ERRORED: ", ce$msg)
} else {
  df1 <- as.data.frame(ce[[1L]])
  say("rows = ", nrow(df1), "  NA in estimate = ",
      sum(is.na(df1[["estimate__"]])), "  NA in lower = ",
      sum(is.na(df1[["lower__"]])))
}

hdr("1f the autoscale pre-fit")
af2 <- try_fit(frm(bf(yo ~ x + cs(f)), family = sratio(), data = d,
                   control = frmtmb_control(autoscale = "always")))
if (inherits(af2, "revfail")) {
  say("autoscale fit ERRORED: ", af2$msg)
} else {
  say("autoscale logLik = ", sprintf("%.6f", as.numeric(logLik(af2))))
}

hdr("1g frm_multiple() over imputations, one missing level c's row")
imps <- list(d, d, d)
imps[[2L]] <- d[-which(d$f == "c"), , drop = FALSE]      # level declared,
imps[[3L]] <- droplevels(d[-which(d$f == "c"), , drop = FALSE])  # dropped
say("imp2 levels ", paste(levels(imps[[2L]]$f), collapse = ","),
    "  imp3 levels ", paste(levels(imps[[3L]]$f), collapse = ","))
fm <- try_fit(frm_multiple(bf(yo ~ x + cs(f)), family = sratio(),
                           data = imps))
if (inherits(fm, "revfail")) {
  say("frm_multiple ERRORED: ", fm$msg)
} else {
  say("frm_multiple fits = ", length(fm$fits %||% list()))
}
# and the same over three imputations that all keep every level
imps2 <- lapply(1:3, function(i) { dd <- d; dd$x <- d$x + rnorm(n, 0, .01); dd })
fm2 <- try_fit(frm_multiple(bf(yo ~ x + cs(f)), family = sratio(),
                            data = imps2))
if (inherits(fm2, "revfail")) {
  say("frm_multiple (all levels) ERRORED: ", fm2$msg)
} else {
  say("frm_multiple (all levels) fits = ", length(fm2$fits %||% list()))
}

hdr("1h influence(groups = ) where a group holds the only 'c' row")
d$g <- factor(rep(1:12, each = 10))
say("group of the 'c' row = ", as.character(d$g[d$f == "c"]))
fitg <- try_fit(frm(bf(yo ~ x + cs(f)), family = sratio(), data = d))
ing <- try_fit(influence(fitg, groups = "g"))
if (inherits(ing, "revfail")) {
  say("influence(groups) ERRORED: ", ing$msg)
} else {
  naro <- which(apply(ing$fixed, 1L, function(r) all(is.na(r))))
  say("group units = ", nrow(ing$fixed), "  all-NA units = ", length(naro),
      "  ids ", paste(names(naro), collapse = ","))
}
cat("\nDONE ", TAG, "\n")
