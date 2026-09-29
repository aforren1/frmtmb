# Re-check 2 and 3: the new influence() warning, and the branched alias
# advice. Seed 90291 is the round-1 construction, unchanged.
#   Rscript dev/csfactor-rev2-inf.R <lib> <tag>
args <- commandArgs(TRUE)
LIB <- args[[1L]]; TAG <- args[[2L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("== build", TAG, ":", as.character(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), "\n")
say <- function(...) cat(..., "\n", sep = "")
short <- function(e) gsub("[\r\n]+", " ", conditionMessage(e))

# ------------------------------------------------ round 1's own data
set.seed(90291)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$f <- factor(c(rep("a", n - 41L), rep("b", 40L), "c"),
              levels = c("a", "b", "c"))
eta <- 0.6 * d$x + c(a = 0, b = 0.8, c = -0.5)[as.character(d$f)]
pr <- function(e) {
  p1 <- plogis(-0.6 - e); p2 <- plogis(0.9 - e)
  cbind(p1, p2 - p1, 1 - p2)
}
P <- pr(eta)
d$yo <- apply(P, 1L, function(p) sample.int(3L, 1L, prob = pmax(p, 1e-9)))
say("table(f) = ", paste(table(d$f), collapse = " "),
    "   the single 'c' row is ", which(d$f == "c"))

# one warning, or none: captured rather than suppressed
run_inf <- function(lab, fit) {
  ws <- character(0)
  inf <- withCallingHandlers(
    suppressMessages(influence(fit, force = TRUE)),
    warning = function(w) {
      ws <<- c(ws, short(w)); invokeRestart("muffleWarning")
    })
  cd <- suppressWarnings(cooks.distance(inf))
  say("\n[", lab, "]")
  say("  influence-table columns: ",
      paste(colnames(inf$fixed), collapse = ", "))
  say("  warnings: ", length(ws))
  for (w in ws) say("    W: ", substr(w, 1L, 420L))
  say("  cooks.distance NA count = ", sum(is.na(cd)), " of ", length(cd))
  say("  rows with any NA in the coefficient table = ",
      sum(apply(inf$fixed, 1L, anyNA)))
  invisible(list(ws = ws, cd = cd, inf = inf))
}

fit_cs <- frm(bf(yo ~ x + cs(f)), family = sratio(), data = d)
say("cs fit coef names (get_coef): ",
    paste(names(frmtmb:::get_coef.frmtmb_fit(fit_cs)), collapse = ", "))
say("cs par components: ", paste(names(fit_cs$estimates), collapse = ", "))
a <- run_inf("cs(f) with a singleton level: the warning case", fit_cs)

# control 1: no cs() at all, same data and the same singleton level
fit_no <- frm(bf(yo ~ x + f), family = sratio(), data = d)
b <- run_inf("no cs(): f as an ordinary factor", fit_no)

# control 2: cs() on a factor whose every level has many rows, so no
# deletion loses a coefficient
set.seed(90291)
d2 <- d
d2$f2 <- factor(rep(c("a", "b", "c"), length.out = n))
fit_ok <- frm(bf(yo ~ x + cs(f2)), family = sratio(), data = d2)
say("\ntable(f2) = ", paste(table(d2$f2), collapse = " "))
cc <- run_inf("cs(f2), every level with many rows", fit_ok)

# control 3: cs() on a numeric column, no levels to lose
dd <- run_inf("cs(z) numeric", frm(bf(yo ~ x + cs(z)), family = sratio(),
                                   data = d))

say("\n=== the direct mapping helper, called on both spellings ===")
say("influence_coef_labels(fit, c('fcc[1]','bcs3_1','x')) = ",
    paste(frmtmb:::influence_coef_labels(
      fit_cs, c("fcc[1]", "bcs3_1", "x")), collapse = " | "))
say("the cs par names in this fit: ",
    paste(vapply(fit_cs$frame$linpreds[["yo.mu"]][["cs"]], `[[`, "",
                 "par"), collapse = ", "))
say("the cs labels in this fit: ",
    paste(vapply(fit_cs$frame$linpreds[["yo.mu"]][["cs"]], `[[`, "",
                 "label"), collapse = ", "))

say("\n=== 3  the alias advice, one line per formula ===")
set.seed(1907)
n2 <- 300
e <- data.frame(x = rnorm(n2), z = rnorm(n2))
e$m <- factor(sample(0:3, n2, TRUE), levels = 0:3, ordered = TRUE)
et <- 0.5 * e$x
q1 <- plogis(-0.7 - et); q2 <- plogis(0.8 - et)
Q <- cbind(q1, q2 - q1, 1 - q2)
e$yo <- apply(Q, 1L, function(p) sample.int(3L, 1L, prob = pmax(p, 1e-9)))
shapes <- list(
  "x + cs(x)"            = quote(bf(yo ~ x + cs(x))),
  "poly(x,2) + cs(x)"    = quote(bf(yo ~ poly(x, 2) + cs(x))),
  "s(x) + cs(x)"         = quote(bf(yo ~ s(x) + cs(x))),
  "x + z + cs(I(x+z))"   = quote(bf(yo ~ x + z + cs(I(x + z)))),
  "mo(m) + cs(m)"        = quote(bf(yo ~ mo(m) + cs(m))),
  "z + mo(m):z + cs(m)"  = quote(bf(yo ~ z + mo(m):z + cs(m))),
  "I(x^2) + cs(x)"       = quote(bf(yo ~ I(x^2) + cs(x)))
)
for (nmx in names(shapes)) {
  r <- tryCatch(suppressWarnings(suppressMessages(
         frm(eval(shapes[[nmx]]), family = sratio(), data = e))),
        error = function(ex) structure(list(m = short(ex)),
                                       class = "revfail"))
  say("\n-- ", nmx, " --")
  if (inherits(r, "revfail")) {
    cat(strwrap(paste("REFUSED:", r$m), 74), sep = "\n  ")
    cat("\n")
  } else {
    say("  FITTED, df = ", attr(logLik(r), "df"), ", logLik = ",
        sprintf("%.9f", as.numeric(logLik(r))),
        ", non-finite se = ",
        sum(!is.finite(sqrt(diag(suppressWarnings(vcov(r)))))))
  }
}
say("\nDONE ", TAG)
