# Reviewer, punch round 1: the cmc = FALSE refusal on level-read
# covariance structures. Where does it fire, and on which correct models?
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages(expr)),
           error = function(e) structure(conditionMessage(e), class = "ERR"))
}
set.seed(41)
G <- 30; Tn <- 4
d <- data.frame(g = factor(rep(seq_len(G), each = Tn)),
                t = factor(rep(seq_len(Tn), G)),
                x1 = rnorm(G * Tn), x2 = rnorm(G * Tn),
                f = factor(sample(c("a", "b", "c"), G * Tn, TRUE)))
d$tim <- num_factor(rep(seq_len(Tn), G))
d$pos <- num_factor(rep(c(0, 1, 0, 1), G), rep(c(0, 0, 1, 1), G))
d$y <- rnorm(G * Tn, 1 + d$x1) + rnorm(G, 0, 0.7)[d$g]
res <- function(r) {
  if (inherits(r, "ERR")) paste("REFUSED:", substr(r, 1, 110)) else {
    cn <- unlist(lapply(r$re_blocks, function(b)
      lapply(b$components, `[[`, "cnms")))
    paste("frame ok, cnms", paste(cn, collapse = ","))
  }
}
fr <- function(f) res(q(frm(f, data = d, dry_run = "frame")))
lvl <- c("ar1", "hetar1", "cs", "homcs", "toep", "homtoep")
dist <- c("ou", "exp", "gau", "mat")
cat("== factor LHS: cmc default, TRUE, FALSE\n")
for (cs in c(lvl, dist)) {
  lhs <- if (cs == "ou") "0 + tim" else if (cs %in% dist) "0 + pos" else
    "0 + t"
  rhs <- sprintf("y ~ x1 + %s(%s | g)", cs, lhs)
  for (cm in list(NULL, TRUE, FALSE)) {
    cat(sprintf("%-32s cmc=%-5s %s\n", rhs, format(cm %||% "dflt"),
                fr(bf(as.formula(rhs), cmc = cm))))
  }
}
cat("== numeric LHS, cmc = FALSE (the same columns either way)\n")
for (cs in lvl) {
  rhs <- sprintf("y ~ x1 + %s(0 + x1 + x2 | g)", cs)
  cat(sprintf("%-38s TRUE : %s\n", rhs, fr(bf(as.formula(rhs)))))
  cat(sprintf("%-38s FALSE: %s\n", rhs,
              fr(bf(as.formula(rhs), cmc = FALSE))))
}
cat("== the structure with an intercept, cmc = FALSE\n")
cat("cs(1 + x1 | g):", fr(bf(y ~ x1 + cs(1 + x1 | g), cmc = FALSE)), "\n")
cat("== cmc = FALSE on another formula than the structure's\n")
cat("mu ar1(0 + t | g), lf(sigma ~ 0 + f, cmc = FALSE):",
    fr(bf(y ~ x1 + ar1(0 + t | g)) + lf(sigma ~ 0 + f, cmc = FALSE)), "\n")
cat("mu 0 + f cmc = FALSE, lf(sigma ~ 1 + ar1(0 + t | g)):",
    fr(bf(y ~ 0 + f, cmc = FALSE) + lf(sigma ~ 1 + ar1(0 + t | g))), "\n")
cat("== unstructured and diagonal, cmc = FALSE\n")
cat("us(0 + f | g):", fr(bf(y ~ x1 + us(0 + f | g), cmc = FALSE)), "\n")
cat("diag(0 + f | g):", fr(bf(y ~ x1 + diag(0 + f | g), cmc = FALSE)), "\n")
cat("(0 + f | g):", fr(bf(y ~ x1 + (0 + f | g), cmc = FALSE)), "\n")
cat("== the autocorrelation term ar(), cmc = FALSE\n")
cat("frm 0 + f + ar(t, g):",
    fr(bf(y ~ 0 + f + ar(t, g), cmc = FALSE)), "\n")
d$tn <- as.integer(d$t)
sd <- q(brms::standata(brms::bf(y ~ 0 + f + ar(tn, g), cmc = FALSE), d))
cat("brms 0 + f + ar(tn, g):", if (inherits(sd, "ERR")) sd else
  paste("ok, X", paste(colnames(sd$X), collapse = ",")), "\n")
cat("== brms on the structures themselves\n")
for (cs in c("ar1", "cs", "toep", "ou")) {
  f <- as.formula(sprintf("y ~ x1 + %s(0 + t | g)", cs))
  r <- q(brms::standata(brms::bf(f, cmc = FALSE), d))
  cat(cs, "brms cmc = FALSE:", if (inherits(r, "ERR")) substr(r, 1, 110) else
    paste("ok", paste(grep("^Z_", names(r), value = TRUE), collapse = ",")),
    "\n")
}
cat("== full fits where the guard is absent\n")
f1 <- q(frm(bf(y ~ x1 + ar1(0 + t | g)), data = d))
cat("ar1 cmc default fit:", if (inherits(f1, "ERR")) f1 else
  as.numeric(logLik(f1)), "\n")
f2 <- q(frm(bf(y ~ x1 + ar1(0 + t | g), cmc = TRUE), data = d))
cat("ar1 cmc = TRUE identical logLik:", identical(logLik(f1), logLik(f2)),
    "\n")
cat("DONE\n")
