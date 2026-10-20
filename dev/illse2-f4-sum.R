# Lane illse2, fix round 1 resumed: every table of the findings from
# the logs of dev/illse2-f4-jobs.sh (dev/illse2-log/f4/), and the base
# rows from the logs of 0.71.0 (rellib-r9) the lane wrote before.
#   Rscript dev/illse2-f4-sum.R > dev/illse2-log/f4/summary.txt
D <- "C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-log"
F4 <- file.path(D, "f4")
rd <- function(f) readLines(f, warn = FALSE)
num <- function(x, re) as.numeric(sub(re, "\\1", regmatches(x, regexpr(re, x))))
e1 <- function(v) if (length(v)) sprintf("%.1e", v) else "-"
rng <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return("-")
  if (length(unique(signif(v, 2))) == 1L) return(e1(v[1]))
  paste0(e1(min(v)), "-", e1(max(v)))
}

cat("== 1. cooks.distance, D1, D2 (largest |ratio - 1|, seeds 1:3)\n")
for (re in 0:1) for (deg in 5:8) {
  x <- rd(file.path(F4, sprintf("p1-lane-%d-%d.txt", re, deg)))
  one <- function(tag) {
    l <- grep(tag, x, value = TRUE)
    nerr <- sum(grepl("ERROR", l))
    ok <- l[!grepl("ERROR", l)]
    # the first "stat" is the ratio; the later "exact stat" the value
    st <- regmatches(ok, regexpr("stat [0-9.]+", ok))
    v <- c(num(ok, ".*max dev ([0-9.e+-]+).*"),
           abs(as.numeric(sub("stat ", "", st)) - 1))
    paste0(if (length(v)) e1(max(v)) else "",
           if (nerr) sprintf(" %d of %d error", nerr, length(l)) else "")
  }
  warn <- sum(grepl("warn [^-]", grep("^  (cooks|D1 all|D2 all)", x,
                                       value = TRUE)))
  cat(sprintf("lane %-6s deg %d | cooks %-8s | D1 all %-8s | D2 all %-8s | warned lines %d\n",
              if (re) "(1|g)" else "no RE", deg, one("^  cooks"),
              one("^  D1 all"), one("^  D2 all"), warn))
}

cat("\n== 2. degree 7 seed 4\n")
x <- rd(file.path(F4, "p2-lane-7-4.txt"))
cat(grep("^frmtmb:|lost:|profile max|vcov[(]full|exact SE", x, value = TRUE),
    sep = "\n")

cat("\n== 2. the 80 sd-at-zero fits (short = exact max - fit)\n")
e4 <- function(lines) {
  sh <- num(lines, ".*short ([+-][0-9.e+-]+).*")
  v <- sub(".*theta_1 ([a-z:]+|SE).*", "\\1", lines)
  v <- sub("lost:", "", v)
  c(fits = length(lines), gt1e3 = sum(sh > 1e-3),
    mid = sum(sh > 1e-6 & sh <= 1e-3), le1e6 = sum(sh <= 1e-6),
    table(factor(v, c("SE", "boundary", "nonfinite", "short", "flat"))),
    escape = sum(!grepl("escape -", lines)),
    warned = sum(!grepl("warn -$", lines)))
}
lane <- unlist(lapply(5:8, function(d) {
  grep("^deg", rd(file.path(F4, sprintf("edge4-ref-%d.txt", d))), value = TRUE)
}))
base <- grep("^deg", rd(file.path(D, "edge4-base-ref.txt")), value = TRUE)
print(rbind(base = e4(base), lane = e4(lane)))
key <- function(l) sub(" [|].*", "", l)
bm <- setNames(base, key(base)); lm_ <- setNames(lane, key(lane))
cat("fits whose theta_1 verdict, escape or short bucket differ:\n")
bk <- function(l) {
  s <- num(l, ".*short ([+-][0-9.e+-]+).*")
  paste(sub(".*(theta_1 [^|]*).*", "\\1", l),
        cut(s, c(-Inf, 1e-6, 1e-3, Inf)))
}
for (k in names(lm_)) {
  if (!identical(bk(bm[[k]]), bk(lm_[[k]]))) {
    cat("  base", bm[[k]], "\n  lane", lm_[[k]], "\n")
  }
}
cat("lane fits short by 1e-6 to 1e-3:\n")
s <- num(lane, ".*short ([+-][0-9.e+-]+).*")
cat(paste(" ", lane[s > 1e-6 & s <= 1e-3]), sep = "\n")

cat("\n== 3. nonlinear bodies\n")
for (c in c("poly", "sin", "rat")) {
  cat(grep("seed", rd(file.path(F4, sprintf("p3nl-lane-%s.txt", c))),
           value = TRUE), sep = "\n")
}
cat("\n== 3. rational body seed 3\n")
cat(rd(file.path(F4, "p3ratdbg-lane-3.txt"))[-1], sep = "\n")
cat("\n== 3/4. asymmetry over 27 fits\n")
cat(tail(rd(file.path(F4, "asym2.txt")), 1), "\n")

cat("\n== 3. REML (gaussian, (1 | g), seeds 1:3)\n")
for (deg in 5:8) {
  x <- rd(file.path(F4, sprintf("p3reml-lane-%d.txt", deg)))
  h <- grep("^seed", x, value = TRUE)
  se <- grep("frmtmb vs exact", x, value = TRUE)
  ft <- grep("fitted[(]re_formula", x, value = TRUE)
  th <- grep("share of the variance", x, value = TRUE)
  cat(sprintf(paste0("lane deg %d | short %s | SE dev %s | cor %s | ",
                     "fitted SE dev %s | theta share %s | warn %s | vcov warn %s\n"),
              deg, rng(num(h, ".*REML short ([0-9.e+-]+).*")),
              rng(num(se, ".*max dev ([0-9.e+-]+)[)].*")),
              rng(num(se, ".*cor ([0-9.e+-]+).*")),
              rng(num(ft, ".*max dev ([0-9.e+-]+)[)].*")),
              paste(range(c(num(th, ".*adds: ([0-9.e+-]+) to.*"),
                            num(th, ".*to ([0-9.e+-]+)$"))), collapse = "-"),
              paste(unique(sub(".*\\| warn ([^|]*) \\|.*", "\\1", h)),
                    collapse = ","),
              paste(unique(sub(".*vcov warn (.*)$", "\\1", h)),
                    collapse = ",")))
  l4 <- grep("^  lme4", x, value = TRUE)
  if (length(l4)) cat("  lme4:", paste(sub(".*max dev ([0-9.e+-]+)[)].*",
                                           "\\1", l4), collapse = " "),
                      "| lme4 short", paste(sub(".*REML short ([0-9.e+-]+).*",
                                                "\\1", l4), collapse = " "),
                      "\n")
  else cat("  lme4:", unique(sub("^ *", "", grep("lme4", x, value = TRUE))),
           "\n")
}
cat("\n== 3. REML log-likelihood\n")
cat(rd(file.path(F4, "remlll-lane.txt"))[-1], sep = "\n")
cat("\n== 3. non-gaussian REML\n")
for (fam in c("bernoulli", "poisson")) {
  cat(grep("seed", rd(file.path(F4, sprintf("p3remlglm-lane-%s.txt", fam))),
           value = TRUE), sep = "\n")
}

cat("\n== 4. bounds (20 fits)\n")
x <- grep("^deg", rd(file.path(F4, "p4bound-lane.txt")), value = TRUE)
ll <- num(x, ".*logLik - lm ([+-][0-9.e+-]+).*")
sd <- suppressWarnings(num(x, ".*max dev ([0-9.e+-]+)[)].*"))
at <- abs(ll) < 1e-4
cat(sprintf("at the constrained max (|logLik - lm| < 1e-4): %d; warned and short: %d\n",
            sum(at), sum(!at)))
cat(sprintf("  largest short among them %.1e, SE dev up to %.1e\n",
            max(-ll[at]), max(sd[is.finite(sd)])))
cat(paste(" ", x[!at | grepl("code 1", x)]), sep = "\n")
cat("\n== 4. wide fits\n")
for (K in c(150, 292)) {
  cat(grep("^seed|^cap", rd(file.path(F4, sprintf("p4wide-%d.txt", K))),
           value = TRUE), sep = "\n")
}

cat("\n== 5. accuracy of the scheme (largest cor over seeds 1:3)\n")
for (job in c("gaussian-0", "gaussian-1", "bernoulli-0", "poisson-0")) {
  x <- rd(file.path(F4, sprintf("acc-%s.txt", job)))
  r <- vapply(6:8, function(dg) {
    l <- grep(sprintf(" deg %d seed", dg), x, value = TRUE)
    e1(max(num(l, ".*cor ([0-9.e+-]+).*")))
  }, "")
  w <- sum(!grepl("warn -$", grep("seed", x, value = TRUE)))
  cat(sprintf("%-12s deg6 %s deg7 %s deg8 %s | warned %d\n", job, r[1],
              r[2], r[3], w))
}
cat("\n== 3. reviewer's cases\n")
for (f in c("asym-1-lane", "asym-0.2-lane", "asym-0.05-lane",
            "asym-0.2-lane-off", "asym-0.05-lane-off", "esc2-lane",
            "p2tier-lane", "hess-poisson6-lane", "hess-bernoulli7-lane",
            "hess-poisson8-lane")) {
  cat("--", f, "\n")
  cat(rd(file.path(F4, paste0(f, ".txt")))[-1], sep = "\n")
}
for (f in c("hess-poisson6-base", "hess-bernoulli7-base",
            "hess-poisson8-base")) {
  cat("--", f, "(0.71.0, written 2026-10-14)\n")
  cat(rd(file.path(D, "fx", paste0(f, ".txt")))[-1], sep = "\n")
}
